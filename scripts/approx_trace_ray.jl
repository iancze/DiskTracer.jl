#!/usr/bin/env julia

using ArgParse

s = ArgParseSettings()
@add_arg_table s begin
    "-M"
    help = "Print out the list of files this program will create."
    action = :store_true
    "config"
    help = "a YAML configuration file"
    default = "config.yaml"
end

parsed_args = parse_args(ARGS, s)

import YAML
config = YAML.load(open(parsed_args["config"]))

using DiskTracer.model
using DiskTracer.constants
using DiskTracer.geometry
using SpecialFunctions

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]
pars = convert_dict(config["parameters"], config["model"])

mol = mols[species*transition]

nr=128
rmax = 500 * AU
nz=128
zmax = 200 * AU
# Get the interpolator object for this particular setup.
interp = DiskTracer.model.get_grids(pars, mol, mol.nu_0, nr, rmax, nz, zmax)


using LaTeXStrings
using Plots
pyplot()

# pick a pixel (xprime, yprime), and trace a full ray through the disk

function plot_ray(xprime, yprime)

    v0 = 1.0 # km/s

    # Get a decent estimate of DeltaVmax
    # Calculate DeltaVmax for this ray
    rcyl_min = abs(xprime)
    DeltaVmax = sqrt(2 * kB * temperature(rcyl_min, pars)/mol.mol_weight + (pars.ksi * 1e5)^2) * 1e-5 # km/s

    # DeltaVmax = 0.5 # km/s
    rmax = 600 * AU
    nzp = 512
    zps = linspace(-600 * AU, 600 * AU, nzp)

    rcyls = Array{Float64}(nzp)
    zs = Array{Float64}(nzp)
    Ss = Array{Float64}(nzp)
    Deltavs = Array{Float64}(nzp)
    tildeDeltavs = Array{Float64}(nzp)
    DeltaVs = Array{Float64}(nzp)
    Us = Array{Float64}(nzp)
    tildeUs = Array{Float64}(nzp)
    Deltas = Array{Float64}(nzp)
    alphas = Array{Float64}(nzp)
    tildealphas = Array{Float64}(nzp)
    tildealphasG = Array{Float64}(nzp)

    z0s = DiskTracer.geometry.get_zps(xprime, yprime, pars, v0)
    println("z0s ", z0s)

    # Get the minimum rcyl
    rcyl_min = xprime

    # Occurs at what value of zprime?
    zprime_rcyl_min = cosd(pars.incl)/sind(pars.incl) * yprime

    # Actually, we want the zprime that corresponds to the midplane crossing
    zprime_midplane = -tand(pars.incl) * yprime
    # Calculate the rcyl_mid here

    # Get amplitude at midplane crossing
    # Corresponds to r_cyl_mid, 0.0
    r_cyl_mid = sqrt(xprime^2 + (cosd(pars.incl) * yprime + sind(pars.incl) * zprime_midplane)^2)

    # println("rcyl min ", rcyl_min/AU, " zprime_rcyl_min ", zprime_rcyl_min/AU, " zprime_midplane ", zprime_midplane/AU, " rcyl_mid ", r_cyl_mid/AU)

    # Get scale height at rcyl_min (will be an underestimate)
    H_rcyl_min = DiskTracer.model.Hp(rcyl_min, pars)
    H_rcyl_mid = DiskTracer.model.Hp(r_cyl_mid, pars)

    DeltaV = sqrt(2 * kB * DiskTracer.model.temperature(rcyl_min, pars)/mol.mol_weight + (pars.ksi * 1e5)^2)

    sigma_Upsilon = 10 * H_rcyl_mid / cosd(pars.incl)
    a_Upsilon = Upsilon_nu(r_cyl_mid, 0.0, mol.nu_0, pars, mol)
    mu_Upsilon = zprime_midplane

    DeltaV = sqrt(2 * kB * DiskTracer.model.temperature(rcyl_min, pars)/mol.mol_weight + (pars.ksi * 1e5)^2)

    mu_upsilon1, mu_upsilon2 = z0s
    sigma_upsilon = DeltaV * sqrt(2)/3 * (xprime * sqrt(G * pars.M_star * M_sun) * sind(pars.incl))^(4./3) / (mu_upsilon1 * (v0 * 1.e5)^(7/3))
    # a_upsilon = sqrt(2pi) * sigma_upsilon

    sigma_alpha2 = 1/(sigma_Upsilon^(-2) + sigma_upsilon^(-2)) # actually sigma^2, hence the 2 in the name
    mu_alpha1 = sigma_alpha2 * (mu_Upsilon/sigma_Upsilon^2 + mu_upsilon1/sigma_upsilon^2)
    mu_alpha2 = sigma_alpha2 * (mu_Upsilon/sigma_Upsilon^2 + mu_upsilon2/sigma_upsilon^2)

    a_alpha1 = a_Upsilon * exp(- (mu_Upsilon - mu_upsilon1)^2/(2 * (sigma_Upsilon^2 + sigma_upsilon^2)))
    a_alpha2 = a_Upsilon * exp(- (mu_Upsilon - mu_upsilon2)^2/(2 * (sigma_Upsilon^2 + sigma_upsilon^2)))


    for i=1:nzp
        zprime = zps[i]
        rcyl = DiskTracer.geometry.get_r_cyl(xprime, yprime, zprime, pars.incl)
        z = DiskTracer.geometry.get_z(xprime, yprime, zprime, pars.incl)

        rcyls[i] = rcyl
        zs[i] = z

        Deltav = v0 * 1e5 - DiskTracer.geometry.get_vlos(xprime, yprime, zprime, pars)
        Deltavs[i] = Deltav

        if zprime > 0
            tildeDeltav = DiskTracer.geometry.get_Delta_v_approx(xprime, yprime, zprime, pars, v0, z0s[1])
        else
            tildeDeltav = DiskTracer.geometry.get_Delta_v_approx(xprime, yprime, zprime, pars, v0, z0s[2])
        end
        tildeDeltavs[i] = tildeDeltav * 1e-5 #km/s

        DeltaV = sqrt(2 * kB * DiskTracer.model.temperature(rcyl, pars)/mol.mol_weight + (pars.ksi * 1e5)^2)
        DeltaVs[i] = DeltaV * 1e-5 # km/s

        Deltas[i] = exp(-Deltav^2/DeltaV^2)

        Us[i] = DiskTracer.model.Upsilon_nu(rcyl, z, mol.nu_0, pars, mol)
        tildeUs[i] = a_Upsilon * exp(-0.5 * (zprime - zprime_midplane)^2/sigma_Upsilon^2)

        alphas[i] = Us[i] * exp(-Deltav^2/DeltaV^2)
        # tildealphas[i] = tildeUs[i] * exp(-tildeDeltavs[i]^2 / DeltaVs[i]^2)

        if zprime > 0
            tildealphasG[i] = a_alpha1 * exp(-0.5 * (zprime - mu_alpha1)^2/sigma_alpha2)
        else
            tildealphasG[i] = a_alpha2 * exp(-0.5 * (zprime - mu_alpha2)^2/sigma_alpha2)
        end
    end

    # get bounding zps
    zbounds = DiskTracer.geometry.get_bounding_zps(xprime, yprime, pars, v0, DeltaVmax, rmax)
    if length(zbounds) == 2
        zstart, zend = zbounds
        println("H bound verified merged: ", DiskTracer.geometry.verify_h_bound(xprime, yprime, zstart, zend, pars, 0.5))
    else
        z1start, z1end, z2start, z2end = zbounds
        println("H bound verified front-side: ", DiskTracer.geometry.verify_h_bound(xprime, yprime, z1start, z1end, pars, 0.5))
        println("H bound verified rear-side: ", DiskTracer.geometry.verify_h_bound(xprime, yprime, z2start, z2end, pars, 0.5))
    end

    # Call the routine that actually traces the ray, and returns zp, tau, and I quantities
    sol = DiskTracer.trace.trace_pixel(xprime, yprime, v0, mol, pars, interp)
    zpts = sol.t


    # Also call our custom integrator to get tau
    zpts_c, tau_c = DiskTracer.trace_simple.trace_pixel(xprime, yprime, v0, mol, pars, interp)


    # tildetau = Array{Float64}(length(zpts))
    # for i=1:length(zpts)
    #     z = zpts[i]
    #     tildetau[i] = a_alpha1 * erfc((z - mu_c1)/(sqrt(2 * sigma_c2)))/(2 * sqrt(pi))
    # end
    #
    #
    # "Given tilde_tau, invert to find zprime"
    # function tilde_tau_inv(tilde_tau)
    #     zprime = sqrt(2 * sigma_c2) * erfcinv(tilde_tau * 2 * sqrt(pi)/a_alpha1) + mu_c1
    # end
    #
    # zpts_inv = Array{Float64}(length(zpts))
    # for i=length(zpts)
    #     zpts_inv[i] = tilde_tau_inv(tildetau[i])
    # end

    p1 = plot(zps./AU, Us, xlabel=L"$z^\prime$", ylabel=L"$\Upsilon$") #, yaxis=(scale=:log10))

    p2 = plot(zps./AU, tildeUs, ylabel=L"$\tilde{\Upsilon}$") #, yaxis=(scale=:log10))
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p3 = plot(zps./AU, alphas, xlabel=L"$z^\prime$", ylabel=L"$\alpha$")
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))
    # plot!(zps./AU, tildealphas, ylabel=L"$\tilde{\alpha}$")
    plot!(zps./AU, tildealphasG, ylabel=L"$\tilde{\alpha}$")

    p4 = plot(zps./AU, alphas, xlabel=L"$z^\prime$", ylabel=L"$\alpha$", yaxis=(scale=:log10))
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))
    # plot!(zps./AU, tildealphas, ylabel=L"$\tilde{\alpha}$")
    plot!(zps./AU, tildealphasG, ylabel=L"$\tilde{\alpha}$", yaxis=(scale=:log10))

    # p5 = plot(zpts./AU, sol[1,:], xlim=(-600, 600), ylabel=L"$\tau$")
    p55 = plot(zpts_c./AU, tau_c, xlim=(-600, 600), ylabel=L"$\tau$")
    # plot!(zpts./AU, tildetau, xlim=(-600, 600))

    # p6 = plot(zpts./AU, sol[2,:], xlim=(-600, 600), ylabel=L"$I$")
    # p6 = plot(zpts./AU, zpts_inv./AU, xlim=(-600, 600))

    # println("final ", sol[2,end])

    p = plot(p1,p2,p3,p4,p55, layout=(7,1), legend=false, size=(800,1200))


    savefig(p, "approx_ray_trace.png")

    # println(DiskTracer.trace.trace_pixel_tilde(xprime, yprime, v0, mol, pars, interp))

end


plot_ray(100.*AU, -100.0*AU)
