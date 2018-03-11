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

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]
pars = convert_dict(config["parameters"], config["model"])

mol = mols[species*transition]

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

    z0s = DiskTracer.geometry.get_zps(xprime, yprime, pars, v0)
    println("z0s ", z0s)

    # Get the minimum rcyl
    rcyl_min = xprime

    # Occurs at what value of zprime?
    zprime_rcyl_min = cosd(pars.incl)/sind(pars.incl) * yprime

    # Actually, we want the zprime that corresponds to the midplane crossing
    zprime_midplane = -tand(pars.incl) * yprime

    println("rcyl min ", rcyl_min/AU, " zprime_rcyl_min ", zprime_rcyl_min/AU, " zprime_midplane ", zprime_midplane/AU)

    # Get scale height at rcyl_min
    H_rcyl_min = DiskTracer.model.Hp(rcyl_min, pars)

    # The amplitude of tilde Upsilon is just the value at the midplane
    Upsilon_amp = DiskTracer.model.Upsilon_nu(rcyl_min, 0.0, mol.nu_0, pars, mol)

    # The sigma of the Gaussian for tilde Upsilon is this lengthened by the inclination
    Upsilon_sigma = H_rcyl_min / cosd(pars.incl)

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


        # Ss[i] = DiskTracer.model.S_nu(rcyl, z, mol.nu_0, pars)

        Us[i] = DiskTracer.model.Upsilon_nu(rcyl, z, mol.nu_0, pars, mol)
        tildeUs[i] = Upsilon_amp * exp(-0.5 * (zprime - zprime_midplane)^2/Upsilon_sigma^2)

        # println(Us[i], " ", tildeUs[i])

        alphas[i] = Us[i] * exp(-Deltav^2/DeltaV^2)
        # tilde_alphas[i] =

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

    # hbounds, derivs = DiskTracer.geometry.get_bounding_scale_heights(xprime, yprime, pars, 1.0)


    p1 = plot(zps./AU, rcyls./AU, ylabel=L"$r_\mathrm{cyl}$ [AU]")
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p2 = plot(zps./AU, zs./AU, ylabel=L"$z$ [AU]")
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p3 = plot(zps./AU, Deltavs .* 1e-5, ylabel=L"$\Delta v$")
    plot!(zps./AU, tildeDeltavs, ylims=(-1., 1))

    p4 = plot(zps./AU, Deltas, xlabel=L"$z^\prime$", ylabel=L"$\exp(-\Delta v^2/\Delta V^2)$")
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))
    vline!((z0s./AU)', line=([:black :black]))

    println(minimum(Us))
    println(minimum(tildeUs))

    p5 = plot(zps./AU, Us, xlabel=L"$z^\prime$", ylabel=L"$\Upsilon$") #, yaxis=(scale=:log10))
    p55 = plot(zps./AU, tildeUs, ylabel=L"$\tilde{\Upsilon}$") #, yaxis=(scale=:log10))
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    # p45 = plot(zps./AU, DeltaVs, xlabel=L"$z^\prime$", ylabel=L"$\Delta V$")

    p6 = plot(zps./AU, alphas, xlabel=L"$z^\prime$", ylabel=L"$\alpha$")
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p = plot(p1,p2,p3,p4,p5,p55,p6,layout=(7,1), legend=false, size=(800,1200))


    savefig(p, "approx_ray_plot.png")


end

plot_ray(100.*AU, 100.0*AU)
