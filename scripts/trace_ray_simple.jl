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
    Us = Array{Float64}(nzp)
    Deltas = Array{Float64}(nzp)
    alphas = Array{Float64}(nzp)

    z0s = DiskTracer.geometry.get_zps(xprime, yprime, pars, v0)
    println("z0s ", z0s)


    for i=1:nzp
        zprime = zps[i]
        rcyl = DiskTracer.geometry.get_r_cyl(xprime, yprime, zprime, pars.incl)
        z = DiskTracer.geometry.get_z(xprime, yprime, zprime, pars.incl)

        rcyls[i] = rcyl
        zs[i] = z

        Deltav = v0 * 1e5 - DiskTracer.geometry.get_vlos(xprime, yprime, zprime, pars)
        DeltaV = sqrt(2 * kB * DiskTracer.model.temperature(rcyl, pars)/mol.mol_weight + (pars.ksi * 1e5)^2)

        Deltas[i] = exp(-Deltav^2/DeltaV^2)

        Us[i] = DiskTracer.model.Upsilon_nu(rcyl, z, mol.nu_0, pars, mol)
        alphas[i] = Us[i] * exp(-Deltav^2/DeltaV^2)

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
    zpts, sol = DiskTracer.trace.trace_pixel(xprime, yprime, v0, mol, pars, interp)
    println("end tau, I: ", sol.u[end])
    println("fancy steps: ", length(zpts))

    # Also call our custom integrator to get tau
    zpts_c, tau_c, tot_i = DiskTracer.trace_simple.trace_pixel(xprime, yprime, v0, mol, pars, interp)
    println("num iterations: ", length(zpts_c))
    println("tot intensity: ", tot_i)
    # println("zpts_c: ", zpts_c)
    # println("tau_c: ", tau_c)

    p1 = plot(zps./AU, Us, xlabel=L"$z^\prime$", ylabel=L"$\Upsilon$") #, yaxis=(scale=:log10))


    p2 = plot(zps./AU, alphas, xlabel=L"$z^\prime$", ylabel=L"$\alpha$")
    # vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p3 = plot(zps./AU, alphas, xlabel=L"$z^\prime$", ylabel=L"$\alpha$", yaxis=(scale=:log10))
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))


    p4 = plot(zpts_c./AU, tau_c, xlim=(200, 400), ylabel=L"$\tau$")
    plot!(zpts./AU, sol[1,:], xlim=(200, 400), ylabel=L"$\tau$")

    # p6 = plot(zpts./AU, sol[2,:], xlim=(-600, 600), ylabel=L"$I$")
    # p6 = plot(zpts./AU, zpts_inv./AU, xlim=(-600, 600))

    # println("final ", sol[2,end])

    p = plot(p1,p2,p3,p4, layout=(4,1), legend=false, size=(800,1200))


    savefig(p, "ray_trace_simple.png")


end


plot_ray(100.*AU, -100.0*AU)
