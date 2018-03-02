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

    vlos = 0.5 # km/s
    DeltaVmax = 0.2 # km/s
    rmax = 600 * AU
    nzp = 512
    zps = linspace(-400 * AU, 400 * AU, nzp)

    rcyls = Array{Float64}(nzp)
    zs = Array{Float64}(nzp)
    Ss = Array{Float64}(nzp)
    Us = Array{Float64}(nzp)
    Deltas = Array{Float64}(nzp)
    alphas = Array{Float64}(nzp)

    for i=1:nzp
        zprime = zps[i]
        rcyl = DiskTracer.geometry.get_r_cyl(xprime, yprime, zprime, pars.incl)
        z = DiskTracer.geometry.get_z(xprime, yprime, zprime, pars.incl)

        rcyls[i] = rcyl
        zs[i] = z

        Ss[i] = DiskTracer.model.S_nu(rcyl, z, mol.nu_0, pars)
        Us[i] = DiskTracer.model.Upsilon_nu(rcyl, z, mol.nu_0, pars, mol)

        # Calculate Delta v at this position
        Deltav = vlos * 1e5 - DiskTracer.geometry.get_vlos(xprime, yprime, zprime, pars)

        # Look up DeltaV2 from nearest-neighbor interp.
        DeltaV = sqrt(2 * kB * DiskTracer.model.temperature(rcyl, pars)/mol.mol_weight + (pars.ksi * 1e5)^2)

        Deltas[i] = exp(-Deltav^2/DeltaV^2)

        alphas[i] = Us[i] * exp(-Deltav^2/DeltaV^2)

    end

    # get bounding zps
    zbounds = DiskTracer.geometry.get_bounding_zps(xprime, yprime, pars, vlos, DeltaVmax, rmax)

    p1 = plot(zps./AU, rcyls./AU, ylabel=L"$r_\mathrm{cyl}$ [AU]")
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p2 = plot(zps./AU, zs./AU, ylabel=L"$z$ [AU]")
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p25 = plot(zps./AU, abs.(zs./rcyls), ylabel=L"$|z/r|$")
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p3 = plot(zps./AU, Ss, ylabel=L"$S_\nu$", yaxis=(scale=:log10))
    # vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p4 = plot(zps./AU, Us, xlabel=L"$z^\prime$", ylabel=L"$\Upsilon$")# , yaxis=(scale=:log10))
    # vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p5 = plot(zps./AU, Deltas, xlabel=L"$z^\prime$", ylabel=L"$\exp(\Delta v^2/\Delta V^2)$")# , yaxis=(scale=:log10))
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p6 = plot(zps./AU, alphas, xlabel=L"$z^\prime$", ylabel=L"$\alpha$")
    vline!((zbounds./AU)', line=([:red :green :blue :cyan]))

    p = plot(p1,p2,p25,p3,p4,p5,p6,layout=(7,1), legend=false, size=(800,1000)) #, dpi=200)


    savefig(p, "ray_plot.png")


end

plot_ray(10.*AU, 90.0*AU)
