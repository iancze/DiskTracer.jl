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

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]
pars = convert_dict(config["parameters"], config["model"])

mol = mols[species*transition]

import PyPlot.plt
using LaTeXStrings


# Make a 2D plot of the source function and Upsilon
# density structure
function plot_S(pars::AbstractParameters)

    nr = 128
    nz = 128
    rs = linspace(0.1, 600 * AU, nr)
    zs = linspace(0.1, 600 * AU, nz)

    rr = rs./AU
    zz = zs./AU

    xx = Array{Float64}(nz, nr)
    yy = Array{Float64}(nz, nr)

    Ss = Array{Float64}(nz, nr)

    for i=1:nz
        xx[i, :] = rr
    end

    for j=1:nr
        yy[:, j] = zz
    end

    for i=1:nz
        for j=1:nr
            Ss[i,j] = DiskTracer.model.S_nu(rs[j], zs[i], mol.nu_0, pars)
        end
    end


    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    ax[:set_ylabel](L"$z$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")


    img = ax[:contourf](xx, yy, log10.(Ss)) #
    # img = ax[:contourf](xx, yy, Ss, levels=levels)

    # ax[:set_xscale]("log")


    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.77)

    cax = fig[:add_axes]([0.82, 0.22, 0.03, 0.65])
    cb = fig[:colorbar](img, cax=cax)

    plt[:savefig]("source_function.png")
end

function plot_Upsilon(pars::AbstractParameters)

    nr = 128
    nz = 128
    rs = linspace(0.1, 600 * AU, nr)
    zs = linspace(0.1, 600 * AU, nz)

    rr = rs./AU
    zz = zs./AU

    xx = Array{Float64}(nz, nr)
    yy = Array{Float64}(nz, nr)

    Us = Array{Float64}(nz, nr)

    for i=1:nz
        xx[i, :] = rr
    end

    for j=1:nr
        yy[:, j] = zz
    end

    for i=1:nz
        for j=1:nr
            Us[i,j] = DiskTracer.model.Upsilon(rs[j], zs[i], mol.nu_0, pars, mol)
        end
    end


    fig = plt[:figure]()
    ax = fig[:add_subplot](111)

    ax[:set_ylabel](L"$z$ [AU]")
    ax[:set_xlabel](L"$r$ [AU]")

    ext = (rr[1], rr[end], zz[1], zz[end])
    img = ax[:imshow](log10.(Us), vmin=-20, origin="lower", extent=ext) #
    # img = ax[:contourf](xx, yy, Ss, levels=levels)

    # ax[:set_xscale]("log")


    fig[:subplots_adjust](left=0.15, bottom=0.15, right=0.77)

    cax = fig[:add_axes]([0.82, 0.22, 0.03, 0.65])
    cb = fig[:colorbar](img, cax=cax)

    plt[:savefig]("Upsilon.png")
end


# plot_S(pars)
plot_Upsilon(pars)
