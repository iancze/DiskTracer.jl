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

using LaTeXStrings
using Plots
# gr()
# plotlyjs()
pyplot()

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
interp_S = DiskTracer.model.get_S_nu_grid(pars, mol, mol.nu_0, nr, rmax, nz, zmax)
interp_U = DiskTracer.model.get_Upsilon_nu_grid(pars, mol, mol.nu_0, nr, rmax, nz, zmax)
interp_D2 = DiskTracer.model.get_DeltaV2_grid(pars, mol, nr, rmax, nz, zmax)

# Make a 2D plot of the source function and Upsilon
# density structure
function plot_S(pars::AbstractParameters)

    nr = 128
    nz = 128
    rs = linspace(0.1, 600 * AU, nr)
    zs = linspace(0.1, 300 * AU, nz)

    rr = rs./AU
    zz = zs./AU

    Ss = Array{Float64}(nz, nr)

    for i=1:nz
        for j=1:nr
            Ss[i,j] = interp_S(rs[j], zs[i])
        end
    end

    p = heatmap(Ss, xlabel=L"$r$ [AU]", ylabel=L"$z$ [AU]") #, size=(800,800))
    savefig(p, "source_function_interp.png")


end

function plot_Upsilon(pars::AbstractParameters)

    nr = 128
    nz = 128
    rs = linspace(0.1, 600 * AU, nr)
    zs = linspace(0.1, 300 * AU, nz)


    Us = Array{Float64}(nz, nr)

    for i=1:nz
        for j=1:nr
            Us[i,j] = interp_U(rs[j], zs[i])
        end
    end

    p = heatmap(Us, xlabel=L"$r$ [AU]", ylabel=L"$z$ [AU]") #, size=(800,800))
    savefig(p, "upsilon_interp.png")


end

function plot_DeltaV2(pars::AbstractParameters)

    nr = 128
    nz = 128
    rs = linspace(0.1, 600 * AU, nr)
    zs = linspace(0.1, 300 * AU, nz)


    D2s = Array{Float64}(nz, nr)

    for i=1:nz
        for j=1:nr
            D2s[i,j] = interp_D2(rs[j], zs[i])
        end
    end

    p = heatmap(D2s, xlabel=L"$r$ [AU]", ylabel=L"$z$ [AU]") #, size=(800,800))
    savefig(p, "DeltaV2_interp.png")


end

plot_S(pars)
plot_Upsilon(pars)
plot_DeltaV2(pars)
