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

using LaTeXStrings
using Plots
pyplot()

# Make a 2D plot of the source function and Upsilon
# density structure
function plot_S(pars::AbstractParameters)

    nr = 128
    nz = 128
    rs = linspace(0.1, 600 * AU, nr)
    zs = linspace(0.1, 600 * AU, nz)

    rr = rs./AU
    zz = zs./AU

    Ss = Array{Float64}(nz, nr)

    for i=1:nz
        for j=1:nr
            Ss[i,j] = DiskTracer.model.S_nu(rs[j], zs[i], mol.nu_0, pars)
        end
    end

    p = heatmap(log10.(Ss), xlabel=L"$r$ [AU]", ylabel=L"$z$ [AU]") #, size=(800,800))
    savefig(p, "source_function.png")
end

function plot_Upsilon(pars::AbstractParameters)

    nr = 128
    nz = 128
    rs = linspace(0.1, 200 * AU, nr)
    zs = linspace(0.1, 200 * AU, nz)

    rr = rs./AU
    zz = zs./AU

    Us = Array{Float64}(nz, nr)

    for i=1:nz
        for j=1:nr
            Us[i,j] = DiskTracer.model.Upsilon_nu(rs[j], zs[i], mol.nu_0, pars, mol)
        end
    end

    # show 5 decades of scale
    p = heatmap(log10.(Us+ 1e-5 * maximum(Us)), xlabel=L"$r$ [AU]", ylabel=L"$z$ [AU]")
    savefig(p, "Upsilon.png")

end


plot_S(pars)
plot_Upsilon(pars)
