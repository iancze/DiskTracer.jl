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



# DeltaV max
println("Delta V max, ", sqrt(DeltaV2(0.1 * AU, 0.0, pars, mol)) * 1e-5)

# For a set of velocities, make channel maps to test whether we can mask appropriately.

# V4046 Sgr has: -3 to + 3 kms, in 0.4 km/s chunks.
nv = 22
vs = linspace(-4.0, 4.0, nv)

npix = 128
xs = linspace(-500 * AU, 500 * AU, npix)
ys = linspace(-500 * AU, 500 * AU, npix)


arr = Array{Bool}(npix, npix, nv)

# Fill the boolean mask
for k=1:nv
    for j=1:npix
        for i=1:npix
            arr[j,i,k] = verify_pixel(xs[i], ys[j], pars, vs[k], 0.5, 400*AU)
        end
    end
end

println("Fraction traced: ",sum(arr)/(npix * npix * nv))

function plot_channel_maps(arr)
    parr = [heatmap(arr[:,:,k]) for k=1:nv]
    p = plot(parr..., layout=(2, 11), size=(11 * 250,500))
    savefig(p, "mask.png")
end

plot_channel_maps(arr)
