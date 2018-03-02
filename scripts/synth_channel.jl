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
using DiskTracer.trace

species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]
pars = convert_dict(config["parameters"], config["model"])

mol = mols[species*transition]

using NPZ

# Set up the grid interpolators
nr=128
rmax = 500 * AU
nz=128
zmax = 200 * AU

# Get the interpolator object for this particular setup.
interp = DiskTracer.model.get_grids(pars, mol, mol.nu_0, nr, rmax, nz, zmax)

npix = 128
xs = linspace(-500 * AU, 500 * AU, npix)
ys = linspace(-500 * AU, 500 * AU, npix)

nv = 9
vs = linspace(-4.0, 4.0, nv)

img = zeros(Float64, npix, npix, nv)
mask = Array{Bool}(npix, npix, nv)


tic()
# Fill the boolean mask and trace image
for k=1:nv
    for j=1:npix
        for i=1:npix
            # verify_pixel(xs[i], ys[j], pars, vs[k], 0.5, 400*AU)
            yes = verify_pixel(xs[i], ys[j], pars, vs[k], 0.5, 400*AU)
            mask[j,i,k] = yes
            if yes
                img[j,i,k] = trace_pixel(xs[i], ys[j], vs[k], mol, pars, interp)
            end
        end
    end
end
toc()

# Save the channel maps and mask

npzwrite("img.npy", img)
npzwrite("mask.npy", mask)
