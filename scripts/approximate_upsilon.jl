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

using LaTeXStrings
using Plots
pyplot()


species = config["species"]
transition = config["transition"]
lam0 = lam0s[species*transition]
model = config["model"]
pars = convert_dict(config["parameters"], config["model"])

mol = mols[species*transition]


# Plot the vertical profile of Upsilon for a bunch of radii.

# Set up the grid interpolators
nr=128
rmax = 500 * AU
nz=128
zmax = 200 * AU


# choose some radii we want to sample
zs = linspace(0., 4 * AU, nz)
Upsilons = Float64[DiskTracer.model.Upsilon_nu(1. * AU, z, mol.nu_0, pars, mol) for z in zs]
p1 = plot(zs./AU, Upsilons, ylabel=L"$\Upsilon$", xlabel=L"$r_\mathrm{cyl}$")

zs = linspace(0., 10 * AU, nz)
Upsilons = Float64[DiskTracer.model.Upsilon_nu(10. * AU, z, mol.nu_0, pars, mol) for z in zs]
p2 = plot(zs./AU, Upsilons, ylabel=L"$\Upsilon$", xlabel=L"$r_\mathrm{cyl}$")

zs = linspace(0., 40 * AU, nz)
Upsilons = Float64[DiskTracer.model.Upsilon_nu(100. * AU, z, mol.nu_0, pars, mol) for z in zs]
p3 = plot(zs./AU, Upsilons, ylabel=L"$\Upsilon$", xlabel=L"$r_\mathrm{cyl}$")

zs = linspace(0., 150 * AU, nz)
Upsilons = Float64[DiskTracer.model.Upsilon_nu(300. * AU, z, mol.nu_0, pars, mol) for z in zs]
p4 = plot(zs./AU, Upsilons, ylabel=L"$\Upsilon$", xlabel=L"$r_\mathrm{cyl}$")

p = plot(p1,p2,p3,p4,layout=(4,1), legend=false, size=(800,1200))

savefig(p, "Upsilon-traces.png")
