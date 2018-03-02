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


using BenchmarkTools
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


nr=128
rmax = 500 * AU
nz=128
zmax = 200 * AU

# Get the interpolator object for this particular setup.
interp = DiskTracer.model.get_grids(pars, mol, mol.nu_0, nr, rmax, nz, zmax)

# pick a pixel (xprime, yprime), and trace a full ray through the disk

v0 = 2.0 # km/s


I = trace_pixel(100*AU, 80*AU, v0, mol, pars, interp)
println(I)

quit()
# @btime trace_pixel(100*AU, 80*AU, v0, mol, pars, grids)
# @profile trace_pixel(100*AU, 80*AU, v0, mol, pars, grids)
# Profile.print()

DeltaVmax = 0.5 # km/s
rmax = 600*AU #
xprime = 100*AU
yprime = 80*AU
z1start, z1end = geometry.get_bounding_zps(xprime, yprime, pars, v0, DeltaVmax, rmax)
gpre = precalc_geo(xprime, yprime, pars)
args = (z1start, v0, pars, interp, gpre)

zprime = 2*AU

rcyl = DiskTracer.geometry.get_r_cyl(xprime, yprime, zprime, pars.incl)
z = DiskTracer.geometry.get_z(xprime, yprime, zprime, pars.incl)

# Try evaluating f
println("Get vlos")
@btime DiskTracer.geometry.get_vlos(xprime, yprime, zprime, pars)

println("Get r_cyl")
@btime DiskTracer.geometry.get_r_cyl(xprime, yprime, zprime, pars.incl)

println("Get z")
@btime DiskTracer.geometry.get_z(xprime, yprime, zprime, pars.incl)

println("Interp")
@btime interp(rcyl, z)

println("f")
@btime DiskTracer.trace.f(Float64[0.0, 0.0], Float64[0.0, 0.0], args, 2*AU)

# DiskTracer.trace.f(Float64[0.0, 0.0], Float64[0.0, 0.0], args, 2*AU)
# # Profile.init(delay=1e-16)
# @profile [DiskTracer.trace.f(Float64[0.0, 0.0], Float64[0.0, 0.0], args, 2*AU) for i=1:100]
# Profile.print()
