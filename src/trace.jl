"Module for tracing a ray."
module trace

export trace_pixel

using ..model
using ..constants
using ..geometry

using DifferentialEquations
using SpecialFunctions

# Some of setups for the integrator.
# const dt = 1e-3 * AU
alg = Tsit5()
# const rt=1e-6
# const at=1e-8

const rt=Float64[1e-8, 1e-12]
const at=Float64[1e-8, 1e-12]
# const dtmax = 5 * AU

const tau_thresh = 8.0

# Set up the terminator callback, stop ODE integration when this condition is satisfied.
function condition(u, t, integrator)
    u[1] > tau_thresh # when tau = 8
end
affect!(integrator) = terminate!(integrator)

cb = DiscreteCallback(condition, affect!)


# The integrand. p encapsulates the passed parameters.
function f(du, u, args, t)

    zstart, v0, pars, interp, gpre = args

    zprime = zstart - t

    # Write in coordinate conversions for querying alpha and sfunc
    rcyl, z, vlos = get_coords(zprime, gpre)

    # Calculate Delta v at this position
    Deltav = v0 * 1.e5 - vlos # get_vlos(xprime, yprime, zprime, pars)

    # Look up RT quantities from nearest-neighbor interp.
    DV2, sfunc, Upsilon = interp(rcyl, z)

    # Calculate alpha
    alpha = Upsilon * exp(-Deltav^2/DV2)

    # Update dtau/ds and dI/ds
    du[1] = alpha
    du[2] = exp(-u[1]) * sfunc * alpha
end


"v0 is the central velocity of the channel, relative to the disk systemic velocity."
function trace_pixel(xprime, yprime, v0, mol::Molecule, pars::AbstractParameters, interp; dt=0.5 * AU, dtmax=5 * AU)

    rmax = 700.0 * AU

    # Calculate DeltaVmax for this ray
    rcyl_min = abs(xprime)
    DeltaVmax = sqrt(2 * kB * temperature(rcyl_min, pars)/mol.mol_weight + (pars.ksi * 1e5)^2) * 1e-5 # km/s

    # Calculate the bounding positions for emission
    # If it's two, just trace it.
    # If it's four, break it up into two integrals.

    bounds = get_bounding_zps(xprime, yprime, pars, v0, DeltaVmax, rmax)
    if length(bounds) == 4
        z1start, z1end, z2start, z2end = bounds
    elseif length(bounds) == 2
        z1start, z1end = bounds
    end

    # Start the ode to integrate the first path.
    u0 = Float64[0.0; 0.0]

    gpre = precalc_geo(xprime, yprime, pars)

    args = (z1start, v0, pars, interp, gpre)

    tspan = (0.0, z1start - z1end)

    prob = ODEProblem(f, u0, tspan, args)

    sol = solve(prob, alg, callback=cb, reltol=rt, abstol=at, dense=true, save_everystep=true, dtmax=dtmax, dt=dt)
    # println(sol)

    tau_final, int_final = sol.u[end]dtmax

    zps = z1start - sol.t

    # println(sol.t ./ AU)
    # println("Number of evaluations ", length(sol.t))

    if (tau_final < tau_thresh) & (length(bounds) == 4)
        #We didn't reach the threshold in the first bounding region, so continue in the second.
        tspan = (0, z2start - z2end)

        # Start the ode to integrate the first path from the end of the second.
        u0 = Float64[tau_final; int_final]
        args = (z2start, v0, pars, interp, gpre)
        prob = ODEProblem(f, u0, tspan, args)
        sol = solve(prob, alg, callback=cb, reltol=rt, abstol=at, dense=true, save_everystep=true, dtmax=dtmax, dt=dt)
        tau_final, int_final = sol.u[end]
        zps = z2start - sol.t
        # return int_final
        return (zps, sol)

    else
        # return int_final
        return (zps, sol)
    end
end



end # module
