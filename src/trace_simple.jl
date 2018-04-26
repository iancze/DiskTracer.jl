"Module for tracing a ray."
module trace_simple

export trace_pixel

using ..model
using ..constants
using ..geometry

# Simply use midpoint formula and adaptive steps to measure tau


const TAU_THRESH = 8.0

"Adaptive stepper using Midpoint Method. https://en.wikipedia.org/wiki/Midpoint_method "
function integrate_tau(zstart::Real, v0::Real, pars::AbstractParameters, interp, gpre; h_tau::Real=0.1, nsteps::Int=1000, tau_start::Real=0.0, max_ds::Real=-5*AU)
    # Initialize an array of taus
    zps = Array{Float64}(nsteps)
    alpha_mps = Array{Float64}(nsteps) # alpha evaluated at the midpoints
    taus = Array{Float64}(nsteps)

    tot_intensity = 0.0

    zps[1] = zstart
    taus[1] = tau_start

    # Write in coordinate conversions for querying alpha and sfunc
    rcyl, z, vlos = get_coords(zstart, gpre)

    # Calculate Delta v at this position
    Deltav = v0 * 1.e5 - vlos # get_vlos(xprime, yprime, zprime, pars)

    # Look up RT quantities from nearest-neighbor interp.
    DV2, sfunc, Upsilon = interp(rcyl, z)

    # Evaluate alpha at current zstart position
    alpha = Upsilon * exp(-Deltav^2/DV2)
    alpha_mps[1] = alpha
    zp2 = 0.0

    println("Max ds: ", max_ds)

    for i=2:nsteps
        # Based upon current alpha value, calculate how much of a dz we would need to get dz * alpha = h_tau
        # Because these are negative numbers, the smaller step is actually the maximum (closer to 0)
        dzp = maximum([-h_tau / alpha_mps[i-1], max_ds])

        # Update current z position
        zps[i] = dzp + zps[i-1]

        # Midpoint step position
        zp2 = zps[i-1] + 0.5 * dzp

        # Write in coordinate conversions for querying alpha and sfunc
        rcyl, z, vlos = get_coords(zp2, gpre)

        # Calculate Delta v at this position
        Deltav = v0 * 1.e5 - vlos # get_vlos(xprime, yprime, zprime, pars)

        # Look up RT quantities from nearest-neighbor interp.
        # sfunc is evaluated at the midpoint
        DV2, sfunc, Upsilon = interp(rcyl, z)

        # Calculate alpha at new midpoint location
        alpha_mps[i] = Upsilon * exp(-Deltav^2/DV2)
        taus[i] = taus[i-1] - dzp * alpha_mps[i]

        # dtau = - alpha * dzp
        # tau_(n+1) = tau_n + dzp * alpha(zp_n + h/2)
        # since we have a simple ODE, explicit and implicit techniques are the same.

        # Calculate S at the new midpoint location
        # tot_intensity += -exp(-0.5(taus[i] + taus[i-1])) * sfunc * alpha_mps[i] * dzp # why is this the most accurate?
        # tot_intensity += -exp(-taus[i]) * sfunc * alpha_mps[i] * dzp
        # the reason it wasn't working before is because  sfunc is actually already evaluated at the midpoint.
        tot_intensity += (exp(-taus[i-1]) - exp(-taus[i])) * sfunc # this is also very accurate

        # Update trailing sfuncs
        # sfunc_ = sfunc

        # Query alpha at current z position
        if taus[i] > TAU_THRESH
            # total intensity
            return (zps[1:i], taus[1:i], tot_intensity)
        end

    end

    return (zps, taus, tot_intensity)
end



"v0 is the central velocity of the channel, relative to the disk systemic velocity."
function trace_pixel(xprime, yprime, v0, mol::Molecule, pars::AbstractParameters, interp)

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

    gpre = precalc_geo(xprime, yprime, pars)

    ans = integrate_tau(z1start, v0, pars, interp, gpre)

    return ans

    # if (tau_final < TAU_THRESH) & (length(bounds) == 4)
    #     #We didn't reach the threshold in the first bounding region, so continue in the second.
    #     tspan = (0, z2start - z2end)
    #
    #     # Start the ode to integrate the first path from the end of the second.
    #     u0 = Float64[tau_final; int_final]
    #     args = (z2start, v0, pars, interp, gpre)
    #     prob = ODEProblem(f, u0, tspan, args)
    #     sol = solve(prob, alg, callback=cb, reltol=rt, abstol=at, dense=false, save_everystep=false, dtmax=dtmax, dt=dt)
    #     tau_final, int_final = sol.u[end]
    #     return int_final
    #
    # else
    #     return int_final
    # end
end



end # module
