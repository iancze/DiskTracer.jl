"Module for tracing a ray."
module trace

export trace_pixel, trace_pixel_debug

using ..model
using ..constants
using ..geometry

using DifferentialEquations
using SpecialFunctions

# Some of setups for the integrator.
const dt = 1e-3 * AU
alg = Tsit5()
const rt=1e-8
const at=1e-10

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
function trace_pixel(xprime, yprime, v0, mol::Molecule, pars::AbstractParameters, interp; get_tau=false)

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

    sol = solve(prob, alg, callback=cb, reltol=rt, abstol=at, dense=false, save_everystep=false, dtmax=5*AU, dt=dt)
    # println(sol)

    tau_final, int_final = sol.u[end]

    if (tau_final < tau_thresh) & (length(bounds) == 4)
        #We didn't reach the threshold in the first bounding region, so continue in the second.
        tspan = (0, z2start - z2end)

        # Start the ode to integrate the first path from the end of the second.
        u0 = Float64[tau_final; int_final]
        args = (z2start, v0, pars, interp, gpre)
        prob = ODEProblem(f, u0, tspan, args)
        sol = solve(prob, alg, callback=cb, reltol=rt, abstol=at, dense=false, save_everystep=false, dtmax=3*AU, dt=dt)
        tau_final, int_final = sol.u[end]

        if get_tau
            return int_final, tau_final
        else
            return int_final
        end
    else
        if get_tau
            return int_final, tau_final
        else
            return int_final
        end
    end
end

# The integrand. p encapsulates the passed parameters.
"The idea here is to have a function which steps in dtau, not ds."
function f_tilde(du, u, args, t)

    # tilde_alpha is a function
    zstart, v0, pars, interp, gpre, tilde_alpha_func, tilde_tau_inv = args

    # t is now tilde_tau, not s
    zprime = tilde_tau_inv(t)
    println("tilde tau ", t, " tau ", u[1], " zprime ", zprime/AU)

    # Write in coordinate conversions for querying alpha and sfunc
    rcyl, z, vlos = get_coords(zprime, gpre)

    # Calculate Delta v at this position
    Deltav = v0 * 1.e5 - vlos # get_vlos(xprime, yprime, zprime, pars)

    # Look up RT quantities from nearest-neighbor interp.
    DV2, sfunc, Upsilon = interp(rcyl, z)

    # Calculate alpha
    alpha = Upsilon * exp(-Deltav^2/DV2)

    tilde_alpha = tilde_alpha_func(zprime)

    # Update dtau/ds and dI/ds
    du[1] = alpha/tilde_alpha
    du[2] = exp(-u[1]) * sfunc * alpha/tilde_alpha
end


"Use estimates of the disk structure to take better steps in tau."
function trace_pixel_tilde(xprime, yprime, v0, mol::Molecule, pars::AbstractParameters, interp)
    rmax = 700.0 * AU

    # Calculate DeltaVmax for this ray
    rcyl_min = abs(xprime)
    DeltaVmax = sqrt(2 * kB * temperature(rcyl_min, pars)/mol.mol_weight + (pars.ksi * 1e5)^2) * 1e-5 # km/s


    z0s = geometry.get_zps(xprime, yprime, pars, v0)

    # create a tilde_alpha_func
    # Get the minimum rcyl
    rcyl_min = xprime

    # Occurs at what value of zprime?
    zprime_rcyl_min = cosd(pars.incl)/sind(pars.incl) * yprime

    # Actually, we want the zprime that corresponds to the midplane crossing
    zprime_midplane = -tand(pars.incl) * yprime
    # Calculate the rcyl_mid here

    # Get amplitude at midplane crossing
    # Corresponds to r_cyl_mid, 0.0
    r_cyl_mid = sqrt(xprime^2 + (cosd(pars.incl) * yprime + sind(pars.incl) * zprime_midplane)^2)


    # Get scale height at rcyl_min (will be an underestimate)
    H_rcyl_min = model.Hp(rcyl_min, pars)
    H_rcyl_mid = model.Hp(r_cyl_mid, pars)


    sigma_Upsilon = H_rcyl_mid / cosd(pars.incl)
    a_Upsilon = Upsilon_nu(r_cyl_mid, 0.0, mol.nu_0, pars, mol) * sqrt(2pi) * sigma_Upsilon
    mu_Upsilon = zprime_midplane

    DeltaV = sqrt(2 * kB * model.temperature(rcyl_min, pars)/mol.mol_weight + (pars.ksi * 1e5)^2)

    mu_upsilon1, mu_upsilon2 = z0s
    sigma_upsilon = DeltaV * sqrt(2)/3 * (xprime * sqrt(G * pars.M_star * M_sun) * sind(pars.incl))^(4./3) / (mu_upsilon1 * (v0 * 1.e5)^(7/3))
    a_upsilon = sqrt(2pi) * sigma_upsilon

    sigma_c2 = 1/(sigma_Upsilon^(-2) + sigma_upsilon^(-2))
    mu_c1 = sigma_c2 * (mu_Upsilon/sigma_Upsilon^2 + mu_upsilon1/sigma_upsilon^2)
    mu_c2 = sigma_c2 * (mu_Upsilon/sigma_Upsilon^2 + mu_upsilon2/sigma_upsilon^2)

    a_alpha1 = a_Upsilon * a_upsilon / (sqrt(2pi) * sqrt(sigma_Upsilon^2 + sigma_upsilon^2)) * exp(- (mu_Upsilon - mu_upsilon1)^2/(2 * (sigma_Upsilon^2 + sigma_upsilon^2)))
    a_alpha2 = a_Upsilon * a_upsilon / (sqrt(2pi) * sqrt(sigma_Upsilon^2 + sigma_upsilon^2)) * exp(- (mu_Upsilon - mu_upsilon2)^2/(2 * (sigma_Upsilon^2 + sigma_upsilon^2)))

    function tilde_alpha_func(zprime)
        if zprime > 0
            return a_alpha1 / sqrt(2pi * sigma_c2) * exp(-0.5 * (zprime - mu_c1)^2/sigma_c2)
        else
            return a_alpha2 / sqrt(2pi * sigma_c2) * exp(-0.5 * (zprime - mu_c2)^2/sigma_c2)
        end

    end

    "Given tilde_tau, invert to find zprime"
    function tilde_tau_inv(tilde_tau)
        zprime = sqrt(2 * sigma_c2) * erfcinv(tilde_tau * 2 * sqrt(pi)/a_alpha1) + mu_c1
    end

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

    args = (z1start, v0, pars, interp, gpre, tilde_alpha_func, tilde_tau_inv)

    println("z1start ", z1start/AU)
    # stop when we get to dtau of a large value (will trigger completion before)
    tspan = (0.001, 100.)

    prob = ODEProblem(f_tilde, u0, tspan, args)

    # sol = solve(prob, alg, callback=cb, reltol=rt, abstol=at, dtmax=5*AU, dt=dt)
    sol = solve(prob, alg, dtmax=0.5, dt=0.2, callback=cb)
    # println(sol)


    tau_final, int_final = sol.u[end]
    return tau_final, int_final

    # if (tau_final < tau_thresh) & (length(bounds) == 4)
    #     #We didn't reach the threshold in the first bounding region, so continue in the second.
    #     tspan = (0, z2start - z2end)
    #
    #     # Start the ode to integrate the first path from the end of the second.
    #     u0 = Float64[tau_final; int_final]
    #     args = (z2start, v0, pars, interp, gpre)
    #     prob = ODEProblem(f, u0, tspan, args)
    #     # sol = solve(prob, alg, callback=cb, reltol=rt, abstol=at, dense=false, save_everystep=false, dtmax=3*AU, dt=dt)
    #     sol = solve(prob, alg, reltol=rt, abstol=at, dtmax=3*AU, dt=dt)
    #     tau_final, int_final = sol.u[end]
    #
    #     zps = z2start - sol.t
    #
    #     return zps, sol
    # else
    #     return zps, sol
    # end

end

"v0 is the central velocity of the channel, relative to the disk systemic velocity."
function trace_pixel_debug(xprime, yprime, v0, mol::Molecule, pars::AbstractParameters, interp)

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

    # sol = solve(prob, alg, callback=cb, reltol=rt, abstol=at, dtmax=5*AU, dt=dt)
    sol = solve(prob, alg, reltol=rt, abstol=at, dtmax=3*AU, dt=dt)
    # println(sol)

    zps = z1start - sol.t

    tau_final, int_final = sol.u[end]

    if (tau_final < tau_thresh) & (length(bounds) == 4)
        #We didn't reach the threshold in the first bounding region, so continue in the second.
        tspan = (0, z2start - z2end)

        # Start the ode to integrate the first path from the end of the second.
        u0 = Float64[tau_final; int_final]
        args = (z2start, v0, pars, interp, gpre)
        prob = ODEProblem(f, u0, tspan, args)
        # sol = solve(prob, alg, callback=cb, reltol=rt, abstol=at, dense=false, save_everystep=false, dtmax=3*AU, dt=dt)
        sol = solve(prob, alg, reltol=rt, abstol=at, dtmax=3*AU, dt=dt)
        tau_final, int_final = sol.u[end]

        zps = z2start - sol.t

        return zps, sol
    else
        return zps, sol
    end

end


end # module
