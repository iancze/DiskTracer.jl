
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
    println("alpha / tilde_alpha, ", alpha/tilde_alpha)

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
    a_Upsilon = Upsilon_nu(r_cyl_mid, 0.0, mol.nu_0, pars, mol)
    mu_Upsilon = zprime_midplane

    DeltaV = sqrt(2 * kB * model.temperature(rcyl_min, pars)/mol.mol_weight + (pars.ksi * 1e5)^2)

    mu_upsilon1, mu_upsilon2 = z0s
    sigma_upsilon = DeltaV * sqrt(2)/3 * (xprime * sqrt(G * pars.M_star * M_sun) * sind(pars.incl))^(4./3) / (mu_upsilon1 * (v0 * 1.e5)^(7/3))
    # a_upsilon = sqrt(2pi) * sigma_upsilon

    sigma_alpha2 = 1/(sigma_Upsilon^(-2) + sigma_upsilon^(-2)) # actually sigma^2, hence the 2 in the name
    mu_alpha1 = sigma_alpha2 * (mu_Upsilon/sigma_Upsilon^2 + mu_upsilon1/sigma_upsilon^2)
    mu_alpha2 = sigma_alpha2 * (mu_Upsilon/sigma_Upsilon^2 + mu_upsilon2/sigma_upsilon^2)

    a_alpha1 = a_Upsilon * exp(- (mu_Upsilon - mu_upsilon1)^2/(2 * (sigma_Upsilon^2 + sigma_upsilon^2)))
    a_alpha2 = a_Upsilon * exp(- (mu_Upsilon - mu_upsilon2)^2/(2 * (sigma_Upsilon^2 + sigma_upsilon^2)))

    function tilde_alpha_func(zprime)
        if zprime > 0
            return a_alpha1 * exp(-0.5 * (zprime - mu_alpha1)^2/sigma_alpha2)
        else
            return a_alpha2 * exp(-0.5 * (zprime - mu_alpha2)^2/sigma_alpha2)
        end

    end

    "Given tilde_tau, invert to find zprime"
    function tilde_tau_inv(tilde_tau)
        zprime = mu_alpha1 + sqrt(2 * sigma_alpha2) * erfcinv(tilde_tau * sqrt(2)/(a_alpha1 * sqrt(sigma_alpha2 * pi)))
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
    tilde_tau_start = tilde_alpha_func(z1start)
    println("Starting at tilde_tau ", tilde_tau_start)

    # stop when we get to dtau of a large value (will trigger completion before)
    tspan = (tilde_tau_start, 100.)

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
    sol = solve(prob, alg, reltol=rt, abstol=at, dt=dt)
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
