module model

export size_au
export AbstractParameters, ParametersStandard, ParametersNuker, convert_vector, convert_dict, registered_params
export Upsilon_nu, S_nu, temperature, DeltaV2, get_grids

# The double dot is because we are now inside the model module, and we want to import the
# constants module, which is part of the enclosing DiskJockey package.
using ..constants

# Define an abstract Parameters type, then subset of parameters that can be used to dispatch specifics

# In this file, we need some way of specifying the RT constants necessary
# E.g., frequency (or rel velocity), transition, isotope
# E.g., 0.1 km/s, J=2-1, 12CO


abstract type AbstractParameters end

"Parameters for the standard model."
mutable struct ParametersStandard <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    T_10::Float64 # [K] temperature at 10 AU
    q::Float64 # temperature gradient exponent
    gamma::Float64 # surface density gradient exponent
    Sigma_c::Float64 # [g/cm^2] surface density at characteristic radius
    ksi::Float64 # [cm s^{-1}] microturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end

"Parameters for the NUKER model."
mutable struct ParametersNuker <: AbstractParameters
    M_star::Float64 # [M_sun] stellar mass
    r_c::Float64 # [AU] characteristic radius
    T_10::Float64 # [K] temperature at 10 AU
    q::Float64 # temperature gradient exponent
    gamma::Float64 # surface density gradient within r_c (negative values yield holes)
    alpha::Float64 # sharpness of transition (2 = smooth, 16 = sharp)
    beta::Float64 # gradient power law outside r_c (~7)
    Sigma_c::Float64 # [g/cm^2] surface density at characteristic radius
    ksi::Float64 # [cm s^{-1}] microturbulence
    dpc::Float64 # [pc] distance to system
    incl::Float64 # [degrees] inclination 0 deg = face on, 90 = edge on.
    PA::Float64 # [degrees] position angle (East of North)
    vel::Float64 # [km/s] systemic velocity (positive is redshift/receeding)
    mu_RA::Float64 # [arcsec] central offset in RA
    mu_DEC::Float64 # [arcsec] central offset in DEC
end


"A dictionary of parameter lists for conversion."
registered_params = Dict([("standard", ["M_star", "r_c", "T_10", "q", "gamma", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"]),
("nuker", ["M_star", "r_c", "T_10", "q", "gamma", "alpha", "beta", "Sigma_c", "ksi", "dpc", "incl", "PA", "vel", "mu_RA", "mu_DEC"])])

registered_types = Dict([("standard", ParametersStandard), ("nuker", ParametersNuker)])

"Unroll a vector of parameter values into a parameter type."
function convert_vector(p::Vector{Float64}, model::AbstractString, fix_params::Vector; args...)
    args = Dict{Symbol}{Float64}(args)

    # The goal is to assemble a length-nparam vector that can be unrolled into a parameter type
    # e.g., ParametersStandard(M_star, r_c, T_10, q, gamma, Sigma_c, ksi, dpc, incl, PA, vel, mu_RA, mu_DEC)

    # Select the registerd parameters corresponding to this model
    # These are the names listed in this file (model.jl)
    reg_params = registered_params[model]

    # fix_params is a list of strings from a config file that list the names of parameters to be fixed to the config
    # value. These names correspond to the registered names.

    # fit_params are the ones in reg_params that are not in (∉) fix_params
    fit_params = filter(x->∉(x,fix_params), reg_params)

    nparams = length(reg_params)

    # Make an empty vector of this same length
    par_vec = Array{Float64}(nparams)

    # This requires assigning p to fit_params
    # Find the indexes that correspond to fit_params
    par_indexes = findin(reg_params, fit_params)
    # Stuff p directly into these
    par_vec[par_indexes] = p

    # Then reading fix_params from args
    # First, create an array of fixed values analogous to p
    p_fixed = Float64[args[convert(Symbol, par)] for par in fix_params]
    par_indexes = findin(reg_params, fix_params)
    par_vec[par_indexes] = p_fixed

    # Now that we are sampling for log10M_gas for the verticalEta model, this part gets tricky.

    # Find the location of logSigma_c and make it Sigma_c
    # Even if we are using verticalEta, this will still be in here because it is in reg_params
    # Only though it currently corresponds to log10M_gas instead of log10Sigma_c
    indSigma_c = findin(reg_params, ["Sigma_c"])
    @assert length(indSigma_c) == 1 "Could not find Sigma_c in order to convert from logSigma_c or logM_gas."

    if model == "standard"
      # Convert from log10M_gas to Sigma_c
      M_gas = 10.^par_vec[indSigma_c] * M_sun # [g]

      # Find gamma and r_c
      r_c = par_vec[2] * AU # [cm]
      gamma = par_vec[5]

      Sigma_c = M_gas * (2 - gamma) / (2 * pi * r_c^2)
      par_vec[indSigma_c] = Sigma_c

    elseif model == "nuker"
      indalpha = findin(reg_params, ["alpha"])
      # println("Converting logalpha at index: ", indalpha)
      par_vec[indalpha] = 10.^par_vec[indalpha]
      par_vec[indSigma_c] = 10.^par_vec[indSigma_c]
    else
      par_vec[indSigma_c] = 10.^par_vec[indSigma_c]
    end
    # Then assembling these in the same orignial order as registered_params, into the parameter
    # type corresponding to the model.
    return registered_types[model](par_vec...)

end

"Used to turn a dictionary of parameter values (from config.yaml) directly into a parameter type. Generally used for synthesis and plotting command line scripts."
function convert_dict(p::Dict, model::AbstractString)
    # Select the registerd parameters corresponding to this model
    reg_params = registered_params[model]
    nparams = length(reg_params)

    if model == "standard"

      M_gas = 10.^p["logM_gas"] * M_sun # [g]

      # Find gamma and r_c
      r_c = p["r_c"] * AU # [cm]
      gamma = p["gamma"]

      Sigma_c = M_gas * (2 - gamma) / (2 * pi * r_c^2)
      p["Sigma_c"] = Sigma_c

    elseif model == "nuker"
      # add a new field, which is the conversion of logSigma_c to Sigma_c
      p["alpha"] = 10^p["logalpha"]
      p["Sigma_c"] = 10^p["logSigma_c"]
    else
      # add a new field, which is the conversion of logSigma_c to Sigma_c
      p["Sigma_c"] = 10^p["logSigma_c"]
    end

    # Using this order of parameters, unpack the dictionary p into a vector
    # reg_params reads Sigma_c, not logSigma_c; alpha, not logalpha
    par_vec = Float64[p[par_name] for par_name in reg_params]

    return registered_types[model](par_vec...)

end


# Assume all inputs to these functions are in CGS units and in *cylindrical* coordinates.
# Parametric type T allows passing individual Float64 or Vectors.
# # Alternate functions accept pars passed around, where pars is in M_star, AU, etc...
"The Keplerian velocity assuming non-zero thickness to the disk."
function velocity(r::Float64, z::Float64, M_star::Float64)
    sqrt.(G * M_star / (r^2 + z^2)^(3./2)) * r
end
velocity{T}(r::T, z::T, pars::AbstractParameters) = velocity(r, z, pars.M_star * M_sun)

"Calculate temperature in cylindrical coordinates."
function temperature{T}(r::T, T_10::Float64, q::Float64)
    T_10 * (r ./ (10. * AU)).^(-q)
end
temperature{T}(r::T, pars::AbstractParameters) = temperature(r, pars.T_10, pars.q)

"Scale height, calculate in cylindrical coordinates."
function Hp{T}(r::T, M_star::Float64, T_10::Float64, q::Float64)
    temp = temperature(r, T_10, q)
    sqrt.(kB * temp .* r.^3./(mu_gas * m_H * G * M_star))
end
Hp{T}(r::T,  pars::AbstractParameters) = Hp(r, pars.M_star * M_sun, pars.T_10, pars.q)


"Calculate the gas surface density using cylindrical coordinates."
function Sigma(r::Float64, pars::AbstractParameters)
    r_c = pars.r_c * AU

    gamma = pars.gamma
    Sigma_c = pars.Sigma_c

    S = Sigma_c * (r/r_c)^(-gamma) * exp(-(r/r_c)^(2 - gamma))

    return S
end

"
    Sigma(r::Float64, pars::ParametersNuker)

Calculate the gas surface density using the Nuker profile."
function Sigma(r::Float64, pars::ParametersNuker)
    r_c = pars.r_c * AU

    gamma = pars.gamma
    Sigma_c = pars.Sigma_c

    alpha = pars.alpha
    beta = pars.beta

    S = Sigma_c * (r/r_c)^(-gamma) * (1 + (r/r_c)^alpha)^((gamma - beta)/alpha)

end


# Delivers a gas density in g/cm^3
function rho(r::Float64, z::Float64, pars::AbstractParameters)
    H = Hp(r, pars)
    S = Sigma(r, pars)

    # Calculate the density
    rho = S/(sqrt(2. * pi) * H) * exp(-0.5 * (z/H)^2)

    return rho
end

# Ksi is microturbulent broadining width in units of km/s. Output of this function
# is in cm/s for RADMC (RADMC manual, eqn 7.12)
function microturbulence(ksi::Float64)
    return ksi * 1.e5 # convert from km/s to cm/s
end

microturbulence(pars::AbstractParameters) = microturbulence(pars.ksi)

# Calculate the partition function for the temperature
# uses Mangum and Shirley expansion
# assumes B0 in Hz, depends on molecule
function Z_partition(T::Float64, mol::Molecule)
    nugget = (h * mol.B0) / (kB * T)
    return 1/nugget + 1/3. + 1/15. * nugget + 4 * nugget^2 / 315 + nugget^3 / 315
end

# Calculate the source function
"Calculate the source function. nu in Hz."
function S_nu(r::Real, z::Real, nu::Real, pars::ParametersStandard)
    T = temperature(r, pars)
    return (2 * h * nu^3)/cc^2 / (exp(h * nu / (kB * T)) - 1)
end

# Calculate Upsilon
"Calculate Upsilon, the pre-factor."
function Upsilon_nu(r::Real, z::Real, nu::Real, pars::ParametersStandard, mol::Molecule)

    T = temperature(r, pars)

    # Many of the terms in Upsilon could be pre-calculated for the molecule.
    g_l = 2 * mol.l + 1

    n_l = mol.X_mol * rho(r, z, pars) * g_l * exp(-mol.T_L / T) / Z_partition(T, mol)

    # sigma_0
    sigma_0 = 16 * pi^3 / (h^2 * cc) * (kB * mol.T_L) * (mol.l + 1) / (mol.l * (2 * mol.l + 1)) * mol.mu^2

    # calculate Delta V
    DeltaV = sqrt(2 * kB * T / mol.mol_weight + microturbulence(pars)^2)

    return n_l * sigma_0 * cc / (DeltaV * mol.nu_0 * sqrt(pi)) * (1 - exp(- h * nu / (kB * T)))
end

"Calculate (Delta V)^2. Returns in cm/s."
function DeltaV2(r::Real, z::Real, pars::ParametersStandard, mol::Molecule)
    T = temperature(r, pars)
    return 2 * kB * T/mol.mol_weight + (pars.ksi * 1e5)^2
end

# Functions to compute grid objects for fast interpolation.
# Because we will always be querying for the same rcyl, z values,
# Return the three quantities at once.

function get_grids(pars, mol, nu, nr, rmax, nz, zmax)
    # Fill out 2D array spaced in zmax, rmax
    rs = linspace(0.0, rmax, nr)
    zs = linspace(0.0, zmax, nz)

    DeltaV2_grid = Array{Float64}(nz, nr)
    S_nu_grid = Array{Float64}(nz, nr)
    Upsilon_nu_grid = Array{Float64}(nz, nr)

    # Don't fill out the first r column, since this blows up to infty due to temperature.
    # Simply copy the second value.
    for i=1:nz
        DeltaV2_grid[i,1] = DeltaV2(rs[2], zs[i], pars, mol)
        S_nu_grid[i,1] = S_nu(rs[2], zs[i], nu, pars)
        Upsilon_nu_grid[i,1] = Upsilon_nu(rs[2], zs[i], nu, pars, mol)
        for j=2:nr
            DeltaV2_grid[i,j] = DeltaV2(rs[j], zs[i], pars, mol)
            S_nu_grid[i,j] = S_nu(rs[j], zs[i], nu, pars)
            Upsilon_nu_grid[i,j] = Upsilon_nu(rs[j], zs[i], nu, pars, mol)
        end
    end

    # Closure which works on this grid.
    function interp(rcyl, z)
        # Look up what index rcyl is in rs
        ind_r = searchsortedfirst(rs, rcyl)

        # Look up what index z is in zs (assume -z = z)
        ind_z = searchsortedfirst(zs, abs(z))

        # If either of these errored, return a large value for Delta V, and 0.0 for S_nu and Upsilon
        if (ind_r > nr) | (ind_z > nz)
            return 1.0, 0.0, 0.0
        else
            return DeltaV2_grid[ind_z, ind_r], S_nu_grid[ind_z, ind_r], Upsilon_nu_grid[ind_z, ind_r]
        end

    end

    return interp
end

# Calculate a grid/interpolation object
function get_S_nu_grid(pars, mol, nu, nr, rmax, nz, zmax)

    # Fill out 2D array spaced in zmax, rmax
    rs = linspace(0.0, rmax, nr)
    zs = linspace(0.0, zmax, nz)

    S_nu_grid = Array{Float64}(nz, nr)

    # Don't fill out the first r column, since this blows up to infty due to temperature.
    # Simply copy the second value.
    for i=1:nz
        S_nu_grid[i,1] = S_nu(rs[2], zs[i], nu, pars)
        for j=2:nr
            S_nu_grid[i,j] = S_nu(rs[j], zs[i], nu, pars)
        end
    end

    # Closure which works on this grid.
    function interp(rcyl, z)
        # Look up what index rcyl is in rs
        ind_r = searchsortedfirst(rs, rcyl)

        # Look up what index z is in zs (assume -z = z)
        ind_z = searchsortedfirst(zs, abs(z))

        # If either of these errored, return 0.0
        if (ind_r > nr) | (ind_z > nz)
            return 0.0
        else
            return S_nu_grid[ind_z, ind_r]
        end

    end

    return interp
end

function get_Upsilon_nu_grid(pars, mol, nu, nr, rmax, nz, zmax)

    # Fill out 2D array spaced in zmax, rmax
    rs = linspace(0.0, rmax, nr)
    zs = linspace(0.0, zmax, nz)

    Upsilon_nu_grid = Array{Float64}(nz, nr)

    # Don't fill out the first r column, since this blows up to infty due to temperature.
    # Simply copy the second value.
    for i=1:nz
        Upsilon_nu_grid[i,1] = Upsilon_nu(rs[2], zs[i], nu, pars, mol)
        for j=2:nr
            Upsilon_nu_grid[i,j] = Upsilon_nu(rs[j], zs[i], nu, pars, mol)
        end
    end

    # Closure which works on this grid.
    function interp(rcyl, z)
        # Look up what index rcyl is in rs
        ind_r = searchsortedfirst(rs, rcyl)

        # Look up what index z is in zs (assume -z = z)
        ind_z = searchsortedfirst(zs, abs(z))

        # If either of these errored, return 0.0
        if (ind_r > nr) | (ind_z > nz)
            return 0.0
        else
            return Upsilon_nu_grid[ind_z, ind_r]
        end

    end

    return interp
end

function get_DeltaV2_grid(pars, mol, nr, rmax, nz, zmax)

    # Fill out 2D array spaced in zmax, rmax
    rs = linspace(0.0, rmax, nr)
    zs = linspace(0.0, zmax, nz)

    DeltaV2_grid = Array{Float64}(nz, nr)

    # Don't fill out the first r column, since this blows up to infty due to temperature.
    # Simply copy the second value.
    for i=1:nz
        DeltaV2_grid[i,1] = DeltaV2(rs[2], zs[i], pars, mol)
        for j=2:nr
            DeltaV2_grid[i,j] = DeltaV2(rs[j], zs[i], pars, mol)
        end
    end

    # Closure which works on this grid.
    function interp(rcyl, z)
        # Look up what index rcyl is in rs
        ind_r = searchsortedfirst(rs, rcyl)

        # Look up what index z is in zs (assume -z = z)
        ind_z = searchsortedfirst(zs, abs(z))

        # If either of these errored, return 0.0
        if (ind_r > nr) | (ind_z > nz)
            return 0.0
        else
            return DeltaV2_grid[ind_z, ind_r]
        end

    end

    return interp
end


# Various projection matrices which may or may not be correct.

function P_x(var)
    mat =  Float64[[1, 0, 0]  [0, cos(var), -sin(var)] [0, sin(var), cos(var)]]
    return mat
end

function P_z(var)
    mat =   Float64[[cos(var), -sin(var), 0] [sin(var), cos(var),   0] [0,        0,          1]]
    return mat
end

# function to transform spherical to cartesian
function P_project(theta, phi)
    mat = Float64[ [sin(theta) * cos(phi), sin(theta) * sin(phi), cos(theta)] [cos(theta) * cos(phi), cos(theta) * sin(phi), -sin(theta)] [-sin(phi), cos(phi), 0]]
    return mat
end

function convert_position(pars, r, theta, phi)
    x = r * sin(theta) * cos(phi)
    y = r * sin(theta) * sin(phi)
    z = r * cos(theta)

    sphere = Float64[r, theta, phi]
    cart = Float64[x, y, z]

    Px = P_x(pars.incl_inner * pi/180)
    Pz = P_z(pars.Omega * pi/180)

    cart_prime = Px.' * Pz.' * cart

    x_prime, y_prime, z_prime = cart_prime

    # now convert to spherical coordinates in the primed frame
    r_prime = r
    theta_prime = acos(z_prime / r)
    phi_prime = atan2(y_prime, x_prime)

    if phi_prime < 0
        phi_prime += 2pi
    end

    sphere_prime = Float64[r_prime, theta_prime, phi_prime]

    # Assert that the magnitudes of all of these vectors are the same
    @assert isapprox(norm(cart), r) "Cartesian norm doesn't match "
    @assert isapprox(norm(cart_prime), r) "Rotated cartesion norm doesn't match"

    return (cart, cart_prime, sphere_prime)
end


# # project it to the primed cartesian unit vectors
# vel_cart_prime = Pproj * vel_sphere_prime
#
# # now rotate this from the primed cartesian frame to the RADMC-3D cartesian frame
# vel_cart = Pz * Px * vel_cart_prime
#
# # now project it on to the spherical unit vectors in the RADMC-3D frame
# vel_sphere = Pproj.' * vel_cart



end # module
