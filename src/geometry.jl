module geometry

using ..constants
using ..model

export get_r_cyl, get_z, get_vlos, get_bounding_zps, verify_pixel, verify_pixel_height, GeoPrecalc, precalc_geo, get_coords, get_bounding_scale_heights

struct GeoPrecalc
    xp2::Float64
    a1::Float64
    a2::Float64
    b1::Float64
    b2::Float64
    c1::Float64
    c2::Float64
end

# assumes i is in degrees
function get_r_cyl(xprime::Float64, yprime::Float64, zprime::Float64, i::Float64)
    sqrt(xprime^2 + (cosd(i) * yprime - sind(i) * zprime)^2)
end

# assumes i is in degrees
function get_z(xprime::Float64, yprime::Float64, zprime::Float64, i::Float64)
    sind(i) * yprime + cosd(i) * zprime
end

function get_vlos(xprime::Float64, yprime::Float64, zprime::Float64, pars::AbstractParameters)
    sqrt(G * pars.M_star * M_sun) * sind(pars.incl) * xprime /(xprime^2 + yprime^2 + zprime^2)^(3/4.)
end

"Precalculate the quantities necessary to make repeated calls to `get_r_cyl`, `get_z`, and `get_vlos` efficient."
function precalc_geo(xprime::Float64, yprime::Float64, pars::AbstractParameters)
    # For get_r_cyl
    xp2 = xprime^2
    a1 = cosd(pars.incl) * yprime
    a2 = sind(pars.incl)
    # r_cyl = sqrt(xp2 + (a1 - a2 * zprime)^2)

    # For get_z
    b1 = sind(pars.incl) * yprime
    b2 = cosd(pars.incl)
    # z = b1 + b2 * zprime

    # For get_vlos
    c1 = sqrt(G * pars.M_star * M_sun) * sind(pars.incl) * xprime
    c2 = xprime^2 + yprime^2
    # vlos = c1 /(c2 + zprime^2)^(3/4.)

    return GeoPrecalc(xp2, a1, a2, b1, b2, c1, c2)
end

function get_coords(zprime::Float64, gcalc::GeoPrecalc)
    rcyl = sqrt(gcalc.xp2 + (gcalc.a1 - gcalc.a2 * zprime)^2)
    z = gcalc.b1 + gcalc.b2 * zprime
    vlos = gcalc.c1 / (gcalc.c2 + zprime^2)^(3./4)
    return (rcyl, z, vlos)
end


"Verify whether a given pixel has any any emitting regions. Assumes xprime, yprime, and rmax are in cm. The velocity equivalent to the frequency to be traced is given by v0 [km/s]. DeltaVmax is also in km/s.

Returns true/false"
function verify_pixel(xprime::Float64, yprime::Float64, pars::AbstractParameters, v0::Float64, DeltaVmax::Float64, rmax::Float64)

    vb_min = (v0 - 3 * DeltaVmax) * 1e5 # convert from km/s to cm/s
    vb_max = (v0 + 3 * DeltaVmax) * 1e5 # convert from km/s to cm/s

    rho2 = xprime^2 + yprime^2

    if rho2 > rmax^2
        return false
    end

    overlap = (vb_min < 0) & (vb_max > 0)

    if xprime > 0
        if vb_max < 0
            return false
        elseif overlap
            return true
        elseif rho2 > (xprime * sqrt(G * pars.M_star * M_sun) * sind(pars.incl) / vb_min)^(4/3)
            return false
        else
            return true
        end
    elseif xprime < 0
        if vb_min > 0
            return false
        elseif overlap
            return true
        elseif rho2 > (xprime * sqrt(G * pars.M_star * M_sun) * sind(pars.incl) / vb_max)^(4/3)
            return false
        else
            return true
        end
    else
        return true
    end

end

"Verify whether a given pixel has any any emitting regions, based upon a scale height argument too. Assumes xprime, yprime, and rmax are in cm. The velocity equivalent to the frequency to be traced is given by v0 [km/s]. DeltaVmax is also in km/s.

hlim is the threshold for scale height. If no scale heights at or below this value exist in the ray, don't trace.

Returns true/false"
function verify_pixel_height(xprime::Float64, yprime::Float64, pars::AbstractParameters, v0::Float64, DeltaVmax::Float64, rmax::Float64, hlim::Float64)

    vb_min = (v0 - 3 * DeltaVmax) * 1e5 # convert from km/s to cm/s
    vb_max = (v0 + 3 * DeltaVmax) * 1e5 # convert from km/s to cm/s

    rho2 = xprime^2 + yprime^2

    if rho2 > rmax^2
        return false
    end

    overlap = (vb_min < 0) & (vb_max > 0)

    if xprime > 0
        if vb_max < 0
            return false
        elseif !overlap  && (rho2 > (xprime * sqrt(G * pars.M_star * M_sun) * sind(pars.incl) / vb_min)^(4/3))
            return false
        else
            # Just get the bounds and evaluate
            bounds = get_bounding_zps(xprime, yprime, pars, v0, DeltaVmax, rmax)
            if length(bounds) == 2
                z1start, z1end = bounds
                return verify_h_bound(xprime, yprime, z1start, z1end, pars, hlim)
            else
                z1start, z1end, z2start, z2end = bounds
                return verify_h_bound(xprime, yprime, z1start, z1end, pars, hlim) || verify_h_bound(xprime, yprime, z2start, z2end, pars, hlim)
            end

        end
    elseif xprime < 0
        if vb_min > 0
            return false
        elseif !overlap && (rho2 > (xprime * sqrt(G * pars.M_star * M_sun) * sind(pars.incl) / vb_max)^(4/3))
            return false
        else
            # Just get the bounds and evaluate
            bounds = get_bounding_zps(xprime, yprime, pars, v0, DeltaVmax, rmax)
            if length(bounds) == 2
                z1start, z1end = bounds
                return verify_h_bound(xprime, yprime, z1start, z1end, pars, hlim)
            else
                z1start, z1end, z2start, z2end = bounds
                return verify_h_bound(xprime, yprime, z1start, z1end, pars, hlim) || verify_h_bound(xprime, yprime, z2start, z2end, pars, hlim)
            end
        end
    else
        # Just get the bounds and evaluate
        bounds = get_bounding_zps(xprime, yprime, pars, v0, DeltaVmax, rmax)
        if length(bounds) == 2
            z1start, z1end = bounds
            return verify_h_bound(xprime, yprime, z1start, z1end, pars, hlim)
        else
            z1start, z1end, z2start, z2end = bounds
            return verify_h_bound(xprime, yprime, z1start, z1end, pars, hlim) || verify_h_bound(xprime, yprime, z2start, z2end, pars, hlim)
        end

    end
    println("Reached end?")
end


"Get the starting and ending bounding regions on zprime, based only upon the kinematic/geometrical constraints. Assumes xprime, yprime, and rmax are in cm. The velocity equivalent to the frequency to be traced is given by v0 [km/s]. DeltaVmax is also in km/s."
function get_bounding_zps(xprime::Float64, yprime::Float64, pars::AbstractParameters, v0::Float64, DeltaVmax::Float64, rmax::Float64)

    vb_min = (v0 - 3. * DeltaVmax) * 1.e5 # convert from km/s to cm/s
    vb_max = (v0 + 3. * DeltaVmax) * 1.e5 # convert from km/s to cm/s

    vb_crit = xprime * sqrt(G * pars.M_star * M_sun) * sind(pars.incl)/(xprime^2 + yprime^2)^(3./4)

    # We want to assert that this pixel has already fulfilled the zeroth order check that there will be emission
    @assert (((xprime > 0) & (vb_max > 0)) | ((xprime < 0) & (vb_min < 0))) "Pixel will have no emission."

    # Does the range overlap zero?
    overlap = (vb_min < 0) & (vb_max > 0)
    overlap_crit = (vb_crit > vb_min) & (vb_crit < vb_max)

    if (xprime > 0) & (vb_max > 0)
        if overlap
            z1start = rmax
            z2end = -rmax
        elseif vb_min > 0
            z1start = sqrt((xprime  * sqrt(G  * pars.M_star * M_sun) * sind(pars.incl)/vb_min)^(4/3.) - xprime^2 - yprime^2)
            z2end = -sqrt((xprime  * sqrt(G  * pars.M_star * M_sun) * sind(pars.incl)/vb_min)^(4/3.) - xprime^2 - yprime^2)
        end

        if overlap_crit
            # There exists a vb within the range of vbs which yields zprime = 0, so the two regions merge.
            # z1end = 0.0
            # z2start = 0.0
            return Float64[z1start, z2end]
        else
            z1end = sqrt((xprime  * sqrt(G  * pars.M_star * M_sun) * sind(pars.incl)/vb_max)^(4/3.) - xprime^2 - yprime^2)
            z2start = -sqrt((xprime  * sqrt(G  * pars.M_star * M_sun) * sind(pars.incl)/vb_max)^(4/3.) - xprime^2 - yprime^2)
        end

    elseif (xprime < 0) & (vb_min < 0)
        if overlap
            z1start = rmax
            z2end = -rmax
        elseif vb_max < 0
            z1start = sqrt((xprime  * sqrt(G  * pars.M_star * M_sun) * sind(pars.incl)/vb_max)^(4/3.) - xprime^2 - yprime^2)
            z2end = -sqrt((xprime  * sqrt(G  * pars.M_star * M_sun) * sind(pars.incl)/vb_max)^(4/3.) - xprime^2 - yprime^2)
        end

        if overlap_crit
            # z1end = 0.0
            # z2start = 0.0
            return Float64[z1start, z2end]
        else
            z1end = sqrt((xprime  * sqrt(G  * pars.M_star * M_sun) * sind(pars.incl)/vb_min)^(4/3.) - xprime^2 - yprime^2)
            z2start = -sqrt((xprime  * sqrt(G  * pars.M_star * M_sun) * sind(pars.incl)/vb_min)^(4/3.) - xprime^2 - yprime^2)
        end
    end

    # TODO could be faster with pre-allocated array?
    return Float64[z1start, z1end, z2start, z2end]

end


function get_h(xprime::Float64, yprime::Float64, zprime::Float64, pars::AbstractParameters)
    return (sind(pars.incl) * yprime + cosd(pars.incl) * zprime) / sqrt(xprime^2 + (cosd(pars.incl) * yprime - sind(pars.incl) * zprime)^2)
end


"Given a pair of z start and end values, calculate the value of h at these points and the derivative. Make an assesment about whether the path between zstart and zend contains points that are below the threshold value of hlim.

Return true if region should be traced, false if no."
function verify_h_bound(xprime::Float64, yprime::Float64, zstart::Float64, zend::Float64, pars::AbstractParameters, hlim::Float64)
    # Calculate h
    hstart = get_h(xprime, yprime, zstart, pars)
    hend = get_h(xprime, yprime, zend, pars)

    # If either of these values are below the threshold, we know that we're going to need to trace the ray.
    if (abs(hstart) < hlim) || (abs(hend) < hlim)
        return true
    end

    # If hstart and hend are of opposite sign, then we know a zero-crossing occurred.
    if hstart/hend < 0.0
        return true
    end


    # If neither of these values are below the threshold, then the only way we'll want to trace the ray is if
    # the scale height is decreasing as you move from the edges towards the middle of the region.

    # Note that for some rays near the systemic velocity, this isn't actually a good critereon.

    # Due to the way we're tracing (in the -z direction) assume that zend < zstart, which means that
    # the derivative at zstart must be positive and the derivative at zend must be negative for
    # it to even be possible to have a lower value
    dhstart, dhend = get_h_derivative(xprime, yprime, pars, Float64[zstart, zend])

    # Also, the ray is always guaranteed to cross the midplane somewhere, and if these bounds really do satisfy this
    # derivative criterion, it seems like it must cross the midplane in between zstart and zend.

    return (dhstart == true) && (dhend == false)

end

function get_h_derivative{T}(xprime::Float64, yprime::Float64, pars::AbstractParameters, zprime::T)
    return (cosd(pars.incl) * sind(pars.incl) * (sind(pars.incl) * yprime + cosd(pars.incl) .* zprime) .* (cosd(pars.incl) * yprime - sind(pars.incl) .* zprime) .> 0.0)
end

function get_bounding_scale_heights(xprime::Float64, yprime::Float64, pars::AbstractParameters, h::Float64)
    a = h^2 * sind(pars.incl)^2 - cosd(pars.incl)^2
    b = -2 * yprime * cosd(pars.incl) * sind(pars.incl) * (h^2 + 1)
    c = h^2 * xprime^2 - yprime^2 * (sind(pars.incl)^2 - h^2 * cosd(pars.incl)^2)
    sol1 = (-b + sqrt(b^2 - 4 * a * c))/(2 * a)
    sol2 = (-b - sqrt(b^2 - 4 * a * c))/(2 * a)

    sol = Float64[sol1, sol2]
    # Also, get whether the derivative is positive or negative at each one--don't actually care about the value.
    # true is positive, false is negative
    der = get_h_derivative(xprime, yprime, pars, sol)

    return sol, der
end

end  # module
