module geometry

using ..constants
using ..model

export get_r_cyl, get_z, get_vlos, get_bounding_zps, verify_pixel, GeoPrecalc, precalc_geo, get_coords

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

Returns true/false"
function verify_pixel_height(xprime::Float64, yprime::Float64, pars::AbstractParameters, v0::Float64, DeltaVmax::Float64, rmax::Float64)
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


end  # module
