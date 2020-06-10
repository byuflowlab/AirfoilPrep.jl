import Xfoil

"""
Uses Xfoil.jl to create a ::AirfoilPrep.Polar object for the specified airfoil.

Arguments:

* Re : Reynolds number at which a ::AirfoilPrep.Polar object is desired
* xs : vector containing x coordinates of the airfoil contour
* ys : vector containing y coordinates of the airfoil contour

Optional arguments:

* αs                  : angles of attack (radians) at which to evaluate XFOIL for Polar construction
* M                   : desired Mach number
* `iter::Int`         : XFOIL convergence iterations, default to 100.
* `alpha_iter::Int`   : Max AOA convergence iterations, default to 10.
* `verbose::Bool`     : Verbose.

"""
function runXFOIL(Re, xs, ys; αs=range(-25,stop=25,length=101)*pi/180.0, M=0.0, iter::Int=200, alpha_ite::Int=10, verbose=false)

    # run XFOIL
    cls, cds, cdps, cms, convs = Xfoil.xfoilsweep(xs, ys, αs.*180/pi, Re; iter=iter, npan=160, mach=0, percussive_maintenance=true, printdata=verbose, zeroinit=true, clmaxstop=true, clminstop=true)
    # remove unconverged points
    αs_converged, cls_converged, cds_converged, cms_converged = removeSingularities(αs, cls, cds, cdps, cms, convs)
    # use AirfoilPrep to wrap airfoilpreppy
    polar = Polar(Re, αs .* 180/pi, cls_converged, cds_converged, cms_converged, xs, ys)

    return polar
end

"""
Returns:

* αs_converged  : vector of converged coefficients
* cls_converged : vector of converged coefficients
* cds_converged : vector of converged coefficients
* cdps_converged: vector of converged coefficients
* cms_converged : vector of converged coefficients

"""
function removeSingularities(αs, cls, cds, cdps, cms, convs)
    iconverged = Bool.(convs)
    αs_converged = αs[iconverged]
    cls_converged = cls[iconverged]
    cds_converged = cds[iconverged]
    cdps_converged= cdps[iconverged]
    cms_converged = cms[iconverged]

    return αs_converged, cls_converged, cds_converged, cms_converged
end
