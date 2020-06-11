

"""
Uses Xfoil.jl to create a ::AirfoilPrep.Polar object for the specified airfoil.

Arguments:

* Re : Reynolds number at which a ::AirfoilPrep.Polar object is desired
* xs : vector containing x coordinates of the airfoil contour
* ys : vector containing y coordinates of the airfoil contour

Optional arguments:

* alphas              : angles of attack (radians) at which to evaluate XFOIL for Polar construction
* M                   : desired Mach number
* `iter::Int`         : XFOIL convergence iterations, default to 100.
* `alpha_iter::Int`   : Max AOA convergence iterations, default to 10.
* `verbose::Bool`     : Verbose.

"""
function runXFOIL(Re, xs, ys; alphas=range(-25,stop=25,length=101)*pi/180.0, M=0.0, iter::Int=200, alpha_ite::Int=10, verbose=false)

    # run XFOIL
    cls, cds, cdps, cms, convs = Xfoil.xfoilsweep(xs, ys, alphas.*180/pi, Re; iter=iter, npan=160, mach=0, percussive_maintenance=true, printdata=verbose, zeroinit=true, clmaxstop=true, clminstop=true)
    # remove unconverged points
    alphas_converged, cls_converged, cds_converged, cms_converged = removeSingularities(alphas, cls, cds, cdps, cms, convs)
    # use AirfoilPrep to wrap airfoilpreppy
    polar = Polar(Re, alphas .* 180/pi, cls_converged, cds_converged, cms_converged, xs, ys)

    return polar
end

"""
Returns:

* alphas_converged  : vector of converged coefficients
* cls_converged : vector of converged coefficients
* cds_converged : vector of converged coefficients
* cdps_converged: vector of converged coefficients
* cms_converged : vector of converged coefficients

"""
function removeSingularities(alphas, cls, cds, cdps, cms, convs)
    iconverged = Bool.(convs)
    alphas_converged = alphas[iconverged]
    cls_converged = cls[iconverged]
    cds_converged = cds[iconverged]
    cdps_converged= cdps[iconverged]
    cms_converged = cms[iconverged]

    return alphas_converged, cls_converged, cds_converged, cms_converged
end
