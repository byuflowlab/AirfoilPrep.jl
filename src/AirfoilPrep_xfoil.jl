

"""
Uses Xfoil.jl to create a ::AirfoilPrep.Polar object for the specified airfoil.

Arguments:

* Re : Reynolds number at which a ::AirfoilPrep.Polar object is desired
* xs : vector containing x coordinates of the airfoil contour
* ys : vector containing y coordinates of the airfoil contour

Optional arguments:

* alphas              : angles of attack (DEGREES) at which to evaluate XFOIL for Polar construction
* M                   : desired Mach number
* `iter::Int`         : XFOIL convergence iterations, default to 100.
* `alpha_iter::Int`   : Max AOA convergence iterations, default to 10.
* `verbose::Bool`     : Verbose.

**NOTE: Airfoil points must go from trailing edge around the top, then the
bottom and end back at the trailing edge.**
"""
function runXFOIL(xs, ys, Re; alphas=range(-10, stop=25, step=0.5), Mach=0.0,
                                iter::Int=100, verbose=false,
                                npanels=160, ncrit=9, optargs...)

    xseps = zeros(2, length(alphas))

    # run XFOIL
    cls, cds, cdps, cms, convs = Xfoil.xfoilsweep(xs, ys, alphas, Re; iter=iter,
                                                    npan=npanels, mach=Mach,
                                                    ncrit=ncrit,
                                                    xsep=xseps,
                                                    printdata=verbose,
                                                    optargs...)

    # remove unconverged points
    (alphas_converged, cls_converged,
     cds_converged, cms_converged,
        xseps_converged) = removeSingularities(alphas, cls, cds, cdps, cms, xseps, convs)

    # use AirfoilPrep to wrap airfoilpreppy
    polar = Polar(Re, alphas_converged, cls_converged, cds_converged, cms_converged;
                        Ma=Float64(Mach), npanels=npanels, ncrit=Float64(ncrit),
                            x=xs, y=ys, xsepup=xseps[1, :], xseplo=xseps[2, :])

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
function removeSingularities(alphas, cls, cds, cdps, cms, xseps, convs)
    iconverged = Bool.(convs)
    alphas_converged = alphas[iconverged]
    cls_converged = cls[iconverged]
    cds_converged = cds[iconverged]
    cdps_converged= cdps[iconverged]
    cms_converged = cms[iconverged]
    xseps_converged = xseps[:, iconverged]

    return alphas_converged, cls_converged, cds_converged, cms_converged, xseps_converged
end
