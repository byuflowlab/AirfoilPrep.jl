module AirfoilPrep



# ------------ GENERIC MODULES -------------------------------------------------

#-------------Sub Routines --------------------------------------------
include("AirfoilPreppy_Wrapper.jl")
include("NDtools.jl")

function NDTable_correction3D_extrap(NDtable)
    #for each aoa sweep
    map_in = genMapInput(NDtable.var_input[2:end]) #everything except AOA

    grid_alphas=[i for i in -180:1.0:180]

    coord = zeros(10,2) #TODO? airfoil not included here
    function f2 = (Re,M)
        polar = Polar(Re, alpha, cl, cd, cm, coord[:,1], coord[:,2])
        # 3D corrected Polar
        #TODO: Use xfoil finding of min and max linear?
        newpolar = correction3D(polar, r_over_R, c_over_r, TSR,
        alpha_linear_min=AoAclmin, alpha_linear_max=AoAclmax2, alpha_max_corr=maximum(alpha))
        # Extrapolated polar
        extrap_polar = extrapolate(newpolar, CDmax;nalpha = 40)

        cl = extrap_polar.init_cl
        cd = extrap_polar.init_cd
        cm = extrap_polar.init_cm
        alpha_extrap = extrap_polar.init_alpha

        clspl = Dierckx.Spline1D(alpha_extrap,cl,bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
        clgrid_alphas = clspl(grid_alphas)

        cdspl = Dierckx.Spline1D(alpha_extrap,cd,bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
        cdgrid_alphas = cdspl(grid_alphas)

        cmspl = Dierckx.Spline1D(alpha_extrap,cm,bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
        cmgrid_alphas = cmspl(grid_alphas)
    end

    response_values = map(f2,map_in...) #TODO: Replace with pmap for parallelization




end



end # module
