module AirfoilPrep

# ------------ GENERIC MODULES -------------------------------------------------
using Dierckx
using LiftProps
#-------------Sub Routines --------------------------------------------
include("AirfoilPreppy_wrapper.jl")
include("NDtools.jl")
include("xfoil_wrapper.jl")


"""Internal: Break up functions calling airfoilpreppy for clarity,
handles aoa extraction out of the ND array to feed into airfoilpreppy
then puts it back into the correct ND array format for b-splining"""
function afpreppy_wrap3(NDtable,coord,grid_alphas,r_over_R,c_over_r,TSR,CDmax,var_indices)
    # var_indices for airfoil data should be Re, Mach, then whatever
    var_indices = round.(Int,var_indices)
    idx_Re = var_indices[1]
    idx_M = var_indices[2]
    #TODO: match input variable names to correct index instead of assuming
    Re = NDtable.var_input[2][idx_Re] #Re assumed to be second variable
    M = NDtable.var_input[3][idx_M]
    alpha = NDtable.var_input[1] #AOA assumed to be in first variable


    cl = NDtable.response_values[1][:,var_indices...] #Assumed to be first response
    cd = NDtable.response_values[2][:,var_indices...] #Assumed to be second response
    cm = NDtable.response_values[3][:,var_indices...] #Assumed to be 3rd response
    conv = NDtable.response_values[4][:,var_indices...] #Assumed to be 3rd response

    cl2 = []
    cd2 = []
    cm2 = []
    alpha2 = []

    # build good data
    i2 = 1
    for i = 1:length(cl)
        if conv[i]==true
            push!(cl2,cl[i])
            push!(cd2,cd[i])
            push!(cm2,cm[i])
            push!(alpha2,alpha[i])
            i2 += 1
        end
    end

    # Warn if the number of converged points is less than half
    if i2<(length(cl)/2)
        warn("Percent of converged solutions is $(i2/length(cl)*100)%")
    end
    # PyPlot.figure("$Re $M")
    # PyPlot.plot(alpha2,cl2)
    polar = Polar(Re, alpha2, cl2, cd2, cm2, coord[:,1], coord[:,2])
    # 3D corrected Polar
    #especially too high or low aoa for the conditions, possibly pop out data, then spline and resample?
    liftslope,zeroliftangle,aoafit,clfit = fitliftslope(alpha2,cl2)
    aoaclmaxlinear,_ = LiftProps.findclmaxlinear(alpha2,cl2,liftslope,zeroliftangle;tol=0.1,interpolate=true)
    aoaclminlinear,_ = LiftProps.findclmaxlinear(alpha2,-cl2,liftslope,zeroliftangle;tol=0.1,interpolate=true)

    if aoaclmaxlinear<aoaclminlinear
        warn("$aoaclmaxlinear max <$aoaclminlinear min degrees for linear region, correcting with min and max converged aoas (deg)")
        aoaclmaxlinear = maximum(alpha2)
        aoaclminlinear = minimum(alpha2)
    end

    # println("min: $aoaclminlinear")
    # println("max: $aoaclmaxlinear")


    newpolar = correction3D(polar, r_over_R, c_over_r, TSR,
    alpha_linear_min=aoaclminlinear, alpha_linear_max=aoaclmaxlinear, alpha_max_corr=maximum(alpha2))

    # Extrapolated polar
    extrap_polar = APextrapolate(newpolar,CDmax; nalpha = 40)

    cl = extrap_polar.init_cl
    cd = extrap_polar.init_cd
    cm = extrap_polar.init_cm
    alpha_extrap = extrap_polar.init_alpha
    alpha_extrap2 = copy(extrap_polar.init_alpha)

    #Remove elements if duplicated aoa
    j=1
    for i = 2:length(alpha_extrap2)
        if alpha_extrap2[i]==alpha_extrap2[i-j]
            deleteat!(alpha_extrap,i-j)
            deleteat!(cl,i-j)
            deleteat!(cd,i-j)
            deleteat!(cm,i-j)
            j+=1 #update index for deleting in array
        end
    end

    clspl = Dierckx.Spline1D(alpha_extrap,cl,bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
    clgrid_alphas = clspl(grid_alphas)

    cdspl = Dierckx.Spline1D(alpha_extrap,cd,bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
    cdgrid_alphas = cdspl(grid_alphas)

    cmspl = Dierckx.Spline1D(alpha_extrap,cm,bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
    cmgrid_alphas = cmspl(grid_alphas)

    return clgrid_alphas,cdgrid_alphas,cmgrid_alphas
end

"""
    `NDTable_correction3D_extrap(NDtable,r_over_R,c_over_r,TSR;grid_alphas=[i for i in -180:1.0:180])`

runs airfoilprepy for each aoa vs cl,cd,cm,etc contained in a TableND object

        Parameters
        ----------
        NDtable: TableND
            a TableND object
        r_over_R: Float
            blade radial position
        c_over_r: Float
            chord to radius ratio at radial position
        TSR: Float
            tip speed ratio: wind turbine analysis definition
        grid_alphas: Array
            linearly space alphas for response value splines

        Returns
        -------
        tableND : TableND
            a new TableND object

        Notes
        -----
        does include some robust measures to get rid of duplicate points, non-converged airfoil data

"""
function NDTable_correction3D_extrap(NDtable,r_over_R,c_over_r,TSR;grid_alphas=[i for i in -180:1.0:180],CDmax = 1.3)
    #Create array of indices to correctly access each cl curve
    var_indices = []# Array{Array{Int64,length(var_input)},1}(length(var_input))
    for i = 2:length(NDtable.var_input) #skip first dimension, assumed to be AOA
        push!(var_indices,range(1,stop=length(NDtable.var_input[i]),length=length(NDtable.var_input[i])))
    end
    map_in = genMapInput(var_indices) #everything except AOA

    coord = zeros(10,2) #TODO? airfoil not included here

    function afpreppy_wrap2(var_indices...)
        cl,cd,cm = afpreppy_wrap3(NDtable,coord,grid_alphas,r_over_R,c_over_r,TSR,CDmax,var_indices)
        return (cl,cd,cm)
    end

    response_values = map(afpreppy_wrap2,map_in...) #TODO: Replace with pmap for parallelization

    cl = [tup[1] for tup in response_values]
    cd = [tup[2] for tup in response_values]
    cm = [tup[3] for tup in response_values]

    cl2 = zeros(length(grid_alphas),prod(size(NDtable.response_values[1])[2:end]))
    cd2 = similar(cl2) * 0.0
    cm2 = similar(cl2) * 0.0
    response_values_extrap = []

    for i = 1:length(grid_alphas)
        cl2[i,:] = [aoa[i] for aoa in cl]
        cd2[i,:] = [aoa[i] for aoa in cd]
        cm2[i,:] = [aoa[i] for aoa in cm]
    end
    cl2 = reshape(cl2,length(grid_alphas),size(NDtable.response_values[1])[2:end]...)
    cd2 = reshape(cd2,length(grid_alphas),size(NDtable.response_values[1])[2:end]...)
    cm2 = reshape(cm2,length(grid_alphas),size(NDtable.response_values[1])[2:end]...)

    # Reshape to be a float of floats(size(variable_inputs))
    response_values2 = Array{Array{Float64,length(NDtable.var_input)},1}(undef,length(NDtable.response_values))
    response_values2[1] = cl2
    response_values2[2] = cd2
    response_values2[3] = cm2
    response_values2[4] = ones(size(cm2))

    #update the variable input for the extrapolated aoa
    var_input = (grid_alphas,NDtable.var_input[2:end]...)

    return TableND(response_values2,NDtable.response_names,var_input,NDtable.var_names)
end #NDTable_correction3D_extrap

end # module
