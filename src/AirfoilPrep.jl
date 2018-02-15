module AirfoilPrep



# ------------ GENERIC MODULES -------------------------------------------------

#-------------Sub Routines --------------------------------------------
include("AirfoilPreppy_Wrapper.jl")
include("NDtools.jl")

function afpreppy_wrap3(NDtable,coord,grid_alphas,r_over_R,c_over_r,TSR,indices)
    # indices for airfoil data should be Re, Mach, then whatever
    indices = round.(Int,indices)
    idx_Re = indices[1]
    #TODO: match input variable names to correct index instead of assuming
    Re = round(Int,NDtable.var_input[2][idx_Re]) #Re assumed to be second variable
    alpha = NDtable.var_input[1] #AOA assumed to be in first variable


    cl = NDtable.response_values[1][:,indices...] #Assumed to be first response
    cd = NDtable.response_values[2][:,indices...] #Assumed to be second response
    cm = NDtable.response_values[3][:,indices...] #Assumed to be 3rd response
    conv = NDtable.response_values[4][:,indices...] #Assumed to be 3rd response

    cl2 = []
    cd2 = []
    cm2 = []

    # build good data
    for i = 1:length(cl)
        if conv[i]==true
            push!(cl2,cl[i])
            #TODO: TODO: TODO: HERE
        end
    end

    #TODO: Use xfoil finding of min and max linear?
    #TODO: include logic to filter out bad data, or non-converged solutions,

    polar = Polar(Re, alpha, cl, cd, cm, coord[:,1], coord[:,2])
    # 3D corrected Polar
    #especially too high or low aoa for the conditions, possibly pop out data, then spline and resample?
    AoAclmax = 5.0
    AoAclmin = -2.0
    CDmax = 1.3
    newpolar = correction3D(polar, r_over_R, c_over_r, TSR,
    alpha_linear_min=AoAclmin, alpha_linear_max=AoAclmax, alpha_max_corr=5.0)#maximum(alpha))
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

    return clgrid_alphas,cdgrid_alphas,cmgrid_alphas
end

function NDTable_correction3D_extrap(NDtable,r_over_R,c_over_r,TSR;grid_alphas=[i for i in -180:1.0:180])
    #Create array of indices to correctly access each cl curve
    var_indices = []# Array{Array{Int64,length(var_input)},1}(length(var_input))
    for i = 2:length(NDtable.var_input) #skip first dimension, assumed to be AOA
        push!(var_indices,linspace(1,length(NDtable.var_input[i]),length(NDtable.var_input[i])))
    end
    map_in = genMapInput(var_indices) #everything except AOA

    coord = zeros(10,2) #TODO? airfoil not included here

    function afpreppy_wrap2(indices...)
        cl,cd,cm = afpreppy_wrap3(NDtable,coord,grid_alphas,r_over_R,c_over_r,TSR,indices)
        return (cl,cd,cm)
    end

    response_values = map(afpreppy_wrap2,map_in...) #TODO: Replace with pmap for parallelization

    cl = [tup[1] for tup in response_values]
    cd = [tup[2] for tup in response_values]
    cm = [tup[3] for tup in response_values]

    cl2 = zeros(length(grid_alphas),size(NDtable.response_values[1])[2:end]...)
    cd2 = zeros(cl2)
    cm2 = zeros(cl2)
    response_values_extrap = []
    for i = 1:length(grid_alphas)
        cl2[i,:] = [aoa[i] for aoa in cl]
        cd2[i,:] = [aoa[i] for aoa in cd]
        cm2[i,:] = [aoa[i] for aoa in cm]
    end

    # Reshape to be a float of floats(size(variable_inputs))
    response_values2 = Array{Array{Float64,length(NDtable.var_input)},1}(length(NDtable.response_values))
    response_values2[1] = cl2
    response_values2[2] = cd2
    response_values2[3] = cm2
    response_values2[4] = ones(cm2)

    #update the variable input for the extrapolated aoa
    var_input = (grid_alphas,NDtable.var_input[2:end]...)

    return TableND(response_values2,NDtable.response_names,var_input,NDtable.var_names)
end #NDTable_correction3D_extrap

end # module
