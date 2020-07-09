"""
`TableND(response_name, spl_response, m_vars, b_vars, 
            var_input, var_names)`

Defines an arbitrary input arbitrary output ND spline with
the following properties:

# Arguments
* response_values::ND Array              : Values of response being splined
* response_name::String               : Name of the response being splined
* var_input::Tuple(Arrays)    : Input variable values for reference
* var_names::String     : Names of input variables
"""
mutable struct TableND
  # Initialization variables (USER INPUT)

  response_values
  response_names

  # Keep track of what it was run at
  var_input
  var_names


  TableND(response_values, response_names, var_input, var_names) = new(response_values, 
      response_names, var_input, var_names)

end

"""
`SplineND(response_name, spl_response, m_vars, b_vars, 
var_input, var_names)`

Defines an arbitrary input arbitrary output ND spline with
the following properties:

# Arguments
* response_name::String               : Name of the response being splined
* spl_response::ND Spline via interpolation     : ND Spline of response
* m_vars::Tuple    : Spline access linear function slopes for input variables
* b_vars::Tuple    : Spline access linear function offset for input variables
* var_input::Tuple(Arrays)    : Input variable values for reference
* var_names::String     : Names of input variables


NOTE:
"""
mutable struct SplineND
  # Initialization variables (USER INPUT)
  response_name
  spl_response

  m_vars
  b_vars

  var_input
  var_names


  SplineND(response_name, spl_response, m_vars, b_vars, 
      var_input, var_names) = new(response_name, spl_response, m_vars, 
                                    b_vars, var_input, var_names)

end


"""
`genMapInput(var_input)`

Defines the arrays needed for the julia map function to run all of the combinations of input variables

  Parameters
  ----------
  var_input::Tuple(Arrays)    : Input variable values for reference

  Returns
  ----------
  map_in: Array(Arrays)     :Arrays used by the julia map function to run the airfoil/whatever generation function

  NOTE:
  """
function genMapInput(var_input)
    # Set up map input based on variable inputs
    map_in = Array{Array{Float64,length(var_input)},1}(undef,length(var_input))
    # Iterate through each variable to get N matrices for the map function combinations
    for i = 1:length(var_input)
      # Check if var_inputs are linear

      if (length(var_input[i])>2) #Catch single variable inputs before accessing them
        if (!isapprox((var_input[i][1]-var_input[i][2]), (var_input[i][2]-var_input[i][3])))
          error("Variable inputs MUST be linear for gridded interpolation
          $(var_input[i][1]-var_input[i][2]) != $(var_input[i][2]-var_input[i][3])")
        end
      end

      # Create array that specifies dimensions to reshape function
      tmp = ones(Int64,length(var_input))
      tmp[i] = length(var_input[i])
      # println(tmp)

      # Put values into an array of the correct N
      tmp2 = reshape(var_input[i],tmp...)

      # Create array that specifies dimensions to repeat function
      tmp3 = ones(Int64,length(var_input))
      for j = 1:length(var_input)
        if i != j
          tmp3[j] = length(var_input[j])
        end
      end

      # Repeat the values into the full array and add to the map input tuple
      map_in[i] = repeat(tmp2,outer = tmp3)
    end
    return map_in
end

"""
`genNDarray(f,response_names,var_input,var_names)`
generates ND array of responses for M outputs and N inputs with
a call to an external function (like xfoil).  See example for how
to correctly wrap a function

Parameters
----------
response_names: array string
names of the output values
var_input: tuple array(floats)
input values that airfoil will be run at
var_names: array string
names of the input variables

Returns
-------
response_values : Array(Arrays(Tuples))
values from the data function

Notes
-----
input variables MUST be linearly spaced and accessing the spline outside of
the defined variable inputs could return undefined behavior

"""
function genNDarray(f, response_names, var_input, var_names)

    # Iterate through each variable to get N matrices for the map function combinations
    map_in = genMapInput(var_input)
    response_values = map(f, map_in...) #TODO: Replace with pmap for parallelization

    return response_values #tableND
end

"""
`SplineND_from_tableND(tableND)`

defines the

Parameters
----------
var_input::Tuple(Arrays)    : Input variable values for reference

Returns
----------
map_in: Array(Arrays)     :Arrays used by the julia map function to run the airfoil/whatever generation function

NOTE:
"""
function SplineND_from_tableND(tableND)

    #create the 3D splines
    splND = []
    for i = 1:length(tableND.response_names)
        response = tableND.response_values[i]

        spl_response = interpolate(response, BSpline(Cubic(Line(OnCell()))))

        #create splines for the reverse mapping of variables
        #y = mx+b or # = m*value+b (put in value, get out it's number)
        m_vars = zeros(length(tableND.var_names))
        b_vars = similar(m_vars)
        for j = 1:length(tableND.var_names)
            m_vars[j] = (length(tableND.var_input[j]) - 1) / (maximum(tableND.var_input[j]) - minimum(tableND.var_input[j]))
            b_vars[j] = Float64(1-m_vars[j] * minimum(tableND.var_input[j]))
        end

        splND = push!(splND, SplineND(tableND.response_names[i], spl_response,
                                    m_vars, b_vars, tableND.var_input, tableND.var_names))
    end

    return splND
end

"""helper function for accessing ND spline"""
function interpND(splND, vars)

    var_nums = zeros(length(splND.var_names))

    for i = 1:length(vars)

        if vars[i] > maximum(splND.var_input[i])
            # warn("Accessing spline in extrapolated area, capping at max. Variable $i value: $(vars[i])")
            vars[i] = maximum(splND.var_input[i])
        elseif vars[i]<minimum(splND.var_input[i])
            # warn("Accessing spline in extrapolated area, capping at min. Variable $i value: $(vars[i])")
            vars[i] = minimum(splND.var_input[i])
        end
        #Map to the spline index
        var_nums[i] = splND.m_vars[i] * vars[i] + splND.b_vars[i]
    end
    response = splND.spl_response[var_nums...]
    # g = gradient(splND, var_nums...)

    return response#, g
end

"""Not ready yet"""
function plotNDtable(tableND)

    for i = 1:length(tableND.response_names)
        PyPlot.figure("$(tableND.response_names[i])_vs_$(tableND.var_names[1]),$(tableND.var_names[j])")

        for j = 1:length(tableND.var_names)
            # plot aoa and response against each other variable
        end
    end

end



"""Internal: Break up functions calling airfoilpreppy for clarity,
handles aoa extraction out of the ND array to feed into airfoilpreppy
then puts it back into the correct ND array format for b-splining"""
function afpreppy_wrap3(NDtable, coord, grid_alphas, r_over_R, c_over_r, TSR, CDmax, var_indices)
    # var_indices for airfoil data should be Re, Mach, then whatever
    var_indices = round.(Int, var_indices)
    idx_Re = var_indices[1]
    idx_M = var_indices[2]
    #TODO: match input variable names to correct index instead of assuming
    Re = NDtable.var_input[2][idx_Re] #Re assumed to be second variable
    M = NDtable.var_input[3][idx_M]
    alpha = NDtable.var_input[1] #AOA assumed to be in first variable


    cl = NDtable.response_values[1][:, var_indices...] #Assumed to be first response
    cd = NDtable.response_values[2][:, var_indices...] #Assumed to be second response
    cm = NDtable.response_values[3][:, var_indices...] #Assumed to be 3rd response
    conv = NDtable.response_values[4][:, var_indices...] #Assumed to be 3rd response

    cl2 = []
    cd2 = []
    cm2 = []
    alpha2 = []

    # build good data
    i2 = 1
    for i = 1:length(cl)
        if conv[i]==true
            push!(cl2, cl[i])
            push!(cd2, cd[i])
            push!(cm2, cm[i])
            push!(alpha2, alpha[i])
            i2 += 1
        end
    end

    # Warn if the number of converged points is less than half
    if i2 < (length(cl) / 2)
        @warn("Percent of converged solutions is $(i2 / length(cl) * 100)%")
    end
    # PyPlot.figure("$Re $M")
    # PyPlot.plot(alpha2, cl2)
    polar = Polar(Re, alpha2, cl2, cd2, cm2, x=coord[:,1], y=coord[:,2])
    # 3D corrected Polar
    #especially too high or low aoa for the conditions, possibly pop out data, then spline and resample?
    liftslope, zeroliftangle, aoafit, clfit = fitliftslope(alpha2, cl2)
    aoaclmaxlinear, _ = findclmaxlinear(alpha2, cl2, liftslope, zeroliftangle; tol=0.1, interpolate=true)
    aoaclminlinear, _ = findclmaxlinear(alpha2, -cl2, liftslope, zeroliftangle; tol=0.1, interpolate=true)

    if aoaclmaxlinear<aoaclminlinear
        @warn("$aoaclmaxlinear max <$aoaclminlinear min degrees for linear region, correcting with min and max converged aoas (deg)")
        aoaclmaxlinear = maximum(alpha2)
        aoaclminlinear = minimum(alpha2)
    end

    # println("min: $aoaclminlinear")
    # println("max: $aoaclmaxlinear")


    newpolar = correction3D(polar, r_over_R, c_over_r, TSR,
                                alpha_linear_min=aoaclminlinear, 
                                alpha_linear_max=aoaclmaxlinear, 
                                alpha_max_corr=maximum(alpha2))

    # Extrapolated polar
    extrap_polar = extrapolate(newpolar, CDmax; nalpha = 40)

    cl = extrap_polar.init_cl
    cd = extrap_polar.init_cd
    cm = extrap_polar.init_cm
    alpha_extrap = extrap_polar.init_alpha
    alpha_extrap2 = copy(extrap_polar.init_alpha)

    #Remove elements if duplicated aoa
    j = 1
    for i = 2:length(alpha_extrap2)
        if alpha_extrap2[i] == alpha_extrap2[i-j]
            deleteat!(alpha_extrap, i - j)
            deleteat!(cl, i - j)
            deleteat!(cd, i - j)
            deleteat!(cm, i - j)
            j += 1 #update index for deleting in array
        end
    end

    clspl = Dierckx.Spline1D(alpha_extrap, cl, bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
    clgrid_alphas = clspl(grid_alphas)

    cdspl = Dierckx.Spline1D(alpha_extrap, cd, bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
    cdgrid_alphas = cdspl(grid_alphas)

    cmspl = Dierckx.Spline1D(alpha_extrap, cm, bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
    cmgrid_alphas = cmspl(grid_alphas)

    return clgrid_alphas, cdgrid_alphas, cmgrid_alphas
end

"""
    `NDTable_correction3D_extrap(NDtable, r_over_R, c_over_r, TSR;grid_alphas=[i for i in -180:1.0:180])`

runs airfoilprepy for each aoa vs cl, cd, cm, etc contained in a TableND object

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
function NDTable_correction3D_extrap(NDtable, r_over_R, c_over_r, TSR;grid_alphas=[i for i in -180:1.0:180], CDmax = 1.3)
    #Create array of indices to correctly access each cl curve
    var_indices = [] # Array{Array{Int64, length(var_input)}, 1}(length(var_input))
    for i = 2:length(NDtable.var_input) #skip first dimension,  assumed to be AOA
        push!(var_indices, range(1, stop=length(NDtable.var_input[i]), length=length(NDtable.var_input[i])))
    end
    map_in = genMapInput(var_indices) #everything except AOA

    coord = zeros(10, 2) #TODO? airfoil not included here

    function afpreppy_wrap2(var_indices...)
        cl, cd, cm = afpreppy_wrap3(NDtable, coord, grid_alphas, r_over_R, c_over_r, TSR, CDmax, var_indices)
        return (cl, cd, cm)
    end

    response_values = map(afpreppy_wrap2, map_in...) #TODO: Replace with pmap for parallelization

    cl = [tup[1] for tup in response_values]
    cd = [tup[2] for tup in response_values]
    cm = [tup[3] for tup in response_values]

    cl2 = zeros(length(grid_alphas), prod(size(NDtable.response_values[1])[2:end]))
    cd2 = similar(cl2) * 0.0
    cm2 = similar(cl2) * 0.0
    response_values_extrap = []

    for i = 1:length(grid_alphas)
        cl2[i, :] = [aoa[i] for aoa in cl]
        cd2[i, :] = [aoa[i] for aoa in cd]
        cm2[i, :] = [aoa[i] for aoa in cm]
    end
    cl2 = reshape(cl2, length(grid_alphas), size(NDtable.response_values[1])[2:end]...)
    cd2 = reshape(cd2, length(grid_alphas), size(NDtable.response_values[1])[2:end]...)
    cm2 = reshape(cm2, length(grid_alphas), size(NDtable.response_values[1])[2:end]...)

    # Reshape to be a float of floats(size(variable_inputs))
    response_values2 = Array{Array{Float64, length(NDtable.var_input)}, 1}(undef, length(NDtable.response_values))
    response_values2[1] = cl2
    response_values2[2] = cd2
    response_values2[3] = cm2
    response_values2[4] = ones(size(cm2))

    #update the variable input for the extrapolated aoa
    var_input = (grid_alphas, NDtable.var_input[2:end]...)

    return TableND(response_values2, NDtable.response_names, var_input, NDtable.var_names)
end #NDTable_correction3D_extrap
