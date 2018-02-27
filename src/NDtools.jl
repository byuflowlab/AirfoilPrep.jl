using JLD
using Interpolations

"""
  `TableND(response_name,spl_response,m_vars,b_vars,
  var_input,var_names)`

Defines an arbitrary input arbitrary output ND spline with
the following properties:

  # Arguments
  * response_values::ND Array              : Values of response being splined
  * response_name::String               : Name of the response being splined
  * var_input::Tuple(Arrays)    : Input variable values for reference
  * var_names::String     : Names of input variables

NOTE:
"""
type TableND
  # Initialization variables (USER INPUT)

  response_values
  response_names

  # Keep track of what it was run at
  var_input
  var_names


  TableND(response_values,response_names,var_input,var_names) = new(response_values,
  response_names,var_input,var_names)

end

"""
  `SplineND(response_name,spl_response,m_vars,b_vars,
  var_input,var_names)`

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
type SplineND
  # Initialization variables (USER INPUT)
  response_name
  spl_response

  m_vars
  b_vars

  var_input
  var_names


  SplineND(response_name,spl_response,m_vars,b_vars,
  var_input,var_names) = new(response_name,spl_response,m_vars,
  b_vars,var_input,var_names)

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
    map_in = Array{Array{Float64,length(var_input)},1}(length(var_input))
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
function genNDarray(f,response_names,var_input,var_names)

    # Iterate through each variable to get N matrices for the map function combinations
    map_in = genMapInput(var_input)
    response_values = map(f,map_in...) #TODO: Replace with pmap for parallelization

    return response_values#tableND
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

      spl_response = interpolate(response, BSpline(Cubic(Line())), OnCell())

      #create splines for the reverse mapping of variables
      #y = mx+b or # = m*value+b (put in value, get out it's number)
      m_vars = zeros(length(tableND.var_names))
      b_vars = zeros(m_vars)
      for j = 1:length(tableND.var_names)
          m_vars[j] = (length(tableND.var_input[j])-1)/(maximum(tableND.var_input[j])-minimum(tableND.var_input[j]))
          b_vars[j] = Float64(1-m_vars[j]*minimum(tableND.var_input[j]))
      end

      splND=push!(splND,SplineND(tableND.response_names[i],spl_response,
      m_vars,b_vars,tableND.var_input,tableND.var_names))
  end
  return splND
end

"""helper function for accessing ND spline"""
function interpND(splND,vars)

    var_nums = zeros(length(splND.var_names))
    for i = 1:length(splND.var_names)
        var_nums[i] = splND.m_vars[i]*vars[i]+splND.b_vars[i]
        if vars[i]>maximum(splND.var_input[i]) || vars[i]<minimum(splND.var_input[i])
            warn("Accessing spline in extrapolated area, undefined behavior")
        end
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
