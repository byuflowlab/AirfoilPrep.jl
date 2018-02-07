using JLD
using Interpolations

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

function genNDarray(f,response_names,var_input,var_names;
    savefile=false,tablename="tableND")

    # var_input = (aoa,Re,M,tc)

    # Set up map input based on variable inputs
    map_in = Array{Array{Float64,length(var_input)},1}(length(var_input))

    # Iterate through each variable to get N matrices for the map function combinations
    for i = 1:length(var_input)
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

    response_values = map(f,map_in...)

    tableND = TableND(response_values,response_names,var_input,var_names)

    if savefile
        JLD.save("$tablename.txt",tablename,tableND)
    end

    return tableND
end


function SplineND_from_tableND(tableND)
  #create the 3D splines
  splND = []
  for i = 1:length(tableND.response_names)
      response = [tup[1] for tup in tableND.response_values]
      println(response)

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

function interpND(splND,vars)

    var_nums = zeros(length(splND.var_names))
    for i = 1:length(splND.var_names)
        var_nums[i] = splND.m_vars[i]*vars[i]+splND.b_vars[i]
        if vars[i]>maximum(splND.var_input[i]) || vars[i]<minimum(splND.var_input[i])
            error("Accessing spline in extrapolated area, undefined behavior")
        end
    end
    response = splND.spl_response[var_nums...]
    # g = gradient(splND, var_nums...)
    return response#, g
end

function plotNDtable(tableND)
    for i = 1:length(tableND.response_names)
        PyPlot.figure("$(tableND.response_names[i])_vs_$(tableND.var_names[1]),$(tableND.var_names[j])")
        for j = 1:length(tableND.var_names)
            # plot aoa and response against each other variable
        end
    end


end



#TEST CODE


function xfoil(airfoil,argsxfoil,aoa,Re,M)
    cl = aoa+Re+M
    return cl,cl/100,cl/10
end

aoa = collect(linspace(1,3,3))#linspace(-10,10,20)
Re = collect(linspace(1,3,3))#linspace(1000,1E8,20)
M = collect(linspace(1,3,4))#linspace(.001,1,20)
tc = collect(linspace(1,3,3))#linspace(.001,1,20)

airfoil = "test"
argsxfoil = (true)
function f(aoa,Re,M)
    return xfoil(airfoil,argsxfoil,aoa,Re,M)
end


var_input = (aoa,Re,M)
var_names = ["aoa","Re","M"]
response_names = ["cl","cd","cm"]
tableout = genNDarray(f,response_names,var_input,var_names;
    savefile=false,tablename="tableND")
# response = [tup[1] for tup in tableout.response_values]
splout = SplineND_from_tableND(tableout)

vars = (1,2,1)
outputWORKS = interpND(splout[1],vars)
