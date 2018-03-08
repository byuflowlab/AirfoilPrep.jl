using AirfoilPrep
using Xfoil
using PyPlot
PyPlot.close("all")
using JLD
using AirfoilManip


# Define operating conditions
aoas = collect(linspace(-5,20,51))#linspace(-10,10,20)
Res = collect(linspace(1000,1000000,2))
Ms = collect(linspace(0,.75,2))
r_over_R = 0.75
c_over_r = 0.3
TSR = 8.0
grid_alphas=[i for i in -180:1.0:180]

#**NOTE 2: Airfoil points x,y must go from trailing edge around the top, then the
#bottom and end back at the trailing edge.**

#Load in airfoil x and y points
folder,_ = splitdir(@__FILE__)
airfoil_file = joinpath(folder,"data","ClarkY_xy_top.txt")
headerlines = 0
open(airfoil_file,"r") do f
    global  xtop = Array{Float64,1}(0)
    global  ytop = Array{Float64,1}(0)
    for (i,line) in enumerate(eachline(f))
        if i>headerlines
            xtop = append!(xtop,parse(split(chomp(line))[1]))
            ytop = append!(ytop,parse(split(chomp(line))[2]))
        else
        end
    end
end

folder,_ = splitdir(@__FILE__)
airfoil_file = joinpath(folder,"data","ClarkY_xy_bot.txt")
headerlines = 0
open(airfoil_file,"r") do f
    global  xbot = Array{Float64,1}(0)
    global  ybot = Array{Float64,1}(0)
    for (i,line) in enumerate(eachline(f))
        if i>headerlines
            xbot = append!(xbot,parse(split(chomp(line))[1]))
            ybot = append!(ybot,parse(split(chomp(line))[2]))
        else
        end
    end
end


reverse!(xtop)
reverse!(ytop)
x = [xtop;xbot]
y = [ytop;ybot]
PyPlot.close("all")
PyPlot.figure()
PyPlot.plot(x,y)


#Wapper function for my analysis code: Xfoil
println("Running Xfoil as Airfoil Data generator")
function f(Re,M)
    # Note that aoas is inherited for more efficient xfoil operation
    cls,cds,cdps,cms,convs =Xfoil.xfoilsweep(x,y,aoas,Re;iter=100,npan=140,mach=M,
    percussive_maintenance=true,printdata=false,zeroinit=true,clmaxstop=true,clminstop=true)

    println("Re:$(round(Int,Re)), M:$(round(M,3))")
    return cls,cds+cdps,cms,convs
end


var_input = (aoas,Res,Ms)
var_names = ["aoa","Re","M"]
response_names = ["cl","cd","cm","convs"]

#Since the version of Xfoil being used is more efficient if it handles the aoa sweep, we'll not generate a table with it yet.
response_values = AirfoilPrep.genNDarray(f,response_names,var_input[2:end],var_names[2:end])

# Reformat to get the ND array in the right format for my specific problem, ie
# a table of responses lining up to aoa, re, mach
cls = zeros(length(aoas),length(Res),length(Ms))
cds = zeros(cls)
cms = zeros(cls)
convs = zeros(cls)

for i = 1:length(aoas)
    for j = 1:length(Res)
        for k = 1:length(Ms)
            cls[i,j,k] = response_values[j,k][1][i]
            cds[i,j,k] = response_values[j,k][2][i]
            cms[i,j,k] = response_values[j,k][3][i]
            convs[i,j,k] = response_values[j,k][4][i]
        end
    end
end

#Put the response values into the format required by NDtools
response_values2 = [cls,cds,cms,convs]
NDtable = AirfoilPrep.TableND(response_values2,response_names,var_input,var_names)

#Save the table
fileLoc,_ = splitdir(@__FILE__)
JLD.save("$(fileLoc)/Data/af_prop_ClarkY.jld", "NDtable", NDtable)

# After saving your airfoil data, you would need to load it and do the following for your in-the-loop analysis

# #------- Example Extrapolation and Splining --------#
fileLoc,_ = splitdir(@__FILE__)
NDtable = JLD.load("$(fileLoc)/Data/af_prop_ClarkY.jld")
NDtable = NDtable["NDtable"]
## Includes convergence checking for airfoil data
NDextrap3D_3Dtable = AirfoilPrep.NDTable_correction3D_extrap(NDtable,r_over_R,c_over_r,TSR;grid_alphas=grid_alphas)



# Spline the table
splout_extrap = AirfoilPrep.SplineND_from_tableND(NDextrap3D_3Dtable)

#Plot cl
extrap_cl =zeros(length(grid_alphas),length(Res))
PyPlot.figure("Verify_cl_Re")
for i = 1:length(Res) #length of the airfoiltools data
    for j = 1:length(grid_alphas)
        vars = [grid_alphas[j],Res[i],0.0]
        extrap_cl[j,i] = AirfoilPrep.interpND(splout_extrap[1],vars)

    end
    PyPlot.plot(grid_alphas,extrap_cl[:,i],label = "Re: $(round(Int,Res[i]))")
end
PyPlot.xlabel("aoa (deg)")
PyPlot.ylabel("cl")
PyPlot.legend(loc = "best")

extrap_cl =zeros(length(grid_alphas),length(Ms))
PyPlot.figure("Verify_cl_M")
for i = 1:length(Ms) #length of the airfoiltools data
    for j = 1:length(grid_alphas)
        vars = [grid_alphas[j],Res[end],Ms[i]]
        extrap_cl[j,i] = AirfoilPrep.interpND(splout_extrap[1],vars)

    end
    PyPlot.plot(grid_alphas,extrap_cl[:,i],label = "M: $(Ms[i])")
end
PyPlot.xlabel("aoa (deg)")
PyPlot.ylabel("cl")
PyPlot.legend(loc = "best")

#Plot cd
extrap_cd =zeros(length(grid_alphas),length(Res))
PyPlot.figure("Verify_cd_Re")
for i = 1:length(Res) #length of the airfoiltools data
    for j = 1:length(grid_alphas)
        vars = [grid_alphas[j],Res[i],0.0]
        extrap_cd[j,i] = AirfoilPrep.interpND(splout_extrap[2],vars)

    end
    PyPlot.plot(grid_alphas,extrap_cd[:,i],label = "Re: $(round(Int,Res[i]))")
end
PyPlot.xlabel("aoa (deg)")
PyPlot.ylabel("cd")
PyPlot.legend(loc = "best")

extrap_cd =zeros(length(grid_alphas),length(Ms))
PyPlot.figure("Verify_cd_M")
for i = 1:length(Ms) #length of the airfoiltools data
    for j = 1:length(grid_alphas)
        vars = [grid_alphas[j],Res[1],Ms[i]]
        extrap_cd[j,i] = AirfoilPrep.interpND(splout_extrap[2],vars)

    end
    PyPlot.plot(grid_alphas,extrap_cd[:,i],label = "M: $(Ms[i])")
end
PyPlot.xlabel("aoa (deg)")
PyPlot.ylabel("cd")
PyPlot.legend(loc = "best")
