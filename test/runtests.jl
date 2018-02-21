# using AirfoilPrep
include("/Users/kmoore/.julia/v0.6/AirfoilPrep/src/AirfoilPrep.jl")
using Base.Test
using Xfoil

using CSV
using PyPlot

modulepath,_ = splitdir(@__FILE__)
modulepath = modulepath*"/"

AP = AirfoilPrep

data_path = "data/"

#-------- TEST AIRFOILPREPY WRAPPER --------#
"Verify that 3D corrections (stall delay) to 2D airfoil
properties is done correctly. The stall delay is implemented in `airfoilprep.py`
as explained in Du, Z., & Selig, M. (1998), *A 3-D stall-delay model for
horizontal axis wind turbine performance prediction*, and is verified by
replicating the results Du reports in Fig. 5 over the FFA airfoil in 5WPX."
function verify_correction3D_2()
    Re = Int(0.5*10^5) # Reynolds number
    Vinf = 8.8     # (m/s) Wind speed
    RPM = 158
    D = 5.35        # (m) Rotor diameter
    TSR = (2*pi*RPM/60)*(D/2) / (Vinf)

    fig1 = figure("correction3D")
    title("FFA 5WPX airfoil lift coefficient")
    xlabel(L"Angle of attack $\alpha (^\circ)$")
    ylabel(L"$C_l$")
    grid(true, color="0.8", linestyle="--")

    # 2D curve
    data = CSV.read(modulepath*data_path*"FFA-2Dcalc.csv";
    header=["angle", "Cl"], datarow=1)
    plot(data[:,1], data[:,2], "-k", label="2D calculated")

    polar = AP.Polar(Re, data[:,1], data[:,2],
    zeros(data[:,1]), zeros(data[:,1]),
    Float64[], Float64[])


    # 3D curves
    for (mrkr, style, c_over_r, r_over_R) in [
        ("^r","--r",16,55),("Pg","-.g",37,30)]
        data = CSV.read(modulepath*data_path*"FFA-3D0$c_over_r.csv";
        header=["angle", "Cl"], datarow=1)
        plot(data[:,1], data[:,2], mrkr,
        label="3D wind tunnel c/r=0.$c_over_r")

        newpolar = AP.correction3D(polar, r_over_R/100, c_over_r/100, TSR,
        alpha_linear_min=0, alpha_linear_max=7,
        alpha_max_corr=30)

        plot(AP.get_cl(newpolar)[1], AP.get_cl(newpolar)[2], style,
        label="3D-airfoilprep.py corrected c/r=0.$c_over_r")
    end
    legend(loc="best")
end


"Tests (but no verifies, since there is nothing to compare to) that the
extrapolation method from `airfoilprep.py` is being called correctly. This
method takes airfoil curves an extrapolates them all around 360 degrees.
NOTE: CDmax in that function is the drag coeff at alpha=90deg or
CDmax=1.11+0.18AR for a blade of aspect ratio AR<50."
function test_extrapolation(; alphas=[i for i in -10:1.0:20], iter=100,
    verbose=true)
    Re=1*10^5
    CDmax = 1.3

    Vinf = 10     # (m/s) Wind speed
    RPM = 72
    D = 10        # (m) Rotor diameter
    TSR = (2*pi*RPM/60)*(D/2) / (Vinf)
    r_over_R = 30
    c_over_r = 18

    data = CSV.read(modulepath*data_path*"FFA-2Dcalc.csv";
    header=["angle", "Cl"], datarow=1)

    polar = AP.Polar(Re, data[:,1], data[:,2],
    zeros(data[:,1]), zeros(data[:,1]),
    Float64[], Float64[])

    newpolar = AP.correction3D(polar, r_over_R/100, c_over_r/100, TSR,
    alpha_linear_min=0, alpha_linear_max=7,
    alpha_max_corr=30)

    # Extrapolated polar
    # extrap_polar1 = AP.extrapolate(newpolar, CDmax)
    extrap_polar = AP.extrapolate(polar, CDmax)

    # AP.plot(extrap_polar1; cdpolar=false)
    AP.plot(extrap_polar; cdpolar=false)

    return extrap_polar
end


#-------- TEST NDtools --------#

#Use S809 NREL airfoil:
#1) Verify Cl, Cd, Cm
#2) Verify 3d correction
#3) Verify extrapolation

#--- Load in Airfoiltools.com S809 Data ---#
S809_Re2E5 = CSV.read(modulepath*data_path*"xf-s809-nr-200000.csv";
header=["alpha" ,"CL","CD","CDp","CM","Top_Xtr","Bot_Xtr"],delim = ",", datarow=13)

S809_Re5E5 = CSV.read(modulepath*data_path*"xf-s809-nr-500000.csv";
header=["alpha" ,"CL","CD","CDp","CM","Top_Xtr","Bot_Xtr"],delim = ",", datarow=13)

S809_Re1E6 = CSV.read(modulepath*data_path*"xf-s809-nr-1000000.csv";
header=["alpha" ,"CL","CD","CDp","CM","Top_Xtr","Bot_Xtr"],delim = ",", datarow=13)

AirfoilToolsData = (S809_Re2E5,S809_Re5E5,S809_Re1E6)

folder,_ = splitdir(@__FILE__)
airfoil_file = joinpath(folder,"data","S809.txt")
headerlines = 2
open(airfoil_file,"r") do f
    global  x = Array{Float64,1}(0)
    global  y = Array{Float64,1}(0)
    for (i,line) in enumerate(eachline(f))
        if i>headerlines
            x = append!(x,parse(split(chomp(line))[1]))
            y = append!(y,parse(split(chomp(line))[2]))
        else
        end
    end
end

aoas = collect(linspace(-15,25,41))#linspace(-10,10,20)
Res = [2E5,3E5,4E5,5E5,6E5,7E5,8E5,9E5,1E6]
Ms = [0.0,0.01]

#Wapper function for my analysis code: happens to be Xfoil
function f(Re,M)
    cls,cds,cdps,cms,convs =Xfoil.xfoilsweep(x,y,aoas,Re;iter=100,npan=140,mach=M,
    percussive_maintenance=true,printdata=true,zeroinit=true,clmaxstop=true,clminstop=true)

    return cls,cds+cdps,cms,convs
end


var_input = (aoas,Res,Ms)
var_names = ["aoa","Re","M"]
response_names = ["cl","cd","cm","convs"]

#Since the version of Xfoil being used is more efficient if it handles the aoa sweep, we'll not generate a table with it yet.
response_values = AirfoilPrep.genNDarray(f,response_names,var_input[2:end],var_names[2:end];
savefile=false,tablename="tableND")

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

#Access the table example
myindices = (1,1,1)
cl = NDtable.response_values[1][myindices...] #Assumed cl to be first response

# Spline the table
#Warning, the airfoil data here has non-converged points
splout_non_extrap = AirfoilPrep.SplineND_from_tableND(NDtable)
test_cl = AirfoilPrep.interpND(splout_non_extrap[1],(0.0,2E5,0.0))
#Plot the results
Re_airfoiltools = [2E5,5E5,1E6]
XfoilData_cl = zeros(length(aoas),length(Re_airfoiltools))
vars = []
for i = 1:length(Re_airfoiltools) #length of the airfoiltools data #TODO TODO TODO
    for j = 1:length(aoas)
        vars = (aoas[j],Re_airfoiltools[i],0.0)
        XfoilData_cl[j,i] = AirfoilPrep.interpND(splout_non_extrap[1],vars)
    end
end

#Plot the AirfoilTools Data vs the newly generated data
rc("figure", figsize=(4.5, 2.6))
rc("font", size=10.0)
rc("lines", linewidth=1.5)
rc("lines", markersize=3.0)
rc("legend", frameon=false)
rc("axes.spines", right=false, top=false)
rc("figure.subplot", left=0.18, bottom=0.18, top=0.97, right=0.92)
rc("axes", color_cycle=["348ABD", "A60628", "009E73", "7A68A6", "D55E00", "CC79A7"])
color_cycle=["#348ABD", "#A60628", "#009E73", "#7A68A6", "#D55E00", "#CC79A7", "#000486","#700000","#006907","#4C0099"]

PyPlot.close("all")
PyPlot.figure("Verify_NDspline")
for i = 1:length(Re_airfoiltools)
    PyPlot.plot(aoas,XfoilData_cl[:,i],".-",color = color_cycle[i],label = "NDtools Re $(round(Int,Re_airfoiltools[i]))")
    PyPlot.plot(AirfoilToolsData[i][:,1],AirfoilToolsData[i][:,2],"--",color = color_cycle[i],label = "Aifoiltools.com Re $(round(Int,Re_airfoiltools[i]))")
end
PyPlot.xlabel("AOA")
PyPlot.ylabel("cl")
PyPlot.legend(loc = "best")

# Test airfoilpreppy on the ND table
r_over_R = 0.1
c_over_r = 0.3
TSR = 10.0

grid_alphas=[i for i in -180:1.0:180]

coord = (x,y)

NDextrap3D_3Dtable = AirfoilPrep.NDTable_correction3D_extrap(NDtable,r_over_R,c_over_r,TSR;grid_alphas=grid_alphas)
#
splout_extrap = AirfoilPrep.SplineND_from_tableND(NDextrap3D_3Dtable)
splout_non_extrap = AirfoilPrep.SplineND_from_tableND(NDtable)

vars = (0,1e6,0.01)
outputWORKS = AirfoilPrep.interpND(splout_extrap[1],vars)
outputWORKS2 = AirfoilPrep.interpND(splout_non_extrap[1],vars)
#
# verify_correction3D_2()
# test_extrapolation()
#
# @test 1 == 2
