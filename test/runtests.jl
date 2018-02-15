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

# airfoil="./data/R1_smoothed.dat"
# # airfoil="$fileloc/data/NACA0006.dat"
# Re_array = round.(Int,(linspace(1e3,3e6,2)))
# M_array = (linspace(.1,.9,2))
# alphas = (linspace(-10,20,60))
# spline_3D = AirfoilSpline3D.gen3Dspl(airfoil,Re_array,M_array;alphas = alphas,
#   TSR=0.7,r_over_R=0.4,c_over_r=0.1, verify=true, returntable=false)

#-------- TEST NDtools --------#

folder,_ = splitdir(@__FILE__)
airfoil_file = joinpath(folder,"data","naca0006.dat")
open(airfoil_file,"r") do f
    global  x = Array{Float64,1}(0)
    global  y = Array{Float64,1}(0)
    for line in eachline(f)
        x = append!(x,parse(split(chomp(line))[1]))
        y = append!(y,parse(split(chomp(line))[2]))
    end
end

aoas = collect(linspace(-2,5,5))#linspace(-10,10,20)
Res = collect(linspace(1e4,1e6,4))#linspace(1000,1E8,20)
Ms = collect(linspace(0.001,0.1,3))#linspace(.001,1,20)

#Wapper function for my analysis code: happens to be Xfoil
function f(Re,M)
    cls,cds,cdps,cms,convs =Xfoil.xfoilsweep(x,y,aoas,Re;iter=100,npan=140,mach=M,
    percussive_maintenance=true,printdata=true,zeroinit=true,clmaxstop=true,clminstop=true)


    nonconv_iter = 0.0
    for i = 1:length(convs)
        if convs[i]==false
            # cls[i] = 0.0
            # cds[i] = 0.0
            # cdps[i] = 0.0
            # cms[i] = 0.0
            nonconv_iter += 1.0
        end
    end

    if nonconv_iter>length(aoas)*2/3
        warn("more than 2/3 of the airfoil data did not converge")
    end

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
indices = (1,1,2)
cl = NDtable.response_values[1][indices...] #Assumed cl to be first response

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

vars = (0,1e5,0.01)
outputWORKS = AirfoilPrep.interpND(splout_extrap[1],vars)
outputWORKS2 = AirfoilPrep.interpND(splout_non_extrap[1],vars)
#
# verify_correction3D_2()
# test_extrapolation()
#
# @test 1 == 2
