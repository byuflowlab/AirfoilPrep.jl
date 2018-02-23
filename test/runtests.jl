using AirfoilPrep
# include("/Users/kmoore/.julia/v0.6/AirfoilPrep/src/AirfoilPrep.jl")
using Base.Test
using Xfoil

using CSV
using PyPlot

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

    # Hard code old solutions for error checking
    old_cl_cr_16 = [0.0299769, 0.0473358, 0.0641102, 0.0826366, 0.0996613, 0.116936, 0.132459, 0.149484, 0.166759, 0.183533, 0.199056, 0.216706, 0.23443, 0.2518, 0.266208, 0.283751, 0.304505, 0.328773, 0.350701, 0.369912, 0.390744, 0.411294, 0.429392, 0.447877, 0.465695, 0.485271, 0.500772, 0.517415, 0.542212, 0.565233, 0.583722, 0.599897, 0.619426, 0.634199, 0.653888, 0.673884, 0.695413, 0.716739, 0.744017, 0.772024, 0.798671, 0.823187, 0.847241, 0.873737, 0.886028, 0.884262, 0.878168, 0.873115, 0.866446, 0.860485, 0.848607, 0.838469, 0.826784, 0.815457, 0.804488, 0.794233, 0.783979, 0.772321, 0.760691, 0.750102, 0.759824, 0.773619, 0.785396]
    old_cl_cr_37 = [0.0337801, 0.0513294, 0.0677512, 0.0849978, 0.101537, 0.118195, 0.134027, 0.150567, 0.167224, 0.183646, 0.199479, 0.21546, 0.23513, 0.251588, 0.262161, 0.281975, 0.30213, 0.326092, 0.346911, 0.366866, 0.386748, 0.406807, 0.425557, 0.443637, 0.461402, 0.484257, 0.501784, 0.514734, 0.538474, 0.56623, 0.585163, 0.603007, 0.623283, 0.64132, 0.662523, 0.683871, 0.713161, 0.74411, 0.802981, 0.866457, 0.922474, 0.969816, 1.01609, 1.08993, 1.17158, 1.21422, 1.23692, 1.26182, 1.2851, 1.30616, 1.31761, 1.32563, 1.33206, 1.33866, 1.34542, 1.35253, 1.35963, 1.36693, 1.37509, 1.39523, 1.45718, 1.54739, 1.62453]
    old_cl = (old_cl_cr_16,old_cl_cr_37)
    correction3D_error_max = 0.0


    # 3D curves
    i = 1
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
        correction3D_error_max = max(NDTable_error_max,maximum(abs.(AP.get_cl(newpolar)[2]-old_cl[i])))
        i+=1
    end
    legend(loc="best")

    return correction3D_error_max
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


# """#-------- TEST NDtools --------#
#
# #Use S809 NREL airfoil:
# #1) Verify Cl, Cd, Cm
# #2) Verify 3d correction
# #3) Verify extrapolation
# """
function ValidateNDtools_from_Xfoil()

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
    aoas21 = []#linspace(-15,25,100)
    aoas2 = []
    XfoilData_cl1 =[] #zeros(length(aoas2),length(Re_airfoiltools))
    XfoilData_cl =[] #zeros(length(aoas2),length(Re_airfoiltools))
    vars = []
    for i = 1:length(Re_airfoiltools) #length of the airfoiltools data
        for j = 1:length(aoas)
            if NDtable.response_values[4][j,i,1]==true
                vars = (aoas[j],Re_airfoiltools[i],0.0)
                push!(XfoilData_cl1,AirfoilPrep.interpND(splout_non_extrap[1],vars))
                push!(aoas21,aoas[j])
            end
        end
        push!(XfoilData_cl,XfoilData_cl1)
        push!(aoas2,aoas21)
        XfoilData_cl1 =[]
        aoas21 = []
    end

    #Plot the AirfoilTools Data vs the newly generated data
    PyPlot.figure("Verify_NDspline")
    NDTable_error_max = 0.0
    OldXfoilData_s809_Re2E5 = [-0.473515, -0.648893, -0.528673, -0.492563, -0.357591, -0.209986, -0.110082, -0.0481071, 0.00891502, 0.113733, 0.17017, 0.235873, 0.327773, 0.41861, 0.527213, 0.633993, 0.746845, 0.87758, 0.977806, 0.954624, 0.926883]
    OldXfoilData_s809_Re5E5 = [-0.746883, -0.650123, -0.584852, -0.551076, -0.524587, -0.512261, -0.439827, -0.318644, -0.202955, -0.0871805, 0.0295255, 0.146385, 0.263179, 0.37841, 0.492835, 0.60978, 0.726073, 0.837511, 0.943375, 0.960604, 0.983048, 1.01456, 1.03991, 1.0668, 1.09374, 1.1144, 1.14403, 1.15379, 1.1609, 1.1431]
    OldXfoilData_s809_Re1E6 = [-0.809019, -0.867028, -0.75136, -0.688495, -0.653602, -0.622135, -0.57855, -0.519902, -0.448765, -0.327558, -0.208416, -0.0891631, 0.0295293, 0.148423, 0.267439, 0.386354, 0.50429, 0.621041, 0.736629, 0.851222, 0.914147, 0.97305, 1.02751, 1.07725, 1.12474, 1.17066, 1.20861, 1.23612, 1.2648, 1.29016, 1.29755, 1.31023]
    XfoilData_cl_old = (OldXfoilData_s809_Re2E5,OldXfoilData_s809_Re5E5,OldXfoilData_s809_Re1E6)
    for i = 1:length(Re_airfoiltools)
        PyPlot.plot(aoas2[i],XfoilData_cl[i],".-",color = color_cycle[i],label = "NDtools Re $(round(Int,Re_airfoiltools[i]))")
        PyPlot.plot(AirfoilToolsData[i][:,1],AirfoilToolsData[i][:,2],"--",color = color_cycle[i],label = "Aifoiltools.com Re $(round(Int,Re_airfoiltools[i]))")
        NDTable_error_max = max(NDTable_error_max,maximum(abs.(XfoilData_cl[i]-XfoilData_cl_old[i])))
    end
    PyPlot.xlabel("AOA")
    PyPlot.ylabel("cl")
    PyPlot.legend(loc = "best")



    return NDTable_error_max,NDtable
end

# NDTable_error_max,NDtable = ValidateNDtools_from_Xfoil()
# Test airfoilpreppy on the ND table
r_over_R = 0.1
c_over_r = 0.3
TSR = 10.0
grid_alphas=[i for i in -180:1.0:180]
# function verifyNDtable_extrap(NDtable)
# Includes convergence checking
# using Gallium
# @enter AirfoilPrep.NDTable_correction3D_extrap(NDtable,r_over_R,c_over_r,TSR;grid_alphas=grid_alphas)
NDextrap3D_3Dtable = AirfoilPrep.NDTable_correction3D_extrap(NDtable,r_over_R,c_over_r,TSR;grid_alphas=grid_alphas)

# Spline the new table
splout_extrap = AirfoilPrep.SplineND_from_tableND(NDextrap3D_3Dtable)


#Plot the output and get the error
Re_airfoiltools = [2E5,5E5,1E6]
aoas = collect(linspace(-15,25,41))#linspace(-10,10,20)
extrap_aoas = grid_alphas
extrap_cl =zeros(length(extrap_aoas),length(Re_airfoiltools))
vars = []
ND_corr3Dextr_maxerror = 0.0
#START:Debug output of splined cl 
for i = 1:length(Re_airfoiltools) #length of the airfoiltools data
    for j = 1:length(aoas)
        vars = (extrap_aoas[j],Re_airfoiltools[i],0.0)
        extrap_cl[j,i] = AirfoilPrep.interpND(splout_extrap[1],vars)

    end
    println(extrap_cl[:,i])
    # ND_corr3Dextr_maxerror = max(ND_corr3Dextr_maxerror,maximum(abs.(extrap_cl[:,i]-extrap_cl)old[i])))
end

vars = (0,1e6,0.01)
outputWORKS = AirfoilPrep.interpND(splout_extrap[1],vars)


#     return NDTable_error_max,NDTable_correction3D_extrap_error_max
# end

ERROR_TOL = 1E-4
@test NDTable_error_max <=ERROR_TOL

correction3D_error_max = verify_correction3D_2()
@test correction3D_error_max <=ERROR_TOL
