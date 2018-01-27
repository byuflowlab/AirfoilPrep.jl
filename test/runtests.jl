using AirfoilPrep
using Base.Test

using CSV
using PyPlot

modulepath,_ = splitdir(@__FILE__)
modulepath = modulepath*"/"

AP = AirfoilPrep

data_path = "data/"


"Verify that both XFOIL and `airfoilprep.py` are running well by comparing
calculated curves on a NACA0006 with data from airfoiltools.com"
function verify_xfoil(; alphas=[i for i in -10:1.0:20], iter=100, verbose=true)
  Re=1*10^5
  airfoil = data_path*"naca0006.dat"

  # Reads data from airfoiltools.com
  data = CSV.read(modulepath*data_path*"xf-naca0006-il-100000.csv";
                      datarow=12,
                      header=split("Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr", ','))

  println("****** Reading geometry file...")
  @time x,y = AP.readcontour(modulepath*airfoil)

  println("****** Running XFOIL...")
  @time polar = AP.runXFOIL(x, y, Re; verbose=verbose, iter=iter,
                                      alphas=alphas)

  # Plots data
  data_polar = AP.Polar(Re, data[:,1], data[:,2], data[:,3],
                        data[:,5], x, y)
  AP.plot(data_polar; geometry=false, label="XFOIL airfoiltools.com",cdpolar=false)

  # Plots calculation
  AP.plot(polar; style="^--r", label="XFOIL jlxlight",cdpolar=false)

  # Determines error of each curve
  tol = 0.15
  errs = AP.compare(data_polar, polar)
  test = [0<=err<tol for err in errs]
  labels = ["Cl", "Cd", "Cm"]
  ignore = ["Cm"]
  println("RESULTS")
  res = true
  for (i,label) in enumerate(labels)
    res = res*(label in ignore ? true : test[i])
    println("\t$label curve test: $(test[i]) (val=$(errs[i]), tol=$(tol))")
  end
  println("OVERALL RESULT: $(res ? "Passed" : "Failed")")

  return res
end

"Verify that 3D corrections (stall delay) to 2D airfoil
properties is done correctly. The stall delay is implemented in `airfoilprep.py`
as explained in Du, Z., & Selig, M. (1998), *A 3-D stall-delay model for
horizontal axis wind turbine performance prediction*, and is verified by
replicating the results Du reports in Fig. 7 over an S809 airfoil.
RESULTS: Failed. It doesn't match what Du reports."
function verify_correction3D(; alphas=[i for i in -10:1.0:30], iter=100,
                                verbose=true)
  Re = 3*10^5 # Reynolds number
  # Re = Int(0.5*10^5) # Reynolds number
  Vinf = 10     # (m/s) Wind speed
  RPM = 72
  D = 10        # (m) Rotor diameter
  TSR = (2*pi*RPM/60)*(D/2) / (Vinf)

  airfoil = data_path*"S809.txt"
  println("\n****** Reading geometry file...")
  @time x,y = AP.readcontour(modulepath*airfoil, header_len=2)
  # Reverts the order of point to go from TE around the bottom to LE over top
  # reverse!(x)
  # reverse!(y)

  println("\n****** Running XFOIL...")
  @time polar = AP.runXFOIL(x, y, Re; verbose=verbose, iter=iter,
                                      alphas=alphas)

  fig1 = figure("correction3D")
  title("S809 airfoil lift coefficient")
  xlabel(L"Angle of attack $\alpha (^\circ)$")
  ylabel(L"$C_l$")
  grid(true, color="0.8", linestyle="--")

  # 2D curve
  data = CSV.read(modulepath*data_path*"S809-2DWindtunnel.csv";
                    header=["angle", "Cl"], datarow=1)
  plot(data[:,1], data[:,2], ".k", label="2D wind tunnel data")
  plot(AP.get_cl(polar)[1], AP.get_cl(polar)[2], "-k", label="2D XFOIL")


  # 3D curves
  for (mrkr, style, c_over_r, r_over_R) in [
                  ("^r","--r",11,80),("Pg","-.g",18,47),("sb",":b",30,30)]
    data = CSV.read(modulepath*data_path*"S809-3D0$c_over_r.csv";
                      header=["angle", "Cl"], datarow=1)
    plot(data[:,1], data[:,2], mrkr,
                    label="3D wind tunnel c/r=0.$c_over_r")

    println("\n****** Running 3D corrections (airfoilprep.py)...")
    @time newpolar = AP.correction3D(polar, r_over_R/100, c_over_r/100, TSR)

    plot(AP.get_cl(newpolar)[1], AP.get_cl(newpolar)[2], style,
                      label="3D-corrected XFOIL c/r=0.$c_over_r")
  end
  legend(loc="best")
end


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
  # airfoil = data_path*"S809.txt"
  airfoil = data_path*"naca0006.dat"
  x,y = AP.readcontour(modulepath*airfoil; header_len=1)


  println("\n****** Running XFOIL...")
  @time polar = AP.runXFOIL(x, y, Re; verbose=verbose, iter=iter,
                                      alphas=alphas)

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

verify_xfoil()
verify_correction3D()
verify_correction3D_2()
test_extrapolation()

@test 1 == 2
