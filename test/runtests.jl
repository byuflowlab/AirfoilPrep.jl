using AirfoilPrep
using Base.Test

using CSV
using PyPlot

modulepath,_ = splitdir(@__FILE__)
modulepath = modulepath*"/"

AP = AirfoilPrep

data_path = "data/"



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

verify_correction3D_2()
test_extrapolation()

@test 1 == 2
