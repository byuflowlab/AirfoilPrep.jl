# ------------ GENERIC MODULES -------------------------------------------------
using Dierckx
using Roots
using PyCall
using LaTeXStrings
using PyPlot
# using QuadGK

# Wrap airfoilprep.py
filepath,_ = splitdir(@__FILE__)
# @pyimport imp # syntax deprecated 

const prepy = PyNULL()

function __init__()
    imp = pyimport("imp")
    (file, filename, data) = imp.find_module("airfoilprep", ["$filepath/../../AirfoilPreppy/src/"])
    copy!(prepy, imp.load_module("airfoilprep", file, filename, data))
end

include("airfoilpreppy_wrapper_misc.jl")



################################################################################
# airfoilprep.py WRAPPER CLASS
################################################################################
"""
  `Polar(Re, alpha, cl, cd, cm, x, y)`

Defines section lift, drag, and pitching moment coefficients as a
function of angle of attack at a particular Reynolds number. This object acts
as an interface for the Python object `Polar` from `airfoilprep.py`, with
the following properties:

  * `self.PyPolar.Re`     : (float) Reynolds number.
  * `self.PyPolar.alpha`  : (array) Angles of attacks.
  * `self.PyPolar.cl`     : (array) Lift coefficient at each AOA.
  * `self.PyPolar.cd`     : (array) Drag coefficient at each AOA.
  * `self.PyPolar.cm`     : (array) Moment coefficient at each AOA.

  # Arguments
  * Re::Int64               : Reynolds number of this polar
  * alpha::Array{Float64,1} : Angles of attack (deg)
  * cl::Array{Float64,1}    : Lift coefficient at each AOA.
  * cd::Array{Float64,1}    : Drag coefficient at each AOA.
  * cm::Array{Float64,1}    : Moment coefficient at each AOA.
  * x::Array{Float64,1}     : x-coordinates of airfoil geometry
  * y::Array{Float64,1}     : y-coordinates of airfoil geometry

NOTE: `alpha` input and output is always in degrees.
**NOTE 2: Airfoil points x,y must go from trailing edge around the top, then the
bottom and end back at the trailing edge.**
"""
mutable struct Polar
  # Initialization variables (USER INPUT)
  init_Re::Float64                          # Reynolds number of this polar
  init_alpha::Array{Float64,1}            # Angles of attack (deg)
  init_cl::Array{Float64,1}               # Lift coefficient at each AOA.
  init_cd::Array{Float64,1}               # Drag coefficient at each AOA.
  init_cm::Array{Float64,1}               # Moment coefficient at each AOA.
  x::Array{Float64,1}                     # x-coordinates of airfoil geometry
  y::Array{Float64,1}                     # y-coordinates of airfoil geometry

  # Internal variables
  pyPolar::PyCall.PyObject

#   Polar(init_Re, init_alpha, init_cl, init_cd, init_cm, x=Float64[], y=Float64[],
#           pyPolar=prepy[:Polar](init_Re, init_alpha, init_cl, init_cd, init_cm)
#         ) = new(
#         init_Re, init_alpha, init_cl, init_cd, init_cm, x, y,
#           pyPolar)
    Polar(init_Re, init_alpha, init_cl, init_cd, init_cm, x=Float64[], y=Float64[],
    pyPolar=prepy.Polar(init_Re, init_alpha, init_cl, init_cd, init_cm)
    ) = new(
    init_Re, init_alpha, init_cl, init_cd, init_cm, x, y,
    pyPolar)
end

"Returns Re of this Polar"
function get_Re(self::Polar)
  return self.pyPolar.Re
end
"Returns (alpha, cl) points of this Polar"
function get_cl(self::Polar)
  return (self.pyPolar.alpha, self.pyPolar.cl)
end
"Returns (alpha, cd) points of this Polar"
function get_cd(self::Polar)
  return (self.pyPolar.alpha, self.pyPolar.cd)
end
"Returns (alpha, cm) points of this Polar"
function get_cm(self::Polar)
  return (self.pyPolar.alpha, self.pyPolar.cm)
end
"Returns (x, y) points of the airfoil geometry of this Polar"
function get_geometry(self::Polar)
  return (self.x, self.y)
end
"Returns a dummy polar object"
function dummy_polar()
  return Polar(-1, Float64[], Float64[], Float64[], Float64[],
    Float64[], Float64[])
end

"""
  `correction3D(self::Polar, r_over_R::Float64, chord_over_r::Float64,
                        tsr::Float64; alpha_max_corr=30, alpha_linear_min=-5,
                        alpha_linear_max=5)`

Applies 3-D corrections for rotating sections from the 2-D data.

  Parameters
  ----------
  r_over_R : float
      local radial position / rotor radius
  chord_over_r : float
      local chord length / local radial location
  tsr : float
      tip-speed ratio
  alpha_max_corr : float, optional (deg)
      maximum angle of attack to apply full correction
  alpha_linear_min : float, optional (deg)
      angle of attack where linear portion of lift curve slope begins
  alpha_linear_max : float, optional (deg)
      angle of attack where linear portion of lift curve slope ends

  Returns
  -------
  polar : Polar
      A new Polar object corrected for 3-D effects

  Notes
  -----
  The Du-Selig method :cite:`Du1998A-3-D-stall-del` is used to correct lift, and
  the Eggers method :cite:`Eggers-Jr2003An-assessment-o` is used to correct drag.
"""
function correction3D(self::Polar, r_over_R::Float64, chord_over_r::Float64,
                      tsr::Float64; alpha_max_corr=30, alpha_linear_min=-5,
                      alpha_linear_max=5)
  new_pyPolar = self.pyPolar.correction3D(r_over_R, chord_over_r, tsr,
                  alpha_linear_min=alpha_linear_min,
                  alpha_linear_max=alpha_linear_max,
                  alpha_max_corr=alpha_max_corr)
  new_polar = _pyPolar2Polar(new_pyPolar, self.x, self.y)
  return new_polar
end

"""
  `APextrapolate(self::Polar, cdmax::Float64; AR=nothing, cdmin=0.001,
                    nalpha=15)`
Extrapolates force coefficients up to +/- 180 degrees using Viterna's method
:cite:`Viterna1982Theoretical-and`.

  Parameters
  ----------
  cdmax : float
      maximum drag coefficient
  AR : float, optional
      aspect ratio = (rotor radius / chord_75% radius)
      if provided, cdmax is computed from AR
  cdmin: float, optional
      minimum drag coefficient.  used to prevent negative values that can
      sometimes occur with this extrapolation method
  nalpha: int, optional
      number of points to add in each segment of Viterna method

  Returns
  -------
  polar : Polar
      a new Polar object

  Notes
  -----
  If the current polar already supplies data beyond 90 degrees then
  this method cannot be used in its current form and will just return itself.

  If AR is provided, then the maximum drag coefficient is estimated as

  cdmax = 1.11 + 0.018*AR

"""
function APextrapolate(self::Polar, cdmax::Float64; AR=nothing, cdmin=0.001,
                      nalpha=15)
  new_pyPolar = self.pyPolar.extrapolate(cdmax, AR=AR, cdmin=cdmin,
                                            nalpha=nalpha)
  new_polar = _pyPolar2Polar(new_pyPolar, self.x, self.y)
  return new_polar
end

"Plots this Polar. Give it geometry=(x,y,1) for ploting the airfoil geometry,
where x and y are the points of the airfoil (the third number gives the size
of the plot)"
function plot(self::Polar; geometry::Bool=true, label="", style=".-",
                              cdpolar=true)

  # Geometry
  if geometry
    x,y = get_geometry(self)
    plot_airfoil(x, y; label=label, style="-", figfactor=1.0)
  end

  fig2 = figure("polar_curves", figsize=(7*3,5*1))

  alpha, cl = get_cl(self)
  _, cd = get_cd(self)
  _, cm = get_cm(self)

  # Lift
  subplot(131)
  title("Lift curve at Re=$(self.pyPolar.Re)")
  PyPlot.plot(alpha, cl, style, label=label)
  xlabel(L"Angle of attack $\alpha (^\circ)$")
  ylabel(L"C_l")
  grid(true, color="0.8", linestyle="--")
  if label!=""; legend(loc="best"); end;

  subplot(132)
  title("Drag polar at Re=$(self.pyPolar.Re)")
  PyPlot.plot( cdpolar ? cl : alpha, cd, style, label=label)
  if cdpolar
    xlabel(L"C_l")
  else
    xlabel(L"Angle of attack $\alpha (^\circ)$")
  end
  ylabel(L"C_d")
  grid(true, color="0.8", linestyle="--")
  if label!=""; legend(loc="best"); end;

  subplot(133)
  title("Moment curve at Re=$(self.pyPolar.Re)")
  PyPlot.plot(alpha, cm, style, label=label)
  xlabel(L"Angle of attack $\alpha (^\circ)$")
  ylabel(L"C_m")
  grid(true, color="0.8", linestyle="--")
  if label!=""; legend(loc="best"); end;

  return fig2
end

"Compares two polars. Returns [err1, err2, err3] with an error between 0 and 1
of Cl, Cd, and Cm curves, respectively."
function compare(polar1::Polar, polar2::Polar; verbose::Bool=false)
  out = []
  curves = [:cl, :cd, :cm]

  # Determines the range of alphas to compare
  min_alpha1 = minimum(polar1.pyPolar.alpha)
  min_alpha2 = minimum(polar2.pyPolar.alpha)
  max_alpha1 = maximum(polar1.pyPolar.alpha)
  max_alpha2 = maximum(polar2.pyPolar.alpha)
  min_alpha = max(min_alpha1, min_alpha2)
  max_alpha = min(max_alpha1, max_alpha2)
  if max_alpha < min_alpha # Case that ranges don't intercept
    warn("Requested comparison of polar that don't coincide. Returning -1.")
    return [-1,-1,-1]
  end

  if verbose; println("min_alpha:$min_alpha\tmax_alpha:$max_alpha"); end;

  # Creates splines of each curve on each polar
  splines = Dict()
  for polar in [polar1, polar2] # each polar
    this_splines = []
    push!(this_splines, Dierckx.Spline1D(polar.pyPolar.alpha, polar.pyPolar.cl))
    push!(this_splines, Dierckx.Spline1D(polar.pyPolar.alpha, polar.pyPolar.cd))
    push!(this_splines, Dierckx.Spline1D(polar.pyPolar.alpha, polar.pyPolar.cm))
    splines[polar] = this_splines
  end

  for (i,curve) in enumerate(curves)
    f1(x) = splines[polar1][i](x)
    f2(x) = splines[polar2][i](x)
    absf1(x) = abs(f1(x))
    absf2(x) = abs(f2(x))
    dif(x) = abs(f2(x) - f1(x))

    # Area of the difference
    num = quadgk(dif, min_alpha, max_alpha)[1]
    # Areas of both
    den = max(quadgk(absf1, min_alpha, max_alpha)[1],
                  quadgk(absf2, min_alpha, max_alpha)[1])
    push!(out, num/den)

    if verbose; println("num=$num\tden=$den"); end;
  end

  return out
end


### INTERNAL FUNCTIONS #########################################################
"Given a Python Polar object, it returns a julia Polar object"
function _pyPolar2Polar(pyPolar::PyCall.PyObject, x, y)
  polar = Polar(pyPolar.Re, pyPolar.alpha, pyPolar.cl, pyPolar.cd,
                pyPolar.cm, x, y)
  return polar
end
### END OF airfoilprep.py CLASS ################################################
