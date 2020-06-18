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
"Reads a polar as downloaded from Airfoiltools.com"
function read_polar(file_name::String; path::String="", x=Float64[],y=Float64[])
  header = ["Alpha","Cl","Cd","Cdp","Cm","Top_Xtr","Bot_Xtr"]
  data = CSV.read(joinpath(path,file_name), datarow=12, header=header)
  polar = Polar(-1, data[:,1], data[:,2], data[:,3], data[:,5], x, y)
  return polar
end
"Reads a polar as saved from a Polar object"
function read_polar2(file_name::String; path::String="", x=Float64[],y=Float64[])
  header = ["Alpha","Cl","Cd","Cm"]
  data = CSV.read(joinpath(path,file_name), datarow=2, header=header)
  polar = Polar(-1, data[:,1], data[:,2], data[:,3], data[:,4], x, y)
  return polar
end
"Saves a polar in Polar format"
function save_polar2(self::Polar, file_name::String; path::String="")
  _file_name = file_name*(occursin(".", file_name) ? "" : ".csv")
  f = open(joinpath(path,_file_name),"w")

  # Header
  header = ["Alpha","Cl","Cd","Cm"]
  for (i,h) in enumerate(header)
    write(f, h)
    write(f, i!=size(header)[1] ? "," : "\n")
  end

  # Data
  alpha, cl = get_cl(self)
  _, cd = get_cd(self)
  _, cm = get_cm(self)
  for (i,a) in enumerate(alpha)
    write(f, "$a,$(cl[i]),$(cd[i]),$(cm[i])\n")
  end

  close(f)
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
                      tsr; alpha_max_corr=30, alpha_linear_min=-5,
                      alpha_linear_max=5)
  new_pyPolar = self.pyPolar.correction3D(r_over_R, chord_over_r, tsr,
                  alpha_linear_min=alpha_linear_min,
                  alpha_linear_max=alpha_linear_max,
                  alpha_max_corr=alpha_max_corr)
  new_polar = _pyPolar2Polar(new_pyPolar, self.x, self.y)
  return new_polar
end

"""
  `extrapolate(self::Polar, cdmax::Float64; AR=nothing, cdmin=0.001,
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
function extrapolate(self::Polar, cdmax::Float64; AR=nothing, cdmin=0.001,
                      nalpha=15)
  new_pyPolar = self.pyPolar.extrapolate(cdmax, AR=AR, cdmin=cdmin,
                                            nalpha=nalpha)
  new_polar = _pyPolar2Polar(new_pyPolar, self.x, self.y)
  return new_polar
end

"Plots this Polar. Give it `geometry=(x,y,1)` for ploting the airfoil geometry,
where x and y are the points of the airfoil (the third number gives the size
of the plot)"
function plot(self::Polar; geometry::Bool=true, label="", style=".-",
                              cdpolar=true, rfl_style="-", fig_id="polar_curves",
                              figsize=[7, 5],
                              title_str="automatic", polar_optargs=[],
                              legend_optargs=[(:loc, "best")],
                              to_plot=["Cl", "Cd", "Cm"])

  # Geometry
  if geometry
    x,y = get_geometry(self)
    plot_airfoil(x, y; label=label, style=rfl_style, figfactor=1.0, fig_id=fig_id*"_rfl")
  end

  fig2 = figure(fig_id, figsize=figsize.*[length(to_plot), 1])

  alpha, cl = get_cl(self)
  _, cd = get_cd(self)
  _, cm = get_cm(self)

  if title_str !== "automatic"
      suptitle(title_str)
  end

  ttl = title_str=="automatic" ? "$(self.pyPolar.Re)" : title_str
  dims = 100 + 10*length(to_plot)

  for (pi, this_plot) in enumerate(to_plot)

      subplot(dims + pi)

      if this_plot=="Cl"
          if title_str=="automatic"; title("Lift curve at Re=$(self.pyPolar.Re)"); end;
          PyPlot.plot(alpha, cl, style, label=label; polar_optargs...)
          xlabel(L"Angle of attack $\alpha (^\circ)$")
          ylabel(L"C_l")
          grid(true, color="0.8", linestyle="--")


      elseif this_plot=="Cd"
          if title_str=="automatic"; title("Drag polar at Re=$(self.pyPolar.Re)"); end;
          PyPlot.plot( cdpolar ? cl : alpha, cd, style, label=label; polar_optargs...)
          if cdpolar
            xlabel(L"C_l")
          else
            xlabel(L"Angle of attack $\alpha (^\circ)$")
          end
          ylabel(L"C_d")
          grid(true, color="0.8", linestyle="--")


      elseif this_plot=="Cm"
          if title_str=="automatic"; title("Moment curve at Re=$(self.pyPolar.Re)"); end;
          PyPlot.plot(alpha, cm, style, label=label; polar_optargs...)
          xlabel(L"Angle of attack $\alpha (^\circ)$")
          ylabel(L"C_m")
          grid(true, color="0.8", linestyle="--")
      end

      if pi==length(to_plot); legend(; legend_optargs...); end;
  end
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
    @warn("Requested comparison of polar that don't coincide. Returning -1.")
    return [-1,-1,-1]
  end

  if verbose; println("min_alpha:$min_alpha\tmax_alpha:$max_alpha"); end;

  # Creates splines of each curve on each polar
  splines = Dict()
  for polar in [polar1, polar2] # each polar
    this_splines = []
    for curve in curves # each curve
      this_spline = Dierckx.Spline1D(polar.pyPolar.alpha, getproperty(polar.pyPolar, curve))
      push!(this_splines, this_spline)
    end
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

"Returns a polar that is injective in alpha (angles are unique)"
function injective(polar::Polar; start_i::Int64=2)
  alpha, cl = get_cl(polar)
  _, cd = get_cd(polar)
  _, cm = get_cm(polar)

  # New polar data
  out_a, out_cl, out_cd, out_cm = Float64[], Float64[], Float64[], Float64[]

  # Iterates over data averaging repeated angles
  ave_val = Dict("a"=>0.0, "cl"=>0.0, "cd"=>0.0, "cm"=>0.0, "cum_i"=>0)
  for (i,this_a) in enumerate(alpha)

    # Accumulates repeated values
    ave_val["a"] += this_a
    ave_val["cl"] += cl[i]
    ave_val["cd"] += cd[i]
    ave_val["cm"] += cm[i]
    ave_val["cum_i"] += 1

    # Averages repeated values
    next_a = i!=size(alpha)[1] ? alpha[i+1] : nothing
    if i<start_i || this_a!=next_a
        push!(out_a, ave_val["a"]/ave_val["cum_i"])
        push!(out_cl, ave_val["cl"]/ave_val["cum_i"])
        push!(out_cd, ave_val["cd"]/ave_val["cum_i"])
        push!(out_cm, ave_val["cm"]/ave_val["cum_i"])
        ave_val = Dict("a"=>0.0, "cl"=>0.0, "cd"=>0.0, "cm"=>0.0, "cum_i"=>0)
    end
  end

  # Injective polar
  x, y = get_geometry(polar)
  new_polar = Polar(get_Re(polar), out_a, out_cl, out_cd, out_cm, x, y)

  return new_polar
end

### INTERNAL FUNCTIONS #########################################################
"Given a Python Polar object, it returns a julia Polar object"
function _pyPolar2Polar(pyPolar::PyCall.PyObject, x, y)
  polar = Polar(pyPolar.Re, pyPolar.alpha, pyPolar.cl, pyPolar.cd,
                pyPolar.cm, x, y)
  return polar
end
### END OF airfoilprep.py CLASS ################################################
