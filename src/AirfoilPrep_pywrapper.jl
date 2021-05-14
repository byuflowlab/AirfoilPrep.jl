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
  init_Re::Float64                        # Reynolds number of this polar
  init_alpha::Array{Float64,1}            # Angles of attack (deg)
  init_cl::Array{Float64,1}               # Lift coefficient at each AOA.
  init_cd::Array{Float64,1}               # Drag coefficient at each AOA.
  init_cm::Array{Float64,1}               # Moment coefficient at each AOA.
  Ma::Float64                             # Mach number
  npanels::Int                            # Number of panels
  ncrit::Float64                          # XFOIL's turbulence transition parameter
  x::Array{Float64,1}                     # x-coordinates of airfoil geometry
  y::Array{Float64,1}                     # y-coordinates of airfoil geometry
  xsepup::Array{Float64,1}                # Separation point over upper surface
  xseplo::Array{Float64,1}                # Separation point over lower surface

  # Internal variables
  pyPolar::PyCall.PyObject

  Polar(init_Re, init_alpha, init_cl, init_cd, init_cm;
            Ma=0.0, npanels=0, ncrit=0,
            x=Float64[], y=Float64[],
            xsepup=Float64[], xseplo=Float64[],
            pyPolar=prepy.Polar(init_Re, init_alpha, init_cl, init_cd, init_cm)
        ) = new(
            init_Re, init_alpha, init_cl, init_cd, init_cm,
            Ma, npanels, ncrit,
            x, y,
            xsepup, xseplo,
            pyPolar)
end

"Returns Re of this Polar"
function get_Re(self::Polar)
    return self.pyPolar.Re
end
"Returns Ma of this Polar"
function get_Ma(self::Polar)
    return self.Ma
end
"Returns number of panels of this Polar"
function get_npanels(self::Polar)
    return self.npanels
end
"Returns ncrit of this Polar"
function get_ncrit(self::Polar)
    return self.ncrit
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
"Returns (alpha, xsepup) points of this Polar"
function get_xsepup(self::Polar)
    return (self.init_alpha, self.xsepup)
end
"Returns (alpha, xseplo) points of this Polar"
function get_xseplo(self::Polar)
    return (self.init_alpha, self.xseplo)
end
"Returns (x, y) points of the airfoil geometry of this Polar"
function get_geometry(self::Polar)
    return (self.x, self.y)
end
"Returns a dummy polar object"
function dummy_polar()
    return Polar(-1, Float64[], Float64[], Float64[], Float64[])
end
"Reads a polar as downloaded from Airfoiltools.com"
function read_polar(file_name::String; path::String="", optargs...)
    header = ["Alpha","Cl","Cd","Cdp","Cm","Top_Xtr","Bot_Xtr"]
    data = DataFrames.DataFrame(CSV.File(joinpath(path,file_name), datarow=12, header=header))
    if data[:,1] == []
        error("Unusable data file. Possibly insufficient data in $file_name.  Data must start on row 12 (XFoil polar output formatting is assumed).")
    end
    polar = Polar(-1, data[:,1], data[:,2], data[:,3], data[:,5]; optargs...)
    return polar
end
"Reads a polar as saved from a Polar object"
function read_polar2(file_name::String; path::String="", optargs...)
    header = ["Alpha","Cl","Cd","Cm"]
    data = DataFrames.DataFrame(CSV.File(joinpath(path,file_name), datarow=2, header=header))
    polar = Polar(-1, data[:,1], data[:,2], data[:,3], data[:,4]; optargs...)
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

    new_polar = _pyPolar2Polar(new_pyPolar; _get_nonpypolar_args(self)...)

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

    new_polar = _pyPolar2Polar(new_pyPolar; _get_nonpypolar_args(self)...)

    return new_polar
end

"Plots this Polar. Give it `geometry=(x,y,1)` for ploting the airfoil geometry,
where x and y are the points of the airfoil (the third number gives the size
of the plot)"
function plot(self::Polar; geometry::Bool=true, label="", style=".-",
                              cdpolar=true, rfl_style="-", fig_id="polar_curves",
                              figsize=[7, 5], rfl_figfactor=2/3,
                              title_str="automatic", polar_optargs=[],
                              legend_optargs=[(:loc, "best")],
                              to_plot=["Cl", "Cd", "Cm"],
                              fig=nothing, axs=nothing)

    # Geometry
    if geometry
        x,y = get_geometry(self)
        plot_airfoil(x, y; label=label, style=rfl_style, figfactor=rfl_figfactor, fig_id=fig_id*"_rfl")
    end

    _fig = fig!=nothing ? fig : figure(fig_id, figsize=figsize.*[length(to_plot), 1])
    _axs = axs!=nothing ? axs : _fig.subplots(1, length(to_plot))

    alpha, cl = get_cl(self)
    _, cd = get_cd(self)
    _, cm = get_cm(self)

    if title_str !== "automatic"
        _fig.suptitle(title_str)
    end

    ttl = title_str=="automatic" ? "$(self.pyPolar.Re)" : title_str

    for (pi, this_plot) in enumerate(to_plot)

        ax = _axs[pi]

        if this_plot=="Cl"
            if title_str=="automatic"; ax.title.set_text("Lift curve at Re=$(self.pyPolar.Re)"); end;
            ax.plot(alpha, cl, style, label=label; polar_optargs...)
            ax.set_xlabel(L"Angle of attack $\alpha (^\circ)$")
            ax.set_ylabel(L"C_l")
            ax.grid(true, color="0.8", linestyle="--")

        elseif this_plot=="Cd"
            if title_str=="automatic"; ax.title.set_text("Drag polar at Re=$(self.pyPolar.Re)"); end;
            ax.plot( cdpolar ? cl : alpha, cd, style, label=label; polar_optargs...)
            if cdpolar
              ax.set_xlabel(L"C_l")
            else
              ax.set_xlabel(L"Angle of attack $\alpha (^\circ)$")
            end
            ax.set_ylabel(L"C_d")
            ax.grid(true, color="0.8", linestyle="--")

        elseif this_plot=="Cm"
            if title_str=="automatic"; ax.title.set_text("Moment curve at Re=$(self.pyPolar.Re)"); end;
            ax.plot(alpha, cm, style, label=label; polar_optargs...)
            ax.set_xlabel(L"Angle of attack $\alpha (^\circ)$")
            ax.set_ylabel(L"C_m")
            ax.grid(true, color="0.8", linestyle="--")
        end

        if pi==length(to_plot) && legend_optargs!=nothing; ax.legend(; legend_optargs...); end;
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
        return [-1, -1, -1]
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

    for (i, curve) in enumerate(curves)
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
    na = length(alpha)

    # New polar data
    out_a, out_cl, out_cd, out_cm = zeros(na), zeros(na), zeros(na), zeros(na)
    # Iterates over data averaging repeated angles
    ave_a, ave_cl, ave_cd, ave_cm, cum_i = 0.0, 0.0, 0.0, 0.0, 0
    count = 0
    for (i,this_a) in enumerate(alpha)

        # Accumulates repeated values
        ave_a += this_a
        ave_cl += cl[i]
        ave_cd += cd[i]
        ave_cm += cm[i]
        cum_i += 1

        # Averages repeated values
        next_a = i!=size(alpha)[1] ? alpha[i+1] : nothing
        if i<start_i || this_a!=next_a
            count += 1
            out_a[count] = ave_a/cum_i
            out_cl[count] = ave_cl/cum_i
            out_cd[count] = ave_cd/cum_i
            out_cm[count] = ave_cm/cum_i
            ave_a, ave_cl, ave_cd, ave_cm, cum_i = 0.0, 0.0, 0.0, 0.0, 0
        end
    end

    # Injective polar
    new_polar = Polar(get_Re(polar), out_a[1:count], out_cl[1:count],
                                    out_cd[1:count], out_cm[1:count];
                                          _get_nonpypolar_args(polar)...)

    return new_polar
end

### INTERNAL FUNCTIONS #########################################################
"Given a Python Polar object, it returns a julia Polar object"
function _pyPolar2Polar(pyPolar::PyCall.PyObject; nonpypolar_args...)

    polar = Polar(pyPolar.Re, pyPolar.alpha, pyPolar.cl, pyPolar.cd,
                  pyPolar.cm; nonpypolar_args...)

    return polar
end

"Return arguments unrelated to airfoilpreppy"
function _get_nonpypolar_args(polar::Polar)

    ma = get_Ma(polar)
    npanels = get_npanels(polar)
    ncrit = get_ncrit(polar)
    x, y = get_geometry(polar)
    xsepup = get_xsepup(polar)[2]
    xseplo = get_xseplo(polar)[2]

    return ((:Ma, ma), (:npanels, npanels), (:ncrit, ncrit), (:x, x), (:y, y),
                                        (:xsepup, xsepup), (:xseplo, xseplo))
end
### END OF airfoilprep.py CLASS ################################################
