################################################################################
# MISCELLANOUS
################################################################################
"Given the output of `airfoil_props()` it generates a file in aerodyn format"
function generate_aerodyn_file(file_name, re, AoA_d,clmin,clmax,AoA_zl_d,
                        lift_slope,cldata,cdfdata,cdpdata,cmdata; numPoints=100)
  ncrit = "-"

  # Generates splines
  Top_Xtr = zeros(AoA_d) # Dummy arrays
  Bot_Xtr = zeros(AoA_d)
  (Cl_spl, Cd_spl, Cdp_spl, Cm_spl, Top_Xtr_spl, Bot_Xtr_spl, aoaZeroCl,
  aoaMinDrag, cdmin, aoamin, aoamax) = _makeSpline(AoA_d, cldata, cdfdata, cdpdata,
  cmdata, Top_Xtr, Bot_Xtr)

  # Creates file
  _file_name = file_name[end-3:end]==".dat" ? file_name : file_name*".dat"

  write(_file_name,"    See filename for airfoil ID
  Cl and Cd values uncorrected
  line
  1        Number of airfoil tables in this file
  $(re/1e6)     Reynolds numbers in millions
  0.0      Control setting
  $(ncrit)     Stall angle (deg)
  $aoaZeroCl   Zero Cn angle of attack (deg)
  0.0   Cn slope for zero lift (dimensionless)
  0.0   Cn extrapolated to value at positive stall angle of attack
  0.0   Cn at stall value for negative angle of attack
  $aoaMinDrag     Angle of attack for minimum CD (deg)
  $cdmin   Minimum CD value\n")

  f = open(_file_name,"a")
  open(_file_name,"a") do x
      for aoa in range(aoamin, aoamax, length=numPoints)
          write(x,"$(aoa)\t $(Cl_spl(aoa))\t $(Cd_spl(aoa))\t $(Cm_spl(aoa)) \n")
      end
      write(x,"EOT")
  end

end

"Plots the contour of an airfoil given in x,y"
function plot_airfoil(x::Array{Float64,1}, y::Array{Float64,1};
                      label="", style="-k", figfactor=1.0,
                      title_str="Airfoil geometry", fig_id="airfoil_geometry")
  # Sizes the figure
  figsize = [7*1.5,5*0.5]*figfactor
  xmin, xmax = -0.05, 1.05
  yrange = (xmax-xmin)/figsize[1] * figsize[2]
  ymin, ymax = -yrange/2, yrange/2
  fig1 = figure(fig_id, figsize=(figsize[1], figsize[2]))
  xlim([xmin, xmax])
  ylim([ymin, ymax])

  PyPlot.plot(x,y, style, label=label)
  xlabel("x")
  ylabel("y")
  grid(true, color="0.8", linestyle="--")
  title(title_str)

  # if label!=""; legend(loc="best"); end;
end

"Receives a .dat file as pulled from airfoiltools.com containing the x and y
contour coordinates of an airfoil, and returns arrays x and y."
function readcontour(file_name; header_len=1, delim=" ", path="", output="arrays")

    x, y = Float64[], Float64[]

    open(joinpath(path,file_name)) do f
        for (i,line) in enumerate(eachline(f))

            # Ignores header
            if i<=header_len
                nothing
                # Parses each line
            else
                this_x, this_y = split(line, delim; keepempty=false)
                push!(x, Meta.parse(Float64, this_x))
                push!(y, Meta.parse(Float64, this_y))
            end

        end
    end

    if output=="arrays"
        return x,y
    elseif output=="matrix"
        return hcat(x,y)
    else
        error("Invalid `output` argument $(output).")
    end
end

"""
  Receives an airfoil contour and splits it up in upper and lower surfaces as
  divided by the chord line. It returns `(upper, lower)` with `upper=(x,y)` the
  points of the upper surface, ditto for `lower`. Both `upper` and `lower` are
  given in increasing order in x (i.e., from leading to trailing edge).
"""
function splitcontour(x,y)
  # ERROR CASES
  if !(x[1] in [0.0, 1.0])
    error("Invalid contour. x[1] must be either 0.0 or 1.0, got $(x[1]).")
  end

  # Flag indicating whether the contour start at the trailing or leading edge
  start_TE = x[1]==1.0

  # Find the opposite end of the contour
  end_i = -1
  for (i, xi) in enumerate(x)
    if i==1
      nothing
    # Case of starting from the trailing edge
    elseif start_TE && xi > x[i-1]
      end_i = i-1
      break
    # Case of leading edge
    elseif !start_TE  && xi < x[i-1]
      end_i = i-1
      break
    end
  end

  # ERROR CASE
  if end_i==-1
    error("Logic error! End of airfoil not found!")
  end

  # Splits them up
  x_sec1, y_sec1 = x[1:end_i], y[1:end_i]
  x_sec2, y_sec2 = x[end_i:end], y[end_i:end]

  # Sorts them from LE to TE
  if x_sec1[1] > 0.5; reverse!(x_sec1); reverse!(y_sec1); end;
  if x_sec2[1] > 0.5; reverse!(x_sec2); reverse!(y_sec2); end;

  # Determines upper and lower surfaces
  if mean(y_sec1) > mean(y_sec2)
    upper = [x_sec1, y_sec1]
    lower = [x_sec2, y_sec2]
  else
    upper = [x_sec2, y_sec2]
    lower = [x_sec1, y_sec1]
  end

  return upper, lower
end

### INTERNAL FUNCTIONS #########################################################
"Receives the output of `airfoil_props()` and returns splines of them"
function _makeSpline(Alpha,Cl,Cd,Cdp,Cm,Top_Xtr,Bot_Xtr)
  # Splines
  Cl_spl = Dierckx.Spline1D(Alpha,Cl)
  Cd_spl = Dierckx.Spline1D(Alpha,Cd+Cdp)
  Cdp_spl = Dierckx.Spline1D(Alpha,Cdp)
  Cm_spl = Dierckx.Spline1D(Alpha,Cm)
  Top_Xtr_spl = Dierckx.Spline1D(Alpha,Top_Xtr)
  Bot_Xtr_spl = Dierckx.Spline1D(Alpha,Bot_Xtr)

  # Finds the zero-lift AOA
  zerolift(aoa) = Cl_spl(aoa)
  aoaZeroCl = Roots.fzero(zerolift, -8,10)

  # Finds AOA of minimum drag
  step = 1e-6
  zeroMinDrag(aoa) = (Cd_spl(aoa) - Cd_spl(aoa+step))/step
  aoaMinDrag = Roots.fzero(zeroMinDrag, -8.0,10.0)
  cdmin = Cd_spl(aoaMinDrag)

  aoamin = minimum(Alpha)
  aoamax = maximum(Alpha)

  return (Cl_spl, Cd_spl, Cdp_spl, Cm_spl, Top_Xtr_spl, Bot_Xtr_spl,
          aoaZeroCl,aoaMinDrag, cdmin, aoamin, aoamax)
end
### END OF TOOLS ###############################################################
