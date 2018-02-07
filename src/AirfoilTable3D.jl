"""
A julia wrapper for creating 3D tables and splines of airfoil data.
"""


# ------------ GENERIC MODULES -------------------------------------------------
using Dierckx
using Interpolations
using PyPlot


# ------------ GLOBAL VARIABLES ------------------------------------------------
modulepath,_ = splitdir(@__FILE__)


# ------------ OTHER MODULES ------------------------------------------------
push!(LOAD_PATH, joinpath(modulepath,"jlxlight/julia/"))
using Xfoil
push!(LOAD_PATH, "$modulepath/")
import airfoilprep
import AirfoilProps
import Airfoil
import Kulfan

################################################################################
# 3D SPLINE WRAPPER FUNCTIONS
################################################################################
"""
  `gen3Dspl()`

Given an airfoil geometry, Mach, Re, it calculates its 2D aerodynamic properties calling
XFOIL and returns a Spline3D object.

  # ARGUMENTS
  * ` `     : .
  #EVERYTHING in DEGREES

**NOTE: Airfoil points must go from trailing edge around the top, then the
bottom and end back at the trailing edge. Also, Mach, Re, alpha must be
LINEAR for reverse lookup (map to 3D spline) to work**
"""
type Spline3D
  # Initialization variables (USER INPUT)

  spline_3D_cl#::Array{Float64,1}               # Lift coefficient at each AOA, Mach, Reynols.
  spline_3D_cd#::Array{Float64,1}               # Drag coefficient at each AOA, Mach, Reynols.
  spline_3D_cm#::Array{Float64,1}               # Moment coefficient at each AOA, Mach, Reynols.

  #y = mx+b or # = m*value+b (put in value, get out it's number that maps to the 3D spline)
  m_Re_array
  m_M_array
  m_alphas

  b_Re_array
  b_M_array
  b_alphas

  airfoil                                   # Aifoil name

  # Keep track of what it was run at
  Re_array
  M_array
  alphas

  Spline3D(spline_3D_cl,spline_3D_cd,spline_3D_cm,m_Re_array,m_M_array,m_alphas,
  b_Re_array,b_M_array,b_alphas,airfoil,Re_array,M_array,alphas) = new(spline_3D_cl,spline_3D_cd,
  spline_3D_cm,m_Re_array,m_M_array,m_alphas,b_Re_array,b_M_array,b_alphas,
  airfoil,Re_array,M_array,alphas)

end

type Table3D
  # Initialization variables (USER INPUT)

  table_3D_cl#::Array{Float64,1}               # Lift coefficient at each AOA, Mach, Reynols.
  table_3D_cd#::Array{Float64,1}               # Drag coefficient at each AOA, Mach, Reynols.
  table_3D_cm#::Array{Float64,1}               # Moment coefficient at each AOA, Mach, Reynols.

  airfoil                                   # Aifoil name

  # Keep track of what it was run at
  Re_array
  M_array
  grid_alphas

  # Include the original data as well
  # alpha_extrap #TODO:
  # cl

  Table3D(table_3D_cl,table_3D_cd,table_3D_cm,airfoil,Re_array,M_array,
  grid_alphas) = new(table_3D_cl,table_3D_cd,table_3D_cm,airfoil,Re_array,
  M_array,grid_alphas)

end

function Spline3D_from_table3D(table3D)
  #create the 3D spline
  table_3D_cl = table3D.table_3D_cl
  table_3D_cd = table3D.table_3D_cd
  table_3D_cm = table3D.table_3D_cm
  Re_array = table3D.Re_array
  M_array = table3D.M_array
  grid_alphas = table3D.grid_alphas
  airfoil = table3D.airfoil

  spline_3D_cl = interpolate(table_3D_cl, BSpline(Cubic(Line())), OnCell())
  spline_3D_cd = interpolate(table_3D_cd, BSpline(Cubic(Line())), OnCell())
  spline_3D_cm = interpolate(table_3D_cm, BSpline(Cubic(Line())), OnCell())

  #create splines for the reverse mapping of variables
  #y = mx+b or # = m*value+b (put in value, get out it's number)
  m_Re_array = (length(Re_array)-1)/(maximum(Re_array)-minimum(Re_array))
  m_M_array = (length(M_array)-1)/(maximum(M_array)-minimum(M_array))
  m_alphas = (length(grid_alphas)-1)/(maximum(grid_alphas)-minimum(grid_alphas))

  b_Re_array = Float64(1-m_Re_array*minimum(Re_array))
  b_M_array = Float64(1-m_M_array*minimum(M_array))
  b_alphas = Float64(1-m_alphas*minimum(grid_alphas))
  return Spline3D(spline_3D_cl,spline_3D_cd,spline_3D_cm,m_Re_array,m_M_array,
  m_alphas,b_Re_array,b_M_array,b_alphas,airfoil,Re_array,M_array,grid_alphas)
end

"""helper function for accessing 3D spline"""
function interp3(spline_3D,alpha,Re,M;case = "cl")
    if case=="cl"
        spl3D = spline_3D.spline_3D_cl
    elseif case=="cd"
        spl3D = spline_3D.spline_3D_cd
    elseif case=="cm"
        spl3D = spline_3D.spline_3D_cm
    else
        error("interpolation case not existent, use cl, cd, or cm")
    end

    num_alpha = spline_3D.m_alphas*alpha+spline_3D.b_alphas
    num_Re = spline_3D.m_Re_array*Re+spline_3D.b_Re_array
    num_M = spline_3D.m_M_array*M+spline_3D.b_M_array

    case_out = spl3D[num_alpha,num_Re,num_M]
    return case_out
end


function gen3Dspl(airfoil_name,Re_array,M_array;alphas=[i for i in -10:1.0:20],
  grid_alphas=[i for i in -180:1.0:180],TSR = 0.5,r_over_R = 0.3,
  c_over_r=0.18,iter=100,verbose=true, CDmax = 1.3,verify = false, returntable = true)

  # Make input robust
  Re_array = round.(Int, Re_array)
  M_array = collect(M_array)
  alphas = collect(alphas)
  # #1 aoa 2 M 3 RE 4
  # airfoil="./data/naca0006.dat"
  # #MUST be in LINEAR for reverse lookup (map to 3D spline) to work
  # Re_array = round.(Int,collect(linspace(1e3,3e6,3)))
  # Re_array = collect(linspace(.01,.5,3))
  # alphas = collect(linspace(-10,20,31))
  # grid_alphas = collect(linspace(-180,180,131))#[linspace(-180,minimum(alphas),50);alphas;linspace(maximum(alphas),180,50)]

  table_3D_cl = zeros(length(grid_alphas),length(Re_array),length(M_array))
  table_3D_cd = zeros(length(grid_alphas),length(Re_array),length(M_array))
  table_3D_cm = zeros(length(grid_alphas),length(Re_array),length(M_array))

  clgrid_alphas = []
  cdgrid_alphas = []
  cmgrid_alphas = []

  #Make this rhobust by applying the last solved for run, will change spline especially if switching from high to low on Re or M #TODO better way
  alpha_solved_last = []
  cl_solved_last = []
  cdf_solved_last = []
  cdp_solved_last = []
  cm_solved_last = []
  liftslope_last = []
  AoAzl_d_last = []
  clmin_last = []
  clmax1_last = []
  clmax2_last = []
  AoAclmin_last = []
  AoAclmax2_last = []

  alpha_solved = []
  cl_solved = []
  cdf_solved = []
  cdp_solved = []
  cm_solved = []
  liftslope = []
  AoAzl_d = []
  clmin = []
  clmax1 = []
  clmax2 = []
  AoAclmin = []
  AoAclmax2 = []

  for i = 1:length(Re_array)
    for j = 1:length(M_array)

      # x,y = airfoilprep.readcontour(airfoil; header_len=0)
      if false
        PyPlot.figure()
        PyPlot.plot(x,y)
        PyPlot.pause(0.001)
      end

      println("\n****** Running XFOIL...")
      # polar = airfoilprep.runXFOIL(x, y, Re_array[i]; verbose=verbose, Mach = M_array[j], iter=iter,alphas=alphas)

      Mach = M_array[j]
      Re = Re_array[i]
      coord = Airfoil.readAirfoil(airfoil_name)
      coord = Airfoil.fixAirfoil(coord) #gets it in the correct format for the kulfan parameterization method
      coord[:,2] -= coord[1,2] #TODO issue with kulfan implementation is that the trailing edge must be at y=0
      Al,Au=Kulfan.CoordtoCST(coord,3)


      try
        alpha_solved,cl_solved,cdf_solved,cdp_solved,cm_solved,liftslope,AoAzl_d,clmin,clmax1,clmax2,AoAclmin,AoAclmax2 = AirfoilProps.airfoilprops(;Al = Al,
                                Au = Au,
                                AoA_d = alphas,
                                re = Re,
                                # iter = 100,
                                # npan = 160,
                                # lincltol = 0.05,
                                # clmax2tol = lincltol,
                                verbose = true,
                                # plots = false,
                                # zeroinit = true,
                                percussive_maintenance = true,
                                # clmaxstop = false,
                                # clminstop = false,
                                # filter = true,
                                # filterwindow = 5,
                                # filtermaxstd = 2.0,
                                Mach = Mach
                                )
      catch
        println("####################################################")
        warn("Xfoil Failed at Mach $Mach Re $Re $airfoil_name \n skipping.")
        println("####################################################")
        alpha_solved = alpha_solved_last
        cl_solved = cl_solved_last
        cdf_solved = cdf_solved_last
        cdp_solved = cdp_solved_last
        cm_solved = cm_solved_last
        liftslope = liftslope_last
        AoAzl_d = AoAzl_d_last
        clmin = clmin_last
        clmax1 = clmax1_last
        clmax2 = clmax2_last
        AoAclmin = AoAclmin_last
        AoAclmax2 = AoAclmax2_last
      end
      #Reset
      alpha_solved_last = alpha_solved
      cl_solved_last = cl_solved
      cdf_solved_last = cdf_solved
      cdp_solved_last = cdp_solved
      cm_solved_last = cm_solved
      liftslope_last = liftslope
      AoAzl_d_last = AoAzl_d
      clmin_last = clmin
      clmax1_last = clmax1
      clmax2_last = clmax2
      AoAclmin_last = AoAclmin
      AoAclmax2_last = AoAclmax2

      cd_solved = cdf_solved+cdp_solved #TODO break these apart?

      polar = airfoilprep.Polar(Re, alpha_solved, cl_solved, cd_solved, cm_solved, coord[:,1], coord[:,2])

      # 3D corrected Polar
      newpolar = airfoilprep.correction3D(polar, r_over_R, c_over_r, TSR,
      alpha_linear_min=AoAclmin, alpha_linear_max=AoAclmax2, alpha_max_corr=maximum(alpha_solved))
      # Extrapolated polar
      extrap_polar = airfoilprep.extrapolate(newpolar, CDmax;nalpha = 40)

      cl2 = extrap_polar.init_cl
      cd2 = extrap_polar.init_cd
      cm2 = extrap_polar.init_cm
      alpha_extrap2 = extrap_polar.init_alpha

      #filter out duplicate aoa
      cl = []
      cd = []
      cm = []
      alpha_extrap = []
      alpha_extrap_prev = alpha_extrap2[1]
      for ii = 2:(length(alpha_extrap2))
        if alpha_extrap_prev==alpha_extrap2[ii]
          #skip it
        else
          push!(alpha_extrap,alpha_extrap2[ii])
          push!(cl,cl2[ii])
          push!(cd,cd2[ii])
          push!(cm,cm2[ii])
        end
        alpha_extrap_prev=alpha_extrap2[ii]
      end

      #Taylor's xfoil doesn't always return the specified aoa (percussive maintenance), and also sometimes drops points that absolutely fail
      clspl = Dierckx.Spline1D(alpha_extrap,cl,bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
      clgrid_alphas = clspl(grid_alphas)

      cdspl = Dierckx.Spline1D(alpha_extrap,cd,bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
      cdgrid_alphas = cdspl(grid_alphas)

      cmspl = Dierckx.Spline1D(alpha_extrap,cm,bc="nearest") # Evaluate at alphas that were specified to keep the 3D gridded spline consistent
      cmgrid_alphas = cmspl(grid_alphas)

      if verify
        PyPlot.figure()
        PyPlot.plot(grid_alphas,clgrid_alphas,".", label = "Extrapolated & Corrected Splined at Specified Alphas")
        PyPlot.plot(alpha_extrap,cl,"+", label = "Extrapolated & Corrected Alphas")
        PyPlot.plot(alpha_solved,cl_solved,"+", label = "XFOIL only Alphas")
        PyPlot.xlabel("AOA (deg)")
        PyPlot.xlabel("cl")
        PyPlot.legend()
      end
      table_3D_cl[:,i,j] = clgrid_alphas
      table_3D_cd[:,i,j] = cdgrid_alphas
      table_3D_cm[:,i,j] = cmgrid_alphas
    end
  end

  #populate table object with the info from the run
  table3D = Table3D(table_3D_cl,table_3D_cd,table_3D_cm,airfoil_name,Re_array,M_array,grid_alphas)
  #create the 3D spline
  spline3D = Spline3D_from_table3D(table3D)

  if verify == true && returntable==false

    #create splines for the reverse mapping of variables
    #y = mx+b or # = m*value+b (put in value, get out it's number)
    m_Re_array = (length(Re_array)-1)/(maximum(Re_array)-minimum(Re_array))
    m_M_array = (length(M_array)-1)/(maximum(M_array)-minimum(M_array))
    m_alphas = (length(grid_alphas)-1)/(maximum(grid_alphas)-minimum(grid_alphas))

    b_Re_array = Float64(1-m_Re_array*minimum(Re_array))
    b_M_array = Float64(1-m_M_array*minimum(M_array))
    b_alphas = Float64(1-m_alphas*minimum(grid_alphas))

    numRe_array = collect(linspace(1,length(Re_array),length(Re_array)))
    nummach = collect(linspace(1,length(M_array),length(M_array)))
    numalpha = collect(linspace(1,length(grid_alphas),length(grid_alphas)))
    PyPlot.figure()
    PyPlot.plot(Re_array,(m_Re_array*Re_array+b_Re_array),label = "Re_array fit")
    PyPlot.plot(Re_array,numRe_array,label = "actual nums")
    PyPlot.xlabel("Alpha")
    PyPlot.ylabel("#")
    PyPlot.legend()

    PyPlot.figure()
    PyPlot.plot(M_array,(m_M_array*M_array+b_M_array),label = "alpha fit")
    PyPlot.plot(M_array,nummach,label = "actual nums")
    PyPlot.xlabel("Alpha")
    PyPlot.ylabel("#")
    PyPlot.legend()

    PyPlot.figure()
    PyPlot.plot(grid_alphas,(m_alphas*grid_alphas+b_alphas),label = "alpha fit")
    PyPlot.plot(grid_alphas,numalpha,label = "actual alpha")
    PyPlot.xlabel("Alpha")
    PyPlot.ylabel("#")
    PyPlot.legend()

    PyPlot.figure()
    for k = 1:length(Re_array)
      for j = 1:length(M_array)
        cdtest = zeros(cdgrid_alphas)
        for i =1:length(cdgrid_alphas)
          # i = 1
          # j = 1
          # k = 1
          #test alpha spline
          alphatest = grid_alphas[i]
          num_alpha = m_alphas*alphatest+b_alphas

          #test re spline
          Retest = Re_array[k]
          num_Re = m_Re_array*Retest+b_Re_array

          #test mach spline
          Mtest = M_array[j]
          num_M = m_M_array*Mtest+b_M_array

          cdtest[i] = spline3D.spline_3D_cd[num_alpha,num_Re,num_M]
        end
        PyPlot.plot(grid_alphas,cdtest,"r.-",label = "interp cd")
        PyPlot.plot(grid_alphas,table_3D_cd[:,k,j],"b-",label = "real cd")
        PyPlot.legend()
        PyPlot.xlabel("AOA")
        PyPlot.ylabel("cd")
      end
    end
    return spline3D
  end


  if returntable==false
    return spline3D
  else
    return table3D
  end
end #gen3Dspl()

### END OF 3D spline ###############################################################
