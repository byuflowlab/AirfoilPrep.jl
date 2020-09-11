"""
A julia wrapper and interface for both XFOIL and airfoilprep.py, and other
functions for airfoil processing.
"""
module AirfoilPrep

# ------------ GENERIC MODULES -------------------------------------------------
using PyCall
using PyPlot
using JLD
using Dierckx
using Interpolations
using Roots
using LaTeXStrings
using Statistics: mean
import CSV
import DataFrames

# ------------ FLOW CODES ------------------------------------------------------
# Xfoil from https://github.com/byuflowlab/Xfoil.jl
import Xfoil

# ------------ GLOBAL VARIABLES ------------------------------------------------
const module_path = splitdir(@__FILE__)[1]          # Path to this module

# ------------ LOAD airfoilprep.py ---------------------------------------------
path_airfoilpreppy = module_path                    # Path to airfoilprep.py
prepy = PyNULL()                                    # airfoilpreppy module

function __init__()
    imp = pyimport("imp")
    (file, filename, data) = imp.find_module("airfoilprep", [path_airfoilpreppy])
    copy!(prepy, imp.load_module("airfoilprep", file, filename, data))
end

# ------------ HEADERS ---------------------------------------------------------
for header_name in ["pywrapper", "xfoil", "ndtools", "liftprops", "misc"]
    include("AirfoilPrep_"*header_name*".jl")
end

end # END OF MODULE
