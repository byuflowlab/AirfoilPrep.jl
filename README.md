# AirfoilPrep.jl

Airfoil pre-processing tools, polar generation, and post-processing (splining of multi-dimensional lookup airfoil tables).
Includes calls to Xfoil, AirfoilPreppy, and the interpolations packages.


### Dependencies
- Python 3
   - matplotlib
   - mpmath
   - scipy
   
   
### Setting up PyCall

The airfoilprep.py package (wrapped by the AirfoilPrep.jl package) is written in Python 3.8, so make sure that the Python version linked to PyCall.jl is 3.8. After installing PyCall (] add PyCall), you can do this by running the following:

``` shell
import Pkg
Pkg.add("PyCall")
ENV["PYTHON"] = "path/to/your/python3"
Pkg.build("PyCall")
```

Then close and reopen the Julia REPL, and run:

``` shell
import PyCall
PyCall.pyversion
```

which should reveal your Python version:

``` shell
v"3.8"
```
