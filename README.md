# AirfoilPrep.jl

Collection of airfoil pre-processing routines for creating splined airfoil data over AOA, RE, Mach.  Specify airfoil geometry and operating conditions (including rotor parameters if a rotor) and receive a 3D table of airfoil lift and drag data.  Includes calls to Xfoil, AirfoilPreppy, and the interpolations packages.  

### To Install

1. Clone the repo:
```julia
Pkg.clone(...)
```

## Run tests

```julia
Pkg.test("AirfoilPrep")
```

## To Use

```julia
using AirfoilPrep
```

See examples in tests.


[![Build Status](https://travis-ci.org/moore54/AirfoilPrep.jl.svg?branch=master)](https://travis-ci.org/moore54/AirfoilPrep.jl)

[![Coverage Status](https://coveralls.io/repos/moore54/AirfoilPrep.jl/badge.svg?branch=master&service=github)](https://coveralls.io/github/moore54/AirfoilPrep.jl?branch=master)

[![codecov.io](http://codecov.io/github/moore54/AirfoilPrep.jl/coverage.svg?branch=master)](http://codecov.io/github/moore54/AirfoilPrep.jl?branch=master)
