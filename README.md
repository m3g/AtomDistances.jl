# AtomDistances.jl

Simple script to illustrate how to compute distance-dependent properties using `Chemfiles.jl`, `PDBTools.jl`, etc.
The idea is that the user will modify the script to compute whatever property he/she wants.

## How to use it

```julia
julia> using Revise # recommended to be in .julia/config/startup.jl

julia> includet("distance.jl") # note the "t"

julia> plt, distances = run_namd()
```

Next, modify the `run_namd` and/or `run_gromacs` functions to compute the properties desired. Also, comment the line
that adds the package to avoid reinstalling or updating the packages on every script execution. Execute again the
functions (the changes will be tracked by the `Revise` package).








