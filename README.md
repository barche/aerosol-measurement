# aerosol-measurement
Course materials for aerosol measurement

## Slides
The slides are available at https://barche.github.io/aerosol-intro-slides

## Notebooks installation

1. Download Julia 1.6 from https://julialang.org/downloads/
2. Start Julia and switch the prompt to package mode by hitting the `]` key, and then run:
```
(@v1.6) pkg> add https://github.com/Hinamoooon/jlmie.git
(@v1.6) pkg> develop https://github.com/barche/aerosol-measurement.git
(@v1.6) pkg> precompile
```
3. Exit package mode by hitting backspace and enter:
```julia
julia> using AerosolMeasurement
julia> cd(aerosoldir())
```
4. Enter Pkg mode again (hit `]`) and run (note the `.` at the end):
```
(@v1.6) pkg> activate .
```
5. Finally hit backspace and run:
```julia
julia> import Pluto
julia> Pluto.run()
```

Steps 1 and 2 are only needed at the initial setup. Your browser should now open, and you can load the notebooks by typing their name in the box, e.g. `mie.jl` (hit tab to get a list, the `.jl` files are Julia notebooks).
