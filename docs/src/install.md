# Install

As usual with Julia you can install the package via the internal package manager. Until I decide if should register the TSC interpolation package you need to install if by hand first (I don't know why).

```@example
julia> ]
pkg> add https://github.com/LudwigBoess/TriangularShapedCloudInterpolation.jl.git
pkg> add https://github.com/LudwigBoess/SPHtoGrid.jl
```

If you want to get the latest version use

```@example
julia> ]
pkg> add https://github.com/LudwigBoess/SPHtoGrid.jl#development
```


Now you should be good to go!