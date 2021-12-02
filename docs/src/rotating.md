# Rotating images

By default `sphMapping` only maps the xy-plane. To change the mapping you have two options, [Project along axis](@ref) and [Define Euler Angle](@ref).

## Project along axis

If you only want to change the axis along which you want to project the data you can use the wrapper function [`project_along_axis`](@ref). 

```@docs
project_along_axis
```

If you for example want to project along the x-axes, so in the yz-plane use

```julia
axis = 1
pos_new = (pos_old, axis)
```

## Define Euler Angle

If projection along one of the principle axis is too crude for you, you can define individual angles α, β and γ corresponding to rotations around the x, y, and z-axis respectively and use the function [`rotate_3D`](@ref).

```@docs
rotate_3D
```

These angles have to be given in degrees.
So to rotate a 3D quantity 45 degrees around the x-axis you can use:

```julia
pos_new = rotate_3D(pos_old, 45.0, 0.0, 0.0)
```
