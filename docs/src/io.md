# IO

If you want to store the constructed images you can save them as FITS files. This way you can open them in any code/program you want to make plots later.

## Saving

To save an image and the relevant fields from the [`mappingParameters`](@ref) struct you can use the function [`write_fits_image`](@ref):

```@docs
write_fits_image
```

The keyword arguments `units` and `snap` are optional and are used to store a unit string and the snapshot number for the image, respectively.


## Reading

To read a mapped image from a FITS file you can use [`read_fits_image`](@ref):

```@docs
read_fits_image
```