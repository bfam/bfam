# ISSUES

## fields_m and fields_p
This storage needs to be removed from bfam_subdomain once dgx_quad is completely
removed

## Legendre-Gauss
Add the LG points to bfam.

## dump lua and include files
dump the lua file to the screen as well as any included lua files

## output handling with lua
Make it so that lua can set up different output types

## material projection stability
Check stability if projected material properties are used.


## communicator prefix
need to pass field prefix through to the get / put functions

## field init
Need to pass treeid through to the initialization function.

Possibly pass subdomain tags through to the field init function

## mesh
handle gmsh and split subdomains based on mesh tags as well as sudomain ID and
order and ???

## revisit output
consider using ADIOS for output

## output at points
Need to be able to output fields at a set of points (as opposed to just volume)

## move mesh points
call back function to allow a user to modify where the mesh points are

## beard fault grids
Currently things are not handled as cleanly in beard as they could with regards
to initial stresses on faults. Perhaps use a call back function or something
like that as in the old version
