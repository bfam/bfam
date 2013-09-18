# ISSUES

## dump lua and include files
dump the lua file to the screen as well as any included lua files

## output handling with lua
Make it so that lua can set up different output types

## material projection stability
Check stability if projected material properties are used.

## boundary conditions
Fix it so that multiple boundary conditions can be handled

## communicator prefix
need to pass field prefix through to the get / put functions

## field init
Need to pass treeid through to the initialization function.

Possibly pass subdomain tags through to the field init function

## mesh
handle gmsh and split subdomains based on mesh tags as well as sudomain ID and
order and ???

## interfaces
Figure out how to handle internal interfaces

## dimensionality
generalize the code to be dimension independent. Make glue grids and subdomains
one uniform type and let all glue grids have many plus sides

## revisit output
consider using ADIOS for output

## output at points
Need to be able to output fields at a set of points (as opposed to just volume)

## move mesh points
call back function to allow a user to modify where the mesh points are
