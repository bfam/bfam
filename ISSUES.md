# ISSUES

## communicator data
pass a void * through communicator initialization that is then always passed to
the subdomain get / put buffers routines to allow for special data handling.

Also will need to get passed through time steppers.

Proposal in subdomain is a struct that contains char ** to a list of strings for
scalars, vectors, and tensors. Vectors and tensors will be rotated to normal and
perpendicular components.

## communicator prefix
need to pass field prefix through to the get / put functions
