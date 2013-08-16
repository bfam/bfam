# ISSUES

## communicator data
Data is now passed through to the communicator, but nobody (i.e. sbp or dgx)
does anything with it.

Still todo: Pass a struct through the pointer as outlined below.

Proposal in subdomain is a struct that contains char ** to a list of strings for
scalars, vectors, and tensors. Vectors and tensors will be rotated to normal and
perpendicular components. Default behavior will communicate all fields on the
minus and plus side if the pointer is NULL.

## communicator prefix
need to pass field prefix through to the get / put functions
