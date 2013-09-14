-- refinement parameters
min_refine_level = 2
max_refine_level = 4

-- connectivity info
connectivity = "brick"
brick = {
  nx = 2,
  ny = 3,
  periodic_x = 1,
  periodic_y = 1,
}

-- time stepper to use
lsrk_method  = "KC54"

