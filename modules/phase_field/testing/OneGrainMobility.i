#
# Test the parsed function free enery Allen-Cahn Bulk kernel
#

[Mesh]
  type = GeneratedMesh
  dim = 2 #dimensions of mesh
  nx = 60  #Number of elements in the x direction
  ny = 60 #number of elements in the Y direction
  xmax = 125 #x dimension of simulation
  ymax = 125
  elem_type = QUAD4
[]

[AuxVariables]
  [./chi]
    [./InitialCondition]
      type = FunctionIC
      function = 'x/48+0.5'
    [../]
  [../]
[]

[Variables]
  [./eta]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = SmoothCircleIC
      x1 = 60.0
      y1 = 60.0
      radius = 20.0
      invalue = 1.0
      outvalue = 0.1
      int_width = 8.0
    [../]
  [../]
[]

[Kernels]
  [./detadt]
    type = TimeDerivative
    variable = eta
  [../]

  [./ACBulk]
    type = AllenCahn
    variable = eta
    f_name = F
  [../]

  [./ACInterface]
    type = ACInterface
    variable = eta
    kappa_name = 1
    variable_L = true
    args = chi
  [../]
[]

[Materials]
  [./L]
    type = DerivativeParsedMaterial
    f_name = L
    args = 'eta chi'
    function = '0.1 * eta^2 + chi^2'
    derivative_order = 2
  [../]

  [./free_energy]
    type = DerivativeParsedMaterial
    f_name = F
    args = 'eta'
    function = '2 * eta^2 * (1-eta)^2 - 0.2*eta'
    derivative_order = 2
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'
  solve_type = 'NEWTON'
  num_steps = 100
  dt = 1
[]

[Outputs]
  exodus = true
[]
