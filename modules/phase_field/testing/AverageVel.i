#
# Test the parsed function free enery Allen-Cahn Bulk kernel
#

[Mesh]
  type = GeneratedMesh
  dim = 2 #dimensions of mesh
  nx = 100  #Number of elements in the x direction
  ny = 100 #number of elements in the Y direction
  xmax = 125 #x dimension of simulation
  ymax = 125
  elem_type = QUAD4
[]
[GlobalParams]
  op_num = 1
  var_name_base = gr
  #_dt = 1.00
[]
[AuxVariables]
  [./chi]
    [./InitialCondition]
      type = FunctionIC
      function = '0.5'
    [../]
  [../]
  [./velocity]
  order = CONSTANT
  family = MONOMIAL
  [../]
[]

[Variables]
#[./PolycrystalVariables]
  [./gr0]
    order = FIRST
    family = LAGRANGE
      [./InitialCondition]
      type = SmoothCircleIC
      x1 = 60.0
      y1 = 60.0
      radius = 30.0
      invalue = 1.0
      outvalue = 0.1
      int_width = 8.0
    [../]
  [../]
[]
[AuxKernels]
  [./velocity]
  type = calcVelocityAux
  variable = velocity
#  grain_num = 1
  [../]
[]
[Kernels]
  [./detadt]
    type = TimeDerivative
    variable = gr0
  [../]

  [./ACBulk]
    type = AllenCahn
    variable = gr0
    f_name = F
  [../]

  [./ACInterface]
    type = ACInterface
    variable = gr0
    kappa_name = 1
    variable_L = true
    args = chi
  [../]
[]

[Materials]
  [./L]
    type = DerivativeParsedMaterial
    f_name = L
    args = 'gr0 chi'
    function = '0.1 * gr0^2 + chi^2'
    derivative_order = 2
    #output = exodus
  [../]

  [./free_energy]
    type = DerivativeParsedMaterial
    f_name = F
    args = 'gr0'
    function = '2 * gr0^2 * (1-gr0)^2 - 0.2*gr0'
    derivative_order = 2
  [../]
[]
[Postprocessors]
  [./AvgVelocity]
    type = AverageGBVelocity
    variable = gr0
  [../]
[]
#  [./drdt]
#    type = GBVelocity
#  [../]
#  []

[Executioner]
  type = Transient
  scheme = 'bdf2'
  solve_type = 'NEWTON'
  num_steps = 100
  dt = 1.0
[]

[Outputs]
  exodus = true
  csv = true
[]
