[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 50
  ny = 50
  nz = 0
  xmin = 0
  xmax = 1000
  ymin = 0
  ymax = 1000
  zmin = 0
  zmax = 0
  elem_type = QUAD4

  uniform_refine = 2
[]

[GlobalParams]
  op_num = 2
  var_name_base = gr
  v = 'gr0 gr1'
[]

[Variables]
  [./gr0]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ThumbIC
      xcoord = 500.0
      height = 600.0
      width = 400.0
      invalue = 0.0
      outvalue = 1.0
    [../]
  [../]

  [./gr1]
    order = FIRST
    family = LAGRANGE
    [./InitialCondition]
      type = ThumbIC
      xcoord = 500.0
      height = 600.0
      width = 400.0
      invalue = 1.0
      outvalue = 0.0
    [../]
  [../]
[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
  [../]
  [./velocity]
  order = CONSTANT
  family = MONOMIAL
  [../]
[]

[Kernels]
  [./PolycrystalKernel]
  [../]
[]

[AuxKernels]
  [./BndsCalc]
    type = BndsCalcAux
    variable = bnds
  [../]
  [./velocity]
    type = calcVelocityAux
    variable = velocity
  [../]
[]

[BCs]
  active = ' '

  [./Periodic]
    [./left_right]
      primary = 0
      secondary = 2
      translation = '0 1000 0'
    [../]

    [./top_bottom]
      primary = 1
      secondary = 3
      translation = '-1000 0 0'
    [../]
  [../]
[]

[Materials]
  [./Copper]
    type = GBEvolution
    T = 500 # K
    wGB = 60 # nm
    GBmob0 = 2.5e-6 #m^4/(Js) from Schoenfelder 1997
    Q = 0.23 #Migration energy in eV
    GBenergy = 0.708 #GB energy in J/m^2
  [../]
[]

[Postprocessors]
  [./gr_area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
  [../]
  [./pos]
    type = FindValueOnLine
    target = 0.5
    v = gr1
    start_point = '500 500 0'
    end_point = '500 1000 0'
    #execute_on = 'initial timestep_end'
  [../]
[]

[Preconditioning]
  [./SMP]
   type = SMP
   full = true
  [../]
[]

[Executioner]
  type = Transient
  scheme = 'bdf2'

  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  l_tol = 1.0e-4
  l_max_its = 30

  nl_max_its = 20
  nl_rel_tol = 1.0e-9

  start_time = 0.0
  num_steps = 50
  dt = 80.0

  [./Adaptivity]
    initial_adaptivity = 2
    refine_fraction = 0.8
    coarsen_fraction = 0.05
    max_h_level = 2
  [../]
[]

[Outputs]
  execute_on = 'timestep_end'
  csv = true
  exodus = true
[]
