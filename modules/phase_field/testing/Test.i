[Mesh]
  type = GeneratedMesh
  dim = 2
  nx = 10
  ny = 10
  nz = 0
  xmax = 1000
  ymax = 1000
  zmax = 0
  elem_type = QUAD4
  uniform_refine = 2
[]

[GlobalParams]
  op_num = 3
  var_name_base = 'gr'
[]

[Variables]
  [./gr0]
    [./InitialCondition]
      type = Tricrystal2CircleGrainsIC
      op_index = 0
      op_num = 3
      variable = gr0
    [../]
  [../]
    [./gr1]
      [./InitialCondition]
        type = Tricrystal2CircleGrainsIC
        op_index = 1
        op_num = 3
        variable = gr1
      [../]
    [../]
      [./gr2]
        [./InitialCondition]
          type = Tricrystal2CircleGrainsIC
          op_index = 2
          op_num = 3
          variable = gr2
        [../]
      [../]
[]
#[UserObjects]
#  [./voronoi]
#    type = PolycrystalVoronoi
#    rand_seed = 102
#    grain_num = 3
#    coloring_algorithm = bt
#  [../]
#[]

[AuxVariables]
  [./bnds]
    order = FIRST
    family = LAGRANGE
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
    execute_on = 'timestep_end'
  [../]
[]

[BCs]
  [./Periodic]
    [./All]
      auto_direction = 'x y'
    [../]
  [../]
[]

[Materials]
  [./Mobil0]
    type = GBEvolution
    time_scale = 1.0
    GBmob0 = 3.986e-6
    T = 500 # K
    wGB = 60 # nm
    Q = 1.0307
    GBenergy = 2.4
  [../]
[]

[Postprocessors]   # grain 0 is theback grain, grain 1 is the left and 2 is right.
  [./gr1area]
    type = ElementIntegralVariablePostprocessor
    variable = gr1
    execute_on = 'initial timestep_end'
  [../]
  [./avg_grain_vol]
    type = AverageGrainVolume
    grain_num = 3
    execute_on = 'initial timestep_end'
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
  scheme = bdf2
  solve_type = 'NEWTON'

  petsc_options_iname = '-pc_type -pc_hypre_type -ksp_gmres_restart'
  petsc_options_value = 'hypre boomeramg 31'

  l_tol = 1.0e-4
  l_max_its = 30
  nl_max_its = 20
  nl_rel_tol = 1.0e-9
  start_time = 0.0
  num_steps = 25
  dt = 2
[]

[Outputs]
  exodus = true
  csv = true
[]
