# UserObject J2 test, with exponential hardening
# plot vm_stress vs intnl to see exponential relationship
# bcs and mesh designed to match j2_hard1_mod & LSH_mod for comparison

[Mesh]
  type = GeneratedMesh
  dim = 3
  nx = 4
  ny = 4
  nz = 4
  xmin = -0.5
  xmax = 0.5
  ymin = -0.5
  ymax = 0.5
  zmin = -0.5
  zmax = 0.5
[]


[Variables]
  [./disp_x]
  [../]
  [./disp_y]
  [../]
  [./disp_z]
  [../]
[]

[Kernels]
  [./TensorMechanics]
    displacements = 'disp_x disp_y disp_z'
  [../]
[]


[BCs]
  [./left]
    type = PresetBC
    variable = disp_x
    boundary = left
    value = 0.0
  [../]
  [./bottom]
    type = PresetBC
    variable = disp_y
    boundary = bottom
    value = 0.0
  [../]
  [./back]
    type = PresetBC
    variable = disp_z
    boundary = back
    value = 0.0
  [../]
  [./z]
    type = FunctionPresetBC
    variable = disp_z
    boundary = front
    function = 't/60'
  [../]
[]

[AuxVariables]
  [./stress_zz]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./intnl]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./vm_stress]
    order = CONSTANT
    family = MONOMIAL
  [../]
  [./eq_pl_strain]
    order = CONSTANT
    family = MONOMIAL
  [../]
[]

[AuxKernels]
  [./stress_zz]
    type = RankTwoAux
    rank_two_tensor = stress
    variable = stress_zz
    index_i = 2
    index_j = 2
  [../]
  [./intnl]
    type = MaterialStdVectorAux
    index = 0
    property = plastic_internal_parameter
    variable = intnl
  [../]
  [./vm_stress]
    type = RankTwoScalarAux
    rank_two_tensor = stress
    scalar_type = VonMisesStress
    variable = vm_stress
  [../]
  [./eq_pl_strain]
    type = RankTwoScalarAux
    rank_two_tensor = total_strain
    scalar_type = EquivalentPlasticStrain
    variable = eq_pl_strain
  [../]
[]

[Postprocessors]
  [./s_zz]
    type = PointValue
    point = '0 0 0'
    variable = stress_zz
  [../]
  [./intnl]
    type = PointValue
    point = '0 0 0'
    variable = intnl
  [../]
  [./vm_stress]
    type = PointValue
    point = '0 0 0'
    variable = vm_stress
  [../]
  [./eq_pl_strain]
    type = PointValue
    point = '0 0 0'
    variable = eq_pl_strain
  [../]
[]

[UserObjects]
  [./str]
    type = TensorMechanicsHardeningExponential
    value_0 = 100
    value_residual = 500
    rate = 100
  [../]
  [./j2]
    type = TensorMechanicsPlasticJ2
    yield_strength = str
    yield_function_tolerance = 1E-5
    internal_constraint_tolerance = 1E-9
    update_yield_strength = true
  [../]
[]

[Materials]
  [./elasticity_tensor]
    type = ComputeElasticityTensor
    block = 0
    fill_method = symmetric_isotropic
    C_ijkl = '0 1E6'
  [../]
  [./strain]
    type = ComputeFiniteStrain
    block = 0
    displacements = 'disp_x disp_y disp_z'
  [../]
  [./mc]
    type = ComputeMultiPlasticityStress
    block = 0
    ep_plastic_tolerance = 1E-9
    plastic_models = j2
  [../]
[]


[Executioner]
  type = Transient

  #Preconditioned JFNK (default)
  solve_type = 'PJFNK'

  petsc_options = '-snes_ksp_ew'
  petsc_options_iname = '-ksp_gmres_restart'
  petsc_options_value = '101'

  line_search = 'none'

  l_max_its = 100
  nl_max_its = 100
  nl_rel_tol = 1e-6
  nl_abs_tol = 1e-10
  l_tol = 1e-4
  start_time = 0.0
  end_time = 3.0
  dt = 0.1
[]


[Outputs]
  output_initial = true
  csv = true
  print_linear_residuals = false
  print_perf_log = true
[]
