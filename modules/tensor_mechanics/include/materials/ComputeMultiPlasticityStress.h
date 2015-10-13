/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef COMPUTEMULTIPLASTICITYSTRESS_H
#define COMPUTEMULTIPLASTICITYSTRESS_H

#include "ComputeStressBase.h"
#include "MultiPlasticityDebugger.h"

class ComputeMultiPlasticityStress;

template<>
InputParameters validParams<ComputeMultiPlasticityStress>();

/**
 * ComputeMultiPlasticityStress performs the return-map
 * algorithm and associated stress updates for plastic
 * models defined by a General User Objects
 *
 * Note that if run in debug mode you might have to use
 * the --no-trap-fpe flag because PETSc-LAPACK-BLAS
 * explicitly compute 0/0 and 1/0, and this causes
 * Libmesh to trap the floating-point exceptions
 */
class ComputeMultiPlasticityStress :
  public ComputeStressBase,
  public MultiPlasticityDebugger
{
public:
  ComputeMultiPlasticityStress(const InputParameters & parameters);

protected:
  virtual void computeQpStress();
  virtual void initQpStatefulProperties();

  /// Maximum number of Newton-Raphson iterations allowed
  unsigned int _max_iter;

  /// Minimum fraction of applied strain that may be applied during adaptive stepsizing
  Real _min_stepsize;

  /// "dumb" deactivation will only be used if the stepsize falls below this quantity
  Real _max_stepsize_for_dumb;

  /// Even if the returnMap fails, return the best values found for stress and internal parameters
  bool _ignore_failures;

  /// Set to true in input file to use radial return map
  bool _radial_return;

  /// Set to true in input file to use radial return map consistent tangent operator
  /// this will eventually disappear --- just here to toggle for time comparison
  bool _radial_return_tan;

  /// The type of tangent operator to return.  tangent operator = d(stress_rate)/d(strain_rate).
  enum TangentOperatorEnum {
    elastic, linear, nonlinear
  } _tangent_operator_type;

  /// Tolerance on the plastic strain increment ("direction") constraint
  Real _epp_tol;

  /*
   * Scheme by which constraints are deactivated.
   * This enum is defined here for computational
   * efficiency.  If you add to this list you must
   * also add to the MooseEnum in the .C file
   */
  enum DeactivationSchemeEnum {
    optimized, safe, dumb, optimized_to_safe, safe_to_dumb,
    optimized_to_safe_to_dumb, optimized_to_dumb
  } _deactivation_scheme;


  /// User supplied the transverse direction vector
  bool _n_supplied;

  /// the supplied transverse direction vector
  RealVectorValue _n_input;

  /// rotation matrix that takes _n to (0, 0, 1)
  RealTensorValue _rot;


  /// plastic strain
  MaterialProperty<RankTwoTensor> & _plastic_strain;

  /// Old value of plastic strain
  MaterialProperty<RankTwoTensor> & _plastic_strain_old;

  /// internal parameters
  MaterialProperty<std::vector<Real> > & _intnl;

  /// old values of internal parameters
  MaterialProperty<std::vector<Real> > & _intnl_old;

  /// yield functions
  MaterialProperty<std::vector<Real> > & _yf;

  /// Number of Newton-Raphson iterations used in the return-map
  MaterialProperty<Real> & _iter;

  /// Whether a line-search was needed in the latest Newton-Raphson process (1 if true, 0 otherwise)
  MaterialProperty<Real> & _linesearch_needed;

  /// Whether linear-dependence was encountered in the latest Newton-Raphson process (1 if true, 0 otherwise)
  MaterialProperty<Real> & _ld_encountered;

  /// Whether constraints were added in during the latest Newton-Raphson process (1 if true, 0 otherwise)
  MaterialProperty<Real> & _constraints_added;

  /// current value of transverse direction
  MaterialProperty<RealVectorValue> & _n;

  /// old value of transverse direction
  MaterialProperty<RealVectorValue> & _n_old;

  /// strain increment (coming from ComputeIncrementalSmallStrain, for example)
  const MaterialProperty<RankTwoTensor> & _strain_increment;

  /// Old value of total strain (coming from ComputeIncrementalSmallStrain, for example)
  const MaterialProperty<RankTwoTensor> & _total_strain_old;

  /// Rotation increment (coming from ComputeIncrementalSmallStrain, for example)
  const MaterialProperty<RankTwoTensor> & _rotation_increment;

  /// Old value of stress
  MaterialProperty<RankTwoTensor> & _stress_old;

  /// Old value of elastic strain
  MaterialProperty<RankTwoTensor> & _elastic_strain_old;

  /// Elasticity tensor that can be rotated by this class (ie, its not const)
  ElasticityTensorR4 _my_elasticity_tensor;

  /// Strain increment that can be rotated by this class, and split into multiple increments (ie, its not const)
  RankTwoTensor _my_strain_increment;

  /**
   * makes all deactivated_due_to_ld false, and if >0 of them were initially true, returns true
   */
  virtual bool reinstateLinearDependentConstraints(std::vector<bool> & deactivated_due_to_ld);


  /**
   * counts the number of active constraints
   */
  virtual unsigned int numberActive(const std::vector<bool> & active);

  /**
   * Implements the return map
   *
   * Note that this algorithm doesn't do any rotations.  In order to find the
   * final stress and plastic_strain must be rotated using _rotation_increment.
   * This is usually done in computeQpStress
   *
   * @param stress_old The value of stress at the previous "time" step
   * @param stress (output) The stress after returning to the yield surface
   * @param intnl_old The internal variables at the previous "time" step
   * @param intnl    (output) All the internal variables after returning to the yield surface
   * @param plastic_strain_old The value of plastic strain at the previous "time" step
   * @param plastic_strain   (output) The value of plastic strain after returning to the yield surface
   * @param E_ijkl   The elasticity tensor.  If no plasiticity then stress = stress_old + E_ijkl*strain_increment
   * @param strain_increment   The applied strain increment
   * @param f  (output) All the yield functions after returning to the yield surface
   * @param iter (output) The number of Newton-Raphson iterations used
   * @param can_revert_to_dumb  If the _deactivation_scheme is set to revert to dumb, it will only be allowed to do so if this parameter is true
   * @param linesearch_needed (output) True if a linesearch was needed at any stage during the Newton-Raphson proceedure
   * @param ld_encountered (output) True if a linear-dependence of the flow directions was encountered at any stage during the Newton-Raphson proceedure
   * @param constraints_added (output) True if constraints were added into the active set at any stage during the Newton-Raphson proceedure
   * @param final_step Each strain increment may be decomposed into a sum of smaller increments if the return-map algorithm fails.  This flag indicates whether this is the last application of incremental strain
   * @param consistent_tangent_operator (output) The consistent tangent operator d(stress_rate)/d(strain_rate).  This is only output if final_step=true, and the return value of returnMap is also true.
   * @param cumulative_pm (input/output) Upon input: the plastic multipliers before the return map.  Upon output: the plastic multipliers after this return map, if the return map was successful
   * @return true if the stress was successfully returned to the yield surface
   */
  virtual bool returnMap(const RankTwoTensor & stress_old, RankTwoTensor & stress, const std::vector<Real> & intnl_old, std::vector<Real> & intnl, const RankTwoTensor & plastic_strain_old, RankTwoTensor & plastic_strain, const RankFourTensor & E_ijkl, const RankTwoTensor & strain_increment, std::vector<Real> & f, unsigned int & iter, const bool & can_revert_to_dumb, bool & linesearch_needed, bool & ld_encountered, bool & constraints_added, const bool & final_step, RankFourTensor & consistent_tangent_operator, std::vector<Real> & cumulative_pm);

  /**
   * Performs a single Newton-Raphson + linesearch step
   * Constraints are deactivated and the step is re-done if
   * deactivation_scheme is set appropriately
   * @param nr_res2 (input/output) Residual-squared that the line-search will reduce
   * @param stress (input/output) stress
   * @param intnl_old (input) old values of the internal parameters
   * @param intnl (input/output) internal parameters
   * @param pm (input/output) plastic multipliers
   * @param delta_dp (input/output) Change in plastic strain from start of "time" step to current configuration (plastic_strain - plastic_strain_old)
   * @param E_inv (input) Inverse of the elasticity tensor
   * @param f (input/output) Yield function(s).  Upon successful exit only the active constraints are contained in f
   * @param epp (input/output) Plastic strain increment constraint
   * @param ic (input/output) Internal constraint.  Upon successful exit only the active constraints are contained in ic
   * @param active The active constraints.  This is may be modified, depending upon deactivation_scheme
   * @param deactivation_scheme The scheme used for deactivating constraints
   * @param linesearch_needed (output) True if a linesearch was employed during this Newton-Raphson step
   * @param ld_encountered (output) True if a linear-dependence of the flow directions was encountered at any stage during the Newton-Raphson proceedure
   * @return true if the step was successful, ie, if the linesearch was successful and the number of constraints wasn't reduced to zero via deactivation
   */
  virtual bool singleStep(Real & nr_res2, RankTwoTensor & stress, const std::vector<Real> & intnl_old,
                          std::vector<Real> & intnl, std::vector<Real> & pm, RankTwoTensor & delta_dp,
                          const RankFourTensor & E_inv, const RankFourTensor & E_ijkl, std::vector<Real> & f,
                          RankTwoTensor & epp, std::vector<Real> & ic, std::vector<bool> & active,
                          DeactivationSchemeEnum deactivation_scheme,
                          bool & linesearch_needed, bool & ld_encountered);

  /**
   * Checks whether the yield functions are in the admissible region
   * @param stress stress
   * @param intnl internal parameters
   * @param all_f (output) the values of all the yield functions
   * @return return false if any yield functions exceed their tolerance
   */
  virtual bool checkAdmissible(const RankTwoTensor & stress, const std::vector<Real> & intnl,
                               std::vector<Real> & all_f);


  /**
   * Builds the order which "dumb" activation will take.
   * @param stress stress to evaluate yield functions and derivatives at
   * @param intnl internal parameters to evaluate yield functions and derivatives at
   * @param dumb_order (output) dumb_order[0] will be the yield surface furthest away from (stress, intnl), dumb_order[1] will be the next yield surface, etc.  The distance measure used is f/|df_dstress|.  This array can then be fed into incrementDumb in order to first try the yield surfaces which are farthest away from the (stress, intnl).
   */
  void buildDumbOrder(const RankTwoTensor & stress, const std::vector<Real> & intnl,
                      std::vector<unsigned int> & dumb_order);

  /**
   * Increments "dumb_iteration" by 1, and sets "act" appropriately
   * (act[alpha] = true iff alpha_th bit of dumb_iteration == 1)
   * @param (input/output) dumb_iteration Used to set act bitwise - the "dumb" scheme tries all possible combinations of act until a successful return
   * @param (output) act active constraints
   */
  virtual void incrementDumb(int & dumb_iteration, const std::vector<unsigned int> & dumb_order,
                             std::vector<bool> & act);

  /**
   * Checks Kuhn-Tucker conditions, and alters "active" if appropriate.
   * Do not let the simplicity of this routine fool you!
   * Explicitly:
   * (1) checks that pm = 0 for all the f < 0.  If not, then active is set to false for that constraint.  This may be triggered if upon exit of the NR loops a constraint got deactivated due to linear dependence, and then f<0 and its pm>0.
   * (2) checks that pm = 0 for all inactive constraints.  This should always be true unless someone has screwed with the code.
   * (3) if any pm < 0, then active is set to false for that constraint.  This may be triggered if _deactivation_scheme!="optimized".
   * @param f values of the active yield functions
   * @param pm values of all the plastic multipliers
   * @param active the active constraints (true if active)
   * @return return false if any of the Kuhn-Tucker conditions were violated (and hence the set of active constraints was changed)
   */
  virtual bool checkKuhnTucker(const std::vector<Real> & f, const std::vector<Real> & pm,
                               const std::vector<bool> & active);

  /**
   * Checks Kuhn-Tucker conditions, and alters "active" if appropriate.
   * Do not let the simplicity of this routine fool you!
   * Explicitly:
   * (1) checks that pm = 0 for all the f < 0.  If not, then active is set to false for that constraint.  This may be triggered if upon exit of the NR loops a constraint got deactivated due to linear dependence, and then f<0 and its pm>0.
   * (2) checks that pm = 0 for all inactive constraints.  This should always be true unless someone has screwed with the code.
   * (3) if any pm < 0, then active is set to false for that constraint.  This may be triggered if _deactivation_scheme!="optimized".
   * @param f values of the active yield functions
   * @param pm values of all the plastic multipliers
   * @param active the active constraints (true if active)
   * @return return false if any of the Kuhn-Tucker conditions were violated (and hence the set of active constraints was changed)
   */
  virtual void applyKuhnTucker(const std::vector<Real> & f, const std::vector<Real> & pm,
                               std::vector<bool> & active);

  // gets called before any return-map
  virtual void preReturnMap();

  // gets called after return-map
  virtual void postReturnMap();

  /*
   * performs an elastic step
   *
   * @param stress_old The value of stress at the previous "time" step
   * @param stress (output) stress = E_ijkl*plastic_strain
   * @param intnl_old The internal variables at the previous "time" step
   * @param intnl  (output) intnl = intnl_old
   * @param plastic_strain_old The value of plastic strain at the previous "time" step
   * @param plastic_strain   (output) plastic_strain = plastic_strain_old
   * @param E_ijkl   The elasticity tensor.
   * @param strain_increment   The applied strain increment
   * @param yf  (output) All the yield functions at (stress, intnl)
   * @param iterations (output) zero
   * @param consistent_tangent_operator (output) The consistent tangent operator d(stress_rate)/d(strain_rate)
   * @return true if the (stress, intnl) are admissible
   */
  virtual bool elasticStep(const RankTwoTensor & stress_old, RankTwoTensor & stress, const std::vector<Real> & intnl_old, std::vector<Real> & intnl, const RankTwoTensor & plastic_strain_old, RankTwoTensor & plastic_strain, const RankFourTensor & E_ijkl, const RankTwoTensor & strain_increment, std::vector<Real> & yf, unsigned int & iterations, RankFourTensor & consistent_tangent_operator);

  /*
   * performs a plastic step
   *
   * @param stress_old The value of stress at the previous "time" step
   * @param stress (output) stress after returning to the yield surface
   * @param intnl_old The internal variables at the previous "time" step
   * @param intnl  (output) internal variables after returning to the yield surface
   * @param plastic_strain_old The value of plastic strain at the previous "time" step
   * @param plastic_strain   (output) plastic_strain after returning to the yield surface
   * @param E_ijkl   The elasticity tensor.
   * @param strain_increment   The applied strain increment
   * @param yf  (output) All the yield functions at (stress, intnl)
   * @param iterations (output) The total number of Newton-Raphson iterations used
   * @param linesearch_needed (output) True if a linesearch was needed at any stage during the Newton-Raphson proceedure
   * @param ld_encountered (output) True if a linear-dependence of the flow directions was encountered at any stage during the Newton-Raphson proceedure
   * @param constraints_added (output) True if constraints were added into the active set at any stage during the Newton-Raphson proceedure
   * @param consistent_tangent_operator (output) The consistent tangent operator d(stress_rate)/d(strain_rate)
   * @return true if the (stress, intnl) are admissible.  Otherwise, if _ignore_failures==true, the output variables will be the best admissible ones found during the return-map.  Otherwise, if _ignore_failures==false, this routine will perform some finite-diference checks and call mooseError
   */
  virtual bool plasticStep(const RankTwoTensor & stress_old, RankTwoTensor & stress, const std::vector<Real> & intnl_old, std::vector<Real> & intnl, const RankTwoTensor & plastic_strain_old, RankTwoTensor & plastic_strain, const RankFourTensor & E_ijkl, const RankTwoTensor & strain_increment, std::vector<Real> & yf, unsigned int & iterations, bool & linesearch_needed, bool & ld_encountered, bool & constraints_added, RankFourTensor & consistent_tangent_operator);


  //  bool checkAndModifyConstraints(const bool & nr_exit_condition, const RankTwoTensor & stress, const std::vector<Real> & intnl, const std::vector<Real> & pm, const std::vector<bool> & initial_act, const bool & can_revert_to_dumb, const RankTwoTensor & initial_stress, const std::vector<Real> & intnl_old, const std::vector<Real> & f, DeactivationSchemeEnum deact_scheme, std::vector<bool> & act, int & dumb_iteration, std::vector<unsigned int> dumb_order, bool & die);

  bool canChangeScheme(DeactivationSchemeEnum current_deactivation_scheme, const bool & can_revert_to_dumb);

  bool canIncrementDumb(const int & dumb_iteration);

  void changeScheme(const std::vector<bool> & initial_act, const bool & can_revert_to_dumb,
                    const RankTwoTensor & initial_stress, const std::vector<Real> & intnl_old,
                    DeactivationSchemeEnum & current_deactivation_scheme, std::vector<bool> & act,
                    int & dumb_iteration, std::vector<unsigned int> & dumb_order);

  bool canAddConstraints(const std::vector<bool> & act, const std::vector<Real> & all_f);

  unsigned int activeCombinationNumber(const std::vector<bool> & act);

  /*
   * Computes the consistent tangent operator
   * (another name for the jacobian = d(stress_rate)/d(strain_rate)
   *
   * The computations performed depend upon _tangent_operator_type
   *
   * @param stress The value of stress after the return map algorithm has converged
   * @param intnl The internal parameters after the return map has converged
   * @param E_ijkl The elasticity tensor (in the case of no plasticity this is the jacobian)
   * @param pm_this_step The plastic multipliers coming from the final strain increment.  In many cases these will be equal to cumulative_pm, but in the case where the returnMap algorithm had to be performed in multiple substeps of smaller applied strain increments, pm_this_step are just the plastic multipliers for the final application of the strain incrment
   * @param cumulative_pm The plastic multipliers needed for this current Return (this is the sum of the plastic multipliers over all substeps if the strain increment was applied in small substeps)
   */
  RankFourTensor consistentTangentOperator(const RankTwoTensor & stress, const std::vector<Real> & intnl, const RankFourTensor & E_ijkl, const std::vector<Real> & pm_this_step, const std::vector<Real> & cumulative_pm);

};

#endif //COMPUTEMULTIPLASTICITYSTRESS_H
