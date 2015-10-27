/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef TENSORMECHANICSPLASTICJ2_H
#define TENSORMECHANICSPLASTICJ2_H

#include "TensorMechanicsPlasticModel.h"
#include "TensorMechanicsHardeningModel.h"


class TensorMechanicsPlasticJ2;


template<>
InputParameters validParams<TensorMechanicsPlasticJ2>();

/**
 * J2 plasticity, associative, with hardning.
 * Yield_function = sqrt(3*J2) - yield_strength
 */
class TensorMechanicsPlasticJ2 : public TensorMechanicsPlasticModel
{
 public:
  TensorMechanicsPlasticJ2(const InputParameters & parameters);

  /// returns the model name (J2)
  virtual std::string modelName() const;

 protected:

  /**
   * The yield function
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @return the yield function
   */
  Real yieldFunction(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of yield function with respect to stress
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @return df_dstress(i, j) = dyieldFunction/dstress(i, j)
   */
  RankTwoTensor dyieldFunction_dstress(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of yield function with respect to the internal parameter
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @return the derivative
   */
  Real dyieldFunction_dintnl(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The flow potential
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @return the flow potential
   */
  RankTwoTensor flowPotential(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of the flow potential with respect to stress
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @return dr_dstress(i, j, k, l) = dr(i, j)/dstress(k, l)
   */
  RankFourTensor dflowPotential_dstress(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of the flow potential with respect to the internal parameter
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @return dr_dintnl(i, j) = dr(i, j)/dintnl
   */
  RankTwoTensor dflowPotential_dintnl(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * YieldStrength.  The yield function is sqrt(3*J2) - yieldStrength.
   * In this class yieldStrength = 1, but this
   * may be over-ridden by derived classes with nontrivial hardning
   */
  virtual Real yieldStrength(const Real & intnl) const;

  /// d(yieldStrength)/d(intnl)
  virtual Real dyieldStrength(const Real & intnl) const;

  /// yield strength, from user input
  const TensorMechanicsHardeningModel & _strength;

  /// user-set flag: custom return map will run in a loop to get an exact solution & check yf < _f_tol
  const bool _update_strength;

  /// max iters for custom return map loop
  unsigned _max_iters;

  /// pre-made identity tensors
  const RankTwoTensor _iden;
  const RankFourTensor _iden4;

  /**
   * Performs a custom return-map.
   * You may choose to over-ride this in your
   * derived TensorMechanicsPlasticXXXX class,
   * and you may implement the return-map
   * algorithm in any way that suits you.  Eg, using
   * a Newton-Raphson approach, or a radial-return,
   * etc.
   * This may also be used as a quick way of ascertaining
   * whether (trial_stress, intnl_old) is in fact admissible.
   *
   * For over-riding this function, please note the
   * following three points:
   *
   * (1) In your implementation, you will want to set:
   *   returned_stress (the returned value of stress)
   *   returned_intnl  (the returned value of the internal variable)
   *   delta_dp   (the change in plastic strain)
   * The default implementation does not alter these quantities!
   *
   * (2) you must correctly signal "successful_return"
   * using the return value of this function.  Don't assume
   * the calling function will do Kuhn-Tucker checking and so forth!
   *
   * (3) if you return successful_return = false, you MUST
   * place the yield function values at (trial_stress, intnl_old)
   * into yf so the calling function can use this information
   * optimally.  You will have already calculated these
   * yield function values, which can be quite expensive, and
   * it's not very optimal for the calling function to have
   * to re-calculate them.
   *
   * @param trial_stress The trial stress
   * @param intnl_old Value of the internal parameter
   * @param E_ijkl Elasticity tensor
   * @param ep_plastic_tolerance Tolerance defined by the user for the plastic strain
   * @param[out] returned_stress Lies on the yield surface after returning and produces the correct plastic strain (normality condition)
   * @param[out] returned_intnl The value of the internal parameter after returning
   * @param[out] delta_dp The change in plastic strain induced by the return process
   * @param[out] The value of the yield function, evaluated at (trial_stress, intnl_old) in the case of return value = false (failure), or at (returned_stress, returned_intnl) for return value = true (success)
   * @param[out] trial_stress_inadmissible Should be set to 0 if the trial_stress is admissible, and 1 if the trial_stress is inadmissible.  This can be used by the calling prorgram
   * @return true if a successful return (or a return-map not needed), false if the trial_stress is inadmissible but the return process failed
   */
  bool returnMap(const RankTwoTensor & trial_stress, const Real & intnl_old, const RankFourTensor & E_ijkl, Real ep_plastic_tolerance,
                        RankTwoTensor & returned_stress, Real & returned_intnl, std::vector<Real> & model_pm, RankTwoTensor & delta_dp,
                        std::vector<Real> & yf, unsigned & trial_stress_inadmissible, const bool & update_pm) const;

  /**
    * Calculates a custom consistent tangent operator.
    * You may choose to over-ride this in your
    * derived TensorMechanicsPlasticXXXX class.
    *
    * @param stress current stress state
    * @param intnl internal parameters
    * @param E_ijkl elasticity tensor
    * @return the consistent tangent operator for J2 (radial return) case
    */
  RankFourTensor consistentTangentOperator(const RankTwoTensor & stress, const std::vector<Real> & intnl, const RankFourTensor & E_ijkl) const;

};

#endif // TENSORMECHANICSPLASTICJ2_H
