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

  const TensorMechanicsHardeningModel & _strength;

   /// pre-made identity tensors
   RankTwoTensor _iden;
   RankFourTensor _iden4;

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
  bool returnMap(const RankTwoTensor & trial_stress, const Real & intnl_old, const RankFourTensor & E_ijkl, Real ep_plastic_tolerance, RankTwoTensor & returned_stress, Real & returned_intnl, RankTwoTensor & delta_dp, std::vector<Real> & yf, unsigned & trial_stress_inadmissible) const;

  void nrStep(const RankTwoTensor & stress, const std::vector<Real> & intnl_old, const std::vector<Real> & intnl, const std::vector<Real> & pm, const RankFourTensor & E_ijkl, RankTwoTensor & delta_dp, RankTwoTensor & dstress, std::vector<Real> & dpm, std::vector<Real> & dintnl, const std::vector<bool> & active, std::vector<bool> & deactivated_due_to_ld) const;

  RankFourTensor consistentTangentOperator(const RankTwoTensor & stress, const std::vector<Real> & intnl, const RankFourTensor & E_ijkl, const std::vector<Real> & pm_this_step, const std::vector<Real> & cumulative_pm) const;

  void calculateJacobian(const RankTwoTensor & stress, const std::vector<Real> & intnl, const std::vector<Real> & pm, const RankFourTensor & E_ijlk, const std::vector<bool> & active, const std::vector<bool> & deactivated_due_to_ld, std::vector<Real> & df_dintnl, std::vector<RankTwoTensor> & df_dstress, std::vector<Real> & jac, const Real & mu) const;

};

#endif // TENSORMECHANICSPLASTICJ2_H
