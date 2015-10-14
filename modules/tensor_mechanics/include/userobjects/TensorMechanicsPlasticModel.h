/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef TENSORMECHANICSPLASTICMODEL_H
#define TENSORMECHANICSPLASTICMODEL_H

#include "GeneralUserObject.h"
#include "RankTwoTensor.h"

class TensorMechanicsPlasticModel;


template<>
InputParameters validParams<TensorMechanicsPlasticModel>();

/**
 * Plastic Model base class
 * The virtual functions written below must be
 * over-ridden in derived classes to provide actual values
 *
 * It is assumed there is only one internal parameter, and
 * that is a function of the plastic multiplier, with rate
 * given by hardPotential
 *
 * For better or worse, I have created two versions of
 * all functions (eg yieldFunction, flowPotential, etc).
 * This is so that for single-surface plasticity you can
 * just override the 'protected' functions:
 * Real yieldFunction(stress, intnl)
 * (and similar), and don't have to worry about all the
 * multi-surface stuff, since in multi-surface yieldFunction
 * (etc) return std::vectors of stuff.  In the case of
 * multi-surface plasticity models you DO need to override the
 * 'public' functions (with a 'V' in their name):
 * void yieldFunctionV(stress, intnl, f)
 * versions.
 */
class TensorMechanicsPlasticModel : public GeneralUserObject
{
 public:
  TensorMechanicsPlasticModel(const InputParameters & parameters);

  void initialize();
  void execute();
  void finalize();

  /// The number of yield surfaces for this plasticity model
  virtual unsigned int numberSurfaces() const;

  /**
   * Calculates the yield functions.  Note that for single-surface plasticity
   * you don't want to override this - override the private yieldFunction below
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @param[out] f the yield functions
   */
  virtual void yieldFunctionV(const RankTwoTensor & stress, const Real & intnl, std::vector<Real> & f) const;

  /**
   * The derivative of yield functions with respect to stress
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @param[out] df_dstress df_dstress[alpha](i, j) = dyieldFunction[alpha]/dstress(i, j)
   */
  virtual void dyieldFunction_dstressV(const RankTwoTensor & stress, const Real & intnl, std::vector<RankTwoTensor> & df_dstress) const;

  /**
   * The derivative of yield functions with respect to the internal parameter
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @param[out] df_dintnl df_dintnl[alpha] = df[alpha]/dintnl
   */
  virtual void dyieldFunction_dintnlV(const RankTwoTensor & stress, const Real & intnl, std::vector<Real> & df_dintnl) const;

  /**
   * The flow potentials
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @param[out] r r[alpha] is the flow potential for the "alpha" yield function
   */
  virtual void flowPotentialV(const RankTwoTensor & stress, const Real & intnl, std::vector<RankTwoTensor> & r) const;

  /**
   * The derivative of the flow potential with respect to stress
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @param[out] dr_dstress dr_dstress[alpha](i, j, k, l) = dr[alpha](i, j)/dstress(k, l)
   */
  virtual void dflowPotential_dstressV(const RankTwoTensor & stress, const Real & intnl, std::vector<RankFourTensor> & dr_dstress) const;

  /**
   * The derivative of the flow potential with respect to the internal parameter
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @param[out] dr_dintnl  dr_dintnl[alpha](i, j) = dr[alpha](i, j)/dintnl
   */
  virtual void dflowPotential_dintnlV(const RankTwoTensor & stress, const Real & intnl, std::vector<RankTwoTensor> & dr_dintnl) const;

  /**
   * The hardening potential
   * @param stress the stress at which to calculate the hardening potential
   * @param intnl internal parameter
   * @param[out] h h[alpha] is the hardening potential for the "alpha" yield function
   */
  virtual void hardPotentialV(const RankTwoTensor & stress, const Real & intnl, std::vector<Real> & h) const;

  /**
   * The derivative of the hardening potential with respect to stress
   * @param stress the stress at which to calculate the hardening potentials
   * @param intnl internal parameter
   * @param[out] dh_dstress dh_dstress[alpha](i, j) = dh[alpha]/dstress(i, j)
   */
  virtual void dhardPotential_dstressV(const RankTwoTensor & stress, const Real & intnl, std::vector<RankTwoTensor> & dh_dstress) const;

  /**
   * The derivative of the hardening potential with respect to the internal parameter
   * @param stress the stress at which to calculate the hardening potentials
   * @param intnl internal parameter
   * @param[out] dh_dintnl dh_dintnl[alpha] = dh[alpha]/dintnl
   */
  virtual void dhardPotential_dintnlV(const RankTwoTensor & stress, const Real & intnl, std::vector<Real> & dh_dintnl) const;

  /**
   * The active yield surfaces, given a vector of yield functions.
   * This is used by FiniteStrainMultiPlasticity to determine the initial
   * set of active constraints at the trial (stress, intnl) configuration.
   * It is up to you (the coder) to determine how accurate you want the
   * returned_stress to be.  Currently it is only used by FiniteStrainMultiPlasticity
   * to estimate a good starting value for the Newton-Rahson procedure,
   * so currently it may not need to be super perfect.
   * @param f values of the yield functions
   * @param stress stress tensor
   * @param intnl internal parameter
   * @param Eijkl elasticity tensor (stress = Eijkl*strain)
   * @param[out] act act[i] = true if the i_th yield function is active
   * @param[out] returned_stress Approximate value of the returned stress
   */
  virtual void activeConstraints(const std::vector<Real> & f, const RankTwoTensor & stress, const Real & intnl, const RankFourTensor & Eijkl, std::vector<bool> & act, RankTwoTensor & returned_stress) const;

  /// Returns the model name (eg "MohrCoulom")
  virtual std::string modelName() const;

  /// Tolerance on yield function
  Real _f_tol;

  /// Tolerance on internal constraint
  Real _ic_tol;

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
  virtual bool returnMap(const RankTwoTensor & trial_stress, const Real & intnl_old, const RankFourTensor & E_ijkl, Real ep_plastic_tolerance, RankTwoTensor & returned_stress, Real & returned_intnl, RankTwoTensor & delta_dp, std::vector<Real> & yf, unsigned & trial_stress_inadmissible) const;

  virtual void nrStep(const RankTwoTensor & stress, const std::vector<Real> & intnl_old, const std::vector<Real> & intnl, const std::vector<Real> & pm, const RankFourTensor & E_ijkl, RankTwoTensor & delta_dp, RankTwoTensor & dstress, std::vector<Real> & dpm, std::vector<Real> & dintnl, const std::vector<bool> & active, std::vector<bool> & deactivated_due_to_ld) const;

  virtual void calculateJacobian(const RankTwoTensor & stress, const std::vector<Real> & intnl, const std::vector<Real> & pm, const RankFourTensor & E_ijlk, const std::vector<bool> & active, const std::vector<bool> & deactivated_due_to_ld, std::vector<Real> & df_dintnl, std::vector<RankTwoTensor> & df_dstress, std::vector<Real> & jac, const Real & mu) const;

  virtual RankFourTensor consistentTangentOperator(const RankTwoTensor & stress, const std::vector<Real> & intnl, const RankFourTensor & E_ijkl, const std::vector<Real> & pm_this_step, const std::vector<Real> & cumulative_pm) const;


 protected:

  /// The following functions are what you should override when building single-plasticity models
  /**
   * The yield function
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @return the yield function
   */
  virtual Real yieldFunction(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of yield function with respect to stress
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @return df_dstress(i, j) = dyieldFunction/dstress(i, j)
   */
  virtual RankTwoTensor dyieldFunction_dstress(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of yield function with respect to the internal parameter
   * @param stress the stress at which to calculate the yield function
   * @param intnl internal parameter
   * @return the derivative
   */
  virtual Real dyieldFunction_dintnl(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The flow potential
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @return the flow potential
   */
  virtual RankTwoTensor flowPotential(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of the flow potential with respect to stress
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @return dr_dstress(i, j, k, l) = dr(i, j)/dstress(k, l)
   */
  virtual RankFourTensor dflowPotential_dstress(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of the flow potential with respect to the internal parameter
   * @param stress the stress at which to calculate the flow potential
   * @param intnl internal parameter
   * @return dr_dintnl(i, j) = dr(i, j)/dintnl
   */
  virtual RankTwoTensor dflowPotential_dintnl(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The hardening potential
   * @param stress the stress at which to calculate the hardening potential
   * @param intnl internal parameter
   * @return the hardening potential
   */
  virtual Real hardPotential(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of the hardening potential with respect to stress
   * @param stress the stress at which to calculate the hardening potentials
   * @param intnl internal parameter
   * @return dh_dstress(i, j) = dh/dstress(i, j)
   */
  virtual RankTwoTensor dhardPotential_dstress(const RankTwoTensor & stress, const Real & intnl) const;

  /**
   * The derivative of the hardening potential with respect to the internal parameter
   * @param stress the stress at which to calculate the hardening potentials
   * @param intnl internal parameter
   * @return the derivative
   */
  virtual Real dhardPotential_dintnl(const RankTwoTensor & stress, const Real & intnl) const;

};

#endif // TENSORMECHANICSPLASTICMODEL_H
