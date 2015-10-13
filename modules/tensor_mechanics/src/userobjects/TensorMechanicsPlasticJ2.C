/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "TensorMechanicsPlasticJ2.h"

template<>
InputParameters validParams<TensorMechanicsPlasticJ2>()
{
  InputParameters params = validParams<TensorMechanicsPlasticModel>();
  params.addRequiredParam<UserObjectName>("yield_strength", "A TensorMechanicsHardening UserObject that defines hardening of the yield strength");
  params.addClassDescription("J2 plasticity, associative, with hardening");

  return params;
}

TensorMechanicsPlasticJ2::TensorMechanicsPlasticJ2(const InputParameters & parameters) :
    TensorMechanicsPlasticModel(parameters),
    _strength(getUserObject<TensorMechanicsHardeningModel>("yield_strength")),
    _iden(RankTwoTensor::initIdentity),
    _iden4(RankFourTensor::initIdentitySymmetricFour)
{
}


Real
TensorMechanicsPlasticJ2::yieldFunction(const RankTwoTensor & stress, const Real & intnl) const
{
  return std::pow(3*stress.secondInvariant(), 0.5) - yieldStrength(intnl);
}

RankTwoTensor
TensorMechanicsPlasticJ2::dyieldFunction_dstress(const RankTwoTensor & stress, const Real & /*intnl*/) const
{
  Real sII = stress.secondInvariant();
  if (sII == 0.0)
    return RankTwoTensor();
  else
    return 0.5*std::pow(3/sII, 0.5)*stress.dsecondInvariant();
}


Real
TensorMechanicsPlasticJ2::dyieldFunction_dintnl(const RankTwoTensor & /*stress*/, const Real & intnl) const
{
  return -dyieldStrength(intnl);
}

RankTwoTensor
TensorMechanicsPlasticJ2::flowPotential(const RankTwoTensor & stress, const Real & intnl) const
{
  return dyieldFunction_dstress(stress, intnl);
}

RankFourTensor
TensorMechanicsPlasticJ2::dflowPotential_dstress(const RankTwoTensor & stress, const Real & /*intnl*/) const
{
  Real sII = stress.secondInvariant();
  if (sII == 0)
    return RankFourTensor();

  RankFourTensor dfp = 0.5*std::pow(3/sII, 0.5)*stress.d2secondInvariant();
  Real pre = -0.25*std::pow(3, 0.5)*std::pow(sII, -1.5);
  RankTwoTensor dII = stress.dsecondInvariant();
  for (unsigned i = 0 ; i < 3 ; ++i)
    for (unsigned j = 0 ; j < 3 ; ++j)
      for (unsigned k = 0 ; k < 3 ; ++k)
        for (unsigned l = 0 ; l < 3 ; ++l)
          dfp(i, j, k, l) += pre*dII(i, j)*dII(k, l);
  return dfp;
}

RankTwoTensor
TensorMechanicsPlasticJ2::dflowPotential_dintnl(const RankTwoTensor & /*stress*/, const Real & /*intnl*/) const
{
  return RankTwoTensor();
}

Real
TensorMechanicsPlasticJ2::yieldStrength(const Real & intnl) const
{
  return _strength.value(intnl);
}

Real
TensorMechanicsPlasticJ2::dyieldStrength(const Real & intnl) const
{
  return _strength.derivative(intnl);
}

std::string
TensorMechanicsPlasticJ2::modelName() const
{
  return "J2";
}

void
TensorMechanicsPlasticJ2::nrStep(const RankTwoTensor & stress, const std::vector<Real> & intnl_old, const std::vector<Real> & intnl, const std::vector<Real> & pm, const RankFourTensor & E_ijkl, RankTwoTensor & delta_dp, RankTwoTensor & dstress, std::vector<Real> & dpm, std::vector<Real> & dintnl, const std::vector<bool> & active, std::vector<bool> & deactivated_due_to_ld) const
{
  // Calculate yield function and Jacobian (jac is scalar in this case)
  Real yf = yieldFunction(stress, intnl[0]);

  Real mu = E_ijkl(0,1,0,1);

  std::vector<Real> jac;
  jac.assign(1,0);
  std::vector<Real> df_dintnl;
  df_dintnl.assign(1,0);
  std::vector<RankTwoTensor> df_dstress;
  df_dstress.resize(1);
  calculateJacobian(stress, intnl, pm, E_ijkl, active, deactivated_due_to_ld, df_dintnl, df_dstress, jac, mu);

  dpm.assign(1,0);
  dpm[0] = yf / jac[0];
  //dpm = f / (3*mu - df_dintnl);

  dstress = -2 * mu * dpm[0] * df_dstress[0];
  dintnl.assign(1,0);
  dintnl[0] = dpm[0]; // only internal variable is hardening; q1 = pm
}

void
TensorMechanicsPlasticJ2::calculateJacobian(const RankTwoTensor & stress, const std::vector<Real> & intnl, const std::vector<Real> & pm, const RankFourTensor & E_ijlk, const std::vector<bool> & active, const std::vector<bool> & deactivated_due_to_ld, std::vector<Real> & df_dintnl, std::vector<RankTwoTensor> & df_dstress, std::vector<Real> & jac, const Real & mu) const
{
  // plastic modulus H = -df/dq = -df_dintnl
  // needed to calc dpm in nrStep()
  df_dintnl[0] = dyieldFunction_dintnl(stress, intnl[0]);

  // df/dstress = sqrt(3/2) * n_hat
  df_dstress[0] = dyieldFunction_dstress(stress, intnl[0]);

  jac[0] = 3*mu - df_dintnl[0];
}


RankFourTensor
TensorMechanicsPlasticJ2::consistentTangentOperator(const RankTwoTensor & stress, const std::vector<Real> & intnl, const RankFourTensor & E_ijkl, const std::vector<Real> & pm_this_step, const std::vector<Real> & cumulative_pm) const
{

  Real mu = E_ijkl(0,1,0,1);
  Real le = E_ijkl(1,1,0,0);
  Real K = le + (2.0*mu)/3.0; // bulk modulus, Hooke's law

  // plastic modulus H = -df/dq = -df_dintnl
  Real df_dintnl;
  df_dintnl = dyieldFunction_dintnl(stress, intnl[0]);

  // df/dstress = sqrt(3/2) * n_hat
  RankTwoTensor df_dstress;
  df_dstress = dyieldFunction_dstress(stress, intnl[0]);

  Real gamma = 1.0/(1.0 - (df_dintnl/(3.0*mu)));
  RankTwoTensor n_hat = (-df_dstress) / (std::sqrt(1.5));

  RankFourTensor iden2 = _iden.outerProduct(_iden);

  RankFourTensor iden_dev = _iden4 - ((1.0/3.0) * iden2);

  RankFourTensor n_hat2 = n_hat.outerProduct(n_hat);

  // need to calculate stress_dev in order to get effective_trial_stress
  RankTwoTensor dev_trial_stress(stress);
  dev_trial_stress.addIa(-dev_trial_stress.trace()/3.0);

  // compute effective trial stress = sqrt ( 3/2 s_dev : s_dev)
  Real dts_squared = dev_trial_stress.doubleContraction(dev_trial_stress);
  Real eff_trial_stress = std::sqrt(1.5 * dts_squared);

  Real beta = yieldStrength(intnl[0]) / eff_trial_stress;

  Real gamma_bar = gamma - (1.0 - beta);

  // consistent tangent operator for radial return map case
  return K * iden2 + 2.0*mu*beta * iden_dev - 2.0*mu*gamma_bar * n_hat2;
}
