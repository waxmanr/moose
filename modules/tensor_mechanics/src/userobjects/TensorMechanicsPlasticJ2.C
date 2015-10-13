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
    _strength(getUserObject<TensorMechanicsHardeningModel>("yield_strength"))
{
}


Real
TensorMechanicsPlasticJ2::yieldFunction(const RankTwoTensor & stress, const std::vector<Real> & intnl) const
{
  return std::pow(3*stress.secondInvariant(), 0.5) - yieldStrength(intnl[0]);
}

RankTwoTensor
TensorMechanicsPlasticJ2::dyieldFunction_dstress(const RankTwoTensor & stress, const std::vector<Real> & /*intnl*/) const
{
  Real sII = stress.secondInvariant();
  if (sII == 0.0)
    return RankTwoTensor();
  else
    return 0.5*std::pow(3/sII, 0.5)*stress.dsecondInvariant();
}


std::vector<Real>
TensorMechanicsPlasticJ2::dyieldFunction_dintnl(const RankTwoTensor & /*stress*/, const std::vector<Real> & intnl) const
{
  std::vector<Real> dyf_dintnl(1,0);
  dyf_dintnl[0] = -dyieldStrength(intnl[0]);
  return dyf_dintnl;
}

RankTwoTensor
TensorMechanicsPlasticJ2::flowPotential(const RankTwoTensor & stress, const std::vector<Real> & intnl) const
{
  return dyieldFunction_dstress(stress, intnl);
}

RankFourTensor
TensorMechanicsPlasticJ2::dflowPotential_dstress(const RankTwoTensor & stress, const std::vector<Real> & /*intnl*/) const
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

std::vector<RankTwoTensor> // not sure why this fcn exists; overrides with same
TensorMechanicsPlasticJ2::dflowPotential_dintnl(const RankTwoTensor & /*stress*/, const std::vector<Real> & /*intnl*/) const
{
  std::vector<RankTwoTensor> dfl_dintnl;
  for (unsigned i = 0 ; numberICs() ; ++i)
    dfl_dintnl[i] = RankTwoTensor();
  return dfl_dintnl;
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
