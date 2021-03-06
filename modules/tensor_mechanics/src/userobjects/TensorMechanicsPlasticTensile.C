/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#include "TensorMechanicsPlasticTensile.h"
#include <math.h> // for M_PI

template<>
InputParameters validParams<TensorMechanicsPlasticTensile>()
{
  InputParameters params = validParams<TensorMechanicsPlasticModel>();
  params.addRequiredParam<UserObjectName>("tensile_strength", "A TensorMechanicsHardening UserObject that defines hardening of the tensile strength");
  params.addRangeCheckedParam<Real>("tensile_edge_smoother", 25.0, "tensile_edge_smoother>=0 & tensile_edge_smoother<=30", "Smoothing parameter: the edges of the cone are smoothed by the given amount.");
  MooseEnum tip_scheme("hyperbolic cap", "hyperbolic");
  params.addParam<MooseEnum>("tip_scheme", tip_scheme, "Scheme by which the pyramid's tip will be smoothed.");
  params.addRequiredRangeCheckedParam<Real>("tensile_tip_smoother", "tensile_tip_smoother>=0", "For the 'hyperbolic' tip_scheme, the pyramid vertex will be smoothed by the given amount.  For the 'cap' tip_scheme, additional smoothing will occur.  Typical value is 0.1*tensile_strength");
  params.addParam<Real>("cap_start", 0.0, "For the 'cap' tip_scheme, smoothing is performed in the stress_mean > cap_start region");
  params.addRangeCheckedParam<Real>("cap_rate", 0.0, "cap_rate>=0", "For the 'cap' tip_scheme, this controls how quickly the cap degenerates to a hemisphere: small values mean a slow degeneration to a hemisphere (and zero means the 'cap' will be totally inactive).  Typical value is 1/tensile_strength");
  params.addParam<Real>("tensile_lode_cutoff", "If the second invariant of stress is less than this amount, the Lode angle is assumed to be zero.  This is to gaurd against precision-loss problems, and this parameter should be set small.  Default = 0.00001*((yield_Function_tolerance)^2)");
  params.addClassDescription("Associative tensile plasticity with hardening/softening, and tensile_strength = 1");

  return params;
}

TensorMechanicsPlasticTensile::TensorMechanicsPlasticTensile(const InputParameters & parameters) :
    TensorMechanicsPlasticModel(parameters),
    _strength(getUserObject<TensorMechanicsHardeningModel>("tensile_strength")),
    _tip_scheme(getParam<MooseEnum>("tip_scheme")),
    _small_smoother2(std::pow(getParam<Real>("tensile_tip_smoother"), 2)),
    _cap_start(getParam<Real>("cap_start")),
    _cap_rate(getParam<Real>("cap_rate")),
    _tt(getParam<Real>("tensile_edge_smoother")*M_PI/180.0),
    _sin3tt(std::sin(3*_tt)),
    _lode_cutoff(parameters.isParamValid("tensile_lode_cutoff") ? getParam<Real>("tensile_lode_cutoff") : 1.0E-5*std::pow(_f_tol, 2))

{
  if (_lode_cutoff < 0)
    mooseError("tensile_lode_cutoff must not be negative");
  _ccc = (-std::cos(3*_tt)*(std::cos(_tt) - std::sin(_tt)/std::sqrt(3.0)) - 3*std::sin(3*_tt)*(std::sin(_tt) + std::cos(_tt)/std::sqrt(3.0)))/(18*std::pow(std::cos(3*_tt), 3));
  _bbb = (std::sin(6*_tt)*(std::cos(_tt) - std::sin(_tt)/std::sqrt(3.0)) - 6*std::cos(6*_tt)*(std::sin(_tt) + std::cos(_tt)/std::sqrt(3.0)))/(18*std::pow(std::cos(3*_tt), 3));
  _aaa = -std::sin(_tt)/std::sqrt(3.0) - _bbb*std::sin(3*_tt) - _ccc*std::pow(std::sin(3*_tt), 2) + std::cos(_tt);
}


Real
TensorMechanicsPlasticTensile::yieldFunction(const RankTwoTensor & stress, const Real & intnl) const
{
  Real mean_stress = stress.trace()/3.0;
  Real sin3Lode = stress.sin3Lode(_lode_cutoff, 0);
  if (sin3Lode <= _sin3tt)
  {
    // the non-edge-smoothed version
    std::vector<Real> eigvals;
    stress.symmetricEigenvalues(eigvals);
    return mean_stress + std::sqrt(smooth(stress) + std::pow(eigvals[2] - mean_stress, 2)) - tensile_strength(intnl);
  }
  else
  {
    // the edge-smoothed version
    Real kk = _aaa + _bbb*sin3Lode + _ccc*std::pow(sin3Lode, 2);
    Real sibar2 = stress.secondInvariant();
    return mean_stress + std::sqrt(smooth(stress) + sibar2*std::pow(kk, 2)) - tensile_strength(intnl);
  }
}


RankTwoTensor
TensorMechanicsPlasticTensile::dyieldFunction_dstress(const RankTwoTensor & stress, const Real & /*intnl*/) const
{
  Real mean_stress = stress.trace()/3.0;
  RankTwoTensor dmean_stress = stress.dtrace()/3.0;
  Real sin3Lode = stress.sin3Lode(_lode_cutoff, 0);
  if (sin3Lode <= _sin3tt)
  {
    // the non-edge-smoothed version
    std::vector<Real> eigvals;
    std::vector<RankTwoTensor> deigvals;
    stress.dsymmetricEigenvalues(eigvals, deigvals);
    Real denom = std::sqrt(smooth(stress) + std::pow(eigvals[2] - mean_stress, 2));
    return dmean_stress + (0.5*dsmooth(stress)*dmean_stress + (eigvals[2] - mean_stress)*(deigvals[2] - dmean_stress))/denom;
  }
  else
  {
    // the edge-smoothed version
    Real kk = _aaa + _bbb*sin3Lode + _ccc*std::pow(sin3Lode, 2);
    RankTwoTensor dkk = (_bbb + 2*_ccc*sin3Lode)*stress.dsin3Lode(_lode_cutoff);
    Real sibar2 = stress.secondInvariant();
    RankTwoTensor dsibar2 = stress.dsecondInvariant();
    Real denom = std::sqrt(smooth(stress) + sibar2*std::pow(kk, 2));
    return dmean_stress + (0.5*dsmooth(stress)*dmean_stress + 0.5*dsibar2*std::pow(kk, 2) + sibar2*kk*dkk)/denom;
  }
}


Real
TensorMechanicsPlasticTensile::dyieldFunction_dintnl(const RankTwoTensor & /*stress*/, const Real & intnl) const
{
  return -dtensile_strength(intnl);
}

RankTwoTensor
TensorMechanicsPlasticTensile::flowPotential(const RankTwoTensor & stress, const Real & intnl) const
{
  // This plasticity is associative so
  return dyieldFunction_dstress(stress, intnl);
}

RankFourTensor
TensorMechanicsPlasticTensile::dflowPotential_dstress(const RankTwoTensor & stress, const Real & /*intnl*/) const
{
  Real mean_stress = stress.trace()/3.0;
  RankTwoTensor dmean_stress = stress.dtrace()/3.0;
  Real sin3Lode = stress.sin3Lode(_lode_cutoff, 0);
  if (sin3Lode <= _sin3tt)
  {
    // the non-edge-smoothed version
    std::vector<Real> eigvals;
    std::vector<RankTwoTensor> deigvals;
    std::vector<RankFourTensor> d2eigvals;
    stress.dsymmetricEigenvalues(eigvals, deigvals);
    stress.d2symmetricEigenvalues(d2eigvals);

    Real denom = std::sqrt(smooth(stress) + std::pow(eigvals[2] - mean_stress, 2));
    Real denom3 = std::pow(denom, 3.0);
    RankTwoTensor numer_part = deigvals[2] - dmean_stress;
    RankTwoTensor numer_full = 0.5*dsmooth(stress)*dmean_stress + (eigvals[2] - mean_stress)*numer_part;
    Real d2smooth_over_denom = d2smooth(stress)/denom;

    RankFourTensor dr_dstress = (eigvals[2] - mean_stress)*d2eigvals[2]/denom;
    for (unsigned i = 0 ; i < 3 ; ++i)
      for (unsigned j = 0 ; j < 3 ; ++j)
        for (unsigned k = 0 ; k < 3 ; ++k)
          for (unsigned l = 0 ; l < 3 ; ++l)
          {
            dr_dstress(i, j, k, l) += 0.5*d2smooth_over_denom*dmean_stress(i, j)*dmean_stress(k, l);
            dr_dstress(i, j, k, l) += numer_part(i, j)*numer_part(k, l)/denom;
            dr_dstress(i, j, k, l) -= numer_full(i, j)*numer_full(k, l)/denom3;
          }
    return dr_dstress;
  }
  else
  {
    // the edge-smoothed version
    RankTwoTensor dsin3Lode = stress.dsin3Lode(_lode_cutoff);
    Real kk = _aaa + _bbb*sin3Lode + _ccc*std::pow(sin3Lode, 2);
    RankTwoTensor dkk = (_bbb + 2*_ccc*sin3Lode)*dsin3Lode;
    RankFourTensor d2kk = (_bbb + 2*_ccc*sin3Lode)*stress.d2sin3Lode(_lode_cutoff);
    for (unsigned i = 0 ; i < 3 ; ++i)
      for (unsigned j = 0 ; j < 3 ; ++j)
        for (unsigned k = 0 ; k < 3 ; ++k)
          for (unsigned l = 0 ; l < 3 ; ++l)
            d2kk(i, j, k, l) += 2*_ccc*dsin3Lode(i, j)*dsin3Lode(k, l);

    Real sibar2 = stress.secondInvariant();
    RankTwoTensor dsibar2 = stress.dsecondInvariant();
    RankFourTensor d2sibar2 = stress.d2secondInvariant();

    Real denom = std::sqrt(smooth(stress) + sibar2*std::pow(kk, 2));
    Real denom3 = std::pow(denom, 3.0);
    Real d2smooth_over_denom = d2smooth(stress)/denom;
    RankTwoTensor numer_full = 0.5*dsmooth(stress)*dmean_stress + 0.5*dsibar2*kk*kk + sibar2*kk*dkk;

    RankFourTensor dr_dstress = (0.5*d2sibar2*std::pow(kk, 2) + sibar2*kk*d2kk)/denom;
    for (unsigned i = 0 ; i < 3 ; ++i)
      for (unsigned j = 0 ; j < 3 ; ++j)
        for (unsigned k = 0 ; k < 3 ; ++k)
          for (unsigned l = 0 ; l < 3 ; ++l)
          {
            dr_dstress(i, j, k, l) += 0.5*d2smooth_over_denom*dmean_stress(i, j)*dmean_stress(k, l);
            dr_dstress(i, j, k, l) += (dsibar2(i, j)*dkk(k, l)*kk + dkk(i, j)*dsibar2(k, l)*kk + sibar2*dkk(i, j)*dkk(k, l))/denom;
            dr_dstress(i, j, k, l) -= numer_full(i, j)*numer_full(k, l)/denom3;
          }
    return dr_dstress;
  }
}

RankTwoTensor
TensorMechanicsPlasticTensile::dflowPotential_dintnl(const RankTwoTensor & /*stress*/, const Real & /*intnl*/) const
{
  return RankTwoTensor();
}


Real
TensorMechanicsPlasticTensile::tensile_strength(const Real internal_param) const
{
  return _strength.value(internal_param);
}

Real
TensorMechanicsPlasticTensile::dtensile_strength(const Real internal_param) const
{
  return _strength.derivative(internal_param);
}


Real
TensorMechanicsPlasticTensile::smooth(const RankTwoTensor & stress) const
{
  Real smoother2 = _small_smoother2;
  if (_tip_scheme == "cap")
  {
    Real x = stress.trace()/3.0 - _cap_start;
    Real p = 0;
    if (x > 0)
      p = x*(1 - std::exp(-_cap_rate*x));
    smoother2 += std::pow(p, 2);
  }
  return smoother2;
}


Real
TensorMechanicsPlasticTensile::dsmooth(const RankTwoTensor & stress) const
{
  Real dsmoother2 = 0;
  if (_tip_scheme == "cap")
  {
    Real x = stress.trace()/3.0 - _cap_start;
    Real p = 0;
    Real dp_dx = 0;
    if (x > 0)
    {
      p = x*(1 - std::exp(-_cap_rate*x));
      dp_dx = (1 - std::exp(-_cap_rate*x)) + x*_cap_rate*std::exp(-_cap_rate*x);
    }
    dsmoother2 += 2*p*dp_dx;
  }
  return dsmoother2;
}

Real
TensorMechanicsPlasticTensile::d2smooth(const RankTwoTensor & stress) const
{
  Real d2smoother2 = 0;
  if (_tip_scheme == "cap")
  {
    Real x = stress.trace()/3.0 - _cap_start;
    Real p = 0;
    Real dp_dx = 0;
    Real d2p_dx2 = 0;
    if (x > 0)
    {
      p = x*(1 - std::exp(-_cap_rate*x));
      dp_dx = (1 - std::exp(-_cap_rate*x)) + x*_cap_rate*std::exp(-_cap_rate*x);
      d2p_dx2 = 2*_cap_rate*std::exp(-_cap_rate*x) - x*std::pow(_cap_rate, 2)*std::exp(-_cap_rate*x);
    }
    d2smoother2 += 2*std::pow(dp_dx, 2) + 2*p*d2p_dx2;
  }
  return d2smoother2;
}

std::string
TensorMechanicsPlasticTensile::modelName() const
{
  return "Tensile";
}

