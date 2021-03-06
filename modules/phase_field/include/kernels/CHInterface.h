/****************************************************************/
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*          All contents are licensed under LGPL V2.1           */
/*             See LICENSE for full restrictions                */
/****************************************************************/
#ifndef CHINTERFACE_H
#define CHINTERFACE_H

#include "Kernel.h"
#include "Material.h"
#include "JvarMapInterface.h"
#include "DerivativeKernelInterface.h"

//Forward Declarations
class CHInterface;

template<>
InputParameters validParams<CHInterface>();

/**
 * This is the Cahn-Hilliard equation base class that implements the interfacial or gradient energy term of the equation.
 * See M.R. Tonks et al. / Computational Materials Science 51 (2012) 20–29 for more information.
 */
class CHInterface : public DerivativeMaterialInterface<JvarMapInterface<Kernel> >
{
public:
  CHInterface(const InputParameters & parameters);

protected:
  virtual Real computeQpResidual();
  virtual Real computeQpJacobian();
  virtual Real computeQpOffDiagJacobian(unsigned int jvar);

  const MaterialProperty<Real> & _kappa;

  ///@{
  /// Mobility material property value and concentration derivatives
  const MaterialProperty<Real> & _M;
  const MaterialProperty<Real> & _dMdc;
  const MaterialProperty<Real> & _d2Mdc2;
  ///@}

  ///@{
  /// Variables for second order derivatives
  VariableSecond & _second_u;
  VariableTestSecond & _second_test;
  VariablePhiSecond & _second_phi;
  ///@}

  ///Number of variables
  unsigned int _nvar;

  ///@{
  /// Mobility derivatives w.r.t. its dependent variables
  std::vector<const MaterialProperty<Real> *> _dMdarg;
  std::vector<const MaterialProperty<Real> *> _d2Mdcdarg;
  std::vector<std::vector<const MaterialProperty<Real>* > > _d2Mdargdarg;
  ///@}

  /// Coupled variables used in mobility
  std::vector<VariableGradient *> _coupled_grad_vars;
};

#endif //CHINTERFACE_H
