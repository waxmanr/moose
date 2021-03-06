/****************************************************************/
/*               DO NOT MODIFY THIS HEADER                      */
/* MOOSE - Multiphysics Object Oriented Simulation Environment  */
/*                                                              */
/*           (c) 2010 Battelle Energy Alliance, LLC             */
/*                   ALL RIGHTS RESERVED                        */
/*                                                              */
/*          Prepared by Battelle Energy Alliance, LLC           */
/*            Under Contract No. DE-AC07-05ID14517              */
/*            With the U. S. Department of Energy               */
/*                                                              */
/*            See COPYRIGHT for full restrictions               */
/****************************************************************/

#include "ExampleFunction.h"

template<>
InputParameters validParams<ExampleFunction>()
{
  InputParameters params = validParams<Function>();
  params.addParam<Real>("alpha", 1.0, "The value of alpha");
  return params;
}

ExampleFunction::ExampleFunction(const InputParameters & parameters) :
    Function(parameters),
    _alpha(getParam<Real>("alpha"))
{}

Real
ExampleFunction::value(Real /*t*/, const Point & p)
{
  return _alpha*_alpha*libMesh::pi*libMesh::pi*std::sin(_alpha*libMesh::pi*p(0));  // p(0) == x
}
