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

#include "QuadraturePointMarker.h"
#include "FEProblem.h"
#include "MooseEnum.h"

template<>
InputParameters validParams<QuadraturePointMarker>()
{
  InputParameters params = validParams<Marker>();

  MooseEnum third_state("DONT_MARK=-1 COARSEN DO_NOTHING REFINE", "DONT_MARK");
  params.addParam<MooseEnum>("third_state", third_state, "The Marker state to apply to values falling in-between the coarsen and refine thresholds.");
  params.addParam<Real>("coarsen", "The threshold value for coarsening.  Elements with variable values beyond this will be marked for coarsening.");
  params.addParam<Real>("refine", "The threshold value for refinement.  Elements with variable values beyond this will be marked for refinement.");
  params.addParam<bool>("invert", false, "If this is true then values _below_ 'refine' will be refined and _above_ 'coarsen' will be coarsened.");
  params.addRequiredParam<VariableName>("variable", "The values of this variable will be compared to 'refine' and 'coarsen' to see what should be done with the element");
  return params;
}


QuadraturePointMarker::QuadraturePointMarker(const InputParameters & parameters) :
    Marker(parameters),
    Coupleable(parameters, false),
    MaterialPropertyInterface(parameters),
    _qrule(_assembly.qRule()),
    _q_point(_assembly.qPoints()),
    _qp(0)
{
  const std::vector<MooseVariable *> & coupled_vars = getCoupledMooseVars();
  for (unsigned int i=0; i<coupled_vars.size(); i++)
    addMooseVariableDependency(coupled_vars[i]);
}

Marker::MarkerValue
QuadraturePointMarker::computeElementMarker()
{
  MarkerValue current_mark = DONT_MARK;

  for (_qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    MarkerValue new_mark = computeQpMarker();

    current_mark = std::max(current_mark, new_mark);
  }

  return current_mark;
}

