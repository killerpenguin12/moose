//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "AverageGBVelocity.h"
#include "MooseUtils.h"

registerMooseObject("PhaseFieldApp", AverageGBVelocity);

template <>
InputParameters
validParams<AverageGBVelocity>()
{
  InputParameters params = validParams<ElementIntegralPostprocessor>();
  params.addClassDescription("Calculate average grain boundary velocity");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  params.addParam<Real>("",
                        2.0,
                        "Number of order parameters contacting a boundary "
                        "(should be 2.0 in polycrystals and 1.0 for "
                        "dispersed particles)");
  params.addParam<Real>("op_range",
                        1.0,
                        "Range over which order parameters change across an "
                        "interface. By default order parameters are assumed to "
                        "vary from 0 to 1");
  return params;
}

AverageGBVelocity::AverageGBVelocity(const InputParameters & parameters)
  : ElementIntegralPostprocessor(parameters),
    _op_num(coupledComponents("v")),
    _grads(_op_num),
    _factor(getParam<Real>("grains_per_side") * getParam<Real>("op_range"))
{
  // make sure user input is valid
  if (MooseUtils::absoluteFuzzyEqual(_factor, 0.0))
    mooseError("Neither grains_per_side nor op_range may be zero.");

  // Loop over variables (ops)
  for (MooseIndex(_op_num) op_index = 0; op_index < _op_num; ++op_index)
    _grads[op_index] = &coupledGradient("v", op_index);
}

Real
AverageGBVelocity::computeQpIntegral()
{
  Real grad_sum = 0.0;
  for (auto grad : _grads)
    grad_sum += (*grad)[_qp].norm();
  return grad_sum;
}

Real
AverageGBVelocity::getValue()
{
  return ElementIntegralPostprocessor::getValue() / _factor;
}
