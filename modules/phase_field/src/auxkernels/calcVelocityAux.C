//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "calcVelocityAux.h"
#include "FEProblem.h"

registerMooseObject("PhaseFieldApp", calcVelocityAux);

template <>
InputParameters
validParams<calcVelocityAux>()
{
  InputParameters params = validParams<AuxKernel>();
  params.addClassDescription("calculate the velocity of grains while avoiding numerical instability.");
  params.addParam<Real>("eta_bottom",0.20,"start point on area of eta we are interested in.");
  params.addParam<Real>("eta_top",0.80,"end point on area of eta we are interested in.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

calcVelocityAux::calcVelocityAux(const InputParameters & parameters)
  : AuxKernel(parameters),
    _eta_bottom(getParam<Real>("eta_bottom")),
    _eta_top(getParam<Real>("eta_top")),
    _op_num(coupledComponents("v")),
    _eta(_op_num),
    deta_dt(_op_num),
    _grad_eta(_op_num)

{
  for (unsigned int i = 0; i < _op_num; ++i)
  {
    _grad_eta[i] = &coupledGradient("v", i);
    deta_dt[i] = &coupledDot("v",i);
    _eta[i] = &coupledValue("v",i);
  }
}


void
calcVelocityAux::precalculateValue()
{

}

Real
calcVelocityAux::computeValue()
{
  _value = 0;
  new_val = 0;
  for(unsigned int i = 0; i < _op_num; ++i)
  {
    if((*_eta[i])[_qp] > _eta_bottom && (*_eta[i])[_qp] < _eta_top){
      new_val = (1/((*_grad_eta[i])[_qp].norm())) * ((*deta_dt[i])[_qp]);}

    if(new_val > _value)
      _value = new_val;
  }

  return _value;
}
