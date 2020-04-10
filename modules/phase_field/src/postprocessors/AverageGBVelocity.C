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

InputParameters
AverageGBVelocity::validParams()
{
  InputParameters params = ElementIntegralVariablePostprocessor::validParams();
  params.addClassDescription("Compute the average velocity of grain boundaries.");
  params.addParam<Real>("eta_bottom", 0.20, "start point on range of eta we are interested in.");
  params.addParam<Real>("eta_top", 0.80, "end point on range of eta we are interested in.");
  params.addRequiredCoupledVarWithAutoBuild(
      "v", "var_name_base", "op_num", "Array of coupled variables");
  return params;
}

AverageGBVelocity::AverageGBVelocity(const InputParameters & parameters)
  : ElementIntegralVariablePostprocessor(parameters),
    _eta_bottom(getParam<Real>("eta_bottom")),
    _eta_top(getParam<Real>("eta_top")),
    _op_num(coupledComponents("v")),
    _eta(_op_num),
    _deta_dt(_op_num),
    _grad_eta(_op_num),
    _qPoints(0)
{
  for (unsigned int i = 0; i < _op_num; ++i)
  {
    // Grab necessary variables
    _grad_eta[i] = &coupledGradient("v", i);
    _deta_dt[i] = &coupledDot("v", i);
    _eta[i] = &coupledValue("v", i);
  }
}

void
AverageGBVelocity::initialize()
{
  ElementIntegralVariablePostprocessor::initialize();
//  std::cout << "Iniatialize" << std::endl;
  _qPoints = 0;
  _value = 0;
}

void
AverageGBVelocity::execute()
{
  ElementIntegralVariablePostprocessor::execute();
/*  std::cout << "Execute" << std::endl;

  for(unsigned int _qp = 0; _qp < _qrule->n_points(); _qp++)
   for (unsigned int i = 0; i < _op_num; ++i)
   {
     if ((*_eta[i])[_qp] > _eta_bottom && (*_eta[i])[_qp] < _eta_top)
     {
       _qPoints += 1;
     }
   }*/

}

void
AverageGBVelocity::threadJoin(const UserObject & y)
{
  ElementIntegralVariablePostprocessor::threadJoin(y);
  //std::cout << "ThreadJoin" << std::endl;
  const AverageGBVelocity & pps = static_cast<const AverageGBVelocity &>(y);
  _qPoints += pps._qPoints;
  _value += pps._value;
}
Real
AverageGBVelocity::computeIntegral()
{
  //_value = 0;
  //Real _new_value = 0;
//std::cout << "Compute Integral" << std::endl;
for(unsigned int _qp = 0; _qp < _qrule->n_points(); _qp++)
  {
    for (unsigned int i = 0; i < _op_num; ++i)
    {
      _qPoints += 1;
      if ((*_eta[i])[_qp] > _eta_bottom && (*_eta[i])[_qp] < _eta_top)
      {
        _value += (1.0 / (*_grad_eta[i])[_qp].norm() * abs((*_deta_dt[i])[_qp]));

        //if((*_deta_dt[i])[_qp] < 0)
        //  std::cout << "_deta_dt is negative" << std::endl;
        //if(1.0 / (*_grad_eta[i])[_qp].norm() < 0)
        //  std::cout << "_grad_eta is negative" << std::endl;
      }
    }
  }
 //if(_value < 0) {
  //std::cout << " _value: " << _value << " qPoints: " << _qPoints << std::endl;}
  //std::cout << "IS THIS WORKING?? COME ON MAN!, _value: " << _value << std::endl;
  //if(  > 0)
  //  return sum/_qPoints;
  //else
  return _value;
}
