//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralVariablePostprocessor.h"

/**
 * Calculate total grain boundary length in 2D and area in 3D.
 */
class AverageGBVelocity : public ElementIntegralVariablePostprocessor
{
public:
  static InputParameters validParams();

  AverageGBVelocity(const InputParameters & parameters);

protected:

  virtual void initialize() override;
  virtual void execute() override;
  virtual Real computeIntegral() override;
  virtual Real computeQpIntegral(){std::cout << "when does this run?" << std::endl;
    return 0;};

  virtual Real getValue(){
    std::cout << "qPoints: " << _qPoints << " value: " << _value << std::endl;
    return _value/_qPoints;};
  virtual void threadJoin(const UserObject & y) override;


  /// Number of order parameters
  const Real _eta_bottom;
  const Real _eta_top;
  const unsigned int _op_num;

  std::vector<const VariableValue *> _eta;
  std::vector<const VariableValue *> _deta_dt;
  std::vector<const VariableGradient *> _grad_eta;

  Real _value;
  Real _qPoints;
  /// Order parameters
  //std::vector<const VariableGradient *> _grads;

  /// normalization factor, depending on order parameter range and grains per side
//  const Real _factor;
};
