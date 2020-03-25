//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "ElementIntegralPostprocessor.h"

class AverageGBVelocity;

template <>
InputParameters validParams<AverageGBVelocity>();

/**
 * Calculate total grain boundary length in 2D and area in 3D.
 */
class AverageGBVelocity : public ElementIntegralPostprocessor
{
public:
  AverageGBVelocity(const InputParameters & parameters);

  virtual Real getValue() override;

protected:
  virtual Real computeQpIntegral() override;

  /// Number of order parameters
  const unsigned int _op_num;

  /// Order parameters
  std::vector<const VariableGradient *> _grads;

  /// normalization factor, depending on order parameter range and grains per side
  const Real _factor;
};
