//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "AuxKernel.h"
//#include "Kernel.h"
//#include "JvarMapInterface.h"
#include "DerivativeMaterialInterface.h"
//#include "ACInterface.h"

// Forward Declarations
class calcVelocityAux;
class GrainTracker;

template <>
InputParameters validParams<calcVelocityAux>();

/**
 * Output euler angles from user object to an AuxVariable.
 */
class calcVelocityAux : public AuxKernel
{
public:
  calcVelocityAux(const InputParameters & parameters);

protected:
  virtual Real computeValue();
  virtual void precalculateValue();

  /// number of order parameters
  const Real _eta_bottom;
  const Real _eta_top;
  const unsigned int _op_num;



  std::vector<const VariableValue *> _eta;
  std::vector<const VariableValue *> deta_dt;
  std::vector<const VariableGradient *> _grad_eta;

  /// precalculated element value

  Real _value;
  Real new_val;
};
