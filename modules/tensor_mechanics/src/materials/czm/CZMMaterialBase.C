//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "Assembly.h"
#include "CZMMaterialBase.h"
#include "RotationMatrix.h"

template <>
InputParameters
validParams<CZMMaterialBase>()
{
  InputParameters params = validParams<InterfaceMaterial>();

  params.addClassDescription("Base class for cohesive zone mateirla models");
  params.addRequiredCoupledVar("displacements",
                               "The string of displacements suitable for the problem statement");
  return params;
}

CZMMaterialBase::CZMMaterialBase(const InputParameters & parameters)
  : InterfaceMaterial(parameters),
    _normals(_assembly.normals()),
    _ndisp(coupledComponents("displacements")),
    _disp(_ndisp),
    _disp_neighbor(_ndisp),
    _displacement_jump_global(declareProperty<RealVectorValue>("displacement_jump_global")),
    _displacement_jump(declareProperty<RealVectorValue>("displacement_jump")),
    _traction_global(declareProperty<RealVectorValue>("traction_global")),
    _traction(declareProperty<RealVectorValue>("traction")),
    _traction_derivatives_global(declareProperty<RankTwoTensor>("traction_derivatives_global")),
    _traction_derivatives(declareProperty<RankTwoTensor>("traction_derivatives"))
{
  if (_ndisp > 3 || _ndisp < 1)
    mooseError("the CZM material requires 1, 2 or 3 displacement variables");

  if (getParam<bool>("use_displaced_mesh") == true)
    mooseError("This material cannot be used with use_displaced_mesh = true");

  // initializing the displacement vectors
  for (unsigned int i = 0; i < _ndisp; ++i)
  {
    _disp[i] = &coupledValue("displacements", i);
    _disp_neighbor[i] = &coupledNeighborValue("displacements", i);
  }
}

void
CZMMaterialBase::computeQpProperties()
{

  RealTensorValue RotationGlobalToLocal =
      RotationMatrix::rotVec1ToVec2(_normals[_qp], RealVectorValue(1, 0, 0));

  // computing the displacement jump
  for (unsigned int i = 0; i < _ndisp; i++)
    _displacement_jump_global[_qp](i) = (*_disp_neighbor[i])[_qp] - (*_disp[i])[_qp];
  for (unsigned int i = _ndisp; i < 3; i++)
    _displacement_jump_global[_qp](i) = 0;

  // rotate the displacement jump to the local coordiante system
  _displacement_jump[_qp] = RotationGlobalToLocal * _displacement_jump_global[_qp];

  // compute local traction
  _traction[_qp] = computeTraction();

  // compute local traction derivatives wrt the displacement jump
  _traction_derivatives[_qp] = computeTractionDerivatives();

  // rotate local traction and derivatives to the global coordinate system
  _traction_global[_qp] = RotationGlobalToLocal.transpose() * _traction[_qp];
  _traction_derivatives_global[_qp] = _traction_derivatives[_qp];
  _traction_derivatives_global[_qp].rotate(RotationGlobalToLocal.transpose());
}
