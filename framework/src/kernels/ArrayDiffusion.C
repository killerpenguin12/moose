//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ArrayDiffusion.h"

registerMooseObject("MooseApp", ArrayDiffusion);

defineLegacyParams(ArrayDiffusion);

InputParameters
ArrayDiffusion::validParams()
{
  InputParameters params = ArrayKernel::validParams();
  params.addParam<MaterialPropertyName>("diffusion_coefficient",
                                        "The name of the diffusivity, "
                                        "can be scalar, vector, or matrix.");
  params.addClassDescription(
      "The array Laplacian operator ($-\\nabla \\cdot \\nabla u$), with the weak "
      "form of $(\\nabla \\phi_i, \\nabla u_h)$.");
  return params;
}

ArrayDiffusion::ArrayDiffusion(const InputParameters & parameters)
  : ArrayKernel(parameters),
    _d(hasMaterialProperty<Real>("diffusion_coefficient")
           ? &getMaterialProperty<Real>("diffusion_coefficient")
           : nullptr),
    _d_array(hasMaterialProperty<RealEigenVector>("diffusion_coefficient")
                 ? &getMaterialProperty<RealEigenVector>("diffusion_coefficient")
                 : nullptr),
    _d_2d_array(hasMaterialProperty<RealEigenMatrix>("diffusion_coefficient")
                    ? &getMaterialProperty<RealEigenMatrix>("diffusion_coefficient")
                    : nullptr)
{
  if (!_d && !_d_array && !_d_2d_array)
  {
    MaterialPropertyName mat = getParam<MaterialPropertyName>("diffusion_coefficient");
    mooseError("Property " + mat + " is of unsupported type for ArrayDiffusion");
  }
}

RealEigenVector
ArrayDiffusion::computeQpResidual()
{
  if (_d)
    return _grad_u[_qp] * _array_grad_test[_i][_qp] * (*_d)[_qp];

  else if (_d_array)
  {
    mooseAssert((*_d_array)[_qp].size() == _var.count(),
                "diffusion_coefficient size is inconsistent with the number of components of array "
                "variable");
    return (*_d_array)[_qp].cwiseProduct(_grad_u[_qp] * _array_grad_test[_i][_qp]);
  }

  else
  {
    mooseAssert((*_d_2d_array)[_qp].cols() == _var.count(),
                "diffusion_coefficient size is inconsistent with the number of components of array "
                "variable");
    mooseAssert((*_d_2d_array)[_qp].rows() == _var.count(),
                "diffusion_coefficient size is inconsistent with the number of components of array "
                "variable");
    return (*_d_2d_array)[_qp] * (_grad_u[_qp] * _array_grad_test[_i][_qp]);
  }
}

RealEigenVector
ArrayDiffusion::computeQpJacobian()
{
  if (_d)
    return RealEigenVector::Constant(_var.count(),
                                     _grad_phi[_j][_qp] * _grad_test[_i][_qp] * (*_d)[_qp]);
  else if (_d_array)
    return _grad_phi[_j][_qp] * _grad_test[_i][_qp] * (*_d_array)[_qp];
  else
    return _grad_phi[_j][_qp] * _grad_test[_i][_qp] * (*_d_2d_array)[_qp].diagonal();
}

RealEigenMatrix
ArrayDiffusion::computeQpOffDiagJacobian(MooseVariableFEBase & jvar)
{
  if (jvar.number() == _var.number() && _d_2d_array)
    return _grad_phi[_j][_qp] * _grad_test[_i][_qp] * (*_d_2d_array)[_qp];
  else
    return ArrayKernel::computeQpOffDiagJacobian(jvar);
}
