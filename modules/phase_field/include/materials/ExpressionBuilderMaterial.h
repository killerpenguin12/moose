//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#pragma once

#include "Material.h"
#include "DerivativeMaterialInterface.h"
#include "ExpressionBuilder.h"

// Forward Declarations
template <typename T = Real>
class ExpressionBuilderMaterial;

template <>
InputParameters validParams<ExpressionBuilderMaterial<> >();

/**
 * %Material base class central for all Materials that provide a Function as a
 * material property value.
 */
template <typename T>
class ExpressionBuilderMaterial : public DerivativeMaterialInterface<Material>, public ExpressionBuilder
{
public:
  ExpressionBuilderMaterial(const InputParameters & parameters);

  class CoupledVariableDescriptor;

  InputParameters validParams();

  void declareEBMaterial(EBTerm & term);

  EBTerm setEBMaterialPropertyReal(const std::string & name);
  EBVector setEBMaterialPropertyRealVectorValue(const std::string & name);
  EBMatrix setEBMaterialPropertyRankTwoTensor(const std::string & name);

  EBTerm setCoupledVariable(const std::string & name, unsigned int comp = 0);
  std::vector<EBTerm> setCoupledVarVector(const std::string & name, unsigned int _op_num);
  EBVector setCoupledGradient(const std::string & name, unsigned int comp = 0);
  std::vector<EBVector> setCoupledGradVector(const std::string & name, unsigned int _op_num);
  EBMatrix setCoupledSecond(const std::string & name, unsigned int comp = 0);
  std::vector<EBMatrix> setCoupledSecVector(const std::string & name, unsigned int _op_num);

protected:
  /**
   * ExpressionBuilderMaterial keeps an internal list of all the variables the derivatives are taken
   * w.r.t.
   * We provide the MOOSE variable names in _arg_names, the libMesh variable numbers in
   * _arg_numbers, and the
   * input file parameter names in _arg_param_names. All are indexed by the argument index.
   * This method returns the argument index for a given the libMesh variable number.
   *
   * This mapping is necessary for internal classes which maintain lists of derivatives indexed by
   * argument index
   * and need to pull from those lists from the computeDF, computeD2F, and computeD3F methods, which
   * receive
   * libMesh variable numbers as parameters.
   */

  /**
   * Name of the function value material property and used as a base name to
   * concatenate the material property names for the derivatives.
   */
  std::string _F_name;

  /// String vector of all argument names.
  std::map<std::string, std::shared_ptr<EBTermNode> > _coup_vars;

  /// Calculate (and allocate memory for) the third derivatives of the free energy.
  bool _third_derivatives;

  std::vector<std::pair<MaterialProperty<T> *, EBTerm> > _mat_props;

  virtual void computeQpProperties();
};
