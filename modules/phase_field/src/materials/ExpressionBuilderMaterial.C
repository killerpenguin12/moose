//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "ExpressionBuilderMaterial.h"

template <>
InputParameters
validParams<ExpressionBuilderMaterial<> >()
{
  InputParameters params = validParams<Material>();
  params.addClassDescription("Material to provide a function (such as a free energy)");
  params.addParam<std::string>(
      "f_name",
      "F",
      "Base name of the free energy function (used to name the material properties)");
  params.addParam<bool>("third_derivatives", true, "Should we take third derivatives?");
  return params;
}

template <typename T>
ExpressionBuilderMaterial<T>::ExpressionBuilderMaterial(const InputParameters & parameters)
  : DerivativeMaterialInterface<Material>(parameters),
    _F_name(getParam<std::string>("f_name")),
    _third_derivatives(getParam<bool>("third_derivatives"))
{}

template <typename T>
void
ExpressionBuilderMaterial<T>::declareEBMaterial(EBTerm & term)
{
  //MaterialProperty<T> * temp;
  EBTerm temp;
  EBTerm temp2;
  EBTerm temp3;

  std::map<std::string, std::shared_ptr<EBTermNode> >::iterator it;
  std::map<std::string, std::shared_ptr<EBTermNode> >::iterator it2;
  std::map<std::string, std::shared_ptr<EBTermNode> >::iterator it3;
  term.simplify();
  _mat_props.push_back(std::pair<MaterialProperty<T> *, EBTerm>(&declareProperty<T>(_F_name),term));
  for(it = _coup_vars.begin(); it != _coup_vars.end(); ++it)
  {
    derivative_of = it->second;
    temp = term.takeDerivative();
    if(temp.simplify())
    {
      _mat_props.push_back(std::pair<MaterialProperty<T> *, EBTerm>(&declarePropertyDerivative<T>(_F_name, it->first),temp));
      for(it2 = _coup_vars.begin(); it2 != _coup_vars.end(); ++it2)
      {
        derivative_of = it2->second;
        temp2 = temp.takeDerivative();
        if(temp2.simplify())
        {
          _mat_props.push_back(std::pair<MaterialProperty<T> *, EBTerm>(&declarePropertyDerivative<T>(_F_name, it->first, it2->first),temp2));
          if(_third_derivatives)
          {
            for(it3 = _coup_vars.begin(); it3 != _coup_vars.end(); ++it3)
            {
              std::cout << "this happens" << std::endl;
              derivative_of = it3->second;
              temp3 = temp2.takeDerivative();
              if(temp3.simplify())
                _mat_props.push_back(std::pair<MaterialProperty<T> *, EBTerm>(&declarePropertyDerivative<T>(_F_name, it->first, it2->first, it3->first),temp3));
            }
          }
        }
      }
    }
  }
}

template <typename T>
ExpressionBuilder::EBTerm
ExpressionBuilderMaterial<T>::setEBMaterialPropertyReal(const std::string & name)
{
  FunctionMaterialPropertyDescriptor<Real> * mpd = new FunctionMaterialPropertyDescriptor<Real>(name, this);
  mpd->value();
  ExpressionBuilder::EBTerm new_term(std::make_shared<EBMatPropNode<Real> >(name, mpd));
  return new_term;
}

template <typename T>
ExpressionBuilder::EBVector
ExpressionBuilderMaterial<T>::setEBMaterialPropertyRealVectorValue(const std::string & name)
{
  FunctionMaterialPropertyDescriptor<RealVectorValue> * mpd = new FunctionMaterialPropertyDescriptor<RealVectorValue>(name, this);
  mpd->value();
  ExpressionBuilder::EBVector new_vec;
  new_vec[0] = EBTerm(std::make_shared<EBMatPropNode<RealVectorValue> >(name, mpd, 0));
  new_vec[1] = EBTerm(std::make_shared<EBMatPropNode<RealVectorValue> >(name, mpd, 1));
  new_vec[2] = EBTerm(std::make_shared<EBMatPropNode<RealVectorValue> >(name, mpd, 2));
  return new_vec;
}

template <typename T>
ExpressionBuilder::EBMatrix
ExpressionBuilderMaterial<T>::setEBMaterialPropertyRankTwoTensor(const std::string & name)
{
  FunctionMaterialPropertyDescriptor<RankTwoTensor> * mpd = new FunctionMaterialPropertyDescriptor<RankTwoTensor>(name, this);
  mpd->value();
  ExpressionBuilder::EBMatrix new_mat(3,3);
  new_mat[0][0] = EBTerm(std::make_shared<EBMatPropNode<RankTwoTensor> >(name, mpd, 0, 0));
  new_mat[0][1] = EBTerm(std::make_shared<EBMatPropNode<RankTwoTensor> >(name, mpd, 0, 1));
  new_mat[0][2] = EBTerm(std::make_shared<EBMatPropNode<RankTwoTensor> >(name, mpd, 0, 2));
  new_mat[1][0] = EBTerm(std::make_shared<EBMatPropNode<RankTwoTensor> >(name, mpd, 1, 0));
  new_mat[1][1] = EBTerm(std::make_shared<EBMatPropNode<RankTwoTensor> >(name, mpd, 1, 1));
  new_mat[1][2] = EBTerm(std::make_shared<EBMatPropNode<RankTwoTensor> >(name, mpd, 1, 2));
  new_mat[2][0] = EBTerm(std::make_shared<EBMatPropNode<RankTwoTensor> >(name, mpd, 2, 0));
  new_mat[2][1] = EBTerm(std::make_shared<EBMatPropNode<RankTwoTensor> >(name, mpd, 2, 1));
  new_mat[2][2] = EBTerm(std::make_shared<EBMatPropNode<RankTwoTensor> >(name, mpd, 2, 2));
  return new_mat;
}

template <typename T>
ExpressionBuilder::EBTerm
ExpressionBuilderMaterial<T>::setCoupledVariable(const std::string & name, unsigned int comp)
{
  ExpressionBuilder::EBTerm new_term;
  const VariableValue * value = &coupledValue(name, comp);
  std::map<std::string, std::shared_ptr<EBTermNode> >::iterator it = _coup_vars.find(name + std::to_string(comp));
  if(it == _coup_vars.end())
  {
    std::shared_ptr<EBTermNode> new_node = std::make_shared<EBCoupledVarNode<VariableValue> >(name + std::to_string(comp), value);
    new_term = EBTerm(new_node);
    _coup_vars[name + std::to_string(comp)] = new_node;
  }
  else
    new_term = EBTerm(it->second);
  return new_term;
}

template <typename T>
std::vector<ExpressionBuilder::EBTerm>
ExpressionBuilderMaterial<T>::setCoupledVarVector(const std::string & name, unsigned int _op_num)
{
  std::vector<ExpressionBuilder::EBTerm> new_terms;
  for(unsigned int i = 0; i < _op_num; ++i)
  {
    new_terms.push_back(setCoupledVariable(name, i));
  }
  return new_terms;
}

template <typename T>
ExpressionBuilder::EBVector
ExpressionBuilderMaterial<T>::setCoupledGradient(const std::string & name, unsigned int comp)
{
  ExpressionBuilder::EBVector new_vec;
  const VariableGradient * value = &coupledGradient(name, comp);
  std::map<std::string, std::shared_ptr<EBTermNode> >::iterator it;
  std::shared_ptr<EBTermNode> new_node;
  for(unsigned int i = 0; i < 3; ++i)
  {
    it = _coup_vars.find(name + std::to_string(comp) + std::to_string(i));
    if(it == _coup_vars.end())
    {
      new_node = std::make_shared<EBCoupledVarNode<VariableGradient> >(name + std::to_string(comp) + std::to_string(i), value, i);
      new_vec[i] = EBTerm(new_node);
      _coup_vars[name + std::to_string(comp) + std::to_string(i)] = new_node;
    }
    else
      new_vec[i] = EBTerm(it->second);
  }
  return new_vec;
}

template <typename T>
std::vector<ExpressionBuilder::EBVector>
ExpressionBuilderMaterial<T>::setCoupledGradVector(const std::string & name, unsigned int _op_num)
{
  std::vector<ExpressionBuilder::EBVector> new_vecs;
  for(unsigned int i = 0; i < _op_num; ++i)
  {
    new_vecs.push_back(setCoupledGradient(name, i));
  }
  return new_vecs;
}

template <typename T>
ExpressionBuilder::EBMatrix
ExpressionBuilderMaterial<T>::setCoupledSecond(const std::string & name, unsigned int comp)
{
  ExpressionBuilder::EBMatrix new_mat(3,3);
  const VariableSecond * value = &coupledSecond(name, comp);
  std::map<std::string, std::shared_ptr<EBTermNode> >::iterator it;
  std::shared_ptr<EBTermNode> new_node;
  for(unsigned int i = 0; i < 3; ++i)
    for(unsigned int j = 0; j < 3; ++j)
    {
      it = _coup_vars.find(name + std::to_string(comp) + std::to_string(i) + std::to_string(j));
      if(it == _coup_vars.end())
      {
        new_node = std::make_shared<EBCoupledVarNode<VariableSecond> >(name + std::to_string(comp) + std::to_string(i) + std::to_string(j), value, i, j);
        new_mat[i][j] = EBTerm(new_node);
        _coup_vars[name + std::to_string(comp) + std::to_string(i) + std::to_string(j)] = new_node;
      }
      else
        new_mat[i][j] = EBTerm(it->second);
    }
  return new_mat;
}

template <typename T>
std::vector<ExpressionBuilder::EBMatrix>
ExpressionBuilderMaterial<T>::setCoupledSecVector(const std::string & name, unsigned int _op_num)
{
  std::vector<ExpressionBuilder::EBMatrix> new_mats;
  for(unsigned int i = 0; i < _op_num; ++i)
  {
    new_mats.push_back(setCoupledSecond(name, i));
  }
  return new_mats;
}

template <typename T>
void
ExpressionBuilderMaterial<T>::computeQpProperties()
{
  current_quad_point = _qp;
  typename std::vector<std::pair<MaterialProperty<T> *, EBTerm> >::iterator it;
  for(it = _mat_props.begin(); it != _mat_props.end(); ++it)
  {
    it->second.cleanUp();
    (*it->first)[_qp] = it->second.evaluate();
  }
}
