//* This file is part of the MOOSE framework
//* https://www.mooseframework.org
//*
//* All rights reserved, see COPYRIGHT for full restrictions
//* https://github.com/idaholab/moose/blob/master/COPYRIGHT
//*
//* Licensed under LGPL 2.1, please see LICENSE for details
//* https://www.gnu.org/licenses/lgpl-2.1.html

#include "SuwaMobility.h"

registerMooseObject("PhaseFieldApp", SuwaMobility);

template <>
InputParameters
validParams<SuwaMobility>()
{
  InputParameters params = validParams<SuwaMobilitybase>();
  params.addRequiredParam<Real>("GBenergy", "Grain boundary energy in J/m^2");
  return params;
}

SuwaMobility::SuwaMobility(const InputParameters & parameters)
  : SuwaMobilitybase(parameters), _GBEnergy(getParam<Real>("GBenergy"))
{
}

void
SuwaMobility::computeQpProperties()
{
  // eV/nm^2
  _sigma[_qp] = _GBEnergy * _JtoeV * (_length_scale * _length_scale);

  SuwaMobilitybase::computeQpProperties();
}
