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

// Forward Declarations
class SuwaMobilitybase;

template <>
InputParameters validParams<SuwaMobilitybase>();

class SuwaMobilitybase : public DerivativeMaterialInterface<Material>
{
public:
  SuwaMobilitybase(const InputParameters & parameters);

protected:
  virtual void computeQpProperties();

  Real _f0s;          //GB energy constant
  Real _wGB;          //diffuse GB width in proper length scale
  Real _length_scale;   //for converting to proper scale
  Real _time_scale;    //for converting to proper time scale
  Real _GBmob0;        //have to provide this and Q or _GBMobility, this is for
                      //temperature dependant calcs
  Real _Q;            //activation energy in eV
  Real _GBMobility;   //provided mobility for GB
  Real _molar_vol;
  Real _Mis;

  const VariableValue & _T;

  MaterialProperty<Real> & _sigma;
  MaterialProperty<Real> & _M_GB;
  MaterialProperty<Real> & _kappa;
  MaterialProperty<Real> & _gamma;
  MaterialProperty<Real> & _L;
  MaterialProperty<Real> * _dLdT;
  MaterialProperty<Real> & _l_GB;
  MaterialProperty<Real> & _mu;
  MaterialProperty<Real> & _entropy_diff;
  MaterialProperty<Real> & _molar_volume;
  MaterialProperty<Real> & _act_wGB;

  const Real _thetaRad;
  const Real _Lo;
  const Real _kb;
  const Real _JtoeV;
};
