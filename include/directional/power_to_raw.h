// This file is part of libdirectional, a library for directional field processing.
// Copyright (C) 2018 Amir Vaxman <avaxman@gmail.com>
//
// This Source Code Form is subject to the terms of the Mozilla Public License
// v. 2.0. If a copy of the MPL was not distributed with this file, You can
// obtain one at http://mozilla.org/MPL/2.0/.
#ifndef DIRECTIONAL_POWER_TO_RAW_H
#define DIRECTIONAL_POWER_TO_RAW_H
#include <directional/rotation_to_representative.h>
#include <directional/representative_to_raw.h>
#include <directional/power_to_representative.h>
#include <igl/local_basis.h>


namespace directional
{
  // Converts the power complex representation to raw representation.
  // Inputs:
  //  B1, B2:
  //  B3: #F normals for each face/B3 from igl::local_base.
  //  N: the degree of the field.
  //  powerField: Representation of the field as complex double
  // Outputs:
  //  rawField: #F by 3*N matrix with all N explicit vectors of each directional in the order X,Y,Z,X,Y,Z, ...
  IGL_INLINE void power_to_raw(const Eigen::MatrixXd& B1,
                               const Eigen::MatrixXd& B2,
                               const Eigen::MatrixXd& B3,
                               const Eigen::MatrixXcd& powerField,
                               int N,
                               Eigen::MatrixXd& rawField,
                               bool normalize=false)
  {
    Eigen::MatrixXd representative;
    power_to_representative(B1, B2, powerField, N, representative);
    if (normalize)
      representative.rowwise().normalize();
    representative_to_raw(B3, representative, N, rawField);
  }
  
  // version without auxiliary data
  IGL_INLINE void power_to_raw(const Eigen::MatrixXd& V,
                               const Eigen::MatrixXi& F,
                               const Eigen::MatrixXcd& powerField,
                               int N,
                               Eigen::MatrixXd& rawField,
                               bool normalize=false)
  {
    Eigen::MatrixXd B1, B2, B3;
    igl::local_basis(V, F, B1, B2, B3);
    power_to_raw(B1, B2, B3, powerField, N, rawField,normalize);
  }
}

#endif
