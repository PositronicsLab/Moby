/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_CONSTANTS_H
#define _MOBY_CONSTANTS_H

#include <limits>
#include <cmath>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/Vector3.h>
#include <Moby/SVector6.h>

namespace Moby {

// constants
const Real NEAR_ZERO = std::sqrt(std::numeric_limits<Real>::epsilon());
const Vector3 ZEROS_3 = Vector3::zero();
const Matrix3 ZEROS_3x3 = Matrix3::zero();
const Matrix3 IDENTITY_3x3 = Matrix3::identity();
const Matrix4 IDENTITY_4x4 = Matrix4::identity();
const SVector6 ZEROS_6 = SVector6::zero();
const VectorN EMPTY_VEC(0);

// debugging bits
const unsigned LOG_SIMULATOR = 1;
const unsigned LOG_CONTACT = 2;
const unsigned LOG_DYNAMICS = 4;
const unsigned LOG_BV = 8;
const unsigned LOG_ADF = 16;
const unsigned LOG_COLDET = 32;
const unsigned LOG_COMPGEOM = 64;
const unsigned LOG_LINALG = 128;
const unsigned LOG_OPT = 256;
const unsigned LOG_DEFORM = 512;

}

#endif
