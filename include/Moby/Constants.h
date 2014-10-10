/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_CONSTANTS_H
#define _MOBY_CONSTANTS_H

#include <limits>
#include <cmath>
#include <Ravelin/Vector3d.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Pose3d.h>

namespace Moby {

// constants
const boost::shared_ptr<const Ravelin::Pose3d> GLOBAL;
const double NEAR_ZERO = std::sqrt(std::numeric_limits<double>::epsilon());
const Ravelin::Matrix3d ZEROS_3x3 = Ravelin::Matrix3d::zero();
const Ravelin::Matrix3d IDENTITY_3x3 = Ravelin::Matrix3d::identity();
const Ravelin::VectorNd EMPTY_VEC(0);

// debugging bits
const unsigned LOG_SIMULATOR = 1;
const unsigned LOG_CONSTRAINT = 2;
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
