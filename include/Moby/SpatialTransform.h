/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SPATIAL_TRANSFORM_H
#define _SPATIAL_TRANSFORM_H

#include <Moby/SVector6.h>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/SpatialRBInertia.h>
#include <Moby/SpatialABInertia.h>

namespace Moby {

/// A spatial transformation matrix, equivalent to a 6x6 matrix 
class SpatialTransform
{
  public:
    SpatialTransform();
    SpatialTransform(const SpatialTransform& source) { operator=(source); }
    SpatialTransform(const Matrix4& source, const Matrix4& target);
    SpatialTransform(const Matrix3& sourceR, const Vector3& sourceX, const Matrix3& targetR, const Vector3& targetX); 
    SpatialTransform(const Matrix3& sourceR, const Vector3& sourceX, const Matrix4& target); 
    SpatialTransform(const Matrix4& source, const Matrix3& targetR, const Vector3& targetX); 
    void set_identity();
    static SpatialTransform identity() { SpatialTransform m; m.set_identity(); return m; }
/*
    MatrixN& mult(const SMatrix6N* m, MatrixN& result) const;
*/
    SpatialTransform& operator=(const SpatialTransform& source);
    SpatialABInertia transform(const SpatialABInertia& m) const;
    SpatialRBInertia transform(const SpatialRBInertia& m) const;
    SpatialRBInertia transform(Real m, const Matrix3& J) const;
    SpatialTransform concatenate(const SpatialTransform& m) const;
    SVector6 transform(const SVector6& v) const;
    SMatrix6N transform(const SMatrix6N& m) const;
    SMatrix6N& transform(const SMatrix6N& m, SMatrix6N& result) const;
    static SpatialTransform to_global(const Matrix4& T);
    static SpatialTransform from_global(const Matrix4& T);

    /// The rotation component
    Matrix3 E;

    /// The translation component
    Vector3 r;
}; // end class

std::ostream& operator<<(std::ostream& out, const SpatialTransform& m);

} // end namespace

#endif

