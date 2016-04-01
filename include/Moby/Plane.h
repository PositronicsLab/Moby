/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_PLANE_H_
#define _MOBY_PLANE_H_

#include <Moby/Constants.h>
#include <Moby/Triangle.h>
#include <Ravelin/Vector3d.h>

namespace Moby {

/// Defines a plane using the equation <n, x> = d
class Plane
{
  public:
    Plane() { _normal.set_zero(); offset = 0; }
    Plane(double nx, double ny, double nz, double d) { _normal = Ravelin::Vector3d(nx, ny, nz); _normal.normalize(); offset = d; }
    Plane(const Triangle& t) { _normal = t.calc_normal(); offset = t.calc_offset(_normal); }
    Plane(const Ravelin::Vector3d& n, double d) { _normal = Ravelin::Vector3d::normalize(n); offset = d; }
    Plane(const Ravelin::Vector3d& normal, const Point3d& point) { _normal = Ravelin::Vector3d::normalize(normal); offset = normal.dot(point); }
    Plane(const Plane& p) { operator=(p); }
    void operator=(const Plane& p) { _normal = p._normal; offset = p.offset; }
    double calc_signed_distance(const Point3d& p) const { return _normal.dot(p) - offset; }
    bool on_plane(const Point3d& p) { return std::fabs(calc_signed_distance(p)) < NEAR_ZERO * std::max((double) 1.0, std::max(p.norm_inf(), std::fabs(offset))); }
    bool operator==(const Plane& p) const { double dot = _normal.dot(p.get_normal()); return (std::fabs(dot - 1.0) < NEAR_ZERO && std::fabs(offset - p.offset) < NEAR_ZERO); }
    bool operator<(const Plane& p) const;
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return _normal.pose; }
    Point3d project(const Point3d& p) const;
    Ravelin::Origin2d to_2D(const Point3d& p) const;
    Plane transform(const Ravelin::Transform3d& T) const;

    /// Assuming this plane represents a half-space, negates the half-space
    Plane operator-() const { return Plane(-_normal, -offset); }

    /// Gets the outward pointing normal to the plane
    const Ravelin::Vector3d& get_normal() const { return _normal; }

    /// Sets the outward pointing normal to the plane
    void set_normal(const Ravelin::Vector3d& n) { _normal = Ravelin::Vector3d::normalize(n); }

    /// The plane offset such that the plane equation is given by <n, x> = offset
    double offset;

  private:
    Ravelin::Vector3d _normal;
}; // end class

/// Transforms a plane
inline Plane Plane::transform(const Ravelin::Transform3d& T) const
{
  // setup the new plane
  Plane p;

  // get the new normal
  p._normal = T.transform_vector(_normal);

  // compute the new offset using a point on the old plane
  // NOTE: n' * (n*offset) - offset = 0
  Point3d plane_point = _normal * offset;

  // transform that point
  Point3d new_plane_point = T.transform_point(plane_point);

  // now compute the offset
  p.offset = p._normal.dot(new_plane_point);

  return p;
}

/// Transforms a point to 2D
inline Ravelin::Origin2d Plane::to_2D(const Point3d& p) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup the projection matrix
  Ravelin::Vector3d v1, v2;
  Ravelin::Vector3d::determine_orthonormal_basis(_normal, v1, v2);
  Ravelin::Matrix3d R;
  R.set_row(X, Ravelin::Origin3d(v1));
  R.set_row(Y, Ravelin::Origin3d(v2));
  R.set_row(Z, Ravelin::Origin3d(_normal));

  // multiply
  Ravelin::Origin3d result = R * Ravelin::Origin3d(p);

  return Ravelin::Origin2d(result[X], result[Y]);
}

/// Projects a point onto the plane
inline Point3d Plane::project(const Point3d& p) const
{
  double d = _normal.dot(p);
  return p - _normal*d;
}

inline bool Plane::operator<(const Plane& p) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  if (_normal[X] < p.get_normal()[X])
    return true;
  else if (_normal[X] == p.get_normal()[X])
  {
    if (_normal[Y] < p.get_normal()[Y])
      return true;
    else if (_normal[Y] == p.get_normal()[Y])
    {
      if (_normal[Z] < p.get_normal()[Z])
        return true;
      else if (_normal[Z] == p.get_normal()[Z])
        return offset < p.offset;
    }
  }

  return false;
} 

inline std::ostream& operator<<(std::ostream& out, const Plane& p)
{
  out << "normal: " << p.get_normal() << " offset: " << p.offset;
  return out;
}

} // end namespace

#endif

