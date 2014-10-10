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
    Plane operator-() const { return Plane(-_normal, -offset); }
    bool operator==(const Plane& p) const { double dot = _normal.dot(p.get_normal()); return (std::fabs(dot - 1.0) < NEAR_ZERO && std::fabs(offset - p.offset) < NEAR_ZERO); }
    bool operator<(const Plane& p) const;
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return _normal.pose; }

    /// Gets the outward pointing normal to the plane
    const Ravelin::Vector3d& get_normal() const { return _normal; }

    /// Sets the outward pointing normal to the plane
    void set_normal(const Ravelin::Vector3d& n) { _normal = Ravelin::Vector3d::normalize(n); }

    /// The plane offset such that the plane equation is given by <n, x> = offset
    double offset;

  private:
    Ravelin::Vector3d _normal;
}; // end class

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

