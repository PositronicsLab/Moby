/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_PLANE_H_
#define _MOBY_PLANE_H_

#include <Moby/Constants.h>
#include <Moby/Triangle.h>
#include <Moby/Vector3.h>

namespace Moby {

/// Defines a plane using the equation <n, x> = d
class Plane
{
  public:
    Plane() { _normal = ZEROS_3; offset = 0; }
    Plane(Real nx, Real ny, Real nz, Real d) { _normal = Vector3(nx, ny, nz); _normal.normalize(); offset = d; }
    Plane(const Triangle& t) { _normal = t.calc_normal(); offset = t.calc_offset(_normal); }
    Plane(const Vector3& n, Real d) { _normal = Vector3::normalize(n); offset = d; }
    Plane(const Plane& p) { operator=(p); }
    void operator=(const Plane& p) { _normal = p._normal; offset = p.offset; }
    Real calc_signed_distance(const Vector3& p) const { return _normal.dot(p) - offset; }
    bool on_plane(const Vector3& p) { return std::fabs(calc_signed_distance(p)) < NEAR_ZERO; }
    Plane operator-() const { return Plane(-_normal, -offset); }
    bool operator==(const Plane& p) const { Real dot = _normal.dot(p.get_normal()); return (std::fabs(dot - 1.0) < NEAR_ZERO && std::fabs(offset - p.offset) < NEAR_ZERO); }
    bool operator<(const Plane& p) const;

    /// Gets the outward pointing normal to the plane
    const Vector3& get_normal() const { return _normal; }

    /// Sets the outward pointing normal to the plane
    void set_normal(const Vector3& n) { _normal = Vector3::normalize(n); }

    /// The plane offset such that the plane equation is given by <n, x> = offset
    Real offset;

  private:
    Vector3 _normal;
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

