/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/Tetrahedron.h>
#include <Moby/Triangle.h>
#include <Moby/CompGeom.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/Log.h>
#include <Moby/GJK.h>

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

std::ostream& GJK::SVertex::output(std::ostream& out) const
{
  out << " " << v << std::endl;
  out << "   A vertex: " << vA << std::endl;
  out << "   B vertex: " << vB << std::endl;

  return out;
}

std::ostream& GJK::Simplex::output(std::ostream& out) const
{
  out << " type: ";
  switch (_type)
  {
    case ePoint:       out << "point" << std::endl; break;
    case eSegment:     out << "segment" << std::endl; break;
    case eTriangle:    out << "triangle" << std::endl; break;
    case eTetra:       out << "tetrahedron" << std::endl; break;
    default:           assert(false);
  }

  for (unsigned i=0; i< num_vertices(); i++)
    get_vertex(i).output(out) << std::endl;

  return out;
}

/// Gets the desired vertex
const GJK::SVertex& GJK::Simplex::get_vertex(unsigned i) const
{
  switch (i)
  {
    case 0: return _v1;
    case 1: return _v2;
    case 2: return _v3;
    case 3: return _v4;
    default: assert(false); return _v1;
  }
}

/// Gets the number of vertices in the simplex
unsigned GJK::Simplex::num_vertices() const
{
  switch (_type)
  {
    case ePoint:    return 1;
    case eSegment:  return 2;
    case eTriangle: return 3;
    case eTetra:    return 4;
    default: assert(false); return 0;
  }
}

/// Finds the closest point on the simplex (debugging function)
Point3d GJK::Simplex::find_closest()
{
  if (_type == ePoint)
    return _v1.v;
  else if (_type == eSegment)
  {
    double t;
    CompGeom::calc_dist(LineSeg3(_v1.v, _v2.v), Point3d(0,0,0,_v1.v.pose), t);
    if (t < NEAR_ZERO) 
      return _v1.v;
    else
      return _v1.v*t + _v2.v*(1.0-t);
  }
  else if (_type == eTriangle)
  {
    Point3d cp;
    Triangle tri(_v1.v, _v2.v, _v3.v);
    Triangle::calc_sq_dist(tri, Point3d(0,0,0,_v1.v.pose), cp);
    return cp;
  }
  else 
  {
    Point3d cp;
    Tetrahedron tetra(_v1.v, _v2.v, _v3.v, _v4.v);
    tetra.calc_signed_dist(Point3d(0,0,0,_v1.v.pose), cp);
    double u, v, w;
    tetra.determine_barycentric_coords(cp, u, v, w);
    return cp;
  }
}

/// Finds the closest point to the origin and simplifies the simplex (if possible)
Point3d GJK::Simplex::find_closest_and_simplify()
{
  if (_type == ePoint)
    return _v1.v;
  else if (_type == eSegment)
  {
    double t;
    CompGeom::calc_dist(LineSeg3(_v1.v, _v2.v), Point3d(0,0,0,_v1.v.pose), t);
    if (t < NEAR_ZERO) 
    {
      // new point is the closest; remove the old point
      _v1 = _v2;
      _type = ePoint;
      assert((find_closest()-_v1.v).norm() < NEAR_ZERO);
      return _v1.v;
    }
    else
    {
      // setup the closest point
      return _v1.v*t + _v2.v*(1.0-t);
    }
  }
  else if (_type == eTriangle)
  {
    Point3d cp;
    Triangle tri(_v1.v, _v2.v, _v3.v);
    Triangle::calc_sq_dist(tri, Point3d(0,0,0,_v1.v.pose), cp);
    double s, t;
    tri.determine_barycentric_coords(cp, s, t);
    assert(s >= -NEAR_ZERO && t >= -NEAR_ZERO && s+t <= 1.0+NEAR_ZERO);
    if (s < NEAR_ZERO && t < NEAR_ZERO)
    {
      // new simplex is the newest vertex
      _v1 = _v3;
      _type = ePoint;
      assert((find_closest()-_v1.v).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(t-1.0) < NEAR_ZERO)
    {
      // simplex is second point
      _v1 = _v2;
      _type = ePoint;
      assert((find_closest()-_v1.v).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(s-1.0) < NEAR_ZERO)
    {
      // simplex is first point (shouldn't be here?)
      _type = ePoint;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(s+t-1.0) < NEAR_ZERO)
    {
      // simplex is edge ab
      _type = eSegment;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (s < NEAR_ZERO)
    {
      // simplex is edge bc
      _type = eSegment;
      _v1 = _v2;
      _v2 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (t < NEAR_ZERO)
    {
      // simplex is edge ca
      _type = eSegment;
      _v2 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else
    {
      // simplex is entire triangle
      return cp;
    }
  }
  else if (_type == eTetra)
  {
    Point3d cp;
    Tetrahedron tetra(_v1.v, _v2.v, _v3.v, _v4.v);
    tetra.calc_signed_dist(Point3d(0,0,0,_v1.v.pose), cp);
    double u, v, w;
    tetra.determine_barycentric_coords(cp, u, v, w);
    assert(u >= -NEAR_ZERO && v >= -NEAR_ZERO && w >= -NEAR_ZERO && u+v+w <= 1.0+NEAR_ZERO);

    // handle case of individual vertices first
    if (u < NEAR_ZERO && v < NEAR_ZERO && w < NEAR_ZERO)
    {
      // simplex simplified to newest vertex
      _type = ePoint;
      _v1 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(u-1.0) < NEAR_ZERO)
    {
      _type = ePoint;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(v-1.0) < NEAR_ZERO)
    {
      _type = ePoint;
      _v1 = _v2;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(w-1.0) < NEAR_ZERO)
    {
      _type = ePoint;
      _v1 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(u+v-1.0) < NEAR_ZERO)
    {
      // edge ab
      _type = eSegment;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(v+w-1.0) < NEAR_ZERO)
    {
      // edge bc
      _type = eSegment;
      _v1 = _v2;
      _v2 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(u+w-1.0) < NEAR_ZERO)
    {
      // edge ac
      _type = eSegment;
      _v2 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (v < NEAR_ZERO && w < NEAR_ZERO)
    {
      // edge ad
      _type = eSegment;
      _v2 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (u < NEAR_ZERO && w < NEAR_ZERO)
    {
      // edge bd
      _type = eSegment;
      _v1 = _v2;
      _v2 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (u < NEAR_ZERO && v < NEAR_ZERO)
    {
      // edge cd
      _type = eSegment;
      _v1 = _v3;
      _v2 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (std::fabs(u+v+w-1.0) < NEAR_ZERO)
    {
      // triangle abc
      _type = eTriangle;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (w < NEAR_ZERO)
    {
      // triangle abd
      _type = eTriangle;
      _v3 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (u < NEAR_ZERO)
    {
      // triangle bcd
      _type = eTriangle;
      _v1 = _v2;
      _v2 = _v3;
      _v3 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else if (v < NEAR_ZERO)
    {
      // triangle cad
      _type = eTriangle;
      _v2 = _v1;
      _v1 = _v3;
      _v3 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      return cp;
    }
    else
    {
      // simplex is tetrahedron
      return cp;
    }
  }
  else
  {
    // should never be here...
    assert(false);
    return Point3d(0,0,0,_v1.v.pose);
  }
}

/// Adds a vertex to the simplex
void GJK::Simplex::add(const SVertex& v)
{
  if (_type == ePoint)
  {
    _type = eSegment;
    _v2 = v;
  }
  else if (_type == eSegment)
  {
    _type = eTriangle;
    _v3 = v;
    if (std::fabs(Triangle(_v1.v, _v2.v, _v3.v).calc_area()) < NEAR_ZERO)
      _type = eSegment;
  }
  else if (_type == eTriangle)
  {
    _type = eTetra;
    _v4 = v;
    double vol = Tetrahedron(_v1.v, _v2.v, _v3.v, _v4.v).calc_volume();
    if (std::fabs(vol) < NEAR_ZERO)
      _type = eTriangle;
    else if (vol < 0)
      std::swap(_v2.v, _v3.v);
  }
  else if (_type == eTetra)
    assert(false);
}

double GJK::do_gjk(CollisionGeometryPtr A, CollisionGeometryPtr B, Point3d& closestA, Point3d& closestB, unsigned max_iter)
{
  const double INF = std::numeric_limits<double>::max();

  // setup a random direction
  Point3d rdir((double) rand() / RAND_MAX * 2.0 - 1.0,(double) rand() / RAND_MAX * 2.0 - 1.0, (double) rand() / RAND_MAX * 2.0 - 1.0, GLOBAL);
  Point3d pA = A->get_supporting_point(-rdir);
  Point3d pB = B->get_supporting_point(rdir); 

  // setup the initial support (a point)
  SVertex p(pA, pB);
  Simplex S = p;
  if (LOGGING(LOG_COLDET))
  {
    std::ostringstream oss;
    S.output(oss); 
    FILE_LOG(LOG_COLDET) << "GJK::do_gjk() entered" << std::endl;
    FILE_LOG(LOG_COLDET) << " -- initial simplex: " << oss.str() << std::endl;
  }

  // setup the minimum dot
  double min_dot = std::numeric_limits<double>::max();
  double min_dist = std::numeric_limits<double>::max();

  // GJK loop
  for (unsigned i=0; i< max_iter; i++)
  {
    // find the closest point in the simplex to the origin
    Point3d p = S.find_closest_and_simplify();
    if (LOGGING(LOG_COLDET))
    {
      std::ostringstream oss;
      S.output(oss); 
      FILE_LOG(LOG_COLDET) << " -- closest point on simplex to origin: " << p << std::endl;
      FILE_LOG(LOG_COLDET) << " -- distance to origin: " << p.norm() << std::endl;
      FILE_LOG(LOG_COLDET) << " -- new simplex: " << oss.str() << std::endl;
    }

    // look and see whether the origin is contained in the simplex
    double pnorm = p.norm();
    if (pnorm < NEAR_ZERO)
    {
      // A and B are intersecting
      // determine the interpenetration distance
      double pen_dist = INF;
      const unsigned NV = S.num_vertices();
      for (unsigned i=0; i< NV; i++)
      {
        double dA = A->calc_signed_dist(S.get_vertex(i).vB);
        if (dA < 0.0) 
          pen_dist = std::min(pen_dist, -dA);
        double dB = B->calc_signed_dist(S.get_vertex(i).vA);
        if (dB < 0.0) 
          pen_dist = std::min(pen_dist, -dB);
      }

      return -pen_dist;
    }

    // get the new supporting points and determine the new vertex
    Point3d pA = A->get_supporting_point(-p);
    Point3d pB = B->get_supporting_point(p); 
    SVertex V(pA, pB);
    if (LOGGING(LOG_COLDET))
    {
      std::ostringstream oss;
      V.output(oss); 
      FILE_LOG(LOG_COLDET) << " -- new vertex: " << oss.str() << std::endl;
    }

    // get the minimum distance  
    min_dist = std::min(min_dist, pnorm);

    // look to see whether no intersection
    double vdotd = V.v.dot(-p);
    FILE_LOG(LOG_COLDET) << " -- <new vertex, direction> : " << vdotd << std::endl;
    if (false && vdotd < min_dot)
    {
      min_dot = vdotd;
      closestA = pA;
      closestB = pB;
      if (vdotd < 0.0)
        return min_dist;
    }
    else
    {
      // add the new vertex to the simplex
      S.add(V);
      if (LOGGING(LOG_COLDET))
      {
        std::ostringstream oss;
        S.output(oss); 
        FILE_LOG(LOG_COLDET) << "GJK::do_gjk() entered" << std::endl;
        FILE_LOG(LOG_COLDET) << "added new point to simplex, now: " << oss.str() << std::endl;
      }
    }
  }

  FILE_LOG(LOG_COLDET) << "GJK::do_gjk() exited" << std::endl;
  return min_dist;
}

