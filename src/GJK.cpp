/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
    Point3d closest;
    CompGeom::calc_dist(LineSeg3(_v1.v, _v2.v), Point3d(0,0,0,_v1.v.pose), t, closest);
    if (t < NEAR_ZERO) 
      return _v1.v;
    else
      return closest; 
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
  FILE_LOG(LOG_COLDET) << "Simplex::find_closest_and_simplify() entered" << std::endl;

  if (_type == ePoint)
    return _v1.v;
  else if (_type == eSegment)
  {
    FILE_LOG(LOG_COLDET) << " -- current simplex is segment" << std::endl;

    double t;
    Point3d closest;
    CompGeom::calc_dist(LineSeg3(_v1.v, _v2.v), Point3d(0,0,0,_v1.v.pose), t, closest);
    if (t >= 1.0-NEAR_ZERO) 
    {
      // new point is the closest; remove the old point
      _v1 = _v2;
      _type = ePoint;
      assert((find_closest()-_v1.v).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is v2 (removing v1)" << std::endl;
      return _v1.v;
    }
    else
    {
      FILE_LOG(LOG_COLDET) << " -- closest point is on line segment" << std::endl;

      // setup the closest point
      return closest; 
    }
  }
  else if (_type == eTriangle)
  {
    FILE_LOG(LOG_COLDET) << " -- current simplex is on triangle" << std::endl;
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
      FILE_LOG(LOG_COLDET) << " -- closest point is to newest vertex (v3); reverting to point simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(t-1.0) < NEAR_ZERO)
    {
      // simplex is second point
      _v1 = _v2;
      _type = ePoint;
      assert((find_closest()-_v1.v).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to vertex v2; reverting to point simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(s-1.0) < NEAR_ZERO)
    {
      // simplex is first point (shouldn't be here?)
      _type = ePoint;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to vertex v1 (shouldn't be here!); reverting to point simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(s+t-1.0) < NEAR_ZERO)
    {
      // simplex is edge ab
      _type = eSegment;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to edge (v1,v2); reverting to edge simplex" << std::endl;
      return cp;
    }
    else if (s < NEAR_ZERO)
    {
      // simplex is edge bc
      _type = eSegment;
      _v1 = _v2;
      _v2 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to edge (v2,v3); reverting to edge simplex" << std::endl;
      return cp;
    }
    else if (t < NEAR_ZERO)
    {
      // simplex is edge ca
      _type = eSegment;
      _v2 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to edge (v3,v1); reverting to edge simplex" << std::endl;
      return cp;
    }
    else
    {
      // simplex is entire triangle
      FILE_LOG(LOG_COLDET) << " -- closest point is interior of triangle (not reducing simplex)" << std::endl;
      return cp;
    }
  }
  else if (_type == eTetra)
  {
    FILE_LOG(LOG_COLDET) << " -- current simplex is tetrahedron" << std::endl;
    Point3d cp;
    Tetrahedron tetra(_v1.v, _v2.v, _v3.v, _v4.v);
    tetra.calc_signed_dist(Point3d(0,0,0,_v1.v.pose), cp);
    double u, v, w;
    tetra.determine_barycentric_coords(cp, u, v, w);
    FILE_LOG(LOG_COLDET) << " u: " << u << " v: " << v << " w: " << w << std::endl;
    assert(u >= -NEAR_ZERO && v >= -NEAR_ZERO && w >= -NEAR_ZERO && u+v+w <= 1.0+NEAR_ZERO);

    // handle case of individual vertices first
    if (u < NEAR_ZERO && v < NEAR_ZERO && w < NEAR_ZERO)
    {
      // simplex simplified to newest vertex
      _type = ePoint;
      _v1 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to newest vertex (v4); reverting to point simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(u-1.0) < NEAR_ZERO)
    {
      _type = ePoint;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to vertex (v1)- shouldn't be here? reverting to point simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(v-1.0) < NEAR_ZERO)
    {
      _type = ePoint;
      _v1 = _v2;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to vertex (v2)- shouldn't be here?; reverting to point simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(w-1.0) < NEAR_ZERO)
    {
      _type = ePoint;
      _v1 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to vertex (v3)- shouldn't be here?; reverting to point simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(u+v-1.0) < NEAR_ZERO)
    {
      // edge ab
      _type = eSegment;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to edge (v1,v2)- shouldn't be here?; reverting to edge simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(v+w-1.0) < NEAR_ZERO)
    {
      // edge bc
      _type = eSegment;
      _v1 = _v2;
      _v2 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to edge (v2,v3)- shouldn't be here?; reverting to edge simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(u+w-1.0) < NEAR_ZERO)
    {
      // edge ac
      _type = eSegment;
      _v2 = _v3;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to edge (v1,v3)- shouldn't be here?; reverting to edge simplex" << std::endl;
      return cp;
    }
    else if (v < NEAR_ZERO && w < NEAR_ZERO)
    {
      // edge ad
      _type = eSegment;
      _v2 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to edge (v1,v4); reverting to edge simplex" << std::endl;
      return cp;
    }
    else if (u < NEAR_ZERO && w < NEAR_ZERO)
    {
      // edge bd
      _type = eSegment;
      _v1 = _v2;
      _v2 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to edge (v2,v4); reverting to edge simplex" << std::endl;
      return cp;
    }
    else if (u < NEAR_ZERO && v < NEAR_ZERO)
    {
      // edge cd
      _type = eSegment;
      _v1 = _v3;
      _v2 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to edge (v3,v4); reverting to edge simplex" << std::endl;
      return cp;
    }
    else if (std::fabs(u+v+w-1.0) < NEAR_ZERO)
    {
      // triangle abc
      _type = eTriangle;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to triangle (v1,v2,v3)- shouldn't be here? reverting to triangle simplex" << std::endl;
      return cp;
    }
    else if (w < NEAR_ZERO)
    {
      // triangle abd
      _type = eTriangle;
      _v3 = _v4;
      assert((find_closest()-cp).norm() < NEAR_ZERO);
      FILE_LOG(LOG_COLDET) << " -- closest point is to triangle (v1,v2,v4); reverting to triangle simplex" << std::endl;
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
      FILE_LOG(LOG_COLDET) << " -- closest point is to triangle (v2,v3,v4); reverting to triangle simplex" << std::endl;
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
      FILE_LOG(LOG_COLDET) << " -- closest point is to triangle (v1,v3,v4); reverting to triangle simplex" << std::endl;
      return cp;
    }
    else
    {
      // simplex is tetrahedron
      FILE_LOG(LOG_COLDET) << " -- closest point is to interior of tetrahedron" << std::endl;
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
    {
      std::swap(_v2.v, _v3.v);
      assert(Tetrahedron(_v1.v, _v2.v, _v3.v, _v4.v).calc_volume() > 0.0);
    }
  }
  else if (_type == eTetra)
    assert(false);
}

// TODO: remove this
/*
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
    FILE_LOG(LOG_COLDET) << "GJK::do_gjk() iteration: " << i << std::endl;

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
      FILE_LOG(LOG_COLDET) << "GJK::do_gjk() shapes are intersecting"  << std::endl;

      // A and B are intersecting
      // determine the interpenetration distance
      // THE CORRECT WAY to compute interpenetration for convex bodies
      // using GJK is to construct the minkowski difference (NOTE: not sure
      // if actual construction is necessary); the interpenetration distance
      // (corresponds to minimum translation necessary to separate the bodies)
      // will be the minimum l2-norm of a point on the boundary of the
      // Minkowski difference. Note that there may be several such points with
      // minimum norm (which would correspond to several points of deepest
      // interpenetration) 
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
    if (vdotd < min_dot)
    {
      min_dot = vdotd;
      closestA = pA;
      closestB = pB;
      if (vdotd < 0.0)
      {
        FILE_LOG(LOG_COLDET) << "GJK::do_gjk() dist=" << min_dist << ", exiting" << std::endl;
        return min_dist;
      }
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

  FILE_LOG(LOG_COLDET) << "GJK::do_gjk() [max iterations exceeded] dist=" << min_dist << ", exiting" << std::endl;

  return min_dist;
}
*/

/// Does GJK using primitives and poses defined in collision geometry frames
double GJK::do_gjk(shared_ptr<const Primitive> A, shared_ptr<const Primitive> B, shared_ptr<const Pose3d> PA, shared_ptr<const Pose3d> PB, Point3d& closestA, Point3d& closestB, unsigned max_iter)
{
  const double INF = std::numeric_limits<double>::max();

  // setup a random direction
  Point3d rdir((double) rand() / RAND_MAX * 2.0 - 1.0,(double) rand() / RAND_MAX * 2.0 - 1.0, (double) rand() / RAND_MAX * 2.0 - 1.0, GLOBAL);
  Point3d pA = A->get_supporting_point(-Pose3d::transform_vector(PA, rdir));
  Point3d pB = B->get_supporting_point(Pose3d::transform_vector(PB, rdir)); 

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
      FILE_LOG(LOG_COLDET) << "GJK::do_gjk() shapes are intersecting"  << std::endl;

      // A and B are intersecting
      // determine the interpenetration distance
      // THE CORRECT WAY to compute interpenetration for convex bodies
      // using GJK is to construct the minkowski difference (NOTE: not sure
      // if actual construction is necessary); the interpenetration distance
      // (corresponds to minimum translation necessary to separate the bodies)
      // will be the minimum l2-norm of a point on the boundary of the
      // Minkowski difference. Note that there may be several such points with
      // minimum norm (which would correspond to several points of deepest
      // interpenetration) 
      double pen_dist = INF;
      const unsigned NV = S.num_vertices();
      for (unsigned i=0; i< NV; i++)
      {
        double dA = A->calc_signed_dist(Pose3d::transform_point(PA, S.get_vertex(i).vB));
        if (dA < 0.0) 
          pen_dist = std::min(pen_dist, -dA);
        double dB = B->calc_signed_dist(Pose3d::transform_point(PB, S.get_vertex(i).vA));
        if (dB < 0.0) 
          pen_dist = std::min(pen_dist, -dB);
      }

      return -pen_dist;
    }
    // look for no progress
    else if (pnorm > min_dist-NEAR_ZERO)
    {
      FILE_LOG(LOG_COLDET) << "GJK::do_gjk() unable to progress!"  << std::endl;
      return pnorm;
    }

    // get the new supporting points and determine the new vertex
    Point3d pA = A->get_supporting_point(-Pose3d::transform_vector(PA, p));
    Point3d pB = B->get_supporting_point(Pose3d::transform_vector(PB, p)); 
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
    if (vdotd < 0.0)
    {
      FILE_LOG(LOG_COLDET) << "GJK::do_gjk() dist=" << min_dist << ", exiting" << std::endl;
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
        FILE_LOG(LOG_COLDET) << "GJK::do_gjk() added new point to simplex, now: " << oss.str() << std::endl;
      }
    }
  }

  throw std::runtime_error("maximum GJK iterations exceeded");
  return min_dist;
}

