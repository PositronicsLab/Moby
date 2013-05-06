/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/BV.h>
#include <Moby/OBB.h>
#include <Moby/BoundingSphere.h>
#include <Moby/AABB.h>
#include <Moby/SSR.h>
#include <Moby/SSL.h>
#include <Moby/Constants.h>
#include <Moby/DummyBV.h>

using std::pair;
using std::endl;
using namespace Ravelin;
using namespace Moby;

/// Computes the inverse of a relative pose
pair<Quatd, Origin3d> BV::inverse(const pair<Quatd, Origin3d>& p)
{
  pair<Quatd, Origin3d> invp;
  invp.first = Quatd::invert(p.first);
  invp.second = invp.first * -p.second;
  return invp;
}

/// Computes the distance between two abstract bounding volumes and stores the closest points
/**
 * \param cp1 the closest point on a to b
 * \param cp2 the closest point on b to a
 * \return the distance between the bounding volumes
 */
double BV::calc_distance(const BV* a, const BV* b, Point3d& cp1, Point3d& cp2)
{
  throw std::runtime_error("This method not implemented!");
  return -1.0;
}

/// Computes the distance between two abstract bounding volumes
/**
 * \param aTb the relative transformation from b's frame to a's frame
 * \param cp1 the closest point on a to b
 * \param cp2 the closest point on b to a
 * \return the distance between the bounding volumes
 */
double BV::calc_distance(const BV* a, const BV* b, const pair<Quatd, Origin3d>& aTb, Point3d& cp1, Point3d& cp2)
{
  throw std::runtime_error("This method not implemented!");
  return -1.0;
}

/// Computes whether two abstract bounding volumes intersect
bool BV::intersects(const BV* a, const BV* b)
{
  // look for dummy type
  if (dynamic_cast<const DummyBV*>(a) || dynamic_cast<const DummyBV*>(b))
    return true; 

  // look for OBB type
  if (dynamic_cast<const OBB*>(a))
  {
    // OBB / OBB intersection
    if (dynamic_cast<const OBB*>(b))
      return OBB::intersects(*((const OBB*) a), *((const OBB*) b));
    // OBB / AABB intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((const OBB*) a, (const AABB*) b);
    // OBB / SSR intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const OBB*) a, (const SSR*) b);
    // OBB / SSL intersection
    if (dynamic_cast<const SSL*>(b))
      return intersects((const OBB*) a, (const SSL*) b);
    // OBB / bounding sphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return intersects((const OBB*) a, (const BoundingSphere*) b);
  }
  // look for SSR type
  if (dynamic_cast<const SSR*>(a))
  {
    // OBB / SSR intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const SSR*) a);
    // SSR / SSR intersection
    if (dynamic_cast<const SSR*>(b))
      return SSR::intersects(*((const SSR*) a), *((const SSR*) b));
    // SSR / SSL intersection
    if (dynamic_cast<const SSL*>(b))
      return intersects((const SSR*) a, (const SSL*) b);
    // SSR / bounding sphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return intersects((const SSR*) a, (const BoundingSphere*) b);
    // SSR / AABB intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((const SSR*) a, (const AABB*) b);
  }
  // look for SSL type
  else if (dynamic_cast<const SSL*>(a))
  {
    // OBB / SSL intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const SSL*) a);
    // AABB / SSL intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((AABB*) b, (SSL*) a);
    // SSR / SSL intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const SSR*) b, (const SSL*) a);
    // SSL / SSL intersection
    if (dynamic_cast<const SSL*>(b))
      return SSL::intersects(*((const SSL*) a), *((const SSL*) b));
    // SSL / bounding sphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return intersects((const SSL*) a, (const BoundingSphere*) b);
  }
  // look for bounding sphere type
  else if (dynamic_cast<const BoundingSphere*>(a))
  {
    // OBB / bounding sphere intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const BoundingSphere*) a);
    // AABB / bounding sphere intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((AABB*) b, (BoundingSphere*) a);
    // SSR / bounding sphere intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const SSR*) b, (const BoundingSphere*) a);
    // SSL / bounding sphere intersection
    if (dynamic_cast<const SSL*>(b))
      return intersects((const SSL*) b, (const BoundingSphere*) a);
    // bounding sphere / bounding sphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return BoundingSphere::intersects(*((const BoundingSphere*) a), *((const BoundingSphere*) b));
  }
  else if (dynamic_cast<const AABB*>(a))
  {
    // AABB / OBB intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const AABB*) a);
    if (dynamic_cast<const SSR*>(b))
      return intersects((const SSR*) b, (const AABB*) a);
    if (dynamic_cast<const SSL*>(b))
      return intersects((const AABB*) a, (const SSL*) b);
    if (dynamic_cast<const BoundingSphere*>(b))
      return intersects((const AABB*) a, (const BoundingSphere*) b);
    if (dynamic_cast<const AABB*>(b))
      return AABB::intersects(*((const AABB*) a), *((const AABB*) b));
  }

  // shouldn't be here
  assert(false);
  return true;
}

/// Computes whether two abstract bounding volumes intersect
/**
 * \param aTb the relative transformation from b to a
 */
bool BV::intersects(const BV* a, const BV* b, const pair<Quatd, Origin3d>& aTb)
{
  // look for dummy type
  if (dynamic_cast<const DummyBV*>(a) || dynamic_cast<const DummyBV*>(b))
    return true; 

  // look for OBB type
  if (dynamic_cast<const OBB*>(a))
  {
    // OBB / OBB intersection
    if (dynamic_cast<const OBB*>(b))
      return OBB::intersects(*((const OBB*) a), *((const OBB*) b), aTb);
    // OBB / SSR intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const OBB*) a, (const SSR*) b, aTb);
    // OBB / SSL intersection
    if (dynamic_cast<const SSL*>(b))
      return intersects((const OBB*) a, (const SSL*) b, aTb);
    // OBB / bounding sphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return intersects((const OBB*) a, (const BoundingSphere*) b, aTb);
    // OBB / AABB intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((const OBB*) a, (const AABB*) b, aTb);
  }
  // look for SSR type
  else if (dynamic_cast<const SSR*>(a))
  {
    // SSR / SSR intersection
    if (dynamic_cast<const SSR*>(b))
      return SSR::intersects(*((const SSR*) a), *((const SSR*) b), aTb);
    // SSR / SSL intersection
    if (dynamic_cast<const SSL*>(b))
      return intersects((const SSR*) a, (const SSL*) b, aTb);
    // SSR / bounding sphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return intersects((const SSR*) a, (const BoundingSphere*) b, aTb);
    // OBB / SSR intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const SSR*) a, inverse(aTb));
    // AABB / SSR intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((const SSR*) a, (const AABB*) b, aTb);
  }
  else if (dynamic_cast<const SSL*>(a))
  {
    // OBB / SSL intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const SSL*) a, inverse(aTb));
    // SSR / SSL intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const SSR*) b, (const SSL*) a, inverse(aTb));
    // SSL / SSL intersection
    if (dynamic_cast<const SSL*>(b))
      return SSL::intersects(*((const SSL*) a), *((const SSL*) b), aTb);
    // SSL / BoundingSphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return intersects((const SSL*) a, (const BoundingSphere*) b, inverse(aTb));
    // SSL / AABB intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((const AABB*) b, (const SSL*) a, inverse(aTb));
  }
  else if (dynamic_cast<const BoundingSphere*>(a))
  {
    // OBB / bounding sphere intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const BoundingSphere*) a, inverse(aTb));
    // SSR / bounding sphere intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const SSR*) b, (const BoundingSphere*) a, inverse(aTb));
    // SSL / bounding sphere intersection
    if (dynamic_cast<const SSL*>(b))
      return intersects((const SSL*) b, (const BoundingSphere*) a, inverse(aTb));
    // bounding sphere / bounding sphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return BoundingSphere::intersects(*((const BoundingSphere*) a), *((const BoundingSphere*) b), aTb);
    // AABB / bounding sphere intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((const AABB*) b, (const BoundingSphere*) a, inverse(aTb));
  }
  else if (dynamic_cast<const AABB*>(a))
  {
    // OBB / AABB intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const AABB*) a, inverse(aTb));
    // SSR / AABB intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const SSR*) b, (const AABB*) a, inverse(aTb));
    // SSL / AABB intersection
    if (dynamic_cast<const SSL*>(b))
      return intersects((const AABB*) a, (const SSL*) b, aTb);
    // BoundingSphere / AABB intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return intersects((const AABB*) a, (const BoundingSphere*) b, aTb);
    // AABB / AABB intersection
    if (dynamic_cast<const AABB*>(b))
      return AABB::intersects(*((const AABB*) a), *((const AABB*) b), aTb);
  }

  // shouldn't be here
  assert(false);
  return true;
}

/// Tests intersection between an OBB and an AABB
bool BV::intersects(const OBB* O, const AABB* A, const pair<Quatd, Origin3d>& OTA)
{
  OBB Ao = A->get_OBB();

  return OBB::intersects(*O, Ao, OTA);
}

/// Tests intersection between an OBB and an AABB
bool BV::intersects(const OBB* O, const AABB* A)
{
  // make an OBB from the AABB
  OBB Ao = A->get_OBB();

  return OBB::intersects(*O, Ao);
}

/// Tests intersection between an OBB and a bounding sphere
bool BV::intersects(const OBB* O, const BoundingSphere* S)
{
  const unsigned THREE_D = 3;

  // transform the sphere center to OBB space
  Point3d center = O->R.transpose_mult(S->center - O->center);

  FILE_LOG(LOG_COLDET) << "  -- sphere center: " << S->center << endl;
  FILE_LOG(LOG_COLDET) << "  -- sphere center: " << center << " (OBB frame)" << endl;

  // find the square of the distance from the sphere to the box
  double d = 0;
  for (unsigned i=0; i< THREE_D; i++)
  {
    if (center[i] < -O->l[i])
    {
      double s = center[i] + O->l[i];
      d += s*s;
    }
    else if (center[i] > O->l[i])
    {
      double s = center[i] - O->l[i];
      d += s*s;
    }
  }

  FILE_LOG(LOG_COLDET) << " -- squared distance (negative indicates interpenetration): " << (d - S->radius*S->radius) << endl;
  return d <= S->radius*S->radius;
}

/// Tests intersection between an OBB and a bounding sphere
/**
 * \param OTS the matrix transforming S's frame to O's frame
 */
bool BV::intersects(const OBB* O, const BoundingSphere* S, const pair<Quatd, Origin3d>& OTS)
{
  // create a new bounding sphere in O's frame
  BoundingSphere s = *S;
  s.center = OTS.first * s.center + OTS.second;
  return intersects(O, &s);
}

/// Checks for intersection between a AABB and a bounding sphere
bool BV::intersects(const AABB* A, const BoundingSphere* S)
{
  const unsigned THREE_D = 3;

  FILE_LOG(LOG_COLDET) << "BV::intersects() [AABB/sphere] entered" << endl;

  // transform the sphere center to OBB space
  Point3d center = S->center - A->minp + A->maxp;

  FILE_LOG(LOG_COLDET) << "  -- sphere center: " << S->center << endl;
  FILE_LOG(LOG_COLDET) << "  -- sphere center: " << center << " (OBB frame)" << endl;

  // get the half-lengths of the AABB
  Point3d l = A->maxp*0.5 - A->minp*0.5;

  // find the square of the distance from the sphere to the box
  double d = 0;
  for (unsigned i=0; i< THREE_D; i++)
  {
    if (center[i] < -l[i])
    {
      double s = -center[i] + l[i];
      d += s*s;
    }
    else if (center[i] > l[i])
    {
      double s = center[i] - l[i];
      d += s*s;
    }
  }

  FILE_LOG(LOG_COLDET) << " -- dist: " << d << endl;

  return d <= S->radius*S->radius;
}

/// Tests intersection between a AABB and a bounding sphere
/**
 * \param ATS the matrix transforming S's frame to A's frame
 */
bool BV::intersects(const AABB* A, const BoundingSphere* S, const pair<Quatd, Origin3d>& ATS)
{
  // create a new bounding sphere in O's frame
  BoundingSphere s = *S;
  s.center = ATS.first * s.center + ATS.second;
  return intersects(A, &s);
}

/// Tests intersection between a SSR and a bounding sphere
/**
 * \param S the sphere-swept rectangle
 * \param B the bounding sphere
 */
bool BV::intersects(const SSR* S, const BoundingSphere* B)
{
  // determine the distance between S and center of the bounding sphere
  double dist = SSR::calc_dist(*S, B->center);

  // check whether the distance is within the radius of the bounding sphere
  return dist - B->radius <= (double) 0.0;
}

/// Tests intersection between a SSR and a bounding sphere
/**
 * \param S the sphere-swept rectangle
 * \param B the bounding sphere
 * \param STB transformation from B's frame to S's frame
 */
bool BV::intersects(const SSR* S, const BoundingSphere* B, const pair<Quatd, Origin3d>& STB)
{
  // transform the center of the bounding sphere
  Point3d xc = STB.first * B->center + STB.second;

  // determine the distance between S and xformed center of the bounding sphere
  double dist = SSR::calc_dist(*S, xc);

  // check whether the distance is within the radius of the bounding sphere
  return dist - B->radius <= (double) 0.0;
}

/// Tests intersection between an OBB and a SSR
/**
 * \param O the oriented bounding box
 * \param S the sphere-swept rectangle
 */
bool BV::intersects(const OBB* O, const SSR* S)
{
  // create a AABB around the SSR
  AABB Sx;
  Sx.minp = ((SSR*) S)->get_lower_bounds(Pose3d::identity());
  Sx.maxp = ((SSR*) S)->get_upper_bounds(Pose3d::identity());
  return intersects(O, &Sx);
}

/// Tests intersection between an OBB and a SSR
/**
 * \param O the oriented bounding box
 * \param S the sphere-swept rectangle
 * \param OTS the transformation from S's frame to O's frame
 */
bool BV::intersects(const OBB* O, const SSR* S, const pair<Quatd, Origin3d>& OTS)
{
  AABB Sx;
  Pose3d T(OTS.first, OTS.second);
  Sx.minp = ((SSR*) S)->get_lower_bounds(T);
  Sx.maxp = ((SSR*) S)->get_upper_bounds(T);
  return intersects(O, &Sx);
}

/// Tests intersection between a SSR and a AABB
/**
 * \param S the sphere-swept rectangle
 * \param A the axis-aligned bounding box
 * \param STA the transformation from A's frame to S's frame
 */
bool BV::intersects(const SSR* S, const AABB* A, const pair<Quatd, Origin3d>& STA)
{
  pair<Quatd, Origin3d> ATS = inverse(STA);
  AABB Sx;
  Pose3d T(ATS.first, ATS.second);
  Sx.minp = ((SSR*) S)->get_lower_bounds(T);
  Sx.maxp = ((SSR*) S)->get_upper_bounds(T);
  return AABB::intersects(*A, Sx);
}

/// Tests intersection between a SSR and a AABB
/**
 * \param S the sphere-swept rectangle
 * \param A the axis-aligned bounding box
 */
bool BV::intersects(const SSR* S, const AABB* A)
{
  AABB Sx;
  Sx.minp = ((SSR*) S)->get_lower_bounds(Pose3d::identity());
  Sx.maxp = ((SSR*) S)->get_upper_bounds(Pose3d::identity());
  return AABB::intersects(*A, Sx);
}

/// Tests intersection between a SSL and a AABB
bool BV::intersects(const AABB* A, const SSL* B)
{
  AABB Bx;
  Bx.minp = ((SSL*) B)->get_lower_bounds(Pose3d::identity());
  Bx.maxp = ((SSL*) B)->get_upper_bounds(Pose3d::identity());
  return AABB::intersects(*A, Bx);
}

/// Tests intersection between a SSL and a AABB
bool BV::intersects(const AABB* A, const SSL* B, const pair<Quatd, Origin3d>& aTb)
{
  AABB Bx;
  Pose3d T(aTb.first, aTb.second);
  Bx.minp = ((SSL*) B)->get_lower_bounds(T);
  Bx.maxp = ((SSL*) B)->get_upper_bounds(T);
  return AABB::intersects(*A, Bx);
}

/// Tests intersection between a SSL and an OBB
bool BV::intersects(const OBB* A, const SSL* B)
{
  AABB Bx;
  Bx.minp = ((SSL*) B)->get_lower_bounds(Pose3d::identity());
  Bx.maxp = ((SSL*) B)->get_upper_bounds(Pose3d::identity());
  return intersects(A, &Bx);
}

/// Tests intersection between a SSL and an OBB
bool BV::intersects(const OBB* A, const SSL* B, const pair<Quatd, Origin3d>& aTb)
{
  AABB Bx;
  Pose3d T(aTb.first, aTb.second);
  Bx.minp = ((SSL*) B)->get_lower_bounds(T);
  Bx.maxp = ((SSL*) B)->get_upper_bounds(T);
  return intersects(A, &Bx);
}

/// Tests intersection between a SSL and a bounding sphere
bool BV::intersects(const SSL* A, const BoundingSphere* B)
{
  double dist = SSL::calc_dist(*A, B->center);
  return dist <= B->radius;
}

/// Tests intersection between a SSL and a bounding sphere
bool BV::intersects(const SSL* A, const BoundingSphere* B, const pair<Quatd, Origin3d>& aTb)
{
  double dist = SSL::calc_dist(*A, B->center);
  return dist <= B->radius;
}

/// Tests intersection between a SSR and a SSL
bool BV::intersects(const SSR* A, const SSL* B)
{
  double dist = SSR::calc_dist(*A, LineSeg3(B->p1, B->p2));
  return dist <= B->radius;
}

/// Tests intersection between a SSR and a SSL
bool BV::intersects(const SSR* A, const SSL* B, const pair<Quatd, Origin3d>& aTb)
{
  Point3d Bp1 = aTb.first * B->p1 + aTb.second;
  Point3d Bp2 = aTb.first * B->p2 + aTb.second;
  double dist = SSR::calc_dist(*A, LineSeg3(Bp1, Bp2));
  return dist <= B->radius;
}


