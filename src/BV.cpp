/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
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
double BV::calc_distance(const BV* a, const BV* b, const Transform3d& aTb, Point3d& cp1, Point3d& cp2)
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
bool BV::intersects(const BV* a, const BV* b, const Transform3d& aTb)
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
      return intersects((const OBB*) b, (const SSR*) a, aTb.inverse());
    // AABB / SSR intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((const SSR*) a, (const AABB*) b, aTb);
  }
  else if (dynamic_cast<const SSL*>(a))
  {
    // OBB / SSL intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const SSL*) a, aTb.inverse());
    // SSR / SSL intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const SSR*) b, (const SSL*) a, aTb.inverse());
    // SSL / SSL intersection
    if (dynamic_cast<const SSL*>(b))
      return SSL::intersects(*((const SSL*) a), *((const SSL*) b), aTb);
    // SSL / BoundingSphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return intersects((const SSL*) a, (const BoundingSphere*) b, aTb.inverse());
    // SSL / AABB intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((const AABB*) b, (const SSL*) a, aTb.inverse());
  }
  else if (dynamic_cast<const BoundingSphere*>(a))
  {
    // OBB / bounding sphere intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const BoundingSphere*) a, aTb.inverse());
    // SSR / bounding sphere intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const SSR*) b, (const BoundingSphere*) a, aTb.inverse());
    // SSL / bounding sphere intersection
    if (dynamic_cast<const SSL*>(b))
      return intersects((const SSL*) b, (const BoundingSphere*) a, aTb.inverse());
    // bounding sphere / bounding sphere intersection
    if (dynamic_cast<const BoundingSphere*>(b))
      return BoundingSphere::intersects(*((const BoundingSphere*) a), *((const BoundingSphere*) b), aTb);
    // AABB / bounding sphere intersection
    if (dynamic_cast<const AABB*>(b))
      return intersects((const AABB*) b, (const BoundingSphere*) a, aTb.inverse());
  }
  else if (dynamic_cast<const AABB*>(a))
  {
    // OBB / AABB intersection
    if (dynamic_cast<const OBB*>(b))
      return intersects((const OBB*) b, (const AABB*) a, aTb.inverse());
    // SSR / AABB intersection
    if (dynamic_cast<const SSR*>(b))
      return intersects((const SSR*) b, (const AABB*) a, aTb.inverse());
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
bool BV::intersects(const OBB* O, const AABB* A, const Transform3d& OTA)
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

  // verify that the relative poses are identical
  assert(O->get_relative_pose() == S->get_relative_pose());

  // transform the sphere center to OBB space
  Origin3d center = O->R.transpose_mult(Origin3d(S->center - O->center));

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
bool BV::intersects(const OBB* O, const BoundingSphere* S, const Transform3d& OTS)
{
  // create a new bounding sphere in O's frame
  BoundingSphere s = *S;
  s.center = OTS.transform_point(S->center);
  return intersects(O, &s);
}

/// Checks for intersection between a AABB and a bounding sphere
bool BV::intersects(const AABB* A, const BoundingSphere* S)
{
  const unsigned THREE_D = 3;

  FILE_LOG(LOG_COLDET) << "BV::intersects() [AABB/sphere] entered" << endl;

  // transform the sphere center to OBB space
  assert(A->get_relative_pose() == S->get_relative_pose());
  Point3d center = S->center;
  center.pose = A->get_relative_pose();
  center += A->maxp - A->minp;

  FILE_LOG(LOG_COLDET) << "  -- sphere center: " << center << " (AABB frame)" << endl;

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
bool BV::intersects(const AABB* A, const BoundingSphere* S, const Transform3d& ATS)
{
  // create a new bounding sphere in O's frame
  BoundingSphere s = *S;
  s.center = ATS.transform_point(S->center);
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
bool BV::intersects(const SSR* S, const BoundingSphere* B, const Transform3d& STB)
{
  // transform the center of the bounding sphere
  Point3d xc = STB.transform_point(B->center);

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
  AABB S_aabb;
  S_aabb.minp = S->get_lower_bounds();
  S_aabb.maxp = S->get_upper_bounds();
  return intersects(O, &S_aabb);
}

/// Tests intersection between an OBB and a SSR
/**
 * \param O the oriented bounding box
 * \param S the sphere-swept rectangle
 * \param OTS the transformation from S's frame to O's frame
 */
bool BV::intersects(const OBB* O, const SSR* S, const Transform3d& OTS)
{
  SSR Sx;
  AABB S_aabb;
  S->transform(OTS, &Sx);
  S_aabb.minp = Sx.get_lower_bounds();
  S_aabb.maxp = Sx.get_upper_bounds();
  return intersects(O, &S_aabb);
}

/// Tests intersection between a SSR and a AABB
/**
 * \param S the sphere-swept rectangle
 * \param A the axis-aligned bounding box
 * \param STA the transformation from A's frame to S's frame
 */
bool BV::intersects(const SSR* S, const AABB* A, const Transform3d& STA)
{
  Transform3d ATS = STA.inverse();
  SSR Sx;
  AABB S_aabb;
  S->transform(ATS, &Sx);
  S_aabb.minp = Sx.get_lower_bounds();
  S_aabb.maxp = Sx.get_upper_bounds();
  return AABB::intersects(*A, S_aabb);
}

/// Tests intersection between a SSR and a AABB
/**
 * \param S the sphere-swept rectangle
 * \param A the axis-aligned bounding box
 */
bool BV::intersects(const SSR* S, const AABB* A)
{
  AABB Sx;
  Sx.minp = S->get_lower_bounds();
  Sx.maxp = S->get_upper_bounds();
  return AABB::intersects(*A, Sx);
}

/// Tests intersection between a SSL and a AABB
bool BV::intersects(const AABB* A, const SSL* B)
{
  AABB Bx;
  Bx.minp = B->get_lower_bounds();
  Bx.maxp = B->get_upper_bounds();
  return AABB::intersects(*A, Bx);
}

/// Tests intersection between a SSL and a AABB
bool BV::intersects(const AABB* A, const SSL* B, const Transform3d& aTb)
{
  AABB B_aabb;
  SSL Bx;

  B->transform(aTb, &Bx);
  B_aabb.minp = Bx.get_lower_bounds();
  B_aabb.maxp = Bx.get_upper_bounds();

  return AABB::intersects(*A, B_aabb);
}

/// Tests intersection between a SSL and an OBB
bool BV::intersects(const OBB* A, const SSL* B)
{
  AABB Bx;
  Bx.minp = B->get_lower_bounds();
  Bx.maxp = B->get_upper_bounds();
  return intersects(A, &Bx);
}

/// Tests intersection between a SSL and an OBB
bool BV::intersects(const OBB* A, const SSL* B, const Transform3d& aTb)
{
  SSL Bx;
  AABB B_aabb;
  B->transform(aTb, &Bx);
  B_aabb.minp = Bx.get_lower_bounds();
  B_aabb.maxp = Bx.get_upper_bounds();
  return intersects(A, &B_aabb);
}

/// Tests intersection between a SSL and a bounding sphere
bool BV::intersects(const SSL* A, const BoundingSphere* B)
{
  double dist = SSL::calc_dist(*A, B->center);
  return dist <= B->radius;
}

/// Tests intersection between a SSL and a bounding sphere
bool BV::intersects(const SSL* A, const BoundingSphere* B, const Transform3d& aTb)
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
bool BV::intersects(const SSR* A, const SSL* B, const Transform3d& aTb)
{
  Point3d Bp1 = aTb.transform_point(B->p1);
  Point3d Bp2 = aTb.transform_point(B->p2);
  double dist = SSR::calc_dist(*A, LineSeg3(Bp1, Bp2));
  return dist <= B->radius;
}

