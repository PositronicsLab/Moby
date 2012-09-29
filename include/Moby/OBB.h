/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_OBB_H_
#define _MOBY_OBB_H_

#include <stack>
#include <iostream>
#include <queue>
#include <boost/tuple/tuple.hpp>
#include <boost/enable_shared_from_this.hpp>
#include <Moby/AAngle.h>
#include <Moby/BV.h>
#include <Moby/Constants.h>
#include <Moby/Types.h>
#include <Moby/FastThreadable.h>
#include <Moby/CompGeom.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix4.h>
//#include <Moby/TriangleMeshPrimitive.h>

namespace Moby {

/// An oriented bounding box that optionally allows building an OBB tree
/**
 * \note the OBB is generally constructed such that its frame is aligned with
 *       that of the underlying rigid body (or collision geometry).  
 *       Therefore, the center of the OBB is computed relative to the 
 *       center-of-mass of the body (or the center of the geometry).  The
 *       orientation of the OBB will always be identical to the orientation of
 *       the rigid body (or collision geometry).
 */
class OBB : public BV 
{
  public:
    OBB();
    OBB(const OBB& obb) { operator=(obb); }
    OBB(const Vector3& center, const Matrix3& R, const Vector3& l);
    OBB(const OBB& o, const Vector3& v);
    void operator=(const OBB& obb);
    virtual BVPtr calc_vel_exp_BV(CollisionGeometryPtr g, Real dt, const Vector3& lv, const Vector3& av) const;
    static Real calc_sq_dist(const OBB& o, const Vector3& p);
    static Real calc_dist(const OBB& a, const OBB& b, Vector3& cpa, Vector3& cpb);
    static Real calc_dist(const OBB& a, const OBB& b, const Matrix4& aTb, Vector3& cpa, Vector3& cpb);
    static bool intersects(const OBB& a, const OBB& b);
    static bool intersects(const OBB& a, const OBB& b, const Matrix4& T);
    static bool intersects(const OBB& a, const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q);
    virtual bool intersects(const LineSeg3& seg, Real& tmin, Real tmax, Vector3& q) const { return OBB::intersects(*this, seg, tmin, tmax, q); }
    static bool outside(const OBB& a, const Vector3& point, Real tol = NEAR_ZERO);
    virtual bool outside(const Vector3& point, Real tol = NEAR_ZERO) const { return OBB::outside(*this, point, tol); }
    OBBPtr get_this() { return boost::dynamic_pointer_cast<OBB>(shared_from_this()); }
    boost::shared_ptr<const OBB> get_this() const { return boost::dynamic_pointer_cast<const OBB>(shared_from_this()); }
    virtual std::ostream& to_vrml(std::ostream& out, const Matrix4& T) const;
    unsigned calc_size() const;
    XMLTreePtr save_to_xml_tree() const;
    virtual Vector3 get_lower_bounds(const Matrix4& T);
    virtual Vector3 get_upper_bounds(const Matrix4& T);
    static OBBPtr load_from_xml(XMLTreeConstPtr root);

    template <class ForwardIterator>
    void expand_to_fit(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    OBB(ForwardIterator begin, ForwardIterator end);
 
    template <class ForwardIterator>
    static OBB calc_min_volume_OBB(ForwardIterator begin, ForwardIterator end);

    template <class OutputIterator>
    OutputIterator get_vertices(OutputIterator begin) const;

    /// Calculates 1/8th of the volume of the bounding box
    virtual Real calc_volume() const { return l[0] * l[1] * l[2]; }

    /// Center of the bounding box
    Vector3 center;

    /// Half-lengths of the OBBs in the three axes directions (i.e., in the box frame)
    Vector3 l;

    /// Orientation of this OBB
    Matrix3 R;

  private:
    template <class ForwardIterator>
    static OBB calc_low_dim_OBB(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    static void calc_lengths(const Vector3& d1, const Vector3& d2, const Vector3& d3, const Vector3& center, ForwardIterator begin, ForwardIterator end, Real lengths[3]);

    template <class ForwardIterator>
    static void align(ForwardIterator begin, ForwardIterator end, const Vector3& d1, Vector3& d2);
}; // end class

// include inline functions
#include "OBB.inl"

} // end namespace

#endif

