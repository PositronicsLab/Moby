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
#include <Ravelin/Point3d.h>
#include <Ravelin/Vector3d.h>
#include <Ravelin/AAngled.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/LinAlgd.h>
#include <Moby/BV.h>
#include <Moby/Constants.h>
#include <Moby/Types.h>
#include <Moby/FastThreadable.h>
#include <Moby/CompGeom.h>
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
    OBB(const Ravelin::Point3d& center, const Ravelin::Matrix3d& R, const Ravelin::Vector3d& l);
    OBB(const OBB& o, const Ravelin::Vector3d& v);
    void operator=(const OBB& obb);
    virtual BVPtr calc_vel_exp_BV(CollisionGeometryPtr g, double dt, const Ravelin::Vector3d& lv, const Ravelin::Vector3d& av) const;
    static double calc_sq_dist(const OBB& o, const Ravelin::Point3d& p);
    static double calc_dist(const OBB& a, const OBB& b, Ravelin::Point3d& cpa, Ravelin::Point3d& cpb);
    static double calc_dist(const OBB& a, const OBB& b, const Ravelin::Pose3d& aTb, Ravelin::Point3d& cpa, Ravelin::Point3d& cpb);
    static bool intersects(const OBB& a, const OBB& b);
    static bool intersects(const OBB& a, const OBB& b, const Ravelin::Pose3d& T);
    static bool intersects(const OBB& a, const LineSeg3& seg, double& tmin, double tmax, Ravelin::Point3d& q);
    virtual bool intersects(const LineSeg3& seg, double& tmin, double tmax, Ravelin::Point3d& q) const { return OBB::intersects(*this, seg, tmin, tmax, q); }
    static bool outside(const OBB& a, const Ravelin::Point3d& point, double tol = NEAR_ZERO);
    virtual bool outside(const Ravelin::Point3d& point, double tol = NEAR_ZERO) const { return OBB::outside(*this, point, tol); }
    OBBPtr get_this() { return boost::dynamic_pointer_cast<OBB>(shared_from_this()); }
    boost::shared_ptr<const OBB> get_this() const { return boost::dynamic_pointer_cast<const OBB>(shared_from_this()); }
    virtual std::ostream& to_vrml(std::ostream& out, const Ravelin::Pose3d& T) const;
    unsigned calc_size() const;
    XMLTreePtr save_to_xml_tree() const;
    virtual Ravelin::Point3d get_lower_bounds(const Ravelin::Pose3d& T);
    virtual Ravelin::Point3d get_upper_bounds(const Ravelin::Pose3d& T);
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
    virtual double calc_volume() const { return l[0] * l[1] * l[2]; }

    /// Center of the bounding box
    Ravelin::Point3d center;

    /// Half-lengths of the OBBs in the three axes directions (i.e., in the box frame)
    Ravelin::Vector3d l;

    /// Orientation of this OBB
    Ravelin::Matrix3d R;

  private:
    template <class ForwardIterator>
    static OBB calc_low_dim_OBB(ForwardIterator begin, ForwardIterator end);

    template <class ForwardIterator>
    static void calc_lengths(const Ravelin::Vector3d& d1, const Ravelin::Vector3d& d2, const Ravelin::Vector3d& d3, const Ravelin::Point3d& center, ForwardIterator begin, ForwardIterator end, double lengths[3]);

    template <class ForwardIterator>
    static void align(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& d1, Ravelin::Vector3d& d2);
}; // end class

// include inline functions
#include "OBB.inl"

} // end namespace

#endif

