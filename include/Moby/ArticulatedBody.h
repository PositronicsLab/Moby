/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _ARTICULATED_BODY_H
#define _ARTICULATED_BODY_H

#include <vector>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/SMomentumd.h>
#include <Ravelin/sorted_pair>
#include <Moby/UnilateralConstraint.h>
#include <Moby/ControlledBody.h>
#include <Moby/Joint.h>
#include <Ravelin/ArticulatedBodyd.h>

namespace Moby {

class RigidBody;

/// Abstract class for articulated bodies
class ArticulatedBody : public virtual Ravelin::ArticulatedBodyd, public virtual ControlledBody
{
  public:
    ArticulatedBody();
    virtual ~ArticulatedBody() {}
    virtual void update_visualization();
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    void update_joint_constraint_violations();
    bool is_joint_constraint_violated() const;
    virtual void ode_noexcept(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data, Ravelin::SharedVectorNd& dx);
    virtual void prepare_to_calc_ode(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data);
    virtual void prepare_to_calc_ode_sustained_constraints(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data);
    virtual void ode(double t, double dt, void* data, Ravelin::SharedVectorNd& dx);
//    virtual unsigned num_generalized_coordinates(Ravelin::DynamicBodyd::GeneralizedCoordinateType gctype) const;

    /// Finds (joint) limit constraints 
    template <class OutputIterator>
    OutputIterator find_limit_constraints(OutputIterator begin) const;

    /// Gets shared pointer to this object as type ArticulatedBody
    ArticulatedBodyPtr get_this() { return boost::dynamic_pointer_cast<ArticulatedBody>(Ravelin::ArticulatedBodyd::shared_from_this()); }

    /// Gets shared pointer to this object as type const ArticulateBody
    boost::shared_ptr<const ArticulatedBody> get_this() const { return boost::dynamic_pointer_cast<const ArticulatedBody>(Ravelin::ArticulatedBodyd::shared_from_this()); }

  private:
    ArticulatedBody(const ArticulatedBody& ab) {}

    // joint constraint violation
    std::vector<double> _cvio;

    // joint velocity tolerances (for joints at constraints)
    std::vector<double> _cvel_vio;

    // temporary variables
    Ravelin::VectorNd _dq;

}; // end class

#include "ArticulatedBody.inl"

} // end namespace Moby

#endif
