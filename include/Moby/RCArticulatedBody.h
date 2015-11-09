/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _RC_ARTICULATED_BODY_H
#define _RC_ARTICULATED_BODY_H

#include <pthread.h>
#include <map>
#include <list>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/SForced.h>
#include <Moby/Constants.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Ravelin/RCArticulatedBodyd.h>

namespace Moby {

class Joint;
class UnilateralConstraintProblemData;

/// Defines an articulated body for use with reduced-coordinate dynamics algorithms
/**
 * Reduced-coordinate articulated bodies cannot rely upon the integrator to automatically update
 * the states (i.e., positions, velocities) of the links, as is done with maximal-coordinate 
 * articulated bodies.  Rather, the integrator updates the joint positions and velocities; the
 * states are obtained from this reduced-coordinate representation.
 * Notes about concurrency: <br /><br />
 *
 * It is generally desirable to be able to run forward dynamics and inverse 
 * dynamics algorithms concurrently to simulate actual robotic systems.   In 
 * general, derived classes should not operate on state variables
 * (joint positions, velocities, accelerations and floating base positions, 
 * velocites, and accelerations) directly during execution of the algorithm.  
 * Rather, derived classes should operate on copies of the state
 * variables, updating the state variables on conclusion of the algorithms.  
 */
class RCArticulatedBody : public virtual ArticulatedBody, public virtual Ravelin::RCArticulatedBodyd
{
  public:
    RCArticulatedBody();
    virtual ~RCArticulatedBody() {}
/*
    virtual Ravelin::MatrixNd& calc_jacobian(const Point3d& point, RigidBodyPtr link, Ravelin::MatrixNd& J);
    virtual Ravelin::MatrixNd& calc_jacobian(const Point3d& point, const Ravelin::Pose3d& base_pose, const std::map<JointPtr, Ravelin::VectorNd>& q, RigidBodyPtr link, Ravelin::MatrixNd& J);
*/
    RCArticulatedBodyPtr clone() const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void set_links_and_joints(const std::vector<RigidBodyPtr>& links, const std::vector<JointPtr>& joints);
    virtual void apply_generalized_impulse(const Ravelin::SharedVectorNd& gj);

  protected:
     virtual void compile();

  private:
    RCArticulatedBody(const RCArticulatedBody& rcab) {}
}; // end class

} // end namespace
#endif

