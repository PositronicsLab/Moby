/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _GEARS_H
#define _GEARS_H

#include <boost/weak_ptr.hpp>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Pose3d.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/SAcceld.h>
#include <Ravelin/LinAlgd.h>
#include <Ravelin/DynamicBodyd.h>
#include <Moby/Base.h>
#include <Moby/RigidBody.h>
#include <Moby/Visualizable.h>
#include <Moby/Joint.h>

namespace Moby {

class VisualizationData;

/// Defines a gear ("joint") constraint 
/**
 * \todo implement a rest position for q?
 */
class Gears : public Joint 
{
  public:

    Gears();
    virtual ~Gears() {}
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void set_inboard_pose(boost::shared_ptr<const Ravelin::Pose3d> link, bool update_joint_pose);
    virtual void set_outboard_pose(boost::shared_ptr<const Ravelin::Pose3d> link, bool update_joint_pose);
    virtual void evaluate_constraints(double C[]);
    virtual void evaluate_constraints_dot(double C[]);
    virtual void calc_constraint_jacobian(bool inboard, Ravelin::MatrixNd& Cq);
    virtual void calc_constraint_jacobian_dot(bool inboard, Ravelin::MatrixNd& Cq);
    virtual unsigned num_dof() const { return 0; }
    virtual unsigned num_constraint_eqns() const { return 1; }
    virtual void determine_q(Ravelin::VectorNd& q) {}
    virtual boost::shared_ptr<const Ravelin::Pose3d> get_induced_pose() { return boost::shared_ptr<const Ravelin::Pose3d>(); }
    virtual bool is_singular_config() const { return false; }
    virtual const std::vector<Ravelin::SVelocityd>& get_spatial_axes_dot() { return _sdot; }

  private:
    double _ratio; // the input:output ratio for the gears

    std::vector<Ravelin::SVelocityd> _sdot; // this will be unused
}; // end class
} // end namespace

#endif

