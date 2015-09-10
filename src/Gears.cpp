/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <cmath>
#include <iostream>
#include <Ravelin/RevoluteJointd.h>
#include <Ravelin/ArticulatedBodyd.h>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/XMLTree.h>
#include <Moby/UndefinedAxisException.h>
#include <Moby/Gears.h>

using boost::const_pointer_cast;
using boost::dynamic_pointer_cast;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

/// Initializes the joint
/**
 * The gear ratio is set to 1:1 
 */
Gears::Gears() : Joint()
{
  init_data();  
  _ratio = 1.0;
}

/// Evaluates the joint constraint
void Gears::evaluate_constraints(double C[])
{
  // constraint always evaluates to zero
  C[0] = 0.0;
}

/// Evaluates the time derivative of the constraint
void Gears::evaluate_constraint_dot(double C[])
{
  // constraint always evaluates to zero
  C[0] = 0.0;
}

/// Sets the inboard pose 
void Gears::set_inboard_pose(shared_ptr<const Pose3d> inboard_pose, bool update_joint_pose)
{
  // call parent method since it does all of the work
  shared_ptr<Pose3d> ib_nc = const_pointer_cast<Pose3d>(inboard_pose);
  Jointd::set_inboard_pose(ib_nc, update_joint_pose);
}

/// Sets the outboard pose 
void Gears::set_outboard_pose(shared_ptr<const Pose3d> outboard_pose, bool update_joint_pose)
{
  // call parent method since it does all of the work
  shared_ptr<Pose3d> ob_nc = const_pointer_cast<Pose3d>(outboard_pose);
  Jointd::set_outboard_pose(ob_nc, update_joint_pose);
}

/// Computes the constraint Jacobian
void Gears::calc_constraint_jacobian(bool inboard, MatrixNd& Cq)
{
  // get the number of coordinates in the articulated body
  const unsigned NGC = get_inboard_link()->get_articulated_body()->num_generalized_coordinates(DynamicBodyd::eSpatial);

  // clear the matrix 
  Cq.set_zero(1, NGC);

  // set the appropriate entry of the Jacobian
  if (inboard)
  {
    // get the revolute joint below the inboard link
    shared_ptr<Jointd> rji = get_inboard_link()->get_inner_joint_explicit();
    assert(dynamic_pointer_cast<RevoluteJointd>(rji));

    // get the index
    unsigned index = rji->get_coord_index();

    // set the coordinate index for this joint
    *(Cq.column_iterator_begin()+index) = 1.0; 
  }
  else
  {
    // get the revolute joint below the outboard link
    shared_ptr<Jointd> rjo = get_outboard_link()->get_inner_joint_explicit();
    assert(dynamic_pointer_cast<RevoluteJointd>(rjo));

    // get the index
    unsigned index = rjo->get_coord_index();

    // set the coordinate index for this joint
    *(Cq.column_iterator_begin()+index) = -_ratio; 
  }
}

/// Computes the time derivative of the constraint Jacobian
void Gears::calc_constraint_jacobian_dot(bool inboard, MatrixNd& Cq)
{
  // this will not be implemented
  assert(false); 
}

/// Implements Base::load_from_xml()
void Gears::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // read the information from the articulated body joint
  Joint::load_from_xml(node, id_map);

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "Gears") == 0);

  // read the gear ratio, if given
  XMLAttrib* gear_ratio_attrib = node->get_attrib("gear-ratio");
  if (gear_ratio_attrib)
    _ratio = gear_ratio_attrib->get_real_value();
}

/// Implements Base::save_to_xml()
void Gears::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // get info from Joint::save_to_xml()
  Joint::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "Gears";

  // write the gear ratio 
  node->attribs.insert(XMLAttrib("gear-ratio", _ratio));
}

