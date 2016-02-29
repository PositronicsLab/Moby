/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
#include <iostream>
#include <boost/shared_ptr.hpp>
//#include <Moby/Simulator.h>
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/GravityForce.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/PrismaticJoint.h>
#include <Moby/Log.h>

#include "viewer.h"

//----------------------------------------------------------------------------
// the push controller applies an impulse to the link once
/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
Ravelin::VectorNd& push_controller( Moby::ControlledBodyPtr dbp, Ravelin::VectorNd& u, double t, void* ) {
  static bool pushed = false;

  // only apply the force once
  if( pushed ) return u;

  // cast the dynamic body to an articulated body
  Moby::ArticulatedBodyPtr ab = boost::dynamic_pointer_cast<Moby::ArticulatedBody>( dbp );
  if( !ab ) {
    std::cout << "Failed to cast DynamicBody as ArticulatedBody" << std::endl;
    std::cout << "Failed to push link" << std::endl;
    pushed = true;  // disable the controller
    return u;
  } 
  
  const std::vector<boost::shared_ptr<Ravelin::Jointd> >& js = ab->get_joints();
  boost::shared_ptr<Ravelin::Jointd> joint;
  for(std::vector<boost::shared_ptr<Ravelin::Jointd> >::const_iterator it = js.begin(); it != js.end(); it++) {
    if( (*it)->joint_id == "joint" ) joint = *it;
  }
  if( !link ) {
    std::cout << "Failed to find joint" << std::endl;
    std::cout << "Failed to push joint" << std::endl;
    pushed = true;  // disable the controller
    return u;
  } 

  u.resize(1);
  u[0] = 1000;

  // disable the controller
  pushed = true; 
  return u;
}

//----------------------------------------------------------------------------
int main( void ) {

// uncomment to log dynamics
  Moby::Log<Moby::OutputToFile>::reporting_level = 7;

  boost::shared_ptr<Moby::TimeSteppingSimulator> sim( new Moby::TimeSteppingSimulator() );

  boost::shared_ptr<Moby::GravityForce> g( new Moby::GravityForce() );
  g->gravity = Ravelin::Vector3d( 0, 0, -9.8 );

  Moby::RCArticulatedBodyPtr ab( new Moby::RCArticulatedBody() );
  ab->id = "gripper";
  ab->algorithm_type = Moby::RCArticulatedBody::eCRB;

  std::vector< Moby::RigidBodyPtr > links;
  Moby::RigidBodyPtr base( new Moby::RigidBody() );
  {
    Moby::PrimitivePtr box( new Moby::BoxPrimitive(0.05,1.0,0.05) );
    box->set_mass( 1 );

    // static
    base->id = "base";
    base->set_visualization_data( box->create_visualization() );
    base->set_inertia( box->get_inertia() );
    base->set_enabled( false );
  
    base->set_pose( Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(0,0,0,1)), Ravelin::Origin3d(0,0,0) ) ); 
    links.push_back( base );
  }

  Moby::RigidBodyPtr link( new Moby::RigidBody() );
  {
    Moby::PrimitivePtr box( new Moby::BoxPrimitive(0.05,0.05,0.5) );
    box->set_mass( 1 );

    link->id = "link";
    link->set_visualization_data( box->create_visualization() );
    link->set_inertia( box->get_inertia() );
    link->set_enabled( true );
    link->get_recurrent_forces().push_back( g );
  
    link->set_pose( Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(0,0,0,1)), Ravelin::Origin3d(0,0,-0.275) ) ); 
    links.push_back( link ); 
  }

  std::vector< Moby::JointPtr > joints;
  boost::shared_ptr<Moby::PrismaticJoint> joint( new Moby::PrismaticJoint() );
  {
    joint->id = "joint";
    joint->set_location( Ravelin::Vector3d(0,0,0,base->get_pose()), base, link );
    joint->set_axis( Ravelin::Vector3d(0,1,0,Moby::GLOBAL) );
    joint->lolimit = -0.5;
    joint->hilimit = 0.5;

    joints.push_back( joint );
  }

  ab->set_links_and_joints( links, joints );
  ab->get_recurrent_forces().push_back( g );
  ab->set_floating_base(false);
  ab->controller = &push_controller;

  sim->add_dynamic_body( ab );

  Viewer viewer( sim, Ravelin::Origin3d(-5,0,-1), Ravelin::Origin3d(0,0,-1), Ravelin::Origin3d(0,0,1) );

  while(true) {
    if( !viewer.update() ) break;

    sim->step( 0.001 );
    Ravelin::Pose3d pose = *link->get_pose();
    pose.update_relative_pose(Moby::GLOBAL);
    std::cout << "t: " << sim->current_time << " x: " << pose.x << std::endl;
  }

  return 0;
}
//----------------------------------------------------------------------------

