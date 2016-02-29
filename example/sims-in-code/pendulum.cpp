/*----------------------------------------------------------------------------
----------------------------------------------------------------------------*/
#include <iostream>
#include <boost/shared_ptr.hpp>
#include <Moby/Simulator.h>
#include <Moby/GravityForce.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/Log.h>

#include "viewer.h"

//----------------------------------------------------------------------------
int main( void ) {

// uncomment to log dynamics
//  Moby::Log<Moby::OutputToFile>::reporting_level = 7;

  // create a simulator
  boost::shared_ptr<Moby::Simulator> sim( new Moby::Simulator() );

  // create a gravity vector
  boost::shared_ptr<Moby::GravityForce> g( new Moby::GravityForce() );
  g->gravity = Ravelin::Vector3d( 0, 0, -9.8 );

  // create an articulated body for a kinematic chain
  Moby::RCArticulatedBodyPtr pendulum( new Moby::RCArticulatedBody() );
  pendulum->id = "pendulum";
  pendulum->algorithm_type = Moby::RCArticulatedBody::eCRB;

  // vectors for references used in defining the body
  std::vector< Moby::RigidBodyPtr > links;
  std::vector< Moby::JointPtr > joints;

  // create the static base;
  Moby::RigidBodyPtr base( new Moby::RigidBody() );
  {
    // create a primitive box for inertia and rendering prototypes
    Moby::PrimitivePtr box( new Moby::BoxPrimitive(0.1,0.1,0.1) );
    box->set_mass( 1 );

    // assign the static base parameters
    base->id = "base";                        // identify the link
    base->set_visualization_data( box->create_visualization() );  // attach a visualization for osg
    base->set_inertia( box->get_inertia() );  // use the inertia of the box as a model
    base->set_enabled( false );               // disable physics for the base (this makes it static)
 
    // compute the pose of the base 
    Ravelin::Quatd rotation(Ravelin::Quatd::normalize(Ravelin::Quatd(0,0,0,1)));
    Ravelin::Origin3d position(0,0,0);
    Ravelin::Pose3d pose( rotation, position );
    base->set_pose( pose ); 

    // add the base to the set of links
    links.push_back( base );
  }

  // create the dynamic arm link
  Moby::RigidBodyPtr arm( new Moby::RigidBody() );
  {
    // create a primitive cylinder for inertia and rendering prototypes
    Moby::PrimitivePtr cylinder( new Moby::CylinderPrimitive(0.025,1) );
    cylinder->set_mass( 1 );

    // assign the arm parameters
    arm->id = "arm";                              // identify the link
    arm->set_visualization_data( cylinder->create_visualization() );  // attach a visualization for osg
    arm->set_inertia( cylinder->get_inertia() );  // use the inertia of the cylinder as a model
    arm->set_enabled( true );                     // enable physics for the arm
    arm->get_recurrent_forces().push_back( g );   // add the gravity force to the arm
 
    // compute the pose of the arm 
    Ravelin::Quatd rotation(Ravelin::Quatd::normalize(Ravelin::Quatd(0,0,0,1)));
    Ravelin::Origin3d position(0,-0.5,0);
    Ravelin::Pose3d pose( rotation, position );
    arm->set_pose( pose ); 

    // add the arm to the set of links
    links.push_back( arm ); 
  }

  // create the pivot joint
  boost::shared_ptr<Moby::RevoluteJoint> pivot( new Moby::RevoluteJoint() );
  {
    // compute the position of the joint... center of the base w.r.t global frame
    Ravelin::Pose3d pose = *base->get_pose();
    Ravelin::Vector3d position(pose.x.x(), pose.x.y(), pose.x.z(), Moby::GLOBAL);

    // compute the axis of rotation
    Ravelin::Vector3d axis(1,0,0,Moby::GLOBAL);   // actuates around global x-axis

    // assign the pivot parameters
    // Note: set_location(...) must precede set_axis(...)
    pivot->id = "pivot";
    pivot->set_location( position, base, arm );
    pivot->set_axis( axis ); 

    // add the pivot to the set of joints
    joints.push_back( pivot );
  }

  // construct the pendulum articulated body from the set of links and joints created above
  pendulum->set_links_and_joints( links, joints );
  // add gravity as force to affect the pendulum
  pendulum->get_recurrent_forces().push_back( g );
  // pendulum has a fixed base
  pendulum->set_floating_base(false);
 
  // add the pendulum to the simulaiton
  sim->add_dynamic_body( pendulum );

  // create a viewer to visualize the simulation (Must have built with OSG support to see)
  Viewer viewer( sim, Ravelin::Origin3d(-5,0,-1), Ravelin::Origin3d(0,0,-1), Ravelin::Origin3d(0,0,1) );

  // start the main simulation loop
  while(true) {
    // update the viewer, if the viewer fails to update then it has been closed so exit
    if( !viewer.update() ) { 
      break;
    }
    // step the simulation forward one millisecond
    sim->step( 0.001 );

    // get an updated arm pose w.r.t the global frame and print it to the console 
    Ravelin::Pose3d pose = *arm->get_pose();
    pose.update_relative_pose(Moby::GLOBAL);
    std::cout << "t: " << sim->current_time << ", pose: " << pose << std::endl;
  }

  return 0;
}

