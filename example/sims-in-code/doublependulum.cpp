#include <iostream>
#include <boost/shared_ptr.hpp>
#include <Moby/Simulator.h>
#include <Moby/EulerIntegrator.h>
#include <Moby/GravityForce.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/Log.h>

#include "viewer.h"

int main( void ) {

// uncomment to log dynamics
//  Moby::Log<Moby::OutputToFile>::reporting_level = 7;

  boost::shared_ptr<Moby::Simulator> sim( new Moby::Simulator() );
  sim->integrator = boost::shared_ptr<Moby::Integrator>( new Moby::EulerIntegrator() );

  boost::shared_ptr<Moby::GravityForce> g( new Moby::GravityForce() );
  g->gravity = Ravelin::Vector3d( 0, 0, -9.8 );

  Moby::RCArticulatedBodyPtr ab( new Moby::RCArticulatedBody() );
  ab->id = "doublependulum";
  ab->algorithm_type = Moby::RCArticulatedBody::eCRB;

  std::vector< Moby::RigidBodyPtr > links;
  std::vector< Moby::JointPtr > joints;

  Moby::RigidBodyPtr base( new Moby::RigidBody() );
  {
    Moby::PrimitivePtr box( new Moby::BoxPrimitive(0.1,0.1,0.1) );
    box->set_mass( 1000 );

    // static
    base->id = "base";
    base->set_visualization_data( box->create_visualization() );
    base->set_inertia( box->get_inertia() );
    base->set_enabled( false );
  
    base->set_pose( Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(0,0,0,1)), Ravelin::Origin3d(0,0,0) ) ); 
    links.push_back( base );
  }

  Moby::PrimitivePtr cylinder( new Moby::CylinderPrimitive(0.025,1) );
  cylinder->set_mass( 1 );

  Moby::RigidBodyPtr link1( new Moby::RigidBody() );
  {
    link1->id = "link1";
    link1->set_visualization_data( cylinder->create_visualization() );
    link1->set_inertia( cylinder->get_inertia() );
    link1->set_enabled( true );
    link1->get_recurrent_forces().push_back( g );
  
    //link1->set_pose( Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(1,0,0,1)), Ravelin::Origin3d(0,0,-0.5) ) ); 
    link1->set_pose( Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(0,0,0,1)), Ravelin::Origin3d(0,-0.5,0) ) ); 
    links.push_back( link1 ); 
  }

  boost::shared_ptr<Moby::RevoluteJoint> joint1( new Moby::RevoluteJoint() );
  {
    joint1->id = "joint1";
    joint1->set_location( Ravelin::Vector3d(0,0,0,base->get_inertial_pose()), base, link1 );
    joint1->set_axis( Ravelin::Vector3d(1,0,0,Moby::GLOBAL) );

    joints.push_back( joint1 );
  }

  Moby::RigidBodyPtr link2( new Moby::RigidBody() );
  {
    link2->id = "link2";
    link2->set_visualization_data( cylinder->create_visualization() );
    link2->set_inertia( cylinder->get_inertia() );
    link2->set_enabled( true );
    link2->get_recurrent_forces().push_back( g );
  
    //link2->set_pose( Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(1,0,0,1)), Ravelin::Origin3d(0,0,-1.5) ) ); 
    link2->set_pose( Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(0,0,0,1)), Ravelin::Origin3d(0,-1.5,0) ) ); 
    links.push_back( link2 ); 
  }

  boost::shared_ptr<Moby::RevoluteJoint> joint2( new Moby::RevoluteJoint() );
  {
    joint2->id = "joint2";
    joint2->set_location( Ravelin::Vector3d(0,-0.5,0,link1->get_inertial_pose()), link1, link2 );
    joint2->set_axis( Ravelin::Vector3d(1,0,0,Moby::GLOBAL) );

    joints.push_back( joint2 );
  }

  ab->set_links_and_joints( links, joints );
  ab->get_recurrent_forces().push_back( g );
  ab->set_floating_base(false);

  sim->add_dynamic_body( ab );

  Viewer viewer( sim, Ravelin::Origin3d(-10,0,-1), Ravelin::Origin3d(0,0,-1), Ravelin::Origin3d(0,0,1) );

/*
  boost::shared_ptr<Ravelin::Pose3d> impulse_pos( new Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(0,1,0,1)), Ravelin::Origin3d(0,0,-0.5)) );
  Ravelin::SForced impulse(0,100,0,0,0,0,impulse_pos);
  link->add_force( impulse );
*/
  
  while(true) {
    if( !viewer.update() ) break;

    sim->step( 0.001 );
    Ravelin::Pose3d pose = *link2->get_pose();
    pose.update_relative_pose(Moby::GLOBAL);
    std::cout << "t: " << sim->current_time << " x: " << pose.x << std::endl;

  }

  return 0;
}

