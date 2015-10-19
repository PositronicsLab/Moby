#include <iostream>
#include <boost/shared_ptr.hpp>
#include <Moby/Simulator.h>
#include <Moby/EulerIntegrator.h>
#include <Moby/GravityForce.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/PrismaticJoint.h>
#include <Moby/Log.h>

#include "viewer.h"

int main( void ) {

// uncomment to log dynamics
  Moby::Log<Moby::OutputToFile>::reporting_level = 7;

  boost::shared_ptr<Moby::Simulator> sim( new Moby::Simulator() );
  sim->integrator = boost::shared_ptr<Moby::Integrator>( new Moby::EulerIntegrator() );

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
  
    //link->set_pose( Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(1,0,0,1)), Ravelin::Origin3d(0,0,-0.5) ) ); 
    link->set_pose( Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(0,0,0,1)), Ravelin::Origin3d(0,0,-0.275) ) ); 
    links.push_back( link ); 
  }

  std::vector< Moby::JointPtr > joints;
  boost::shared_ptr<Moby::PrismaticJoint> joint( new Moby::PrismaticJoint() );
  {
    joint->id = "joint";
    joint->set_location( Ravelin::Vector3d(0,0,0,base->get_pose()), base, link );
    joint->set_axis( Ravelin::Vector3d(1,0,0,Moby::GLOBAL) );
    joint->set_constraint_type(Moby::Joint::eExplicit);
    //joint->lolimit = -0.5;
    //joint->hilimit = 0.5;

    std::cout << "joint: { lolimit[" << joint->lolimit << "], hilimit[" << joint->hilimit << "], maxforce[" << joint->maxforce <<"] }" << std::endl;
    //joint->maxforce = Ravelin::VectorNd(1);
    //maxforce[0]
    
    joints.push_back( joint );
  }

  ab->set_links_and_joints( links, joints );
  ab->get_recurrent_forces().push_back( g );
  ab->set_floating_base(false);

  sim->add_dynamic_body( ab );

  Viewer viewer( sim, Ravelin::Origin3d(-5,0,-1), Ravelin::Origin3d(0,0,-1), Ravelin::Origin3d(0,0,1) );

/*
  boost::shared_ptr<Ravelin::Pose3d> impulse_pos( new Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(1,0,0,1)), Ravelin::Origin3d(0,0,-0.5)) );
  Ravelin::SForced impulse(100,0,0,0,0,0,impulse_pos);
  link->add_force( impulse );
*/

  Ravelin::VectorNd f(1);
  f[0] = 1000;
  joint->add_force(f);

  while(true) {
    if( !viewer.update() ) break;

    joint->add_force(f);

    sim->step( 0.001 );
    Ravelin::Pose3d pose = *link->get_pose();
    pose.update_relative_pose(Moby::GLOBAL);
    std::cout << "t: " << sim->current_time << " x: " << pose.x << std::endl;

  }

  return 0;
}

