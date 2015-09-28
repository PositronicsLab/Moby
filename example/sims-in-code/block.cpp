#include <iostream>
#include <boost/shared_ptr.hpp>
#include <Moby/Simulator.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/EulerIntegrator.h>
#include <Moby/GravityForce.h>
#include <Moby/Log.h>

#include "viewer.h"

int main( void ) {

  boost::shared_ptr<Moby::Simulator> sim( new Moby::Simulator() );
  sim->integrator = boost::shared_ptr<Moby::Integrator>( new Moby::EulerIntegrator() );

  boost::shared_ptr<Moby::GravityForce> g( new Moby::GravityForce() );
  g->gravity = Ravelin::Vector3d( 0, 0, -9.8 );

  Moby::PrimitivePtr box( new Moby::BoxPrimitive(1,1,1) );
  box->set_mass( 1000 );

  Moby::RigidBodyPtr rb( new Moby::RigidBody() );
  rb->id = "block";
  rb->set_visualization_data( box->create_visualization() );
  rb->set_inertia( box->get_inertia() );
  rb->set_enabled( true );
  rb->get_recurrent_forces().push_back( g );

  rb->set_pose( Ravelin::Pose3d( Ravelin::Quatd(0,0,0,1), Ravelin::Origin3d(0,0,0) ) );

  sim->add_dynamic_body( rb );

  Viewer viewer( sim, Ravelin::Origin3d(0,0,4), Ravelin::Origin3d(0,0,0), Ravelin::Origin3d(1,0,0) );  

  while(true) {
    if( !viewer.update() ) break;

    sim->step( 0.001 );
    boost::shared_ptr<const Ravelin::Pose3d> pose = rb->get_pose();
    std::cout << "t: " << sim->current_time << " x: " << pose->x << std::endl;
  }

  return 0;
}

