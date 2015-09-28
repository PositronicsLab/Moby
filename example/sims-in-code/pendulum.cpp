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

//#define USE_OSG

#ifdef USE_OSG
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Geode>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgGA/TrackballManipulator>
#include <osgGA/StateSetManipulator>
#endif


#ifdef USE_OSG
/// The OpenInventor group node for Moby
osg::Group* MOBY_GROUP;
/// The OpenInventor root group node for this application
osg::Group* MAIN_GROUP;
/// Pointer to the viewer
osgViewer::Viewer* viewer_pointer;
#endif

int main( void ) {
  // setup osg vectors
  #ifdef USE_OSG
  osgViewer::Viewer viewer;
  viewer.setThreadingModel(osgViewer::Viewer::SingleThreaded);
  viewer_pointer = &viewer;

  // setup the window
  viewer_pointer->setUpViewInWindow(0, 0, 1024, 1024);

  osg::Vec3d position_osg(-4,0,-1);
  osg::Vec3d target_osg(0, 0, -1);
  osg::Vec3d up_osg(0, 0, 1);
  // set the camera view
  if (viewer_pointer && viewer_pointer->getCameraManipulator())
  {
  osg::Camera* camera = viewer_pointer->getCamera();
  camera->setViewMatrixAsLookAt(position_osg, target_osg, up_osg);
  // setup the manipulator using the camera, if necessary
  viewer_pointer->getCameraManipulator()->setHomePosition(position_osg, target_osg, up_osg);
  }
  #endif

  boost::shared_ptr<Moby::Simulator> sim( new Moby::Simulator() );
  sim->integrator = boost::shared_ptr<Moby::Integrator>( new Moby::EulerIntegrator() );

  boost::shared_ptr<Moby::GravityForce> g( new Moby::GravityForce() );
  g->gravity = Ravelin::Vector3d( 0, 0, -9.8 );

  Moby::RCArticulatedBodyPtr ab( new Moby::RCArticulatedBody() );
  ab->id = "pendulum";

  std::vector< Moby::RigidBodyPtr > links;
  Moby::RigidBodyPtr base( new Moby::RigidBody() );
  {
    Moby::PrimitivePtr box( new Moby::BoxPrimitive(1,1,1) );
    box->set_mass( 1000 );

    // static
    base->id = "base";
    base->set_visualization_data( box->create_visualization() );
    base->set_inertia( box->get_inertia() );
    base->set_enabled( true );
    //base->get_recurrent_forces().push_back( g );
  
    base->set_pose( Ravelin::Pose3d( Ravelin::Quatd(0,0,0,1), Ravelin::Origin3d(0,0,0) ) ); 
    links.push_back( base );
  }

  Moby::RigidBodyPtr link( new Moby::RigidBody() );
  {
    Moby::PrimitivePtr cylinder( new Moby::CylinderPrimitive(0.025,1) );
    cylinder->set_mass( 1 );

    link->id = "link";
    link->set_visualization_data( cylinder->create_visualization() );
    link->set_inertia( cylinder->get_inertia() );
    link->set_enabled( true );
    link->get_recurrent_forces().push_back( g );
  
    link->set_pose( Ravelin::Pose3d( Ravelin::Quatd(0,0,0,1), Ravelin::Origin3d(0,0,-0.5) ) ); 
    links.push_back( link ); 
  }

  std::vector< Moby::JointPtr > joints;
  boost::shared_ptr<Moby::RevoluteJoint> joint( new Moby::RevoluteJoint() );
  {
    joint->id = "joint";
    joint->set_location( Ravelin::Vector3d(0,0,0,Moby::GLOBAL), base, link );
    joint->set_axis( Ravelin::Vector3d(1,0,0,Moby::GLOBAL) );

    joints.push_back( joint );
  }

  ab->set_links_and_joints( links, joints );

  sim->add_dynamic_body( ab );

/*
  boost::shared_ptr<Ravelin::Pose3d> impulse_pos( new Ravelin::Pose3d( Ravelin::Quatd::normalize(Ravelin::Quatd(0,1,0,1)), Ravelin::Origin3d(0,0,-0.5)) );
  Ravelin::SForced impulse(0,100,0,0,0,0,impulse_pos);
  link->add_force( impulse );
*/
  for( unsigned i = 0; i < 10; i++ ) {
    sim->step( 0.001 );
    boost::shared_ptr<const Ravelin::Pose3d> pose = link->get_pose();
    std::cout << "step: " << i << " x: " << pose->x << std::endl;

    sim->update_visualization();
  }

  return 0;
}

