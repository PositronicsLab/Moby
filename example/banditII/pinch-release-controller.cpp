#include <cmath>
#include <sys/types.h>
#include <sys/times.h>
#include <fstream>
#include <Moby/XMLReader.h>
#include <Moby/Simulator.h>
#include <Moby/EventDrivenSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/DeformableBody.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/GravityForce.h>
#include <Moby/LinAlg.h>
#include <Moby/Constants.h>
#include <Moby/RNEAlgorithm.h>
#include <osgDB/WriteFile>

#include <iomanip>

using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

/*****************************************************************************

 ****************************************************************************/

/// plugin must be "extern C"

extern "C" {

int ITER;

// The Simulator
boost::shared_ptr<EventDrivenSimulator> sim;

// The ball
DeformableBodyPtr ball;

// The pincer(s)
ArticulatedBodyPtr pincers;
RigidBodyPtr left_block;
RigidBodyPtr right_block;

// The position of the ball in the pincers' coordinate systems
bool first = true;
Vector3 _ball_pincer_left;
Vector3 _ball_pincer_right;

const Real PINCH_FORCE_MAGNITUDE = 10.0;
const Real PINCH_DURATION = 0.25;

// pincer controller
void pincer_controller(DynamicBodyPtr robot, Real time, void* data)
{

  if( time < PINCH_DURATION ) {
    // close pincers
    Vector3 l_pinch_force = Vector3( -1.0, 0.0, 0.0 ) * PINCH_FORCE_MAGNITUDE;
    left_block->add_force( l_pinch_force );
    Vector3 r_pinch_force = Vector3( 1.0, 0.0, 0.0 ) * PINCH_FORCE_MAGNITUDE;
    right_block->add_force( r_pinch_force );
  } else {
    // open pincers
    Vector3 l_pinch_force = Vector3( 1.0, 0.0, 0.0 ) * PINCH_FORCE_MAGNITUDE;
    left_block->add_force( l_pinch_force );
    Vector3 r_pinch_force = Vector3( -1.0, 0.0, 0.0 ) * PINCH_FORCE_MAGNITUDE;
    right_block->add_force( r_pinch_force );
  }

}

static int step = 0;
static int event_iter = 0;

// Event Callback Function for EventDrivenSimulator
void event_callback_fn( std::vector<Event>& events, boost::shared_ptr<void> p ) {

  if( event_iter == ITER ) {
    step++;  
  } else {
    step = 1;
  }

  std::stringstream ss_state_trans;
  ss_state_trans << "state_trans_" << std::setfill('0') << std::setw(4) << ITER << "_" << std::setfill('0') << std::setw(4) << step << ".stl";
  std::string trans_state_file_name = ss_state_trans.str();

  osg::Node* node_trans = sim->get_transient_vdata();
  osgDB::writeNodeFile(*node_trans, trans_state_file_name);

  std::stringstream ss_state_persist;
  ss_state_persist << "state_persist_" << std::setfill('0') << std::setw(4) << ITER << "_" << std::setfill('0') << std::setw(4) << step << ".stl";
  std::string persist_state_file_name = ss_state_persist.str();

  osg::Node* node_persist = sim->get_persistent_vdata();
  osgDB::writeNodeFile(*node_persist, persist_state_file_name);

  event_iter = ITER;
}

void init(void* separator, const std::map<std::string, BasePtr>& read_map, Real time)
{
  // ---  Sim  ---
  if (read_map.find("simulator") == read_map.end())
    throw std::runtime_error("pinch-release-controller.cpp:init()- unable to find simulator!");
  // get a reference to the EventDrivenSimulator instance
  sim = dynamic_pointer_cast<EventDrivenSimulator>( read_map.find( "simulator" )->second );
  // register the event callback function
  sim->event_callback_fn = &event_callback_fn;

  // ---  Right Pincer  ---
  // find the right_pincer
  if (read_map.find("pincers") == read_map.end())
    throw std::runtime_error("pinch-release-controller.cpp:init()- unable to find pincers!");
  pincers = dynamic_pointer_cast<ArticulatedBody>(read_map.find("pincers")->second);
  if (!pincers)
    throw std::runtime_error("pinch-release-controller.cpp:init()- unable to cast pincers to type ArticulatedBody");

  // set up the pincer controller
  pincers->controller = pincer_controller;

  // ---  Ball  ---
  // find the ball
  if (read_map.find("ball") == read_map.end())
    throw std::runtime_error("pinch-release-controller.cpp:init()- unable to find ball!");
  ball = dynamic_pointer_cast<DeformableBody>(read_map.find("ball")->second);
  if (!ball)
    throw std::runtime_error("pinch-release-controller.cpp:init()- unable to cast ball to type DeformableBody");

  // setup the pincer rigid body links
  left_block = pincers->get_links()[1]; 
  right_block = pincers->get_links()[2]; 
}

} // end extern C

