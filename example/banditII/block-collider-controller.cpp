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

using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

/*****************************************************************************
  Demonstrates/Validates applying a force to a block and pushing the block
  into another block resulting in the two blocks moving due to the force
  applied on the active block.

  To be used in conjunction with pinch_ball/block-collider.xml
 ****************************************************************************/

/// plugin must be "extern C"

extern "C" {

int ITER;

// The Simulator
boost::shared_ptr<EventDrivenSimulator> sim;

RigidBodyPtr left_block;

/// left_pincer control loop
void left_pincer_controller(DynamicBodyPtr robot, Real time, void* data)
{
  left_block->add_force( Vector3( -100.0, 0.0, 0.0 ) );
}

static int step = 0;

// callback for EventDrivenSimulator
void event_callback_fn( std::vector<Event>& events, boost::shared_ptr<void> p ) {

  std::stringstream ss_state_trans;
  ss_state_trans << "state_trans_" << ITER << "_" << step << ".stl";
  std::string trans_state_file_name = ss_state_trans.str();

  osg::Node* node_trans = sim->get_transient_vdata();
  osgDB::writeNodeFile(*node_trans, trans_state_file_name);

  std::stringstream ss_state_persist;
  ss_state_persist << "state_persist_" << ITER << "_" << step << ".stl";
  std::string persist_state_file_name = ss_state_persist.str();

  osg::Node* node_persist = sim->get_persistent_vdata();
  osgDB::writeNodeFile(*node_persist, persist_state_file_name);

  step++;  
}

void init(void* separator, const std::map<std::string, BasePtr>& read_map, Real time)
{
  // ---  Sim  ---
  if (read_map.find("simulator") == read_map.end())
    throw std::runtime_error("block-collider-controller.cpp:init()- unable to find simulator!");
  // get a reference to the EventDrivenSimulator instance
  sim = dynamic_pointer_cast<EventDrivenSimulator>( read_map.find( "simulator" )->second );
  // register the event callback function
  sim->event_callback_fn = &event_callback_fn;

  // ---  Left Pincer  ---
  // find the left_pincer
  if (read_map.find("left_pincer") == read_map.end())
    throw std::runtime_error("block-collider-controller.cpp:init()- unable to find left_pincer!");
  left_block = dynamic_pointer_cast<RigidBody>(read_map.find("left_pincer")->second);
  if (!left_block)
    throw std::runtime_error("block-collider-controller.cpp:init()- unable to cast left_pincer to type RigidBody");

  // set up the left_pincer controller
  DynamicBodyPtr left_pincer = dynamic_pointer_cast<DynamicBody>( left_block );
  if (!left_pincer)
    throw std::runtime_error("block-collider-controller.cpp:init()- unable to cast left_block to type DynamicBody");
  left_pincer->controller = left_pincer_controller;
}

} // end extern C

