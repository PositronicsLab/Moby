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

 ****************************************************************************/

/// plugin must be "extern C"

extern "C" {

int ITER;

// The Simulator
boost::shared_ptr<EventDrivenSimulator> sim;

// The ball
DeformableBodyPtr ball;

// ball_push_force_controller constants
const Real PUSH_FORCE_MAGNITUDE = 0.001;	// How strong a push
const int FIRST_PUSH_AT_ITERATION = 20;		// When to start pushing in iterations
const Real PUSH_DURATION = 0.01;		// How long to push for in time
const Real REST_DURATION = 0.01;		// How long to rest for in time
const unsigned int MAX_PUSHES = 2;		// The maximum number of times to push

// ball_push_force_controller variables and initial conditions
bool is_push_mode = false;			// is_push_mode == true -> push else rest
Real push_timer_start_time = 0.0;		// When to begin a mode
Real push_timer_end_time = 0.0;			// When to end a mode
Vector3 push_force = Vector3( 0.0, 0.0, 0.0 );	// The push force
Real push_force_time = 0.0;			// The time the controller last executed
unsigned int push_count = 0;			// The number of times a push has been applied

void ball_push_force_controller(DynamicBodyPtr robot, Real time, void* data)
{
  // Delay adding force until the pincer is in proper configuration
  if( ITER < FIRST_PUSH_AT_ITERATION ) return;

  if( time > push_timer_end_time ) {
    // change mode
    is_push_mode = !is_push_mode;

    // recalculate the timer interval
    push_timer_start_time = time;
    if( is_push_mode )
      push_timer_end_time = push_timer_start_time + PUSH_DURATION;
    else
      push_timer_end_time = push_timer_start_time + REST_DURATION;

    if( is_push_mode ) {
      // generate a new push force
/*
      // Generate a force either parallel to pincers (pushes through pincers)
      Real x = 0.0;
      Real y = 0.0;
      Real z = (Real)rand() / ( (Real)RAND_MAX / 2.0 ) - 1.0;	// [-1.0, 1.0]
*//*
      // Generate a force in a random direction of known magnitude
      Real x = (Real)rand() / ( (Real)RAND_MAX / 2.0 ) - 1.0;	// [-1.0, 1.0]
      Real y = (Real)rand() / ( (Real)RAND_MAX / 2.0 ) - 1.0;	// [-1.0, 1.0]
      Real z = (Real)rand() / ( (Real)RAND_MAX / 2.0 ) - 1.0;	// [-1.0, 1.0]
*/
///*
      // Generate a force along the z-axis (toward the viewer in default xml configurations)
      Real x = 0.0;
      Real y = 0.0;
      Real z = -1.0;
//*/
      push_force = Vector3( x, y, z );

      push_force.normalize();

      // Exceedingly unlikely, but throw anyway if arises
      if( push_force[0] == 0.0 && push_force[1] == 0.0 && push_force[2] == 0.0 )
        throw std::runtime_error("push-controller.cpp:ball_force_controller()- generated a zero push force!");

      push_force *= PUSH_FORCE_MAGNITUDE;
      push_count++;
    }
  }

  if( push_count > MAX_PUSHES ) return;

  if( is_push_mode ) {
    // -- the controller is pushing --
    // Apply the force
    // Note: only do this once per time slice so make sure that time has actually advanced
    if( time > push_force_time ) {
      //std::cout << "applying push: " << push_force[0] << ", " << push_force[1] << ", " << push_force[2] << "\n";
      ball->add_force( push_force );
    }
  }

  // update the time the controller last executed
  push_force_time = time;
}

static int step = 0;

// Event Callback Function for EventDrivenSimulator
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
    throw std::runtime_error("push-controller.cpp:init()- unable to find simulator!");
  // get a reference to the EventDrivenSimulator instance
  sim = dynamic_pointer_cast<EventDrivenSimulator>( read_map.find( "simulator" )->second );
  // register the event callback function
  sim->event_callback_fn = &event_callback_fn;

  // ---  Ball  ---
  // find the ball
  if (read_map.find("ball") == read_map.end())
    throw std::runtime_error("push-controller.cpp:init()- unable to find ball!");
  ball = dynamic_pointer_cast<DeformableBody>(read_map.find("ball")->second);
  if (!ball)
    throw std::runtime_error("push-controller.cpp:init()- unable to cast ball to type DeformableBody");

  DynamicBodyPtr ball_pusher = dynamic_pointer_cast<DynamicBody>(ball);
  if (!ball_pusher)
    throw std::runtime_error("push-controller.cpp:init()- unable to cast ball to type DynamicBody");
  // setup the push controller
  ball_pusher->controller = ball_push_force_controller;

}

} // end extern C

