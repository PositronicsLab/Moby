#include <Moby/XMLReader.h>
#include <Moby/InventorOutput.h>
#include <Moby/Simulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/RCArticulatedBodyJoint.h>
#include <Moby/GravityForce.h>
#include <Moby/RNEAlgorithm.h>
#include <Moby/VRMLWriter.h>
#include <Inventor/Qt/SoQt.h>
#include <Inventor/Qt/viewers/SoQtExaminerViewer.h>
#include <Inventor/So.h>
#include <iostream>

// test for reading simulation from XML (FS Articulated bodies), running it,
// and displaying the results in VRML

using namespace Moby;

boost::shared_ptr<InventorOutput> viz;
const double STEP_SIZE = 0.01;
unsigned iteration = 0;
RNEAlgorithm rne;

/// Runs the simulator and updates the transforms
void render(void* arg, SoSensor* sensor)
{
	std::cout << "beginning iteration: " << iteration << std::endl;

	// get the pointer to the simulator
	boost::shared_ptr<Simulator> s = *(boost::shared_ptr<Simulator>*) arg;

	// get the articulated body and the links
	boost::shared_ptr<DynamicBody> db = s->get_dynamic_bodies().back();
	boost::shared_ptr<RCArticulatedBody> rcab = boost::dynamic_pointer_cast<RCArticulatedBody>(db);
	const std::vector<boost::shared_ptr<RCArticulatedBodyLink> >& links = rcab->get_links();

	// setup data for inverse dynamics
		std::map<boost::shared_ptr<RCArticulatedBodyLink>, RCArticulatedBodyInvDynData> inv_dyn_data;
	for (unsigned i=0; i< links.size(); i++)
	{
		RCArticulatedBodyInvDynData iddata;
		iddata._fext =  Vector3(0,-9.81,0) * links[i]->get_mass();
		iddata._text = ZEROS_3;
		iddata._qdd = VectorN::zero(1);
		inv_dyn_data[links[i]] = iddata;
	}

	// calculate inverse dynamics
	std::map<boost::shared_ptr<RCArticulatedBodyJoint>, VectorN> act_forces = rne.calc_inv_dyn(rcab, inv_dyn_data);
 
	// add joint forces
	for (unsigned i=1; i< links.size(); i++)
	{
		boost::shared_ptr<RCArticulatedBodyJoint> joint = links[i]->get_inner_joint();
		joint->add_force(act_forces[joint]);
	}

	// step the simulator
	s->step(STEP_SIZE);
	viz->update();
	iteration++;

	// output the scene
	char fname[30];
	sprintf(fname,"chain%d.wrl", iteration);
//	VRMLWriter::write_scene(std::string(fname), s);
}

int main(int argc, char** argv)
{
	// verify that argument is given
	if (argc < 2)
	{
		std::cerr << "syntax: chain <XML file>" << std::endl;
		return -1;
	}
				
	// setup the simulation
	XMLReader::read(std::string(argv[1]));

	// get the (only) simulation object
	boost::shared_ptr<Simulator> s = XMLReader::_sim_objs.front();

	// setup the visualization
	QWidget* mainwin = SoQt::init(argc, argv, argv[0]);
	SoQtExaminerViewer* viewer = new SoQtExaminerViewer(mainwin);
	SoSeparator* main_sep = new SoSeparator;
	main_sep->ref();

	// add a camera
	SoPerspectiveCamera* camera = new SoPerspectiveCamera;
	camera->position = SbVec3f(0, 0, 150);
	camera->pointAt(SbVec3f(0,0,1));
	main_sep->addChild(camera);

	// add lights
	SoDirectionalLight* light = new SoDirectionalLight;
	light->direction = SbVec3f(0,0,1);
	SoDirectionalLight* light2 = new SoDirectionalLight;
	light2->direction = SbVec3f(0,0,-1);
	SoDirectionalLight* light3 = new SoDirectionalLight;
	light3->direction = SbVec3f(-10,-5,-1);
	SoDirectionalLight* light4 = new SoDirectionalLight;
	light4->direction = SbVec3f(-10,-5,1);
	SoDirectionalLight* light5 = new SoDirectionalLight;
	light5->direction = SbVec3f(0,-10,1);
	main_sep->addChild(light);
	main_sep->addChild(light2);
	main_sep->addChild(light3);
	main_sep->addChild(light4);
	main_sep->addChild(light5);
	
	// setup the simulator visualization
	viz = boost::shared_ptr<InventorOutput>(new InventorOutput(s));
	SoSeparator* sep = viz->get_root();
	main_sep->addChild(sep);

	// start rendering
	viewer->setSceneGraph(sep);
	viewer->show();

	// add a timer sensor that runs the simulation
	SoTimerSensor* timer = new SoTimerSensor(&render, &s);
	SbTime interval(STEP_SIZE);
	timer->setInterval(interval);
	timer->schedule();

	// popup the main window
	SoQt::show(mainwin);
	SoQt::mainLoop();
	delete viewer;
}

