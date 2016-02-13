#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
#include <unistd.h>
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Group>
#include <osg/Geode>
#include <osg/ShapeDrawable>
#include <osg/Material>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgGA/TrackballManipulator>
#include <osgGA/StateSetManipulator>

int frame_rate = 0;
int frame = 1;  // we start at 1 so that the last frame will be removed
std::vector<osg::Node*> nodes;
osg::Group* mainroot;

void cycle()
{
  osg::Node* last = (frame > 0) ? nodes[frame-1] : nodes.back();
  osg::Node* n = nodes[frame++];
  if (frame == nodes.size())
    frame = 0;
  mainroot->removeChild(last);
  mainroot->addChild(n);
}

int main(int argc, char ** argv)
{
  if (argc < 3) 
  {
    std::cerr << "syntax: view <frame rate> <input1> [input2] [input3] ... [inputn]" << std::endl << std::endl;
    std::cerr << "Views 3D files as an animation" << std::endl;
    std::cerr << " (if frame rate <= 0, then files will be output all at same time)" << std::endl;
    exit(1);
  }

  // get the frame rate
  frame_rate = std::atoi(argv[1]);

  // setup the osg Viewer
  osgViewer::Viewer viewer;
  viewer.setCameraManipulator(new osgGA::TrackballManipulator());
  viewer.setThreadingModel(osgViewer::Viewer::SingleThreaded);
  viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));
  viewer.addEventHandler(new osgViewer::WindowSizeHandler);
  viewer.addEventHandler(new osgViewer::StatsHandler);
  viewer.setLightingMode(osg::View::HEADLIGHT);

  // setup the view window size
  viewer.setUpViewInWindow(0, 0, 640, 480);

  // setup the camera
  viewer.getCamera()->setViewMatrixAsLookAt(osg::Vec3d(0,0,10), osg::Vec3d(0,0,0), osg::Vec3d(0,1,0));

  // create the group
  mainroot = new osg::Group;
  mainroot->ref();

  // attempt to open the scenery file
  std::string scenery_fname("scene.osg");
  std::ifstream in(scenery_fname.c_str());
  if (!in.fail())
  {
    std::cout << "...using scenery file" << std::endl;
    mainroot->addChild(osgDB::readNodeFile(scenery_fname));
    in.close();
  }

  // Open the argument file(s)..
  for (int i=2; i< argc; i++)
  {
    // do different things based on the file extension
    std::string fname(argv[i]);

    // read the file 
    nodes.push_back(osgDB::readNodeFile(fname));
    nodes.back()->ref();
  }
  mainroot->addChild(nodes.front());

  // go
  viewer.setSceneData(mainroot);
  viewer.realize();
  while (true)
  {
    if (viewer.done())
      break;
    viewer.frame();
    cycle();
    usleep(1.0/frame_rate);
  }

  return 0;
}

