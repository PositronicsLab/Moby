#ifndef _VIEWER_H_
#define _VIEWER_H_

#include <boost/shared_ptr.hpp>
#include <Moby/Simulator.h>

#ifdef USE_OSG
#include <osgViewer/Viewer>
#include <osgViewer/ViewerEventHandlers>
#include <osg/Geode>
#include <osgDB/ReadFile>
#include <osgDB/WriteFile>
#include <osgGA/TrackballManipulator>
#include <osgGA/StateSetManipulator>
#endif

class Viewer {
public:

  Viewer( boost::shared_ptr<Moby::Simulator> sim, Ravelin::Origin3d camera_position, Ravelin::Origin3d camera_viewpoint, Ravelin::Origin3d camera_up ) {
    _sim = sim;

    #ifdef USE_OSG
    _viewer.setThreadingModel(osgViewer::Viewer::SingleThreaded);
    _viewer_pointer = &_viewer;

    // setup any necessary handlers here
    _viewer.setCameraManipulator(new osgGA::TrackballManipulator());
    _viewer.addEventHandler(new osgGA::StateSetManipulator(_viewer.getCamera()->getOrCreateStateSet()));
    _viewer.addEventHandler(new osgViewer::WindowSizeHandler);
    _viewer.addEventHandler(new osgViewer::StatsHandler);
 
   // setup the window
    _viewer_pointer->setUpViewInWindow(0, 0, 1024, 1024);

    osg::Vec3d osg_eye(camera_position[0],camera_position[1],camera_position[2]);
    osg::Vec3d osg_center(camera_viewpoint[0], camera_viewpoint[1], camera_viewpoint[2]);
    osg::Vec3d osg_up(camera_up[0], camera_up[1], camera_up[2]);

    // set the camera view
    if(_viewer_pointer && _viewer_pointer->getCameraManipulator()) {
      osg::Camera* camera = _viewer_pointer->getCamera();
      camera->setViewMatrixAsLookAt(osg_eye, osg_center, osg_up);

       // setup the manipulator using the camera, if necessary
      _viewer_pointer->getCameraManipulator()->setHomePosition(osg_eye, osg_center, osg_up);
    }


    _MAIN_GROUP = new osg::Group;

    _MAIN_GROUP->addChild(sim->get_persistent_vdata());
    _MAIN_GROUP->addChild(sim->get_transient_vdata());

    _viewer.setSceneData(_MAIN_GROUP);
    _viewer.realize();
    #endif // USE_OSG
  }

  virtual ~Viewer( void ) {

  }

  bool update( void ) {
    _sim->update_visualization();

    #ifdef USE_OSG
    if( _viewer.done() ) return false;
    _viewer.frame();
    #endif // USE_OSG
    return true;
  }

private:
  #ifdef USE_OSG
  /// The OpenInventor root group node for this application
  osg::Group* _MAIN_GROUP;
  /// Pointer to the viewer
  osgViewer::Viewer* _viewer_pointer;
  osgViewer::Viewer _viewer;
  #endif // USE_OSG

  boost::shared_ptr<Moby::Simulator> _sim;
};


#endif // _VIEWER_H_
