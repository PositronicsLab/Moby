#ifndef _VIEWER_H_
#define _VIEWER_H_

#include <boost/shared_ptr.hpp>
#include <Moby/Simulator.h>

#ifdef USE_OSG
  #include <osgViewer/Viewer>
  #include <osgViewer/ViewerEventHandlers>
  #include <osgGA/TrackballManipulator>
  #include <osgGA/StateSetManipulator>
#endif

class Viewer {
public:

  Viewer( boost::shared_ptr<Moby::Simulator> sim, Ravelin::Origin3d camera_position, Ravelin::Origin3d camera_viewpoint, Ravelin::Origin3d camera_up ) : _sim(sim) {

    #ifdef USE_OSG
    // convert parameters to compatible osg types
    osg::Vec3d osg_eye(camera_position[0],camera_position[1],camera_position[2]);
    osg::Vec3d osg_center(camera_viewpoint[0], camera_viewpoint[1], camera_viewpoint[2]);
    osg::Vec3d osg_up(camera_up[0], camera_up[1], camera_up[2]);

    // initialize viewer
    _viewer.setThreadingModel(osgViewer::Viewer::SingleThreaded);

    // initialize handlers
    _viewer.setCameraManipulator(new osgGA::TrackballManipulator());
    _viewer.addEventHandler(new osgGA::StateSetManipulator(_viewer.getCamera()->getOrCreateStateSet()));
    _viewer.addEventHandler(new osgViewer::WindowSizeHandler);
    _viewer.addEventHandler(new osgViewer::StatsHandler);
 
    // initialize window
    _viewer.setUpViewInWindow(0, 0, 1024, 1024);

    // initialize camera
    osg::Camera* camera = _viewer.getCamera();
    camera->setViewMatrixAsLookAt(osg_eye, osg_center, osg_up);
    _viewer.getCameraManipulator()->setHomePosition(osg_eye, osg_center, osg_up);

    // establish the scene data
    _root = new osg::Group;
    _root->addChild(sim->get_persistent_vdata());
    _root->addChild(sim->get_transient_vdata());
    _viewer.setSceneData(_root);
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
  osg::Group* _root;            // the osg scene graph root node
  osgViewer::Viewer _viewer;    // the osg viewer
  #endif // USE_OSG

  boost::shared_ptr<Moby::Simulator> _sim;   // the simulator
};


#endif // _VIEWER_H_
