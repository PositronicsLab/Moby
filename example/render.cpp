// viewer.cpp
#ifdef USE_OSG
#include <osgDB/WriteFile>
#include <osg/CameraNode>
#include <osgViewer/Viewer>
#include <osg/Node>
#include <osg/PositionAttitudeTransform>
#include <osgDB/ReadFile>
#include <osgGA/TrackballManipulator>
#include <osgGA/StandardManipulator>
#include <osgGA/StateSetManipulator>
#include <osg/MatrixTransform>
#include <osg/LightSource>
#include <osg/Light>
#include <osg/LightModel>
#include <osg/Shader>
#include <osg/ShadowVolumeOccluder>
#include <osg/ShadeModel>
#include <osg/ShaderAttribute>
#include <osg/ShaderComposer>
#include <osg/ShapeDrawable>

std::string sceneFile;
int WIDTH, HEIGHT;

class SnapImageDrawCallback : public ::osg::CameraNode::DrawCallback
{
       public:
       SnapImageDrawCallback()
       {
       _snapImageOnNextFrame = false;
       }

       void setFileName(const std::string& filename) { _filename = filename; }
       const std::string& getFileName() const { return _filename; }
       void setSnapImageOnNextFrame(bool flag) { _snapImageOnNextFrame = flag; }
       bool getSnapImageOnNextFrame() const { return _snapImageOnNextFrame; }
       virtual void operator () (const ::osg::CameraNode& camera) const
       {
         if (!_snapImageOnNextFrame) return;
//         ::osg::notify(::osg::NOTICE) << "Saving screen image to '"<<_filename<<"'"<< std::endl;
         int x,y,width,height;
         x = camera.getViewport()->x();
         y = camera.getViewport()->y();
         width = (WIDTH != 0)? WIDTH : camera.getViewport()->width();
         height = (HEIGHT != 0)? HEIGHT : camera.getViewport()->height();
//         ::osg::notify(::osg::NOTICE) << "Capturing image from: (" << x << ", " << y<< ")    " <<width<< " x "<< height << std::endl;

         ::osg::ref_ptr< ::osg::Image> image = new ::osg::Image();
         image->readPixels(x,y,width,height,GL_RGB,GL_UNSIGNED_BYTE);

         if (::osgDB::writeImageFile(*image,_filename)){
//           ::osg::notify(::osg::NOTICE) << "Saved screen image to '"<<_filename<<"'"<< std::endl;
         } else {
           ::osg::notify(::osg::NOTICE) << "Image not saved"<< std::endl;
         }
         _snapImageOnNextFrame = false;
       }

       protected:

       ::std::string _filename;
       mutable bool _snapImageOnNextFrame;

};

void
renderSceneToImage(::osg::Node* node, const ::std::string& sFileName_,double position[3],double target[3],double up[3])
{
  osg::Group* root = new osg::Group();

  // Declare transform, initialize with defaults.

  osg::PositionAttitudeTransform* nodeXform =
     new osg::PositionAttitudeTransform();

  // Use the 'addChild' method of the osg::Group class to
  // add the transform as a child of the root node and the
  // node node as a child of the transform.

  root->addChild(nodeXform);

  nodeXform->addChild(node);

  if(!sceneFile.empty()){
    ::osg::Node* sceneNode = osgDB::readNodeFile(sceneFile);
    nodeXform->addChild(sceneNode);
  }

  // Declare and initialize a Vec3 instance to change the
  // position of the node model in the scene
  osg::Vec3 nodePosit(0,0,0);
  nodeXform->setPosition( nodePosit );

  // Declare a 'viewer'
  osgViewer::Viewer viewer;

  // Next we will need to assign the scene graph we created
  // above to this viewer:
  viewer.setSceneData( root );

  viewer.setCameraManipulator(new osgGA::TrackballManipulator());
  viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));

  osg::Vec3d position_osg(position[0], position[1], position[2]);
  osg::Vec3d target_osg(target[0], target[1], target[2]);
  osg::Vec3d up_osg(up[0], up[1], up[2]);

  osg::Vec3d view = target_osg - position_osg;

  // compute the up and normal vectors 
  //osg::Quat rot;
  // compute the rotation from the view vector to the world up
  //rot.makeRotate( up_osg, view );    // #unused
  // find the normal vector by crossing the view and world up vectors
  osg::Vec3d n = view^up_osg;
  // find desired up vector by crossing the normal vector with the view vector
  osg::Vec3d up_desired = n^view;

  //osg::Vec3d up_new = rot * up_osg;  // #unused

  // replace the up vector with the desired up
  up_osg = up_desired;

  // set the camera view

  osg::Camera* camera = viewer.getCamera();
  camera->setViewMatrixAsLookAt(position_osg, target_osg, up_osg);

  // setup the manipulator using the camera, if necessary
  viewer.getCameraManipulator()->setHomePosition(position_osg, target_osg, up_osg);

  ::osg::ref_ptr<SnapImageDrawCallback> snapImageDrawCallback = new SnapImageDrawCallback();
  camera->setPostDrawCallback (snapImageDrawCallback.get());

  snapImageDrawCallback->setFileName(sFileName_);
  snapImageDrawCallback->setSnapImageOnNextFrame(true);

  // Add a Light to the scene
  osg::ref_ptr<osg::Group> lightGroup (new osg::Group);
  osg::ref_ptr<osg::StateSet> lightSS (root->getOrCreateStateSet());
  osg::ref_ptr<osg::LightSource> lightSource1 = new osg::LightSource;

  double xCenter = 10, yCenter=10;
  osg::Vec4f lightPosition (osg::Vec4f(xCenter, yCenter,75,1.0f));
  osg::ref_ptr<osg::Light> light = new osg::Light;
  light->setLightNum(1);
  light->setPosition(lightPosition);
  light->setAmbient(osg::Vec4(0.3f,0.3f,0.3f,0.4f));
  light->setDiffuse(osg::Vec4(0.2f,0.2f,0.2f,0.5f));
//  light->setSpecular(osg::Vec4(0.1,0.1,0.1,0.3));
//  light->setConstantAttenuation(0.5f);
  light->setDirection(osg::Vec3(0.1f, 0.1f, -1.0f));

  lightSource1->setLight(light.get());

  lightSource1->setLocalStateSetModes(osg::StateAttribute::ON);
  lightSource1->setStateSetModes(*lightSS,osg::StateAttribute::ON);
  //osg::StateSet* lightSS (lightGroup->getOrCreateStateSet());

  lightGroup->addChild(lightSource1.get());

  //Light markers: small spheres
  osg::ref_ptr<osg::Geode> lightMarkerGeode (new osg::Geode);
  lightMarkerGeode->addDrawable(new osg::ShapeDrawable(new osg::Sphere(osg::Vec3f(xCenter,yCenter,75),10.0f)));


  //Tuto 9: lighting code
//  root->addChild(lightGroup.get());
  //Tuto 9: Adding the light marker geode
//  root->addChild(lightMarkerGeode.get());

  viewer.realize();

  int x,y,width,height;
  x = camera->getViewport()->x();
  y = camera->getViewport()->y();
         width = (WIDTH != 0)? WIDTH : camera->getViewport()->width();
         height = (HEIGHT != 0)? HEIGHT : camera->getViewport()->height();
//    ::osg::notify(::osg::NOTICE) << "Capturing image from: (" << x << ", " << y<< ")    " <<width<< " x "<< height << std::endl;

  // Prevent this from opening a window by making pbuffer context
//  osg::ref_ptr<osg::GraphicsContext::Traits> traits = new osg::GraphicsContext::Traits;
//  traits->x = 0;
//  traits->y = 0;
//  traits->width = width;
//  traits->height = height;
//  traits->red = 8;
//  traits->green = 8;
//  traits->blue = 8;
//  traits->alpha = 8;
//  traits->windowDecoration = false;
//  traits->pbuffer = true;
//  traits->doubleBuffer = true;
//  traits->sharedContext = 0;

//  osg::ref_ptr<osg::GraphicsContext> pbuffer;
//  pbuffer = ::osg::GraphicsContext::createGraphicsContext(traits.get());
//  if (pbuffer.valid())
//  {
//      ::osg::notify(osg::NOTICE)<<"Pixel buffer has been created successfully."<<std::endl;
//  }
//  else
//  {
//      ::osg::notify(osg::NOTICE)<<"Pixel buffer has not been created successfully."<<std::endl;
//  }

//  if (pbuffer.valid())
//  {
//      osg::ref_ptr<osg::Camera> camera = new osg::Camera;
//      camera->setGraphicsContext(pbuffer.get());
//      camera->setViewport(new osg::Viewport(0,0,width,height));
//      GLenum buffer = pbuffer->getTraits()->doubleBuffer ? GL_BACK : GL_FRONT;
//      camera->setDrawBuffer(buffer);
//      camera->setReadBuffer(buffer);
////      camera->setFinalDrawCallback(new WindowCaptureCallback(mode, position, readBuffer));
//      camera->setFinalDrawCallback(snapImageDrawCallback.get());

//      viewer.addSlave(camera.get(), osg::Matrixd(), osg::Matrixd());

      viewer.realize();

      viewer.frame();
//  }

}

void render(::osg::Node* node, const ::std::string& sFileName_){

  osg::Group* root = new osg::Group();

  // Declare transform, initialize with defaults.

  osg::PositionAttitudeTransform* nodeXform =
     new osg::PositionAttitudeTransform();

  // Use the 'addChild' method of the osg::Group class to
  // add the transform as a child of the root node and the
  // node node as a child of the transform.

  root->addChild(nodeXform);

  nodeXform->addChild(node);

  // Declare and initialize a Vec3 instance to change the
  // position of the node model in the scene
  osg::Vec3 nodePosit(0,0,0);
  nodeXform->setPosition( nodePosit );

  // Declare a 'viewer'
  osgViewer::Viewer viewer;

  // Next we will need to assign the scene graph we created
  // above to this viewer:
  viewer.setSceneData( root );

  // attach a trackball manipulator to all user control of the view
  viewer.setCameraManipulator(new osgGA::TrackballManipulator);

  // create the windows and start the required threads.
  viewer.realize();

  // Enter the simulation loop. viewer.done() returns false
  // until the user presses the 'esc' key.
  // (This can be changed by adding your own keyboard/mouse
  // event handler or by changing the settings of the default
  // keyboard/mouse event handler)

  while( !viewer.done() )
  {
     // dispatch the new frame, this wraps up the follow Viewer operations:
     //   advance() to the new frame
     //   eventTraversal() that collects events and passes them on to the event handlers and event callbacks
     //   updateTraversal() to calls the update callbacks
     //   renderingTraversals() that runs syncronizes all the rendering threads (if any) and dispatch cull, draw and swap buffers
     viewer.frame();
  }
}

int main(int argc, char** argv)
{
  ::osg::Node* pRoot;
  std::string sFileName;

  if (argc < 3)
  {
    std::cerr << "syntax: moby-render <nodefile> [options] <output file>" << std::endl;
    std::cerr << "  option -p x y z (position of camera)" << std::endl;
    std::cerr << "  option -t x y z (target position)" << std::endl;
    std::cerr << "  option -u x y z (up vector)" << std::endl;
    std::cerr << "  option -s filename (scene file)" << std::endl;
    std::cerr << "  option -w width height (rendering width and height)" << std::endl;
    return -1;
  }

  // get all options
  pRoot = osgDB::readNodeFile(argv[1]);
  sFileName = std::string(argv[argc-1]);
  double position[3] = {1,1,1};
  double target[3] = {0,0,0};
  double up[3] = {0,0,1};
  WIDTH = 0; HEIGHT = 0;
  for (int i=2; i< argc-1; i++)
  {
    // get the option
    std::string option(argv[i]);

    if (option.find("-p") != std::string::npos){
      position[0] = std::atof(&argv[i+1][0]);
      position[1] = std::atof(&argv[i+2][0]);
      position[2] = std::atof(&argv[i+3][0]);
    } else if (option.find("-t") != std::string::npos){
      target[0] = std::atof(&argv[i+1][0]);
      target[1] = std::atof(&argv[i+2][0]);
      target[2] = std::atof(&argv[i+3][0]);
    } else if (option.find("-u") != std::string::npos){
      up[0] = std::atof(&argv[i+1][0]);
      up[1] = std::atof(&argv[i+2][0]);
      up[2] = std::atof(&argv[i+3][0]);
    } else if (option.find("-s=") != std::string::npos){
      sceneFile = std::string(&argv[i][3]);
    } else if (option.find("-w") != std::string::npos){
      WIDTH = std::atoi(&argv[i+1][0]);
      HEIGHT = std::atoi(&argv[i+2][0]);
    }

  }
  ::osg::notify(::osg::NOTICE) << "Capturing image from: (" << position[0] << ", " << position[1]<<", " << position[2]<< ")  "
                                  "of object at  (" << target[0] << ", " << target[1]<<", " << target[2]<< ")" << std::endl;

/*
  // For tracking a specific target
  ::osg::Group* rt = pRoot->asGroup();
  unsigned count = rt->getNumChildren();
  printf( "root children: %u\n", count );
  for( unsigned i = 0; i < count; i++ ) {
    osg::Node* child = rt->getChild( i );
    printf( "child: %s", child->className() );
    osg::Group* childgroup = child->asGroup();
    if( childgroup != 0 ) {
      unsigned count2 = childgroup->getNumChildren();
      printf( " children: %u\n", count2 );
      for( unsigned j = 0; j < count2; j++ ) {
        osg::Node* grandchild = childgroup->getChild( i );
        printf( "  grandchild: %s\n", grandchild->className() );
        osg::MatrixTransform* mt = (osg::MatrixTransform*) grandchild;
        osg::Matrixd m = mt->getMatrix();
        std::cout << "  " <<  m(3,0) << "," << m(3,1) << "," << m(3,2) << std::endl;
        target[0] = m(3,0);
        target[1] = m(3,1);
        target[2] = m(3,2);
      }
    } else {
      printf( "\n" );
    }
  }
*/

  renderSceneToImage(pRoot,sFileName,position,target,up);
//  render(pRoot,sFileName);
}
#else
int main() {}
#endif

