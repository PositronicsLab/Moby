#include <cstring>
#include <vector>
#include <iostream>
#include <fstream>
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

// reads the OBJ file and creates a separator
osg::Group* create_group_from_obj(const char* filename)
{
  const unsigned BUF_SIZE = 2048;
  char buffer[BUF_SIZE];
  double v1, v2, v3;
  int i1, i2, i3;
  std::vector<double*> vertices;
  std::vector<int*> faces;
  const unsigned A = 0, B = 1, C = 2, X = 0, Y = 1, Z = 2;

  // open the filename for reading
  std::ifstream in(filename);

  // read in the file
  while (true)
  {
    // read in the line identifier
    std::string id;
    in >> id;
    if (in.eof())
      break;

    // determine whether the line describes a vertex
    if (strcasecmp(id.c_str(),"v") == 0)
    {
      // read in the vertex and the rest of the line
      in >> v1;
      in >> v2;
      in >> v3;
      if (in.eof())
        return false;

      // create the vertex
      double* vertex = new double[3];
      vertex[0] = v1;
      vertex[1] = v2;
      vertex[2] = v3;

      // store it in the linked list
      vertices.push_back(vertex);
      continue;
    }

    // determine whether the read line describes a face
    if (strcasecmp(id.c_str(),"f") == 0)
    {
      std::string f1, f2, f3;
      size_t idx;

      // read in the indices
      in >> f1;
      in >> f2;
      in >> f3;
      if (in.eof())
        return false;

      // strip out associated texture and normal vertices for f1
      idx = f1.find_first_of('/');
      if (idx != std::string::npos)
        f1.erase(idx, std::string::npos);

      // strip out associated texture and normal vertices for f2
      idx = f2.find_first_of('/');
      if (idx != std::string::npos)
        f2.erase(idx, std::string::npos);

      // strip out associated texture and normal vertices for f3
      idx = f3.find_first_of('/');
      if (idx != std::string::npos)
        f3.erase(idx, std::string::npos);

      // get the indices
      i1 = atoi(f1.c_str());
      i2 = atoi(f2.c_str());
      i3 = atoi(f3.c_str());

      // if the indices are negative, translate them to positive; otherwise, decrement them
      if (i1 < 0)
        i1 = vertices.size() + i1;
      else
        i1--;
  
      if (i2 < 0)
        i2 = vertices.size() + i2;
      else
        i2--;
    
      if (i3 < 0)
        i3 = vertices.size() + i3;
      else
        i3--;

      // create the index
      int* index = new int[3];
      index[0] = i1;
      index[1] = i2;
      index[2] = i3;

      // store it in the linked list, making new nodes as necessary
      faces.push_back(index);                        
      continue;
    }

    // read the rest of the line
    in.getline(buffer, BUF_SIZE);
  }

  // create the group and add the geometry to it
  osg::Geode* geode = new osg::Geode;
  osg::Group* group = new osg::Group;
  group->addChild(geode);
  osg::Geometry* geom = new osg::Geometry;
  geode->addDrawable(geom);

  // create the vertex array 
  osg::Vec3Array* varray = new osg::Vec3Array(vertices.size());
  for (unsigned i=0; i< vertices.size(); i++)
    (*varray)[i] = osg::Vec3((float) vertices[i][X], (float) vertices[i][Y], (float) vertices[i][Z]);
  geom->setVertexArray(varray);

  // create each face
  for (unsigned i=0; i< faces.size(); i++)
  {
    osg::DrawElementsUInt* face = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
    face->push_back(faces[i][A]);
    face->push_back(faces[i][B]);
    face->push_back(faces[i][C]);
    geom->addPrimitiveSet(face);
  }

  // create a new material with random color
  const float RED = (float) rand() / RAND_MAX;
  const float GREEN = (float) rand() / RAND_MAX;
  const float BLUE = (float) rand() / RAND_MAX;
  osg::Material* mat = new osg::Material;
  mat->setColorMode(osg::Material::DIFFUSE);
  mat->setDiffuse(osg::Material::FRONT, osg::Vec4(RED, GREEN, BLUE, 1.0f)); 
  group->getOrCreateStateSet()->setAttribute(mat);

  // free memory
  for (unsigned i=0; i< vertices.size(); i++)
    delete [] vertices[i];
  for (unsigned i=0; i< faces.size(); i++)
    delete [] faces[i];

  // close the file
  in.close();

  return group;
}

int main(int argc, char ** argv)
{
  if (argc < 2) 
  {
    std::cerr << "syntax: view <input1> [input2] [input3] ... [inputn]" << std::endl << std::endl;
    std::cerr << "Views Open Inventor (.iv) and Wavefront (.obj) files" << std::endl;
    std::cerr << "using Coin and SoQt" << std::endl;
    exit(1);
  }

  // setup the osg Viewer
  osgViewer::Viewer viewer;
  viewer.setCameraManipulator(new osgGA::TrackballManipulator());
  viewer.setThreadingModel(osgViewer::Viewer::SingleThreaded);
  viewer.addEventHandler(new osgGA::StateSetManipulator(viewer.getCamera()->getOrCreateStateSet()));
  viewer.addEventHandler(new osgViewer::WindowSizeHandler);
  viewer.addEventHandler(new osgViewer::StatsHandler);
  viewer.setLightingMode(osg::View::SKY_LIGHT);

  // setup the view window size
  viewer.setUpViewInWindow(0, 0, 640, 480);

  // setup the camera
  viewer.getCamera()->setViewMatrixAsLookAt(osg::Vec3d(0,0,10), osg::Vec3d(0,0,0), osg::Vec3d(0,1,0));

  // Open the argument file(s)..
  osg::Group* mainroot = new osg::Group;
  mainroot->ref();
  for (int i=1; i< argc; i++)
  {
    // do different things based on the file extension
    std::string fname(argv[i]);
    if (fname.rfind(".wrl") == fname.size() - 4 || fname.rfind(".iv") == fname.size() - 3)
    {
      // VRML 97 / 1.0 file
      osg::Node* node = osgDB::readNodeFile(fname);
      if (!node)
        std::cerr << "view: unable to open/read " << fname << std::endl;
      else
        mainroot->addChild(node);
    }
    else if (fname.rfind(".obj") == fname.size() - 4)
    {
      osg::Node* root = create_group_from_obj(argv[i]);
      if (!root) 
        std::cerr << "view: unable to open/read " << fname << std::endl;
      else
        mainroot->addChild(root);
    }
    else
    {
      std::cerr << "unknown filename type for " << fname << std::endl;
      return -1;
    }
  }

  // go
  viewer.setSceneData(mainroot);
  viewer.realize();
  while (true)
  {
    if (viewer.done())
      break;
    viewer.frame();
  }

  return 0;
}

