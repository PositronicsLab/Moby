// intersection.cpp : Defines the entry point for the console application.
#include <iostream>
#include <vector>
#include <cmath>
#include <Moby/Types.h>
#include <Moby/CompGeom.h>
#include <Moby/Polyhedron.h>
#include <Moby/XMLTree.h>
#include <Moby/GaussianMixture.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/OBB.h>
#ifdef USE_OSG
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Material>
#include <osg/LightModel>
#endif

using std::map;
using std::list;
using std::vector;
using std::pair;
using std::make_pair;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

/// Computes the OSG visualization
osg::Node* GaussianMixture::create_visualization()
{
  #ifdef USE_OSG
  const unsigned X = 0, Y = 1, Z = 2;

  // get the pose and compute transform from the global frame to it 
  shared_ptr<const Pose3d> P = get_pose();
  Transform3d T = Pose3d::calc_relative_pose(GLOBAL, P);

  // create necessary OSG elements for visualization
  osg::Group* group = new osg::Group;

  // loop over each Gaussian
  for (unsigned j=0; j< _gauss.size(); j++)
  {
    // create the necessary osg mechanicals
    osg::Group* subgroup = new osg::Group;
    osg::Geode* geode = new osg::Geode;
    osg::Geometry* geom = new osg::Geometry;
    geode->addDrawable(geom);
    subgroup->addChild(geode);
    group->addChild(subgroup); 

    // create a new material with random color
    const float RED = (float) rand() / RAND_MAX;
    const float GREEN = (float) rand() / RAND_MAX;
    const float BLUE = (float) rand() / RAND_MAX;

    // create a new material -- random color
    osg::Material* mat = new osg::Material;
    mat->setColorMode(osg::Material::DIFFUSE);
    mat->setDiffuse(osg::Material::FRONT, osg::Vec4(RED, GREEN, BLUE, 1.0f));
    subgroup->getOrCreateStateSet()->setAttribute(mat);

    // get the vertices for this Gaussian 
    const std::vector<Point3d>& verts = _vertices[j];

    // create the vertex array, backing current transform out
    osg::Vec3Array* varray = new osg::Vec3Array(verts.size());
    for (unsigned i=0; i< verts.size(); i++)
    {
      Point3d v = T.transform_point(verts[i]);
      (*varray)[i] = osg::Vec3((float) v[X], (float) v[Y], (float) v[Z]);
    }
    geom->setVertexArray(varray);

    // create a point cloud
    osg::DrawElementsUInt* cloud = new osg::DrawElementsUInt(osg::PrimitiveSet::POINTS, 0);
    for (unsigned i=0; i< verts.size(); i++)
      cloud->push_back(i);
    geom->addPrimitiveSet(cloud);
  }

  return group;
  #else
  return NULL;
  #endif 
}

/// Computes the distance and normal from a point on the CSG
double GaussianMixture::calc_dist_and_normal(const Point3d& p, Vector3d& normal) const
{
  // TODO: implement this
  assert(false);
  return 0.0;
}

/// Finds the signed distance between the sphere and another primitive
double GaussianMixture::calc_signed_dist(shared_ptr<const Primitive> primitive, shared_ptr<const Pose3d> pose_this, shared_ptr<const Pose3d> pose_primitive, Point3d& pthis, Point3d& pprimitive) const

{
  assert(false);
  return 0.0; 
}

/// Reads the Gaussian parameters
GaussianMixture::Gauss GaussianMixture::read_gauss_node(shared_ptr<const XMLTree> node)
{
  const unsigned X = 0, Y = 1;

  // setup default center and variance
  Vector2d x((double) 0.0, (double) 0.0), sigma((double) 1.0, (double) 1.0);

  // create the structure
  Gauss gauss;

  // read the height, if given 
  XMLAttrib* height = node->get_attrib("height");
  if (height)
    gauss.A = height->get_real_value();
  else
    gauss.A = (double) 1.0;

  // read the center, if given
  XMLAttrib* center = node->get_attrib("center");
  if (center)
    center->get_vector_value(x);
  gauss.x0 = x[X];
  gauss.y0 = x[Y];

  // read the variance, if given
  XMLAttrib* variance = node->get_attrib("variance");
  if (variance)
    variance->get_vector_value(sigma);
  gauss.sigma_x = sigma[X];
  gauss.sigma_y = sigma[Y];

  // read theta, if given
  XMLAttrib* theta = node->get_attrib("theta");
  if (theta)
    gauss.th = theta->get_real_value();
  else
    gauss.th = (double) 0.0;

  return gauss;
}

/// Gets the BVH root
BVPtr GaussianMixture::get_BVH_root(CollisionGeometryPtr geom)
{
  // verify that the geometry matches the stored geometry, if any
  if (_geom)
    assert(geom == _geom);
  else
  {
    // store the geometry
    _geom = geom;

    // build the bounding volume hierarchy
    construct_BVs(geom);
  }

  return _root;
}

/// Loads the primitive from an XML file
void GaussianMixture::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // load data from the Primitive
  Primitive::load_from_xml(node, id_map);

  // check that this node name is correct
  assert(strcasecmp(node->name.c_str(), "GaussianMixture") == 0);

  // read the Gaussian parameters
  vector<Gauss> gauss;
  list<shared_ptr<const XMLTree> > g_nodes = node->find_child_nodes("Gaussian");
  BOOST_FOREACH(shared_ptr<const XMLTree> g_node, g_nodes)
    gauss.push_back(read_gauss_node(g_node));

  // rebuild
  rebuild(gauss);
}

/// Saves the primitive to an XML description
void GaussianMixture::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save primitive data
  Primitive::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "GaussianMixture";

  // write all of the Gaussians
  for (unsigned i=0; i< _gauss.size(); i++)
  {
    // create a new subnode
    XMLTreePtr subnode(new XMLTree("Gaussian"));

    // write all parameters
    subnode->attribs.insert(XMLAttrib("height", _gauss[i].A));
    subnode->attribs.insert(XMLAttrib("center", Vector2d(_gauss[i].x0, _gauss[i].y0)));
    subnode->attribs.insert(XMLAttrib("variance", Vector2d(_gauss[i].sigma_x, _gauss[i].sigma_y)));
    subnode->attribs.insert(XMLAttrib("theta", _gauss[i].th));

    // add the child
    node->add_child(subnode);
  }
}

/// Returns the appropriate submesh
const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& GaussianMixture::get_sub_mesh(BVPtr bv)
{
  // gets the appropriate submesh
  OBBPtr obb = dynamic_pointer_cast<OBB>(bv);
  assert(obb);
  assert(_obbs.find(obb) != _obbs.end());
  unsigned index = _obbs.find(obb)->second;
  return _submesh[index];
}

/// Creates the triangle mesh
/**
 * This method functions by sampling over and then computing the convex hull
 * of each Gaussian. The convex hulls are then just overlaid.
 */
void GaussianMixture::create_mesh()
{
  unsigned j, k;
  double x, y;
  vector<Point3d> vertices;

  // number of samples per Gaussian dimension
  const unsigned NSAMPLES = 10;

  // setup empty initial mesh
  _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray);

  // resize submesh array
  _submesh.resize(_gauss.size());

  for (unsigned i=0; i< _gauss.size(); i++)
  {
    // clear the vector of vertices
    vertices.clear();

    // get x start and y start
    const double X_START = _gauss[i].x0 - _gauss[i].sigma_x*3.0;
    const double X_INC = _gauss[i].sigma_x*6.0/(NSAMPLES-1);
    const double Y_START = _gauss[i].y0 - _gauss[i].sigma_y*3.0;
    const double Y_INC = _gauss[i].sigma_y*6.0/(NSAMPLES-1);

    // iterate
    for (x = X_START, j = 0; j< NSAMPLES; j++, x += X_INC)
      for (y = Y_START, k = 0; k< NSAMPLES; k++, y += Y_INC)
        vertices.push_back(Point3d(x, y, gauss(_gauss[i], x, y)));

    // compute the 3D convex hull
    PolyhedronPtr poly = CompGeom::calc_convex_hull(vertices.begin(), vertices.end());

    // get the number of facets beforehand
    unsigned nbefore = _mesh->num_tris();

    // combine the meshes
    *_mesh = IndexedTriArray::merge(*_mesh, poly->get_mesh()); 

    // get the number of facets afterwards
    unsigned nafter = _mesh->num_tris();

    // setup the list of facets 
    list<unsigned> facets;
    for (unsigned j=nbefore; j< nafter; j++)
      facets.push_back(j);

    // setup the submesh array
    _submesh[i] = make_pair(_mesh, facets);
  }
}

/// Computes the Gaussian height at a point (NOTE: primitive transform is not accounted for!)
double GaussianMixture::gauss(const Gauss& g, double x, double y)
{
  const unsigned X = 0, Y = 1;

  // rotate x and y into the Gaussian space
  Matrix3d R = Matrix3d::rot_Z(g.th);
  Origin3d p = R.transpose_mult(Origin3d(x,y,(double) 0.0));
  x = p[X];
  y = p[Y];

  // precompute some things
  double k1 = (x - g.x0) * (x - g.x0);
  double k2 = (double) 2.0*g.sigma_x*g.sigma_x;
  double k3 = (y - g.y0) * (y - g.y0);
  double k4 = (double) 2.0*g.sigma_y*g.sigma_y;

  return g.A * std::exp(-(k1/k2 + k3/k4));
}

/// Computes the Gaussian gradient at a point on the Gaussian (NOTE: primitive transform is not accounted for!)
Vector3d GaussianMixture::grad(const Gauss& g, double x, double y)
{
  const unsigned X = 0, Y = 1;

  // rotate x and y into the Gaussian space
  Matrix3d R = Matrix3d::rot_Z(g.th);
  Origin3d p = R.transpose_mult(Origin3d(x,y,(double) 0.0));
  x = p[X];
  y = p[Y];

  // precompute some things
  double k1 = (x - g.x0) * (x - g.x0);
  double k2 = (double) 2.0*g.sigma_x*g.sigma_x;
  double k3 = (y - g.y0) * (y - g.y0);
  double k4 = (double) 2.0*g.sigma_y*g.sigma_y;

  // setup the gradient
  Vector3d v;
  v[0] = g.A * std::exp(-(k1/k2 + k3/k4)) * (double) -2.0 * (x - g.x0)/k2;
  v[1] = g.A * std::exp(-(k1/k2 + k3/k4)) * (double) -2.0 * (y - g.y0)/k4;

  // handle divide by zeros
  if (std::isnan(v[0]))
    v[0] = (double) 0.0;
  if (std::isnan(v[1]))
    v[1] = (double) 0.0;
  double nrm = v.norm();
  if (nrm < NEAR_ZERO)
    return Vector3d(0,0,1);
  else
    return v/nrm;
}

void GaussianMixture::get_vertices(vector<Point3d>& vertices)
{
  unsigned j, k;
  double x, y;

  // setup vertex mapping
  vertices.clear();

  // number of samples per Gaussian dimension
  const unsigned NSAMPLES = 100;

  for (unsigned i=0; i< _gauss.size(); i++)
  {
    // get x start and y start
    const double X_START = _gauss[i].x0 - _gauss[i].sigma_x*3.0;
    const double X_INC = _gauss[i].sigma_x*6.0/(NSAMPLES-1);
    const double Y_START = _gauss[i].y0 - _gauss[i].sigma_y*3.0;
    const double Y_INC = _gauss[i].sigma_y*6.0/(NSAMPLES-1);

    // iterate
    for (x = X_START, j = 0; j< NSAMPLES; j++, x += X_INC)
      for (y = Y_START, k = 0; k< NSAMPLES; k++, y += Y_INC)
      {
        // get the height for this Gaussian
        double hi = gauss(_gauss[i], x, y);

        // if the height is zero, we can use it directly
        if (hi > NEAR_ZERO)
        {
          // otherwise, verify that point is not inside any (other) Gaussians
          double maxh = (double) 0.0;
          for (unsigned r=0; r< _gauss.size() && maxh < NEAR_ZERO; r++)
          {
            if (r == i)
              continue;

            // get the height
            double h = gauss(_gauss[r], x, y);
            if (h > maxh)
              maxh = h;
           }

           // see whether maxh is sufficiently large
           if (maxh > hi)
             continue;
         }

        // add the point
        vertices.push_back(Point3d(x, y, hi, get_pose()));
      }
  }
}

/*
void GaussianMixture::get_vertices(BVPtr bv, vector<const Point3d*>& vertices)
{
  typedef map<OBBPtr, unsigned>::const_iterator OBBMapIter;

  // determine which Gaussian the bounding volume corresponds to   
  OBBPtr obb = dynamic_pointer_cast<OBB>(bv);
  assert(obb);

  // if the OBB has no mapping, we need *all* vertices
  OBBMapIter iter = _obbs.find(obb);
  if (iter == _obbs.end())
  {
    vertices.clear();
    for (unsigned j=0; j< _vertices.size(); j++)
    {
      const vector<Point3d>& verts = _vertices[j]; 
      for (unsigned i=0; i< verts.size(); i++)
        vertices.push_back(&verts[i]);
    }

    return;
  }

  // otherwise, get the vector of vertices
  const vector<Point3d>& verts = _vertices[iter->second]; 
  vertices.resize(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    vertices[i] = &verts[i];
}
*/

/// Rebuilds everything from scratch
void GaussianMixture::rebuild(const vector<Gauss>& gauss)
{
  // copy the Gaussians
  _gauss = gauss;

  // construct sets of vertices
  construct_vertices();

  // recreate the mesh
  create_mesh();
}

/// Sets the transform for the primitive
void GaussianMixture::set_pose(const Pose3d& T)
{
  // call parent method
  Primitive::set_pose(T);

  // reconstruct bounding volumes and vertices
  rebuild(_gauss);

  // re-construct the BVs
  if (_geom)
    construct_BVs(_geom);
}

/// Constructs the bounding volumes used for faster intersection testing
void GaussianMixture::construct_BVs(CollisionGeometryPtr geom)
{
  const unsigned X = 0, Y = 1, Z = 2;
  list<BVPtr> children;

  // get the collision geometry pose
  shared_ptr<const Pose3d> gpose = geom->get_pose();

  // get the pose
  shared_ptr<const Pose3d> P = get_pose();
  assert(!P->rpose);

  // setup a transform
  Transform3d T;
  T.source = gpose;
  T.target = gpose;
  T.q = P->q;
  T.x = P->x;

  // iterate over all Gaussians
  for (unsigned i=0; i< _gauss.size(); i++)
  {
    // create an OBB for the Gaussian
    shared_ptr<OBB> obb(new OBB);
    obb->geom = geom;

    // setup the OBB center
    const double HEIGHT = gauss(_gauss[i], _gauss[i].x0, _gauss[i].y0);
    obb->center = T.transform_point(Point3d(_gauss[i].x0, _gauss[i].y0, HEIGHT*0.5, gpose));

    // setup the OBB half-lengths
    obb->l[X] = _gauss[i].sigma_x*3.0;
    obb->l[Y] = _gauss[i].sigma_y*3.0;
    obb->l[Z] = HEIGHT*0.5;

    // setup the R matrix
    obb->R = T.q * Matrix3d::rot_Z(_gauss[i].th);

    // add the OBB as a child of the root
    children.push_back(obb);

    // setup mapping from BV to Gaussian
    _obbs[obb] = i;
  }

  // get the vertices from all children
  vector<Point3d> verts;
  BOOST_FOREACH(BVPtr bv, children)
  {
    shared_ptr<OBB> obb = dynamic_pointer_cast<OBB>(bv);
    obb->get_vertices(std::back_inserter(verts));
  }

  // setup an OBB over *all* Gaussians
  _root = shared_ptr<OBB>(new OBB(verts.begin(), verts.end()));
  _root->geom = geom;
  _root->children = children;
}

/// Constructs the set of vertices
void GaussianMixture::construct_vertices()
{
  unsigned j, k;
  double x, y;

  // clear vector of vertices
  _vertices.resize(_gauss.size());

  // get the transform from the current pose to the global frame
  shared_ptr<const Pose3d> P = get_pose();
  Transform3d T = Pose3d::calc_relative_pose(P, GLOBAL);

  // number of samples per Gaussian dimension
  const unsigned NSAMPLES = 100;

  for (unsigned i=0; i< _gauss.size(); i++)
  {
    // setup vertex mapping
    vector<Point3d>& v = _vertices[i];
    v.clear();

    // get x start and y start
    const double X_START = _gauss[i].x0 - _gauss[i].sigma_x*3.0;
    const double X_INC = _gauss[i].sigma_x*6.0/(NSAMPLES-1);
    const double Y_START = _gauss[i].y0 - _gauss[i].sigma_y*3.0;
    const double Y_INC = _gauss[i].sigma_y*6.0/(NSAMPLES-1);

    // iterate
    for (x = X_START, j = 0; j< NSAMPLES; j++, x += X_INC)
      for (y = Y_START, k = 0; k< NSAMPLES; k++, y += Y_INC)
      {
        // get the height for this Gaussian
        double hi = gauss(_gauss[i], x, y);

        // if the height is zero, we can use it directly
        if (hi > NEAR_ZERO)
        {
          // otherwise, verify that point is not inside any (other) Gaussians
          double maxh = (double) 0.0;
          for (unsigned r=0; r< _gauss.size() && maxh < NEAR_ZERO; r++)
          {
            if (r == i)
              continue;

            // get the height
            double h = gauss(_gauss[r], x, y);
            if (h > maxh)
              maxh = h;
           }

           // see whether maxh is sufficiently large
           if (maxh > hi)
             continue;
         }

        // add the point
        v.push_back(Point3d(x, y, hi, P));
      }

    // transform all points
    for (unsigned i=0; i< v.size(); i++)
      v[i] = T.transform_point(v[i]);
  }
}

/// Determines whether a point is inside one of the Gaussians
bool GaussianMixture::point_inside(BVPtr bv, const Point3d& point, Vector3d& normal) const
{
  double a,b,c,X,Y,Z,X_next,eps;
  static shared_ptr<Pose3d> P;
  
  const double PARAM_BOUND = -5.0;
  int n=0;
  const int NMAX = _gauss.size(); 

  // get the pose for the collision geometry
  shared_ptr<const Pose3d> gpose = bv->geom->get_pose(); 

  // get the pose for this geometry and BV
  shared_ptr<const Pose3d> bpose = get_pose(); 
  assert(!bpose->rpose);

  // setup a new pose
  if (!P)
    P = shared_ptr<Pose3d>(new Pose3d);
  *P = *bpose;
  P->rpose = gpose;

  // get the transform from the point pose to the Gaussian pose
  Transform3d T = Pose3d::calc_relative_pose(point.pose, P); 

  // convert the point to primitive space
  Point3d p = T.transform_point(point);

  // find max(z) of gaussians
  double tempX,tempY,tempZ,tempMax;
  tempX=p[0];
  tempY=p[1];
  tempMax = -std::numeric_limits<double>::max();
  int num = 0;
  for(int i=0;i<_gauss.size();i++){

    // compute some constants
    const double CTH = std::cos(_gauss[i].th);
    const double STH = std::sin(_gauss[i].th); 
    const double S2TH = std::sin((double) 2.0 * _gauss[i].th); 

    a = CTH*CTH/2/_gauss[i].sigma_x*_gauss[i].sigma_x + STH*STH/2/_gauss[i].sigma_y*_gauss[i].sigma_y;
    b = -S2TH/4/_gauss[i].sigma_x*_gauss[i].sigma_x + S2TH/4/_gauss[i].sigma_y*_gauss[i].sigma_y;
    c = STH*STH/2/_gauss[i].sigma_x*_gauss[i].sigma_x + CTH*CTH/2/_gauss[i].sigma_y*_gauss[i].sigma_y;
    tempZ = _gauss[i].A*std::exp( - (a*(tempX-_gauss[i].x0)*a*(tempX-_gauss[i].x0) + 2*b*(tempX-_gauss[i].x0)*(tempY-_gauss[i].y0) + c*(tempY-_gauss[i].y0)*(tempY-_gauss[i].y0)));
    //tempZ=_gauss[i].A/std::exp((tempX - _gauss[i].x0)*(tempX - _gauss[i].x0)/(2*_gauss[i].sigma_x) + ((tempY - _gauss[i].y0)*(tempY - _gauss[i].y0))/(2*_gauss[i].sigma_y)); 
    if(tempZ>tempMax){
      tempMax=tempZ;
      num = i;
    }
  }

  //*************************************************************
  // check p is already inside
  if(p[2] <= tempMax)
  {
    // compute some constants
    const double CTH = std::cos(_gauss[num].th);
    const double STH = std::sin(_gauss[num].th); 
    const double S2TH = std::sin((double) 2.0 * _gauss[num].th); 

    // normal -------------------------------------------------------------
    a = CTH/2/_gauss[num].sigma_x*_gauss[num].sigma_x + STH/2/_gauss[num].sigma_y*_gauss[num].sigma_y;
    b = -S2TH/4/_gauss[num].sigma_x*_gauss[num].sigma_x + S2TH/4/_gauss[num].sigma_y*_gauss[num].sigma_y;
    c = STH*STH/2/_gauss[num].sigma_x*_gauss[num].sigma_x + CTH*CTH/2/_gauss[num].sigma_y*_gauss[num].sigma_y;
      
    X = tempX;
    Y = tempY;
    Z = tempMax;
    // w.r.t F(x,y)
    normal[0] = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(Y - _gauss[num].y0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0))*((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(2*X - 2*_gauss[num].x0));
    normal[1] = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(Y - _gauss[num].y0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0))*((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(2*Y - 2*_gauss[num].y0));
    normal[2] = +1;

    // normalize the normal
    normal.normalize();
    normal = T.inverse_transform_vector(normal);

    if(tempMax-p[2] <= std::exp(PARAM_BOUND*2)) //error bound
    {      
      // cout<<"inside errorbound= "<<endl;
      // cout<<"X= "<<isect[0] <<" Y= "<<isect[1] <<" Z= "<<isect[2] <<endl;
      return true;
    }
    if(a==0 && b==0) // no tangent
    {
      // cout<<"X= "<<isect[0] <<" Y= "<<isect[1] <<" Z= "<<isect[2] <<endl;
      return true;
    }
  }

  // if we're here, point is outside
  return false;
}

/// Evaluates the intersection function (for Newton-Raphson)
double GaussianMixture::f(const Gauss& g, const Point3d& p, const Point3d& q, double t)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // evaluate
  double x = p[X] + (q[X] - p[X])*t;
  double y = p[Y] + (q[Y] - p[Y])*t;
  return p[Z] + (q[Z] - p[Z])*t - gauss(g, x, y);
}

/// Evaluates the intersection function (for Newton-Raphson)
double GaussianMixture::df(const Gauss& g, const Point3d& p, const Point3d& q, double t)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // evaluate
  double dx = (q[X] - p[X]);
  double dy = (q[Y] - p[Y]);
  double x = p[X] + (q[X] - p[X])*t;
  double y = p[Y] + (q[Y] - p[Y])*t;

  // rotate x and y into the Gaussian space
  Matrix3d R = Matrix3d::rot_Z(g.th);
  Origin3d r = R.transpose_mult(Origin3d(x,y,(double) 0.0));
  Origin3d dr = R.transpose_mult(Origin3d(dx, dy, (double) 0.0));
  double xx = r[X];
  double yy = r[Y];
  double dxx = dr[X];
  double dyy = dr[Y];

  // precompute some things
  double k1 = (xx - g.x0) * (xx - g.x0);
  double k2 = (double) 2.0*g.sigma_x*g.sigma_x;
  double k3 = (yy - g.y0) * (yy - g.y0);
  double k4 = (double) 2.0*g.sigma_y*g.sigma_y;

  // compute derivatives of these
  double dk1 = (double) 2.0 * (xx - g.x0) * dxx;
  double dk3 = (double) 2.0 * (yy - g.y0) * dyy;

  return q[Z] - p[Z] - g.A * std::exp(-(k1/k2 + k3/k4)) * -(dk1/k2 + dk3/k4);
}

/// Performs Newton-Raphson to find a point of intersection with a line segment and a Gaussian
double GaussianMixture::newton_raphson(const Gauss& g, const Point3d& p, const Point3d& q)
{
  const unsigned ITER_MAX = 20;
  const double INF = std::numeric_limits<double>::max();
  const double TTOL = NEAR_ZERO;

  double t = (double) 0.5;
  for (unsigned j=0; j< ITER_MAX; j++)
  {
    // do update
    double dt = -f(g, p, q, t)/df(g, p, q, t);
    t += dt;

    // check for convergence
    if (std::fabs(dt) < TTOL)
    {
      // look for bad convergence
      if (t < -NEAR_ZERO)
        return INF;
      else if (t < (double) 0.0)
        return (double) 0.0; 
      else
        return t;
    }
  }
  
  // no convergence? return INF
  return INF;
}

/// Computes line segment intersection (if any)
bool GaussianMixture::intersect_seg(BVPtr bv, const LineSeg3& seg,double& tisect,Point3d& isect, Vector3d& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const double INF = std::numeric_limits<double>::max();
  static shared_ptr<Pose3d> P;

  // setup intersection vectors
  SAFESTATIC vector<double> t, depth;

  // get the pose for the collision geometry
  shared_ptr<const Pose3d> gpose = bv->geom->get_pose(); 

  // get the pose for this geometry and BV
  shared_ptr<const Pose3d> bpose = get_pose(); 
  assert(!bpose->rpose);

  // setup a new pose
  if (!P)
    P = shared_ptr<Pose3d>(new Pose3d);
  *P = *bpose;
  P->rpose = gpose;

  // get the transform from the line segment pose to the primitive pose
  Transform3d T = Pose3d::calc_relative_pose(seg.first.pose, P); 

  // convert the line segment to primitive space
  Point3d p = T.transform_point(seg.first);
  Point3d q = T.transform_point(seg.second);

  // determine whether any starting points are inside the Gaussians
  depth.resize(_gauss.size());
  for (unsigned i=0; i< _gauss.size(); i++)
    depth[i] = p[Z] - gauss(_gauss[i], p[X], p[Y]);

  // see whether any points are inside the Gaussians
  unsigned mini = std::min_element(depth.begin(), depth.end()) - depth.begin();
  if (depth[mini] < (double) 0.0)
  {
    // point is inside, compute and transform the normal
    tisect = (double) 0.0;
    isect = seg.first; 
    
    // compute the transformed normal
    normal = T.inverse_transform_vector(grad(_gauss[mini], p[X], p[Y]));

    return true;
  }

  // no points are inside; make sure that a point on the line segment is
  // inside
  for (unsigned i=0; i< _gauss.size(); i++)
    depth[i] = q[Z] - gauss(_gauss[i], q[X], q[Y]);

  // see whether all points are outside
  mini = std::min_element(depth.begin(), depth.end()) - depth.begin();
  if (depth[mini] > (double) 0.0)
    return false;

  // still here? use Newton-Raphson 
  t.resize(_gauss.size());
  for (unsigned i=0; i< _gauss.size(); i++)
    t[i] = INF;

  // compute t
  for (unsigned i=0; i< _gauss.size(); i++)
  {
    // only apply to appropriate points 
    if (depth[i] > (double) 0.0)
      continue;

    // apply Newton-Raphson
    t[i] = newton_raphson(_gauss[i], p, q);
  }

  // get the first intersection
  mini = std::min_element(t.begin(), t.end()) - t.begin();

  // compute the point of intersection and normal
  tisect = t[mini];
  isect = p + (q-p)*tisect;

  // compute transformed normal
  double x = p[X] + (q[X] - p[X])*tisect;
  double y = p[Y] + (q[Y] - p[Y])*tisect;
  normal = T.inverse_transform_vector(grad(_gauss[mini], x, y));

  // compute and transform intersection point
  isect = seg.first + (seg.second-seg.first)*tisect;

  return true;

/*
  // starting point of line -> q vector
  // Find x1,y1,z1 ->q vector
  x1=q[0];
  y1=q[1];
  z1=q[2];
 
   // seg.second is p vector
  double tempX,tempY,tempZ,tempMax;
  tempX=p[0];
  tempY=p[1];
  tempMax = -std::numeric_limits<double>::max();
  int num = 0;

  // find max(z) of gaussians (we need to find out which Gaussian is tallest
  // at the first point in the line segment)
  for(unsigned i=0; i<_gauss.size();i++){

    // compute some constants
    const double CTH = std::cos(_gauss[i].th);
    const double STH = std::sin(_gauss[i].th); 
    const double S2TH = std::sin((double) 2.0 * _gauss[i].th); 

    a = CTH*CTH/2/_gauss[i].sigma_x*_gauss[i].sigma_x + STH*STH/2/_gauss[i].sigma_y*_gauss[i].sigma_y;
    b = -S2TH/4/_gauss[i].sigma_x*_gauss[i].sigma_x + S2TH/4/_gauss[i].sigma_y*_gauss[i].sigma_y;
    c = STH*STH/2/_gauss[i].sigma_x*_gauss[i].sigma_x + CTH*CTH/2/_gauss[i].sigma_y*_gauss[i].sigma_y;
    tempZ = _gauss[i].A*std::exp( - (a*(tempX-_gauss[i].x0)*a*(tempX-_gauss[i].x0) + 2*b*(tempX-_gauss[i].x0)*(tempY-_gauss[i].y0) + c*(tempY-_gauss[i].y0)*(tempY-_gauss[i].y0)));
    //tempZ=_gauss[i].A/std::exp((tempX - _gauss[i].x0)*(tempX - _gauss[i].x0)/(2*_gauss[i].sigma_x) + ((tempY - _gauss[i].y0)*(tempY - _gauss[i].y0))/(2*_gauss[i].sigma_y)); 
    if(tempZ>tempMax){
      tempMax=tempZ;
      num = i;
    }
  }

  // *************************************************************
  // check whether the first point in the line segment is already inside
  if(p[2] <= tempMax)
  {
    t = (double) 0.0;
    isect = p;

    // compute some constants
    const double CTH = std::cos(_gauss[num].th);
    const double STH = std::sin(_gauss[num].th); 
    const double S2TH = std::sin((double) 2.0 * _gauss[num].th); 

    // normal -------------------------------------------------------------
    a = CTH/2/_gauss[num].sigma_x*_gauss[num].sigma_x + STH/2/_gauss[num].sigma_y*_gauss[num].sigma_y;
    b = -S2TH/4/_gauss[num].sigma_x*_gauss[num].sigma_x + S2TH/4/_gauss[num].sigma_y*_gauss[num].sigma_y;
    c = STH*STH/2/_gauss[num].sigma_x*_gauss[num].sigma_x + CTH*CTH/2/_gauss[num].sigma_y*_gauss[num].sigma_y;
      
    X = tempX;
    Y = tempY;
    Z = tempMax;
    // w.r.t F(x,y)
    normal[0] = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(Y - _gauss[num].y0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0))*((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(2*X - 2*_gauss[num].x0));
    normal[1] = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(Y - _gauss[num].y0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0))*((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(2*Y - 2*_gauss[num].y0));
    normal[2] = +1;

    // normalize the normal
    normal.normalize();

    if(tempMax-p[2] <= std::exp(PARAM_BOUND*2)) //error bound
    {      
      // cout<<"inside errorbound= "<<endl;
      // cout<<"X= "<<isect[0] <<" Y= "<<isect[1] <<" Z= "<<isect[2] <<endl;

      // transform intersection point and normal
      isect = T.transform_point(isect);
      normal = T.transform_vector(normal);

      return true;
    }
    if(a==0 && b==0) // no tangent
    {
      isect[2]=tempMax;      
      // cout<<"X= "<<isect[0] <<" Y= "<<isect[1] <<" Z= "<<isect[2] <<endl;

      // transform intersection point and normal
      isect = T.transform_point(isect);
      normal = T.transform_vector(normal);

      return true;
    }
  }
  else
  {
    // cout<<"no intersection"<<endl;
    return false;
  }

  // otherwise put p is initial value of X and Y.
  X=p[0];
  Y=p[1];
  // Numerical anaysis
  eps=1;num=0;
  while(eps>= std::exp (PARAM_BOUND) && n<=NMAX)
  {
    double tempZ,tempdiffZ;
    double tempMax = -std::numeric_limits<double>::max();
    num = 0;
    // find Max F(x,y)
    for(int i=0;i<_gauss.size();i++){

      // compute some constants
      const double CTH = std::cos(_gauss[i].th);
      const double STH = std::sin(_gauss[i].th); 
      const double S2TH = std::sin((double) 2.0 * _gauss[i].th); 

      a = CTH*CTH/2/_gauss[i].sigma_x*_gauss[i].sigma_x + STH*STH/2/_gauss[i].sigma_y*_gauss[i].sigma_y;
      b = -S2TH/4/_gauss[i].sigma_x*_gauss[i].sigma_x + S2TH/4/_gauss[i].sigma_y*_gauss[i].sigma_y;
      c = STH*STH/2/_gauss[i].sigma_x*_gauss[i].sigma_x + CTH*CTH/2/_gauss[i].sigma_y*_gauss[i].sigma_y;
         
      tempZ = _gauss[i].A/std::exp((a*(X-_gauss[i].x0)*(X-_gauss[i].x0) + 2*b*(X-_gauss[i].x0)*(Y-_gauss[i].y0) + c*(Y-_gauss[i].y0)*(Y-_gauss[i].y0)));      
      if(tempZ>tempMax){
        tempMax=tempZ;
        num = i;
      }
    }

    // compute some constants
    const double CTH = std::cos(_gauss[num].th);
    const double STH = std::sin(_gauss[num].th); 
    const double S2TH = std::sin((double) 2.0 * _gauss[num].th); 

    a = CTH*CTH/2/_gauss[num].sigma_x*_gauss[num].sigma_x + STH*STH/2/_gauss[num].sigma_y*_gauss[num].sigma_y;
    b = -S2TH/4/_gauss[num].sigma_x*_gauss[num].sigma_x + S2TH/4/_gauss[num].sigma_y*_gauss[num].sigma_y;
    c = STH*STH/2/_gauss[num].sigma_x*_gauss[num].sigma_x + CTH*CTH/2/_gauss[num].sigma_y*_gauss[num].sigma_y;

    tempZ = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(y1 - _gauss[num].y0 + (d[1]*(X - x1))/d[0]) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(y1 - _gauss[num].y0 + (d[1]*(X - x1))/d[0])*(y1 - _gauss[num].y0 + (d[1]*(X - x1))/d[0]) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0)) - z1 - (d[2]*(X - x1))/d[0];
    tempdiffZ = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(y1 - _gauss[num].y0 + (d[1]*(X - x1))/d[0])    - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(y1 - _gauss[num].y0 + (d[1]*(X - x1))/d[0])*(y1 - _gauss[num].y0 + (d[1]*(X - x1))/d[0]) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0))*((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(y1 - _gauss[num].y0 + (d[1]*(X - x1))/d[0]) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(2*X - 2*_gauss[num].x0) + (d[1]*(S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0))/d[0] - (2*d[1]*(CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(y1 - _gauss[num].y0 + (d[1]*(X - x1))/d[0]))/d[0]) - d[2]/d[0];

    //cout<<"X= "<<X<<endl;
    X_next = X - tempZ/tempdiffZ;
    eps = std::fabs(X_next-X);
    X= X_next;
    assert(!std::isnan(X));
    Y=d[1]/d[0]*(X-x1)+y1;
    n++;
  }
  if(n>NMAX)
    return false;
  //
  t=(X-p[0])/d[0];

  Y= d[1]/d[0]*(X-x1) +y1;
  Z= d[2]/d[0]*(X-x1) +z1;
  isect[0] = X;
  isect[1] = Y;
  isect[2] = Z;
  // normal -------------------------------------------------------------
  // w.r.t F(x,y)

  // compute some constants
  const double CTH = std::cos(_gauss[num].th);
  const double STH = std::sin(_gauss[num].th); 
  const double S2TH = std::sin((double) 2.0 * _gauss[num].th); 

  normal[0] = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(Y - _gauss[num].y0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0))*((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(2*X - 2*_gauss[num].x0));
  normal[1] = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(Y - _gauss[num].y0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0))*((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(2*Y - 2*_gauss[num].y0));
  normal[2] = +1;

  // normalize the normal
  normal.normalize();

  // transform intersection points and normals
  isect = T.transform_point(isect);
  normal = T.transform_vector(normal);

  cout <<"X= "<<isect[0]<<" Y= "<<isect[1]<<" Z= "<<isect[2]<<endl;
  cout <<"normal[0]= "<<normal[0]<<" normal[1]= "<<normal[1]<<" normal[z]= "<<normal[z]<<endl;

  return true;
  */
}

