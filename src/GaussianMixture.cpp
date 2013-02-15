// intersection.cpp : Defines the entry point for the console application.
#include <iostream>
#include <vector>
#include <cmath>
#include <Moby/Types.h>
#include <Moby/CompGeom.h>
#include <Moby/Polyhedron.h>
#include <Moby/XMLTree.h>
#include <Moby/GaussianMixture.h>
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
using namespace Moby;

/// Computes the OSG visualization
#ifdef USE_OSG
osg::Node* GaussianMixture::create_visualization()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the transform
  const Matrix4& T = get_transform();

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
    const std::vector<Vector3>& verts = _vertices[j];

    // create the vertex array, backing current transform out
    osg::Vec3Array* varray = new osg::Vec3Array(verts.size());
    for (unsigned i=0; i< verts.size(); i++)
    {
      Vector3 v = T.inverse_mult_point(verts[i]);
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
}
#endif 

/// Reads the Gaussian parameters
GaussianMixture::Gauss GaussianMixture::read_gauss_node(XMLTreeConstPtr node)
{
  const unsigned X = 0, Y = 1;

  // setup default center and variance
  Vector2 x((Real) 0.0, (Real) 0.0), sigma((Real) 1.0, (Real) 1.0);

  // create the structure
  Gauss gauss;

  // read the height, if given 
  const XMLAttrib* height = node->get_attrib("height");
  if (height)
    gauss.A = height->get_real_value();
  else
    gauss.A = (Real) 1.0;

  // read the center, if given
  const XMLAttrib* center = node->get_attrib("center");
  if (center)
    center->get_vector_value(x);
  gauss.x0 = x[X];
  gauss.y0 = x[Y];

  // read the variance, if given
  const XMLAttrib* variance = node->get_attrib("variance");
  if (variance)
    variance->get_vector_value(sigma);
  gauss.sigma_x = sigma[X];
  gauss.sigma_y = sigma[Y];

  // read theta, if given
  const XMLAttrib* theta = node->get_attrib("theta");
  if (theta)
    gauss.th = theta->get_real_value();
  else
    gauss.th = (Real) 0.0;

  return gauss;
}

/// Loads the primitive from an XML file
void GaussianMixture::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  // load data from the Primitive
  Primitive::load_from_xml(node, id_map);

  // check that this node name is correct
  assert(strcasecmp(node->name.c_str(), "GaussianMixture") == 0);

  // read the Gaussian parameters
  vector<Gauss> gauss;
  list<XMLTreeConstPtr> g_nodes = node->find_child_nodes("Gaussian");
  BOOST_FOREACH(XMLTreeConstPtr g_node, g_nodes)
    gauss.push_back(read_gauss_node(g_node));

  // rebuild
  rebuild(gauss);
}

/// Saves the primitive to an XML description
void GaussianMixture::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
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
    subnode->attribs.insert(XMLAttrib("center", Vector2(_gauss[i].x0, _gauss[i].y0)));
    subnode->attribs.insert(XMLAttrib("variance", Vector2(_gauss[i].sigma_x, _gauss[i].sigma_y)));
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
  Real x, y;
  vector<Vector3> vertices;

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
    const Real X_START = _gauss[i].x0 - _gauss[i].sigma_x*3.0;
    const Real X_INC = _gauss[i].sigma_x*6.0/(NSAMPLES-1);
    const Real Y_START = _gauss[i].y0 - _gauss[i].sigma_y*3.0;
    const Real Y_INC = _gauss[i].sigma_y*6.0/(NSAMPLES-1);

    // iterate
    for (x = X_START, j = 0; j< NSAMPLES; j++, x += X_INC)
      for (y = Y_START, k = 0; k< NSAMPLES; k++, y += Y_INC)
        vertices.push_back(Vector3(x, y, gauss(_gauss[i], x, y)));

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
Real GaussianMixture::gauss(const Gauss& g, Real x, Real y)
{
  const unsigned X = 0, Y = 1;

  // rotate x and y into the Gaussian space
  Matrix3 R = Matrix3::rot_Z(g.th);
  Vector3 p = R.transpose_mult(Vector3(x,y,(Real) 0.0));
  x = p[X];
  y = p[Y];

  // precompute some things
  Real k1 = (x - g.x0) * (x - g.x0);
  Real k2 = (Real) 2.0*g.sigma_x*g.sigma_x;
  Real k3 = (y - g.y0) * (y - g.y0);
  Real k4 = (Real) 2.0*g.sigma_y*g.sigma_y;

  return g.A * std::exp(-(k1/k2 + k3/k4));
}

void GaussianMixture::get_vertices(BVPtr bv, vector<const Vector3*>& vertices)
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
      const vector<Vector3>& verts = _vertices[j]; 
      for (unsigned i=0; i< verts.size(); i++)
        vertices.push_back(&verts[i]);
    }

    return;
  }

  // otherwise, get the vector of vertices
  const vector<Vector3>& verts = _vertices[iter->second]; 
  vertices.resize(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    vertices[i] = &verts[i];
}

/// Rebuilds everything from scratch
void GaussianMixture::rebuild(const vector<Gauss>& gauss)
{
  // copy the Gaussians
  _gauss = gauss;

  // construct bounding volumes
  construct_BVs();

  // construct sets of vertices
  construct_vertices();

  // recreate the mesh
  create_mesh();
}

/// Sets the transform for the primitive
void GaussianMixture::set_transform(const Matrix4& T)
{
  // call parent method
  Primitive::set_transform(T);

  // reconstruct bounding volumes and vertices
  rebuild(_gauss);
}

/// Constructs the bounding volumes used for faster intersection testing
void GaussianMixture::construct_BVs()
{
  const unsigned X = 0, Y = 1, Z = 2;
  list<BVPtr> children;

  // get the transform
  const Matrix4& T = get_transform();
  
  // iterate over all Gaussians
  for (unsigned i=0; i< _gauss.size(); i++)
  {
    // create an OBB for the Gaussian
    shared_ptr<OBB> obb(new OBB);

    // setup the OBB center
    const Real HEIGHT = gauss(_gauss[i], _gauss[i].x0, _gauss[i].y0);
    obb->center[X] = _gauss[i].x0;
    obb->center[Y] = _gauss[i].y0;
    obb->center[Z] = HEIGHT*0.5;

    // setup the OBB half-lengths
    obb->l[X] = _gauss[i].sigma_x*3.0;
    obb->l[Y] = _gauss[i].sigma_y*3.0;
    obb->l[Z] = HEIGHT*0.5;

    // setup the R matrix
    obb->R = Matrix3::rot_Z(_gauss[i].th);

    // transform the obb
    obb->center = T.mult_point(obb->center);
    obb->R = T.get_rotation() * obb->R;

    // add the OBB as a child of the root
    children.push_back(obb);

    // setup mapping from BV to Gaussian
    _obbs[obb] = i;
  }

  // get the vertices from all children
  vector<Vector3> verts;
  BOOST_FOREACH(BVPtr bv, children)
  {
    shared_ptr<OBB> obb = dynamic_pointer_cast<OBB>(bv);
    obb->get_vertices(std::back_inserter(verts));
  }

  // setup an OBB over *all* Gaussians
  _root = shared_ptr<OBB>(new OBB(verts.begin(), verts.end()));
  _root->children = children;
}

/// Constructs the set of vertices
void GaussianMixture::construct_vertices()
{
  unsigned j, k;
  Real x, y;

  // clear vector of vertices
  _vertices.resize(_gauss.size());

  // get the transform
  const Matrix4& T = get_transform();

  // number of samples per Gaussian dimension
  const unsigned NSAMPLES = 100;

  for (unsigned i=0; i< _gauss.size(); i++)
  {
    // setup vertex mapping
    vector<Vector3>& v = _vertices[i];
    v.clear();

    // get x start and y start
    const Real X_START = _gauss[i].x0 - _gauss[i].sigma_x*3.0;
    const Real X_INC = _gauss[i].sigma_x*6.0/(NSAMPLES-1);
    const Real Y_START = _gauss[i].y0 - _gauss[i].sigma_y*3.0;
    const Real Y_INC = _gauss[i].sigma_y*6.0/(NSAMPLES-1);

    // iterate
    for (x = X_START, j = 0; j< NSAMPLES; j++, x += X_INC)
      for (y = Y_START, k = 0; k< NSAMPLES; k++, y += Y_INC)
      {
        // get the height for this Gaussian
        Real hi = gauss(_gauss[i], x, y);

        // if the height is zero, we can use it directly
        if (hi > NEAR_ZERO)
        {
          // otherwise, verify that point is not inside any (other) Gaussians
          Real maxh = (Real) 0.0;
          for (unsigned r=0; r< _gauss.size() && maxh < NEAR_ZERO; r++)
          {
            if (r == i)
              continue;

            // get the height
            Real h = gauss(_gauss[r], x, y);
            if (h > maxh)
              maxh = h;
           }

           // see whether maxh is sufficiently large
           if (maxh > hi)
             continue;
         }

        // add the point
        v.push_back(Vector3(x, y, hi));
      }

    // transform all points
    for (unsigned i=0; i< v.size(); i++)
      v[i] = T.mult_point(v[i]);
  }
}

/// Determines whether a point is inside one of the Gaussians
bool GaussianMixture::point_inside(BVPtr bv, const Vector3& point, Vector3& normal) const
{
  Real a,b,c,X,Y,Z,X_next,eps;
  
  const Real PARAM_BOUND = -5.0;
  int n=0;
  const int NMAX = _gauss.size(); 

  // get the current transform
  const Matrix4& T = get_transform();
 
  // convert the point to primitive space
  Vector3 p = T.inverse_mult_point(point);

  // find max(z) of gaussians
  Real tempX,tempY,tempZ,tempMax;
  tempX=p[0];
  tempY=p[1];
  tempMax = -std::numeric_limits<Real>::max();
  int num = 0;
  for(int i=0;i<_gauss.size();i++){

    // compute some constants
    const Real CTH = std::cos(_gauss[i].th);
    const Real STH = std::sin(_gauss[i].th); 
    const Real S2TH = std::sin((Real) 2.0 * _gauss[i].th); 

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
    const Real CTH = std::cos(_gauss[num].th);
    const Real STH = std::sin(_gauss[num].th); 
    const Real S2TH = std::sin((Real) 2.0 * _gauss[num].th); 

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

/// Computes line segment intersection (if any)
bool GaussianMixture::intersect_seg(BVPtr bv, const LineSeg3& seg,Real& t,Vector3& isect, Vector3& normal) const
{
  Real a,b,c,X,Y,Z,x1,y1,z1,X_next,eps;
 
  // setup intersections
  SAFESTATIC vector<Real> tx, ty, tz;
  tx.resize(_gauss.size());
  ty.resize(_gauss.size());
  tz.resize(_gauss.size());

 
  const Real PARAM_BOUND = -5.0;
  int n=0;
  const int NMAX = _gauss.size(); 

  // get the current transform
  const Matrix4& T = get_transform();
 
  // convert the line segment to primitive space
  Vector3 p = T.inverse_mult_point(seg.first);
  Vector3 q = T.inverse_mult_point(seg.second);
  Vector3 d = q - p;

  // starting point of line -> q vector
  // Find x1,y1,z1 ->q vector
  x1=q[0];
  y1=q[1];
  z1=q[2];
 
   // seg.second is p vector
  Real tempX,tempY,tempZ,tempMax;
  tempX=p[0];
  tempY=p[1];
  tempMax = -std::numeric_limits<Real>::max();
  int num = 0;

/*
  // find max(z) of gaussians (we need to find out which Gaussian is tallest
  // at the first point in the line segment)
  for(unsigned i=0; i<_gauss.size();i++){

    // compute some constants
    const Real CTH = std::cos(_gauss[i].th);
    const Real STH = std::sin(_gauss[i].th); 
    const Real S2TH = std::sin((Real) 2.0 * _gauss[i].th); 

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
*/

  //*************************************************************
  // check whether the first point in the line segment is already inside
  if(p[2] <= tempMax)
  {
    t = (Real) 0.0;
    isect = p;

    // compute some constants
    const Real CTH = std::cos(_gauss[num].th);
    const Real STH = std::sin(_gauss[num].th); 
    const Real S2TH = std::sin((Real) 2.0 * _gauss[num].th); 

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
      isect = T.mult_point(isect);
      normal = T.mult_vector(normal);

      return true;
    }
    if(a==0 && b==0) // no tangent
    {
      isect[2]=tempMax;      
      // cout<<"X= "<<isect[0] <<" Y= "<<isect[1] <<" Z= "<<isect[2] <<endl;

      // transform intersection point and normal
      isect = T.mult_point(isect);
      normal = T.mult_vector(normal);

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
    Real tempZ,tempdiffZ;
    Real tempMax = -std::numeric_limits<Real>::max();
    num = 0;
    // find Max F(x,y)
    for(int i=0;i<_gauss.size();i++){

      // compute some constants
      const Real CTH = std::cos(_gauss[i].th);
      const Real STH = std::sin(_gauss[i].th); 
      const Real S2TH = std::sin((Real) 2.0 * _gauss[i].th); 

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
    const Real CTH = std::cos(_gauss[num].th);
    const Real STH = std::sin(_gauss[num].th); 
    const Real S2TH = std::sin((Real) 2.0 * _gauss[num].th); 

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
  const Real CTH = std::cos(_gauss[num].th);
  const Real STH = std::sin(_gauss[num].th); 
  const Real S2TH = std::sin((Real) 2.0 * _gauss[num].th); 

  normal[0] = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(Y - _gauss[num].y0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0))*((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(2*X - 2*_gauss[num].x0));
  normal[1] = _gauss[num].A*std::exp((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(Y - _gauss[num].y0)*(Y - _gauss[num].y0) - (CTH*CTH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) + STH*STH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0)*(X - _gauss[num].x0))*((S2TH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x) - S2TH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y))*(X - _gauss[num].x0) - (CTH*CTH/(2*_gauss[num].sigma_y*_gauss[num].sigma_y) + STH*STH/(2*_gauss[num].sigma_x*_gauss[num].sigma_x))*(2*Y - 2*_gauss[num].y0));
  normal[2] = +1;

  // normalize the normal
  normal.normalize();

  // transform intersection points and normals
  isect = T.mult_point(isect);
  normal = T.mult_vector(normal);

  /*
  cout <<"X= "<<isect[0]<<" Y= "<<isect[1]<<" Z= "<<isect[2]<<endl;
  cout <<"normal[0]= "<<normal[0]<<" normal[1]= "<<normal[1]<<" normal[z]= "<<normal[z]<<endl;
  */

  return true;
}
