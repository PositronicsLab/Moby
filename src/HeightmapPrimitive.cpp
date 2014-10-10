/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifdef USE_OSG
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Material>
#include <osg/LightModel>
#endif
#include <boost/algorithm/minmax_element.hpp>
#include <Moby/Constants.h>
#include <Moby/CompGeom.h>
#include <Moby/sorted_pair>
#include <Moby/XMLTree.h>
#include <Moby/BoundingSphere.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/HeightmapPrimitive.h>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr; 
using std::map;
using std::vector;
using std::list;
using std::pair;
using boost::dynamic_pointer_cast;
using boost::const_pointer_cast;
using std::make_pair;

/// Initializes the heightmap primitive
HeightmapPrimitive::HeightmapPrimitive()
{
  _width = _depth = 0.0;
}

/// Initializes the heightmap primitive
HeightmapPrimitive::HeightmapPrimitive(const Ravelin::Pose3d& T) : Primitive(T)
{
  _width = _depth = 0.0;
}

/// Gets the supporting point
Point3d HeightmapPrimitive::get_supporting_point(const Vector3d& d) const 
{
  throw std::runtime_error("HeightmapPrimitive::get_supporting_point(.) - primitive is not convex!");
  return Point3d(0,0,0,GLOBAL);
}

/// Computes the signed distance of the given point from this primitive
double HeightmapPrimitive::calc_signed_dist(const Point3d& p) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end());

  // compute the height
  return calc_height(p);
}

/// Gets the BVH root for the heightmap
BVPtr HeightmapPrimitive::get_BVH_root(CollisionGeometryPtr geom)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the pointer to the geometry
  shared_ptr<OBB>& obb = _obbs[geom];

  // if the OBB hasn't been created, create it and initialize it
  if (!obb)
  {
    // create the OBB 
    obb = shared_ptr<OBB>(new OBB);
    obb->geom = geom;

    // get the pose for the geometry
    shared_ptr<const Pose3d> P = get_pose(geom);

    // setup the rotation matrix for the OBB
    obb->R.set_identity();

    // setup the translation for the OBB
    obb->center.set_zero(P);

    // get the maximum and minimum height
    pair<ColumnIteratord, ColumnIteratord> mmax = boost::minmax_element(_heights.column_iterator_begin(), _heights.column_iterator_end());
    double maxy = *mmax.second;
    double miny = *mmax.first;
    obb->center[Y] = (maxy+miny)*0.5;

    // setup obb dimensions
    obb->l.pose = P;
    obb->l[X] = _width*0.5;
    obb->l[Y] = (maxy-miny)*0.5;
    obb->l[Z] = _depth*0.5;
  }

  return obb;
}

/// Gets the mesh of the heightmap
shared_ptr<const IndexedTriArray> HeightmapPrimitive::get_mesh(shared_ptr<const Pose3d> P)
{
  // TODO: implement this with Delaunay triangulation
  assert(false);
  return shared_ptr<const IndexedTriArray>();
}

/// Gets the vertices of the heightmap that could intersect with a given bounding volume
void HeightmapPrimitive::get_vertices(BVPtr bv, shared_ptr<const Pose3d> P, vector<Point3d>& vertices) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(P)) != _poses.end());
  const unsigned X = 0, Z = 2;

  // clear the vector of vertices
  vertices.clear();

  // get the corners of the bounding box in this frame
  Point3d bv_lo = bv->get_lower_bounds();
  Point3d bv_hi = bv->get_lower_bounds();

  // get the lower i and j indices
  unsigned lowi = (unsigned) ((bv_lo[X]+_width*0.5)*(_heights.rows()-1)/_width);
  unsigned lowj = (unsigned) ((bv_lo[Z]+_depth*0.5)*(_heights.columns()-1)/_depth);

  // get the upper i and j indices
  unsigned upi = (unsigned) ((bv_hi[X]+_width*0.5)*(_heights.rows()-1)/_width)+1;
  unsigned upj = (unsigned) ((bv_hi[Z]+_depth*0.5)*(_heights.columns()-1)/_depth)+1;

  // iterate over all points
  for (unsigned i=lowi; i<= upi; i++)
    for (unsigned j=lowj; j< upj; j++)
    {
      double x = -_width*0.5+_width*i/(_heights.rows()-1);
      double z = -_depth*0.5+_depth*j/(_heights.columns()-1);
      vertices.push_back(Point3d(x, _heights(i,j), z, P));
    }
}

/// Gets the vertices of the heightmap
void HeightmapPrimitive::get_vertices(shared_ptr<const Pose3d> P, vector<Point3d>& vertices) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(P)) != _poses.end() || P == get_pose());

  // clear the vector of vertices
  vertices.clear();

  // iterate over all points
  for (unsigned i=0; i< _heights.rows(); i++)
    for (unsigned j=0; j< _heights.columns(); j++)
    {
      double x = -_width*0.5+_width*i/(_heights.rows()-1);
      double z = -_depth*0.5+_depth*j/(_heights.columns()-1);
      vertices.push_back(Point3d(x, _heights(i,j), z, P));
    }
}

/// Computes the height at a particular point
double HeightmapPrimitive::calc_height(const Point3d& p) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end());
  const unsigned X = 0, Y = 1, Z = 2;

  // get the X and Z query points
  const unsigned qx = p[X];
  const unsigned qz = p[Z];

  // determine the indices 
  unsigned i = (unsigned) ((qx+_width*0.5)*(_heights.rows()-1)/_width);
  unsigned j = (unsigned) ((qz+_depth*0.5)*(_heights.columns()-1)/_depth);
  assert(i < _heights.rows());
  assert(j < _heights.columns());

  // setup inputs
  double x0 = -_width*0.5+_width*i/(_heights.rows()-1);
  double z0 = -_depth*0.5+_depth*j/(_heights.columns()-1);
  double x1 = -_width*0.5+_width*(i+1)/(_heights.rows()-1);
  double z1 = -_depth*0.5+_depth*(j+1)/(_heights.columns()-1);

  // compute s and t
  double s = qx/(x1 - x0);
  double t = qz/(z1 - z0); 

  // get four height values
  const double f00 = _heights(i,j);
  const double f10 = _heights(i+1,j);
  const double f01 = _heights(i,j+1);
  const double f11 = _heights(i+1,j+1);

  return p[Y] - f00*(1.0-s)*(1.0-t) + f10*s*(1.0-t) + f01*(1.0-s)*t + f11*s*t;
}

/// Computes the gradient at a particular point
void HeightmapPrimitive::calc_gradient(const Point3d& p, double& gx, double& gz) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end());
  const unsigned X = 0, Z = 2;

  // get the X and Z query points
  const unsigned qx = p[X];
  const unsigned qz = p[Z];

  // determine the indices 
  unsigned i = (unsigned) (qx*(_heights.rows()-1)/_width + _width*0.5);
  unsigned j = (unsigned) (qz*(_heights.columns()-1)/_depth + _depth*0.5);
  assert(i < _heights.rows());
  assert(j < _heights.columns());

  // setup inputs
  double x0 = -_width*0.5+_width*i/(_heights.rows()-1);
  double z0 = -_depth*0.5+_depth*j/(_heights.columns()-1);
  double x1 = -_width*0.5+_width*(i+1)/(_heights.rows()-1);
  double z1 = -_depth*0.5+_depth*(j+1)/(_heights.columns()-1);

  // compute s and t
  double s = qx/(x1 - x0);
  double t = qz/(z1 - z0); 

  // get four height values
  const double f00 = _heights(i,j);
  const double f10 = _heights(i+1,j);
  const double f01 = _heights(i,j+1);
  const double f11 = _heights(i+1,j+1);

  // compute the x gradient
  gx = f10*(1.0-t) + f11*t;

  // compute the z gradient
  gz = f01*(1.0-s) + f11*s;
}

static void perturb_color(float color[3])
{
  const float EPS = 1e-2;

  for (unsigned i=0; i< 3; i++)
  {
//    color[i] += (float) rand() / RAND_MAX * EPS - 2*EPS;
    color[i] = (float) rand() / RAND_MAX;
    if (color[i] < 0.0f)
      color[i] = 0.0f;
    else if (color[i] > 1.0f)
      color[i] = 1.0f;
  }
}

/// Computes the OSG visualization
osg::Node* HeightmapPrimitive::create_visualization()
{
  #ifdef USE_OSG
  const unsigned X = 0, Y = 1, Z = 2;

  // get the pose and compute transform from the global frame to it 
  shared_ptr<const Pose3d> P = get_pose();
  Transform3d T = Pose3d::calc_relative_pose(P, GLOBAL);

  // create necessary OSG elements for visualization
  osg::Group* group = new osg::Group;

  // get the set of vertices
  std::vector<Point3d> verts;
  get_vertices(P, verts);

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
  osg::Material* mat = new osg::Material;
  mat->setColorMode(osg::Material::DIFFUSE);
  mat->setDiffuse(osg::Material::FRONT, osg::Vec4(1.0f, 1.0f, 1.0f, 1.0f));
  subgroup->getOrCreateStateSet()->setAttribute(mat);

  // create the vertex array
  osg::Vec3Array* varray = new osg::Vec3Array(verts.size()); 
  for (unsigned i=0; i< verts.size(); i++)
    (*varray)[i] = osg::Vec3((float) verts[i][X], (float) verts[i][Y], (float) verts[i][Z]);
  geom->setVertexArray(varray);

  // get the maximum and minimum intensity
  double min_value = std::numeric_limits<double>::max();
  double max_value = -std::numeric_limits<double>::max();
  for (unsigned i=0; i< verts.size(); i++)
  {
    min_value = std::min(min_value, verts[i][Y]);
    max_value = std::max(max_value, verts[i][Y]);
  }

  // create the color array (for vertices)
  osg::Vec4Array* vcolors = new osg::Vec4Array(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
  {
    const float intensity = (verts[i][Y]-min_value)/(max_value-min_value)*0.9;
    (*vcolors)[i] = osg::Vec4(intensity+0.1f, 0.4f, (float) rand()/RAND_MAX*0.1 + 0.2, 1.0f); 
  }
  geom->setColorArray(vcolors);
  geom->setColorBinding(osg::Geometry::BIND_PER_VERTEX);

  // create the faces - we're going to iterate over every grouping of four
  // points
  const unsigned COLS = _heights.columns();
  for (unsigned i=0; i< _heights.rows()-1; i++)
    for (unsigned j=0; j< _heights.columns()-1; j++)
    {
      // get the four indices
      const unsigned V1 =  i*COLS+j; // i, j
      const unsigned V2 =  (i+1)*COLS+j; // i+1, j
      const unsigned V3 =  i*COLS+j+1; // i, j+1
      const unsigned V4 =  (i+1)*COLS+j+1; // i+1, j+1 

      // create the first facet 
      osg::DrawElementsUInt* face1 = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
      face1->push_back(V2);
      face1->push_back(V1);
      face1->push_back(V3);
      geom->addPrimitiveSet(face1);
      
      // create the second facet 
      osg::DrawElementsUInt* face2 = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
      face2->push_back(V3);
      face2->push_back(V4);
      face2->push_back(V2);
      geom->addPrimitiveSet(face2);
    }

  return group;
  #else
  return NULL;
  #endif 
}

/// Computes the distance from a sphere primitive
double HeightmapPrimitive::calc_signed_dist(shared_ptr<const SpherePrimitive> s, Point3d& pthis, Point3d& ps) const
{
  const unsigned X = 0, Z = 2;
  assert(_poses.find(const_pointer_cast<Pose3d>(pthis.pose)) != _poses.end());
  Point3d ps_prime = ps;

  // compute the transform from the sphere to the heightmap
  Transform3d T = Pose3d::calc_relative_pose(ps.pose, pthis.pose);

  // transform the sphere center to the height map space
  Point3d ps_c(0.0, 0.0, 0.0, ps.pose);
  Point3d ps_c_this = T.transform_point(ps_c);

  // get the lowest point on the sphere (toward the heightmap)
  Vector3d vdir(0.0, -1.0*s->get_radius(), 0.0, pthis.pose);

  // get the lowest point on the sphere
  Point3d sphere_lowest = ps_c_this + vdir; 

  // get the height of the sphere center
  double min_dist = calc_height(sphere_lowest);  
  ps = T.inverse_transform_point(sphere_lowest); 

  // get the corners of the bounding box in this frame
  Point3d bv_lo = ps_c_this;
  Point3d bv_hi = ps_c_this;
  bv_lo[X] -= s->get_radius();
  bv_hi[X] += s->get_radius();
  bv_lo[Z] -= s->get_radius();
  bv_hi[Z] += s->get_radius();

  // get the lower i and j indices
  unsigned lowi = (unsigned) ((bv_lo[X]+_width*0.5)*(_heights.rows()-1)/_width);
  unsigned lowj = (unsigned) ((bv_lo[Z]+_depth*0.5)*(_heights.columns()-1)/_depth);

  // get the upper i and j indices
  unsigned upi = (unsigned) ((bv_hi[X]+_width*0.5)*(_heights.rows()-1)/_width)+1;
  unsigned upj = (unsigned) ((bv_hi[Z]+_depth*0.5)*(_heights.columns()-1)/_depth)+1;

  // iterate over all points in the bounding region
  for (unsigned i=lowi; i<= upi; i++)
    for (unsigned j=lowj; j< upj; j++)
    {
      // compute the point on the heightmap
      double x = -_width*0.5+_width*i/(_heights.rows()-1);
      double z = -_depth*0.5+_depth*j/(_heights.columns()-1);
      Point3d p(x, _heights(i,j), z, pthis.pose);
      Point3d ps_prime = Pose3d::transform_point(ps.pose, p);

      // get the distance from the sphere
      double dist = s->calc_signed_dist(ps_prime);

      // see how the distance compares
      if (dist < min_dist)
      {
        min_dist = dist;
        ps = ps_prime;
      }
    }

  return min_dist;
}

/// Transforms the primitive
void HeightmapPrimitive::set_pose(const Pose3d& p)
{
  // convert p to a shared pointer
  shared_ptr<Pose3d> x(new Pose3d(p));

  // determine the transformation from the old pose to the new one 
  Transform3d T = Pose3d::calc_relative_pose(_F, x);

  // go ahead and set the new transform
  Primitive::set_pose(p);
}

/// Finds the signed distance between the sphere and another primitive
double HeightmapPrimitive::calc_signed_dist(shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const
{
  const unsigned Y = 1;

  // first try heightmap/sphere
  shared_ptr<const SpherePrimitive> spherep = dynamic_pointer_cast<const SpherePrimitive>(p);
  if (spherep)
    return calc_signed_dist(spherep, pthis, pp);

  // get p as non-const
  shared_ptr<Primitive> pnc = const_pointer_cast<Primitive>(p);

  // still here? no specialization; get all vertices from other primitive
  vector<Point3d> verts;
  pnc->get_vertices(pp.pose, verts);
  double mindist = std::numeric_limits<double>::max();
  for (unsigned i=0; i< verts.size(); i++)
  {
     Point3d pt = Pose3d::transform_point(pthis.pose, verts[i]);
     const double HEIGHT = calc_height(pt);
     if (HEIGHT < mindist)
     {
       mindist = HEIGHT;
       pp = verts[i];
       pthis = pt;
       pthis[Y] = -(HEIGHT - pt[Y]);
     }
  }

  return mindist;
}

/// Finds the signed distance betwen the heightmap and a point
double HeightmapPrimitive::calc_dist_and_normal(const Point3d& p, std::vector<Vector3d>& normals) const
{
  // compute the distance
  double d = calc_height(p);

  // setup the normal
  normals.push_back(Vector3d());
  Vector3d& normal = normals.back();

  // setup the normal
  if (d >= 0.0)
  {
    double gx, gz;
    calc_gradient(p, gx, gz);
    normal = Vector3d::normalize(Vector3d(gx, 1, gz, p.pose));
  }
  else
    normal = Vector3d(0.0, 1.0, 0.0, p.pose);

  // compute the distance
  return d;
}

/// Implements Base::load_from_xml() for serialization
void HeightmapPrimitive::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // verify that the node type is heightmap
  assert(strcasecmp(node->name.c_str(), "Heightmap") == 0);

  // load the parent data
  Primitive::load_from_xml(node, id_map);

  // read in the height map
  XMLAttrib* file_attr = node->get_attrib("filename");
  if (file_attr)
  {
    std::ifstream in(file_attr->get_string_value().c_str());
    if (!in.fail())
    {
      unsigned rows, cols;
      in >> rows;
      in >> cols;
      _heights.resize(rows, cols);
      for (unsigned i=0; i< rows; i++)
        for (unsigned j=0; j< cols; j++)
          in >> _heights(i,j);
      in.close();
    }
    else
    {
      std::cerr << "HeightmapPrimitive::load_from_xml() - unable to read heightmap!" << std::endl;
      _heights.set_zero(1,1);
    }
  }

  // read in the width, if specified
  XMLAttrib* width_attr = node->get_attrib("width");
  if (width_attr)
    _width = width_attr->get_unsigned_value();

  // read in the depth, if specified
  XMLAttrib* depth_attr = node->get_attrib("depth");
  if (depth_attr)
    _depth = depth_attr->get_unsigned_value();
}

/// Implements Base::save_to_xml() for serialization
void HeightmapPrimitive::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data
  Primitive::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "Heightmap";

  // save the width 
  node->attribs.insert(XMLAttrib("width", _width));

  // save the depth 
  node->attribs.insert(XMLAttrib("depth", _depth));

  // write out the height map
  const unsigned MAX_DIGITS = 28;
  char buffer[MAX_DIGITS+1];
  sprintf(buffer, "%p", this);
  std::string filename = "heightmap" + std::string(buffer) + ".dat";

  // add the filename as an attribute
  node->attribs.insert(XMLAttrib("filename", filename));

  // write the heightmap
  std::ofstream out(filename.c_str());
  if (out.fail())
  {
    std::cerr << "HeightmapPrimitive::save_to_xml() - unexpectedly unable to write heightmap!" << std::endl;
    return;
  }
  out << _heights.rows() << " " << _heights.columns() << std::endl;
  for (unsigned i=0; i< _heights.rows(); i++)
    for (unsigned j=0; j< _heights.columns(); j++)
      out << _heights(i,j);
  out.close();
}


