/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifdef USE_OSG
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#endif
#include <Moby/XMLTree.h>
#include <Moby/OBB.h>
#include <Moby/Constants.h>
#include <Moby/BoxPrimitive.h>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast; 
using std::list;
using std::vector;
using std::make_pair;
using std::pair;
using std::queue;
using std::endl;

/// Constructs a unit cube centered at the origin
BoxPrimitive::BoxPrimitive()
{
  _xlen = 1;
  _ylen = 1;
  _zlen = 1;
  _edge_sample_length = std::numeric_limits<double>::max();
  calc_mass_properties();
}

/// Constructs a cube of the specified size
BoxPrimitive::BoxPrimitive(double xlen, double ylen, double zlen)
{
  _xlen = xlen;
  _ylen = ylen;
  _zlen = zlen;
  _edge_sample_length = std::numeric_limits<double>::max();
  calc_mass_properties();
}

/// Constructs a unit cube transformed by the given matrix
BoxPrimitive::BoxPrimitive(const Pose3d& T) : Primitive(T)
{
  _xlen = 1;
  _ylen = 1;
  _zlen = 1;
  _edge_sample_length = std::numeric_limits<double>::max();
  calc_mass_properties();
}  

/// Constructs a cube of the specified size transformed by the given matrix
BoxPrimitive::BoxPrimitive(double xlen, double ylen, double zlen, const Pose3d& T) : Primitive(T)
{
  _xlen = xlen;
  _ylen = ylen;
  _zlen = zlen;
  _edge_sample_length = std::numeric_limits<double>::max();
  calc_mass_properties();
}

/// Gets a sub-mesh for the primitive
const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& BoxPrimitive::get_sub_mesh(BVPtr bv)
{
  if (!_smesh.first)
    get_mesh(); 
  return _smesh;
}

/// Sets the intersection tolerance
void BoxPrimitive::set_intersection_tolerance(double tol)
{
  Primitive::set_intersection_tolerance(tol);

  // vertices are no longer valid
  _vertices = shared_ptr<vector<Point3d> >();
  _invalidated = true;
}

/// Sets the size of this box
/**
 * \note forces recomputation of the mesh
 */
void BoxPrimitive::set_size(double xlen, double ylen, double zlen)
{
  _xlen = xlen;
  _ylen = ylen;
  _zlen = zlen;
  if (_xlen < 0.0 && _ylen < 0.0 && _zlen < 0.0)
    throw std::runtime_error("Attempting to pass negative box dimensions to BoxPrimitive::set_size()");

  // mesh, vertices are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices = shared_ptr<vector<Point3d> >();
  _smesh = pair<shared_ptr<IndexedTriArray>, list<unsigned> >();
  _invalidated = true;

  // recalculate the mass properties
  calc_mass_properties();

  // need to update visualization
  update_visualization();
}

/// Sets the edge sample length for this box
void BoxPrimitive::set_edge_sample_length(double len)
{
  _edge_sample_length = len;

  // vertices are no longer valid
  _vertices = shared_ptr<vector<Point3d> >();
  _invalidated = true;
}

/// Transforms the primitive
void BoxPrimitive::set_pose(const Pose3d& p)
{
  // convert p to a shared pointer
  shared_ptr<Pose3d> x(new Pose3d(p));

  // determine the transformation from the global frame to the old pose
  Transform3d cTg = Pose3d::calc_relative_pose(GLOBAL, _F);

  // determine the transformation from the old to the new pose
  Transform3d xTc = Pose3d::calc_relative_pose(_F, x);

  // determine the transformation from the new pose to the global frame 
  Transform3d gTx = Pose3d::calc_relative_pose(x, GLOBAL);

  // compute the transformation
  Transform3d T = gTx * xTc * cTg;

  // go ahead and set the new transform
  Primitive::set_pose(p);

  // transform the mesh
  if (_mesh)
  {
    _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(_mesh->transform(T)));
    _smesh.first = _mesh;
  }

  // invalidate this primitive (in case it is part of a CSG)
  _invalidated = true;

  // transform vertices
  if (_vertices)
    for (unsigned i=0; i< _vertices->size(); i++)
      (*_vertices)[i] = T.transform_point((*_vertices)[i]);

  // recalculate the mass properties
  calc_mass_properties();
}

/// Gets the set of vertices for the BoxPrimitive (constructing, if necessary)
void BoxPrimitive::get_vertices(BVPtr bv, vector<const Point3d*>& vertices)
{
  if (!_vertices)
  {
    if (_xlen <= 0.0 || _ylen <= 0.0 || _zlen <= 0.0)
    {
      vertices.clear(); 
      return;
    }

    // get the transform from the primitive to global coordinates 
    shared_ptr<const Pose3d> P = get_pose();
    Transform3d T = Pose3d::calc_relative_pose(P, GLOBAL);

    // determine the vertices in the mesh
    _vertices = shared_ptr<vector<Point3d> >(new vector<Point3d>());

    // setup half-lengths
    const double XLEN = _xlen*(double) 0.5 + _intersection_tolerance;
    const double YLEN = _ylen*(double) 0.5 + _intersection_tolerance;
    const double ZLEN = _zlen*(double) 0.5 + _intersection_tolerance;

    // add the vertices 
    _vertices->push_back(T.transform_point(Point3d(XLEN,YLEN,ZLEN,P)));
    _vertices->push_back(T.transform_point(Point3d(XLEN,YLEN,-ZLEN,P)));
    _vertices->push_back(T.transform_point(Point3d(XLEN,-YLEN,ZLEN,P)));
    _vertices->push_back(T.transform_point(Point3d(XLEN,-YLEN,-ZLEN,P)));
    _vertices->push_back(T.transform_point(Point3d(-XLEN,YLEN,ZLEN,P)));
    _vertices->push_back(T.transform_point(Point3d(-XLEN,YLEN,-ZLEN,P)));
    _vertices->push_back(T.transform_point(Point3d(-XLEN,-YLEN,ZLEN,P)));
    _vertices->push_back(T.transform_point(Point3d(-XLEN,-YLEN,-ZLEN,P)));
    
    // now we want to add vertices by subdividing edges
    // note: these edges come from facets in get_mesh()
    queue<sorted_pair<unsigned> > q;
    q.push(make_sorted_pair(3, 1));
    q.push(make_sorted_pair(0, 2));
    q.push(make_sorted_pair(2, 3));
    q.push(make_sorted_pair(0, 1));
    q.push(make_sorted_pair(1, 5));
    q.push(make_sorted_pair(5, 4));
    q.push(make_sorted_pair(4, 0));
    q.push(make_sorted_pair(6, 2));
    q.push(make_sorted_pair(4, 6));
    q.push(make_sorted_pair(7, 5));
    q.push(make_sorted_pair(3, 7));
    q.push(make_sorted_pair(6, 7));

    // process until queue is empty
    while (!q.empty())
    {
      // get the edge off of the front of the queue and get the corresponding
      // vertices
      sorted_pair<unsigned> sp = q.front();
      q.pop();
      const Point3d& v1 = (*_vertices)[sp.first];
      const Point3d& v2 = (*_vertices)[sp.second];

      // if edge is sufficiently small, don't need to do any more processing
      if ((v1-v2).norm() <= _edge_sample_length)
        continue;

      // subdivide the edge by first creating a new vertex
      unsigned vidx = _vertices->size();
      _vertices->push_back((v1+v2) * (double) 0.5);
      q.push(make_sorted_pair(sp.first, vidx)); 
      q.push(make_sorted_pair(sp.second, vidx)); 
    }
  }

  // copy the addresses of the computed vertices into 'vertices' 
  vertices.resize(_vertices->size());
  for (unsigned i=0; i< _vertices->size(); i++)
    vertices[i] = &(*_vertices)[i];
}

/// Recomputes the triangle mesh for the BoxPrimitive
shared_ptr<const IndexedTriArray> BoxPrimitive::get_mesh()
{
  // create the mesh if necessary
  if (!_mesh)
  {
    const unsigned BOX_VERTS = 8, BOX_FACETS = 12;

    // setup lengths
    double XLEN = _xlen * (double) 0.5;
    double YLEN = _ylen * (double) 0.5;
    double ZLEN = _zlen * (double) 0.5;

    // need eight vertices
    Point3d verts[BOX_VERTS];
    verts[0] = Point3d(XLEN,YLEN,ZLEN);
    verts[1] = Point3d(XLEN,YLEN,-ZLEN);
    verts[2] = Point3d(XLEN,-YLEN,ZLEN);
    verts[3] = Point3d(XLEN,-YLEN,-ZLEN);
    verts[4] = Point3d(-XLEN,YLEN,ZLEN);
    verts[5] = Point3d(-XLEN,YLEN,-ZLEN);
    verts[6] = Point3d(-XLEN,-YLEN,ZLEN);
    verts[7] = Point3d(-XLEN,-YLEN,-ZLEN);

    // create 12 new triangles, making sure to do so ccw
    IndexedTri facets[BOX_FACETS];
    facets[0] = IndexedTri(3, 1, 0);
    facets[1] = IndexedTri(0, 2, 3);
    facets[2] = IndexedTri(0, 1, 5);
    facets[3] = IndexedTri(5, 4, 0);
    facets[4] = IndexedTri(6, 2, 0);
    facets[5] = IndexedTri(0, 4, 6);
    facets[6] = IndexedTri(7, 5, 1);
    facets[7] = IndexedTri(1, 3, 7);
    facets[8] = IndexedTri(6, 7, 3);
    facets[9] = IndexedTri(3, 2, 6);
    facets[10] = IndexedTri(7, 6, 4);
    facets[11] = IndexedTri(4, 5, 7);
  
    // setup the triangle mesh
    _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(verts, verts+BOX_VERTS, facets, facets+BOX_FACETS));

    // setup sub mesh (it will be just the standard mesh)
    list<unsigned> all_tris;
    for (unsigned i=0; i< _mesh->num_tris(); i++)
      all_tris.push_back(i);
    _smesh = make_pair(_mesh, all_tris);
  }

  return _mesh;
}

/// Creates the visualization for this primitive
osg::Node* BoxPrimitive::create_visualization()
{
  #ifdef USE_OSG
  const double HALF = (double) 0.5;
  osg::Box* box = new osg::Box;
  osg::Vec3 half_lens(_xlen * HALF, _ylen * HALF, _zlen * HALF);
  box->setHalfLengths(half_lens);
  osg::Geode* geode = new osg::Geode;
  geode->addDrawable(new osg::ShapeDrawable(box));
  return geode;
  #else
  return NULL;
  #endif
}

/// Implements Base::load_from_xml() for serialization
void BoxPrimitive::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // verify that the node type is BoxPrimitive
  assert(strcasecmp(node->name.c_str(), "Box") == 0);

  // load the parent data
  Primitive::load_from_xml(node, id_map);

  // get the length attributes, if specified
  XMLAttrib* xlen_attr = node->get_attrib("xlen");
  XMLAttrib* ylen_attr = node->get_attrib("ylen");
  XMLAttrib* zlen_attr = node->get_attrib("zlen");

  // set lengths to zero initially
  double xlen = 0, ylen = 0, zlen = 0;

  // get the lengths
  if (xlen_attr) xlen = xlen_attr->get_real_value();
  if (ylen_attr) ylen = ylen_attr->get_real_value();
  if (zlen_attr) zlen = zlen_attr->get_real_value();

  // read in the edge sample length
  XMLAttrib* esl_attr = node->get_attrib("edge-sample-length");
  if (esl_attr)
    set_edge_sample_length(esl_attr->get_real_value());

  // set the lengths of the box
  set_size(xlen, ylen, zlen);

  // recompute mass properties
  calc_mass_properties();
}

/// Implements Base::save_to_xml() for serialization
void BoxPrimitive::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data
  Primitive::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "Box";

  // save the lengths
  node->attribs.insert(XMLAttrib("xlen", _xlen));
  node->attribs.insert(XMLAttrib("ylen", _ylen));
  node->attribs.insert(XMLAttrib("zlen", _zlen));

  // save the edge sample length
  node->attribs.insert(XMLAttrib("edge-sample-length", _edge_sample_length));
}

/// Calculates mass properties for this primitive
void BoxPrimitive::calc_mass_properties()
{
  // get the current transform
  shared_ptr<const Pose3d> T = get_pose();

  // compute the mass if necessary 
  if (_density)
  {
    const double volume = _xlen * _ylen * _zlen;
    _J.m = *_density * volume;
  } 

  // compute some constants
  const double XSQ = _xlen * _xlen;
  const double YSQ = _ylen * _ylen;
  const double ZSQ = _zlen * _zlen;
  const double M = _J.m/12.0;

  // compute the inertia matrix
  _J.J = Matrix3d(M*(YSQ+ZSQ), 0, 0, 0, M*(XSQ+ZSQ), 0, 0, 0, M*(XSQ+YSQ));
}

/// Gets the bounding volume for this plane
BVPtr BoxPrimitive::get_BVH_root()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // box not applicable for deformable bodies 
  if (is_deformable())
    throw std::runtime_error("BoxPrimitive::get_BVH_root() - primitive unusable for deformable bodies!");

  // create the bounding box, if necessary
  if (!_obb)
    _obb = shared_ptr<OBB>(new OBB);

  // get the transform
  shared_ptr<const Pose3d> T = get_pose();

  // setup the center of the OBB 
  _obb->center = Point3d(T->x, T);
  
  // setup the orientation of the OBB
  _obb->R = T->q;

  // setup OBB half-lengths
  _obb->l[X] = _xlen * (double) 0.5;
  _obb->l[Y] = _ylen * (double) 0.5;
  _obb->l[Z] = _zlen * (double) 0.5;

  return _obb;
}

/// Tests whether a point is inside or on the box
bool BoxPrimitive::point_inside(BVPtr bv, const Point3d& point, Vector3d& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // form the transformation from box space to query space (and back) 
  shared_ptr<const Pose3d> P = get_pose();
  Transform3d T = Pose3d::calc_relative_pose(point.pose, P);

  // setup lengths
  Vector3d l(P);
  l[X] = _xlen * (double) 0.5;
  l[Y] = _ylen * (double) 0.5;
  l[Z] = _zlen * (double) 0.5;

  // convert the point to primitive space
  Point3d p = T.transform_point(point);

  FILE_LOG(LOG_COLDET) << "BoxPrimitive::point_inside() entered" << endl; 
  FILE_LOG(LOG_COLDET) << "  -- querying point " << p << " and box: " << l << endl;

  // check whether p is inside the box 
  if (p[X] < -l[X] || p[X] > l[X] ||
      p[Y] < -l[Y] || p[Y] > l[Y] ||
      p[Z] < -l[Z] || p[Z] > l[Z])
  {
    FILE_LOG(LOG_COLDET) << "  ** point is outside" << endl;
    return false;
  }

  // p is inside the box; determine the normal
  double absPX = l[X] - std::fabs(p[X]);
  double absPY = l[Y] - std::fabs(p[Y]);
  double absPZ = l[Z] - std::fabs(p[Z]);
  if (absPZ < absPX - NEAR_ZERO && absPZ < absPY - NEAR_ZERO)
  {
    if (p[Z] < (double) 0.0)
      normal = T.inverse_transform_vector(Vector3d(0,0,-1,P));
    else
      normal = T.inverse_transform_vector(Vector3d(0,0,1,P));
  }
  else if (absPY < absPZ - NEAR_ZERO && absPY < absPX - NEAR_ZERO)
  {
    if (p[Y] < (double) 0.0)
      normal = T.inverse_transform_vector(Vector3d(0,-1,0,P));
    else
      normal = T.inverse_transform_vector(Vector3d(0,1,0,P));
  }
  else if (absPX < absPY - NEAR_ZERO && absPX < absPZ - NEAR_ZERO)
  {
    if (p[X] < (double) 0.0)
      normal = T.inverse_transform_vector(Vector3d(-1,0,0,P));
    else
      normal = T.inverse_transform_vector(Vector3d(1,0,0,P));
  }
  else
  {
    // degenerate normal
    normal .set_zero();
    normal.pose = p.pose;
  }

  FILE_LOG(LOG_COLDET) << "  ** point is inside" << endl;
  return true;
}

/// Computes the intersection of a line segment with the box
/**
 * \note for line segments that are partially or fully inside the box, the
 *       method only returns intersection if the second endpoint of the segment
 *       is farther inside than the first
 */
bool BoxPrimitive::intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  double tmin = (double) 0.0;
  double tmax = (double) 2.0;

  // form a transform between segment space and box space 
  shared_ptr<const Pose3d> P = get_pose();
  Transform3d T = Pose3d::calc_relative_pose(seg.first.pose, P);

  // setup lengths
  Vector3d l(P);
  l[X] = _xlen * (double) 0.5;
  l[Y] = _ylen * (double) 0.5;
  l[Z] = _zlen * (double) 0.5;

  // convert the line segment to primitive space
  Point3d p = T.transform_point(seg.first);
  Point3d q = T.transform_point(seg.second);
  Vector3d d = q - p;

  FILE_LOG(LOG_COLDET) << "BoxPrimitive::intersect_seg() entered" << endl; 
  FILE_LOG(LOG_COLDET) << "  -- checking intersection between line segment " << p << " / " << (d + p) << " and box: " << l << endl;

  // check whether p is already inside the box
  if (p[X] >= -l[X] - NEAR_ZERO && p[X] <= l[X] + NEAR_ZERO &&
      p[Y] >= -l[Y] - NEAR_ZERO && p[Y] <= l[Y] + NEAR_ZERO &&
      p[Z] >= -l[Z] - NEAR_ZERO && p[Z] <= l[Z] + NEAR_ZERO)
  {
    // set time and point of intersection
    t = (double) 0.0;
    isect = seg.first;

    // determine normal
    double absPX = l[X] - std::fabs(p[X]);
    double absPY = l[Y] - std::fabs(p[Y]);
    double absPZ = l[Z] - std::fabs(p[Z]);
    if (absPZ < absPX && absPZ < absPY)
    {
      if (p[Z] < (double) 0.0)
        normal = T.inverse_transform_vector(Vector3d(0,0,-1,P));
      else
        normal = T.inverse_transform_vector(Vector3d(0,0,1,P));
    }
    else if (absPY < absPZ && absPY < absPX)
    {
      if (p[Y] < (double) 0.0)
        normal = T.inverse_transform_vector(Vector3d(0,-1,0,P));
      else
        normal = T.inverse_transform_vector(Vector3d(0,1,0,P));
    }
    else
    {
      if (p[X] < (double) 0.0)
        normal = T.inverse_transform_vector(Vector3d(-1,0,0,P));
      else
        normal = T.inverse_transform_vector(Vector3d(1,0,0,P));
    }

    FILE_LOG(LOG_COLDET) << " -- point is already inside the box..." << endl;
    return true;
  }

  // for all three slabs
  for (unsigned i=0; i< 3; i++)
  {
    if (std::fabs(d[i]) < NEAR_ZERO)
    {
      // line is parallel to slab; no hit if origin not within slab
      if (p[i] < -l[i] || p[i] > l[i])
      {
        FILE_LOG(LOG_COLDET) << "  -- seg parallel to slab " << i << " and origin not w/in slab = no intersection" << endl; 
        FILE_LOG(LOG_COLDET) << "    amin: " << -l[i] << "  amax: " << l[i] << "  p: " << p << endl;
        return false;
      }
    }
    else
    {
      // compute intersection value of line with near and far plane of slab
      double ood = 1.0 / d[i];
      double t1 = (-l[i] - p[i]) * ood;
      double t2 = (l[i] - p[i]) * ood;

      FILE_LOG(LOG_COLDET) << "  ood: " << ood << " t1: " << t1 << "  t2: " << t2 << endl; 

      // make t1 be the intersection with the near plane, t2 with far plane
      if (t1 > t2)
        std::swap(t1, t2);

      // compute the intersection of slab intersection intervals
      tmin = std::max(tmin, t1);
      tmax = std::min(tmax, t2);

      // exit with no collision as soon as slab intersection becomes empty
      if (tmin > tmax + NEAR_ZERO || tmin > (double) 1.0)
      {
        FILE_LOG(LOG_COLDET) << "  tmin (" << tmin << ") > 1.0 or tmax (" << tmax << ") -- seg and box do not intersect" << endl;

        return false;
      }
    }
  }

  // ray intersects all three slabs
  // determine normal; if it is degenerate, report no intersection
  isect = p + d * tmin;
  if (!determine_normal_abs(l, isect, normal))
    return false;

  // transform the normal back out of box space
  normal = T.inverse_transform_vector(normal);
  assert(std::fabs(normal.norm() - (double) 1.0) < NEAR_ZERO);

  // transform intersection point out of box space 
  isect = T.inverse_transform_point(isect);

  FILE_LOG(LOG_COLDET) << "BoxPrimitive::intersects() - seg and box intersect; first intersection: " << tmin << "(" << isect << ")" << endl; 
  FILE_LOG(LOG_COLDET) << "BoxPrimitive::intersects() exiting" << endl; 

  // store the intersection parameter
  t = tmin;
  if (t < 0.0)
    t = 0.0;

  return true;
}

/// Determines the normal to a point inside or on the box
/**
 * \return <b>true</b> if the normal is valid, <b>false</b> if degenerate
 */
void BoxPrimitive::determine_normal(const Vector3d& lengths, const Point3d& p, Vector3d& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const unsigned NFACES = 6;
  pair<double, BoxPrimitive::FaceID> distances[NFACES];

  // get the pose of the primitive
  shared_ptr<const Pose3d> P = get_pose();

  distances[0] = make_pair(p[X] - lengths[X], ePOSX);
  distances[1] = make_pair(p[X] + lengths[X], eNEGX);
  distances[2] = make_pair(p[Y] - lengths[Y], ePOSY);
  distances[3] = make_pair(p[Y] + lengths[Y], eNEGY);
  distances[4] = make_pair(p[Z] - lengths[Z], ePOSZ);
  distances[5] = make_pair(p[Z] + lengths[Z], eNEGZ);
  std::sort(distances, distances+NFACES);

  // if we have an equal distance from (at least) two planes, we'll just use
  // the normal from the first plane regardless..  Hopefully, it will still
  // work ok (normal is degenerate in this case)
  switch (distances[0].second)
  {
    case ePOSX:  normal = Vector3d(1,0,0,P); break;
    case eNEGX:  normal = Vector3d(-1,0,0,P); break;
    case ePOSY:  normal = Vector3d(0,1,0,P); break;
    case eNEGY:  normal = Vector3d(0,-1,0,P); break;
    case ePOSZ:  normal = Vector3d(0,0,1,P); break;
    case eNEGZ:  normal = Vector3d(0,0,-1,P); break;
  }
}

/// Determines the normal to a point on the box
/**
 * \return <b>true</b> if the normal is valid, <b>false</b> if degenerate
 */
bool BoxPrimitive::determine_normal_abs(const Vector3d& lengths, const Point3d& p, Vector3d& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const unsigned NFACES = 6;
  pair<double, BoxPrimitive::FaceID> distances[NFACES];

   // get the pose of the primitive
  shared_ptr<const Pose3d> P = get_pose();

  distances[0] = make_pair(std::fabs(p[X] - lengths[X]), ePOSX);
  distances[1] = make_pair(std::fabs(p[X] + lengths[X]), eNEGX);
  distances[2] = make_pair(std::fabs(p[Y] - lengths[Y]), ePOSY);
  distances[3] = make_pair(std::fabs(p[Y] + lengths[Y]), eNEGY);
  distances[4] = make_pair(std::fabs(p[Z] - lengths[Z]), ePOSZ);
  distances[5] = make_pair(std::fabs(p[Z] + lengths[Z]), eNEGZ);
  std::sort(distances, distances+NFACES);

  // if we have an equal distance from (at least) two planes, we'll just use
  // the normal from the first plane regardless..  Hopefully, it will still
  // work ok (normal is degenerate in this case)
  switch (distances[0].second)
  {
    case ePOSX:  normal = Vector3d(1,0,0,P); break;
    case eNEGX:  normal = Vector3d(-1,0,0,P); break;
    case ePOSY:  normal = Vector3d(0,1,0,P); break;
    case eNEGY:  normal = Vector3d(0,-1,0,P); break;
    case ePOSZ:  normal = Vector3d(0,0,1,P); break;
    case eNEGZ:  normal = Vector3d(0,0,-1,P); break;
  }

  return !CompGeom::rel_equal(distances[0].second, distances[1].second, NEAR_ZERO);
}

