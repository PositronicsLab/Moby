/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <iostream>
#ifdef USE_OSG
#include <osg/Geode>
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/MatrixTransform>
#endif
#include <Moby/Constants.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CompGeom.h>
#include <Moby/XMLTree.h>
#include <Moby/OBB.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/GJK.h>
#include <Moby/ConePrimitive.h>

using namespace Ravelin;
using namespace Moby;
using std::map;
using std::pair;
using boost::shared_ptr;
using boost::const_pointer_cast;
using boost::dynamic_pointer_cast;
using std::list;
using std::vector;

/// Constructs a cone centered at the origin, with the longitudinal axis aligned with the y-axis, radius 1.0, height 1.0, 1 ring and 10 circle points 
ConePrimitive::ConePrimitive()
{
  _radius = 1.0;
  _height = 1.0;
  _npoints = 10;
  _nrings = 1;
  calc_mass_properties();
}

/// Constructs a cone along the y-axis with specified radius and height, centered at the origin, with 1 ring and 10 circle points 
ConePrimitive::ConePrimitive(double radius, double height)
{
  if (height < (double) 0.0)
    throw std::runtime_error("Attempting to set negative height in ConePrimitive (constructor)");
  if (radius < (double) 0.0)
    throw std::runtime_error("Attempting to set negative radius in ConePrimitive (constructor)");

  _radius = radius;
  _height = height;
  _npoints = 10;
  _nrings = 1;
  calc_mass_properties();
}

/// Constructs a cone along the y-axis with specified radius and height, centered at the origin, with 1 ring and 10 circle points 
ConePrimitive::ConePrimitive(double radius, double height, const Pose3d& T) : Primitive(T)
{
  if (height < (double) 0.0)
    throw std::runtime_error("Attempting to set negative height in ConePrimitive (constructor)");
  if (radius < (double) 0.0)
    throw std::runtime_error("Attempting to set negative radius in ConePrimitive (constructor)");

  _radius = radius;
  _height = height;
  _npoints = 10;
  _nrings = 1;
  calc_mass_properties();
}

/// Constructs a cone along the y-axis and centered at the origin with specified, radius, height, points and rings
ConePrimitive::ConePrimitive(double radius, double height, unsigned npoints, unsigned nrings, const Pose3d& T) : Primitive(T)
{
  if (height < (double) 0.0)
    throw std::runtime_error("Attempting to set negative height in ConePrimitive (constructor)");
  if (radius < (double) 0.0)
    throw std::runtime_error("Attempting to set negative radius in ConePrimitive (constructor)");
  if (npoints < 3)
    throw std::runtime_error("Attempting to set number of circle points < 3 in ConePrimitive (constructor)");
  if (nrings < 1)
    throw std::runtime_error("Attempting to set number of rings < 1 in ConePrimitive (constructor)");

  _radius = radius;
  _height = height;
  _npoints = npoints;
  _nrings = nrings;
  calc_mass_properties();
}

/// Gets a supporting point in a particular direction
Point3d ConePrimitive::get_supporting_point(const Vector3d& d) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(d.pose)) != _poses.end());

  // scale the vector
  Vector3d dscal = Vector3d::normalize(d) * (_radius + _height) * 2.0;

  // setup the zero vector
  Vector3d zero(0.0, 0.0, 0.0, dscal.pose);

  // setup a line segment
  Point3d isect;
  Vector3d normal;
  double t;
  bool intersects = intersect_seg(LineSeg3(dscal, zero), t, isect, normal);
  assert(intersects);

  // return the point of intersection
  return isect;
}

/// Computes the distance and normal from a point on the primitive 
double ConePrimitive::calc_dist_and_normal(const Point3d& p, std::vector<Vector3d>& normals) const
{
  const unsigned X = 0, Y = 1, Z = 2;

/*
  // compute distance from point (projected onto plane) to circle
  double dist = std::sqrt(p[X]*p[X] + p[Z]*p[Z]) - _radius;

  // setup closest point on cone
  pcyl.pose = get_pose();
  pcyl[X] = p[X];
  pcyl[Y] = p[Y];
  pcyl[Z] = p[Z];
*/

  // initialize the closest point on the cone
  Point3d pcon;

  // setup cone pose 
  pcon.pose = p.pose;

  // setup the normal
  normals.push_back(Vector3d());
  Vector3d& normal = normals.back();
  normal.pose = p.pose;

  // setup distance
  double dist;

  // see whether the point is above the height extent
  double ht_ext = _height*0.5;
  if (p[Y] > ht_ext)
  {
    // closest point is the height extent
    pcon[X] = pcon[Z] = 0.0;
    pcon[Y] = ht_ext;

    // setup the normal
    normal[X] = normal[Z] = 0.0;
    normal[Y] = 1.0;

    // return the distance
    return (pcon - p).norm();
  }
  else
  {
    // see whether the point is below the negative height extent
    if (p[Y] < -ht_ext)
    {
      // compute distance from point (projected onto plane) to circle
      double d = std::sqrt(p[X]*p[X] + p[Z]*p[Z]) - _radius;

      // setup the closest point on the cone
      if (d > 0.0)
      {
        pcon[Y] = 0.0;
        pcon[X] = p[X];
        pcon[Z] = p[Z];
        pcon *= _radius;
      }
      pcon[Y] = -ht_ext;

      // setup the normal
      normal[X] = normal[Z] = 0.0;
      normal[Y] = -1.0;

      // compute distance
      return (p - pcon).norm();
    }
    else
    {
      // compute the radius at the given height
      double b =  _radius*0.5;
      double m = -b/ht_ext;
      double rad = m*p[Y] + b;
      double d = std::sqrt(p[X]*p[X] + p[Z]*p[Z]) - rad;

      // setup the closest point on the cone
      if (d > 0.0)
      {
        pcon[Y] = 0.0;
        pcon[X] = p[X];
        pcon[Z] = p[Z];
        pcon *= rad;
      }
      pcon[Y] = p[Y];

      // setup the normal
      normal = Vector3d::normalize(pcon);

      // compute distance
      return (p - pcon).norm();
    }
  }
}

double ConePrimitive::calc_dist(const SpherePrimitive* s, Point3d& pcon, Point3d& psph) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the sphere center and put it into the cone's coordinate system
  Point3d sph_c(0.0, 0.0, 0.0, psph.pose);
  Point3d p = Pose3d::transform_point(pcon.pose, sph_c);

/*
  // compute distance from point (projected onto plane) to circle
  double dist = std::sqrt(p[X]*p[X] + p[Z]*p[Z]) - _radius;

  // setup closest point on cone
  pcyl.pose = get_pose();
  pcyl[X] = p[X];
  pcyl[Y] = p[Y];
  pcyl[Z] = p[Z];
*/

  // setup distance
  double dist;

  // see whether the point is above the height extent
  double ht_ext = _height*0.5;
  if (p[Y] > ht_ext)
  {
    // closest point is the height extent
    pcon[X] = pcon[Z] = 0.0;
    pcon[Y] = ht_ext;
  }
  else
  {
    // see whether the point is below the negative height extent
    if (p[Y] < -ht_ext)
    {
      // compute distance from point (projected onto plane) to circle
      double d = std::sqrt(p[X]*p[X] + p[Z]*p[Z]) - _radius;

      // setup the closest point on the cone
      if (d > 0.0)
      {
        pcon[Y] = 0.0;
        pcon[X] = p[X];
        pcon[Z] = p[Z];
        pcon *= _radius;
      }
      pcon[Y] = -ht_ext;

      // compute distance
      dist = (p - pcon).norm() - s->get_radius();
    }
    else
    {
      // compute the radius at the given height
      double b =  _radius*0.5;
      double m = -b/ht_ext;
      double rad = m*p[Y] + b;
      double d = std::sqrt(p[X]*p[X] + p[Z]*p[Z]) - rad;

      // setup the closest point on the cone
      if (d > 0.0)
      {
        pcon[Y] = 0.0;
        pcon[X] = p[X];
        pcon[Z] = p[Z];
        pcon *= rad;
      }
      pcon[Y] = p[Y];

      // compute distance
      dist = (p - pcon).norm() - s->get_radius();
    }
  }

  // determine closest point on the sphere
  if (dist < 0.0)
  {
    // compute farthest interpenetration of cone inside sphere
    Point3d con_c(0.0, 0.0, 0.0, pcon.pose);
    psph = sph_c - Pose3d::transform_point(psph.pose, con_c);
    psph.normalize();
    psph *= -dist;
  }
  else
  {
    // determine the closest point on the sphere using the vector from the
    // sphere center to the closest point on the cone 
    psph = Pose3d::transform_point(psph.pose, pcon) - sph_c;
    psph.normalize();
    psph *= s->get_radius();
  }

  return dist;
}

/// Sets the radius for this cone
void ConePrimitive::set_radius(double radius)
{
  const unsigned X = 0, Y = 1, Z = 2;

  _radius = radius;
  if (_radius < 0.0)
    throw std::runtime_error("Attempting to pass negative radius to ConePrimitive::set_radius()");

  // need to recalculate mass properties
  calc_mass_properties();

  // need to update visualization
  update_visualization();

  // set lengths on each OBB
  for (map<CollisionGeometryPtr, OBBPtr>::iterator i = _obbs.begin(); i != _obbs.end(); i++)
  {
    // setup OBB half-lengths
    i->second->l[X] = _radius;
    i->second->l[Y] = _height*0.5;
    i->second->l[Z] = _radius;
  }
}

/// Sets the height for this cone
void ConePrimitive::set_height(double height)
{
  const unsigned X = 0, Y = 1, Z = 2;

  _height = height;
  if (_height < 0.0)
    throw std::runtime_error("Attempting to pass negative height to ConePrimitive::set_height()");

  // need to recalculate mass properties
  calc_mass_properties();

  // need to update visualization
  update_visualization();

  // set lengths on each OBB
  for (map<CollisionGeometryPtr, OBBPtr>::iterator i = _obbs.begin(); i != _obbs.end(); i++)
  {
    // setup OBB half-lengths
    i->second->l[X] = _radius;
    i->second->l[Y] = _height*0.5;
    i->second->l[Z] = _radius;
  }
}

/// Sets the number of points in the rings of the cone 
void ConePrimitive::set_circle_points(unsigned n)
{
  _npoints = n;
  if (_npoints < 4)
    throw std::runtime_error("Too few points to represent a circle in ConePrimitive::set_circle_points()");
}

/// Sets the number of rings in the cone 
void ConePrimitive::set_num_rings(unsigned n)
{
  _nrings = n;
  if (_npoints < 1)
    throw std::runtime_error("Too few rings in ConePrimitive::set_num_rings()");
}

/// Transforms the primitive
void ConePrimitive::set_pose(const Pose3d& p)
{
  // convert p to a shared pointer
  shared_ptr<Pose3d> x(new Pose3d(p));

  // determine the transformation from the old pose to the new one 
  Transform3d T = Pose3d::calc_relative_pose(_F, x);

  // go ahead and set the new transform
  Primitive::set_pose(p);

  // recalculate the mass properties
  calc_mass_properties();
}

/// Gets the triangle mesh for the cone, computing it if necessary
shared_ptr<const IndexedTriArray> ConePrimitive::get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P)
{
  // verify that the primitive knows about this pose 
  assert(_poses.find(const_pointer_cast<Pose3d>(P)) != _poses.end());

  const double R = _radius;
  const double HH = _height * (double) 0.5;

  // if radius or height = 0, # of rings < 1 or circle points < 3, return now
  if (_radius <= 0.0 || _height <= 0.0 || _nrings < 1 || _npoints < 3)
    return shared_ptr<IndexedTriArray>(new IndexedTriArray());

  // determine the vertices in the mesh
  // NOTE: we wish for points to be defined in the GLOBAL frame
  list<Point3d> points; 
  for (unsigned i=0; i< _npoints; i++)
  {
    const double THETA = i * (M_PI * (double) 2.0/_npoints);
    const double CT = std::cos(THETA);
    const double ST = std::sin(THETA);
    points.push_back(Point3d(CT*R, -HH, ST*R, P));
  }

  // create one more vertex for the tip of the cone
  points.push_back(Point3d(0.0, HH, 0.0, P));

  // compute the convex hull
  PolyhedronPtr hull = CompGeom::calc_convex_hull(points.begin(), points.end());

  // set the mesh
  const std::vector<Origin3d>& v = hull->get_vertices();
  const std::vector<IndexedTri>& f = hull->get_facets();
  return shared_ptr<IndexedTriArray>(new IndexedTriArray(v.begin(), v.end(), f.begin(), f.end()));
}

/// Creates the visualization for this primitive
osg::Node* ConePrimitive::create_visualization()
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  #ifdef USE_OSG
  // setup the cone
  osg::Cone* cone = new osg::Cone;
  cone->setRadius((float) _radius);
  cone->setHeight((float) _height);
  osg::Geode* geode = new osg::Geode;
  geode->addDrawable(new osg::ShapeDrawable(cone));

  // setup a transform by -90 deg around x
  osg::Matrixd T;
  T(X,X) = 1.0; T(Y,X) = 0.0;  T(Z,X) = 0.0;  T(X,W) = 0.0;
  T(X,Y) = 0.0; T(Y,Y) = 0.0;  T(Z,Y) = 1.0; T(Y,W) = 0.0;
  T(X,Z) = 0.0; T(Y,Z) = -1.0;  T(Z,Z) = 0.0;  T(Z,W) = 0.0;
  T(W,X) = T(W,Z) = 0.0;
  T(W,Y) = -_height*0.25;
  T(W,W) = 1.0;
  
  // setup the matrix transform
  osg::MatrixTransform* G = new osg::MatrixTransform;
  G->setMatrix(T);
  G->addChild(geode);

  return G;
  #else
  return NULL;
  #endif
}

/// Implements Base::load_from_xml() for serialization
void ConePrimitive::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // verify that the node type is cone
  assert(strcasecmp(node->name.c_str(), "Cone") == 0);

  // load the parent data
  Primitive::load_from_xml(node, id_map);

  // read in the radius, if specified
  XMLAttrib* radius_attr = node->get_attrib("radius");
  if (radius_attr)
    _radius = radius_attr->get_real_value();

  // read in the height, if specified
  XMLAttrib* height_attr = node->get_attrib("height");
  if (height_attr)
    _height = height_attr->get_real_value();

  // read in the number of circle points, if specified
  XMLAttrib* npoints_attr = node->get_attrib("num-circle-points");
  if (npoints_attr)
    _npoints = npoints_attr->get_unsigned_value();

  // read in the number of rings of the cylinder, if specified
  XMLAttrib* nrings_attr = node->get_attrib("num-rings");
  if (nrings_attr)
    _nrings = nrings_attr->get_unsigned_value();

  // recompute mass properties
  calc_mass_properties();
}

/// Implements Base::save_to_xml() for serialization
void ConePrimitive::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data
  Primitive::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "Cone";

  // save the radius
  node->attribs.insert(XMLAttrib("radius", _radius));

  // save the height
  node->attribs.insert(XMLAttrib("height", _height));

  // save the number of circle points 
  node->attribs.insert(XMLAttrib("num-circle-points", _npoints));

  // save the number of rings
  node->attribs.insert(XMLAttrib("num-rings", _nrings));
}

/// Calculates mass properties for the cone
void ConePrimitive::calc_mass_properties()
{
  // get the current transform
  shared_ptr<const Pose3d> T = get_pose();

  // determine the radius squared (we'll need this) 
  const double RSQ = _radius * _radius;

  // compute the mass if density is given
  if (_density)
  {
    const double volume = 1.0/3.0 * M_PI * RSQ * _height;
    _J.m = *_density * volume;
  }

  // compute the non-longitudinal elements
  const double HSQ = _height * _height;
  const double NL_ELM = 0.1 * _J.m * HSQ + (3.0/20.0) * _J.m * RSQ;
  const double LONG_ELM = (1.0/3.0) * _J.m * RSQ;

  // compute the inertia matrix
  _J.J = Matrix3d(NL_ELM, 0, 0, 0, LONG_ELM, 0, 0, 0, NL_ELM);
}

/// Finds the signed distance between the cylinder and another primitive
double ConePrimitive::calc_signed_dist(shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const
{
  shared_ptr<const HeightmapPrimitive> hmp = dynamic_pointer_cast<const HeightmapPrimitive>(p);
  if (hmp)
    return hmp->calc_signed_dist(dynamic_pointer_cast<const Primitive>(shared_from_this()), pp, pthis);

  // if the primitive is convex, can use GJK
  if (p->is_convex())
  {
    shared_ptr<const Pose3d> Pbox = pthis.pose;
    shared_ptr<const Pose3d> Pgeneric = pp.pose;
    shared_ptr<const Primitive> bthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return GJK::do_gjk(bthis, p, Pbox, Pgeneric, pthis, pp);
  }

  // try cone/(non-convex) trimesh
  shared_ptr<const TriangleMeshPrimitive> trip = dynamic_pointer_cast<const TriangleMeshPrimitive>(p);
  if (trip)
    return trip->calc_signed_dist(dynamic_pointer_cast<const Primitive>(shared_from_this()), pp, pthis);

  assert(false);
  return 0.0;
}

/// Gets vertices from the primitive
void ConePrimitive::get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& verts) const
{
  // clear the vector of vertices
  verts.clear();

  // setup constant for the expanded radius
  const double H = _height;

  // if the radius, height, or rings = 0 or circle points < 3, return now
  if (_radius == 0.0 || _height == 0.0 || _nrings == 0 || _npoints < 3)
    return; 

  // create vertices
  for (unsigned j=0; j< _nrings; j++)
  {
    const double HEIGHT = -(H * (double) 0.5) + (j*H)/_nrings;
    const double R = _radius * (double) (_nrings - j)/_nrings;
    for (unsigned i=0; i< _npoints; i++)
    {
      const double THETA = i*(M_PI * (double) 2.0/_npoints);
      const double CT = std::cos(THETA);
      const double ST = std::sin(THETA);
      verts.push_back(Point3d(CT*R, HEIGHT, ST*R, P));
    }

    // create one more vertex for the tip of the cone
    verts.push_back(Point3d(0.0, H * (double) 0.5, 0.0, P));
  }
}

/// Gets vertices from the primitive
/*
void ConePrimitive::get_vertices(CollisionGeometryPtr geom, std::vector<const Point3d*>& vertices)
{
  // get the vertices for the geometry
  vector<Point3d>& verts = _vertices[geom];

  // create the vector of vertices if necessary
  if (verts.empty())
  {
    // setup constant for the expanded radius
    const double H = _height;

    // if the radius, height, or rings = 0 or circle points < 3, return now
    if (_radius == 0.0 || _height == 0.0 || _nrings == 0 || _npoints < 3)
    {
      vertices.clear();
      return; 
    }

    // get the pose for the geometry
    shared_ptr<const Pose3d> gpose = geom->get_pose();

    // get the pose for this geometry
    shared_ptr<const Pose3d> P = get_pose(); 
    assert(!P->rpose);

    // setup transform
    Transform3d T;
    T.source = gpose;
    T.target = gpose;
    T.x = P->x;
    T.q = P->q;

    // create vertices
    for (unsigned j=0; j< _nrings; j++)
    {
      const double HEIGHT = -(H * (double) 0.5) + (j*H)/_nrings;
      const double R = _radius * (double) (_nrings - j)/_nrings;
      for (unsigned i=0; i< _npoints; i++)
      {
        const double THETA = i*(M_PI * (double) 2.0/_npoints);
        const double CT = std::cos(THETA);
        const double ST = std::sin(THETA);
        verts.push_back(T.transform_point(Point3d(CT*R, HEIGHT, ST*R, gpose)));
      }
    }

    // create one more vertex for the tip of the cone
    verts.push_back(T.transform_point(Point3d(0.0, H * (double) 0.5, 0.0, gpose)));
  }

  // copy the addresses of the computed vertices into 'vertices' 
  vertices.resize(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    vertices[i] = &verts[i];
}
*/

/// Gets the OBB
BVPtr ConePrimitive::get_BVH_root(CollisionGeometryPtr geom)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the pointer to the bounding box
  OBBPtr& obb = _obbs[geom];

  // create the bounding box, if necessary
  if (!obb)
  {
    // create the bounding box
    obb = shared_ptr<OBB>(new OBB);
    obb->geom = geom;

    // get the pose for this geometry
    shared_ptr<const Pose3d> P = get_pose(geom); 

    // setup the obb center and orientation
    obb->center = Point3d(0.0, 0.0, 0.0, P);
    obb->R.set_identity();

    // setup OBB half-lengths
    obb->l[X] = _radius;
    obb->l[Y] = _height*0.5;
    obb->l[Z] = _radius;
  }

  return obb;
}

/// The signum function
namespace Moby {
static double sgn(double x) { return x / std::fabs(x); }
}

/// Determines whether a point is inside the cone; if so, determines the normal
/**
 * Derived/adapted from Eberly. D.  "Intersection of a Line and a Cone"
 */
bool ConePrimitive::point_inside(const Point3d& p, Vector3d& normal) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end());

  // determine the angle theta
  double theta = std::atan(_radius/_height);

  // setup A
  Vector3d A(0,-1,0, p.pose);

  // setup v
  Vector3d V(0,_height,0, p.pose);

  // see whether the point is outside
  if (A.dot(Vector3d::normalize(p - V)) < std::cos(theta))
    return false;

  // determine the normal (the gradient: derived using Mathematica)
  normal = determine_normal(p);

  return true;
}

/// Computes the signed distance from the cylinder
double ConePrimitive::calc_signed_dist(const Point3d& p) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end());

  // determine the angle theta
  double theta = std::atan(_radius/_height);

  // check for point being above cylinder
  if (p[Y] > _height*0.5)
    return (p[Y] - _height*0.5);
  else if (p[Y] < -_height*0.5)
  {
    double d1 = -_height*0.5 - p[Y];
    double d2 = sqr(p[X]) + sqr(p[Z]) - _radius*_radius;
    if (d2 < 0.0)
      return -d1;
    else
      return -std::sqrt(sqr(d1) + d2);
  }

  // get the radius of the cone at the vertical location of the point
  // radius at +1/2 height = 0
  // radius at -1/2 height = R
  const double RR = -_radius * (p[Y] / _height) + (double) 0.5 * _radius ;

  // get the distance from the horizontal part of the cone
  double dcone = RR - std::sqrt(sqr(p[X]) + sqr(p[Z]));

  // check whether point is outside the cone
  if (dcone > 0.0)
    return std::sqrt(dcone);

  // get the distance from the vertical parts of the cone
  double dv1 = (double) 0.5 * _height - p[Y];
  double dv2 = (double) 0.5 * _height + p[Y];

  return std::min(std::min(dv1, dv2), dcone);
}

/// Computes the penetration depth of a point inside the cone
/**
 * \return the penetration depth, or -INF if the point is outside the cone
 */
double ConePrimitive::calc_penetration_depth(const Point3d& p) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const double INF = std::numeric_limits<double>::max();

  // determine the angle theta
  double theta = std::atan(_radius/_height);

  // setup A
  Vector3d A(0,-1,0,p.pose);

  // setup v
  Vector3d V(0,_height,0,p.pose);

  // see whether the point is outside
  if (A.dot(Vector3d::normalize(p - V)) < std::cos(theta))
    return -INF;

  // get the radius of the cone at the vertical location of the point
  // radius at +1/2 height = 0
  // radius at -1/2 height = R
  const double RR = -_radius * (p[Y] / _height) + (double) 0.5 * _radius ;

  // get the distance from the horizontal part of the cone
  double dcone = RR - std::sqrt(sqr(p[X]) + sqr(p[Z]));

  // get the distance from the vertical parts of the cone
  double dv1 = (double) 0.5 * _height - p[Y];
  double dv2 = (double) 0.5 * _height + p[Y];

  return std::min(std::min(dv1, dv2), dcone);
}

/// Determines the normal to a point on the cone
/**
 * \note the normal may be degenerate (NaN)
 */
Vector3d ConePrimitive::determine_normal(const Point3d& query) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(query.pose)) != _poses.end());
  const unsigned X = 0, Y = 1, Z = 2;

  // determine the normal (the gradient: derived using Mathematica)
  const double xx = query[X];
  const double xy = query[Y];
  const double xz = query[Z];
  const double usqrt = xx*xx + (xy-_height)*(xy-_height) + xz*xz;
  const double D = std::pow(usqrt, (double) 1.5);

  // setup the normal
  Vector3d normal(query.pose);
  normal[X] = (xy - _height)*std::fabs(xx)*Moby::sgn(xx)/D;
  normal[Y] = -1.0/std::sqrt(usqrt) + (xy - _height)*std::fabs(xy-_height)*Moby::sgn(xy-_height)/D;
  normal[Z] = (xy - _height)*std::fabs(xz)*Moby::sgn(xz)/D;

  return normal;
}

/// Intersects a cone with a line segment
/**
 * Adapted from Geometric Tools cone/line segment intersection
 * \note for line segments that are partially or fully inside the cone, the
 *       method only returns intersection if the second endpoint of the segment
 *       is farther inside than the first
 */
bool ConePrimitive::intersect_seg(const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(seg.first.pose)) != _poses.end());
  const unsigned Y = 1;

  // we'll need a second intersection point for temporary storage
  Point3d isects[3];

  // first check whether the first point is inside/on the cone
  double dp = calc_penetration_depth(seg.first);
  if (dp >= (double) 0.0)
  {
    // first point is inside/on the cone
    // get the normal
    if (!point_inside(seg.first, normal))
      return false;

    // indicate intersection
    t = (double) 0.0;
    isect = seg.first;
    return true;
  }

  // alias the segments
  const Point3d& p = seg.first;
  const Point3d& q = seg.second;

  // determine the angle theta
  double theta = std::atan(_radius/_height);

  // determine the unit-length line direction vector
  Vector3d dir = Vector3d::normalize(q - p);
  double fAdD = -dir[Y];
  double fCosSqr = std::cos(theta);
  fCosSqr *= fCosSqr;
  Vector3d kE = p - Vector3d(0,_height,0,p.pose);
  double fAdE = -kE[Y];
  double fDdE = dir.dot(kE);
  double fEdE = kE.dot(kE);
  double fC2 = fAdD*fAdD - fCosSqr;
  double fC1 = fAdD*fAdE - fCosSqr*fDdE;
  double fC0 = fAdE*fAdE - fCosSqr*fEdE;

  // solve the quadratic; keep only those points for which A dot X-V >= 0
  if (std::fabs(fC2) >= NEAR_ZERO)
  {
    // c2 != 0
    double fDiscr = fC1*fC1 - fC0*fC2;
    if (fDiscr < (double) 0.0)
    {
      // Q(t) = 0 has no real-valued roots; the line does not intersect
      // the double sided cone
      return false;
    }
    else if (fDiscr > NEAR_ZERO)
    {
      // Q(t) = 0 has two distinct real-valued roots.  However, one or
      // both o fthem might intersect the portion of the double-sided cone
      // behind the vertex.  We are interested only in those intersections in
      // "front" of the vertex
      double fRoot = std::sqrt(fDiscr);
      double fInvC2 = ((double) 1.0)/fC2;
      unsigned nisect = 0;

      double fT = (-fC1 - fRoot)*fInvC2;
      isects[nisect] = p + dir*fT;
      kE = isects[nisect] - Vector3d(0,_height,0,p.pose);          
      double fDot = -kE[Y];
      if (fDot > (double) 0.0)
        nisect++;

      fT = (-fC1 + fRoot)*fInvC2;
      isects[nisect] = p + dir*fT;
      kE = isects[nisect] - Vector3d(0,_height,0,p.pose);
      fDot = -kE[Y];
      if (fDot > (double) 0.0)
        nisect++;

      if (nisect == 2)
      {
        // line intersects the single sided cone in front of the vertex twice
        // determine which point is closer to the first
        if ((isects[0]-p).norm() > (isects[1]-p).norm())
          std::swap(isects[0],isects[1]);

        // determine the intersection parameter
        t = std::sqrt((isects[0] - p).norm_sq()/dir.norm_sq());

        // transform that point from the cone frame
        isect = isects[0];

        // determine the normal to the point
        normal = determine_normal(isects[0]);

        return true;
      }
      else if (nisect == 1)
      {
        // line intersects the single-sided cone in front of the vertex once.
        // the other intersection is with the single-sided cone behind the
        // vertex

        // determine the intersection parameter
        t = std::sqrt((isects[0] - p).norm_sq()/dir.norm_sq());

        // transform the intersection point from the cone frame
        isect = isects[0];

        // determine the normal to the point
        normal = determine_normal(isects[0]);

        return true;
      }
      else
      {
        // line intersects single-sided cone behind the vertex twice
        return false;
      }
    }
    else
    {
      // one repeated real root (line is tangent to the cone)
      isects[0] = p - dir*(fC1/fC2);
      kE = isects[0] - Vector3d(0,_height,0,p.pose);
      if (-kE[Y] > (double) 0.0)
      {
        // have a single point intersection
        // determine the intersection parameter
        t = std::sqrt((isects[0] - p).norm_sq()/dir.norm_sq());

        // transform the intersection point from the cone frame
        isect = isects[0];

        // determine the normal to the point
        normal = determine_normal(isects[0]);

        return true;
      }
      else
      {
        // no intersection
        return false;
      }
    }
  }
  else if (std::fabs(fC1) > NEAR_ZERO)
  {
    // c2 = 0, c1 != 0 (D is a direction vector on the cone boundary)
    isects[0] = p - (((double) 0.5) * fC0/fC1)*dir;
    kE = isects[0] - Vector3d(0,_height,0,p.pose);
    double fDot = -kE[Y];
    if (fDot > (double) 0.0)
    {
      // ray intersection
      // transform point from cone frame
      isect = isects[0];

      // determine the intersection parameter
      t = std::sqrt((isects[0] - p).norm_sq()/dir.norm_sq());
      
      // determine the normal to the point
      normal = determine_normal(isects[0]);

      return true;
    }
  }
  else if (std::fabs(fC0) > NEAR_ZERO)
  {
    // c2 = c1 = 0, c0 != 0
    return false;
  }
  else
  {
    // c2 = c1 = c0 = 0, cone contains ray V+t*D, where V is cone vertex
    // and D is the line direction; NOTE: normal is degenerate, so we return
    // no intersection
    return false;
  }

  // shouldn't still be here
  return false;
}

