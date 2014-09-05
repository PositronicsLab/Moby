/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifdef USE_OSG
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#endif
#include <Moby/Constants.h>
#include <Moby/CompGeom.h>
#include <Moby/sorted_pair>
#include <Moby/XMLTree.h>
#include <Moby/BoundingSphere.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/PlanePrimitive.h>
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/GJK.h>
#include <Moby/SpherePrimitive.h>

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

/// Creates a sphere with radius 1.0 and 100 points 
SpherePrimitive::SpherePrimitive()
{
  _radius = 1.0;
  _npoints = 100;
  calc_mass_properties();
}

/// Creates a sphere with radius 1.0 and 100 points at the given transform
SpherePrimitive::SpherePrimitive(const Pose3d& T) : Primitive(T)
{
  _radius = 1.0;
  _npoints = 100;
  calc_mass_properties();
}

/// Creates a sphere with the specified radius and 100 points 
SpherePrimitive::SpherePrimitive(double radius)
{
  _radius = radius;
  _npoints = 100;
  calc_mass_properties();
}

/// Creates a sphere with the specified radius and number of points
SpherePrimitive::SpherePrimitive(double radius, unsigned n)
{
  _radius = radius;
  _npoints = n;
  calc_mass_properties();
}

/// Creates a sphere with the specified radius and transform
/**
 * The sphere is composed of 100 points.
 */
SpherePrimitive::SpherePrimitive(double radius, const Pose3d& T) : Primitive(T)
{
  _radius = radius;
  _npoints = 100;
  calc_mass_properties();
}

/// Creates a sphere with the specified radius, transform, and number of points 
SpherePrimitive::SpherePrimitive(double radius, unsigned n, const Pose3d& T) : Primitive(T)
{
  _radius = radius;
  _npoints = n;  
  calc_mass_properties();
}

/// Gets the supporting point
Point3d SpherePrimitive::get_supporting_point(const Vector3d& d) const 
{
  assert(_poses.find(const_pointer_cast<Pose3d>(d.pose)) != _poses.end());

  return Vector3d::normalize(d)*_radius;
}

/// Computes the signed distance of the given point from this primitive
double SpherePrimitive::calc_signed_dist(const Point3d& p) const
{
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end());

  return p.norm() - _radius;
}

/// Computes the distance from another sphere primitive
double SpherePrimitive::calc_signed_dist(shared_ptr<const SpherePrimitive> s, Point3d& pthis, Point3d& ps) const
{
  // get the transform from s to this
  Transform3d T = Pose3d::calc_relative_pose(ps.pose, pthis.pose);

  // compute the distance
  double d = T.x.norm() - _radius - s->_radius;

  // setup sphere centers in alternate frames
  Point3d ps_c(0.0, 0.0, 0.0, ps.pose);
  Point3d pthis_c(0.0, 0.0, 0.0, pthis.pose);

  // setup closest points
  pthis = T.transform_point(ps_c);
  pthis.normalize();
  ps = T.inverse_transform_point(pthis_c);
  ps.normalize();

  // scale closest points appropriately
  if (d > 0.0)
  {
    pthis *= _radius;
    ps *= s->_radius;
  }
  else
  {
    pthis *= _radius+d;
    ps *= s->_radius+d;
  }

  return d;
}

/// Calculates mass properties for the sphere 
void SpherePrimitive::calc_mass_properties() 
{
  // get the current transform
  shared_ptr<const Pose3d> T = get_pose();

  // compute the mass if density is given
  if (_density)
  {
    const double volume = M_PI * _radius * _radius * _radius * 4.0/3.0; 
    _J.m = *_density * volume;
  }

  // compute the diagonal element of the untransformed inertia matrix
  const double diag = _radius * _radius * _J.m * 2.0/5.0;

  // form the inertia matrix (untransformed)
  _J.J = Matrix3d(diag, 0, 0, 0, diag, 0, 0, 0, diag);
}

/// Sets the radius for this sphere (forces redetermination of the mesh)
void SpherePrimitive::set_radius(double radius)
{
  _radius = radius;
  if (_radius < 0.0)
    throw std::runtime_error("Attempting to set negative radius in call to SpherePrimitive::set_radius()");

  // recalculate mass properties
  calc_mass_properties();

  // need to update the visualization
  update_visualization();

  // set radius on each bounding sphere
  for (map<CollisionGeometryPtr, shared_ptr<BoundingSphere> >::iterator i = _bsphs.begin(); i != _bsphs.end(); i++)
    i->second->radius = _radius;
}

/// Sets the number of points used in this sphere 
/**
 * \param n the number of points
 * \note forces redetermination of the mesh
 */
void SpherePrimitive::set_num_points(unsigned n)
{
  _npoints = n;
}

/// Transforms the primitive
void SpherePrimitive::set_pose(const Pose3d& p)
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

/// Gets the mesh, computing it if necessary
shared_ptr<const IndexedTriArray> SpherePrimitive::get_mesh(shared_ptr<const Pose3d> T)
{
  // if the radius is zero, create an empty mesh 
  if (_radius == 0.0)
  {
    return shared_ptr<IndexedTriArray>(new IndexedTriArray());
  }

  // verify that the primitive knows about this pose 
  assert(_poses.find(const_pointer_cast<Pose3d>(T)) != _poses.end());

  // determine the vertices in the mesh
  // NOTE: they will all be defined in the global frame
  list<Point3d> points;
  const double INC = (double) M_PI * ((double) 3.0 - std::sqrt((double) 5.0));
  const double OFF = (double) 2.0 / _npoints;
  for (unsigned k=0; k< _npoints; k++)
  {
    const double Y = k * OFF - (double) 1.0 + (OFF * (double) 0.5);
    const double R = std::sqrt((double) 1.0 - Y*Y);
    const double PHI = k * INC;
    Vector3d unit(std::cos(PHI)*R, Y, std::sin(PHI)*R);
    points.push_back(T->transform_point(unit*_radius));
  }

  // add points at the extents
  points.push_back(T->transform_point(Vector3d(+1*_radius,0,0)));
  points.push_back(T->transform_point(Vector3d(-1*_radius,0,0)));
  points.push_back(T->transform_point(Vector3d(0,+1*_radius,0)));
  points.push_back(T->transform_point(Vector3d(0,-1*_radius,0)));
  points.push_back(T->transform_point(Vector3d(0,0,+1*_radius)));
  points.push_back(T->transform_point(Vector3d(0,0,-1*_radius)));

  // compute the convex hull
  PolyhedronPtr hull = CompGeom::calc_convex_hull(points.begin(), points.end());

  // set the mesh
  const vector<Origin3d>& v = hull->get_vertices();
  const vector<IndexedTri>& f = hull->get_facets();
  return boost::shared_ptr<IndexedTriArray>(new IndexedTriArray(v.begin(), v.end(), f.begin(), f.end()));
}

/// Gets vertices for the primitive
void SpherePrimitive::get_vertices(shared_ptr<const Pose3d> P, std::vector<Point3d>& vertices) const
{
  // verify that the primitive knows about this mesh
  assert(_poses.find(const_pointer_cast<Pose3d>(P)) != _poses.end());

  // clear the vertices at first 
  vertices.clear(); 

  // look for no vertices 
  if (_radius == 0.0)
    return;

  // resize number of points
  vertices.resize(_npoints+6);

  // determine the vertices in the mesh
  // NOTE: they will all be defined in the global frame
  const double INC = (double) M_PI * ((double) 3.0 - std::sqrt((double) 5.0));
  const double OFF = (double) 2.0 / _npoints;
  for (unsigned k=0; k< _npoints; k++)
  {
    const double Y = k * OFF - (double) 1.0 + (OFF * (double) 0.5);
    const double R = std::sqrt((double) 1.0 - Y*Y);
    const double PHI = k * INC;
    vertices[k] = Vector3d(std::cos(PHI)*R, Y, std::sin(PHI)*R, P)*_radius;
  }

  // setup vertices at extents of each axis
  vertices[_npoints+0] = Vector3d(1.0*_radius, 0.0, 0.0, P);
  vertices[_npoints+1] = Vector3d(-1.0*_radius, 0.0, 0.0, P);
  vertices[_npoints+2] = Vector3d(0.0, 1.0*_radius, 0.0, P);
  vertices[_npoints+3] = Vector3d(0.0, -1.0*_radius, 0.0, P);
  vertices[_npoints+4] = Vector3d(0.0, 0.0, 1.0*_radius, P);
  vertices[_npoints+5] = Vector3d(0.0, 0.0, -1.0*_radius, P);
}

/// Finds the signed distance between the sphere and another primitive
double SpherePrimitive::calc_signed_dist(shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const
{
  // first try box/sphere
  shared_ptr<const BoxPrimitive> boxp = dynamic_pointer_cast<const BoxPrimitive>(p);
  if (boxp)
  {
    shared_ptr<const SpherePrimitive> thisp = dynamic_pointer_cast<const SpherePrimitive>(shared_from_this());
    return boxp->calc_signed_dist(thisp, pp, pthis);
  }

  // now try sphere/sphere
  shared_ptr<const SpherePrimitive> spherep = dynamic_pointer_cast<const SpherePrimitive>(p);
  if (spherep)
    return calc_signed_dist(spherep, pthis, pp);

  // now try plane/sphere
  shared_ptr<const PlanePrimitive> planep = dynamic_pointer_cast<const PlanePrimitive>(p);
  if (planep)
  {
    shared_ptr<const SpherePrimitive> thisp = dynamic_pointer_cast<const SpherePrimitive>(shared_from_this());
    return planep->calc_signed_dist(thisp, pp, pthis);
  }

  // now try heightmap/sphere
  shared_ptr<const HeightmapPrimitive> hmp = dynamic_pointer_cast<const HeightmapPrimitive>(p);
  if (hmp)
  {
    shared_ptr<const SpherePrimitive> thisp = dynamic_pointer_cast<const SpherePrimitive>(shared_from_this());
    return hmp->calc_signed_dist(thisp, pp, pthis);
  }

  // if the primitive is convex, can use GJK
  if (p->is_convex())
  {
    shared_ptr<const Pose3d> Psph = pthis.pose;
    shared_ptr<const Pose3d> Pgeneric = pp.pose;
    shared_ptr<const Primitive> bthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return GJK::do_gjk(bthis, p, Psph, Pgeneric, pthis, pp);
  }

  // try sphere/(non-convex) trimesh
  shared_ptr<const TriangleMeshPrimitive> trip = dynamic_pointer_cast<const TriangleMeshPrimitive>(p);
  if (trip)
    return trip->calc_signed_dist(dynamic_pointer_cast<const Primitive>(shared_from_this()), pp, pthis);

  assert(false);
  return 0.0;
}

/// Finds the signed distance betwen the sphere and a point
double SpherePrimitive::calc_dist_and_normal(const Point3d& p, std::vector<Vector3d>& normals) const
{
  // setup the normal
  normals.push_back(p);
  Vector3d& normal = normals.back();
  double pnorm = p.norm();
  normal /= pnorm;

  // compute the distance
  return pnorm - _radius;
}

/// Creates the visualization for this primitive
osg::Node* SpherePrimitive::create_visualization()
{
  #ifdef USE_OSG
  osg::Sphere* sphere = new osg::Sphere;
  sphere->setRadius((float) _radius);
  osg::Geode* geode = new osg::Geode;
  geode->addDrawable(new osg::ShapeDrawable(sphere));
  return geode;
  #else
  return NULL;
  #endif
}  

/// Implements Base::load_from_xml() for serialization
void SpherePrimitive::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // verify that the node type is sphere
  assert(strcasecmp(node->name.c_str(), "Sphere") == 0);

  // load the parent data
  Primitive::load_from_xml(node, id_map);

  // read in the radius, if specified
  XMLAttrib* radius_attr = node->get_attrib("radius");
  if (radius_attr)
    set_radius(radius_attr->get_real_value());

  // read in the number of points, if specified
  XMLAttrib* npoints_attr = node->get_attrib("num-points");
  if (npoints_attr)
    set_num_points(npoints_attr->get_unsigned_value());

  // recompute mass properties
  calc_mass_properties();
}

/// Implements Base::save_to_xml() for serialization
void SpherePrimitive::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data
  Primitive::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "Sphere";

  // save the radius
  node->attribs.insert(XMLAttrib("radius", _radius));

  // save the number of points 
  node->attribs.insert(XMLAttrib("num-points", _npoints));
}

/// Gets the root bounding volume
BVPtr SpherePrimitive::get_BVH_root(CollisionGeometryPtr geom) 
{
  // get the pointer to the bounding sphere
  shared_ptr<BoundingSphere>& bsph = _bsphs[geom];

  // if the bounding sphere hasn't been created, create it and initialize it
  if (!bsph)
  {
    // create the sphere
    bsph = shared_ptr<BoundingSphere>(new BoundingSphere);
    bsph->geom = geom;

    // get the pose for the geometry
    shared_ptr<const Pose3d> P = get_pose(geom);

    // setup the bounding sphere center; we're assuming that the primitive
    // pose is defined relative to the geometry frame
    bsph->center = Point3d(0.0, 0.0, 0.0, P);
    bsph->radius = _radius;
  }

  return bsph;
}


