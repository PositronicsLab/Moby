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
#include <Moby/Constants.h>
#include <Moby/CompGeom.h>
#include <Moby/sorted_pair>
#include <Moby/XMLTree.h>
#include <Moby/BoundingSphere.h>
#include <Moby/SpherePrimitive.h>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr; 
using std::vector;
using std::list;
using std::pair;
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

  // mesh, vertices are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices = shared_ptr<vector<Point3d> >();
  _smesh = pair<shared_ptr<IndexedTriArray>, list<unsigned> >();
  _invalidated = true;

  // recalculate mass properties
  calc_mass_properties();

  // need to update the visualization
  update_visualization();
}

/// Sets the number of points used in this sphere 
/**
 * \param n the number of points; must be greater than or equal to six
 * \note forces redetermination of the mesh
 */
void SpherePrimitive::set_num_points(unsigned n)
{
  _npoints = n;
  if (n < 5)
    throw std::runtime_error("Attempting to call SpherePrimitive::set_num_points() with n < 5");

  // vertices are no longer valid
  _vertices = shared_ptr<vector<Point3d> >();
  _invalidated = true;
}

/// Sets the intersection tolerance
void SpherePrimitive::set_intersection_tolerance(double tol)
{
  Primitive::set_intersection_tolerance(tol);

  // vertices are no longer valid
  _vertices = shared_ptr<vector<Point3d> >();
  _invalidated = true;
}

/// Transforms the primitive
void SpherePrimitive::set_pose(const Pose3d& p)
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

  // transform mesh
  if (_mesh)
  {
    _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(_mesh->transform(T)));
    _smesh.first = _mesh;
  }

  // transform vertices
  if (_vertices)
    for (unsigned i=0; i< _vertices->size(); i++)
      (*_vertices)[i] = T.transform_point((*_vertices)[i]);

  // invalidate this primitive
  _invalidated = true;

  // recalculate the mass properties
  calc_mass_properties();
}

/// Gets the mesh, computing it if necessary
shared_ptr<const IndexedTriArray> SpherePrimitive::get_mesh()
{
  if (!_mesh)
  {
    // if the radius is zero or the number of points is less than six, create an
    // empty mesh 
    if (_radius == 0.0 || _npoints < 6)
    {
      _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray());
      _smesh = make_pair(_mesh, list<unsigned>());
      return _mesh;
    }

    // get the translation for the transform
    shared_ptr<const Pose3d> T = get_pose();
    const Origin3d& xlat = T->x;

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
      points.push_back(xlat + unit*_radius);
    }

    // compute the convex hull
    PolyhedronPtr hull = CompGeom::calc_convex_hull(points.begin(), points.end());

    // set the mesh
    const vector<Point3d>& v = hull->get_vertices();
    const vector<IndexedTri>& f = hull->get_facets();
    _mesh = boost::shared_ptr<IndexedTriArray>(new IndexedTriArray(v.begin(), v.end(), f.begin(), f.end()));

    // setup sub mesh (it will be just the standard mesh)
    list<unsigned> all_tris;
    for (unsigned i=0; i< _mesh->num_tris(); i++)
      all_tris.push_back(i);
    _smesh = make_pair(_mesh, all_tris);
  }

  return _mesh;
}

/// Gets a sub-mesh for the primitive
const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& SpherePrimitive::get_sub_mesh(BVPtr bv)
{
  if (!_smesh.first)
    get_mesh(); 
  return _smesh;
}

/// Gets vertices for the primitive
void SpherePrimitive::get_vertices(BVPtr bv, std::vector<const Point3d*>& vertices)
{
  // create the vector of vertices if necessary
  if (!_vertices)
  {
    if (_radius == 0.0 || _npoints < 6)
    {
      vertices.clear(); 
      return;
    }

    // get the translation for the transform
    shared_ptr<const Pose3d> T = get_pose();
    const Origin3d& xlat = T->x;

    // determine the vertices in the mesh
    // NOTE: they will all be defined in the global frame
    _vertices = shared_ptr<vector<Point3d> >(new vector<Point3d>(_npoints));
    const double INC = (double) M_PI * ((double) 3.0 - std::sqrt((double) 5.0));
    const double OFF = (double) 2.0 / _npoints;
    for (unsigned k=0; k< _npoints; k++)
    {
      const double Y = k * OFF - (double) 1.0 + (OFF * (double) 0.5);
      const double R = std::sqrt((double) 1.0 - Y*Y);
      const double PHI = k * INC;
      Vector3d unit(std::cos(PHI)*R, Y, std::sin(PHI)*R);
      (*_vertices)[k] = xlat + unit*(_radius + _intersection_tolerance);
    }
  }

  // copy the addresses of the computed vertices into 'vertices' 
  vertices.resize(_vertices->size());
  for (unsigned i=0; i< _vertices->size(); i++)
    vertices[i] = &(*_vertices)[i];
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
BVPtr SpherePrimitive::get_BVH_root() 
{
  // sphere not applicable for deformable bodies 
  if (is_deformable())
    throw std::runtime_error("SpherePrimitive::get_BVH_root() - primitive unusable for deformable bodies!");

  // set the radius and center
  if (!_bsph)
    _bsph = shared_ptr<BoundingSphere>(new BoundingSphere);

  // set the radius and center
  _bsph->radius = _radius + _intersection_tolerance;
  _bsph->center = get_pose()->x;
  _bsph->center.pose = get_pose();

  return _bsph;
}

/// Determines whether there is an intersection between the primitive and a line segment
/**
 * Determines the normal to the primitive if there is an intersection.
 * \note for line segments that are partially or fully inside the sphere, the
 *       method only returns intersection if the second endpoint of the segment
 *       is farther inside than the first
 */
bool SpherePrimitive::intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // account for sphere center in translation
  const Origin3d& center = get_pose()->x;
  Vector3d p = seg.first - center;
  Vector3d q = seg.second - center;

  // get the radius plus the intersection tolerance
  const double R = _radius;

  // determine whether p is already within the sphere
  double pp = p.dot(p);
  if (pp <= R*R)
  {
    // set the intersection
    t = (double) 0.0;
    isect = p + center;
    double pnorm = std::sqrt(pp);
    normal.pose = get_pose();
    if (pnorm > NEAR_ZERO)
      normal = p/pnorm;
    else
      normal.set_zero();
    return true; 
  }

  // look for:
  // (seg.first*t + seg.second*(1-t))^2 = R^2

  // use quadratic formula
  const double px = p[X];
  const double py = p[Y];
  const double pz = p[Z];
  const double qx = q[X];
  const double qy = q[Y];
  const double qz = q[Z];
  const double a = px*px + py*py + pz*pz - 2*px*qx + qx*qx - 2*py*qy + qy*qy -
                 2*pz*qz + qz*qz;
  const double b = 2*px*qx - 2*qx*qx + 2*py*qy - 2*qy*qy + 2*pz*qz - 2*qz*qz;
  const double c = qx*qx + qy*qy + qz*qz - R*R;

  // check for no solution
  if (a == 0.0)
    return false;
  double disc = b*b - 4*a*c;
  if (disc < 0.0)
    return false;

  // compute solutions
  disc = std::sqrt(disc);
  double t1 = (-b + disc)/(2*a);
  double t2 = (-b - disc)/(2*a);

  // look for lowest solution in [0, 1]
  if (t1 < 0.0)
    t1 = 2.0;
  if (t2 < 0.0)
    t2 = 2.0;
  if (t2 < t1)
    std::swap(t1, t2);
  if (t1 < 0.0 || t1 > 1.0)
    return false;

  // compute the point of intersection and normal
  isect = seg.first*t1 + seg.second*(1-t1);
  normal = Vector3d::normalize(Vector3d(isect));

  t = t1;
  return true;
}

/// Determines whether a point is inside or on the sphere
bool SpherePrimitive::point_inside(BVPtr bv, const Point3d& p, Vector3d& normal) const
{
  // get the sphere translation
  const Origin3d& center = get_pose()->x;

  // subtract sphere translation from the query point
  Vector3d query = p - center;

  // check whether query outside of radius
  if (query.norm_sq() > _radius * _radius) 
    return false;

  // determine normal
  normal = Vector3d::normalize(query);

  return true;
}

