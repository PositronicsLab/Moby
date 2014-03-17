/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iostream>
#ifdef USE_OSG
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#endif
#include <Moby/Constants.h>
#include <Moby/CompGeom.h>
#include <Moby/XMLTree.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/OBB.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CylinderPrimitive.h>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using std::map;
using std::list;
using std::vector;
using std::pair;
using std::endl;

/// Constructs a cylinder centered at the origin, with the longitudinal axis aligned with the y-axis, radius 1.0, height 1.0, 10 circle points, and 2 rings
CylinderPrimitive::CylinderPrimitive()
{
  _radius = 1.0;
  _height = 1.0;
  _npoints = 10;
  _nrings = 2;
  calc_mass_properties();
}

/// Constructs a cylinder along the y-axis with specified radius and height, centered at the origin, with 10 circle points and 2 rings
CylinderPrimitive::CylinderPrimitive(double radius, double height)
{
  if (height < (double) 0.0)
    throw std::runtime_error("Attempting to set negative height in CylinderPrimitive (constructor)");
  if (radius < (double) 0.0)
    throw std::runtime_error("Attempting to set negative radius in CylinderPrimitive (constructor)");

  _radius = radius;
  _height = height;
  _npoints = 10;
  _nrings = 2;
  calc_mass_properties();
}

/// Constructs a cylinder along the y-axis with specified radius and height, centered at the origin, with 10 circle points and 2 rings
CylinderPrimitive::CylinderPrimitive(double radius, double height, const Pose3d& T) : Primitive(T)
{
  if (height < (double) 0.0)
    throw std::runtime_error("Attempting to set negative height in CylinderPrimitive (constructor)");
  if (radius < (double) 0.0)
    throw std::runtime_error("Attempting to set negative radius in CylinderPrimitive (constructor)");

  _radius = radius;
  _height = height;
  _npoints = 10;
  _nrings = 2;
  calc_mass_properties();
}

/// Constructs a cylinder along the y-axis and centered at the origin with specified, radius, height, number of points and number of rings
CylinderPrimitive::CylinderPrimitive(double radius, double height, unsigned n, unsigned nrings, const Pose3d& T) : Primitive(T)
{
  if (height < (double) 0.0)
    throw std::runtime_error("Attempting to set negative height in CylinderPrimitive (constructor)");
  if (radius < (double) 0.0)
    throw std::runtime_error("Attempting to set negative radius in CylinderPrimitive (constructor)");
  if (n < 3)
    throw std::runtime_error("Attempting to set number of circle points < 3 in CylinderPrimitive (constructor)");
  if (nrings < 2)
    throw std::runtime_error("Attempting to set number of rings < 2 in CylinderPrimitive (constructor)");

  _radius = radius;
  _height = height;
  _npoints = n;
  _nrings = nrings;
  calc_mass_properties();
}

/// Finds the signed distance between the sphere and another primitive
double CylinderPrimitive::calc_signed_dist(shared_ptr<const Primitive> primitive, shared_ptr<const Pose3d> pose_this, shared_ptr<const Pose3d> pose_p, Point3d& pthis, Point3d& pprimitive) const
{
  // TODO: implement this
  assert(false);
  return 0.0; 
}

/// Gets the supporting point in a particular direction
Point3d CylinderPrimitive::get_supporting_point(const Vector3d& d) 
{
  return Primitive::get_supporting_point(d);
}

/// Computes the signed distance between the cylinder and a point
double CylinderPrimitive::calc_dist_and_normal(const Point3d& p, Vector3d& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // TODO: implement this properly

  // compute distance from point (projected onto plane) to circle
  double dist = std::sqrt(p[X]*p[X] + p[Z]*p[Z]) - _radius;

  // see whether the point is above/below the height extent
  double ht_ext = _height*0.5;
  if (p[Y] > ht_ext)
  {
    // compute the distance from the plane
    double pdist = p[Y] - ht_ext;

    // setup the closest point on the cylinder
    if (dist < 0.0)
      normal = Vector3d(0,1,0,get_pose());
/*
else
    {
      pcyl[Y] = 0.0;
      pcyl.normalize();  
      pcyl *= _radius;
    }
    pcyl[Y] = ht_ext;
*/

    // compute the distance from the cylinder 
    dist = std::min(pdist, dist);
  }
  else if (p[Y] < -ht_ext)
  {
    // compute the distance from the plane
    double pdist = p[Y] + ht_ext;

    // setup the closest point on the cylinder
    if (dist < 0.0)
      normal = Vector3d(0,-1,0,get_pose());
/*
else
    {
      pcyl[Y] = 0.0;
      pcyl.normalize();  
      pcyl *= _radius;
    }
    pcyl[Y] = -ht_ext;
*/

    // compute the distance from the cylinder 
    dist = std::min(pdist, dist);
  }
  else  // point is inside height extents
  {
/*
    // setup the closest point on the cylinder
    if (dist > 0.0)
    {
      pcyl[Y] = 0.0;
      pcyl.normalize();
      pcyl *= _radius; 
    }
    pcyl[Y] = p[Y];
  }
*/
  }

  return dist;
}

/// Gets the distance of this cylinder from a sphere
double CylinderPrimitive::calc_dist(const SpherePrimitive* s, Point3d& pcyl, Point3d& psph) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the sphere center and put it into the cylinder's coordinate system
  Point3d sph_c(0.0, 0.0, 0.0, s->get_pose());
  Point3d p = Pose3d::transform_point(get_pose(), sph_c);

  // compute distance from point (projected onto plane) to circle
  double dist = std::sqrt(p[X]*p[X] + p[Z]*p[Z]) - _radius;

  // setup closest point on cylinder
  pcyl.pose = get_pose();
  pcyl[X] = p[X];
  pcyl[Y] = p[Y];
  pcyl[Z] = p[Z];

  // see whether the point is above/below the height extent
  double ht_ext = _height*0.5;
  if (p[Y] > ht_ext)
  {
    // compute the distance from the plane
    double pdist = p[Y] - ht_ext;

    // setup the closest point on the cylinder
    if (dist > 0.0)
    {
      pcyl[Y] = 0.0;
      pcyl.normalize();  
      pcyl *= _radius;
    }
    pcyl[Y] = ht_ext;

    // compute the distance from the cylinder 
    dist = std::min(pdist, dist);
  }
  else if (p[Y] < -ht_ext)
  {
    // compute the distance from the plane
    double pdist = p[Y] + ht_ext;

    // setup the closest point on the cylinder
    if (dist > 0.0)
    {
      pcyl[Y] = 0.0;
      pcyl.normalize();  
      pcyl *= _radius;
    }
    pcyl[Y] = -ht_ext;

    // compute the distance from the cylinder 
    dist = std::min(pdist, dist);
  }
  else  // point is inside height extents
  {
    // setup the closest point on the cylinder
    if (dist > 0.0)
    {
      pcyl[Y] = 0.0;
      pcyl.normalize();
      pcyl *= _radius; 
    }
    pcyl[Y] = p[Y];
  }

  // now compute distance from sphere
  dist -= s->get_radius();

  // determine closest point on the sphere
  if (dist < 0.0)
  {
    // compute farthest interpenetration of cylinder inside sphere
    Point3d cyl_c(0.0, 0.0, 0.0, get_pose());
    psph = sph_c - Pose3d::transform_point(s->get_pose(), cyl_c);
    psph.normalize();
    psph *= -dist;
  }
  else
  {
    // determine the closest point on the sphere using the vector from the
    // sphere center to the closest point on the cylinder
    psph = Pose3d::transform_point(s->get_pose(), pcyl) - sph_c;
    psph.normalize();
    psph *= s->get_radius();
  }

  return dist;
}

/// Sets the radius for this cylinder
void CylinderPrimitive::set_radius(double radius)
{
  const unsigned X = 0, Y = 1, Z = 2;

  _radius = radius;
  if (_radius < 0.0)
    throw std::runtime_error("Attempting to set negative radius on call to CylinderPrimitive::set_radius()");

  // mesh, vertices are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices.clear();
  _smesh = pair<shared_ptr<IndexedTriArray>, list<unsigned> >();
  _invalidated = true;

  // need to recompute the mass properties of the cylinder
  calc_mass_properties();

  // need up update visualization
  update_visualization();

  // re-set lengths on each OBB
  for (map<CollisionGeometryPtr, OBBPtr>::iterator i = _obbs.begin(); i != _obbs.end(); i++)
  {
    i->second->l[X] = _radius;
    i->second->l[Y] = _height*0.5;
    i->second->l[Z] = _radius;
  }
}

/// Sets the height for this cylinder
void CylinderPrimitive::set_height(double height)
{
  const unsigned X = 0, Y = 1, Z = 2;

  _height = height;
  if (_height < 0.0)
    throw std::runtime_error("Attempting to set negative height on call to CylinderPrimitive::set_height()");

  // mesh, vertices are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices.clear();
  _smesh = pair<shared_ptr<IndexedTriArray>, list<unsigned> >();
  _invalidated = true;
  
  // need to recompute the mass properties of the cylinder
  calc_mass_properties();

  // need up update visualization
  update_visualization();

  // re-set lengths on each OBB
  for (map<CollisionGeometryPtr, OBBPtr>::iterator i = _obbs.begin(); i != _obbs.end(); i++)
  {
    i->second->l[X] = _radius;
    i->second->l[Y] = _height*0.5;
    i->second->l[Z] = _radius;
  }
}

/// Sets the number of points in the circles of the cylinder 
void CylinderPrimitive::set_num_circle_points(unsigned n)
{
  _npoints = n;
  if (n < 3)
    throw std::runtime_error("Attempting to call CylinderPrimitive::set_circle_points() with n < 3");

  // vertices are no longer valid
  _vertices.clear();
  _invalidated = true;
}

/// Sets the number of rings for determining "virtual points" (points on the tube of the cylinder)
void CylinderPrimitive::set_num_rings(unsigned n)
{
  _nrings = n;
  if (_nrings < 2)
    throw std::runtime_error("Attempting to call CylinderPrimitive::set_num_rings() with n < 2");

  // vertices are no longer valid
  _vertices.clear();
  _invalidated = true;
}

/// Gets the triangle mesh for the cylinder, computing it if necessary 
shared_ptr<const IndexedTriArray> CylinderPrimitive::get_mesh()
{
  // compute the mesh if necessary
  if (!_mesh)
  {
    const double R = _radius;
    const double HH = _height;

    // if radius or height = 0, or circle points < 3, return now
    if (_radius <= 0.0 || _height <= 0.0 || _npoints < 3)
    {
      _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray());
      _smesh = make_pair(_mesh, list<unsigned>());
      return _mesh;
    }

    // setup transform so that points are in global frame
    shared_ptr<const Pose3d> P = get_pose();
    Transform3d T = Pose3d::calc_relative_pose(P, GLOBAL); 

    // create vertices evenly spaced in 2D
    std::list<Point3d> points;
    for (unsigned i=0; i< _npoints; i++)
    {
      const double THETA = i * (M_PI * (double) 2.0/_npoints);
      const double CT = std::cos(THETA);
      const double ST = std::sin(THETA);
      points.push_back(T.transform_point(Point3d(CT*R, HH, ST*R, P)));
      points.push_back(T.transform_point(Point3d(CT*R, -HH, ST*R, P)));
    }

    // compute the convex hull
    PolyhedronPtr hull = CompGeom::calc_convex_hull(points.begin(), points.end());

    // get the mesh
    const std::vector<Origin3d>& v = hull->get_vertices();
    const std::vector<IndexedTri>& f = hull->get_facets();
    _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(v.begin(), v.end(), f.begin(), f.end()));

    // setup sub mesh (it will be just the standard mesh)
    list<unsigned> all_tris;
    for (unsigned i=0; i< _mesh->num_tris(); i++)
      all_tris.push_back(i);
    _smesh = make_pair(_mesh, all_tris);
  }

  return _mesh;
}

/// Creates the visualization for this primitive
osg::Node* CylinderPrimitive::create_visualization()
{
  #ifdef USE_OSG
  osg::Cylinder* cyl = new osg::Cylinder;
  cyl->setRadius((float) _radius);
  cyl->setHeight((float) _height);
//  cyl->setRotation(osg::Quat(0, 0, 1, 0));
  cyl->setRotation(osg::Quat(M_PI/2.0f, osg::Vec3f(1.0f, 0.0f, 0.0f)));
  osg::Geode* geode = new osg::Geode;
  geode->addDrawable(new osg::ShapeDrawable(cyl));
  return geode;
  #else
  return NULL;
  #endif
}

/// Implements Base::load_from_xml() for serialization
void CylinderPrimitive::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // verify that the node type is cylinder
  assert(strcasecmp(node->name.c_str(), "Cylinder") == 0);

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
void CylinderPrimitive::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data
  Primitive::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "Cylinder";

  // save the radius
  node->attribs.insert(XMLAttrib("radius", _radius));

  // save the height
  node->attribs.insert(XMLAttrib("height", _height));

  // save the number of circle points 
  node->attribs.insert(XMLAttrib("num-circle-points", _npoints));

  // save the number of rings
  node->attribs.insert(XMLAttrib("num-rings", _nrings));
}

/// Transforms the primitive
void CylinderPrimitive::set_pose(const Pose3d& p)
{
  // convert p to a shared pointer
  shared_ptr<Pose3d> x(new Pose3d(p));

  // determine the transformation from the old pose to the new one 
  Transform3d T = Pose3d::calc_relative_pose(_F, x);

  // go ahead and set the new transform
  Primitive::set_pose(p);

  // invalidate the vertices and mesh
  _mesh.reset();
  _smesh.first.reset();
  _smesh.second.clear();
  _vertices.clear();

  // invalidate this primitive
  _invalidated = true;

  // recalculate the mass properties
  calc_mass_properties();

 // transform bounding volumes
  for (map<CollisionGeometryPtr, OBBPtr>::iterator i = _obbs.begin(); i != _obbs.end(); i++)
  {
    // get the pose for the geometry
    shared_ptr<const Pose3d> gpose = i->first->get_pose();

    // get the pose for this geometry
    shared_ptr<const Pose3d> P = get_pose(); 
    assert(!P->rpose);

    // setup the obb center and orientation
    i->second->center = Point3d(P->x, gpose);
    i->second->R = P->q;
  }
}

/// Calculates mass properties for the cylinder
void CylinderPrimitive::calc_mass_properties()
{
  // get the current transform
  shared_ptr<const Pose3d> T = get_pose();

  // compute the radius squared (we'll need this)
  const double RSQ = _radius * _radius;

  // compute the mass if density is given
  if (_density)
  {
    const double volume = M_PI * RSQ * _height;
    _J.m = *_density * volume;
  }

  // compute the non-longitudinal elements
  const double HSQ = _height * _height;
  const double ONE_TWELFTH = (double) 1.0/12.0;
  const double NL_ELM = ONE_TWELFTH * _J.m * (HSQ + (double) 3.0 * RSQ);
  const double LONG_ELM = 0.5 * _J.m * RSQ;

  // compute the inertia matrix
  _J.J = Matrix3d(NL_ELM, 0, 0, 0, LONG_ELM, 0, 0, 0, NL_ELM);
}

/// Gets the OBB
BVPtr CylinderPrimitive::get_BVH_root(CollisionGeometryPtr geom)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // cylinder not applicable for deformable bodies 
  if (is_deformable())
    throw std::runtime_error("CylinderPrimitive::get_BVH_root(CollisionGeometryPtr geom) - primitive unusable for deformable bodies!");

  // get the pointer to the bounding box
  OBBPtr& obb = _obbs[geom];

  // create the bounding box, if necessary
  if (!obb)
  {
    // create the bounding box
    obb = shared_ptr<OBB>(new OBB);
    obb->geom = geom;

    // get the pose for the geometry
    shared_ptr<const Pose3d> gpose = geom->get_pose();

    // get the pose for this geometry
    shared_ptr<const Pose3d> P = get_pose(); 
    assert(!P->rpose);

    // setup the obb center and orientation
    obb->center = Point3d(P->x, gpose);
    obb->R = P->q;

    // setup OBB half-lengths
    obb->l[X] = _radius;
    obb->l[Y] = _height*0.5;
    obb->l[Z] = _radius;
  }

  return obb;
}

/// Gets a sub-mesh for the primitive
const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& CylinderPrimitive::get_sub_mesh(BVPtr bv)
{
  if (!_smesh.first)
    get_mesh(); 
  return _smesh;
}

/// Gets vertices from the primitive
void CylinderPrimitive::get_vertices(std::vector<Point3d>& verts)
{
  const double R = _radius;
  const double H = _height;

  // clear the vertices
  verts.clear();

  // if radius or height <= 0, num rings < 2, or circle points < 3, return now
  if (_radius <= 0.0 || _height <= 0.0 || _nrings < 2 || _npoints < 3)
    return; 

  // create vertices evenly spaced in 2D
  // NOTE: we wish to create the vertices in the global frame
  for (unsigned j=0; j< _nrings; j++)
  {
    const double HEIGHT = -(H * (double) 0.5) + (j*H)/(_nrings-1); 
    for (unsigned i=0; i< _npoints; i++)
    {
      double THETA = i*(M_PI * (double) 2.0/_npoints);
      const double CT = std::cos(THETA);
      const double ST = std::sin(THETA);
      verts.push_back(Point3d(CT*R, HEIGHT, ST*R, get_pose()));
    }
  }
}

/// Gets vertices from the primitive
/*
void CylinderPrimitive::get_vertices(CollisionGeometryPtr geom, std::vector<const Point3d*>& vertices)
{
  // get the vertices for the geometry
  vector<Point3d>& verts = _vertices[geom];

  // create the vector of vertices if necessary
  if (verts.empty())
  {
    const double R = _radius;
    const double H = _height;

    // if radius or height <= 0, num rings < 2, or circle points < 3, return now
    if (_radius <= 0.0 || _height <= 0.0 || _nrings < 2 || _npoints < 3)
    {
      vertices.clear();
      return; 
    }

    // get the pose for the geometry
    shared_ptr<const Pose3d> gpose = geom->get_pose();

    // get the pose of this primitive
    shared_ptr<const Pose3d> P = get_pose();

    // setup transform
    Transform3d T;
    T.source = gpose;
    T.target = gpose;
    T.x = P->x;
    T.q = P->q;

    // create vertices evenly spaced in 2D
    // NOTE: we wish to create the vertices in the global frame
    for (unsigned j=0; j< _nrings; j++)
    {
      const double HEIGHT = -(H * (double) 0.5) + (j*H)/(_nrings-1); 
      for (unsigned i=0; i< _npoints; i++)
      {
        double THETA = i*(M_PI * (double) 2.0/_npoints);
        const double CT = std::cos(THETA);
        const double ST = std::sin(THETA);
        verts.push_back(T.transform_point(Point3d(CT*R, HEIGHT, ST*R, gpose)));
      }
    }
  }

  // copy the addresses of the computed vertices into 'vertices' 
  vertices.resize(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    vertices[i] = &verts[i];
}
*/

/// Determines whether a point is inside the cylinder; if so, determines the normal
/*
bool CylinderPrimitive::point_inside(CollisionGeometryPtr geom, const Point3d& p, Vector3d& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const double R = _radius;
  const double halfheight = _height*0.5;
  static shared_ptr<Pose3d> P;

  // get the pose for the collision geometry
  shared_ptr<const Pose3d> gpose = geom->get_pose(); 

  // get the pose for this geometry and BV
  shared_ptr<const Pose3d> bpose = get_pose(); 
  assert(!bpose->rpose);

  // setup a new pose
  if (!P)
    P = shared_ptr<Pose3d>(new Pose3d);
  *P = *bpose;
  P->rpose = gpose;

  // transform the point to cylinder space
  Transform3d T = Pose3d::calc_relative_pose(p.pose, P);
  Point3d query = T.transform_point(p);

  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::point_inside() entered" << std::endl;
  FILE_LOG(LOG_COLDET) << "  cylinder radius: " << R << "  half height: " << halfheight << std::endl;
  FILE_LOG(LOG_COLDET) << "query point: " << p << std::endl;
  FILE_LOG(LOG_COLDET) << "query point (cylinder space): " << query << std::endl;

  // determine whether the point is within the height of the cylinder
  double dcaptop = query[Y] - halfheight;
  double dcapbot = -halfheight - query[Y];
  if (dcaptop > 0.0 || dcapbot > 0.0)
  {
    FILE_LOG(LOG_COLDET) << "point outside of cylinder endcaps" << std::endl;
    return false;
  }

  // determine whether the point is within the circular region of the cylinder
  double cdist_sq = query[X]*query[X] + query[Z]*query[Z] - R*R;
  if (cdist_sq > 0.0)
  {
    FILE_LOG(LOG_COLDET) << "point outside circular region of cylinder" << std::endl;
    return false;
  }

  // see whether normal is degenerate and we should reject
  if (std::sqrt(-cdist_sq) < NEAR_ZERO && 
      (dcaptop >= -NEAR_ZERO || dcapbot >= -NEAR_ZERO))
  {
    FILE_LOG(LOG_COLDET) << "can't decide between normal from cap or from annular region (degenerate)" << std::endl;
    return false;
  }

  // determine the normal
  if (-dcaptop < -dcapbot)
  {
    if (-dcaptop < std::sqrt(-cdist_sq))
      normal = Vector3d(0,1,0,P);
    else
      normal = Vector3d::normalize(Vector3d(query[X],(double) 0.0, query[Z],P));
  }
  else
  {
    if (-dcapbot < std::sqrt(-cdist_sq))
      normal = Vector3d(0,-1,0,P);
    else
      normal = Vector3d::normalize(Vector3d(query[X],(double) 0.0, query[Z],P));
  }

//  if (dcaptop >= dcapbot && dcaptop >= -NEAR_ZERO)
//    normal = Vector3d(0,1,0);
//  else if (dcapbot >= dcaptop && dcapbot >= -NEAR_ZERO)
//    normal = Vector3d(0,-1,0);
//  else
//    normal = Vector3d::normalize(Vector3d(query[X],0,query[Z]));

  // transform the normal
  normal = T.inverse_transform_vector(normal);

  FILE_LOG(LOG_COLDET) << " point is inside cylinder" << std::endl;
  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::point_inside() exited" << std::endl;

  return true;
}
*/

/// Computes the signed distance of a point from the cylinder
double CylinderPrimitive::calc_signed_dist(const Point3d& p)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const double INF = std::numeric_limits<double>::max();
  const double R = _radius;
  const double halfheight = _height*0.5;

  // transform the point to cylinder space
  assert(p.pose == get_pose());

  // compute distances from top and bottom of cylinder and main axis 
  double dcaptop = p[Y] - halfheight;
  double dcapbot = -halfheight - p[Y];
  double cdist_sq = p[X]*p[X] + p[Z]*p[Z] - R*R;

  // check whether point is above or below endcaps
  if (dcaptop > 0.0)
  {
    if (cdist_sq > 0.0)
      return std::sqrt(dcaptop*dcaptop + cdist_sq);
    else
      return dcaptop;
  }
  else if (dcapbot > 0.0)
  {
    if (cdist_sq > 0.0)
      return std::sqrt(dcapbot*dcapbot + cdist_sq);
    else
      return dcapbot;
  }

  // point is within endcaps; check to see whether it is outside cylinder 
  if (cdist_sq > 0.0)
    return std::sqrt(cdist_sq);

  // point is within cylinder: find minimum point of interpenetration 
  double dist = std::min(std::min(-dcaptop, -dcapbot), std::sqrt(-cdist_sq));

  return -dist;
}

/// Gets the distance of a point within the cylinder
/**
 * Returns -INF for points outside the cylinder.
 */
double CylinderPrimitive::calc_penetration_depth(shared_ptr<const Pose3d> P, const Point3d& p) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const double INF = std::numeric_limits<double>::max();
  const double R = _radius;
  const double halfheight = _height*0.5;

  // transform the point to cylinder space
  Point3d query = Pose3d::transform_point(P, p);

  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::calc_penetration_depth() entered" << std::endl;
  FILE_LOG(LOG_COLDET) << "  cylinder radius: " << R << "  half height: " << halfheight << std::endl;
  FILE_LOG(LOG_COLDET) << "query point: " << p << std::endl;
  FILE_LOG(LOG_COLDET) << "query point (cylinder space): " << query << std::endl;

  // determine whether the point is within the height of the cylinder
  double dcaptop = query[Y] - halfheight;
  double dcapbot = -halfheight - query[Y];
  if (dcaptop > 0.0 || dcapbot > 0.0)
  {
    FILE_LOG(LOG_COLDET) << "point outside of cylinder endcaps" << std::endl;
    return -INF;
  }

  // determine whether the point is within the circular region of the cylinder
  double cdist_sq = query[X]*query[X] + query[Z]*query[Z] - R*R;
  if (cdist_sq > 0.0)
  {
    FILE_LOG(LOG_COLDET) << "point outside circular region of cylinder" << std::endl;
    return -INF;
  }

  // return the minimum distance
  double dist = std::min(std::min(-dcaptop, -dcapbot), std::sqrt(-cdist_sq));
  FILE_LOG(LOG_COLDET) << "computed penetration depth: " << dist << std::endl;

  return dist;
}

/// Determines the number of intersections between a line and this cylinder
unsigned CylinderPrimitive::intersect_line(shared_ptr<const Pose3d> Px, const Point3d& origin, const Vector3d& dir, double& t0, double& t1) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // create a coordinate system for the cylinder.  In this system, the 
  // cylinder segment center C is the origin and the cylinder axis direction
  // W is the z-axis.  U and V are the other coordinate axis directions.
  // If P = x*U + y*V + z*W, the cylinder is x^2 + y^2 = r^2, where r is the
  // cylinder radius.  The end caps are |z| = h/2, where h is the cylinder
  // height.
  Vector3d U, V, W = Vector3d(0,1,0,Px);
  Vector3d::determine_orthonormal_basis(W, U, V);
  double halfheight = (double) 0.5 * _height;
  double rsqr = _radius*_radius;

  // convert line origin to cylinder coordinates
  Vector3d diff = origin;
  Vector3d P(U.dot(diff), V.dot(diff), W.dot(diff), Px);

  // get the z-value, in cylinder coordinates, of the incoming line's unit-length direction
  double dz = W.dot(dir);

  // check for line parallel to cylinder axis
  if (std::fabs(dz) >= (double) 1.0 - NEAR_ZERO)
  {
    // determine whether line intersects cylinder end disks
    double radialsqrdist = rsqr - P[X]*P[X] - P[Y]*P[Y];
    if (radialsqrdist < (double) 0.0)
      return 0;

    // line intersects the cylinder end disks
    if (dz > (double) 0.0)
    {
      t0 = -P[Z] - halfheight;
      t1 = -P[Z] + halfheight;
    }
    else
    {
      t0 = P[Z] - halfheight;
      t1 = P[Z] + halfheight;
    }

    return 2; 
  }

  // convert incoming line unit-length direction to cylinder coordinates
  Vector3d D(U.dot(dir), V.dot(dir), dz, Px);

  // check whether line is perpendicular to the cylinder axis
  if (std::fabs(D[Z]) <= NEAR_ZERO)
  {
    // check whether line is outside the planes of the cylinder end disks
    if (std::fabs(P[Z]) > halfheight)
      return 0;

    // test intersection of line P+t*D with infinite cylinder x^2 + y^2 = r^2
    // This reduces to computing the roots of a quadratic equation.  If 
    // P = (px,py,pz) and D = (dx,dy,dz), then the quadratic equation is:
    // (dx^2 + dy^2)*t^2 + 2*(px*dx+py*dy)*t + (px^2+py^2-r^2) = 0
    double a0 = P[X]*P[X] + P[Y]*P[Y] - rsqr;
    double a1 = P[X]*D[X] + P[Y]*D[Y];
    double a2 = D[X]*D[X] + D[Y]*D[Y];
    double disc = a1*a1 - a0*a2;

    // quick check for line does not intersect cylinder
    if (disc < (double) 0.0)
      return 0;
    else if (disc > NEAR_ZERO)
    {
      // line intersects cylinder in two places; determine them
      double root = std::sqrt(disc);
      double inv = (double) 1.0 / a2;
      t0 = (-a1 - root)*inv;
      t1 = (-a1 + root)*inv;
      return 2; 
    }
    else
    {
      // line is tangent to the cylinder
      t0 = (-a1/a2);
      return 1;
    }
  }

  // test plane intersections first
  unsigned quantity = 0;
  double inv = (double) 1.0 / D[Z];
  double s0 = (-halfheight - P[Z])*inv;
  double xtmp = P[X] + D[X]*s0;
  double ytmp = P[Y] + D[Y]*s0;
  if (xtmp*xtmp + ytmp*ytmp <= rsqr)
  {
    t0 = s0;
    quantity++;
  }

  double s1 = (halfheight - P[Z])*inv;
  xtmp = P[X] + D[X]*s1;
  ytmp = P[Y] + D[Y]*s1;
  if (xtmp*xtmp + ytmp*ytmp <= rsqr)
  {
    // planar intersection inside the bottom cylinder end disk
    if (quantity == 0)
      t0 = s1;
    else
      t1 = s1;
    quantity++;
  }

  if (quantity == 2)
  {
    // line intersects both top and bottom cylinder end disks
    if (t0 > t1)
      std::swap(t0, t1);
    return 2;
  }

  // if quantity == 1, then the line must intersect the cylinder wall in a 
  // single point somewhere between the end disks.  This case is detected in
  // the following code that tests for intersection between line and cylinder
  // wall.
  double a0 = P[X]*P[X] + P[Y]*P[Y] - rsqr;
  double a1 = P[X]*D[X] + P[Y]*D[Y];
  double a2 = D[X]*D[X] + D[Y]*D[Y];
  double disc = a1*a1 - a0*a2;
  if (disc < (double) 0.0)
    disc = (double) 0.0;

  if (disc > NEAR_ZERO)
  {
    double root = std::sqrt(disc);
    double inv = (double) 1.0/a2;
    double tValue = (-a1 - root)*inv;
    if (s0 <= s1)
    {
      if (s0 <= tValue && tValue <= s1)
      {
        if (quantity == 0)
          t0 = tValue;
        else
          t1 = tValue;
        quantity++;
      }
    }
    else
    {
      if (s1 <= tValue && tValue <= s0)
      {
        if (quantity == 0)
          t0 = tValue;
        else
          t1 = tValue;
        quantity++;
      }
    }

    if (quantity == 2)
    {
      // line intersects one of the cylinder end disks and once on the cylinder
      // wall
      if (t0 > t1)
        std::swap(t0, t1);
      return 2;
    }

    tValue = (-a1 + root)*inv;
    if (s0 <= s1)
    {
      if (s0 <= tValue && tValue <= s1)
      {
        if (quantity == 0)
          t0 = tValue;
        else
          t1 = tValue;
        quantity++;
      }
    }
    else
    {
      if (s1 <= tValue && tValue <= s0)
      {
        if (quantity == 0)
          t0 = tValue;
        else
          t1 = tValue;
        quantity++;
      }
    }
  }
  else
  {
    double tValue = -a1/a2;
    if (s0 <= s1)
    {
      if (s0 <= tValue && tValue <= s1)
      {
        if (quantity == 0)
          t0 = tValue;
        else
          t1 = tValue;
        quantity++;
      }
    }
    else
    {
      if (s1 <= tValue && tValue <= s1)
      {
        if (quantity == 0)
          t0 = tValue;
        else
          t1 = tValue;
        quantity++;
      }
    }
  }

  if (quantity == 2 && t0 > t1)
    std::swap(t0, t1);

  return quantity;  
}

/**
/// Computes the intersection between a cylinder and a line segment
 * Algorithm adapted from www.geometrictools.com 
 * \note for line segments that are partially or fully inside the cylinder, the
 *       method only returns intersection if the second endpoint of the segment
 *       is farther inside than the first
bool CylinderPrimitive::intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  const unsigned Y = 1;
  const double R = _radius;
  const double halfheight = _height*0.5;
  static shared_ptr<Pose3d> P;

  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() entered" << std::endl;
  FILE_LOG(LOG_COLDET) << "  cylinder radius: " << R << "  half height: " << halfheight << std::endl;
  FILE_LOG(LOG_COLDET) << "  segment: " << seg.first << " / " << seg.second << std::endl;

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

  // transform segment to primitive frame 
  Transform3d T = Pose3d::calc_relative_pose(seg.first.pose, P);

  // compute transformed segment
  LineSeg3 nseg(T.transform_point(seg.first), T.transform_point(seg.second));

  // check whether the first point is already within the cylinder
  double p_depth = calc_penetration_depth(P, nseg.first);

  if (p_depth >= (double) 0.0)
  {
    FILE_LOG(LOG_COLDET) << "-- first point is within the cylinder" << endl;
    FILE_LOG(LOG_COLDET) << "  -- depth of first point: " << p_depth << endl;

    // compute the normal
    if (!point_inside(bv, nseg.first, normal))
      return false;

    // set point and time of intersection
    isect = nseg.first;
    t = (double) 0.0;
    normal = T.inverse_transform_vector(normal);

    FILE_LOG(LOG_COLDET) << "  -- point is inside!" << std::endl;
    FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() exited" << std::endl;
    return true;
  }

  // if line segment is a point- we know it's not within the cylinder- and 
  // we have an easy exit 
  Vector3d dir = nseg.second - nseg.first;
  double seg_len = dir.norm();
  if (seg_len < NEAR_ZERO)
  {
    FILE_LOG(LOG_COLDET) << " -- line segment is a point" << std::endl;
    FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() exited" << std::endl;
    return false;
  }

  // do line / cylinder intersection
  double t0, t1;
  unsigned nisects = intersect_line(P, nseg.first, dir/seg_len, t0, t1);

  // if no intersections or not a ray intersection, quit now
  if (nisects == 0 || (nisects == 1 && t0 < (double) 0.0) || 
      (nisects == 2 && t1 < (double) 0.0))
  {
    FILE_LOG(LOG_COLDET) << " -- line segment does not intersect cylinder" << std::endl;
    FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() exited" << std::endl;
    return false;
  }

  // determine the first point of intersection
  if (nisects == 1)
    t = t0;
  else
    t = std::max((double) 0.0, t0);

  // check whether t is within bounds of the segment
  if (t > seg_len)
  {
    FILE_LOG(LOG_COLDET) << " -- line segment does not intersect cylinder" << std::endl;
    FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() exited" << std::endl;
    return false;
  }

  // convert t from line parameter to line segment parameter
  t /= seg_len;
  
  // determine the point of intersection
  isect = seg.first + (seg.second - seg.first)*t;

  // setup cylinder height
  const double HALFH = _height * (double) 0.5;

  // determine the normal in cylinder coordinates
  Point3d isect_cyl = T.transform_point(isect);
  if (isect_cyl[Y] >= HALFH - NEAR_ZERO)
    normal = Vector3d(0,1,0,P);
  else if (isect_cyl[Y] <= -HALFH + NEAR_ZERO)
    normal = Vector3d(0,-1,0,P);
  else
  {
    isect_cyl[Y] = (double) 0.0;
    normal = Vector3d::normalize(isect_cyl);
  }

  // transform normal to global coords
  normal = T.inverse_transform_vector(normal);

  FILE_LOG(LOG_COLDET) << " -- line segment intersects cylinder" << std::endl;
  FILE_LOG(LOG_COLDET) << "    point of intersection (transformed): " << isect << std::endl;
  FILE_LOG(LOG_COLDET) << "    normal (transformed): " << normal << std::endl;
  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() exited" << std::endl;

  return true;  
}
 */

