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
#include <Moby/OBB.h>
#include <Moby/CylinderPrimitive.h>

using namespace Moby;
using boost::shared_ptr;
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
CylinderPrimitive::CylinderPrimitive(Real radius, Real height)
{
  if (height < (Real) 0.0)
    throw std::runtime_error("Attempting to set negative height in CylinderPrimitive (constructor)");
  if (radius < (Real) 0.0)
    throw std::runtime_error("Attempting to set negative radius in CylinderPrimitive (constructor)");

  _radius = radius;
  _height = height;
  _npoints = 10;
  _nrings = 2;
  calc_mass_properties();
}

/// Constructs a cylinder along the y-axis with specified radius and height, centered at the origin, with 10 circle points and 2 rings
CylinderPrimitive::CylinderPrimitive(Real radius, Real height, const Matrix4& T) : Primitive(T)
{
  if (height < (Real) 0.0)
    throw std::runtime_error("Attempting to set negative height in CylinderPrimitive (constructor)");
  if (radius < (Real) 0.0)
    throw std::runtime_error("Attempting to set negative radius in CylinderPrimitive (constructor)");

  _radius = radius;
  _height = height;
  _npoints = 10;
  _nrings = 2;
  calc_mass_properties();
}

/// Constructs a cylinder along the y-axis and centered at the origin with specified, radius, height, number of points and number of rings
CylinderPrimitive::CylinderPrimitive(Real radius, Real height, unsigned n, unsigned nrings, const Matrix4& T) : Primitive(T)
{
  if (height < (Real) 0.0)
    throw std::runtime_error("Attempting to set negative height in CylinderPrimitive (constructor)");
  if (radius < (Real) 0.0)
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

/// Sets the radius for this cylinder
void CylinderPrimitive::set_radius(Real radius)
{
  _radius = radius;
  if (_radius < 0.0)
    throw std::runtime_error("Attempting to set negative radius on call to CylinderPrimitive::set_radius()");

  // mesh, vertices are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices = shared_ptr<vector<Vector3> >();
  _smesh = pair<shared_ptr<IndexedTriArray>, list<unsigned> >();
  _invalidated = true;

  // need to recompute the mass properties of the cylinder
  calc_mass_properties();

  // need up update visualization
  update_visualization();
}

/// Sets the height for this cylinder
void CylinderPrimitive::set_height(Real height)
{
  _height = height;
  if (_height < 0.0)
    throw std::runtime_error("Attempting to set negative height on call to CylinderPrimitive::set_height()");

  // mesh, vertices are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices = shared_ptr<vector<Vector3> >();
  _smesh = pair<shared_ptr<IndexedTriArray>, list<unsigned> >();
  _invalidated = true;
  
  // need to recompute the mass properties of the cylinder
  calc_mass_properties();

  // need up update visualization
  update_visualization();
}

/// Sets the number of points in the circles of the cylinder 
void CylinderPrimitive::set_num_circle_points(unsigned n)
{
  _npoints = n;
  if (n < 3)
    throw std::runtime_error("Attempting to call CylinderPrimitive::set_circle_points() with n < 3");

  // vertices are no longer valid
  _vertices = shared_ptr<vector<Vector3> >();
  _invalidated = true;
}

/// Sets the number of rings for determining "virtual points" (points on the tube of the cylinder)
void CylinderPrimitive::set_num_rings(unsigned n)
{
  _nrings = n;
  if (_nrings < 2)
    throw std::runtime_error("Attempting to call CylinderPrimitive::set_num_rings() with n < 2");

  // vertices are no longer valid
  _vertices = shared_ptr<vector<Vector3> >();
  _invalidated = true;
}

/// Gets the triangle mesh for the cylinder, computing it if necessary 
shared_ptr<const IndexedTriArray> CylinderPrimitive::get_mesh()
{
  // compute the mesh if necessary
  if (!_mesh)
  {
    const Real R = _radius;
    const Real HH = _height;

    // if radius or height = 0, or circle points < 3, return now
    if (_radius <= 0.0 || _height <= 0.0 || _npoints < 3)
    {
      _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray());
      _smesh = make_pair(_mesh, list<unsigned>());
      return _mesh;
    }

    // get the current transform
    const Matrix4& T = get_transform();

    // create vertices evenly spaced in 2D
    std::list<Vector3> points;
    for (unsigned i=0; i< _npoints; i++)
    {
      const Real THETA = i * (M_PI * (Real) 2.0/_npoints);
      const Real CT = std::cos(THETA);
      const Real ST = std::sin(THETA);
      points.push_back(T.mult_point(Vector3(CT*R, HH, ST*R)));
      points.push_back(T.mult_point(Vector3(CT*R, -HH, ST*R)));
    }

    // compute the convex hull
    PolyhedronPtr hull = CompGeom::calc_convex_hull(points.begin(), points.end());

    // get the mesh
    const std::vector<Vector3>& v = hull->get_vertices();
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
#ifdef USE_OSG
osg::Node* CylinderPrimitive::create_visualization()
{
  osg::Cylinder* cyl = new osg::Cylinder;
  cyl->setRadius((float) _radius);
  cyl->setHeight((float) _height);
//  cyl->setRotation(osg::Quat(0, 0, 1, 0));
  cyl->setRotation(osg::Quat(M_PI/2.0f, osg::Vec3f(1.0f, 0.0f, 0.0f)));
  osg::Geode* geode = new osg::Geode;
  geode->addDrawable(new osg::ShapeDrawable(cyl));
  return geode;
}
#endif

/// Implements Base::load_from_xml() for serialization
void CylinderPrimitive::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  // verify that the node type is cylinder
  assert(strcasecmp(node->name.c_str(), "Cylinder") == 0);

  // load the parent data
  Primitive::load_from_xml(node, id_map);

  // read in the radius, if specified
  const XMLAttrib* radius_attr = node->get_attrib("radius");
  if (radius_attr)
    _radius = radius_attr->get_real_value();

  // read in the height, if specified
  const XMLAttrib* height_attr = node->get_attrib("height");
  if (height_attr)
    _height = height_attr->get_real_value();

  // read in the number of circle points, if specified
  const XMLAttrib* npoints_attr = node->get_attrib("num-circle-points");
  if (npoints_attr)
    _npoints = npoints_attr->get_unsigned_value();

  // read in the number of rings of the cylinder, if specified
  const XMLAttrib* nrings_attr = node->get_attrib("num-rings");
  if (nrings_attr)
    _nrings = nrings_attr->get_unsigned_value();

  // recompute mass properties
  calc_mass_properties();
}

/// Implements Base::save_to_xml() for serialization
void CylinderPrimitive::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
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
void CylinderPrimitive::set_transform(const Matrix4& T)
{
  // determine the transformation from the old to the new transform 
  Matrix4 Trel = T * Matrix4::inverse_transform(_T);

  // go ahead and set the new transform
  Primitive::set_transform(T);

  // OBB is no longer valid
  _obb = OBBPtr();

  // transform mesh
  if (_mesh)
  {
    _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(_mesh->transform(Trel)));
    _smesh.first = _mesh;
  }

  // transform vertices
  if (_vertices)
    for (unsigned i=0; i< _vertices->size(); i++)
      (*_vertices)[i] = Trel.mult_point((*_vertices)[i]);

  // invalidate this primitive
  _invalidated = true;

  // recalculate the mass properties
  calc_mass_properties();
}

/// Calculates mass properties for the cylinder
void CylinderPrimitive::calc_mass_properties()
{
  // get the current transform
  const Matrix4& T = get_transform();

  // compute the radius squared (we'll need this)
  const Real RSQ = _radius * _radius;

  // compute the mass if density is given
  if (_density)
  {
    const Real volume = M_PI * RSQ * _height;
    _mass = *_density * volume;
  }

  // compute the non-longitudinal elements
  const Real HSQ = _height * _height;
  const Real ONE_TWELFTH = (Real) 1.0/12.0;
  const Real NL_ELM = ONE_TWELFTH * _mass * (HSQ + (Real) 3.0 * RSQ);
  const Real LONG_ELM = 0.5 * _mass * RSQ;

  // compute the inertia matrix
  Matrix3 J(NL_ELM, 0, 0, 0, LONG_ELM, 0, 0, 0, NL_ELM);

  // transform the inertia matrix using the current transform
  transform_inertia(_mass, J, ZEROS_3, T, _J, _com);
}

/// Gets the OBB
BVPtr CylinderPrimitive::get_BVH_root()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // cylinder not applicable for deformable bodies 
  if (is_deformable())
    throw std::runtime_error("CylinderPrimitive::get_BVH_root() - primitive unusable for deformable bodies!");

  // create the OBB if necessary
  if (!_obb)
    _obb = shared_ptr<OBB>(new OBB);

  // setup the center of the OBB 
  _obb->center = get_transform().get_translation();
  
  // setup the orientation of the OBB
  _obb->R(X,X) = get_transform()(X,X);
  _obb->R(X,Y) = get_transform()(X,Y);
  _obb->R(X,Z) = get_transform()(X,Z);
  _obb->R(Y,X) = get_transform()(Y,X);
  _obb->R(Y,Y) = get_transform()(Y,Y);
  _obb->R(Y,Z) = get_transform()(Y,Z);
  _obb->R(Z,X) = get_transform()(Z,X);
  _obb->R(Z,Y) = get_transform()(Z,Y);
  _obb->R(Z,Z) = get_transform()(Z,Z);

  // must orthonormalize OBB orientation, b/c T may have scaling applied
  _obb->R.orthonormalize();

  // cylinder nominally points upward
  _obb->l[X] = _radius;
  _obb->l[Y] = _height*0.5;
  _obb->l[Z] = _radius;

  return _obb;
}

/// Gets a sub-mesh for the primitive
const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& CylinderPrimitive::get_sub_mesh(BVPtr bv)
{
  if (!_smesh.first)
    get_mesh(); 
  return _smesh;
}

/// Gets vertices from the primitive
void CylinderPrimitive::get_vertices(BVPtr bv, std::vector<const Vector3*>& vertices)
{
  // create the vector of vertices if necessary
  if (!_vertices)
  {
    const Real R = _radius + _intersection_tolerance;
    const Real H = _height + (_intersection_tolerance * (Real) 2.0);

    // if radius or height <= 0, num rings < 2, or circle points < 3, return now
    if (_radius <= 0.0 || _height <= 0.0 || _nrings < 2 || _npoints < 3)
    {
      vertices.clear();
      return; 
    }

    // get the current primitive transform
    const Matrix4& T = get_transform();

    // create the vector of vertices
    _vertices = shared_ptr<vector<Vector3> >(new vector<Vector3>());

    // create vertices evenly spaced in 2D
    for (unsigned j=0; j< _nrings; j++)
    {
      const Real HEIGHT = -(H * (Real) 0.5) + (j*H)/(_nrings-1); 
      for (unsigned i=0; i< _npoints; i++)
      {
        Real THETA = i*(M_PI * (Real) 2.0/_npoints);
        const Real CT = std::cos(THETA);
        const Real ST = std::sin(THETA);
        _vertices->push_back(T.mult_point(Vector3(CT*R, HEIGHT, ST*R)));
      }
    }
  }

  // copy the addresses of the computed vertices into 'vertices' 
  vertices.resize(_vertices->size());
  for (unsigned i=0; i< _vertices->size(); i++)
    vertices[i] = &(*_vertices)[i];
}

/// Determines whether a point is inside the cylinder; if so, determines the normal
bool CylinderPrimitive::point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Real R = _radius;
  const Real halfheight = _height*0.5;

  // transform the point to cylinder space
  const Matrix4& T = get_transform();
  Vector3 query = T.inverse_mult_point(p);

  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::point_inside() entered" << std::endl;
  FILE_LOG(LOG_COLDET) << "  cylinder radius: " << R << "  half height: " << halfheight << std::endl;
  FILE_LOG(LOG_COLDET) << "query point: " << p << std::endl;
  FILE_LOG(LOG_COLDET) << "query point (cylinder space): " << query << std::endl;

  // determine whether the point is within the height of the cylinder
  Real dcaptop = query[Y] - halfheight;
  Real dcapbot = -halfheight - query[Y];
  if (dcaptop > 0.0 || dcapbot > 0.0)
  {
    FILE_LOG(LOG_COLDET) << "point outside of cylinder endcaps" << std::endl;
    return false;
  }

  // determine whether the point is within the circular region of the cylinder
  Real cdist_sq = query[X]*query[X] + query[Z]*query[Z] - R*R;
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
      normal = Vector3(0,1,0);
    else
      normal = Vector3::normalize(Vector3(query[X],(Real) 0.0, query[Z]));
  }
  else
  {
    if (-dcapbot < std::sqrt(-cdist_sq))
      normal = Vector3(0,-1,0);
    else
      normal = Vector3::normalize(Vector3(query[X],(Real) 0.0, query[Z]));
  }
/*
  if (dcaptop >= dcapbot && dcaptop >= -NEAR_ZERO)
    normal = Vector3(0,1,0);
  else if (dcapbot >= dcaptop && dcapbot >= -NEAR_ZERO)
    normal = Vector3(0,-1,0);
  else
    normal = Vector3::normalize(Vector3(query[X],0,query[Z]));
*/
  // transform the normal
  normal = T.mult_vector(normal);

  FILE_LOG(LOG_COLDET) << " point is inside cylinder" << std::endl;
  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::point_inside() exited" << std::endl;

  return true;
}

/// Gets the distance of a point within the cylinder
/**
 * Returns -INF for points outside the cylinder.
 */
Real CylinderPrimitive::calc_penetration_depth(const Vector3& p) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Real INF = std::numeric_limits<Real>::max();
  const Real R = _radius;
  const Real halfheight = _height*0.5;

  // transform the point to cylinder space
  const Matrix4& T = get_transform();
  Vector3 query = T.inverse_mult_point(p);

  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::calc_penetration_depth() entered" << std::endl;
  FILE_LOG(LOG_COLDET) << "  cylinder radius: " << R << "  half height: " << halfheight << std::endl;
  FILE_LOG(LOG_COLDET) << "query point: " << p << std::endl;
  FILE_LOG(LOG_COLDET) << "query point (cylinder space): " << query << std::endl;

  // determine whether the point is within the height of the cylinder
  Real dcaptop = query[Y] - halfheight;
  Real dcapbot = -halfheight - query[Y];
  if (dcaptop > 0.0 || dcapbot > 0.0)
  {
    FILE_LOG(LOG_COLDET) << "point outside of cylinder endcaps" << std::endl;
    return -INF;
  }

  // determine whether the point is within the circular region of the cylinder
  Real cdist_sq = query[X]*query[X] + query[Z]*query[Z] - R*R;
  if (cdist_sq > 0.0)
  {
    FILE_LOG(LOG_COLDET) << "point outside circular region of cylinder" << std::endl;
    return -INF;
  }

  // return the minimum distance
  Real dist = std::min(std::min(-dcaptop, -dcapbot), std::sqrt(-cdist_sq));
  FILE_LOG(LOG_COLDET) << "computed penetration depth: " << dist << std::endl;

  return dist;

}

/// Determines the number of intersections between a line and this cylinder
unsigned CylinderPrimitive::intersect_line(const Vector3& origin, const Vector3& dir, Real& t0, Real& t1) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // form a Matrix4 from the transform
  const Matrix4& T = get_transform();

  // create a coordinate system for the cylinder.  In this system, the 
  // cylinder segment center C is the origin and the cylinder axis direction
  // W is the z-axis.  U and V are the other coordinate axis directions.
  // If P = x*U + y*V + z*W, the cylinder is x^2 + y^2 = r^2, where r is the
  // cylinder radius.  The end caps are |z| = h/2, where h is the cylinder
  // height.
  Vector3 U, V, W = T.mult_vector(Vector3(0,1,0));
  Vector3 C = T.get_translation();
  Vector3::determine_orthonormal_basis(W, U, V);
  Real halfheight = (Real) 0.5 * _height;
  Real rsqr = _radius*_radius;

  // convert line origin to cylinder coordinates
  Vector3 diff = origin - C;
  Vector3 P(U.dot(diff), V.dot(diff), W.dot(diff));

  // get the z-value, in cylinder coordinates, of the incoming line's unit-length direction
  Real dz = W.dot(dir);

  // check for line parallel to cylinder axis
  if (std::fabs(dz) >= (Real) 1.0 - NEAR_ZERO)
  {
    // determine whether line intersects cylinder end disks
    Real radialsqrdist = rsqr - P[X]*P[X] - P[Y]*P[Y];
    if (radialsqrdist < (Real) 0.0)
      return 0;

    // line intersects the cylinder end disks
    if (dz > (Real) 0.0)
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
  Vector3 D(U.dot(dir), V.dot(dir), dz);

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
    Real a0 = P[X]*P[X] + P[Y]*P[Y] - rsqr;
    Real a1 = P[X]*D[X] + P[Y]*D[Y];
    Real a2 = D[X]*D[X] + D[Y]*D[Y];
    Real disc = a1*a1 - a0*a2;

    // quick check for line does not intersect cylinder
    if (disc < (Real) 0.0)
      return 0;
    else if (disc > NEAR_ZERO)
    {
      // line intersects cylinder in two places; determine them
      Real root = std::sqrt(disc);
      Real inv = (Real) 1.0 / a2;
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
  Real inv = (Real) 1.0 / D[Z];
  Real s0 = (-halfheight - P[Z])*inv;
  Real xtmp = P[X] + D[X]*s0;
  Real ytmp = P[Y] + D[Y]*s0;
  if (xtmp*xtmp + ytmp*ytmp <= rsqr)
  {
    t0 = s0;
    quantity++;
  }

  Real s1 = (halfheight - P[Z])*inv;
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
  Real a0 = P[X]*P[X] + P[Y]*P[Y] - rsqr;
  Real a1 = P[X]*D[X] + P[Y]*D[Y];
  Real a2 = D[X]*D[X] + D[Y]*D[Y];
  Real disc = a1*a1 - a0*a2;
  if (disc < (Real) 0.0)
    disc = (Real) 0.0;

  if (disc > NEAR_ZERO)
  {
    Real root = std::sqrt(disc);
    Real inv = (Real) 1.0/a2;
    Real tValue = (-a1 - root)*inv;
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
    Real tValue = -a1/a2;
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

/// Sets the intersection tolerance
void CylinderPrimitive::set_intersection_tolerance(Real tol)
{
  Primitive::set_intersection_tolerance(tol);

  // vertices are no longer valid
  _vertices = shared_ptr<vector<Vector3> >();
}

/// Computes the intersection between a cylinder and a line segment
/**
 * Algorithm adapted from www.geometrictools.com 
 * \note for line segments that are partially or fully inside the cylinder, the
 *       method only returns intersection if the second endpoint of the segment
 *       is farther inside than the first
 */
bool CylinderPrimitive::intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const
{
  const unsigned Y = 1;
  const Real R = _radius;
  const Real halfheight = _height*0.5;

  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() entered" << std::endl;
  FILE_LOG(LOG_COLDET) << "  cylinder radius: " << R << "  half height: " << halfheight << std::endl;
  FILE_LOG(LOG_COLDET) << "  segment: " << seg.first << " / " << seg.second << std::endl;


  // check whether the first point is already within the cylinder
  Real p_depth = calc_penetration_depth(seg.first);

  if (p_depth >= (Real) 0.0)
  {
    FILE_LOG(LOG_COLDET) << "-- first point is within the cylinder" << endl;
    FILE_LOG(LOG_COLDET) << "  -- depth of first point: " << p_depth << endl;

    // compute the normal
    if (!point_inside(bv, seg.first, normal))
      return false;

    // set point and time of intersection
    isect = seg.first;
    t = (Real) 0.0;

    FILE_LOG(LOG_COLDET) << "  -- point is inside!" << std::endl;
    FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() exited" << std::endl;
    return true;
  }

  // form a Matrix4 from the transform
  const Matrix4& T = get_transform();

  // if line segment is a point- we know it's not within the cylinder- and 
  // we have an easy exit 
  Vector3 dir = seg.second - seg.first;
  Real seg_len = dir.norm();
  if (seg_len < NEAR_ZERO)
  {
    FILE_LOG(LOG_COLDET) << " -- line segment is a point" << std::endl;
    FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() exited" << std::endl;
    return false;
  }

  // do line / cylinder intersection
  Real t0, t1;
  unsigned nisects = intersect_line(seg.first, dir/seg_len, t0, t1);

  // if no intersections or not a ray intersection, quit now
  if (nisects == 0 || (nisects == 1 && t0 < (Real) 0.0) || 
      (nisects == 2 && t1 < (Real) 0.0))
  {
    FILE_LOG(LOG_COLDET) << " -- line segment does not intersect cylinder" << std::endl;
    FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() exited" << std::endl;
    return false;
  }

  // determine the first point of intersection
  if (nisects == 1)
    t = t0;
  else
    t = std::max((Real) 0.0, t0);

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
  const Real HALFH = _height * (Real) 0.5;

  // determine the normal in cylinder coordinates
  Vector3 isect_cyl = T.inverse_mult_point(isect);
  if (isect_cyl[Y] >= HALFH - NEAR_ZERO)
    normal = Vector3(0,1,0);
  else if (isect_cyl[Y] <= -HALFH + NEAR_ZERO)
    normal = Vector3(0,-1,0);
  else
  {
    isect_cyl[Y] = (Real) 0.0;
    normal = Vector3::normalize(isect_cyl);
  }

  // transform normal to global coords
  normal = T.mult_vector(normal);

  FILE_LOG(LOG_COLDET) << " -- line segment intersects cylinder" << std::endl;
  FILE_LOG(LOG_COLDET) << "    point of intersection (transformed): " << isect << std::endl;
  FILE_LOG(LOG_COLDET) << "    normal (transformed): " << normal << std::endl;
  FILE_LOG(LOG_COLDET) << "CylinderPrimitive::intersect_seg() exited" << std::endl;

  return true;  
}

