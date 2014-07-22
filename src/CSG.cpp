/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifdef USE_OSG
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#include <osg/Geometry>
#endif
#include <boost/algorithm/string.hpp>
#include <Ravelin/LinAlgd.h>
#include <Moby/XMLTree.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/CSG.h>

using std::endl;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;
using std::make_pair;
using std::list;
using std::vector;
using std::pair;
using std::string;
using std::cerr;
using boost::dynamic_pointer_cast;

/// Constructs a union CSG by default
CSG::CSG()
{
  _op = eUnion;
  calc_mass_properties();
}

/// Constructs a union CSG with the given transform
CSG::CSG(const Pose3d& T) : Primitive(T)
{
  _op = eUnion;
  calc_mass_properties();
}

void CSG::set_operand1(PrimitivePtr p)
{
  _op1 = p;

  // computed mesh and vertices are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices.clear();
  _smesh = pair<shared_ptr<IndexedTriArray>, list<unsigned> >();
  _invalidated = true;

  // need to recalculate mass properties
  calc_mass_properties();
}

void CSG::set_operand2(PrimitivePtr p)
{
  _op2 = p;

  // computed mesh and vertices are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices.clear();
  _smesh = pair<shared_ptr<IndexedTriArray>, list<unsigned> >();
  _invalidated = true;

  // need to recalculate mass properties
  calc_mass_properties();
}

void CSG::set_operator(BooleanOperation op)
{
  _op = op;

  // computed mesh and vertices are no longer valid
  _mesh = shared_ptr<IndexedTriArray>();
  _vertices.clear();
  _smesh = pair<shared_ptr<IndexedTriArray>, list<unsigned> >();
  _invalidated = true;

  // need to recalculate mass properties
  calc_mass_properties();
}

/// Computes the distance and normal from a point on the CSG
double CSG::calc_dist_and_normal(const Point3d& p, Vector3d& normal) const
{
  // TODO: implement this
  assert(false);
  return 0.0;
}

/// Finds the signed distance between the sphere and another primitive
double CSG::calc_signed_dist(shared_ptr<const Primitive> primitive, shared_ptr<const Pose3d> pose_this, shared_ptr<const Pose3d> pose_primitive, Point3d& pthis, Point3d& pprimitive) const
{
  // TODO: implement this
  assert(false);
  return 0.0; 
}

/*
/// Determines whether a point is inside or on the CSG
bool CSG::point_inside(CollisionGeometryPtr geom, const Point3d& p, Vector3d& normal) const
{
  Vector3d normal1, normal2;

  if (!_op1 || !_op2)
    return false;

  // see whether the point is inside or on one or both primitives
  bool in1 = _op1->point_inside(geom, p, normal1);
  bool in2 = _op2->point_inside(geom, p, normal2);

  // see what to return
  if (_op == eUnion)
  {
    if (!in1 && !in2)
      return false;
    if (!in1)
      normal = normal2;
    else if (!in2)
      normal = normal1;
    else
      // pick a normal arbitrarily
      normal = normal1;

    return true;
  }
  else if (_op == eIntersection)
  {
    if (!in1 || !in2)
      return false;

    // pick a normal arbitrarily
    normal = normal1;
    return true;
  }
  else
  {
    assert(_op == eDifference);
    if (!(in1 && !in2))
      return false;
    normal = normal1;
    return true;
  }
}
*/

/// Transforms the CSG
void CSG::set_pose(const Pose3d& p)
{
  // convert p to a shared pointer
  shared_ptr<Pose3d> x(new Pose3d(p));

  // determine the transformation from the old pose to the new one 
  Transform3d T = Pose3d::calc_relative_pose(_F, x);

  // call the primitive transform
  Primitive::set_pose(p);

  // invalidate vertices and the mesh 
  _vertices.clear();
  _mesh.reset();

  // invalidate this primitive
  _invalidated = true;

  // recalculate the mass properties
  calc_mass_properties();
}

/*
/// Computes distance between a CSG and a primitive
double CSG::calc_dist(Primitive* p, Point3d& pcsg, Point3d& pp)
{
  switch (_op)
  {
    case eUnion:
      // compute distance from each operand to p

      // if both distances >= 0, return minimum distance

      // else return negative distance closest to zero

    case eIntersection:

    case eDifference:
  }
}
*/

/// Gets the bounding volume
BVPtr CSG::get_BVH_root(CollisionGeometryPtr geom)
{
  const unsigned THREE_D = 3;

  // CSG not applicable for deformable bodies 
  if (is_deformable())
    throw std::runtime_error("CSG::get_BVH_root(CollisionGeometryPtr geom) - primitive unusable for deformable bodies!");

  // return no BVH if either operand is not set
  if (!_op1 || !_op2)
    return BVPtr();

  // we always recalc the AABB, as the underlying primitives can change
  if (!_aabb)
    _aabb = shared_ptr<AABB>(new AABB);

// TODO: pose should be taken into account
  // get the bounding volumes for the operands
  shared_ptr<BV> bv1 = _op1->get_BVH_root(geom);
  shared_ptr<BV> bv2 = _op2->get_BVH_root(geom);

  // get the lower and upper bounds for the operands
  Point3d low1 = bv1->get_lower_bounds();
  Point3d low2 = bv2->get_lower_bounds();
  Point3d high1 = bv1->get_upper_bounds();
  Point3d high2 = bv2->get_upper_bounds();

  // determine the overall upper and lower bounds
  Point3d low, high;
  
  if (_op == eUnion)
  {
    // AABB is set to the union of the AABB's
    for (unsigned i=0; i< THREE_D; i++) 
    {
      low[i] = std::min(low1[i], low2[i]);
      high[i] = std::max(high1[i], high2[i]);
    }
  }
  else if (_op == eDifference)
  {
    // AABB is set to the difference of the AABB's
    // NOTE: this is conservative
   low = low1;
   high = high1;
  }
  else if (_op == eIntersection)
  {
    // AABB is set to the intersection of the AABB's
    for (unsigned i=0; i< THREE_D; i++)
    {
      low[i] = std::max(low1[i], low2[i]);
      high[i] = std::min(high1[i], high2[i]);
    }
  }

  // set the bounds
  _aabb->minp = low;
  _aabb->maxp = high;
  _aabb->geom = geom;

  return _aabb;
}

/*
bool CSG::intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  static shared_ptr<Pose3d> P;

  if (!_op1 || !_op2)
    throw std::runtime_error("One or more CSG operands are missing!");

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

  // compute transform from line segment pose to primitive pose
  Transform3d T = Pose3d::calc_relative_pose(seg.first.pose, P);

  // compute updated line segment
  LineSeg3 nseg(T.transform_point(seg.first), T.transform_point(seg.second));

  // reset intersection tolerances (if necessary)
  if (std::fabs(_op1->get_intersection_tolerance() - _intersection_tolerance) > NEAR_ZERO)
    _op1->set_intersection_tolerance(_intersection_tolerance);
  if (std::fabs(_op2->get_intersection_tolerance() - _intersection_tolerance) > NEAR_ZERO)
    _op2->set_intersection_tolerance(_intersection_tolerance);

  bool flag = false;
  switch (_op)
  {
    case eUnion: 
      flag = intersect_seg_union(bv, nseg, t, isect, normal);
      break;

    case eIntersection: 
      flag = intersect_seg_intersect(bv, nseg, t, isect, normal);
      break;

    case eDifference: 
      flag = intersect_seg_diff(bv, nseg, t, isect, normal);
      break;
  }

  // transform intersection point and normal (if necessary)
  if (flag)
  {
    isect = Pose3d::transform_point(seg.first.pose, isect); 
    normal = Pose3d::transform_vector(seg.first.pose, normal);
  }

  return flag;
}

/// Does line segment intersection on a union CSG
bool CSG::intersect_seg_union(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  // compute the first time of intersection for each of the two primitives
  double t1, t2;
  Point3d isect1, isect2;
  Vector3d normal1, normal2;

  // get the root bounding volumes
  BVPtr bv1 = _op1->get_BVH_root(bv->geom);
  BVPtr bv2 = _op2->get_BVH_root(bv->geom);

  // perform the individual intersections
  bool i1 = _op1->intersect_seg(bv1, seg, t1, isect1, normal1);
  bool i2 = _op2->intersect_seg(bv2, seg, t2, isect2, normal2);
  if (!i1 && !i2)
    return false;
  else if (i1 && !i2)
  {
    t = t1;
    isect = isect1;
    normal = normal1;
    return true;
  }
  else if (!i1 && i2)
  {
    t = t2;
    isect = isect2;
    normal = normal2;
    return true;
  }
  else
  {
    if (t1 < t2)
    {
      t = t1;
      isect = isect1;
      normal = normal1;
      return true;
    }
    else
    {
      t = t2;
      isect = isect2;
      normal = normal2;
      return true;
    }
  }
}

/// Does line segment intersection on an intersection 
bool CSG::intersect_seg_intersect(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  // compute the first time of intersection for each of the two primitives
  double t1, t2;
  Vector3d normal1, normal2;
  Point3d isect1, isect2;

  // get the root bounding volumes
  BVPtr bv1 = _op1->get_BVH_root(bv->geom);
  BVPtr bv2 = _op2->get_BVH_root(bv->geom);

  // perform the individual intersections
  bool i1 = _op1->intersect_seg(bv1, seg, t1, isect1, normal1);
  if (!i1)
    return false;

  // perform the individual intersection
  bool i2 = _op2->intersect_seg(bv2, seg, t2, isect2, normal2);
  if (!i2)
    return false;

  // need to return the maximum of the t's
  if (t1 < t2)
  {
    t = t2;
    isect = isect2;
    normal = normal2;
  }
  else
  {
    t = t1;
    isect = isect1;
    normal = normal1;
  }

  return true;
}

/// Does line segment intersection on a difference 
 * \note the approach that I am using to do line segment intersection will not
 *       work properly if one of the operands is not convex!  Proper treatment
 *       would require changing the signature for Primitive::intersect_seg()
 *       to return the exit point (if any) of the segment.
bool CSG::intersect_seg_diff(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  // compute the first time of intersection for each of the two primitives
  double t1, t2;
  Point3d isect1, isect2;
  Vector3d normal1, normal2;

  // get the root bounding volumes
  BVPtr bv1 = _op1->get_BVH_root(bv->geom);
  BVPtr bv2 = _op2->get_BVH_root(bv->geom);

  FILE_LOG(LOG_COLDET) << "CSG::intersect_seg_diff() entered" << endl;
  FILE_LOG(LOG_COLDET) << "  segment (CSG space): " << seg.first << " / " << seg.second << endl;

  // setup the intersection tolerance for both primitives
  _op1->set_intersection_tolerance((double) 0.0);
  _op2->set_intersection_tolerance((double) 0.0);

  // special case: look for interpenetration through primitive _op2 
  // specifically, we know that if the point is inside the fattened op2 and
  // is *not* inside the non-fattened op2 and the point is inside op1 then
  // it is interpenetrating through primitive _op2
  if (_op1->point_inside(bv, seg.first, normal) && !_op2->point_inside(bv, seg.first, normal))
  {
    _op2->set_intersection_tolerance(_intersection_tolerance);
    if (_op2->point_inside(bv, seg.first, normal))
    {
      FILE_LOG(LOG_COLDET) << "  point is in special case!" << endl;
      t = (double) 0.0;
      normal = -normal;
      isect = seg.first;
      return true;
    }
  }

  // setup the intersection tolerance for both primitives
  _op1->set_intersection_tolerance(_intersection_tolerance);
  _op2->set_intersection_tolerance(_intersection_tolerance);

  // do intersection with the first primitive 
  bool i1 = _op1->intersect_seg(bv1, seg, t1, isect1, normal1);
  if (!i1)
  {
    FILE_LOG(LOG_COLDET) << "  segment does not intersect first primitive" << endl;
    FILE_LOG(LOG_COLDET) << "  -- no intersection" << endl;
    return false;
  }

  // perform the individual intersection
  bool i2 = _op2->intersect_seg(bv2, seg, t2, isect2, normal2);

  // if first intersection on first primitive, return that
  if (!i2 || t1 < t2)
  {
    FILE_LOG(LOG_COLDET) << "  segment intersects first primitive; i2=" << i2 << " t1=" << t1 << " t2=" << t2 << endl;
    FILE_LOG(LOG_COLDET) << "  -- reporting intersection" << endl;
    t = t1;
    normal = normal1;
    isect = isect1;
    return true;
  }

  // determine first intersections with reversed segment
  double t1R, t2R;
  LineSeg3 qp(seg.second, seg.first);
  i1 = _op1->intersect_seg(bv1, qp, t1R, isect1, normal1);
  i2 = _op2->intersect_seg(bv2, qp, t2R, isect2, normal2);

  FILE_LOG(LOG_COLDET) << "  segment intersects first primitive, but second primitive intersected first" << endl;
  FILE_LOG(LOG_COLDET) << "  -- intersecting reversed segment" << endl;
  FILE_LOG(LOG_COLDET) << "  -- i1=" << i1 << "  t1R=" << t1R << "  isect point: " << isect1 << endl;
  FILE_LOG(LOG_COLDET) << "  -- i2=" << i2 << "  t2R=" << t2R << "  isect point: " << isect2 << endl;

  // i1 and i2 *should* be both true, but we are dealing with numerical
  // algorithms...
  if (!(i1 && i2))
  {
    FILE_LOG(LOG_COLDET) << "  reversed segment did not intersect both primitives!" << endl;
    return false;
  }

  // invert t1R and t2R
  assert(t1R <= (double) 1.0);
  assert(t2R <= (double) 1.0);
  t1R = (double) 1.0 - t1R;
  t2R = (double) 1.0 - t2R;
  if (t1R < (double) 0.0)
    t1R = (double) 0.0;
  if (t2R < (double) 0.0)
    t2R = (double) 0.0;

  // if t2R happens at/after t1R, the line segment never intersects the CSG
  if (t2R + NEAR_ZERO >= t1R)
  {
    FILE_LOG(LOG_COLDET) << "  reversed segment intersects second primitive last; no intersection" << endl;
    return false;
  }

  // at this point, we know 
  // 1. the line segment intersects both (expanded) operands
  // 2. the line segment intersects the second primitive no later than it
  //    intersects the first primitive
  // 3. the line segment exited the first primitive after it exited the
  //    the second primitive (NOTE: this statement has additional ramifications
  //    if both primitives are convex)
  // 

  // if q is outside of the expanded op2, then we can just project the
  // point of intersection back by epsilon units in the direction of [q, p]
  if (t2R < (double) 1.0 - NEAR_ZERO)
  {
    // determine the new point of intersection
    Vector3d qp = seg.second - seg.first;
    const double qp_len = qp.norm();
    assert(qp_len > NEAR_ZERO);
    isect2 -= qp * (_intersection_tolerance / qp_len);
    
    // determine the new time of intersection
    t = (isect2 - seg.first).norm() / qp_len;

    // if t < 0, set t = 0 (and correct isect2 as well)
    if ((isect2-seg.first).dot(qp) < (double) 0.0)
    {
      t = (double) 0.0;
      isect2 = seg.first;
    }

    FILE_LOG(LOG_COLDET) << "  second endpoint of line segment is outside of expanded op2" << endl;
    FILE_LOG(LOG_COLDET) << "  new point of intersection: " << isect2 << endl;
    FILE_LOG(LOG_COLDET) << "  new time of intersection: " << t << endl;
   }
  else
  {
    // q is within the expanded op2; we need to reset op2 intersection 
    // tolerance to zero for next steps
    _op2->set_intersection_tolerance((double) 0.0);

    // if p is outside of _op2 (non-expanded), then point of intersection is p
    Vector3d dummy_normal;
    if (!_op2->point_inside(bv2, seg.first, dummy_normal))
    {
      isect2 = seg.first;
      t = (double) 0.0;
      // NOTE: normal will not have changed

      FILE_LOG(LOG_COLDET) << "  line segment wholly within epsilon zone" << endl;
      FILE_LOG(LOG_COLDET) << "  new point of intersection: " << isect2 << endl;
      FILE_LOG(LOG_COLDET) << "  new time of intersection: " << t << endl;
    }
    else
    {
      // need to find intersection between non-expanded _op2 and pq
      i2 = _op2->intersect_seg(bv2, seg, t, isect2, dummy_normal);

      // due to cases above, we should always have an intersection
      assert(i2);

      // normal from both expanded and base primitives should be the same
      assert(dummy_normal.dot(normal2) < NEAR_ZERO);

      FILE_LOG(LOG_COLDET) << "  line segment crosses epsilon zone" << endl;
      FILE_LOG(LOG_COLDET) << "  new point of intersection: " << isect2 << endl;
      FILE_LOG(LOG_COLDET) << "  new time of intersection: " << t << endl;
    }

    // reset intersection tolerance for op2
    _op2->set_intersection_tolerance(_intersection_tolerance);
  }

  // setup point of intersection and normal
  isect = isect2;
  normal = -normal2;

  FILE_LOG(LOG_COLDET) << "  reversed segment intersects second primitive first" << endl;
  FILE_LOG(LOG_COLDET) << "  -- reporting intersection at " << t << endl;

  return true;
}
*/

/// Sets the mesh for this primitive
/**
 * Using a precomputed mesh is much faster than my algorithms for computing
 * mesh for a CSG.
 */
void CSG::set_mesh(shared_ptr<const IndexedTriArray> mesh)
{
  // transform the mesh using the current transform
  shared_ptr<const Pose3d> P = get_pose();
  Transform3d T = Pose3d::calc_relative_pose(GLOBAL, P);
  _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(mesh->transform(T)));

  // setup sub mesh (it will be just the standard mesh)
  list<unsigned> all_tris;
  for (unsigned i=0; i< _mesh->num_tris(); i++)
    all_tris.push_back(i);
  _smesh = make_pair(_mesh, all_tris);

  // update the visualization
  update_visualization();

  // (re)calculate the mass properties
  calc_mass_properties();
}

/// Creates the visualization for this primitive
osg::Node* CSG::create_visualization()
{
  #ifdef USE_OSG
  const unsigned X = 0, Y = 1, Z = 2;

  // compute the transform from the global frame to the current pose
  shared_ptr<const Pose3d> P = get_pose();
  Transform3d T = Pose3d::calc_relative_pose(GLOBAL, P);

  // back the transformation out of the mesh (if any); NOTE: we have to do this
  // b/c the base Primitive class uses the transform in the visualization
  IndexedTriArray mesh;
  if (_mesh)
    mesh = _mesh->transform(T);

  // create a new group to hold the geometry
  osg::Group* group = new osg::Group;
  osg::Geode* geode = new osg::Geode;
  osg::Geometry* geom = new osg::Geometry;
  geode->addDrawable(geom);
  group->addChild(geode);

  // get the vertices and facets
  const std::vector<Origin3d>& verts = mesh.get_vertices();
  const std::vector<IndexedTri>& facets = mesh.get_facets();
  
  // setup the vertices 
  osg::Vec3Array* varray = new osg::Vec3Array(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    (*varray)[i] = osg::Vec3((float) verts[i][X], (float) verts[i][Y], (float) verts[i][Z]);
  geom->setVertexArray(varray);

  // create the faces
  for (unsigned i=0; i< facets.size(); i++)
  {
    osg::DrawElementsUInt* face = new osg::DrawElementsUInt(osg::PrimitiveSet::TRIANGLES, 0);
    face->push_back(facets[i].a);
    face->push_back(facets[i].b);
    face->push_back(facets[i].c);
    geom->addPrimitiveSet(face);
  }

  return group;
  #else
  return NULL;
  #endif
}

/// Gets a sub-mesh for the primitive
const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& CSG::get_sub_mesh(BVPtr bv)
{
  if (!_smesh.first)
    get_mesh(); 
  return _smesh;
}

/// Returns the mesh (performs a Boolean operation on two CSG operands, if necessary)
shared_ptr<const IndexedTriArray> CSG::get_mesh()
{
  // compute the mesh if necessary 
  if (!_mesh)
  {
    // if one of the operands is not set, exit now
    if (!_op1 || !_op2)
      return shared_ptr<const IndexedTriArray>();

    // ugh...  we have to compute the mesh.  This will be slow.
    // get the meshes from the two operands
    shared_ptr<const IndexedTriArray> m1 = _op1->get_mesh();
    shared_ptr<const IndexedTriArray> m2 = _op2->get_mesh();

    // construct Polyhedra from the meshes
    Polyhedron p1(*m1);
    Polyhedron p2(*m2);

    IndexedTriArray result;
    if (_op == eIntersection)
      result = Polyhedron::construct_intersection(p1, p2);
    else if (_op == eUnion)
      result = Polyhedron::construct_union(p1, p2);
    else if (_op == eDifference)
      result = Polyhedron::construct_difference(p1, p2);
    else
      assert(false);

    // merge m1 and m2; set the mesh
    set_mesh(shared_ptr<IndexedTriArray>(new IndexedTriArray(result)));
  }
    
  return _mesh;
}

/// Centers the underlying triangle mesh
void CSG::center_mesh()
{
  // if there is no mesh, quit now
  if (!_mesh)
    return;

  // we want to transform the mesh so that the inertial pose is coincident
  // with the primitive pose; first, prepare the transformation
  Transform3d T = Pose3d::calc_relative_pose(_jF, _F);

  // do the transformation
  _mesh = shared_ptr<IndexedTriArray>(new IndexedTriArray(_mesh->transform(T)));

  // re-calculate mass properties
  calc_mass_properties();

  // verify that the inertial frame is near zero
  assert(Point3d(_jF->x, GLOBAL).norm() < NEAR_ZERO);

  // update the visualization
  update_visualization();
}

/// Implements Base::load_from_xml() for serialization
void CSG::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  std::map<std::string, BasePtr>::const_iterator id_iter;

  // verify that the node type is sphere
  assert(strcasecmp(node->name.c_str(), "CSG") == 0);

  // load the parent data
  Primitive::load_from_xml(node, id_map);

  // get operand 1, if present
  XMLAttrib* op1_attrib = node->get_attrib("op1-id");
  if (op1_attrib)
  {
    // get the ID of operand 1
    const std::string& id = op1_attrib->get_string_value();

    // find operand 1
    if ((id_iter = id_map.find(id)) == id_map.end())
    {
      std::cerr << "CSG::load_from_xml() - could not find primitive w/ID " << id << std::endl;
      std::cerr << "  from offending node: " << std::endl << *node;
    }
    else
    {
      // make sure that it is castable to a Primitive
      PrimitivePtr p = dynamic_pointer_cast<Primitive>(id_iter->second);
      if (p)
        _op1 = p;
    }
  }

  // get operand 2, if present
  XMLAttrib* op2_attrib = node->get_attrib("op2-id");
  if (op2_attrib)
  {
    // get the ID of operand 2
    const std::string& id = op2_attrib->get_string_value();

    // find operand 2
    if ((id_iter = id_map.find(id)) == id_map.end())
    {
      std::cerr << "CSG::load_from_xml() - could not find primitive w/ID " << id << std::endl;
      std::cerr << "  from offending node: " << std::endl << *node;
    }
    else
    {
      // make sure that it is castable to a Primitive
      PrimitivePtr p = dynamic_pointer_cast<Primitive>(id_iter->second);
      if (p)
        _op2 = p;
    }
  }

  // get operator, if present
  XMLAttrib* op_attrib = node->get_attrib("operator");
  if (op_attrib)
  {
    // get the string
    std::string op_string = op_attrib->get_string_value();

    // look for union
    if (boost::iequals(op_string, std::string("union")))
      _op = eUnion;
    else if (boost::iequals(op_string, std::string("intersection")))
      _op = eIntersection;
    else if (boost::iequals(op_string, std::string("difference")))
      _op = eDifference;
    else
      std::cerr << "CSG::load_from_xml() - unrecognized Boolean operation type '" << op_string << "' in offending node: " << std::endl << *node;
  }

  // determine whether we are using a loaded mesh
  XMLAttrib* fname_attrib = node->get_attrib("triangle-mesh-filename");
  if (fname_attrib)
  {
    // setup the file extensions
    const char* OBJ_EXT = ".obj";

    // get the filename
    string fname(fname_attrib->get_string_value());

    // get the lowercase version of the filename
    string fname_lower = fname;
    std::transform(fname_lower.begin(), fname_lower.end(), fname_lower.begin(), (int(*)(int)) std::tolower);

    // get the type of file and construct the triangle mesh appropriately
    if (fname_lower.find(string(OBJ_EXT)) == fname_lower.size() - strlen(OBJ_EXT))
      set_mesh(shared_ptr<IndexedTriArray>(new IndexedTriArray(IndexedTriArray::read_from_obj(fname))));
    else
    {
      cerr << "CSG::load_from_xml() - unrecognized filename extension" << endl;
      cerr << "  for attribute 'filename'.  Valid extensions are '.obj' (Wavefront OBJ)" << endl;
    }
  
    // see whether to center the mesh
    XMLAttrib* center_attr = node->get_attrib("center");
    if (center_attr && center_attr->get_bool_value())
      this->center_mesh();
  }

  // (re) calculate mass properties
  calc_mass_properties();

  // update visualization
  update_visualization();
}

/// Implements Base::save_to_xml() for serialization
void CSG::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data
  Primitive::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "CSG";

  // save the operator
  if (_op == eUnion)
    node->attribs.insert(XMLAttrib("operator", "union"));
  else if (_op == eIntersection)
    node->attribs.insert(XMLAttrib("operator", "intersection"));
  else if (_op == eDifference)
    node->attribs.insert(XMLAttrib("operator", "difference"));
  else
    assert(false);

  // save the operands
  node->attribs.insert(XMLAttrib("op1-id", _op1->id));
  node->attribs.insert(XMLAttrib("op2-id", _op2->id));

  // save the mesh if it is present
  if (_mesh)
  {
    // make a filename using "this"
    unsigned max_digits = sizeof(void*) * 2;
    char buffer[max_digits+1];
    sprintf(buffer, "%p", this);
    string filename = "triarray" + string(buffer) + ".obj";

    // add the filename as an attribute
    node->attribs.insert(XMLAttrib("triangle-mesh-filename", filename));

    // do not save the array to the OBJ file if it already exists (which we
    // crudely check for using std::ifstream to avoid OS-specific calls -- note
    // that it is possible that opening a file may fails for other reasons than
    // the file does not exist)
    std::ifstream in(filename.c_str());
    if (in.fail())
    {
      // transform the mesh to the CSG frame
      shared_ptr<const Pose3d> P = get_pose();
      Transform3d T = Pose3d::calc_relative_pose(GLOBAL, P);
      IndexedTriArray mesh_xform = _mesh->transform(T);

      // write the mesh
      mesh_xform.write_to_obj(filename);
    }
    else
      in.close();
  }
}

/// Calculates the mass properties for this CSG
void CSG::calc_mass_properties()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // don't calculate mass properties if one or more operand is not set
  if (!_op1 || !_op2)
    return;

  // get the mesh
  const IndexedTriArray& mesh = *get_mesh();

  // get triangles
  std::list<Triangle> tris;
  mesh.get_tris(std::back_inserter(tris), GLOBAL);

  // compute the centroid of the triangle mesh
  _jF->x = Origin3d(CompGeom::calc_centroid_3D(tris.begin(), tris.end()));

  // calculate volume integrals
  double volume_ints[10];
  mesh.calc_volume_ints(volume_ints);

  // we'll need the volume
  const double volume = volume_ints[0];

  // compute the mass if density is given
  if (_density)
    _J.m = *_density * volume;

// NOTE: we no longer transform the inertial components, since we recompute
  // compute the center-of-mass
  // NOTE: we use this instead of the COM calculated from calc_volume_ints()
  // b/c the latter seems to be less accurate; NOTE: we need to check to see
  // why that is the case
  // Point3d com(_volume_ints[1], _volume_ints[2], _volume_ints[3]);
  // com *= (1.0/volume);
  //  Point3d com = calc_com();

  // compute the inertia tensor relative to world origin
  _J.J(X,X) = (_J.m/volume) * (volume_ints[5] + volume_ints[6]);
  _J.J(Y,Y) = (_J.m/volume) * (volume_ints[4] + volume_ints[6]);
  _J.J(Z,Z) = (_J.m/volume) * (volume_ints[4] + volume_ints[5]);
  _J.J(X,Y) = (-_J.m/volume) * volume_ints[7];
  _J.J(Y,Z) = (-_J.m/volume) * volume_ints[8];
  _J.J(X,Z) = (-_J.m/volume) * volume_ints[9];

  // set the symmetric values for J
  _J.J(Y,X) = _J.J(X,Y);
  _J.J(Z,Y) = _J.J(Y,Z);
  _J.J(Z,X) = _J.J(X,Z);

  // if one or more values of J is NaN, don't verify anything
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=i; j<= Z; j++)
      if (std::isnan(_J.J(i,j)))
      {
        _J.J.set_zero();
        return;
      }

  // verify that the inertial properties are correct
  LinAlgd LA;
  if (!LA.is_SPSD(_J.J, std::sqrt(std::numeric_limits<double>::epsilon())))
  {
    std::cerr << "CSG::calc_mass_properties() warning() - inertial properties ";
    std::cerr << "are bad" << std::endl;
    std::cerr << "  inertia tensor: " << std::endl << _J;
    std::cerr << "  likely due to non-manifold mesh" << std::endl;
    _J.J.set_zero();
  }
}

/// Determines the set of vertices on the surface of the CSG
void CSG::get_vertices(vector<Point3d>& verts)
{
  Vector3d dummy;
  vector<Point3d> v1, v2;

  // TODO: finish implementing this

/*
  // if one of the operands is not set, quit
  if (!_op1 || !_op2)
    return;

  // if union, add vertices from both
  if (_op == eUnion)
  {
    // get the vertices from the two operands
    _op1->get_vertices(v1);
    _op2->get_vertices(v2);
  }
  // intersection: return vertices of v1 inside v2 and vice versa
  else if (_op == eIntersection)
  {
    // get the vertices from the two operands
    _op1->get_vertices(v1);
    _op2->get_vertices(v2);
  }
  else
  {
    assert(_op == eDifference);

    // get the vertices from op1 first 
    _op1->get_vertices(v1);

    // look for vertices of v1 that are not inside v2
    for (unsigned i=0; i< v1.size(); i++)
      if (!_op2->point_inside(v1[i], dummy))
        verts.push_back(v1[i]);

    // get the vertices from op2 
    _op2->get_vertices(v2);

    // look for vertices of v2 that are inside v1
    for (unsigned i=0; i< v2.size(); i++)
      if (_op1->point_inside(v2[i], dummy))
        verts.push_back(v2[i]);
  }
*/
}

/*
void CSG::get_vertices(CollisionGeometryPtr geom, vector<const Point3d*>& vertices)
{
  Vector3d dummy;
  vector<const Point3d*> v1, v2;

// TODO: verify this works with poses
  // if one of the operands is not set, quit
  if (!_op1 || !_op2)
    return;

  // get the vertices for this geometry
  vector<Point3d>& verts = _vertices[geom];

  // if either operand is invalidated, clear the vector of vertices
  if (_op1->_invalidated || _op2->_invalidated)
    verts.clear();

  // if the vector of vertices is NULL, compute it
  if (verts.empty())
  {
    // get the pose for the geometry
    shared_ptr<const Pose3d> gpose = geom->get_pose();

    // create a new pose, which will be defined relative to the geometry's base
    shared_ptr<Pose3d> T(new Pose3d);
    *T = *get_pose();
    T->update_relative_pose(gpose);

    // get the BVHs
    BVPtr bv1 = _op1->get_BVH_root(geom);
    BVPtr bv2 = _op2->get_BVH_root(geom);

    // if union, add vertices from both
    if (_op == eUnion)
    {
      // get the vertices from the two operands
      _op1->get_vertices(bv1, v1);
      _op2->get_vertices(bv2, v2);

      // transform the vertices
      BOOST_FOREACH(const Point3d* v, v1)
        verts.push_back(*v);
      BOOST_FOREACH(const Point3d* v, v2)
        verts.push_back(*v);
    }
    // intersection: return vertices of v1 inside v2 and vice versa
    else if (_op == eIntersection)
    {
      // get the vertices from the two operands
      _op1->get_vertices(bv1, v1);
      _op2->get_vertices(bv2, v2);

      // look for vertices of v2 inside op1
      for (unsigned i=0; i< v2.size(); i++)
        if (_op1->point_inside(bv1, *v2[i], dummy))
          verts.push_back(*v2[i]);

      // look for vertices of v1 inside op2
      for (unsigned i=0; i< v1.size(); i++)
        if (_op2->point_inside(bv2, *v1[i], dummy))
          verts.push_back(*v1[i]);
    }
    else
    {
      assert(_op == eDifference);

      // get the vertices from op1 first 
      _op1->get_vertices(bv1, v1);

      // look for vertices of v1 that are not inside v2
      for (unsigned i=0; i< v1.size(); i++)
        if (!_op2->point_inside(bv2, *v1[i], dummy))
          verts.push_back(*v1[i]);

      // get the vertices from op2 
      _op2->set_intersection_tolerance((double) 0.0);
      _op2->get_vertices(bv2, v2);

      // look for vertices of v2 that are inside v1
      for (unsigned i=0; i< v2.size(); i++)
        if (_op1->point_inside(bv1, *v2[i], dummy))
          verts.push_back(*v2[i]);

      // reset the intersection tolerance
      _op2->set_intersection_tolerance(_intersection_tolerance);
    }

    // transform all vertices using the current transform of the CSG
    for (unsigned i=0; i< verts.size(); i++)
      verts[i] = T->transform_point(verts[i]);

    // clear validation flags for both operands
    _op1->_invalidated = _op2->_invalidated = false;
  }

  // copy _vertices to the passed vector
  vertices.resize(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    vertices[i] = &verts[i];
}
*/
