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
#include <Moby/XMLTree.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/PlanePrimitive.h>
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/OBB.h>
#include <Moby/Constants.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/CP.h>
#include <Moby/GJK.h>
#include <Moby/QP.h>
#include <Moby/BoxPrimitive.h>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using boost::const_pointer_cast; 
using boost::dynamic_pointer_cast; 
using std::map;
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
BoxPrimitive::BoxPrimitive(const Pose3d& T) : PolyhedralPrimitive(T)
{
  _xlen = 1;
  _ylen = 1;
  _zlen = 1;
  _edge_sample_length = std::numeric_limits<double>::max();
  calc_mass_properties();
}  

/// Constructs a cube of the specified size transformed by the given matrix
BoxPrimitive::BoxPrimitive(double xlen, double ylen, double zlen, const Pose3d& T) : PolyhedralPrimitive(T)
{
  _xlen = xlen;
  _ylen = ylen;
  _zlen = zlen;
  _edge_sample_length = std::numeric_limits<double>::max();
  calc_mass_properties();
}

/// Sets up the facets for the box primitive
void BoxPrimitive::get_facets(shared_ptr<const Pose3d> P, MatrixNd& M, VectorNd& q) const
{
  const unsigned N_FACETS = 6;

  // verify that the primitive knows about this pose 
  assert(_poses.find(const_pointer_cast<Pose3d>(P)) != _poses.end());

  // setup the normals to the pose
  Vector3d normals[N_FACETS];
  normals[0] = Vector3d(+0,+1,+0,P);
  normals[1] = Vector3d(+0,-1,+0,P);
  normals[2] = Vector3d(+0,+0,+1,P);
  normals[3] = Vector3d(+0,+0,-1,P);
  normals[4] = Vector3d(+1,+0,+0,P);
  normals[5] = Vector3d(-1,+0,+0,P);

  // setup the points on each plane
  Point3d points[N_FACETS];
  points[0] = Point3d(0.0,+_ylen*0.5,0.0,P);
  points[1] = Point3d(0.0,-_ylen*0.5,0.0,P);
  points[2] = Point3d(0.0,0.0,+_zlen*0.5,P);
  points[3] = Point3d(0.0,0.0,-_zlen*0.5,P);
  points[4] = Point3d(+_xlen*0.5,0.0,0.0,P);
  points[5] = Point3d(-_xlen*0.5,0.0,0.0,P);

  // get the transform to the global frame
  Transform3d T = Pose3d::calc_relative_pose(P, GLOBAL);

  // transform all points and normals
  for (unsigned i=0; i< N_FACETS; i++)
  {
    normals[i] = T.transform_vector(normals[i]);
    points[i] = T.transform_point(points[i]);
  }

  // setup M
  M.resize(N_FACETS,3);
  for (unsigned i=0; i< N_FACETS; i++)
    M.set_row(i, normals[i]);

  // setup q
  q.resize(N_FACETS);
  for (unsigned i=0; i< N_FACETS; i++)
    q[i] = normals[i].dot(points[i]);
assert(q.norm_inf() < 1e5);
}

/// Computes the signed distance from the box to a primitive
double BoxPrimitive::calc_signed_dist(shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const
{
  // now try box/sphere
  shared_ptr<const SpherePrimitive> spherep = dynamic_pointer_cast<const SpherePrimitive>(p);
  if (spherep)
    return calc_signed_dist(spherep, pthis, pp);

  // now try box/plane
  shared_ptr<const PlanePrimitive> planep = dynamic_pointer_cast<const PlanePrimitive>(p);
  if (planep)
  {
    shared_ptr<const Primitive> bthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return planep->calc_signed_dist(bthis, pp, pthis);
  }

  // now try box/heightmap
  shared_ptr<const HeightmapPrimitive> hmp = dynamic_pointer_cast<const HeightmapPrimitive>(p);
  if (hmp)
  {
    shared_ptr<const Primitive> bthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return hmp->calc_signed_dist(bthis, pp, pthis);
  }

  // try box/box
  shared_ptr<const BoxPrimitive> bp = dynamic_pointer_cast<const BoxPrimitive>(p);
  if (bp)
  {
    shared_ptr<const BoxPrimitive> bthis = dynamic_pointer_cast<const BoxPrimitive>(shared_from_this());
    return CP::find_cpoint(bthis, bp, pthis.pose, pp.pose, pthis, pp); 
  }

  // if the primitive is convex, can use GJK
  if (p->is_convex())
  {
    shared_ptr<const Pose3d> Pbox = pthis.pose;
    shared_ptr<const Pose3d> Pgeneric = pp.pose;
    shared_ptr<const Primitive> bthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return GJK::do_gjk(bthis, p, Pbox, Pgeneric, pthis, pp);
  }

  // try box/(non-convex) trimesh
  shared_ptr<const TriangleMeshPrimitive> trip = dynamic_pointer_cast<const TriangleMeshPrimitive>(p);
  if (trip)
  {
    shared_ptr<const Primitive> bthis = dynamic_pointer_cast<const Primitive>(shared_from_this());
    return trip->calc_signed_dist(bthis, pp, pthis);
  }
 
  // should never get here...
  assert(false); 
  return 0.0;
}

/// Finds closest point between a box and a sphere; returns the closest point on/in the box to the center of the sphere
double BoxPrimitive::calc_closest_points(shared_ptr<const SpherePrimitive> s, Point3d& pbox, Point3d& psph) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  Origin3d l, u, c, p;
  Matrix3d G;
  static QP qp;

  // to determine the closest point on/inside the box to the sphere, we
  // 1. compute the sphere center in the box frame
  // 2. find the closest point in the box frame, subject to the length
  //    constraints, to the sphere center
  // This is the following QP:
  // minimize 1/2*||p - c|| w.r.t. p
  // subject to lx < px < ux
  //            ly < py < uy
  //            lz < pz < uz

  // get the sphere center in the box frame
  Point3d sph_c(0.0, 0.0, 0.0, psph.pose);
  Point3d sph_c_A = Pose3d::transform_point(pbox.pose, sph_c);

  // setup the quadratic cost (identity matrix)
  G.set_identity();

  // setup the linear cost (-c, where c is sphere center)
  c[X] = -sph_c_A[X];
  c[Y] = -sph_c_A[Y];
  c[Z] = -sph_c_A[Z];

  // setup the box constraints
  const double HALF_X = _xlen * 0.5;
  const double HALF_Y = _ylen * 0.5;
  const double HALF_Z = _zlen * 0.5;
  l[X] = -HALF_X;  u[X] = HALF_X; 
  l[Y] = -HALF_Y;  u[Y] = HALF_Y; 
  l[Z] = -HALF_Z;  u[Z] = HALF_Z; 
  
  // solve the QP for the point nearest to the sphere center
  qp.qp_gradproj(G, c, l, u, 100, p, NEAR_ZERO);

  // setup the closest point on/in the box 
  pbox[X] = p[X];
  pbox[Y] = p[Y];
  pbox[Z] = p[Z];

  // setup the closest point on the sphere
  psph = Pose3d::transform_point(psph.pose, pbox);

  // get the sphere radius
  const double R = s->get_radius();

  // get the squared distance
  double sq_dist = sqr(p[X]-c[X]) + sqr(p[Y]-c[Y]) + sqr(p[Z]-c[Z]) - R*R;
  if (sq_dist > 0.0)
  {
    // sphere and box are not interpenetrating; find closest point on sphere
    double psph_nrm = psph.norm();
    assert(psph_nrm > NEAR_ZERO);
    psph /= psph_nrm;
    psph *= R;
    return std::sqrt(sq_dist);
  }
  else
    return -std::sqrt(-sq_dist);
}

/// Gets the distance of this box from a sphere
double BoxPrimitive::calc_signed_dist(shared_ptr<const SpherePrimitive> s, Point3d& pbox, Point3d& psph) const
{
  // get the sphere center in this frame 
  Point3d sph_c(0.0, 0.0, 0.0, psph.pose);
  Point3d sph_c_A = Pose3d::transform_point(pbox.pose, sph_c);

  // get the closest point
  double dist = calc_closest_point(sph_c_A, pbox) - s->get_radius();

  // get the vector from the center of the sphere to the box
  Vector3d v = Pose3d::transform_point(psph.pose, pbox);
  double vnorm = v.norm();
  if (vnorm == 0.0)
    psph = sph_c;
  else
    // compute point on sphere 
    psph = sph_c + v*((s->get_radius() + std::min(dist,0.0))/vnorm);

  return dist;
}

/// Sets the size of this box
/**
 * \note forces recomputation of the mesh
 */
void BoxPrimitive::set_size(double xlen, double ylen, double zlen)
{
  const unsigned X = 0, Y = 1, Z = 2;

  _xlen = xlen;
  _ylen = ylen;
  _zlen = zlen;
  if (_xlen < 0.0 && _ylen < 0.0 && _zlen < 0.0)
    throw std::runtime_error("Attempting to pass negative box dimensions to BoxPrimitive::set_size()");

  // recalculate the mass properties
  calc_mass_properties();

  // set lengths on each OBB
  for (map<CollisionGeometryPtr, OBBPtr>::iterator i = _obbs.begin(); i != _obbs.end(); i++)
  {
    // setup OBB half-lengths
    i->second->l[X] = _xlen * (double) 0.5;
    i->second->l[Y] = _ylen * (double) 0.5;
    i->second->l[Z] = _zlen * (double) 0.5;
  }

  // need to update visualization
  update_visualization();
}

/// Sets the edge sample length for this box
void BoxPrimitive::set_edge_sample_length(double len)
{
  _edge_sample_length = len;
}

/// Transforms the primitive
void BoxPrimitive::set_pose(const Pose3d& p)
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

/// Gets the set of vertices for the BoxPrimitive
void BoxPrimitive::get_vertices(shared_ptr<const Pose3d> P, vector<Point3d>& verts) const 
{
  // clear the set of vertices
  verts.clear();

  // verify that the box is setup correctly
  if (_xlen <= 0.0 || _ylen <= 0.0 || _zlen <= 0.0)
    return;

  // verify that the pose is found
  assert(_poses.find(const_pointer_cast<Pose3d>(P)) != _poses.end());

  // setup half-lengths
  const double XLEN = _xlen*(double) 0.5;
  const double YLEN = _ylen*(double) 0.5;
  const double ZLEN = _zlen*(double) 0.5;

  // add the vertices 
  verts.push_back(Point3d(XLEN,YLEN,ZLEN,P));
  verts.push_back(Point3d(XLEN,YLEN,-ZLEN,P));
  verts.push_back(Point3d(XLEN,-YLEN,ZLEN,P));
  verts.push_back(Point3d(XLEN,-YLEN,-ZLEN,P));
  verts.push_back(Point3d(-XLEN,YLEN,ZLEN,P));
  verts.push_back(Point3d(-XLEN,YLEN,-ZLEN,P));
  verts.push_back(Point3d(-XLEN,-YLEN,ZLEN,P));
  verts.push_back(Point3d(-XLEN,-YLEN,-ZLEN,P));

  // *************************************************************************
  // NOTE: following code will not be necessary when edge sample length goes
  // away
  // *************************************************************************

  // don't go any further if the edge sample length is infinity
  if (_edge_sample_length == std::numeric_limits<double>::max())
    return;

  const double ESL_SQ = _edge_sample_length*_edge_sample_length;
  vector<LineSeg3> edges;
  edges.push_back(LineSeg3(verts[0], verts[1]));
  edges.push_back(LineSeg3(verts[0], verts[2]));
  edges.push_back(LineSeg3(verts[0], verts[4]));
  edges.push_back(LineSeg3(verts[1], verts[3]));
  edges.push_back(LineSeg3(verts[1], verts[5]));
  edges.push_back(LineSeg3(verts[2], verts[3]));
  edges.push_back(LineSeg3(verts[2], verts[6]));
  edges.push_back(LineSeg3(verts[3], verts[7]));
  edges.push_back(LineSeg3(verts[4], verts[5]));
  edges.push_back(LineSeg3(verts[4], verts[6]));
  edges.push_back(LineSeg3(verts[5], verts[7]));
  edges.push_back(LineSeg3(verts[6], verts[7]));
  for (unsigned i=0; i< edges.size(); i++)
  {
    if ((edges[i].first - edges[i].second).norm_sq() >= ESL_SQ)
    {
      Point3d midpoint = edges[i].first*0.5 + edges[i].second*0.5;
      verts.push_back(midpoint);
      edges.push_back(LineSeg3(midpoint, edges[i].second));
      edges[i].second = midpoint;
    }
  }
}

/// Gets the set of vertices for the BoxPrimitive (constructing, if necessary)
/*
void BoxPrimitive::get_vertices(CollisionGeometryPtr geom, vector<const Point3d*>& vertices)
{
  // get the vertices for the geometry
  vector<Point3d>& verts = _vertices[geom];

  // see whether the vertices have been defined
  if (verts.empty())
  {
    if (_xlen <= 0.0 || _ylen <= 0.0 || _zlen <= 0.0)
    {
      vertices.clear(); 
      return;
    }

    // get the geometry pose
    shared_ptr<const Pose3d> gpose = geom->get_pose();

    // verify that this pose is defined w.r.t. the global frame
    shared_ptr<const Pose3d> P = get_pose();
    assert(!P->rpose);

    // setup transform
    Transform3d T;
    T.source = gpose;
    T.target = gpose;
    T.x = P->x;
    T.q = P->q;

    // setup half-lengths
    const double XLEN = _xlen*(double) 0.5;
    const double YLEN = _ylen*(double) 0.5;
    const double ZLEN = _zlen*(double) 0.5;

    // add the vertices 
    verts.push_back(T.transform_point(Point3d(XLEN,YLEN,ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(XLEN,YLEN,-ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(XLEN,-YLEN,ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(XLEN,-YLEN,-ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(-XLEN,YLEN,ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(-XLEN,YLEN,-ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(-XLEN,-YLEN,ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(-XLEN,-YLEN,-ZLEN,gpose)));
    
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
      const Point3d& v1 = verts[sp.first];
      const Point3d& v2 = verts[sp.second];

      // if edge is sufficiently small, don't need to do any more processing
      if ((v1-v2).norm() <= _edge_sample_length)
        continue;

      // subdivide the edge by first creating a new vertex
      unsigned vidx = verts.size();
      verts.push_back((v1+v2) * (double) 0.5);
      q.push(make_sorted_pair(sp.first, vidx)); 
      q.push(make_sorted_pair(sp.second, vidx)); 
    }
  }

  // copy the addresses of the computed vertices into 'vertices' 
  vertices.resize(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    vertices[i] = &verts[i];
}

/// Gets the set of vertices for the BoxPrimitive (constructing, if necessary)
void BoxPrimitive::get_vertices(CollisionGeometryPtr geom, vector<const Point3d*>& vertices)
{
  // get the vertices for the geometry
  vector<Point3d>& verts = _vertices[geom];

  // see whether the vertices have been defined
  if (verts.empty())
  {
    if (_xlen <= 0.0 || _ylen <= 0.0 || _zlen <= 0.0)
    {
      vertices.clear(); 
      return;
    }

    // get the geometry pose
    shared_ptr<const Pose3d> gpose = geom->get_pose();

    // verify that this pose is defined w.r.t. the global frame
    shared_ptr<const Pose3d> P = get_pose();
    assert(!P->rpose);

    // setup transform
    Transform3d T;
    T.source = gpose;
    T.target = gpose;
    T.x = P->x;
    T.q = P->q;

    // setup half-lengths
    const double XLEN = _xlen*(double) 0.5;
    const double YLEN = _ylen*(double) 0.5;
    const double ZLEN = _zlen*(double) 0.5;

    // add the vertices 
    verts.push_back(T.transform_point(Point3d(XLEN,YLEN,ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(XLEN,YLEN,-ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(XLEN,-YLEN,ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(XLEN,-YLEN,-ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(-XLEN,YLEN,ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(-XLEN,YLEN,-ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(-XLEN,-YLEN,ZLEN,gpose)));
    verts.push_back(T.transform_point(Point3d(-XLEN,-YLEN,-ZLEN,gpose)));
    
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
      const Point3d& v1 = verts[sp.first];
      const Point3d& v2 = verts[sp.second];

      // if edge is sufficiently small, don't need to do any more processing
      if ((v1-v2).norm() <= _edge_sample_length)
        continue;

      // subdivide the edge by first creating a new vertex
      unsigned vidx = verts.size();
      verts.push_back((v1+v2) * (double) 0.5);
      q.push(make_sorted_pair(sp.first, vidx)); 
      q.push(make_sorted_pair(sp.second, vidx)); 
    }
  }

  // copy the addresses of the computed vertices into 'vertices' 
  vertices.resize(verts.size());
  for (unsigned i=0; i< verts.size(); i++)
    vertices[i] = &verts[i];
}
*/

/// Recomputes the triangle mesh for the BoxPrimitive
shared_ptr<const IndexedTriArray> BoxPrimitive::get_mesh(shared_ptr<const Pose3d> P)
{
  // create the mesh
  const unsigned BOX_VERTS = 8, BOX_FACETS = 12;

  // setup lengths
  double XLEN = _xlen * (double) 0.5;
  double YLEN = _ylen * (double) 0.5;
  double ZLEN = _zlen * (double) 0.5;

  // need eight vertices
  Point3d verts[BOX_VERTS];
  verts[0] = Point3d(XLEN,YLEN,ZLEN,P);
  verts[1] = Point3d(XLEN,YLEN,-ZLEN,P);
  verts[2] = Point3d(XLEN,-YLEN,ZLEN,P);
  verts[3] = Point3d(XLEN,-YLEN,-ZLEN,P);
  verts[4] = Point3d(-XLEN,YLEN,ZLEN,P);
  verts[5] = Point3d(-XLEN,YLEN,-ZLEN,P);
  verts[6] = Point3d(-XLEN,-YLEN,ZLEN,P);
  verts[7] = Point3d(-XLEN,-YLEN,-ZLEN,P);

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
  return shared_ptr<IndexedTriArray>(new IndexedTriArray(verts, verts+BOX_VERTS, facets, facets+BOX_FACETS));
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
BVPtr BoxPrimitive::get_BVH_root(CollisionGeometryPtr geom)
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

    // get the pose for the primitive and geometry
    shared_ptr<const Pose3d> P = _cg_poses[geom]; 

    // setup the obb center and orientation
    obb->center = Point3d(0.0, 0.0, 0.0, P);
    obb->R.set_identity();

    // setup OBB half-lengths
    obb->l[X] = _xlen * (double) 0.5;
    obb->l[Y] = _ylen * (double) 0.5;
    obb->l[Z] = _zlen * (double) 0.5;
  }

  return obb;
}

/// Computes the signed distance to a point
double BoxPrimitive::calc_signed_dist(const Point3d& p) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // verify that the primitive knows about this pose 
  assert(_poses.find(const_pointer_cast<Pose3d>(p.pose)) != _poses.end()); 

  // setup extents
  double extents[3] = { _xlen*0.5, _ylen*0.5, _zlen*0.5 };

  // compute the squared distance to the p on the box
  bool inside = true;
  double sqrDist = 0.0;
  double intDist = -std::numeric_limits<double>::max();
  double delta = 0.0;
  for (unsigned i=0; i< 3; i++)
  {
    // see whether this dimension of the p lies below the negative extent
    if (p[i] < -extents[i])
    {
      delta = p[i] + extents[i];
      sqrDist += delta*delta;
      inside = false;
    }
    // see whether this dimension of the p lies above the positive extent
    else if (p[i] > extents[i])
    {
      delta = p[i] - extents[i];
      sqrDist += delta*delta;
      inside = false;
    }
    else if (inside)
    {
      double dist = std::min(std::fabs(p[i] - extents[i]), std::fabs(p[i] + extents[i]));
      intDist = std::max(intDist, dist);
    }
  }

  return (inside) ? -intDist : std::sqrt(sqrDist);
}

/// Computes the closest point on the box to a point (and returns the distance) 
double BoxPrimitive::calc_closest_point(const Point3d& point, Point3d& closest) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // verify that the primitive knows about this pose 
  assert(_poses.find(const_pointer_cast<Pose3d>(point.pose)) != _poses.end()); 

  // setup extents
  double extents[3] = { _xlen*0.5, _ylen*0.5, _zlen*0.5 };

  // setup the closest point
  closest = point;

  FILE_LOG(LOG_COLDET) << "BoxPrimitive::calc_closest_point() entered" << endl; 
  FILE_LOG(LOG_COLDET) << "  -- querying point " << point << " and box: " << extents[0] << ", " << extents[1] << ", " << extents[2] << endl;

  // compute the squared distance to the point on the box
  bool inside = true;
  double sqrDist = 0.0;
  double intDist = -std::numeric_limits<double>::max();
  double delta;
  for (unsigned i=0; i< 3; i++)
  {
    // see whether this dimension of the point lies below the negative extent
    if (point[i] < -extents[i])
    {
      delta = point[i] + extents[i];
      closest[i] = -extents[i];
      sqrDist += delta*delta;
      inside = false;
    }
    // see whether this dimension of the point lies above the positive extent
    else if (point[i] > extents[i])
    {
      delta = point[i] - extents[i];
      closest[i] = extents[i];
      sqrDist += delta*delta;
      inside = false;
    }
    else if (inside)
    {
      double dist = -std::min(std::fabs(extents[i] - point[i]), 
                              std::fabs(point[i] + extents[i]));
      intDist = std::max(intDist, dist);
    }
  }

  return (inside) ? intDist : std::sqrt(sqrDist);
}


/// Tests whether a point is inside or on the box
double BoxPrimitive::calc_dist_and_normal(const Point3d& point, std::vector<Vector3d>& normals) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // verify that the primitive knows about this pose 
  assert(_poses.find(const_pointer_cast<Pose3d>(point.pose)) != _poses.end()); 

  // clear the normals
  normals.clear();

  // setup extents
  double extents[3] = { _xlen*0.5, _ylen*0.5, _zlen*0.5 };

  FILE_LOG(LOG_COLDET) << "BoxPrimitive::calc_dist_and_normal() entered" << endl; 
  FILE_LOG(LOG_COLDET) << "  -- querying point " << point << " and box: " << extents[0] << ", " << extents[1] << ", " << extents[2] << endl;

  // compute the squared distance to the point on the box
  bool inside = true;
  double sqrDist = 0.0;
  double intDist = -std::numeric_limits<double>::max();
  double delta;
  for (unsigned i=0; i< 3; i++)
  {
    // see whether this dimension of the point lies below the negative extent
    if (point[i] < -extents[i])
    {
      delta = point[i] + extents[i];
      sqrDist += delta*delta;
      inside = false;
    }
    // see whether this dimension of the point lies above the positive extent
    else if (point[i] > extents[i])
    {
      delta = point[i] - extents[i];
      sqrDist += delta*delta;
      inside = false;
    }
    else if (inside)
    {
      double dist = -std::min(std::fabs(point[i] - extents[i]), std::fabs(point[i] + extents[i]));
      intDist = std::max(intDist, dist);
    }
  }

  // two cases: point outside the box and point inside the box
  // case 1: point inside the box
  if (inside)
  {
    FILE_LOG(LOG_COLDET) << "  point is inside the box; interior dist = " << intDist << endl;

    // p is inside the box; determine the normal
    double absPX = extents[X] - std::fabs(point[X]);
    double absPY = extents[Y] - std::fabs(point[Y]);
    double absPZ = extents[Z] - std::fabs(point[Z]);
    if (absPZ < absPX - NEAR_ZERO && absPZ < absPY - NEAR_ZERO)
    {
      if (point[Z] < (double) 0.0)
        normals.push_back(Vector3d(0,0,-1,point.pose));
      else
        normals.push_back(Vector3d(0,0,1,point.pose));
    }
    else if (absPY < absPZ - NEAR_ZERO && absPY < absPX - NEAR_ZERO)
    {
      if (point[Y] < (double) 0.0)
        normals.push_back(Vector3d(0,-1,0,point.pose));
      else
        normals.push_back(Vector3d(0,1,0,point.pose));
    }
    else if (absPX < absPY - NEAR_ZERO && absPX < absPZ - NEAR_ZERO)
    {
      if (point[X] < (double) 0.0)
        normals.push_back(Vector3d(-1,0,0,point.pose));
      else
        normals.push_back(Vector3d(1,0,0,point.pose));
    }
    else
    {
      // degenerate normal; do nothing now
      return (inside) ? intDist : std::sqrt(sqrDist);

      // degenerate normal; check for relative equality 
      if (CompGeom::rel_equal(absPX, absPY) && CompGeom::rel_equal(absPX, absPZ))
      {
        if (point[Z] < (double) 0.0)
          normals.push_back(Vector3d(0,0,-1,point.pose));
        else
          normals.push_back(Vector3d(0,0,1,point.pose));
        if (point[Y] < (double) 0.0)
          normals.push_back(Vector3d(0,-1,0,point.pose));
        else
          normals.push_back(Vector3d(0,1,0,point.pose));
        if (point[X] < (double) 0.0)
          normals.push_back(Vector3d(-1,0,0,point.pose));
        else
          normals.push_back(Vector3d(1,0,0,point.pose));
      }
      else if (CompGeom::rel_equal(absPX, absPY))
      {
        if (point[Y] < (double) 0.0)
          normals.push_back(Vector3d(0,-1,0,point.pose));
        else
          normals.push_back(Vector3d(0,1,0,point.pose));
        if (point[X] < (double) 0.0)
          normals.push_back(Vector3d(-1,0,0,point.pose));
        else
          normals.push_back(Vector3d(1,0,0,point.pose));
      }
      else if (CompGeom::rel_equal(absPX, absPZ))
      {
        if (point[Z] < (double) 0.0)
          normals.push_back(Vector3d(0,0,-1,point.pose));
        else
          normals.push_back(Vector3d(0,0,1,point.pose));
        if (point[X] < (double) 0.0)
          normals.push_back(Vector3d(-1,0,0,point.pose));
        else
          normals.push_back(Vector3d(1,0,0,point.pose));
      }
      else
      {
        if (point[Z] < (double) 0.0)
          normals.push_back(Vector3d(0,0,-1,point.pose));
        else
          normals.push_back(Vector3d(0,0,1,point.pose));
        if (point[Y] < (double) 0.0)
          normals.push_back(Vector3d(0,-1,0,point.pose));
        else
          normals.push_back(Vector3d(0,1,0,point.pose));
      }
    }
  }
  else
  {
    FILE_LOG(LOG_COLDET) << "  point is outside the box, sqr dist: " << sqrDist << endl;

    // the contact will be normal to all extents that the contact point
    // lies relatively outside of
    Vector3d normal;
    normal.set_zero();
    normal.pose = point.pose;

    if (std::fabs(point[X]) > extents[X] || CompGeom::rel_equal(std::fabs(point[X]), extents[X]))
      normal[X] = (point[X] > 0.0) ? 1.0 : -1.0;
    if (std::fabs(point[Y]) > extents[Y] || CompGeom::rel_equal(std::fabs(point[Y]), extents[Y]))
      normal[Y] = (point[Y] > 0.0) ? 1.0 : -1.0;
    if (std::fabs(point[Z]) > extents[Z] || CompGeom::rel_equal(std::fabs(point[Z]), extents[Z]))
      normal[Z] = (point[Z] > 0.0) ? 1.0 : -1.0;

    // normalize the normal vector
    normal.normalize();
  }

  return (inside) ? intDist : std::sqrt(sqrDist);
}

/**
/// Computes the intersection of a line segment with the box
 * \note for line segments that are partially or fully inside the box, the
 *       method only returns intersection if the second endpoint of the segment
 *       is farther inside than the first
bool BoxPrimitive::intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Vector3d& normal) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  double tmin = (double) 0.0;
  double tmax = (double) 2.0;
  static shared_ptr<Pose3d> P;

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

  // form a transform between segment space and box space 
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
  if (!determine_normal_abs(l, isect, P, normal))
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
*/


