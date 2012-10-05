/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <fstream>
#include <Moby/CollisionGeometry.h>
#include <Moby/XMLTree.h>
#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/VectorN.h>
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/DeformableBody.h>

using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::vector;
using std::map;
using std::list;
using std::endl;
using std::queue;
using std::cerr;
using std::string;

DeformableBody::DeformableBody()
{
  // ensure that config updates are enabled
  _config_updates_enabled = true;
}

/// Disables configuration updates
void DeformableBody::disable_config_updates()
{
  _config_updates_enabled = false;
}

/// Enables configuration updates
void DeformableBody::enable_config_updates()
{
  // re-enable position updates
  _config_updates_enabled = true;

  // update configuration-related variables
  calc_com_and_vels();
  update_geometries();
}  

/// Applies a generalized impulse to the deformable body
void DeformableBody::apply_generalized_impulse(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gj)
{
  const unsigned THREE_D = 3;

  assert(num_generalized_coordinates(gctype) == gj.size());
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    Vector3 j(gj[i*THREE_D], gj[i*THREE_D+1], gj[i*THREE_D+2]);
    _nodes[i]->xd += j/_nodes[i]->mass;
  }
}

/// Adds a generalized force to the deformable body
void DeformableBody::add_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gf)
{
  const unsigned THREE_D = 3;

  assert(num_generalized_coordinates(gctype) == gf.size());
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    Vector3 f(gf[i*THREE_D], gf[i*THREE_D+1], gf[i*THREE_D+2]);
    _nodes[i]->f += f;
  }
}

/// Gets the number of coordinates in the deformable body
unsigned DeformableBody::num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const
{
  const unsigned THREE_D = 3;
  return THREE_D * _nodes.size();
}

/// Calculates the position of the center-of-mass of the deformable body and the body velocities
/**
 * See wikipedia entry on angular momentum for more information.
 */
void DeformableBody::calc_com_and_vels()
{
  const unsigned X = 0, Y = 1, Z = 2;

  // init the center of mass and mass
  _x = ZEROS_3;
  _mass = (Real) 0.0;

  // determine the position of the com
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    _x += _nodes[i]->x * _nodes[i]->mass;
    _mass += _nodes[i]->mass;
  }
  _x /= _mass;
  FILE_LOG(LOG_DEFORM) << "new center of mass: " << _x << endl;

  // compute the linear velocity
  _xd = calc_point_vel(_x);

  // now determine the moment of inertia matrix
  _J = ZEROS_3x3;
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    Vector3 relpoint = _nodes[i]->x - _x;
    Real xsq = relpoint[X]*relpoint[X];
    Real ysq = relpoint[Y]*relpoint[Y];
    Real zsq = relpoint[Z]*relpoint[Z];
    _J(X,X) += _nodes[i]->mass * (ysq + zsq);
    _J(Y,Y) += _nodes[i]->mass * (xsq + zsq);
    _J(Z,Z) += _nodes[i]->mass * (xsq + ysq);
    _J(X,Y) -= _nodes[i]->mass * relpoint[X] * relpoint[Y];
    _J(X,Z) -= _nodes[i]->mass * relpoint[X] * relpoint[Z];
    _J(Y,Z) -= _nodes[i]->mass * relpoint[Y] * relpoint[Z];
  }

  // make J symmetric
  _J(Y,X) = _J(X,Y);
  _J(Z,X) = _J(X,Z);
  _J(Z,Y) = _J(Y,Z);

  // invert J
  _Jinv = Matrix3::inverse(_J);

  // now determine the angular momentum of the body
  Vector3 P = Vector3::cross(_x, _xd * _mass);
  for (unsigned i=0; i< _nodes.size(); i++)
    P += Vector3::cross(_nodes[i]->x - _x, _nodes[i]->xd * _nodes[i]->mass);

  // set the angular velocity
  _omega = _Jinv * P;
}

/// Gets the generalized coordinates of the deformable body
VectorN& DeformableBody::get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gc)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  // resize the gc vector (if necessary)
  gc.resize(THREE_D * _nodes.size());

  // populate the gc vector
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    const Vector3& x = _nodes[i]->x;
    gc[j++] = x[X];
    gc[j++] = x[Y];
    gc[j++] = x[Z];
  }

  return gc;
}

/// Gets the generalized velocity of the deformable body
VectorN& DeformableBody::get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, VectorN& gv)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  // resize the gv vector (if necessary)
  gv.resize(THREE_D * _nodes.size());

  // populate the gc vector
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    const Vector3& xd = _nodes[i]->xd;
    gv[j++] = xd[X];
    gv[j++] = xd[Y];
    gv[j++] = xd[Z];
  }

  return gv;
}

/// Gets the generalized velocity of the deformable body
VectorN& DeformableBody::get_generalized_acceleration(DynamicBody::GeneralizedCoordinateType gctype, VectorN& ga)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  // resize the ga vector (if necessary)
  ga.resize(THREE_D * _nodes.size());

  // populate the ga vector
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    // get the node
    shared_ptr<Node> node = dynamic_pointer_cast<Node>(_nodes[i]);

    // a = F / m 
    assert(!std::isnan(_nodes[i]->f[X]));
    assert(!std::isnan(_nodes[i]->f[Y]));
    assert(!std::isnan(_nodes[i]->f[Z]));
    ga[j++] = _nodes[i]->f[X] / _nodes[i]->mass;
    ga[j++] = _nodes[i]->f[Y] / _nodes[i]->mass;
    ga[j++] = _nodes[i]->f[Z] / _nodes[i]->mass;
  }

  return ga;
}

/// Sets the generalized coordinates for the body
void DeformableBody::set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gc)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;

  // verify gc vector is proper size
  assert(gc.size() == THREE_D * _nodes.size());

  // populate the state
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    Vector3& x = _nodes[i]->x;
    x[X] = gc[j++];
    x[Y] = gc[j++];
    x[Z] = gc[j++];
  }

  // see whether to update configuration-based variables
  if (_config_updates_enabled)
  {
    // calculate the position of the center-of-mass
    calc_com_and_vels();

    // update the geometries
    update_geometries();
  }
}

/// Sets the generalized velocity for the body
void DeformableBody::set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const VectorN& gv)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  // verify gv vector is proper size
  assert(gv.size() == THREE_D * _nodes.size());

  // populate the state
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    Vector3& xd = _nodes[i]->xd;
    xd[X] = gv[j++];
    xd[Y] = gv[j++];
    xd[Z] = gv[j++];
  }

  // see whether to update configuration-based variables
  if (_config_updates_enabled)
    calc_com_and_vels();
}

/// Gets the generalized inertia for the body
MatrixN& DeformableBody::get_generalized_inertia(DynamicBody::GeneralizedCoordinateType gctype, MatrixN& M)
{
  const unsigned THREE_D = 3;

  // set the inertia matrix to the correct size
  M.set_zero(_nodes.size() * THREE_D, _nodes.size() * THREE_D);

  // set the inertia matrix
  for (unsigned i=0, j=0; i< _nodes.size(); i++, j+= THREE_D)
    M(j, j) = M(j+1,j+1) = M(j+2, j+2) = _nodes[i]->mass;

  return M;
}

/// Gets the generalized forces for the body
VectorN& DeformableBody::get_generalized_forces(DynamicBody::GeneralizedCoordinateType gctype, VectorN& f)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  // resize f (if necessary)
  f.resize(THREE_D * _nodes.size());

  // populate f
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    const Vector3& fext = _nodes[i]->f;
    f[j++] = fext[X];
    f[j++] = fext[Y];
    f[j++] = fext[Z];
  }

  return f;
}

/// Converts a force and torque to a generalized force on the body
VectorN& DeformableBody::convert_to_generalized_force(DynamicBody::GeneralizedCoordinateType gctype, SingleBodyPtr body, const Vector3& f, const Vector3& t, VectorN& gf)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  gf.resize(THREE_D);
  gf[X] = f[X];
  gf[Y] = f[Y];
  gf[Z] = f[Z];

  return gf;
}

/// Resets the force / torque accumulators on the deformable body
void DeformableBody::reset_accumulators()
{
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    // reset the force on the node
    _nodes[i]->f = ZEROS_3;
  }
}

/// Transforms the dynamic body
void DeformableBody::transform(const Matrix4& T)
{
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    // transform the node position
    _nodes[i]->x = T.mult_point(_nodes[i]->x);
  }

  // calculate the position of the center-of-mass
  calc_com_and_vels();

  // update the geometries too
  update_geometries();
}

/// Calculates the kinetic energy of the deformable body
Real DeformableBody::calc_kinetic_energy() const
{
  Real KE = (Real) 0.0;
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    // update KE with linear component of energy
    KE += _nodes[i]->xd.norm_sq() * _nodes[i]->mass;
  }

  // don't forget to halve the K.E.
  KE *= 0.5;

  return KE;
}

/// Gets the deformed tetrahedron
Tetrahedron DeformableBody::get_tetrahedron(unsigned i) const
{
  return Tetrahedron(_nodes[_tetrahedra[i].a]->x, _nodes[_tetrahedra[i].b]->x,
                     _nodes[_tetrahedra[i].c]->x, _nodes[_tetrahedra[i].d]->x);
}

/// Finds the closest tetrahedron to a point; if multiple closest are found, returns one arbitrarily
unsigned DeformableBody::find_closest_tetrahedron(const Vector3& p) const
{
  FILE_LOG(LOG_DEFORM) << "DeformableBody::find_closest_tetrahedron() entered" << endl;

  // get the vector of tetrahedra
  assert(_tetrahedra.size() > 0);

  // closest tetrahedron is none initially
  unsigned closest = std::numeric_limits<unsigned>::max();

  // determine all tetrahedra that p could lie within
  list<unsigned> tet_ids;
  get_tetrahedra(p, std::back_inserter(tet_ids));
  if (tet_ids.empty())
    FILE_LOG(LOG_DEFORM) << " -- warning: point " << p << " not inside any tetrahedra!" << std::endl;

  // if list is empty, determine the closest bounding boxes 
  if (tet_ids.empty())
    for (unsigned i=0; i< _tetrahedra.size(); i++)
      tet_ids.push_back(i);

  FILE_LOG(LOG_DEFORM) << " -- examining " << tet_ids.size() << " / " << _tetra_mesh->num_tetra() << " tetrahedra" << std::endl;

  // setup the minimum distance
  Real min_dist = std::numeric_limits<Real>::max();

  // iterate through the tetrahedra
  BOOST_FOREACH(unsigned id, tet_ids)
  {
    // get the tetrahedron
    Tetrahedron tet(_nodes[_tetrahedra[id].a]->x, _nodes[_tetrahedra[id].b]->x,
                    _nodes[_tetrahedra[id].c]->x, _nodes[_tetrahedra[id].d]->x);

    FILE_LOG(LOG_DEFORM) << "distance of " << p << " to tetrahedron " << id << std::endl;
    FILE_LOG(LOG_DEFORM) << "  " << _nodes[_tetrahedra[id].a]->x << std::endl;
    FILE_LOG(LOG_DEFORM) << "  " << _nodes[_tetrahedra[id].b]->x << std::endl;
    FILE_LOG(LOG_DEFORM) << "  " << _nodes[_tetrahedra[id].c]->x << std::endl;
    FILE_LOG(LOG_DEFORM) << "  " << _nodes[_tetrahedra[id].d]->x << std::endl;

    // determine the distance from the point to the tetrahedron
    Real dist = tet.calc_signed_dist(p);
    FILE_LOG(LOG_DEFORM) << "is " << dist << std::endl;
    if (dist > min_dist)
      continue;

    // store the minimum distance
    min_dist = dist;

    // store the closest tetrahedron
    closest = id;
  }

  FILE_LOG(LOG_DEFORM) << " closest tetrahedron is " << closest << " with distance " << min_dist << endl;
  FILE_LOG(LOG_DEFORM) << "DeformableBody::find_closest_tetrahedron() exited" << endl;

  return closest;
}

/// Outputs the deformable body to VRML
std::ostream& DeformableBody::to_vrml(std::ostream& o, Real scale) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup random colors for each tetrahedron
  std::vector<Vector3> colors(_tetrahedra.size());
  for (unsigned i=0; i< colors.size(); i++)
  {
    colors[i][0] = (Real) rand() / RAND_MAX;
    colors[i][1] = (Real) rand() / RAND_MAX;
    colors[i][2] = (Real) rand() / RAND_MAX;
  }

  // draw each tetrahedron
  for (unsigned i=0; i< _tetrahedra.size(); i++)
  {
    o << "Shape {" << endl;
    o << "  appearance Appearance { material Material { diffuseColor ";
    o << colors[i][0] << " " << colors[i][1] << " " << colors[i][2] << " } }";
    o << endl;
    o << "  geometry IndexedLineSet {" << endl;
    o << "  coord Coordinate { point [ ";
    o << _nodes[_tetrahedra[i].a]->x[X] << " " << _nodes[_tetrahedra[i].a]->x[Y] << " " << _nodes[_tetrahedra[i].a]->x[Z] << ", ";
    o << _nodes[_tetrahedra[i].b]->x[X] << " " << _nodes[_tetrahedra[i].b]->x[Y] << " " << _nodes[_tetrahedra[i].b]->x[Z] << ", ";
    o << _nodes[_tetrahedra[i].c]->x[X] << " " << _nodes[_tetrahedra[i].c]->x[Y] << " " << _nodes[_tetrahedra[i].c]->x[Z] << ", ";
    o << _nodes[_tetrahedra[i].d]->x[X] << " " << _nodes[_tetrahedra[i].d]->x[Y] << " " << _nodes[_tetrahedra[i].d]->x[Z] << " ]}" << endl;
    o << "  coordIndex [ 0 1 -1 0 2 -1 0 3 -1 1 2 -1 1 3 -1 2 3 -1 ] } }" << endl;
  }

  // get the triangle mesh
  shared_ptr<TriangleMeshPrimitive> trimeshp = dynamic_pointer_cast<TriangleMeshPrimitive>(_tri_mesh);
  shared_ptr<const IndexedTriArray> ita = trimeshp->get_mesh();

  // now draw each vertex in the original trimesh
  for (unsigned i=0; i< ita->get_vertices().size(); i++)
  {
    const Vector3& v = ita->get_vertices()[i];
    o << "Transform {" << endl;
    o << "  translation " << v[X] << " " << v[Y] << " " << v[Z] << endl;
    o << "  children Shape {" << endl;
    o << "    appearance Appearance { material Material { diffuseColor ";
    o << colors[i][0] << " " << colors[i][1] << " " << colors[i][2] << " } }";
    o << endl;
    o << "    geometry Box { size " << scale << " " << scale << " " << scale << " } } }" << endl;   
  }

/*
  // now draw each vertex at its new position
  for (unsigned i=0; i< _vertex_map.size(); i++)
  {
    // get the vertex position
    const VertexMap& vmap = _vertex_map[i];
    const IndexedTetra& tetra = _tetrahedra[vmap.tetra];
    Tetrahedron t(_nodes[tetra.a]->x, _nodes[tetra.b]->x, _nodes[tetra.c]->x, _nodes[tetra.d]->x);
    Vector3 p = t.calc_point(vmap.uvw[0], vmap.uvw[1], vmap.uvw[2]);
    
    o << "Transform {" << endl;
    o << "  translation " << p[X] << " " << p[Y] << " " << p[Z] << endl;
    o << "  children Shape {" << endl;
    o << "    appearance Appearance { material Material { diffuseColor ";
    o << colors[i][0] << " " << colors[i][1] << " " << colors[i][2] << " } }";
    o << endl;
    o << "    geometry Sphere { radius " << scale << " } } }" << endl;   
  }
*/
}

/// Sets the mesh for the deformable body
/**
 * \note Also sets the nodes for the body (and thereby modifies the com 
 *       position and velocity.
 */
void DeformableBody::set_mesh(shared_ptr<const IndexedTetraArray> tetra_mesh, shared_ptr<Primitive> tri_mesh)
{
  // set the default mass
  const Real DEFAULT_MASS = (Real) 1.0;

  // store the meshes
  _tetra_mesh = tetra_mesh;
  _tri_mesh = tri_mesh;

  // store the tetrahedra
  _tetrahedra = tetra_mesh->get_tetra();

  // get the vertices of the meshes
  const vector<Vector3>& tet_vertices = tetra_mesh->get_vertices();
  const vector<Vector3>& tri_vertices = tri_mesh->get_mesh()->get_vertices();

  // setup the nodes in the body based on the tetrahedra vertices
  _nodes = vector<shared_ptr<Node> >(tet_vertices.size());
  for (unsigned i=0; i< tet_vertices.size(); i++)
  {
    _nodes[i] = shared_ptr<Node>(new Node);
    _nodes[i]->x = tet_vertices[i];
    _nodes[i]->xd = _nodes[i]->xdd = ZEROS_3;
    _nodes[i]->mass = DEFAULT_MASS;
  }

  // create an AABB hierarchy around the tetrahedral mesh
  _hroot = build_AABB_tree(_aabb_tetra_map);

  // update com and velocity
  calc_com_and_vels();

  // for each vertex in the triangle mesh, determine (a) _a_ tetrahedron that
  // the mesh belongs to (there may be multiple, if the vertex lies coincident
  // with a vertex or edge in the tetrahedral mesh) and (b) the barycentric
  // coordinates of that vertex in the tetrahedron
  _vertex_map = vector<VertexMap>(tri_vertices.size());
  for (unsigned i=0; i< tri_vertices.size(); i++)
  {
    // get the closest tetrahedron
    unsigned closest = find_closest_tetrahedron(tri_vertices[i]);

    // should be one that is closest
    assert(closest < std::numeric_limits<unsigned>::max());

    // mark the tetrahedron in the vertex map
    _vertex_map[i].tetra = closest; 

    // get the tetrahedron
    Tetrahedron tet = get_tetrahedron(closest);

    // determine the barycentric coordinates
    Real u, v, w;
    tet.determine_barycentric_coords(tri_vertices[i], u, v, w);

    // correct barycentric coordinates
    if (u < (Real) 0.0) u = (Real) 0.0;
    else if (u > (Real) 1.0) u = (Real) 1.0;
    if (v < (Real) 0.0) v = (Real) 0.0;
    else if (u + v > (Real) 1.0) v = (Real) 1.0 - u;
    if (w < (Real) 0.0) w = (Real) 0.0;
    else if (u + v + w > (Real) 1.0) w = (Real) 1.0 - u - v;

    // store barycentric coords
    _vertex_map[i].uvw[0] = u;
    _vertex_map[i].uvw[1] = v;
    _vertex_map[i].uvw[2] = w;
  }

  // create a triangle mesh primitive and update the geometries
  _cgeom_primitive = shared_ptr<TriangleMeshPrimitive>(new TriangleMeshPrimitive);
  _cgeom_primitive->set_deformable(true);
  set_visualization_data(_cgeom_primitive->get_visualization());
  if (!_geometry)
    _geometry = CollisionGeometryPtr(new CollisionGeometry);
  _geometry->set_single_body(get_this());
  _geometry->set_geometry(_cgeom_primitive);
  _geometry->set_transform(IDENTITY_4x4, false);
  update_geometries(); 
}

/// Updates collision and visualization geometries using present state of the nodes of the deformable body
void DeformableBody::update_geometries()
{
  // setup the vertices for the triangle mesh
  vector<Vector3> vertices(_tri_mesh->get_mesh()->get_vertices().size());
 
  // get the facets for the triangle mesh
  const vector<IndexedTri> facets = _tri_mesh->get_mesh()->get_facets();

  // sanity check
  assert(_vertex_map.size() == vertices.size());

  // for each vertex in the triangle mesh, get the index of the tetrahedron
  // and barycentric coordinates and determine the new position of the vertex
  for (unsigned i=0; i< vertices.size(); i++)
  {
    // get the vertex map
    VertexMap& vmap = _vertex_map[i];
       
    // get the indexed tetrahedron
    const IndexedTetra& tetra = _tetrahedra[vmap.tetra];

    // form the tetrahedron referenced by the nodes
    Tetrahedron t(_nodes[tetra.a]->x, _nodes[tetra.b]->x, _nodes[tetra.c]->x, _nodes[tetra.d]->x);

    // get the new position of the vertex
    vertices[i] = t.calc_point(vmap.uvw[0], vmap.uvw[1], vmap.uvw[2]);
  }

  // update the collision geometry
  _cgeom_primitive->set_mesh(shared_ptr<IndexedTriArray>(new IndexedTriArray(vertices.begin(), vertices.end(), facets.begin(), facets.end())));

  // rebuild the AABB hierarchy
  _hroot = build_AABB_tree(_aabb_tetra_map);
}

/// Calculates the velocity of a point on the deformable body
Vector3 DeformableBody::calc_point_vel(const Vector3& p) const
{
  // find the closest tetrahedron
  unsigned closest = find_closest_tetrahedron(p);

  // should be a closest 
  assert(closest != std::numeric_limits<unsigned>::max());

  // get the tetrahedron
  Tetrahedron tet = get_tetrahedron(closest);

  // determine the barycentric coordinates
  Real u, v, w;
  tet.determine_barycentric_coords(p, u, v, w);

  // get the velocities at the vertices
  const IndexedTetra& itet = _tetrahedra[closest];
  const Vector3& va = _nodes[itet.a]->xd;
  const Vector3& vb = _nodes[itet.b]->xd;
  const Vector3& vc = _nodes[itet.c]->xd;
  const Vector3& vd = _nodes[itet.d]->xd;

  // determine the velocity at p using the barycentric function
  return vb*u + vc*v + vd*w + va*((Real) 1.0 - u - v - w);
}

/// Adds a force to the body
void DeformableBody::add_force(const Vector3& f)
{
  Vector3 f_per_node = f/_nodes.size();
  for (unsigned i=0; i< _nodes.size(); i++)
    _nodes[i]->f += f_per_node;
}

/// Adds a force to the deformable body at a given point
void DeformableBody::add_force(const Vector3& f, const Vector3& p)
{
  // find the closest tetrahedron
  unsigned closest = find_closest_tetrahedron(p);

  // should be a closest 
  assert(closest != std::numeric_limits<unsigned>::max());

  // get the tetrahedron
  Tetrahedron tet = get_tetrahedron(closest);

  // determine the barycentric coordinates
  Real u, v, w;
  tet.determine_barycentric_coords(p, u, v, w);

  // apply the force using the barycentric coordinates
  const IndexedTetra& itet = _tetrahedra[closest];
  _nodes[itet.a]->f += f * ((Real) 1.0 - u - v - w);
  _nodes[itet.b]->f += f * u;
  _nodes[itet.c]->f += f * v;
  _nodes[itet.c]->f += f * w;
}

/// Implements Base::load_from_xml()
void DeformableBody::load_from_xml(XMLTreeConstPtr node, map<string, BasePtr>& id_map)
{
  map<std::string, BasePtr>::const_iterator id_iter;
  vector<VertexMap> vertex_map;
  shared_ptr<const IndexedTetraArray> tetra_mesh;
  shared_ptr<Primitive> tri_mesh;

  // load parent data
  SingleBody::load_from_xml(node, id_map);

  // ***********************************************************************
  // don't verify the node name -- this class has derived classes
  // ***********************************************************************

  // read nodes- note that this must be done before loading the meshes
  list<XMLTreeConstPtr> nodes_nodes = node->find_child_nodes("Node");
  if (!nodes_nodes.empty())
  {
    // setup a vector of nodes
    vector<shared_ptr<Node> > new_nodes;

    // read in the set of nodes
    for (list<XMLTreeConstPtr>::const_iterator i = nodes_nodes.begin(); i != nodes_nodes.end(); i++)
    {
      // create a new Node
      shared_ptr<Node> n(new Node);

      // read the node position
      const XMLAttrib* node_x_attr = (*i)->get_attrib("position");
      if (node_x_attr)
        node_x_attr->get_vector_value(n->x);

      // read the node velocity
      const XMLAttrib* node_xd_attr = (*i)->get_attrib("velocity");
      if (node_xd_attr)
        node_xd_attr->get_vector_value(n->xd);      

      // read the node acceleration
      const XMLAttrib* node_xdd_attr = (*i)->get_attrib("accel");
      if (node_xdd_attr)
        node_xdd_attr->get_vector_value(n->xdd);

      // read the force on the node
      const XMLAttrib* node_f_attr = (*i)->get_attrib("force");
      if (node_f_attr)
        node_f_attr->get_vector_value(n->f);

      // read the node mass
      const XMLAttrib* node_mass_attr = (*i)->get_attrib("mass");
      if (node_mass_attr)
        n->mass = node_mass_attr->get_real_value();

      // add the node
      new_nodes.push_back(n);
    }

    _nodes = new_nodes;
  }

  // load the geometry (if any)
  list<XMLTreeConstPtr> cg_nodes = node->find_child_nodes("CollisionGeometry");

  // DeformableBody currently only supports one geometry (this is for
  // practicality, not necessity...)
  if (cg_nodes.size() > 1)
    std::cerr << "DeformableBody::load_from_xml() - multiple CollisionGeometry nodes specified" << endl << "(only one such node is supported; reading only one..." << endl;

  // create the geometry and populate its data
  if (!cg_nodes.empty())
  {
    // create the geometry
    _geometry = CollisionGeometryPtr(new CollisionGeometry());

    // point it to the deformable body
    _geometry->set_single_body(get_this());

    // populate its data
    _geometry->load_from_xml(cg_nodes.front(), id_map);
  }

  // read the tetrahedral mesh primitive
  const XMLAttrib* tetmesh_id_attr = node->get_attrib("tetra-mesh-id");
  if (tetmesh_id_attr)
  {
    // get the tetrahedra mesh ID
    const string& id = tetmesh_id_attr->get_string_value();

    // search for it
    map<string, BasePtr>::const_iterator id_iter = id_map.find(id);
    if (id_iter == id_map.end())
    {
      cerr << "DeformableBody::load_from_xml() - tetra mesh with ID '";
      cerr << id << "'" << endl << "  not found in offending node: ";
      cerr << endl << *node;
    }
    else
    {
      tetra_mesh = dynamic_pointer_cast<IndexedTetraArray>(id_iter->second);
      if (!tetra_mesh)
      {
        cerr << "DeformableBody::load_from_xml() - ID '";
        cerr << id << "'" << endl << "  not an indexed tetrahedron array in ";
        cerr << "node: " << endl << *node;
      }
    }
  }

  // read the primitive that contains the triangle mesh
  const XMLAttrib* trimesh_prim_id_attr = node->get_attrib("tri-mesh-primitive-id");
  if (trimesh_prim_id_attr)
  {
    // get the triangle mesh primitive ID
    const string& id = trimesh_prim_id_attr->get_string_value();

    // search for it
    map<string, BasePtr>::const_iterator id_iter = id_map.find(id);
    if (id_iter == id_map.end())
    {
      cerr << "DeformableBody::load_from_xml() - primitive with ID '";
      cerr << id << "'" << endl << "  not found in offending node: ";
      cerr << endl << *node;
    }
    else
    {
      tri_mesh = dynamic_pointer_cast<Primitive>(id_iter->second);
      if (!tri_mesh)
      {
        cerr << "DeformableBody::load_from_xml() - primitive with ID '";
        cerr << id << "'" << endl << "  not a primitive after all in ";
        cerr << "node: " << endl << *node;
      }
    }
  }

  // load the vertex map, if given
  list<XMLTreeConstPtr> vmap_nodes = node->find_child_nodes("VertexMap");
  if (!vmap_nodes.empty())
  {
    assert(vmap_nodes.size() == 1);
    list<XMLTreeConstPtr> map_nodes = vmap_nodes.front()->find_child_nodes("Mapping");
    BOOST_FOREACH(XMLTreeConstPtr mapping_node, map_nodes)
    {
      const XMLAttrib* tetra_id_attr = node->get_attrib("tetra");
      const XMLAttrib* uvw_id_attr = node->get_attrib("uvw");
      if (tetra_id_attr && uvw_id_attr)
      {
        VertexMap vmap;
        Vector3 uvw;
        vmap.tetra = tetra_id_attr->get_unsigned_value();
        uvw_id_attr->get_vector_value(uvw);
        vmap,uvw[0] = uvw[0];
        vmap,uvw[1] = uvw[1];
        vmap,uvw[2] = uvw[2];
        vertex_map.push_back(vmap);
      }
    }
  }

  // only if both pointers are set do we set the mesh
  if (tetra_mesh && tri_mesh)
  {
    // set the default mass
    const Real DEFAULT_MASS = (Real) 1.0;

    // set the mesh
    set_mesh(tetra_mesh, tri_mesh);

/*
    // store the tetrahedron
    _tetrahedra = tetra_mesh->get_tetra();

    // store the meshes
    _tri_mesh = tri_mesh;
    _tetra_mesh = tetra_mesh;

    // get the vertices of the meshes
    const vector<Vector3>& tet_vertices = tetra_mesh->get_vertices();
    const vector<Vector3>& tri_vertices = tri_mesh->get_mesh()->get_vertices();

    // setup the nodes in the body based on the tetrahedra vertices
    // *if they're not already set properly*
    if (_nodes.size() != tet_vertices.size())
    {
      _nodes = vector<shared_ptr<Node> >(tet_vertices.size());
      for (unsigned i=0; i< tet_vertices.size(); i++)
      {
        _nodes[i] = shared_ptr<Node>(new Node);
        _nodes[i]->x = tet_vertices[i];
        _nodes[i]->xd = _nodes[i]->xdd = ZEROS_3;
        _nodes[i]->mass = DEFAULT_MASS;
      }
    }

    // create an AABB hierarchy around the tetrahedral mesh
    _hroot = build_AABB_tree(_aabb_tetra_map);

    // update com and velocity
    calc_com_and_vels();

    // setup the vertex map, if necessary
    if (vertex_map.size() == tri_vertices.size())
      _vertex_map = vertex_map;
    else
    {
      _vertex_map = vector<VertexMap>(tri_vertices.size());
      for (unsigned i=0; i< tri_vertices.size(); i++)
      {
        // get the closest tetrahedron
        unsigned closest = find_closest_tetrahedron(tri_vertices[i]);

        // should be one that is closest
        assert(closest < std::numeric_limits<unsigned>::max());

        // mark the tetrahedron in the vertex map
        _vertex_map[i].tetra = closest; 

        // get the tetrahedron
        Tetrahedron tet = get_tetrahedron(closest);

        // determine the barycentric coordinates
        Real u, v, w;
        tet.determine_barycentric_coords(tri_vertices[i], u, v, w);

        // correct barycentric coordinates
        if (u < (Real) 0.0) u = (Real) 0.0;
        else if (u > (Real) 1.0) u = (Real) 1.0;
        if (v < (Real) 0.0) v = (Real) 0.0;
        else if (u + v > (Real) 1.0) v = (Real) 1.0 - u;
        if (w < (Real) 0.0) w = (Real) 0.0;
        else if (u + v + w > (Real) 1.0) w = (Real) 1.0 - u - v;

        // store barycentric coords
        _vertex_map[i].uvw[0] = u;
        _vertex_map[i].uvw[1] = v;
        _vertex_map[i].uvw[2] = w;
      }
    }

    // create a triangle mesh primitive and update the geometries
    _cgeom_primitive = shared_ptr<TriangleMeshPrimitive>(new TriangleMeshPrimitive);
    _cgeom_primitive->set_deformable(true);
    set_visualization_data(_cgeom_primitive->get_visualization());
    if (!_geometry)
      _geometry = CollisionGeometryPtr(new CollisionGeometry);
    _geometry->set_single_body(get_this());
    _geometry->set_geometry(_cgeom_primitive);
    _geometry->set_transform(IDENTITY_4x4, false);
    update_geometries(); 
*/
  }

  // read a transform to be applied to the body, if provided
  const XMLAttrib* T_attr = node->get_attrib("transform");
  if (T_attr)
  {
    Matrix4 T;
    T_attr->get_matrix_value(T);
    if (!Matrix4::valid_transform(T))
    {
      cerr << "DeformableBody::load_from_xml() warning: invalid transform ";
      cerr << endl << T << " when reading node " << endl;
      cerr << *node << endl;
      cerr << "  --> possibly a floating-point error..." << endl;
    }
    transform(T);
  }
}

/// Implements Base::save_to_xml()
void DeformableBody::save_to_xml(XMLTreePtr node, list<BaseConstPtr>& shared_objects) const
{
  // save parent data
  SingleBody::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "DeformableBody";

  // save the meshes
  if (_tri_mesh)
  {
    node->attribs.insert(XMLAttrib("tri-mesh-primitive-id", _tri_mesh->id));
    shared_objects.push_back(_tri_mesh);
  }
  if (_tetra_mesh)
  {
    node->attribs.insert(XMLAttrib("tetra-mesh-id", _tetra_mesh->id));
    shared_objects.push_back(_tetra_mesh);
  }

  // save the collision geometry
  if (_geometry)
  {
    XMLTreePtr geom_node(new XMLTree("CollisionGeometry"));
    node->add_child(geom_node);
    _geometry->save_to_xml(geom_node, shared_objects);
  }

  // save the vertex map
  XMLTreePtr vertex_map_node(new XMLTree("VertexMap"));
  for (unsigned i=0; i< _vertex_map.size(); i++)
  {
    XMLTreePtr mapping_node(new XMLTree("Mapping"));
    mapping_node->attribs.insert(XMLAttrib("tetra", _vertex_map[i].tetra));
    Vector3 uvw(_vertex_map[i].uvw[0], _vertex_map[i].uvw[1], _vertex_map[i].uvw[2]);
    mapping_node->attribs.insert(XMLAttrib("uvw", uvw));
    vertex_map_node->add_child(mapping_node);
  }
  node->add_child(vertex_map_node);

  // save the nodes
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    XMLTreePtr node_node(new XMLTree("Node"));
    node_node->attribs.insert(XMLAttrib("position", _nodes[i]->x));
    node_node->attribs.insert(XMLAttrib("velocity", _nodes[i]->xd));
    node_node->attribs.insert(XMLAttrib("accel", _nodes[i]->xdd));
    node_node->attribs.insert(XMLAttrib("force", _nodes[i]->f));
    node_node->attribs.insert(XMLAttrib("mass", _nodes[i]->mass));
    node->add_child(node_node);
  }
}

// *************************************************************************
// Methods for creating AABB hierarchy around tetrahedral mesh begin
// *************************************************************************

/// Creates a AABB hierarchy around a tetrahedral mesh
AABBPtr DeformableBody::build_AABB_tree(map<BVPtr, list<unsigned> >& aabb_tetra_map)
{
  const unsigned THREE_D = 3;
  const Real EXP = 0.1;          // box expansion constant
  AABBPtr child1, child2;
  list<unsigned> ptetra, ntetra;
  Vector3 eps(EXP, EXP, EXP);

  FILE_LOG(LOG_BV) << "DeformableBody::build_AABB_tree() entered" << endl;

  // clear the map
  aabb_tetra_map.clear();

  // build a AABB around all nodes 
  AABBPtr root(new AABB);
  if (!_nodes.empty())
  {
    root->minp = _nodes.front()->x;
    root->maxp = _nodes.front()->x;
    for (unsigned i=1; i< _nodes.size(); i++)
      for (unsigned j=0; j< THREE_D; j++)
      {
        if (root->minp[j] > _nodes[i]->x[j])
          root->minp[j] = _nodes[i]->x[j];
        if (root->maxp[j] < _nodes[i]->x[j])
          root->maxp[j] = _nodes[i]->x[j];
      }
  }

  // set root to point to all tetrahedra indices
  list<unsigned> tetra_idx;
  for (unsigned i=0; i< _tetrahedra.size(); i++)
    tetra_idx.push_back(i);

  // setup mapping from AABB to mesh
  aabb_tetra_map[root] = tetra_idx;

  FILE_LOG(LOG_BV) << "  -- created root: " << root << endl;
  FILE_LOG(LOG_BV) << *root << endl;

  // add the root to a splitting queue *only if it has more than tetra*
  queue<AABBPtr> Q;
  if (tetra_idx.size() > 1)
    Q.push(root);

  // split until we can't split further
  while (!Q.empty())
  {
    // get the bounding box off of the top of the queue
    AABBPtr aabb = Q.front();
    Q.pop();

    // get the list of tetrahedra for this bounding box
    list<unsigned>& tlist = aabb_tetra_map[aabb];

    FILE_LOG(LOG_BV) << "  -- splitting AABB: " << aabb << endl;

    // split the bounding box across each of the three axes
    for (unsigned i=0; i< 3; i++)
      if (split(aabb, child1, child2, i, tlist, ptetra, ntetra))
        break;

    // make sure that this AABB was divisible
    if (!child1)
      continue;

    // store tetrahedra present in the children and clear tetrahedra present
    // in the parent
    aabb_tetra_map[aabb].clear();
    aabb_tetra_map[child1] = ptetra;
    aabb_tetra_map[child2] = ntetra;

    // setup child pointers
    aabb->children.push_back(child1);
    aabb->children.push_back(child2);

    // add children to the queue for processing if they have more than one tri
    if (ptetra.size() > 1)
      Q.push(child1);
    if (ntetra.size() > 1)
      Q.push(child2);

    FILE_LOG(LOG_BV) << "  -- built children " << child1 << ", " << child2 << " from " << aabb << endl;
    FILE_LOG(LOG_BV) << "  -- ID: " << child1 << endl << *child1;
    if (LOGGING(LOG_BV))
    {
      std::ostringstream oss;
      for (list<unsigned>::const_iterator i = ptetra.begin(); i != ptetra.end(); i++)
        oss << " " << *i;
      FILE_LOG(LOG_BV) << "    tetrahedra:" << oss.str() << std::endl;
    }
    FILE_LOG(LOG_BV) << "  -- ID: " << child2 << endl << *child2;
    if (LOGGING(LOG_BV))
    {
      std::ostringstream oss;
      for (list<unsigned>::const_iterator i = ntetra.begin(); i != ntetra.end(); i++)
        oss << " " << *i;
      FILE_LOG(LOG_BV) << "    tetrahedra:" << oss.str() << std::endl;
    }
  }

  // now, collapse the tree
  Q.push(root);
  while (!Q.empty())
  {
    // for any children with a greater volume than the aabb in question,
    // remove the grandchildren and add them as children
    AABBPtr aabb = Q.front();
    Real vol = aabb->calc_volume();
    bool erased_one = false;
    for (list<BVPtr>::iterator i = aabb->children.begin(); i != aabb->children.end(); )
    {
      // get the AABB
      AABBPtr child = dynamic_pointer_cast<AABB>(*i);

      // get the volume of this child
      Real voli = child->calc_volume();
      if (!(*i)->is_leaf() && voli > vol + NEAR_ZERO)
      {
        erased_one = true;
        BOOST_FOREACH(BVPtr gchild, (*i)->children)
          aabb->children.push_back(gchild);
        i = aabb->children.erase(i);
      }
      else
        i++;
    }

    if (!erased_one)
    {
      Q.pop();
      BOOST_FOREACH(BVPtr child, aabb->children)
        if (!child->is_leaf())
          Q.push(dynamic_pointer_cast<AABB>(child));
    }
  }

  // expand the boxes in the tree
  Q.push(root);
  while (!Q.empty())
  {
    AABBPtr aabb = Q.front();
    Q.pop();
    aabb->minp -= eps;
    aabb->maxp += eps;
    if (!aabb->is_leaf())
      BOOST_FOREACH(BVPtr child, aabb->children)
        Q.push(dynamic_pointer_cast<AABB>(child));
  }

  // wipe out userdata 
  Q.push(root);
  while (!Q.empty())
  {
    // get the element off of the front of the queue
    AABBPtr aabb = Q.front();
    Q.pop();

    // add all children to the queue
    if (!aabb->is_leaf())
      BOOST_FOREACH(BVPtr child, aabb->children)
        Q.push(dynamic_pointer_cast<AABB>(child));

    // wipe out userdata
    aabb->userdata = shared_ptr<void>();
  }

  if (LOGGING(LOG_BV))
    for (map<BVPtr, list<unsigned> >::const_iterator i = aabb_tetra_map.begin(); i != aabb_tetra_map.end(); i++)
    {
      std::ostringstream oss;
      oss << "Tetrahedra in BV " << i->first << ": ";
      for (list<unsigned>::const_iterator j = i->second.begin(); j != i->second.end(); j++)
        oss << *j << " ";
      FILE_LOG(LOG_BV) << oss.str() << endl;
    } 

  FILE_LOG(LOG_BV) << "DeformableBody::build_AABB_tree() exited" << endl;

  return root;
}

/// Splits a collection of tetrahedra along a splitting plane into 2 new meshes 
void DeformableBody::split_tetra(const Vector3& point, unsigned axis, const list<unsigned>& otetra, list<unsigned>& ptetra, list<unsigned>& ntetra) 
{
  // determine the splitting plane: ax + by + cz = d
  Real offset = point[axis];

  // setup the splitting plane
  Plane plane;
  plane.offset = offset;
  if (axis == 0)
    plane.set_normal(Vector3(1,0,0));
  else if (axis == 1)
    plane.set_normal(Vector3(0,1,0));
  else
    plane.set_normal(Vector3(0,0,1));

  // determine the side of the splitting plane of the triangles
  BOOST_FOREACH(unsigned i, otetra)
  {
    // get the three signed distances
    Real sa = plane.calc_signed_distance(_nodes[_tetrahedra[i].a]->x);
    Real sb = plane.calc_signed_distance(_nodes[_tetrahedra[i].b]->x);
    Real sc = plane.calc_signed_distance(_nodes[_tetrahedra[i].c]->x);
    Real sd = plane.calc_signed_distance(_nodes[_tetrahedra[i].d]->x);
    Real min_s = std::min(sa, std::min(sb, std::min(sc, sd)));
    Real max_s = std::max(sa, std::max(sb, std::max(sc, sd)));    

    // see whether we can cleanly put the triangle into one side
    if (min_s > 0)
      ptetra.push_back(i);
    else if (max_s < 0)
      ntetra.push_back(i);
    else
    {
      // tetrahedron is split down the middle; get its centroid
      Tetrahedron tet(_nodes[_tetrahedra[i].a]->x, _nodes[_tetrahedra[i].b]->x, _nodes[_tetrahedra[i].c]->x, _nodes[_tetrahedra[i].d]->x);
      Vector3 tet_centroid = tet.calc_centroid();
      Real scent = plane.calc_signed_distance(tet_centroid);
      if (scent > 0)
        ptetra.push_back(i);
      else
        ntetra.push_back(i);
    }
  }
}

/// Splits an AABB along a given axis into two new AABBs; returns true if split successful
bool DeformableBody::split(AABBPtr source, AABBPtr& tgt1, AABBPtr& tgt2, unsigned axis, const list<unsigned>& tetra, list<unsigned>& ptetra, list<unsigned>& ntetra) 
{
  // clear the lists of tetrahedra
  ptetra.clear();
  ntetra.clear();

  // clear both targets
  tgt1 = shared_ptr<AABB>();
  tgt2 = shared_ptr<AABB>();

  // make sure that not trying to split a single tetrahedron
  assert(tetra.size() > 1); 

  // determine the centroid of this set of tetrahedra
  Vector3 centroid = ZEROS_3;
  Real total_volume = 0;
  BOOST_FOREACH(unsigned idx, tetra)
  {
    Tetrahedron tet = get_tetrahedron(idx);
    Real volume = tet.calc_volume();
    centroid += tet.calc_centroid()*volume;
    total_volume += volume;
  }
  centroid /= total_volume;

  // get the side of the splitting plane of the triangles
  split_tetra(centroid, axis, tetra, ptetra, ntetra);
  if (ptetra.empty() || ntetra.empty())
    return false;

  // get vertices from both meshes
  vector<Vector3> pverts, nverts;
  get_vertices(ptetra.begin(), ptetra.end(), std::back_inserter(pverts));
  get_vertices(ntetra.begin(), ntetra.end(), std::back_inserter(nverts));

  // create two new AABBs
  tgt1 = AABBPtr(new AABB(pverts.begin(), pverts.end()));
  tgt2 = AABBPtr(new AABB(nverts.begin(), nverts.end()));

  return true;
}

// *************************************************************************
// Methods for creating AABB hierarchy around tetrahedral mesh end
// *************************************************************************


