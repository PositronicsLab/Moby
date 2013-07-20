/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <boost/make_shared.hpp>
#include <fstream>
#include <Moby/CollisionGeometry.h>
#include <Moby/XMLTree.h>
#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/DeformableBody.h>

using boost::shared_ptr;
using boost::make_shared;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::vector;
using std::map;
using std::list;
using std::endl;
using std::queue;
using std::cerr;
using std::string;
using namespace Ravelin;
using namespace Moby;

DeformableBody::DeformableBody()
{
  // setup the pose
  _F = shared_ptr<Pose3d>(new Pose3d);
  _xd.pose = _F;
  _J.pose = _F;

  // setup the visualization pose - note: it's always relative to the global frame
  _vF = shared_ptr<Pose3d>(new Pose3d);
  _vF->rpose = GLOBAL;

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
void DeformableBody::apply_generalized_impulse( const VectorNd& gj)
{
  const unsigned THREE_D = 3;

  assert(num_generalized_coordinates(DynamicBody::eSpatial) == gj.size());
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    Vector3d j(gj[i*THREE_D], gj[i*THREE_D+1], gj[i*THREE_D+2], GLOBAL);
    _nodes[i]->xd += j/_nodes[i]->mass;
  }
}

/// Adds a generalized force to the deformable body
void DeformableBody::add_generalized_force( const VectorNd& gf)
{
  const unsigned THREE_D = 3;

  assert(num_generalized_coordinates(DynamicBody::eSpatial) == gf.size());
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    Vector3d f(gf[i*THREE_D], gf[i*THREE_D+1], gf[i*THREE_D+2], GLOBAL);
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
  _F->x.set_zero()

;
  _J.m = (double) 0.0;

  // determine the position of the com
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    _F->x = Origin3d(_F->x + Vector3d(_nodes[i]->x) * _nodes[i]->mass);
    _J.m += _nodes[i]->mass;
  }
  _F->x /= _J.m;
  FILE_LOG(LOG_DEFORM) << "new center of mass: " << _F->x << endl;

  // compute the linear velocity
  Vector3d xd = calc_point_vel(Point3d(0.0, 0.0, 0.0, _F));

  // set the inertial offset to zero
  _J.h.set_zero();

  // compute the moment of inertia matrix
  _J.J.set_zero();

  for (unsigned i=0; i< _nodes.size(); i++)
  {
    Vector3d relpoint = _nodes[i]->x - _F->x;
    double xsq = relpoint[X]*relpoint[X];
    double ysq = relpoint[Y]*relpoint[Y];
    double zsq = relpoint[Z]*relpoint[Z];
    _J.J(X,X) += _nodes[i]->mass * (ysq + zsq);
    _J.J(Y,Y) += _nodes[i]->mass * (xsq + zsq);
    _J.J(Z,Z) += _nodes[i]->mass * (xsq + ysq);
    _J.J(X,Y) -= _nodes[i]->mass * relpoint[X] * relpoint[Y];
    _J.J(X,Z) -= _nodes[i]->mass * relpoint[X] * relpoint[Z];
    _J.J(Y,Z) -= _nodes[i]->mass * relpoint[Y] * relpoint[Z];
  }

  // make J symmetric
  _J.J(Y,X) = _J.J(X,Y);
  _J.J(Z,X) = _J.J(X,Z);
  _J.J(Z,Y) = _J.J(Y,Z);

  // invert J
  Matrix3d Jinv = Matrix3d::invert(_J.J);

  // now determine the angular momentum of the body
  Vector3d P(0.0, 0.0, 0.0, _F);
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    Vector3d r = Pose3d::transform(_nodes[i]->x.pose, _F, _nodes[i]->x);
    P += Vector3d::cross(r, _nodes[i]->xd * _nodes[i]->mass);
  }

  // set the angular velocity
  Vector3d omega(Jinv * Origin3d(P), _F);

  // setup the twist
  _xd.set_angular(omega);
  _xd.set_linear(xd);
}

/// Gets the generalized coordinates of the deformable body
VectorNd& DeformableBody::get_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, VectorNd& gc)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  // resize the gc vector (if necessary)
  gc.resize(THREE_D * _nodes.size());

  // populate the gc vector
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    const Point3d& x = _nodes[i]->x;
    gc[j++] = x[X];
    gc[j++] = x[Y];
    gc[j++] = x[Z];
  }

  return gc;
}

/// Gets the generalized velocity of the deformable body
VectorNd& DeformableBody::get_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, VectorNd& gv)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  // resize the gv vector (if necessary)
  gv.resize(THREE_D * _nodes.size());

  // populate the gc vector
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    const Vector3d& xd = _nodes[i]->xd;
    gv[j++] = xd[X];
    gv[j++] = xd[Y];
    gv[j++] = xd[Z];
  }

  return gv;
}

/// Gets the generalized velocity of the deformable body
VectorNd& DeformableBody::get_generalized_acceleration( VectorNd& ga)
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
void DeformableBody::set_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype, const VectorNd& gc)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;

  // verify gc vector is proper size
  assert(gc.size() == THREE_D * _nodes.size());

  // populate the state
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    Point3d& x = _nodes[i]->x;
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
void DeformableBody::set_generalized_velocity(DynamicBody::GeneralizedCoordinateType gctype, const VectorNd& gv)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  // verify gv vector is proper size
  assert(gv.size() == THREE_D * _nodes.size());

  // populate the state
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    Vector3d& xd = _nodes[i]->xd;
    xd[X] = gv[j++];
    xd[Y] = gv[j++];
    xd[Z] = gv[j++];
  }

  // see whether to update configuration-based variables
  if (_config_updates_enabled)
    calc_com_and_vels();
}

/// Solves using the generalized inertia
MatrixNd& DeformableBody::solve_generalized_inertia( const MatrixNd& B, MatrixNd& X)
{
  const unsigned THREE_D = 3;
  assert(B.rows() == THREE_D * _nodes.size());

  // setup X
  X = B;

  // perform the multiplication
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    // multiply three rows
    double imass = (double) 1.0/_nodes[i]->mass;
    CBLAS::scal(B.columns(), imass, &X(j++,0), B.rows());
    CBLAS::scal(B.columns(), imass, &X(j++,0), B.rows());
    CBLAS::scal(B.columns(), imass, &X(j++,0), B.rows());
  }

  return X;
}

/// Solves using the generalized inertia
VectorNd& DeformableBody::solve_generalized_inertia( const VectorNd& b, VectorNd& x)
{
  const unsigned THREE_D = 3;
  assert(b.size() == THREE_D * _nodes.size());

  // setup x
  x = b;

  // perform the multiplication
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    // multiply three rows
    double imass = (double) 1.0/_nodes[i]->mass;
    x[j++] *= imass;
    x[j++] *= imass;
    x[j++] *= imass;
  }

  return x;
}


/// Gets the generalized inertia for the body
MatrixNd& DeformableBody::get_generalized_inertia( MatrixNd& M)
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
VectorNd& DeformableBody::get_generalized_forces( VectorNd& f)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  // resize f (if necessary)
  f.resize(THREE_D * _nodes.size());

  // populate f
  for (unsigned i=0, j=0; i< _nodes.size(); i++)
  {
    const Vector3d& fext = _nodes[i]->f;
    f[j++] = fext[X];
    f[j++] = fext[Y];
    f[j++] = fext[Z];
  }

  return f;
}

/// Converts a force and torque to a generalized force on the body
VectorNd& DeformableBody::convert_to_generalized_force( SingleBodyPtr body, const SForced& w, const Point3d& p, VectorNd& gf)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  gf.set_zero(_nodes.size() * THREE_D);

  // convert the force to the global frame
  SForced w0 = Pose3d::transform(w.pose, GLOBAL, w);
  Vector3d f = w0.get_force();

  // determine in what tetrahedron the point lies (or is closest)
  unsigned i = find_closest_tetrahedron(p);

  // get (and correct) the barycentric coordinates for p
  double bu, bv, bw;
  const Tetrahedron& tet = _tetra_mesh->get_tetrahedron(i);
  tet.determine_barycentric_coords(p, bu, bv, bw);

  // correct barycentric coordinates
  if (bu < 0.0) bu = 0.0;
  else if (bu > 1.0) bu = 1.0;
  if (bv < 0.0) bv = 0.0;
  else if (bu + bv > 1.0) bv = 1.0 - bu;
  if (bw < 0.0) bw = 0.0;
  else if (bu + bv + bw > 1.0) bw = 1.0 - bu - bv;

  // now determine forces on the nodes
  gf.set_sub_vec(_tetrahedra[i].a*3, f*bu);
  gf.set_sub_vec(_tetrahedra[i].b*3, f*bv);
  gf.set_sub_vec(_tetrahedra[i].c*3, f*bw);
  gf.set_sub_vec(_tetrahedra[i].d*3, f*(1.0-bu-bv-bw));

  return gf;
}

/// Resets the force / torque accumulators on the deformable body
void DeformableBody::reset_accumulators()
{
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    // reset the force on the node
    _nodes[i]->f .set_zero()

;
  }
}

/// Transforms the dynamic body
void DeformableBody::translate(const Origin3d& x)
{
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    // transform the node position
    _nodes[i]->x += x;
  }

  // calculate the position of the center-of-mass
  calc_com_and_vels();

  // update the geometries too
  update_geometries();
}

/// Transforms the dynamic body
void DeformableBody::rotate(const Quatd& q)
{
  for (unsigned i=0; i< _nodes.size(); i++)
  {
    // transform the node position
    _nodes[i]->x = Point3d(q * Origin3d(_nodes[i]->x), GLOBAL);
  }

  // calculate the position of the center-of-mass
  calc_com_and_vels();

  // update the geometries too
  update_geometries();
}


/// Calculates the kinetic energy of the deformable body
double DeformableBody::calc_kinetic_energy() const
{
  double KE = (double) 0.0;
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
unsigned DeformableBody::find_closest_tetrahedron(const Point3d& p) const
{
  FILE_LOG(LOG_DEFORM) << "DeformableBody::find_closest_tetrahedron() entered" << endl;

  // convert the point to the global frame
  Point3d p0 = Pose3d::transform(p.pose, GLOBAL, p);

  // get the vector of tetrahedra
  assert(_tetrahedra.size() > 0);

  // closest tetrahedron is none initially
  unsigned closest = std::numeric_limits<unsigned>::max();

  // determine all tetrahedra that p could lie within
  list<unsigned> tet_ids;
  get_tetrahedra(p0, std::back_inserter(tet_ids));
  if (tet_ids.empty())
    FILE_LOG(LOG_DEFORM) << " -- warning: point " << p0 << " not inside any tetrahedra!" << std::endl;

  // if list is empty, determine the closest bounding boxes 
  if (tet_ids.empty())
    for (unsigned i=0; i< _tetrahedra.size(); i++)
      tet_ids.push_back(i);

  FILE_LOG(LOG_DEFORM) << " -- examining " << tet_ids.size() << " / " << _tetra_mesh->num_tetra() << " tetrahedra" << std::endl;

  // setup the minimum distance
  double min_dist = std::numeric_limits<double>::max();

  // iterate through the tetrahedra
  BOOST_FOREACH(unsigned id, tet_ids)
  {
    // get the tetrahedron
    Tetrahedron tet(_nodes[_tetrahedra[id].a]->x, _nodes[_tetrahedra[id].b]->x,
                    _nodes[_tetrahedra[id].c]->x, _nodes[_tetrahedra[id].d]->x);

    FILE_LOG(LOG_DEFORM) << "distance of " << p0 << " to tetrahedron " << id << std::endl;
    FILE_LOG(LOG_DEFORM) << "  " << _nodes[_tetrahedra[id].a]->x << std::endl;
    FILE_LOG(LOG_DEFORM) << "  " << _nodes[_tetrahedra[id].b]->x << std::endl;
    FILE_LOG(LOG_DEFORM) << "  " << _nodes[_tetrahedra[id].c]->x << std::endl;
    FILE_LOG(LOG_DEFORM) << "  " << _nodes[_tetrahedra[id].d]->x << std::endl;

    // determine the distance from the point to the tetrahedron
    double dist = tet.calc_signed_dist(p0);
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
std::ostream& DeformableBody::to_vrml(std::ostream& o, double scale) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup random colors for each tetrahedron
  std::vector<Point3d> colors(_tetrahedra.size());
  for (unsigned i=0; i< colors.size(); i++)
  {
    colors[i][0] = (double) rand() / RAND_MAX;
    colors[i][1] = (double) rand() / RAND_MAX;
    colors[i][2] = (double) rand() / RAND_MAX;
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
    const Point3d& v = ita->get_vertices()[i];
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

  return o;
}

/// Sets the mesh for the deformable body
/**
 * \note Also sets the nodes for the body (and thereby modifies the com 
 *       position and velocity.
 */
void DeformableBody::set_mesh(shared_ptr<const IndexedTetraArray> tetra_mesh, shared_ptr<Primitive> tri_mesh)
{
  // set the default mass
  const double DEFAULT_MASS = (double) 1.0;

  // store the triangle mesh
  _tri_mesh = tri_mesh; 

  // store the tetrahedra
  _tetrahedra = tetra_mesh->get_tetra();

  // transform the tetrahedral mesh 
  Transform3d T = Pose3d::calc_relative_pose(tetra_mesh->get_pose(), GLOBAL);
  _tetra_mesh = make_shared<const IndexedTetraArray>(tetra_mesh->transform(T)); 

  // get the vertices of the meshes
  const vector<Point3d>& tet_vertices = tetra_mesh->get_vertices();
  const vector<Point3d>& tri_vertices = tri_mesh->get_mesh()->get_vertices();

  // setup the nodes in the body based on the tetrahedra vertices
  _nodes = vector<shared_ptr<Node> >(tet_vertices.size());
  for (unsigned i=0; i< tet_vertices.size(); i++)
  {
    _nodes[i] = shared_ptr<Node>(new Node);
    _nodes[i]->x = tet_vertices[i];
    _nodes[i]->xd = _nodes[i]->xdd = Vector3d(0.0, 0.0, 0.0, GLOBAL);
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
    // verify that triangle mesh vertex is in global frame (this is my
    // assumption, but I'm not sure of the implication if it's false)
    assert(tri_vertices[i].pose == GLOBAL);

    // get the closest tetrahedron
    unsigned closest = find_closest_tetrahedron(tri_vertices[i]);

    // should be one that is closest
    assert(closest < std::numeric_limits<unsigned>::max());

    // mark the tetrahedron in the vertex map
    _vertex_map[i].tetra = closest; 

    // get the tetrahedron
    Tetrahedron tet = get_tetrahedron(closest);

    // determine the barycentric coordinates
    double u, v, w;
    tet.determine_barycentric_coords(tri_vertices[i], u, v, w);

    // correct barycentric coordinates
    if (u < 0.0) u = 0.0;
    else if (u > 1.0) u = 1.0;
    if (v < 0.0) v = 0.0;
    else if (u + v > 1.0) v = 1.0 - u;
    if (w < 0.0) w = 0.0;
    else if (u + v + w > 1.0) w = 1.0 - u - v;

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
  update_geometries(); 
}

/// Updates collision and visualization geometries using present state of the nodes of the deformable body
void DeformableBody::update_geometries()
{
  // setup the vertices for the triangle mesh
  vector<Point3d> vertices(_tri_mesh->get_mesh()->get_vertices().size());
 
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
Vector3d DeformableBody::calc_point_vel(const Point3d& p) const
{
  // find the closest tetrahedron
  unsigned closest = find_closest_tetrahedron(p);

  // should be a closest 
  assert(closest != std::numeric_limits<unsigned>::max());

  // get the tetrahedron
  Tetrahedron tet = get_tetrahedron(closest);

  // determine the barycentric coordinates
  double u, v, w;
  tet.determine_barycentric_coords(p, u, v, w);

  // get the velocities at the vertices
  const IndexedTetra& itet = _tetrahedra[closest];
  const Vector3d& va = _nodes[itet.a]->xd;
  const Vector3d& vb = _nodes[itet.b]->xd;
  const Vector3d& vc = _nodes[itet.c]->xd;
  const Vector3d& vd = _nodes[itet.d]->xd;

  // determine the velocity at p using the barycentric function
  Vector3d vx = vb*u + vc*v + vd*w + va*(1.0 - u - v - w);

  // convert the velocity to the desired frame
  return Pose3d::transform(vx.pose, p.pose, vx);
}

/// Adds a force to the body
void DeformableBody::add_force(const Vector3d& f)
{
  for (unsigned i=0; i< _nodes.size(); i++)
    _nodes[i]->f += f;
}

/// Adds a force to the deformable body at a given point
void DeformableBody::add_force(const Vector3d& f, const Point3d& p)
{
  // find the closest tetrahedron
  unsigned closest = find_closest_tetrahedron(p);

  // should be a closest 
  assert(closest != std::numeric_limits<unsigned>::max());

  // get the tetrahedron
  Tetrahedron tet = get_tetrahedron(closest);

  // determine the barycentric coordinates
  double u, v, w;
  tet.determine_barycentric_coords(p, u, v, w);

  // apply the force using the barycentric coordinates
  const IndexedTetra& itet = _tetrahedra[closest];
  _nodes[itet.a]->f += f * ((double) 1.0 - u - v - w);
  _nodes[itet.b]->f += f * u;
  _nodes[itet.c]->f += f * v;
  _nodes[itet.c]->f += f * w;
}

/// Implements Base::load_from_xml()
void DeformableBody::load_from_xml(shared_ptr<const XMLTree> node, map<string, BasePtr>& id_map)
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
  list<shared_ptr<const XMLTree> > nodes_nodes = node->find_child_nodes("Node");
  if (!nodes_nodes.empty())
  {
    // setup a vector of nodes
    vector<shared_ptr<Node> > new_nodes;

    // read in the set of nodes
    for (list<shared_ptr<const XMLTree> >::const_iterator i = nodes_nodes.begin(); i != nodes_nodes.end(); i++)
    {
      // create a new Node
      shared_ptr<Node> n(new Node);

      // read the node position
      const XMLAttrib* node_x_attr = (*i)->get_attrib("position");
      if (node_x_attr)
        n->x = node_x_attr->get_point_value();

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
  list<shared_ptr<const XMLTree> > cg_nodes = node->find_child_nodes("CollisionGeometry");

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
  list<shared_ptr<const XMLTree> > vmap_nodes = node->find_child_nodes("VertexMap");
  if (!vmap_nodes.empty())
  {
    assert(vmap_nodes.size() == 1);
    list<shared_ptr<const XMLTree> > map_nodes = vmap_nodes.front()->find_child_nodes("Mapping");
    BOOST_FOREACH(shared_ptr<const XMLTree> mapping_node, map_nodes)
    {
      const XMLAttrib* tetra_id_attr = node->get_attrib("tetra");
      const XMLAttrib* uvw_id_attr = node->get_attrib("uvw");
      if (tetra_id_attr && uvw_id_attr)
      {
        VertexMap vmap;
        Point3d uvw = uvw_id_attr->get_point_value();
        vmap.tetra = tetra_id_attr->get_unsigned_value();
        vmap.uvw[0] = uvw[0];
        vmap.uvw[1] = uvw[1];
        vmap.uvw[2] = uvw[2];
        vertex_map.push_back(vmap);
      }
    }
  }

  // only if both pointers are set do we set the mesh
  if (tetra_mesh && tri_mesh)
  {
    // set the default mass
    const double DEFAULT_MASS = (double) 1.0;

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
        _nodes[i]->xd = _nodes[i]->xdd .set_zero()

;
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
        double u, v, w;
        tet.determine_barycentric_coords(tri_vertices[i], u, v, w);

        // correct barycentric coordinates
        if (u < (double) 0.0) u = (double) 0.0;
        else if (u > (double) 1.0) u = (double) 1.0;
        if (v < (double) 0.0) v = (double) 0.0;
        else if (u + v > (double) 1.0) v = (double) 1.0 - u;
        if (w < (double) 0.0) w = (double) 0.0;
        else if (u + v + w > (double) 1.0) w = (double) 1.0 - u - v;

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
    update_geometries(); 
*/
  }

  // read a transform to be applied to the body, if provided
  const XMLAttrib* xlat_attr = node->get_attrib("translation");
  const XMLAttrib* rpy_attr = node->get_attrib("rpy");
  const XMLAttrib* quat_attr = node->get_attrib("quat");
  if (xlat_attr && rpy_attr)
  {
    rotate(rpy_attr->get_rpy_value());
    translate(xlat_attr->get_origin_value());
  }
  else if (xlat_attr && quat_attr)
  {
    rotate(rpy_attr->get_quat_value());
    translate(xlat_attr->get_origin_value());
  }
  else if (xlat_attr)
    translate(xlat_attr->get_origin_value());
  else if (rpy_attr)
    rotate(rpy_attr->get_rpy_value());
  else if (quat_attr)
    rotate(rpy_attr->get_quat_value());
}

/// Implements Base::save_to_xml()
void DeformableBody::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
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
    Point3d uvw(_vertex_map[i].uvw[0], _vertex_map[i].uvw[1], _vertex_map[i].uvw[2]);
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
  const double EXP = 0.1;          // box expansion constant
  AABBPtr child1, child2;
  list<unsigned> ptetra, ntetra;
  Vector3d eps(EXP, EXP, EXP);

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
    double vol = aabb->calc_volume();
    bool erased_one = false;
    for (list<BVPtr>::iterator i = aabb->children.begin(); i != aabb->children.end(); )
    {
      // get the AABB
      AABBPtr child = dynamic_pointer_cast<AABB>(*i);

      // get the volume of this child
      double voli = child->calc_volume();
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
      {
        if (!child->is_leaf())
          Q.push(dynamic_pointer_cast<AABB>(child));
      }
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
    {
      BOOST_FOREACH(BVPtr child, aabb->children)
        Q.push(dynamic_pointer_cast<AABB>(child));
    }
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
    {
      BOOST_FOREACH(BVPtr child, aabb->children)
        Q.push(dynamic_pointer_cast<AABB>(child));
    }

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
void DeformableBody::split_tetra(const Point3d& point, unsigned axis, const list<unsigned>& otetra, list<unsigned>& ptetra, list<unsigned>& ntetra) 
{
  // determine the splitting plane: ax + by + cz = d
  double offset = point[axis];

  // setup the splitting plane
  Plane plane;
  plane.offset = offset;
  if (axis == 0)
    plane.set_normal(Vector3d(1,0,0));
  else if (axis == 1)
    plane.set_normal(Vector3d(0,1,0));
  else
    plane.set_normal(Vector3d(0,0,1));

  // determine the side of the splitting plane of the triangles
  BOOST_FOREACH(unsigned i, otetra)
  {
    // get the three signed distances
    double sa = plane.calc_signed_distance(_nodes[_tetrahedra[i].a]->x);
    double sb = plane.calc_signed_distance(_nodes[_tetrahedra[i].b]->x);
    double sc = plane.calc_signed_distance(_nodes[_tetrahedra[i].c]->x);
    double sd = plane.calc_signed_distance(_nodes[_tetrahedra[i].d]->x);
    double min_s = std::min(sa, std::min(sb, std::min(sc, sd)));
    double max_s = std::max(sa, std::max(sb, std::max(sc, sd)));    

    // see whether we can cleanly put the triangle into one side
    if (min_s > 0)
      ptetra.push_back(i);
    else if (max_s < 0)
      ntetra.push_back(i);
    else
    {
      // tetrahedron is split down the middle; get its centroid
      Tetrahedron tet(_nodes[_tetrahedra[i].a]->x, _nodes[_tetrahedra[i].b]->x, _nodes[_tetrahedra[i].c]->x, _nodes[_tetrahedra[i].d]->x);
      Point3d tet_centroid = tet.calc_centroid();
      double scent = plane.calc_signed_distance(tet_centroid);
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
  Point3d centroid;
  centroid.set_zero();

  double total_volume = 0;
  BOOST_FOREACH(unsigned idx, tetra)
  {
    Tetrahedron tet = get_tetrahedron(idx);
    double volume = tet.calc_volume();
    centroid += Vector3d(tet.calc_centroid())*volume;
    total_volume += volume;
  }
  centroid /= total_volume;

  // get the side of the splitting plane of the triangles
  split_tetra(centroid, axis, tetra, ptetra, ntetra);
  if (ptetra.empty() || ntetra.empty())
    return false;

  // get vertices from both meshes
  vector<Point3d> pverts, nverts;
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


