/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/XMLTree.h>
#include <Moby/Optimization.h>
#include <Moby/Event.h>
#include <Moby/PSDeformableBody.h>

using namespace Moby;
using std::list;
using std::map;
using std::vector;
using std::string;
using std::endl;
using boost::shared_ptr;

/// Sets up the deformable body
PSDeformableBody::PSDeformableBody()
{
  default_KP = 1000;
  default_KV = 25;
}

/// Sets the mesh for this body
void PSDeformableBody::set_mesh(shared_ptr<const IndexedTetraArray> tetra_mesh, shared_ptr<Primitive> tri_mesh)
{
  const unsigned EDGES_PER_TETRA = 6;

  // call parent method
  DeformableBody::set_mesh(tetra_mesh, tri_mesh);

  // get the tetrahedron vertices
  const vector<Vector3>& verts = tetra_mesh->get_vertices();

  // now, setup springs vector
  _springs.resize(tetra_mesh->num_tetra() * EDGES_PER_TETRA);
  for (unsigned i=0, j=0; i< _springs.size(); j++) 
  {
    // get the j'th indexed tetrahedron
    const IndexedTetra& tet = tetra_mesh->get_tetra()[j];

    // create a new spring
    Spring s1;
    s1.kp = default_KP;
    s1.kv = default_KV;
    s1.node1 = tet.a;
    s1.node2 = tet.b;
    s1.rest_len = (verts[tet.a] - verts[tet.b]).norm();
    _springs[i++] = s1;

    // create a new spring
    Spring s2;
    s2.kp = default_KP;
    s2.kv = default_KV;
    s2.node1 = tet.a;
    s2.node2 = tet.c;
    s2.rest_len = (verts[tet.a] - verts[tet.c]).norm();
    _springs[i++] = s2;

    // create a new spring
    Spring s3;
    s3.kp = default_KP;
    s3.kv = default_KV;
    s3.node1 = tet.a;
    s3.node2 = tet.d;
    s3.rest_len = (verts[tet.a] - verts[tet.d]).norm();
    _springs[i++] = s3;

    // create a new spring
    Spring s4;
    s4.kp = default_KP;
    s4.kv = default_KV;
    s4.node1 = tet.b;
    s4.node2 = tet.c;
    s4.rest_len = (verts[tet.b] - verts[tet.c]).norm();
    _springs[i++] = s4;

    // create a new spring
    Spring s5;
    s5.kp = default_KP;
    s5.kv = default_KV;
    s5.node1 = tet.b;
    s5.node2 = tet.d;
    s5.rest_len = (verts[tet.b] - verts[tet.d]).norm();
    _springs[i++] = s5;

    // create a new spring
    Spring s6;
    s6.kp = default_KP;
    s6.kv = default_KV;
    s6.node1 = tet.c;
    s6.node2 = tet.d;
    s6.rest_len = (verts[tet.c] - verts[tet.d]).norm();
    _springs[i++] = s6;
  }

  // remove duplicate springs?
  for (unsigned i=0; i< _springs.size(); i++)
    for (unsigned j=i+1; j< _springs.size(); j++)
      if ((_springs[i].node1 == _springs[j].node1 &&
           _springs[i].node2 == _springs[j].node2) ||
          (_springs[i].node1 == _springs[j].node2 &&
           _springs[i].node2 == _springs[j].node1))
      {
        _springs[j] = _springs.back();
        _springs.pop_back();
      }
}

/// Integrates the body forward in time
void PSDeformableBody::integrate(Real t, Real h, shared_ptr<Integrator<VectorN> > integrator)
{
  // don't update c.o.m. or geometries
  disable_config_updates();

  // proceed with integration as normal
  DeformableBody::integrate(t, h, integrator);

  // update c.o.m. and geometries now
  enable_config_updates();
} 

/// Calculates the potential energy of the deformable body
Real PSDeformableBody::calc_potential_energy() const
{
  Real PE = (Real) 0.0;

  // evaluate spring lengths 
  for (unsigned i=0; i< _springs.size(); i++)
  {
    // get the spring properties
    const Spring& s = _springs[i];
    unsigned node1 = s.node1;
    unsigned node2 = s.node2;
    assert(node1 < _nodes.size());
    assert(node2 < _nodes.size());
    Real rest_len = s.rest_len;
    Real kp = s.kp;

    // get the node properties
    const Vector3& x1 = _nodes[node1]->x;
    const Vector3& x2 = _nodes[node2]->x;

    // determine spring forces (from Real Time Physics course notes)
    Vector3 x2mx1 = x2 - x1;
    Real len = x2mx1.norm();
    PE += (Real) 0.5 * (len - rest_len) * (len - rest_len) * kp;
  }

  return PE; 
}

/// Calculates forward dynamics for the deformable body
void PSDeformableBody::calc_fwd_dyn(Real dt)
{
  FILE_LOG(LOG_DEFORM) << "PSDeformableBody::calc_fwd_dyn() entered" << endl;

  // evaluate spring and dampening forces
  for (unsigned i=0; i< _springs.size(); i++)
  {
    // get the spring properties
    const Spring& s = _springs[i];
    unsigned node1 = s.node1;
    unsigned node2 = s.node2;
    assert(node1 < _nodes.size());
    assert(node2 < _nodes.size());
    Real rest_len = s.rest_len;
    Real kp = s.kp;
    Real kv = s.kv;

    // get the node properties
    const Vector3& x1 = _nodes[node1]->x;
    const Vector3& v1 = _nodes[node1]->xd;
    const Vector3& x2 = _nodes[node2]->x;
    const Vector3& v2 = _nodes[node2]->xd;

    // determine spring forces (from Real Time Physics course notes)
    Vector3 x2mx1 = x2 - x1;
    Real len = x2mx1.norm();
    Vector3 dir = x2mx1 / len;
    assert(len > (Real) 0.0); 
    Vector3 f = dir * (len - rest_len) * kp + dir * (v2 - v1).dot(dir) * kv;
    _nodes[node1]->f += f;
    _nodes[node2]->f -= f; 

    FILE_LOG(LOG_DEFORM) << " -- determining force on spring: " << i << endl;
    FILE_LOG(LOG_DEFORM) << "   node (" << node1 << ") position: " << x1 << "  velocity: " << v1 << endl;
    FILE_LOG(LOG_DEFORM) << "   node (" << node2 << ") position: " << x2 << "  velocity: " << v2 << endl;
    FILE_LOG(LOG_DEFORM) << "   current spring length: " << len << endl;
    FILE_LOG(LOG_DEFORM) << "   resting spring length: " << rest_len << endl;
    FILE_LOG(LOG_DEFORM) << "   spring velocity: " << (v2 - v1).dot(dir) << endl;
    FILE_LOG(LOG_DEFORM) << "   force on spring: " << f << endl;
  }

  // determine accelerations
  for (unsigned i=0; i< _nodes.size(); i++)
    _nodes[i]->xdd = _nodes[i]->f / _nodes[i]->mass;

  FILE_LOG(LOG_DEFORM) << "PSDeformableBody::calc_fwd_dyn() exited" << endl;
}

/// Applies an impulse to the deformable body at a given point
void PSDeformableBody::apply_impulse(const Vector3& j, const Vector3& p)
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

  // apply the impulse using the barycentric coordinates
  const IndexedTetra& itet = _tetrahedra[closest];
  _nodes[itet.a]->xd += j * ((Real) 1.0 - u - v - w) / _nodes[itet.a]->mass;
  _nodes[itet.b]->xd += j * u / _nodes[itet.b]->mass;
  _nodes[itet.c]->xd += j * v / _nodes[itet.c]->mass;
  _nodes[itet.d]->xd += j * w / _nodes[itet.d]->mass;

  // clear the force accumulators on the nodes
  _nodes[itet.a]->f = ZEROS_3;
  _nodes[itet.b]->f = ZEROS_3;
  _nodes[itet.c]->f = ZEROS_3;
  _nodes[itet.d]->f = ZEROS_3;
FILE_LOG(LOG_COLDET) << "applied impulse " << j << " to " << p << " uvw: " << u << " " << v << " " << w << " (tetra " << closest << ")" << std::endl;
FILE_LOG(LOG_COLDET) << "  node velocities: " << std::endl;
for (unsigned i=0; i< _nodes.size(); i++)
  FILE_LOG(LOG_COLDET) << "    " << i << ": " << _nodes[i]->xd << std::endl;
}

/// Updates the velocity of this body using impulses computed via event data
void PSDeformableBody::update_velocity(const EventProblemData& q)
{
  // check for easy exit
  if (q.N_CONTACTS == 0)
    return;

  // update velocity using contact impulses
  for (unsigned i=0, s=0; i< q.contact_events.size(); i++)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1.get() != this && sb2.get() != this)
    {
      s += 2;
      continue;
    }

    // get contact point
    const Vector3& p = q.contact_events[i]->contact_point;

    // setup impulse to be along normal direction first 
    Vector3 j = q.contact_events[i]->contact_normal * q.alpha_c[i];

    // now update j with tangent impulses
    j += q.contact_events[i]->contact_tan1*q.beta_c[s++];
    j += q.contact_events[i]->contact_tan2*q.beta_c[s++];

    // see whether to negate the impulse
    if (sb2.get() == this)
      j = -j;

    // apply the impulse
    if (j.norm() > NEAR_ZERO)
      apply_impulse(j, p);
  }
}

/// Adds contributions to the event matrices
void PSDeformableBody::update_event_data(EventProblemData& q) 
{
  if (q.N_CONTACTS == 0)
    return;

  // store the current state of this body
  SAFESTATIC VectorN gc, gv;
  get_generalized_coordinates(DynamicBody::eAxisAngle, gc);
  get_generalized_velocity(DynamicBody::eAxisAngle, gv);

  // setup matrices
  MatrixN Jc_iM_JcT, Jc_iM_DcT, Dc_iM_DcT;
  Jc_iM_JcT.set_zero(q.Jc_iM_JcT.rows(), q.Jc_iM_JcT.columns());
  Jc_iM_DcT.set_zero(q.Jc_iM_DcT.rows(), q.Jc_iM_DcT.columns());
  Dc_iM_DcT.set_zero(q.Dc_iM_DcT.rows(), q.Dc_iM_DcT.columns());
  
  // determine Jc_v and Dc_v
  VectorN Jc_v, Dc_v;
  determine_Jc_v(q.contact_events, Jc_v);
  determine_Dc_v(q.contact_events, Dc_v);

  // update event data Jc_v and Dc_v
  q.Jc_v += Jc_v;
  q.Dc_v += Dc_v;

  // apply impulses and check change in velocities
  VectorN Jc_v_new, Dc_v_new, Jc_v_last, Dc_v_last;
  Jc_v_last.copy_from(Jc_v);
  Dc_v_last.copy_from(Dc_v);
  for (unsigned i=0; i< q.N_CONTACTS; i++)
  {
    // if neither of the bodies in the contact is this, keep looping...
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();
    if (sb1.get() != this && sb2.get() != this)
      continue;
    
    // get the normal
    Vector3 n = q.contact_events[i]->contact_normal;
    if (sb2.get() == this)
      n = -n;

    // apply the impulse along the normal
    apply_impulse(n, q.contact_events[i]->contact_point);

    // determine the change in velocity
    determine_Jc_v(q.contact_events, Jc_v_new);
    determine_Dc_v(q.contact_events, Dc_v_new);
    Jc_v_last -= Jc_v_new;
    Jc_v_last.negate();
    Dc_v_last -= Dc_v_new;
    Dc_v_last.negate();

    // update two of the matrices
    Jc_iM_JcT.set_row(i, Jc_v_last);
    Jc_iM_DcT.set_row(i, Dc_v_last);

    // set Jc_v_last and Dc_dv_last
    Jc_v_last.copy_from(Jc_v_new);
    Dc_v_last.copy_from(Dc_v_new); 
  }

  // now do the same, but for tangent directions only
  for (unsigned i=0, k=0; i< q.N_CONTACTS; i++)
  {
    // if neither of the bodies in the contact is this, keep looping...
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();
    if (sb1.get() != this && sb2.get() != this)
    {
      k += 2;
      continue;
    }

    // evaluate both tangent directions
    Vector3 d1 = q.contact_events[i]->contact_tan1;
    Vector3 d2 = q.contact_events[i]->contact_tan2;
    if (sb2.get() == this)
    {
      d1 = -d1;
      d2 = -d2;
    }

    // apply the impulse in the first tangent direction
    apply_impulse(d1, q.contact_events[i]->contact_point);

    // determine the change in velocity
    determine_Dc_v(q.contact_events, Dc_v_new);
    Dc_v_last -= Dc_v_new;
    Dc_v_last.negate();

    // update the matrix 
    Dc_iM_DcT.set_row(k++, Dc_v_last);
    Dc_v_last.copy_from(Dc_v_new); 

    // apply the impulse in the second tangent direction
    apply_impulse(d2, q.contact_events[i]->contact_point);

    // determine the change in velocity
    determine_Dc_v(q.contact_events, Dc_v_new);
    Dc_v_last -= Dc_v_new;
    Dc_v_last.negate();

    // update the matrix 
    Dc_iM_DcT.set_row(k++, Dc_v_last);
    Dc_v_last.copy_from(Dc_v_new); 
  }

  // update q's matrices
  q.Jc_iM_JcT += Jc_iM_JcT;
  q.Jc_iM_DcT += Jc_iM_DcT;
  q.Dc_iM_DcT += Dc_iM_DcT;

  // restore body state
  set_generalized_coordinates(DynamicBody::eAxisAngle, gc);
  set_generalized_velocity(DynamicBody::eAxisAngle, gv);
}

/// Measures the change in velocity in contact normal directions
void PSDeformableBody::determine_Jc_v(const vector<Event*>& contact_events, VectorN& Jc_v) const
{
  Jc_v.resize(contact_events.size());
  for (unsigned i=0; i< contact_events.size(); i++)
  {
    // get the two single bodies
    SingleBodyPtr sb1 = contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = contact_events[i]->contact_geom2->get_single_body();
    if (sb1.get() != this && sb2.get() != this)
    {
      Jc_v[i] = (Real) 0.0;
      continue;
    }

    // get the velocity at the contact point along the contact normal
    const Vector3& n = contact_events[i]->contact_normal;
    const Vector3& p = contact_events[i]->contact_point;
    Jc_v[i] = n.dot(calc_point_vel(p));
    if (sb2.get() == this)
      Jc_v[i] = -Jc_v[i];
  }
}

/// Measures the change in velocity in contact tangent directions
void PSDeformableBody::determine_Dc_v(const vector<Event*>& contact_events, VectorN& Dc_v) const
{
  // resize the vector
  Dc_v.resize(contact_events.size()*2);

  for (unsigned i=0, k=0; i< contact_events.size(); i++)
  {
    // get the two single bodies
    SingleBodyPtr sb1 = contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = contact_events[i]->contact_geom2->get_single_body();
    if (sb1.get() != this && sb2.get() != this)
    {
      Dc_v[k++] = (Real) 0.0;
      Dc_v[k++] = (Real) 0.0;
      continue;
    }

    // get the contact point
    const Vector3& p = contact_events[i]->contact_point;

    // get the velocity at the contact point along the contact tangents
    Vector3 vel = calc_point_vel(p);
    Real dot1 = contact_events[i]->contact_tan1.dot(vel); 
    Real dot2 = contact_events[i]->contact_tan2.dot(vel); 
    if (sb2.get() == this)
    {
      Dc_v[k++] = -dot1; 
      Dc_v[k++] = -dot2;
    }
    else
    {
      Dc_v[k++] = dot1;
      Dc_v[k++] = dot2;
    } 
  }
}

/// Implements Base::load_from_xml()
void PSDeformableBody::load_from_xml(XMLTreeConstPtr node, map<string, BasePtr>& id_map)
{
  map<std::string, BasePtr>::const_iterator id_iter;

  // verify node name
  assert(strcasecmp(node->name.c_str(), "PSDeformableBody") == 0);

  // read the default spring / damper constants 
  const XMLAttrib* kp_attr = node->get_attrib("default-kp");
  const XMLAttrib* kv_attr = node->get_attrib("default-kv");
  if (kp_attr)
    default_KP = kp_attr->get_real_value();
  if (kv_attr)
    default_KV = kv_attr->get_real_value();

  // load parent data after default springs are set 
  DeformableBody::load_from_xml(node, id_map);

  // read in the default mass, if provided
  const XMLAttrib* def_mass_attr = node->get_attrib("default-mass");
  if (def_mass_attr)
  {
    for (unsigned i=0; i< _nodes.size(); i++)
      _nodes[i]->mass = def_mass_attr->get_real_value();
  }

  // read in the springs
  list<XMLTreeConstPtr> spring_nodes = node->find_child_nodes("Spring");
  if (!spring_nodes.empty())
  {
    // setup a vector of springs
    vector<Spring> springs;

    // read in the springs
    for (list<XMLTreeConstPtr>::const_iterator i = spring_nodes.begin(); i != spring_nodes.end(); i++)
    {
      // create a new spring
      Spring s;

      // read in the spring stiffness
      const XMLAttrib* kp_attr = (*i)->get_attrib("kp");
      if (kp_attr)
        s.kp = kp_attr->get_real_value();

      // read in the spring dampening
      const XMLAttrib* kv_attr = (*i)->get_attrib("kv");
      if (kv_attr)
        s.kv = kv_attr->get_real_value();

      // read in the rest length
      const XMLAttrib* rest_len_attr = (*i)->get_attrib("rest-len");
      if (rest_len_attr)
        s.rest_len = rest_len_attr->get_real_value();

      // read in the first node id
      const XMLAttrib* node1_attr = (*i)->get_attrib("node1");
      if (node1_attr)
        s.node1 = node1_attr->get_unsigned_value();
 
      // read in the second node id
      const XMLAttrib* node2_attr = (*i)->get_attrib("node2");
      if (node2_attr)
        s.node2 = node2_attr->get_unsigned_value();
 
      // save the spring
      springs.push_back(s);
    }

    // save the springs
    _springs = springs;
  }
}

/// Implements Base::save_to_xml()
void PSDeformableBody::save_to_xml(XMLTreePtr node, list<BaseConstPtr>& shared_objects) const
{
  // save parent data
  DeformableBody::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "PSDeformableBody";

  // save the spring and damper constants
  node->attribs.insert(XMLAttrib("default-kp", default_KP));
  node->attribs.insert(XMLAttrib("default-kv", default_KV));

  // save the springs
  for (unsigned i=0; i< _springs.size(); i++)
  {
    XMLTreePtr spring_node(new XMLTree("Spring"));
    spring_node->attribs.insert(XMLAttrib("kp", _springs[i].kp));
    spring_node->attribs.insert(XMLAttrib("kv", _springs[i].kv));
    spring_node->attribs.insert(XMLAttrib("rest-len", _springs[i].rest_len));
    spring_node->attribs.insert(XMLAttrib("node1", _springs[i].node1));
    spring_node->attribs.insert(XMLAttrib("node2", _springs[i].node2));
    node->add_child(spring_node);
  }
}


