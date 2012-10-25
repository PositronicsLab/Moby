/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifdef _OPENMP
#include <omp.h>
#endif
#include <fstream>
#include <set>
#include <cmath>
#include <algorithm>
#include <stack>
#include <queue>
#include <boost/tuple/tuple.hpp>
#include <Moby/CompGeom.h>
#include <Moby/CollisionDetection.h>
#include <Moby/LinAlg.h>
#include <Moby/Event.h>
#include <Moby/Constants.h>
#include <Moby/Polyhedron.h>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/CollisionGeometry.h>  
#include <Moby/XMLTree.h>
#include <Moby/Integrator.h>
#include <Moby/OBB.h>
#include <Moby/GeneralizedCCD.h>
#include <Moby/EventDrivenSimulator.h>

// To delete
#include <Moby/SpherePrimitive.h>
#include <Moby/BoxPrimitive.h>

using namespace Moby;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using boost::shared_ptr;
using boost::tuple;
using boost::make_tuple;
using std::cerr;
using std::endl;
using std::set;
using std::queue;
using std::map;
using std::list;
using std::vector;
using std::priority_queue;
using std::pair;
using std::make_pair;

/// Constructs a collision detector with default tolerances
/**
 * eps_tolerance is set to NEAR_ZERO
 */
GeneralizedCCD::GeneralizedCCD()
{
  eps_tolerance = std::sqrt(std::numeric_limits<Real>::epsilon());
  _max_dexp = std::numeric_limits<unsigned>::max();
  pthread_mutex_init(&_contact_mutex, NULL);
  pthread_mutex_init(&_ve_BVs_mutex, NULL);
  _rebuild_bounds_vecs = true;
  return_all_contacts = true;
}

void GeneralizedCCD::add_collision_geometry(CollisionGeometryPtr cg)
{
  CollisionDetection::add_collision_geometry(cg);
  _rebuild_bounds_vecs = true;
}

void GeneralizedCCD::add_rigid_body(RigidBodyPtr rb)
{
  CollisionDetection::add_rigid_body(rb);
  _rebuild_bounds_vecs = true;
}

void GeneralizedCCD::add_articulated_body(ArticulatedBodyPtr abody, bool disable_adjacent)
{
  CollisionDetection::add_articulated_body(abody, disable_adjacent);
  _rebuild_bounds_vecs = true;
}

void GeneralizedCCD::remove_collision_geometry(CollisionGeometryPtr cg)
{
  CollisionDetection::remove_collision_geometry(cg);
  _rebuild_bounds_vecs = true;
}

void GeneralizedCCD::remove_all_collision_geometries()
{
  CollisionDetection::remove_all_collision_geometries();
  _rebuild_bounds_vecs = true;
}

void GeneralizedCCD::remove_rigid_body(RigidBodyPtr rb)
{
  CollisionDetection::remove_rigid_body(rb);
  _rebuild_bounds_vecs = true;
}

void GeneralizedCCD::remove_articulated_body(ArticulatedBodyPtr abody)
{
  CollisionDetection::remove_articulated_body(abody);
  _rebuild_bounds_vecs = true;
}

/// Computes the velocities from states
map<SingleBodyPtr, pair<Vector3, Vector3> > GeneralizedCCD::get_velocities(const vector<pair<DynamicBodyPtr, VectorN> >& q0, const vector<pair<DynamicBodyPtr, VectorN> >& q1, Real dt) const
{
  // first set the generalized velocities
  #ifndef _OPENMP
  VectorN qd;
  for (unsigned i=0; i< q0.size(); i++)
  {
    qd.copy_from(q1[i].second) -= q0[i].second;
    // Q: why do we set generalized coords to q1?
    q1[i].first->set_generalized_coordinates(DynamicBody::eRodrigues, q0[i].second);
    q1[i].first->set_generalized_velocity(DynamicBody::eRodrigues, qd);
  }
  #else
  SAFESTATIC vector<VectorN> qd;
  qd.resize(q0.size());
  #pragma #omp parallel for
  for (unsigned i=0; i< q0.size(); i++)
  {
    qd[i].copy_from(q1[i].second) -= q0[i].second;
    q1[i].first->set_generalized_coordinates(DynamicBody::eRodrigues, q0[i].second);
    q1[i].first->set_generalized_velocity(DynamicBody::eRodrigues, qd[i]);
  }
  #endif

  // setup the map
  map<SingleBodyPtr, pair<Vector3, Vector3> > vels;

  // now get the single body linear and angular velocities
  for (unsigned i=0; i< q0.size(); i++)
  {
    SingleBodyPtr sb = dynamic_pointer_cast<SingleBody>(q0[i].first);
    if (sb)
    {
      pair<Vector3, Vector3>& v = vels[sb];
      v.first = sb->get_lvel();
      v.second = sb->get_avel();
    }
    else
    {
      ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(q0[i].first);
      assert(ab);
      const vector<RigidBodyPtr>& links = ab->get_links();
      for (unsigned j=0; j< links.size(); j++)
      {
        pair<Vector3, Vector3>& v = vels[links[j]];
        v.first = links[j]->get_lvel();
        v.second = links[j]->get_avel();
      }
    } 
  }

  return vels;
}

#ifndef _OPENMP
/// Determines whether there is a contact in the given time interval
/**
 * \pre body states are at time tf
 */
bool GeneralizedCCD::is_contact(Real dt, const vector<pair<DynamicBodyPtr, VectorN> >& q0, const vector<pair<DynamicBodyPtr, VectorN> >& q1, vector<Event>& contacts)
{
  DStruct ds;
  typedef pair<CollisionGeometryPtr, BVPtr> CG_BV;

  // clear the vector of events 
  contacts.clear();

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::is_contact() entered" << endl;

  // get the map of bodies to velocities
  // NOTE: this also sets each body's coordinates and velocities to q0
  map<SingleBodyPtr, pair<Vector3, Vector3> > vels = get_velocities(q0, q1, dt);

  // clear all velocity expanded BVs
  _ve_BVs.clear();
  BOOST_FOREACH(CollisionGeometryPtr cg, _geoms)
    _ve_BVs[cg] = map<BVPtr, BVPtr>();

  // do broad phase; NOTE: broad phase yields updated BVs
  vector<pair<CollisionGeometryPtr, CollisionGeometryPtr> > to_check;
  broad_phase(vels, to_check);

  // check the geometries
  for (unsigned i=0; i< to_check.size(); i++)
  {
    // get the two geometries
    CollisionGeometryPtr a = to_check[i].first;
    CollisionGeometryPtr b = to_check[i].second;

    // get the two rigid bodies
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(a->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(b->get_single_body());

    // verify that the two bodies are rigid
    if (!rba || !rbb)
      throw std::runtime_error("One or more bodies is not rigid; GeneralizedCCD only works with rigid bodies");

    // get the velocities for the two bodies
    const pair<Vector3, Vector3>& a_vel = vels.find(rba)->second;
    const pair<Vector3, Vector3>& b_vel = vels.find(rbb)->second;

    // get the transforms from a to b and back
    Matrix4 aTb = Matrix4::inverse_transform(a->get_transform()) * b->get_transform();
    Matrix4 bTa = Matrix4::inverse_transform(b->get_transform()) * a->get_transform(); 

    // test the geometries for contact
    check_geoms(dt, a, b, aTb, bTa, a_vel, b_vel, contacts);
  } 

  FILE_LOG(LOG_COLDET) << "contacts:" << endl;
  if (contacts.empty())
    FILE_LOG(LOG_COLDET) << " -- no contacts in narrow phase" << endl;
  if (LOGGING(LOG_COLDET))
    for (unsigned i=0; i< contacts.size(); i++)
      FILE_LOG(LOG_COLDET) << contacts[i] << std::endl;

  // sort the vector of events
  std::sort(contacts.begin(), contacts.end());

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::is_contact() exited" << endl << endl;

  // indicate whether impact has occurred
  return !contacts.empty();
}
#else
/// Determines whether there is a contact in the given time interval
/**
 * \pre body states are at time tf
 */
bool GeneralizedCCD::is_contact(Real dt, const vector<pair<DynamicBodyPtr, VectorN> >& q0, const vector<pair<DynamicBodyPtr, VectorN> >& q1, vector<Event>& contacts)
{
  DStruct ds;
  typedef pair<CollisionGeometryPtr, BVPtr> CG_BV;

  // clear the contacts vector 
  contacts.clear();

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::is_contact() entered" << endl;

  // get the map of bodies to velocities
  // NOTE: this also sets each body's coordinates and velocities to q0
  map<SingleBodyPtr, pair<Vector3, Vector3> > vels = get_velocities(q0, q1, dt);

  // clear all velocity expanded BVs
  _ve_BVs.clear();
  BOOST_FOREACH(CollisionGeometryPtr cg, _geoms)
    _ve_BVs[cg] = map<BVPtr, BVPtr>();

  // do broad phase; NOTE: broad phase yields updated BVs
  vector<pair<CollisionGeometryPtr, CollisionGeometryPtr> > to_check;
  broad_phase(vels, to_check);
  vector<vector<Event> > events(to_check.size());

  // check the geometries
  #pragma omp parallel for
  for (int i=0; i< (int) to_check.size(); i++)
  {
    // get the two geometries
    CollisionGeometryPtr a = to_check[i].first;
    CollisionGeometryPtr b = to_check[i].second;

    // get the two rigid bodies
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(a->get_single_body());
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(b->get_single_body());

    // verify that the two bodies are rigid
    if (!rba || !rbb)
      throw std::runtime_error("One or more bodies is not rigid; GeneralizedCCD only works with rigid bodies");

    // get the velocities for the two bodies
    const pair<Vector3, Vector3>& a_vel = vels.find(rba)->second;
    const pair<Vector3, Vector3>& b_vel = vels.find(rbb)->second;

    // get the transforms from a to b and back
    Matrix4 aTb = Matrix4::inverse_transform(a->get_transform()) * b->get_transform();
    Matrix4 bTa = Matrix4::inverse_transform(b->get_transform()) * a->get_transform(); 

    // test the geometries for contact
    check_geoms(dt, a, b, aTb, bTa, a_vel, b_vel, events[i]);
  } 

  // integrate all contacts into a single structure 
  for (unsigned i=0; i< events.size(); i++)
    contacts.insert(contacts.end(), events[i].begin(), events[i].end());

  FILE_LOG(LOG_COLDET) << "contacts:" << endl;
  if (contacts.empty())
    FILE_LOG(LOG_COLDET) << " -- no contacts in narrow phase" << endl;
  if (LOGGING(LOG_COLDET))
    for (unsigned i=0; i< contacts.size(); i++)
      FILE_LOG(LOG_COLDET) << contacts[i] << std::endl;

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::is_contact() exited" << endl << endl;

  // indicate whether impact has occurred
  return !contacts.empty();
}
#endif

/// Does a collision check for a pair of geometries (serial version)
/**
 * \param a the first geometry
 * \param b the second geometry
 * \param aTb the transform from b's frame to a's frame
 * \param bTa the transform from a's frame to b's frame
 * \param vels linear and angular velocities of bodies
 * \param contacts on return
 */
void GeneralizedCCD::check_geoms(Real dt, CollisionGeometryPtr a, CollisionGeometryPtr b, const Matrix4& aTb, const Matrix4& bTa, const pair<Vector3, Vector3>& a_vel, const pair<Vector3, Vector3>& b_vel, vector<Event>& contacts)
{
  map<BVPtr, vector<const Vector3*> > a_to_test, b_to_test;

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::check_geoms() entered" << endl;
  FILE_LOG(LOG_COLDET) << "  checking geometry " << a->id << " for body " << a->get_single_body()->id << std::endl;
  FILE_LOG(LOG_COLDET) << "  against geometry " << b->id << " for body " << b->get_single_body()->id << std::endl;

  // init statistics for geometry pair
  unsigned n_bv_tests = 0;
  unsigned n_verts_tested = 0;

  // set the earliest TOC
  Real earliest = (Real) 1.0; 

  // setup the contact set for these two geometries 
  vector<Event> local_contacts;

  // get the velocities for a and b
  const Vector3& alv = a_vel.first;
  const Vector3& aav = a_vel.second;
  const Vector3& blv = b_vel.first;
  const Vector3& bav = b_vel.second;

  // get the primitives for a and b
  PrimitivePtr aprimitive = a->get_geometry();
  PrimitivePtr bprimitive = b->get_geometry();

  // get the two top-level BVs for a and b
  BVPtr bv_a = aprimitive->get_BVH_root();
  BVPtr bv_b = bprimitive->get_BVH_root(); 

  // add the two top-level BVs to the queue for processing
  queue<BVProcess> q;
  q.push(BVProcess());
  q.back().bva = get_vel_exp_BV(a, bv_a, alv, aav);
  q.back().bvb = get_vel_exp_BV(b, bv_b, blv, bav);
  q.back().nexp = 0;
  q.back().ax = bv_a;
  q.back().bx = bv_b;

  // process until the queue is empty
  while (!q.empty())
  {
    // get the number of expansions  
    unsigned nexp = q.front().nexp;

    // get the BVs off of the queue
    BVPtr bva = q.front().bva;
    BVPtr bvb = q.front().bvb;
    BVPtr ax = q.front().ax;
    BVPtr bx = q.front().bx;

    // increment the # of bounding volume tests
    n_bv_tests++;

    // if the velocity-expanded BVs do not intersect, continue looping
    if (!BV::intersects(bva, bvb, aTb))
    {
      FILE_LOG(LOG_COLDET) << " -- bounding volumes do not intersect" << endl;
      q.pop();
      continue;
    }

    // calculate the volume for the OBBs
    Real ax_vol = ax->calc_volume();
    Real bx_vol = bx->calc_volume();

    // velocity expanded BVs do intersect; if both BVs are leafs OR maximum
    // depth has been reached, intersect
    // the vertices of triangles in one BV against the triangles of the other
    // (and vice versa)
    if ((ax->is_leaf() && bx->is_leaf()) || nexp >= _max_dexp)
    {
      // get testing vertex vectors for ax and bx
      vector<const Vector3*>& ax_test = a_to_test[ax];
      vector<const Vector3*>& bx_test = b_to_test[bx];

      // get the sets of vertices for ax and bx
      aprimitive->get_vertices(ax, bx_test);
      bprimitive->get_vertices(bx, ax_test);
    }
    // descend into BV with greater volume
    else if ((ax_vol > bx_vol && !ax->is_leaf()) || bx->is_leaf())
    {
      // add all children of bva to the queue for processing
      BOOST_FOREACH(BVPtr bv, ax->children)
      {
        q.push(BVProcess());
        q.back().bva = get_vel_exp_BV(a, bv, alv, aav);
        q.back().bvb = bvb; 
        q.back().nexp = nexp+1;
        q.back().ax = bv;
        q.back().bx = bx;
      }
    }
    else
    {
      // add all children of bvb to the queue for processing
      BOOST_FOREACH(BVPtr bv, bx->children)
      {
        q.push(BVProcess());
        q.back().bva = bva; 
        q.back().bvb = get_vel_exp_BV(b, bv, blv, bav);
        q.back().nexp = nexp+1;
        q.back().ax = ax;
        q.back().bx = bv;
      }
    }

    // pop the element off of the front of the queue
    q.pop();
  }

  // make vectors of vertices unique
  for (map<BVPtr, vector<const Vector3*> >::iterator i = a_to_test.begin(); i != a_to_test.end(); i++)
  {
    std::sort(i->second.begin(), i->second.end());
    i->second.erase(std::unique(i->second.begin(), i->second.end()), i->second.end());
  }
  for (map<BVPtr, vector<const Vector3*> >::iterator i = b_to_test.begin(); i != b_to_test.end(); i++)
  {
    std::sort(i->second.begin(), i->second.end());
    i->second.erase(std::unique(i->second.begin(), i->second.end()), i->second.end());
  }

  // call check_vertices on bounding volumes from geometry b
  for (map<BVPtr, vector<const Vector3*> >::iterator i = b_to_test.begin(); i != b_to_test.end(); i++)
  {
    n_verts_tested += i->second.size();
    check_vertices(dt, a, b, i->first, i->second, bTa, a_vel, b_vel, earliest, local_contacts);
  }

  // call check_vertices on bounding volumes from geometry a
  for (map<BVPtr, vector<const Vector3*> >::iterator i = a_to_test.begin(); i != a_to_test.end(); i++)
  {
    n_verts_tested += i->second.size();
    check_vertices(dt, b, a, i->first, i->second, aTb, b_vel, a_vel, earliest, local_contacts);
  }

  if (LOGGING(LOG_COLDET))
  {
    SingleBodyPtr sba = a->get_single_body(); 
    SingleBodyPtr sbb = b->get_single_body(); 
    FILE_LOG(LOG_COLDET) << " -- contact points for bodies " << sba->id << " and " << sbb->id << endl;
    for (unsigned i=0; i< local_contacts.size(); i++)
      FILE_LOG(LOG_COLDET) << local_contacts[i] << std::endl;
  }

  // get the time-of-impact tolerance
  shared_ptr<EventDrivenSimulator> sim(simulator);
  const Real TOI_TOLERANCE = sim->toi_tolerance;

  // sort the vector of contacts
  std::sort(local_contacts.begin(), local_contacts.end());

  // see what points to insert into the global set of contacts 
  // we return all contacts if using an event-driven method to deal with Zeno
  // points 
  if (!return_all_contacts)
  {
    for (unsigned i=0; i< local_contacts.size(); i++)
      if (local_contacts[i].t > earliest + TOI_TOLERANCE/dt)
      {
        contacts.insert(contacts.end(), local_contacts.begin(), local_contacts.begin()+i);
        break;
      }
  }
  else
  {
    #ifdef _OPENMP
    pthread_mutex_lock(&_contact_mutex);
    #endif

    // insert the contacts from these geometries into the contact map of all
    // geometries
    contacts.insert(contacts.end(), local_contacts.begin(), local_contacts.end());

    #ifdef _OPENMP
    pthread_mutex_unlock(&_contact_mutex);
    #endif
  }

  if (LOGGING(LOG_COLDET))
  {
    // determine how many pairs of bounding volumes
    std::vector<BVPtr> all_bvs;
    bv_a->get_all_BVs(std::back_inserter(all_bvs));
    unsigned n_a_bvs = all_bvs.size();
    all_bvs.clear();
    bv_b->get_all_BVs(std::back_inserter(all_bvs));
    unsigned n_b_bvs = all_bvs.size();
    FILE_LOG(LOG_COLDET) << "  number of BV vs. BV tests: " << n_bv_tests << "/" << (n_a_bvs * n_b_bvs) << std::endl;

    // output number of vertices tested 
    FILE_LOG(LOG_COLDET) << "  number of vertices tested: " << n_verts_tested << std::endl;
  }
  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::check_geoms() exited" << endl;
}

/// Checks a set of vertices of geometry a against geometry b
void GeneralizedCCD::check_vertices(Real dt, CollisionGeometryPtr a, CollisionGeometryPtr b, BVPtr bvb, const std::vector<const Vector3*>& a_verts, const Matrix4& bTa, const pair<Vector3, Vector3>& a_vel, const pair<Vector3, Vector3>& b_vel, Real& earliest, vector<Event>& local_contacts) const
{
  Vector3 point, normal;

  // get the time-of-impact tolerance
  shared_ptr<EventDrivenSimulator> sim(simulator);
  const Real TOI_TOLERANCE = sim->toi_tolerance;

  // get the two bodies
  RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(a->get_single_body()); 
  RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(b->get_single_body()); 

  // get the velocities for a and b
  const Vector3& alv = a_vel.first;
  const Vector3& aav = a_vel.second;
  const Vector3& blv = b_vel.first;
  const Vector3& bav = b_vel.second;

  FILE_LOG(LOG_COLDET) << "  -- checking body " << rba->id << " against " << rbb->id << endl;
  FILE_LOG(LOG_COLDET) << "    -- relative transform " << endl << bTa;

  // setup a "queue" for checking vertices
  vector<pair<Real, pair<Vector3, Vector3> > > Q;
  Q.clear();
  BOOST_FOREACH(const Vector3* v, a_verts)
  {
    // compute point at time t0 in b coordinates
    FILE_LOG(LOG_COLDET) << "    -- checking vertex " << *v << " of " << rba->id << " against " << rbb->id << endl;
    Vector3 p0 = bTa.mult_point(*v);

    FILE_LOG(LOG_COLDET) << "     -- p0 (local): " << p0 << endl; 
    FILE_LOG(LOG_COLDET) << "     -- p0 (global): " << (b->get_transform().mult_point(p0)) << endl;

    // we'll sort on inverse distance from the center of mass (origin of b frame) 
    Real dist = 1.0/p0.norm_sq();

    // push the vertices onto the queue
    // NOTE: assumes that center of geometry of body a is its C.O.M.
    //       (if the assumption is wrong, this will only make things slower)
    Q.push_back(make_pair(dist, make_pair(p0, *v)));
  }

  // sort the queue
  std::sort(Q.begin(), Q.end());

  // populate the DStruct for checking vertices of a against b
  DStruct ds;
  populate_dstruct(&ds, a, b, alv, aav, blv, bav, bvb);

  // check all vertices of a against b
  while (!Q.empty())
  {
    // setup the vertex and u in the DStruct
    ds.p0 = Q.back().second.first;
    ds.u = Q.back().second.second;

    FILE_LOG(LOG_COLDET) << "    -- checking vertex p0 (local) of " << rba->id << ": " <<  ds.p0 << "  u: " << ds.u << endl; 
    FILE_LOG(LOG_COLDET) << "     -- global pos (t0): " << b->get_transform().mult_point(ds.p0) << endl;

    // determine TOI, if any
    // NOTE: this has been modified to account for the method of dealing with
    // Zeno points in ImpulseContactSimulator; it is considerably slower as a
    // result
    Real toi;
    if (!return_all_contacts)
      toi = determine_TOI(0.0, earliest, &ds, point, normal);
    else
      toi = determine_TOI(0.0, (Real) 1.0, &ds, point, normal);
  
    // insert into the contacts set if the TOI is finite 
    if (toi < std::numeric_limits<Real>::max())
    {
      // insert the contacts
      local_contacts.push_back(create_contact(toi, a, b, point, normal));
      if (toi < earliest)
        earliest = std::min(toi + TOI_TOLERANCE/dt, (Real) 1.0);
    }

    // move onto the next vertex
    Q.pop_back();
  }
} 

/// Gets the velocity-expanded OBB for a BV 
BVPtr GeneralizedCCD::get_vel_exp_BV(CollisionGeometryPtr cg, BVPtr bv, const Vector3& lv, const Vector3& av)
{
  // verify that the map for the geometry has already been setup
  map<CollisionGeometryPtr, map<BVPtr, BVPtr> >::iterator vi;
  vi = _ve_BVs.find(cg);
  assert(vi != _ve_BVs.end());

  // see whether the velocity-expanded BV has already been calculated
  map<BVPtr, BVPtr>::const_iterator vj;
  #ifdef _OPENMP
  pthread_mutex_lock(&_ve_BVs_mutex);
  #endif
  vj = vi->second.find(bv);
  #ifdef _OPENMP
  pthread_mutex_unlock(&_ve_BVs_mutex);
  #endif
  if (vj != vi->second.end())
    return vj->second;

  // otherwise, calculate it
  OBBPtr obb = boost::dynamic_pointer_cast<OBB>(bv);
  if (LOGGING(LOG_BV) && obb)
  {
    FILE_LOG(LOG_BV) << "calculating velocity-expanded OBB for: " << obb << std::endl;
    FILE_LOG(LOG_BV) << "unexpanded OBB: " << *obb << std::endl;
  }
  BVPtr ve_bv = bv->calc_vel_exp_BV(cg, (Real) 1.0, lv, av);
  FILE_LOG(LOG_BV) << "new OBB: " << ve_bv << std::endl;

  // store the bounding volume
  #ifdef _OPENMP
  pthread_mutex_lock(&_ve_BVs_mutex);
  #endif
  vi->second[bv] = ve_bv;
  #ifdef _OPENMP
  pthread_mutex_unlock(&_ve_BVs_mutex);
  #endif

  return ve_bv;
}

/// Implements Base::load_from_xml()
void GeneralizedCCD::load_from_xml(XMLTreeConstPtr node, map<std::string, BasePtr>& id_map)
{
  map<std::string, BasePtr>::const_iterator id_iter;

  // verify that the node name is correct
  assert(strcasecmp(node->name.c_str(), "GeneralizedCCD") == 0);

  // call parent method
  CollisionDetection::load_from_xml(node, id_map);

  // get the maximum number of OBB expansions
  const XMLAttrib* _max_dexp_attr = node->get_attrib("max-dexp");
  if (_max_dexp_attr)
    this->_max_dexp = _max_dexp_attr->get_unsigned_value();

  // get the eps tolerance, if specified
  const XMLAttrib* eps_attr = node->get_attrib("eps-tolerance");
  if (eps_attr)
    this->eps_tolerance = eps_attr->get_real_value();
}

/// Implements Base::save_to_xml()
/**
 * \note neither the contact cache nor the pairs currently in collision are 
 *       saved
 */
void GeneralizedCCD::save_to_xml(XMLTreePtr node, list<BaseConstPtr>& shared_objects) const
{
  // call parent save_to_xml() method first
  CollisionDetection::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "GeneralizedCCD";

  // save the eps tolerance
  node->attribs.insert(XMLAttrib("eps-tolerance", eps_tolerance));
}

/****************************************************************************
 Methods for Drumwright-Shell algorithm begin 
****************************************************************************/

/// Creates a contact event given the bare-minimum info
Event GeneralizedCCD::create_contact(Real toi, CollisionGeometryPtr a, CollisionGeometryPtr b, const Vector3& point, const Vector3& normal)
{
  Event e;
  e.t = toi;
  e.event_type = Event::eContact;
  e.contact_point = point;
  e.contact_normal = normal;
  e.contact_geom1 = a;  
  e.contact_geom2 = b;  

  // check for valid normal here
  assert(std::fabs(e.contact_normal.norm() - (Real) 1.0) < NEAR_ZERO);

  // make the body first that comes first alphabetically
  if (LOGGING(LOG_COLDET))
  {
    SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
    SingleBodyPtr sb2 = e.contact_geom2->get_single_body();
    if (sb2->id < sb1->id)
    {
      std::swap(e.contact_geom1, e.contact_geom2);
      e.contact_normal = -e.contact_normal;
    }
  }

  return e;
}

/// Populates the DStruct
void GeneralizedCCD::populate_dstruct(DStruct* ds, CollisionGeometryPtr gb, CollisionGeometryPtr gs, const Vector3& bblv, const Vector3& bbav, const Vector3& bs_lvel, const Vector3& bs_avel, BVPtr s_BV)
{
  // save the BV
  ds->s_BV = s_BV;

  // save the geometries 
  ds->gb = gb;
  ds->gs = gs;

  // get the rigid bodies
  RigidBodyPtr bb = dynamic_pointer_cast<RigidBody>(gb->get_single_body());
  RigidBodyPtr bs = dynamic_pointer_cast<RigidBody>(gs->get_single_body());

  // get positions of the centers-of-mass of b and s
  const Vector3& xs = bs->get_position();
  const Vector3& xb = bb->get_position();

  // get velocities of the rigid bodies
  const Matrix3 omegas_hat = Matrix3::skew_symmetric(bs_avel);
  const Matrix3 omegab_hat = Matrix3::skew_symmetric(bbav);

  // determine the relative angular velocity (vector form)
  ds->romega = bbav - bs_avel;

  // determine magnitudes of linear and angular velocities
  ds->thetad = ds->romega.norm();
  ds->lvd = (bblv - bs_lvel).norm();

  // determine Rb and R;
  Matrix3 Rb, Rs;
  gb->get_transform().get_rotation(&Rb);
  gs->get_transform().get_rotation(&Rs);

  // precalculate constants
  /***********************************************************************
  1. We start from the following formula for calculating the position of
     a point in body S's frame:

     p_s = R_S^T ( R_B u + \bar{x}_B ) - R_S^T \bar{x}_S

  2. Upon differentiation w.r.t. time, we get:

     \dot{p}_b = \dot{R}_S^T ( R_B u + \bar{x}_B + 
                 R_S^T ( \dot{R}_B u + \dot{\bar{x}}_S ) -
                 R_S^T \dot{\bar{x}}_S - \dot{R}_S^T \bar{x}_S
  
  3. After rearranging some terms and using the identity 
     \dot{R} = \hat{\omega} R, we get the constants below.
  ***********************************************************************/
  ds->k1 = Rs.transpose_mult(bblv - bs_lvel - omegas_hat * (xb - xs));
  ds->k2 = Rs.transpose_mult(omegab_hat - omegas_hat)*Rb;

  // compute quaternions
  Quat qs(&Rs);
  Quat qb(&Rb);

  // setup quaternion for interpolation and quaternion derivative
  ds->q0 = Quat::conjugate(qs) * qb;

  // setup velocities of bs
  ds->bs_xd = bs->get_lvel();
  ds->bs_omega = bs->get_avel();
}

/// Implements DETERMINE-TOC()
/**
 * \param t0 the time before integration
 * \param tf the time after integration
 * \param dstruct structure passed to determine_TOI of precomputed data
 * \param pt the contact point (if any) in world coordinates
 * \param normal the normal (if any) in world coordinates
 * \return the time of impact [0, 1] such that t0 = 0 and t1 = 1; if the time 
 *         of impact is greater than 1.0, then the time of impact is not in the 
 *         interval [t0, tf] and should be discarded
 * \pre rigid body states are at time t0
 */
Real GeneralizedCCD::determine_TOI(Real t0, Real tf, const DStruct* ds, Vector3& pt, Vector3& normal) const
{
  const Real INF = std::numeric_limits<Real>::max();
  const unsigned X = 0, Y = 1, Z = 2;
  OBB O;
  Vector3 nalpha, nbeta, ngamma;

  // init bisection statistics
  unsigned nbisects = 0;

  // get the two geometries
  CollisionGeometryPtr gs = ds->gs;
  CollisionGeometryPtr gb = ds->gb;

  // get the primitive for gs 
  PrimitiveConstPtr gs_primitive = gs->get_geometry();

  // get the two rigid bodies
  RigidBodyPtr bs = dynamic_pointer_cast<RigidBody>(gs->get_single_body());
  RigidBodyPtr bb = dynamic_pointer_cast<RigidBody>(gb->get_single_body());

  // setup the quaternions used for interpolation
  Matrix3 Rs = gs->get_transform().get_rotation();
  Matrix3 Rb = gb->get_transform().get_rotation();
  Quat qs(&Rs);
  Quat qb(&Rb);

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI() entered" << endl;
  FILE_LOG(LOG_COLDET) << "  time t0: " << t0 << endl;
  FILE_LOG(LOG_COLDET) << "  time tf: " << tf << endl;

  // get p0 object (for simplicity)
  const Vector3& p0 = ds->p0;

  // get the BV for S
  BVPtr gs_BV = ds->s_BV; 

  // setup quaternions for interpolation
  const Quat& q0 = ds->q0; 

  // integrate q0 to dt to get the other quaternion for interpolation
  Quat qf = q0 + Quat::deriv(q0, ds->romega)*(tf-t0);
  qf.normalize();

  // (Euler) integrate the position of the point to dt; note: we do this
  // because the derivative is constant
  Vector3 pdot = ds->k1 + ds->k2*ds->u;
  Vector3 pf = p0 + pdot*(tf-t0);

  // setup the axes of the bounding box
  Vector3 p0pf = pf - p0;
  Real norm_p0pf = p0pf.norm();
  if (norm_p0pf < std::numeric_limits<Real>::epsilon())
  {
    // arbitrary bounding box
    nalpha = Vector3(1,0,0);
    nbeta = Vector3(0,1,0);
    ngamma = Vector3(0,0,1);
    O.R = Matrix3::identity();
  }
  else
  {
    // primary axis of bounding box aligned w/p0pf
    nalpha = p0pf/norm_p0pf;
    Vector3::determine_orthonormal_basis(nalpha, nbeta, ngamma);
    O.R.set_column(X, nalpha);
    O.R.set_column(Y, nbeta);
    O.R.set_column(Z, ngamma);
  }

  // determine whether minimum/maximum deviation is between two planes;
  // if so, determine the interpolation value that yields the minimum
  // and maximum deviation
  Vector3 normal1, normal2;
  Real rho1_max = -1.0, rho1_min = -1.0, rho2_max = -1.0, rho2_min = -1.0, rho3_max = -1.0, rho3_min = -1.0;
  Real max_d1 = -INF, max_d2 = -INF, max_d3 = -INF;
  Real min_d1 = INF, min_d2 = INF, min_d3 = INF;
  if (bound_u(ds->u, q0, qf, normal1, normal2))
  {
    if (nalpha.dot(normal1)*nalpha.dot(normal2) > 0)
    {
      max_d1 = calc_max_dev(ds->u, nalpha, q0, qf, rho1_max);
      min_d1 = calc_min_dev(ds->u, nalpha, q0, qf, rho1_min);

      // scale rho values to (tf-t0)
      rho1_max *= (tf-t0);
      rho1_min *= (tf-t0);
    }
    if (nbeta.dot(normal1)*nbeta.dot(normal2) > 0)
    {
      max_d2 = calc_max_dev(ds->u, nbeta, q0, qf, rho2_max);
      min_d2 = calc_min_dev(ds->u, nbeta, q0, qf, rho2_min);

      // scale rho values to (tf-t0)
      rho2_max *= (tf-t0);
      rho2_min *= (tf-t0);
    }
    if (ngamma.dot(normal1)*ngamma.dot(normal2) > 0)
    {
      max_d3 = calc_max_dev(ds->u, ngamma, q0, qf, rho3_max);
      min_d3 = calc_min_dev(ds->u, ngamma, q0, qf, rho3_min);

      // scale rho values to (tf-t0)
      rho3_max *= (tf-t0);
      rho3_min *= (tf-t0);
    }
  }

  // setup the queue
  typedef pair<Real, Real> RPair;
  priority_queue<RPair, vector<RPair>, std::greater<RPair> > Q;
  Q.push(make_pair(t0, tf));

  // process until error sufficiently low
  while (!Q.empty())
  {
    // get the element off of the top of the queue
    Real ta = Q.top().first;
    Real tb = Q.top().second;
    Q.pop();

    // setup delta t
    const Real dt = tb - ta;

    // bisect if angular velocity * dt > M_PI
    if (ds->thetad*dt > M_PI)
    {
      // if the linear velocity is near zero, no need to check more than 2*pi
      if (ds->lvd*dt < NEAR_ZERO)
      {
        const Real delta = M_PI/ds->thetad;
        Real tmid = ta + delta;
        tb = tmid + delta;
        assert(ta + NEAR_ZERO > t0);
        assert(tf + NEAR_ZERO > tb);
        Q.push(make_pair(ta, tmid));
        Q.push(make_pair(tmid, tb));
        nbisects++;
        continue;
      }
      else
      {
        Real tmid = (ta+tb)*0.5;
        Q.push(make_pair(ta, tmid));
        Q.push(make_pair(tmid, tb));
        nbisects++;
        continue;  
      }
    }

    // determine pa and pb
    Vector3 pa = p0 + pdot*ta;
    Vector3 pb = p0 + pdot*tb;

    FILE_LOG(LOG_COLDET) << " -- checking segment for time [" << ta << ", " << tb << "]" << endl;
    FILE_LOG(LOG_COLDET) << "  p(" << ta << ") = " << pa << "  p(" << tb << ") ~= " << pb << endl;

    // interpolation parameter ranges from [0,1]; 0 corresponds to t0, 1 corresponds
    // to tf.  Determine what ta, tb correspond to
    const Real sa = ta/(tf-t0);
    const Real sb = tb/(tf-t0);

    // init deviation maxima/minima
    Real dp_alpha = -INF, dp_beta = -INF, dp_gamma = -INF;
    Real dn_alpha = INF, dn_beta = INF, dn_gamma = INF;

    // see whether this interval contains a minimum/maximum 
    if (rho1_max >= ta && rho1_max <= tb) dp_alpha = max_d1;
    if (rho1_min >= ta && rho1_min <= tb) dn_alpha = min_d1;
    if (rho2_max >= ta && rho2_max <= tb) dp_beta = max_d2;
    if (rho2_min >= ta && rho2_min <= tb) dn_beta = min_d2;
    if (rho3_max >= ta && rho3_max <= tb) dp_gamma = max_d3;
    if (rho3_min >= ta && rho3_min <= tb) dn_gamma = min_d3;
    
    // calculate deviation at endpoints
    pair<Real, Real> deva = calc_deviations(ds->u, nalpha, q0, qf, sa, sb);
    pair<Real, Real> devb = calc_deviations(ds->u, nbeta, q0, qf, sa, sb);
    pair<Real, Real> devg = calc_deviations(ds->u, ngamma, q0, qf, sa, sb);
  
    // set deviation maxima/minima    
    dn_alpha = std::min(dn_alpha, deva.first);
    dn_beta = std::min(dn_beta, devb.first);
    dn_gamma = std::min(dn_gamma, devg.first);
    dp_alpha = std::max(dp_alpha, deva.second);
    dp_beta = std::max(dp_beta, devb.second);
    dp_gamma = std::max(dp_gamma, devg.second);

    // determine the half deviations in each direction
    Real len_pab = (pb - pa).norm();
    Real d_alpha = std::max(std::fabs(dp_alpha - dn_alpha), len_pab*(Real) 0.5);
    Real d_beta = std::fabs(dp_beta - dn_beta);
    Real d_gamma = std::fabs(dp_gamma - dn_gamma);

    // setup the bounding box
    Real nu = (d_beta + d_gamma + d_alpha)*2.0 - len_pab;
    assert(nu > -NEAR_ZERO);
    FILE_LOG(LOG_COLDET) << " -- dalpha: [" << dn_alpha << ", " << dp_alpha << "]" << endl;
    FILE_LOG(LOG_COLDET) << " -- dbeta: [" << dn_beta << ", " << dp_beta << "]" << endl;
    FILE_LOG(LOG_COLDET) << " -- dgamma: [" << dn_gamma << ", " << dp_gamma << "]" << endl;
    FILE_LOG(LOG_COLDET) << " -- nu contribution along alpha/beta: " << nu << endl;
    Vector3 cab = (pa + pb)*0.5;
    Vector3 mini = cab - nalpha*d_alpha - nbeta*d_beta - ngamma*d_gamma;
    Vector3 maxi = cab + nalpha*d_alpha + nbeta*d_beta + ngamma*d_gamma;
    O.center = (maxi+mini)*0.5;
    O.l = O.R.transpose_mult((maxi-mini)*0.5);
    O.l[0] = std::fabs(O.l[0]);
    O.l[1] = std::fabs(O.l[1]);
    O.l[2] = std::fabs(O.l[2]);

    // determine whether there is an intersection with the bounding box for S
    FILE_LOG(LOG_COLDET) << " -- nu: " << nu << endl;
    FILE_LOG(LOG_COLDET) << " -- checking BV intersection of: " << endl;
    FILE_LOG(LOG_COLDET) << O << " and " << endl << gs_BV;
    if (LOGGING(LOG_COLDET))
    {
      OBBPtr gs_OBB = dynamic_pointer_cast<OBB>(gs_BV);
      if (gs_OBB)
        FILE_LOG(LOG_COLDET) << *gs_OBB << endl;
    }

    // NOTE: O and gs_BV are both in S's frame
    if (!BV::intersects(&O, gs_BV.get()))
    {
      FILE_LOG(LOG_COLDET) << "   -- no intersection; continuing looping..." << endl;
      continue;
    }

    // ************************************************************************
    // there is a bounding box intersection; bisect or intersect with triangles 
    // ************************************************************************

    if (nu > eps_tolerance)
    {
      FILE_LOG(LOG_COLDET) << "   -- intersection detected; trajectory segment must be bisected" << endl;

      // add two elements to the queue
      Real ti = (ta+tb)*0.5;
      Q.push(make_pair(ta, ti));
      Q.push(make_pair(ti, tb));
      nbisects++;
    }
    else
    {
      FILE_LOG(LOG_COLDET) << "   -- intersection detected and nu less than tolerance" << endl;

      // intersect the line segment with the geometry
      Real t;
      if (gs_primitive->intersect_seg(gs_BV, LineSeg3(pa, pb), t, pt, normal))
      {
        FILE_LOG(LOG_COLDET) << "  intersection detected!  time of impact: " << t << " (true: " << (ta + (tb-ta)*t) << ")" << endl;

        // transform time to account for ta and tb
        t = ta + (tb-ta)*t;

        FILE_LOG(LOG_COLDET) << "     -- point intersection (untransformed): " << pt << endl;
        FILE_LOG(LOG_COLDET) << "     -- normal (untransformed): " << normal << endl;

        // since all calculations are in body S's frame; we need to integrate
        // body S's state to time t to transform the contact point and normal
        // first integrate orientation
        Quat qf = qs + Quat::deriv(qs, ds->bs_omega)*t;
        qf.normalize();

        // now integrate position of the c.o.m. 
        const Vector3& xdot = ds->bs_xd;  
        Vector3 x = gs->get_transform().get_translation();
        x += xdot*t;

        // transform contact point and normal to global coords
        pt = qf * pt + x;
        normal = qf * normal;

        // look for degenerate normal
        if (std::fabs(normal.norm() - (Real) 1.0) > NEAR_ZERO)
        {
          FILE_LOG(LOG_COLDET) << "    -- degenerate normal detected! (" << normal << "); not reporting intersection" << endl;
          continue;
        }
        
        
        FILE_LOG(LOG_COLDET) << "     -- point intersection: " << pt << endl;
        FILE_LOG(LOG_COLDET) << "     -- normal: " << normal << endl;
        FILE_LOG(LOG_COLDET) << "# of bisections: " << nbisects << endl;
        FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI() exited" << endl;

        return t;
      }
      else
      {
        FILE_LOG(LOG_COLDET) << "     -- intersection does not occur in interval [t0, tf]" << endl;
        FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI() exited" << endl;
      }
    }
  }

  FILE_LOG(LOG_COLDET) << "# of bisections: " << nbisects << endl;
  FILE_LOG(LOG_COLDET) << "  no impact detected" << endl;
  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI() exited" << endl;

  // still here?  no contact...
  return INF; 
}

/// Determines the two planes that bound a u vector between two rotations
bool GeneralizedCCD::bound_u(const Vector3& u, const Quat& q0, const Quat& qf, Vector3& normal1, Vector3& normal2)
{
  // compute the rotated u
  Vector3 u0 = q0*u;
  Vector3 uf = qf*u;

  // determine a vector perpendicular to both
  Vector3 perp = Vector3::cross(u0, uf);
  Real perp_norm = perp.norm();
  if (perp_norm < NEAR_ZERO)
    return false;
  else
    perp /= perp_norm;

  // determine the two planes
  normal1 = Vector3::normalize(Vector3::cross(u0, perp));
  normal2 = Vector3::normalize(Vector3::cross(uf, perp));

  // make sure that both vectors are on the same side of the plane
  Real dot1 = u0.dot(normal2);
  Real dot2 = uf.dot(normal1);
  if (dot1*dot2 < 0)
    normal2 = -normal2;
 
  return true;
}

/// Function for computing the maximum deviation in a given direction
Real GeneralizedCCD::calc_deviation(Real alpha, void* params)
{
  const DeviationCalc& data = *((const DeviationCalc*) params);

  // linearly interpolate the quaternions
  Quat q = Quat::lerp(data.q1, data.q2, alpha);

  // do the arithmetic
  return Vector3::dot(data.d, q*data.u);
}

/// Computes the minimum deviation over interval [0,1] using bracketing and Brent's method
/**
 * \param u vector from center-of-mass to point (in local frame)
 * \param d the direction vector
 * \param q1 the orientation at t=0
 * \param q2 the orientation at t=1
 * \param t on return, contains the t at which the deviation is minimized
 * \return the minimum deviation
 */
Real GeneralizedCCD::calc_min_dev(const Vector3& u, const Vector3& d, const Quat& q1, const Quat& q2, Real& t)
{
  const Real TOL = std::sqrt(std::sqrt(std::numeric_limits<Real>::epsilon()));

  Real ft;

  // setup data for deviation calculations
  DeviationCalc dc;
  dc.u = u;
  dc.q1 = q1;
  dc.q2 = q2;
  dc.d = d;

  // find suitable starting value for brent's method
  Real f0 = calc_deviation(0.0, &dc);
  Real f1 = calc_deviation(1.0, &dc);
  t = (f0 < f1) ? 0.0 : 1.0;
  
  // call Brent's method -- note: tolerance is a bit high
  if (!Optimization::brent(0.0, 1.0, t, ft, &calc_deviation, TOL, &dc))
    cerr << "GeneralizedCCD::calc_min_dev() - brent's method failed!" << endl;

  return ft;
}

/// Computes the maximum and minimum deviation in a given direction for two points
pair<Real, Real> GeneralizedCCD::calc_deviations(const Vector3& u, const Vector3& d, const Quat& q1, const Quat& q2, Real t1, Real t2)
{
  // compute the deviation in the two directions
  DeviationCalc dc;
  dc.u = u;
  dc.d = d;
  dc.q1 = q1;
  dc.q2 = q2;
  Real dev1 = calc_deviation(t1, &dc);
  Real dev2 = calc_deviation(t2, &dc);
  pair<Real, Real> dev(dev1, dev2);
  if (dev1 > dev2)
    std::swap(dev.first, dev.second);
  return dev;
}

/// Computes the maximum deviation over interval [0,1] using bracketing and Brent's method
/**
 * \param u vector from center-of-mass to point (in local frame)
 * \param d the direction vector
 * \param q1 the orientation at t=0
 * \param q2 the orientation at t=1
 * \param t on return, contains the t at which the deviation is maximized [0,1]
 * \return the maximum deviation
 */
Real GeneralizedCCD::calc_max_dev(const Vector3& u, const Vector3& d, const Quat& q1, const Quat& q2, Real& t)
{
  const Real TOL = 1e-2;

  Real ft;

  // setup data for deviation calculations
  DeviationCalc dc;
  dc.u = u;
  dc.q1 = q1;
  dc.q2 = q2;

  // set d in the deviation data because brent's method is a minimizer and we
  // want to maximize
  dc.d = -d;

  // find suitable starting value for brent's method
  Real f0 = calc_deviation(0.0, &dc);
  Real f1 = calc_deviation(1.0, &dc);
  t = (f0 < f1) ? 0.0 : 1.0;
  
  // call Brent's method -- note: tolerance is a bit high
  if (!Optimization::brent(0.0, 1.0, t, ft, &calc_deviation, TOL, &dc))
    cerr << "GeneralizedCCD::calc_max_dev() - brent's method failed!" << endl;

  return -ft;
}

/****************************************************************************
 Methods for Drumwright-Shell algorithm end 
****************************************************************************/

/****************************************************************************
 Methods for broad phase begin 
****************************************************************************/
void GeneralizedCCD::broad_phase(const map<SingleBodyPtr, pair<Vector3, Vector3> >& vel_map, vector<pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
{
  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::broad_phase() entered" << std::endl;

  // clear the vector of pairs to check
  to_check.clear();

  // sort the AABBs
  sort_AABBs(vel_map);

  // store how many overlaps we have for pairs
  map<sorted_pair<CollisionGeometryPtr>, unsigned> overlaps;

  // set of active bounds
  set<CollisionGeometryPtr> active_bounds;

  // scan through the x-bounds
  for (unsigned i=0; i< _x_bounds.size(); i++)
  {
    // eliminate from the active bounds if at the end of a bound
    if (_x_bounds[i].second.end)
    {
      assert(active_bounds.find(_x_bounds[i].second.geom) != active_bounds.end());
      active_bounds.erase(_x_bounds[i].second.geom);
    }
    else
    {
      // at the start of a bound
      BOOST_FOREACH(CollisionGeometryPtr cg, active_bounds)
        overlaps[make_sorted_pair(cg, _x_bounds[i].second.geom)]++;

      // add the geometry to the active set
      active_bounds.insert(_x_bounds[i].second.geom);
    }
  }

  // scan through the y-bounds
  for (unsigned i=0; i< _y_bounds.size(); i++)
  {
    // eliminate from the active bounds if at the end of a bound
    if (_y_bounds[i].second.end)
    {
      assert(active_bounds.find(_y_bounds[i].second.geom) != active_bounds.end());
      active_bounds.erase(_y_bounds[i].second.geom);
    }
    else
    {
      // at the start of a bound
      BOOST_FOREACH(CollisionGeometryPtr cg, active_bounds)
        overlaps[make_sorted_pair(cg, _y_bounds[i].second.geom)]++;

      // add the geometry to the active set
      active_bounds.insert(_y_bounds[i].second.geom);
    }
  }

  // scan through the z-bounds
  for (unsigned i=0; i< _z_bounds.size(); i++)
  {
    // eliminate from the active bounds if at the end of a bound
    if (_z_bounds[i].second.end)
    {
      assert(active_bounds.find(_z_bounds[i].second.geom) != active_bounds.end());
      active_bounds.erase(_z_bounds[i].second.geom);
    }
    else
    {
      // at the start of a bound
      BOOST_FOREACH(CollisionGeometryPtr cg, active_bounds)
        overlaps[make_sorted_pair(cg, _z_bounds[i].second.geom)]++;

      // add the geometry to the active set
      active_bounds.insert(_z_bounds[i].second.geom);
    }
  }

  // now setup pairs to check
  for (map<sorted_pair<CollisionGeometryPtr>, unsigned>::const_iterator i = overlaps.begin(); i != overlaps.end(); i++)
  {
    FILE_LOG(LOG_COLDET) << i->second << " overlaps between " << i->first.first << " (" << i->first.first->get_single_body()->id << ") and " << i->first.second << " (" << i->first.second->get_single_body()->id << ")" << std::endl;

    // only consider the pair if they overlap in all three dimensions
    if (i->second < 3)
      continue;

    // if the pair is disabled, continue looping
    if (this->disabled_pairs.find(i->first) != this->disabled_pairs.end())
      continue;

    // get the rigid bodies corresponding to the geometries
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(i->first.first->get_single_body());
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(i->first.second->get_single_body());

    // don't check pairs from the same rigid body
    if (rb1 == rb2)
      continue;

    // if both rigid bodies are disabled, don't check
    if (!rb1->is_enabled() && !rb2->is_enabled())
      continue;

    // if we're here, we have a candidate for the narrow phase
    to_check.push_back(make_pair(i->first.first, i->first.second));
    FILE_LOG(LOG_COLDET) << "  ... checking pair" << std::endl;
  }
  
  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::broad_phase() exited" << std::endl;
}

void GeneralizedCCD::sort_AABBs(const map<SingleBodyPtr, pair<Vector3, Vector3> >& vel_map)
{
  // if a geometry was added or removed, rebuild the vector of bounds
  // and copy it to x, y, z dimensions
  if (_rebuild_bounds_vecs)
  {
    build_bv_vector(vel_map, _x_bounds);
    _y_bounds = _x_bounds;
    _z_bounds = _x_bounds;
  }

  // update bounds vectors
  update_bounds_vector(_x_bounds, vel_map, eXAxis);
  update_bounds_vector(_y_bounds, vel_map, eYAxis);
  update_bounds_vector(_z_bounds, vel_map, eZAxis);
  
  // if geometry was added or removed, do standard sorts of vectors
  if (_rebuild_bounds_vecs)
  {
    std::sort(_x_bounds.begin(), _x_bounds.end());
    std::sort(_y_bounds.begin(), _y_bounds.end());
    std::sort(_z_bounds.begin(), _z_bounds.end());

    // now indicate that bounds vectors have been (re)built
    _rebuild_bounds_vecs = false;
  }
  else
  {
    // bounds are nearly sorted; do insertion sort
    insertion_sort(_x_bounds.begin(), _x_bounds.end());
    insertion_sort(_y_bounds.begin(), _y_bounds.end());
    insertion_sort(_z_bounds.begin(), _z_bounds.end());
  }
}

void GeneralizedCCD::update_bounds_vector(vector<pair<Real, BoundsStruct> >& bounds, const map<SingleBodyPtr, pair<Vector3, Vector3> >& vel_map, AxisType axis)
{
  const unsigned X = 0, Y = 1, Z = 2;

  FILE_LOG(LOG_COLDET) << " -- update_bounds_vector() entered (axis=" << axis << ")" << std::endl;

  // iterate over bounds vector
  for (unsigned i=0; i< bounds.size(); i++)
  {
    // get the bounding volume, collision geometry, and rigid body
    BVPtr bv = bounds[i].second.bv;
    CollisionGeometryPtr geom = bounds[i].second.geom;
    RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(geom->get_single_body());

    // update the velocities
    const Vector3& xd = vel_map.find(rb)->second.first;
    const Vector3& omega = vel_map.find(rb)->second.second;
    bounds[i].second.xd = xd; 
    bounds[i].second.omega = omega;

    // get the expanded bounding volume
    BVPtr bv_exp = get_vel_exp_BV(geom, bv, xd, omega);

    // get the transform for the collision geometry
    const Matrix4& T = geom->get_transform();

    // get the bound for the bounding volume
    Vector3 bound = (bounds[i].second.end) ? bv_exp->get_upper_bounds(T) : bv_exp->get_lower_bounds(T);
    FILE_LOG(LOG_COLDET) << "  updating collision geometry: " << geom << "  rigid body: " << geom->get_single_body()->id << std::endl;

    // update the bounds for the given axis
    switch (axis)
    {
      case eXAxis:
        bounds[i].first = bound[X];
        break;
 
      case eYAxis:
        bounds[i].first = bound[Y];
        break;

      case eZAxis:
        bounds[i].first = bound[Z];
        break;

      default:
        assert(false);
    }

    if (bounds[i].second.end)
      FILE_LOG(LOG_COLDET) << "    upper bound: " << bounds[i].first << std::endl; 
    if (!bounds[i].second.end)  
      FILE_LOG(LOG_COLDET) << "    lower bound: " << bounds[i].first << std::endl;
  }

  FILE_LOG(LOG_COLDET) << " -- update_bounds_vector() exited" << std::endl;
}

void GeneralizedCCD::build_bv_vector(const map<SingleBodyPtr, pair<Vector3, Vector3> >& vel_map, vector<pair<Real, BoundsStruct> >& bounds)
{
  const Real INF = std::numeric_limits<Real>::max();

  // clear the vector
  bounds.clear();

  // iterate over all collision geometries
  for (set<CollisionGeometryPtr>::const_iterator i = _geoms.begin(); i != _geoms.end(); i++)
  {
    // if the geometry is disabled, skip the geometry
    if (this->disabled.find(*i) != this->disabled.end())
      continue;

    // get the rigid body and primitive for the geometry
    RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>((*i)->get_single_body());
    PrimitivePtr p = (*i)->get_geometry();

    // get the top-level BV for the geometry
    BVPtr bv = p->get_BVH_root();

    // get the velocities for the rigid body
    assert(vel_map.find(rb) != vel_map.end());
    const Vector3& xd = vel_map.find(rb)->second.first;
    const Vector3& omega = vel_map.find(rb)->second.second;

    // setup the bounds structure
    BoundsStruct bs;
    bs.end = false;
    bs.geom = *i;
    bs.xd = xd;
    bs.omega = omega;
    bs.bv = bv;

    // add the lower bound
    bounds.push_back(make_pair(-INF, bs));

    // modify the bounds structure to indicate the end bound
    bs.end = true;
    bounds.push_back(make_pair(INF, bs));
  }
}

/****************************************************************************
 Methods for broad phase end 
****************************************************************************/

/****************************************************************************
 Methods for static geometry intersection testing begin 
****************************************************************************/

/// Determines whether there is a collision at the current position and orientation of the bodies
/**
 * \note the epsilon parameter is ignored
 */
bool GeneralizedCCD::is_collision(Real epsilon)
{
  // clear the set of colliding pairs and list of colliding triangles
  colliding_pairs.clear();
  colliding_tris.clear();

  // iterate over geometries 
  for (std::set<CollisionGeometryPtr>::const_iterator i = _geoms.begin(); i != _geoms.end(); i++)
  {
    // get the first geometry and its primitive 
    CollisionGeometryPtr g1 = *i;
    PrimitivePtr g1_primitive = g1->get_geometry();

    // get the transform and the inverse transform for this geometry
    const Matrix4& wTg1 = g1->get_transform();
    Matrix4 g1Tw = Matrix4::inverse_transform(wTg1);

    // loop through all other geometries
    std::set<CollisionGeometryPtr>::const_iterator j = i;
    j++;
    for (; j != _geoms.end(); j++)
    {
      // get the second geometry and its primitive 
      CollisionGeometryPtr g2 = *j;
      PrimitivePtr g2_primitive = g2->get_geometry();

      // see whether to check
      if (!is_checked(g1, g2))
        continue; 

      // get the two BV trees
      BVPtr bv1 = g1_primitive->get_BVH_root();
      BVPtr bv2 = g2_primitive->get_BVH_root();

      // get the transform for g2 and its inverse
      const Matrix4& wTg2 = g2->get_transform(); 

      // if intersects, add to colliding pairs
      if (intersect_BV_trees(bv1, bv2, g1Tw * wTg2, g1, g2))
        colliding_pairs.insert(make_sorted_pair(g1, g2));
    } 
  }

  return !colliding_pairs.empty();
}

/// Intersects two BV trees; returns <b>true</b> if one (or more) pair of the underlying triangles intersects
bool GeneralizedCCD::intersect_BV_trees(BVPtr a, BVPtr b, const Matrix4& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b) 
{
  std::queue<tuple<BVPtr, BVPtr, bool> > q;

  // get address of the last colliding triangle pair on the queue
  CollidingTriPair* last = (colliding_tris.empty()) ? NULL : &colliding_tris.back();

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::intersect_BV_trees() entered" << endl;

  // intersect the BVs at the top level
  if (!BV::intersects(a, b, aTb))
  {
    FILE_LOG(LOG_COLDET) << "  no intersection at top-level BVs" << endl;
    FILE_LOG(LOG_COLDET) << "GeneralizedCCD::intersect_BV_trees() exited" << endl;

    return false;
  }

  // add a and b to the queue
  q.push(make_tuple(a, b, false));

  // drill down alternatingly until both trees are exhausted
  while (!q.empty())
  {
    // get the two nodes off of the front of the queue
    BVPtr bv1 = q.front().get<0>();
    BVPtr bv2 = q.front().get<1>();
    bool rev = q.front().get<2>();

    // no longer need the elm on the front of the queue
    q.pop();

    // check for bv1 and bv2 both leafs
    if (bv1->is_leaf() && bv2->is_leaf())
    {
      if (!rev)
        intersect_BV_leafs(bv1, bv2, aTb, geom_a, geom_b, std::back_inserter(colliding_tris));
      else
        intersect_BV_leafs(bv2, bv1, aTb, geom_a, geom_b, std::back_inserter(colliding_tris));

      // see whether we want to exit early
      if (mode == eFirstContact && !colliding_tris.empty() && last != &colliding_tris.back())
        return true;
    }

    // drill down through bv2, if possible
    if (bv2->is_leaf())
    {
      // check the children of o1
      BOOST_FOREACH(BVPtr child, bv1->children)
        if ((!rev && BV::intersects(child, bv2, aTb)) || (rev && BV::intersects(bv2, child, aTb)))
          q.push(make_tuple(child, bv2, rev));
    }
    else
      BOOST_FOREACH(BVPtr child, bv2->children)
      {
        if ((!rev && BV::intersects(bv1, child, aTb)) || (rev && BV::intersects(child, bv1, aTb)))
          q.push(make_tuple(child, bv1, !rev));
      }
  }

  // see whether we have an intersection
  if (!colliding_tris.empty() && last != &colliding_tris.back())
    return true;

  FILE_LOG(LOG_COLDET) << "  -- all intersection checks passed; no intersection" << endl;
  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::intersect_BV_trees() exited" << endl;

  return false;
} 

/****************************************************************************
 Methods for static geometry intersection testing end 
****************************************************************************/

