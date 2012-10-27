/***************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <cmath>
#include <algorithm>
#include <vector>
#include <queue>
#include <map>
#include <fstream>

#ifdef USE_OSG
#include <osg/Geode>
#include <osg/Geometry>
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/MatrixTransform>
#include <osg/Material>
#endif

#include <Moby/Constants.h>
#include <Moby/CompGeom.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/Log.h>
#include <Moby/AAngle.h>
#include <Moby/Event.h>

using namespace Moby;
using std::list;
using std::vector;
using std::map;
using std::multimap;
using std::set;
using std::endl;
using boost::shared_ptr;

/// Creates an empty event 
Event::Event()
{
  t_true = (Real) -1.0;
  event_type = eNone;
  limit_dof = std::numeric_limits<unsigned>::max();
  limit_epsilon = (Real) 0.0;
  limit_upper = false;
  limit_impulse = (Real) 0.0;
  contact_normal = ZEROS_3;
  contact_impulse = ZEROS_3;
  contact_point = ZEROS_3;
  contact_mu_coulomb = (Real) 0.0;
  contact_mu_viscous = (Real) 0.0;
  contact_epsilon = (Real) 0.0;
  contact_NK = 4;
}

Event& Event::operator=(const Event& e)
{
  t_true = e.t_true;
  t = e.t;
  event_type = e.event_type;
  limit_epsilon = e.limit_epsilon;
  limit_dof = e.limit_dof;
  limit_upper = e.limit_upper;
  limit_impulse = e.limit_impulse;
  limit_joint = e.limit_joint;
  contact_normal = e.contact_normal;
  contact_geom1 = e.contact_geom1;
  contact_geom2 = e.contact_geom2;
  contact_point = e.contact_point;
  contact_impulse = e.contact_impulse;
  contact_mu_coulomb = e.contact_mu_coulomb;
  contact_mu_viscous = e.contact_mu_viscous;
  contact_epsilon = e.contact_epsilon;
  contact_NK = e.contact_NK;
  contact_tan1 = e.contact_tan1;
  contact_tan2 = e.contact_tan2;
  constraint_nimpulse.copy_from(e.constraint_nimpulse);
  constraint_fimpulse.copy_from(e.constraint_fimpulse);
  constraint_joint = e.constraint_joint;

  return *this;
}

/// Sets the contact parameters for this event
void Event::set_contact_parameters(const ContactParameters& cparams)
{
  contact_mu_coulomb = cparams.mu_coulomb;
  contact_mu_viscous = cparams.mu_viscous;
  contact_epsilon = cparams.epsilon;
  contact_NK = cparams.NK;
}

/// Computes the velocity of this event
/**
 * Positive velocity indicates separation, negative velocity indicates
 * impact, zero velocity indicates rest.
 */
Real Event::calc_event_vel() const
{
  if (event_type == eContact)
  {
    assert(contact_geom1 && contact_geom2);
    SingleBodyPtr sb1 = contact_geom1->get_single_body();
    SingleBodyPtr sb2 = contact_geom2->get_single_body();
    assert(sb1 && sb2);
    return sb1->calc_point_vel(contact_point, contact_normal) - sb2->calc_point_vel(contact_point, contact_normal);
  }
  else if (event_type == eLimit)
  {
    Real qd = limit_joint->qd[limit_dof];
    return (limit_upper) ? -qd : qd;
  }
  else
    assert(false);
}  

/// Sends the event to the specified stream
std::ostream& Moby::operator<<(std::ostream& o, const Event& e)
{
  o << "TOI: " << e.t << std::endl;

  switch (e.event_type)
  {
    case Event::eNone:
      o << "(event type: none)" << std::endl;
      return o;

    case Event::eLimit:
      o << "(event type: joint limit)" << std::endl;
      break;

    case Event::eContact:
      o << "(event type: contact)" << std::endl;
      break;

    case Event::eConstraint:
      o << "(event type: constraint)" << std::endl;
  }

  if (e.event_type == Event::eContact)
  {
    if (e.contact_geom1)
    {
      SingleBodyPtr sb1(e.contact_geom1->get_single_body());
      if (sb1)
      {
        o << "body1: " << sb1->id << std::endl;
      }
      else
        o << "body1: (undefined)" << std::endl;
    }
    else
      o << "geom1: (undefined)" << std::endl;
  
    if (e.contact_geom2)
    {
      SingleBodyPtr sb2(e.contact_geom2->get_single_body());
      if (sb2)
      {
        o << "body2: " << sb2->id << std::endl;
      }    
      else
        o << "body2: (undefined)" << std::endl;
     }
    else
      o << "geom2: (undefined)" << std::endl;

    o << "contact point: " << e.contact_point << std::endl;
    o << "normal: " << e.contact_normal << std::endl;

    // determine the relative normal velocity at the contact point
    // get the rigid bodies of the contact
    if (e.contact_geom1 && e.contact_geom2)
    {
      SingleBodyPtr sb1(e.contact_geom1->get_single_body());
      SingleBodyPtr sb2(e.contact_geom2->get_single_body());
      if (sb1 && sb2)
      {
        Real cp1 = sb1->calc_point_vel(e.contact_point, e.contact_normal);
        Real cp2 = sb2->calc_point_vel(e.contact_point, e.contact_normal);
        Real rvel = cp1 - cp2; 
        o << "relative normal velocity: " << rvel << std::endl;
      }
    }
  }

  return o;
}

#ifdef USE_OSG
/// Copies this matrix to an OpenSceneGraph Matrixd object
static void to_osg_matrix(const Matrix4& src, osg::Matrixd& tgt)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  for (unsigned i=X; i<= W; i++)
    for (unsigned j=X; j<= Z; j++)
      tgt(j,i) = src(i,j);

  // set constant values of the matrix
  tgt(X,W) = tgt(Y,W) = tgt(Z,W) = (Real) 0.0;
  tgt(W,W) = (Real) 1.0;
}

/// Makes a contact visualizable
osg::Node* Event::to_visualization_data() const
{
  const float CONE_HEIGHT = .2f;
  const float CONE_RADIUS = .2f;
  const unsigned X = 0, Y = 1, Z = 2;

  // setup the transformation matrix for the cone
  Vector3 x_axis, z_axis;
  Vector3::determine_orthonormal_basis(contact_normal, x_axis, z_axis);
  Matrix3 R;
  R.set_column(X, x_axis);
  R.set_column(Y, contact_normal);
  R.set_column(Z, -z_axis);
  Vector3 x = contact_point + contact_normal;
  Matrix4 T(&R, &x);

  // setup the transform node for the cone
  osg::Matrixd m;
  to_osg_matrix(T, m);
  osg::MatrixTransform* transform = new osg::MatrixTransform;
  transform->setMatrix(m);

  // create the new color
  osg::Material* mat = new osg::Material;
  const float RED = (float) rand() / RAND_MAX;
  const float GREEN = (float) rand() / RAND_MAX;
  const float BLUE = (float) rand() / RAND_MAX;
  mat->setColorMode(osg::Material::DIFFUSE);
  mat->setDiffuse(osg::Material::FRONT, osg::Vec4(RED, GREEN, BLUE, 1.0f));
  transform->getOrCreateStateSet()->setAttribute(mat);

  // create the line
  osg::Geometry* linegeom = new osg::Geometry;
  osg::Vec3Array* varray = new osg::Vec3Array;
  linegeom->setVertexArray(varray);  
  varray->push_back(osg::Vec3((float) contact_point[X], (float) contact_point[Y], (float) contact_point[Z]));
  varray->push_back(osg::Vec3((float) contact_point[X] + (float) contact_normal[X], (float) contact_point[Y] + (float) contact_normal[Y], (float) contact_point[Z] + (float) contact_normal[Z]));
  osg::Geode* geode = new osg::Geode;
  geode->addDrawable(linegeom);

  // create the cone
  osg::Cone* cone = new osg::Cone;
  cone->setRadius(CONE_RADIUS);
  cone->setHeight(CONE_HEIGHT);
  geode->addDrawable(new osg::ShapeDrawable(cone));

  // add the geode
  transform->addChild(geode);

  return transform;
}
#endif

/// Given a vector of events, determines all of the sets of connected events
/**
 * A set of connected events is the set of all events such that, for a
 * given event A in the set, there exists another event B for which A
 * and B share at least one rigid body.  
 * \param events the list of events
 * \param groups the islands of connected events on return
 */
void Event::determine_connected_events(const vector<Event>& events, list<list<Event*> >& groups)
{
  FILE_LOG(LOG_CONTACT) << "Event::determine_connected_contacts() entered" << std::endl;

  // clear the groups
  groups.clear();

  // copy the list of events -- only ones with geometry
  list<Event*> events_copy;
  BOOST_FOREACH(const Event& e, events)
    if (e.event_type != Event::eNone)
      events_copy.push_back((Event*) &e);
  
  // The way that we'll determine the event islands is to treat each rigid
  // body present in the events as a node in a graph; nodes will be connected
  // to other nodes if (a) they are both present in event or (b) they are
  // part of the same articulated body.  Nodes will not be created for disabled
  // bodies.
  set<SingleBodyPtr> nodes;
  multimap<SingleBodyPtr, SingleBodyPtr> edges;
  typedef multimap<SingleBodyPtr, SingleBodyPtr>::const_iterator EdgeIter;

  // get all single bodies present in the events
  for (list<Event*>::const_iterator i = events_copy.begin(); i != events_copy.end(); i++)
  {
    if ((*i)->event_type == Event::eContact)
    {
      SingleBodyPtr sb1((*i)->contact_geom1->get_single_body());
      SingleBodyPtr sb2((*i)->contact_geom2->get_single_body());
      if (sb1->is_enabled())
        nodes.insert(sb1);
      if (sb2->is_enabled())
        nodes.insert(sb2);
      if (sb1->is_enabled() && sb2->is_enabled())
      {
        edges.insert(std::make_pair(sb1, sb2));
        edges.insert(std::make_pair(sb2, sb1));
      }
    }
    else if ((*i)->event_type == Event::eLimit)
    {
      RigidBodyPtr inboard = (*i)->limit_joint->get_inboard_link();
      RigidBodyPtr outboard = (*i)->limit_joint->get_outboard_link();
      nodes.insert(inboard);
      nodes.insert(outboard);
    }
    else if ((*i)->event_type == Event::eConstraint)
    {
      RigidBodyPtr inboard = (*i)->constraint_joint->get_inboard_link();
      RigidBodyPtr outboard = (*i)->constraint_joint->get_outboard_link();
      nodes.insert(inboard);
      nodes.insert(outboard);
    }
    else 
      assert(false);
  }

  FILE_LOG(LOG_CONTACT) << " -- single bodies in events:" << std::endl;
  if (LOGGING(LOG_CONTACT))
    for (set<SingleBodyPtr>::const_iterator i = nodes.begin(); i != nodes.end(); i++)
      FILE_LOG(LOG_CONTACT) << "    " << (*i)->id << std::endl;
  FILE_LOG(LOG_CONTACT) << std::endl;

  // add connections between articulated rigid bodies -- NOTE: don't process
  // articulated bodies twice!
  set<ArticulatedBodyPtr> ab_processed;
  BOOST_FOREACH(SingleBodyPtr sb, nodes)
  {
    // if the body is not part of an articulated body, skip it
    ArticulatedBodyPtr abody = sb->get_articulated_body();
    if (!abody)
      continue;

    // see whether it has already been processed
    if (ab_processed.find(abody) != ab_processed.end())
      continue;

    // indicate that the articulated body will now have been processed
    ab_processed.insert(abody);

    // get all links in the articulated body
    const vector<RigidBodyPtr>& links = abody->get_links();

    // add edges between all pairs for which there are links
    vector<RigidBodyPtr>::const_iterator j, k;
    for (j = links.begin(); j != links.end(); j++)
    {
      // no sense iterating over all other links if link pointed to by j is
      // not a node
      if (nodes.find(*j) == nodes.end())
        continue;

      // iterate over all other nodes
      k = j;
      for (k++; k != links.end(); k++)
        if (nodes.find(*k) != nodes.end())
        {
          edges.insert(std::make_pair(*j, *k));
          edges.insert(std::make_pair(*k, *j));
        }
    }      
  }

  // Now, we'll remove nodes from the set until there are no more nodes.
  // For each removed node, we'll get add all events that contain the single 
  // body to the group; all neighboring nodes will then be processed.
  while (!nodes.empty())
  {
    // get the node from the front
    SingleBodyPtr node = *nodes.begin();

    // add a list to the contact groups
    groups.push_back(list<Event*>());
    FILE_LOG(LOG_CONTACT) << " -- events in group: " << std::endl;

    // create a node queue, with this node added
    std::queue<SingleBodyPtr> node_q;
    node_q.push(node);

    // loop until the queue is empty
    while (!node_q.empty())
    {
      // get the node off of the front of the node queue
      node = node_q.front();
      node_q.pop();

      // erase the node from the set of nodes
      nodes.erase(node);

      // add all neighbors of the node that have not been processed already 
      // to the node queue
      std::pair<EdgeIter, EdgeIter> neighbors = edges.equal_range(node);
      for (EdgeIter i = neighbors.first; i != neighbors.second; i++)
        if (nodes.find(i->second) != nodes.end())
          node_q.push(i->second);

      // loop through all remaining events
      for (list<Event*>::iterator i = events_copy.begin(); i != events_copy.end(); )
      {
        if ((*i)->event_type == Event::eContact)
        {
          SingleBodyPtr sb1((*i)->contact_geom1->get_single_body());
          SingleBodyPtr sb2((*i)->contact_geom2->get_single_body());

          // see whether one of the bodies is equal to the node
          if (sb1 == node || sb2 == node)
          {
            groups.back().push_back(*i);
            i = events_copy.erase(i);
            continue;
          }
          else
            i++;
        }
        else if ((*i)->event_type == Event::eLimit)
        {
          RigidBodyPtr inboard = (*i)->limit_joint->get_inboard_link();
          RigidBodyPtr outboard = (*i)->limit_joint->get_outboard_link();
          if (inboard == node || outboard == node)
          {
            groups.back().push_back(*i);
            i = events_copy.erase(i);
            continue;
          }
          else
            i++;
        }
        else if ((*i)->event_type == Event::eConstraint)
        {
          RigidBodyPtr inboard = (*i)->constraint_joint->get_inboard_link();
          RigidBodyPtr outboard = (*i)->constraint_joint->get_outboard_link();
          if (inboard == node || outboard == node)
          {
            groups.back().push_back(*i);
            i = events_copy.erase(i);
            continue;
          }
          else
            i++;
        }
        else
          assert(false);
      }
    }
  }

  // compute minimal event sets
  for (list<list<Event*> >::iterator i = groups.begin(); i != groups.end(); i++)
    determine_minimal_set(*i);

  FILE_LOG(LOG_CONTACT) << "Event::determine_connected_events() exited" << std::endl;
}

/// Modified Gaussian elimination with partial pivoting -- computes half-rank at the same time
unsigned Event::gauss_elim(MatrixN& A)
{
  const unsigned NROWS = A.rows(), NCOLS = A.columns();
  SAFESTATIC vector<unsigned> piv;

  // resize the pivots array
  piv.resize(NCOLS);
  for (unsigned i=0; i< NCOLS; i++)
    piv[i] = i;

  // setup swpi
  unsigned swpi = 0;

  for (unsigned i=0; i< NROWS; i++)
  {
    // get the pointer to the row
    BlockIterator row = A.block_start(i,i+1, 0, NCOLS);

    // find the largest positive element in the row 
    Real largest = row[piv[swpi]];
    unsigned largest_index = swpi;
    for (unsigned k=swpi; k< NCOLS; k++)
      if (row[piv[k]] > largest)
      {
        largest_index = k;
        largest = row[piv[k]];
      }

    // continue processing only if largest positive element is greater than 0 
    if (largest > std::numeric_limits<Real>::epsilon())
    {
      // swap the pivots
      std::swap(piv[swpi], piv[largest_index]);

      // reduce the columns
      const Real* elm = &A(i, piv[swpi]);
      for (unsigned k=swpi+1; k< NCOLS; k++)
      {
        Real* elm_k = &A(i, piv[k]);
        if (*elm_k > std::numeric_limits<Real>::epsilon())
        {
          Real scal = -*elm_k / *elm;
          CBLAS::axpy(NROWS-i, scal, elm, 1, &A(i,piv[k]), 1);
        }
      }

      // update the swap pivot
      swpi++;

      // quit if swpi too large
      if (swpi == NCOLS)
        break;
    }

    // find the largest negative element in the row
    largest = row[piv[swpi]];
    largest_index = swpi;
    for (unsigned k=swpi+1; k< NCOLS; k++)
      if (row[piv[k]] < largest)
      {
        largest_index = k;
        largest = row[piv[k]];
      }

    // only continue processing if largest negative element is less than 0
    if (largest < -std::numeric_limits<Real>::epsilon())
    {
      // swap the pivots
      std::swap(piv[swpi], piv[largest_index]);
      
      // reduce the columns
      const Real* elm = &A(i, piv[swpi]);
      for (unsigned k=swpi+1; k< NCOLS; k++)
      {
        Real* elm_k = &A(i, piv[k]);
        if (*elm_k < -std::numeric_limits<Real>::epsilon())
        {
          Real scal = -*elm_k / *elm;
          CBLAS::axpy(NROWS-i, scal, elm, 1, &A(i,piv[k]), 1);
        }
      }

      // update the swap pivot
      swpi++;

      // quit if swpi too large
      if (swpi == NCOLS)
        break;
    }
  }

  return swpi;
}

/// Computes normal and contact Jacobians for a body
void Event::compute_contact_jacobians(const Event& e, MatrixN& Jc, MatrixN& Dc, const map<DynamicBodyPtr, unsigned>& gc_indices)
{
  map<DynamicBodyPtr, unsigned>::const_iterator miter;
  SAFESTATIC FastThreadable<VectorN> tmpv;

  // zero the matrices
  Jc.set_zero();
  Dc.set_zero();

  // get the two bodies
  SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
  SingleBodyPtr sb2 = e.contact_geom2->get_single_body();

  // get the super bodies
  DynamicBodyPtr ab1 = sb1->get_articulated_body();
  DynamicBodyPtr ab2 = sb2->get_articulated_body();
  DynamicBodyPtr super1 = (ab1) ? ab1 : sb1;
  DynamicBodyPtr super2 = (ab2) ? ab2 : sb2;

  // process the first body
  miter = gc_indices.find(super1);
  if (miter != gc_indices.end())
  {
    const unsigned index = miter->second;

    // compute the 'r' vector
    Vector3 r = e.contact_point - sb1->get_position();

    // convert the normal force to generalized forces
    super1->convert_to_generalized_force(DynamicBody::eAxisAngle, sb1, e.contact_normal, Vector3::cross(r, e.contact_normal), tmpv());
    Jc.set_sub_mat(index, 0, tmpv());

    // convert first tangent direction to generalized forces
    super1->convert_to_generalized_force(DynamicBody::eAxisAngle, sb1, e.contact_tan1, Vector3::cross(r, e.contact_tan1), tmpv());
    Dc.set_sub_mat(index, 0, tmpv());

    // convert second tangent direction to generalized forces
    super1->convert_to_generalized_force(DynamicBody::eAxisAngle, sb1, e.contact_tan2, Vector3::cross(r, e.contact_tan2), tmpv());
    Dc.set_sub_mat(index, 1, tmpv());
  }

  // process the second body
  miter = gc_indices.find(super2);
  if (miter != gc_indices.end())
  {
    const unsigned index = miter->second;

    // compute the 'r' vector
    Vector3 r = e.contact_point - sb2->get_position();

    // convert the normal force to generalized forces
    super2->convert_to_generalized_force(DynamicBody::eAxisAngle, sb2, -e.contact_normal, Vector3::cross(r, -e.contact_normal), tmpv());
    Jc.set_sub_mat(index, 0, tmpv());

    // convert first tangent direction to generalized forces
    super2->convert_to_generalized_force(DynamicBody::eAxisAngle, sb2, -e.contact_tan1, Vector3::cross(r, -e.contact_tan1), tmpv());
    Dc.set_sub_mat(index, 0, tmpv());

    // convert second tangent direction to generalized forces
    super2->convert_to_generalized_force(DynamicBody::eAxisAngle, sb2, -e.contact_tan2, Vector3::cross(r, -e.contact_tan2), tmpv());
    Dc.set_sub_mat(index, 1, tmpv());
  }

  // equilibrate the columns of Jc and Dc 
  const unsigned NGC = Jc.rows();
  Real* cols[3] = { &Jc(0,0), &Dc(0,0), &Dc(0,1) };
  for (unsigned j=0; j< 3; j++)
  {
    Real* col = cols[j];
    Real max_val = (Real) 0.0;
    for (unsigned i=0; i< NGC; i++)
      max_val = std::max(max_val, std::fabs(col[i]));
    assert(max_val > (Real) 0.0);
    CBLAS::scal(NGC, (Real) 1.0/max_val, col, 1);
  }
}

/// Computes a minimal set of contact events
void Event::determine_minimal_set(list<Event*>& group)
{
  FILE_LOG(LOG_CONTACT) << "Event::determine_minimal_set() entered" << std::endl;
  FILE_LOG(LOG_CONTACT) << " -- initial number of events: " << group.size() << std::endl;

  // get the number of contact events and total number of events
  list<Event*>::iterator start = group.begin();
  unsigned NC = 0, NE = 0;
  while (start != group.end() && ++NE)
  {
    if ((*start)->event_type == Event::eContact)
    {
      (*start)->determine_contact_tangents();
      NC++;
    }
    start++;
  }

  // if there is one or no contacts, quit now
  if (NC <= 1 || NE < 4)
  {
    FILE_LOG(LOG_CONTACT) << " -- initial/final number of contacts: " << NC << std::endl;
    FILE_LOG(LOG_CONTACT) << " -- initial/final number of events: " << NE << std::endl;
    return;
  }

  // determine the number of gc's in the group
  unsigned NGC = 0;
  map<DynamicBodyPtr, unsigned> gc_index;
  vector<DynamicBodyPtr> supers;
  for (list<Event*>::const_iterator i = group.begin(); i != group.end(); i++)
  {
    supers.clear();
    (*i)->get_super_bodies(std::back_inserter(supers));
    for (unsigned j=0; j< supers.size(); j++)
      if (gc_index.find(supers[j]) == gc_index.end())
      {
        gc_index[supers[j]] = NGC;
        NGC += supers[j]->num_generalized_coordinates(DynamicBody::eAxisAngle);
      }
  }

  // initialize the Jacobian matrices (as big as may be necessary)
  MatrixN Jc(NGC, NC), Dc(NGC, NC*2);
  MatrixN sub_Jc(NGC, 1), sub_Dc(NGC, 2);
  MatrixN workM(NGC, NC*2);

  // get the pointer to the first contact
  start = group.begin();
  while ((*start)->event_type != Event::eContact)
    start++;

  // get the equilibrated Jacobians (normal & tangential) for the first contact
  Jc.resize(NGC, 1);
  Dc.resize(NGC, 2);
  compute_contact_jacobians(**start, Jc, Dc, gc_index);

  // setup current rank of the tangent Jacobians
  unsigned Jc_rank = 1;
  unsigned Dc_rank = 2;

  // verify that the column rank is two
  #ifndef NDEBUG
  assert(LinAlg::calc_rank(Dc) == 2);
  compute_contact_jacobians(**start, Jc, Dc, gc_index);
  #endif

  // setup the contact index
  unsigned ci = 1;

  // loop through the remainder of contacts
  for (list<Event*>::iterator i = ++start; i != group.end(); )
  {
    // if this isn't a contact event, skip it
    if ((*i)->event_type != Event::eContact)
    {
      i++;
      continue;
    }

    // get the equilibrated Jacobians (normal & tangential) for this contact
    compute_contact_jacobians(**i, sub_Jc, sub_Dc, gc_index);

    // setup candidate Jacobian matrices
    Jc.resize(NGC, ci+1);
    Dc.resize(NGC, (ci+1)*2);
    Jc.set_sub_mat(0,ci,sub_Jc);
    Dc.set_sub_mat(0,ci*2,sub_Dc);

    // try putting normal Jacobian into rref
    unsigned new_Jc_rank = (Jc_rank == NGC*2) ? NGC*2 : gauss_elim(Jc);

    // see whether adding tangential Jacobians increases rank
    unsigned new_Dc_rank = (Dc_rank == NGC) ? NGC : LinAlg::calc_rank(workM.copy_from(Dc));

    // if we don't have a rank increase from add either normal or tangential
    // Jacobians, remove the contact event
    if (!(new_Dc_rank > Dc_rank || new_Jc_rank > Jc_rank))
    {
      // erase the contact event
      i = group.erase(i);
    }
    else
    {
      // update Jc rank and Dc rank
      Jc_rank = new_Jc_rank;
      Dc_rank = new_Dc_rank;

      // advance to the next contact
      ci++;

      // advance to the next event
      i++;
    }
  }

  FILE_LOG(LOG_CONTACT) << " -- final number of events: " << group.size() << std::endl;
}

/// Removes groups of contacts that contain no impacts
void Event::remove_nonimpacting_groups(list<list<Event*> >& groups, Real tol)
{
  typedef list<list<Event*> >::iterator ListIter;

  for (ListIter i = groups.begin(); i != groups.end(); )
  {
    // look for impact in list i
    bool impact_detected = false;
    BOOST_FOREACH(Event* e, *i)
    {
      if (e->is_impacting(tol))
      {
        impact_detected = true;
        break;
      }
    }

    // if no impact in the list, remove the list
    if (!impact_detected)
    {
      ListIter j = i;
      j++;
      groups.erase(i);
      i = j;
    }
    else
      i++;
  }
}

/// Writes an event to the specified filename in VRML format for visualization
/**
 * \todo add a cone onto the arrows
 */
void Event::write_vrml(const std::string& fname, Real sphere_radius, Real normal_length) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  std::ofstream out;
  
  // open the file for writing
  out.open(fname.c_str());
  if (out.fail())
    throw std::runtime_error("Unable to open file for writing in Event::write_vrml()");

  // write the VRML header
  out << "#VRML V2.0 utf8" << std::endl << std::endl;

  // *************************************************
  // first, write the contact point 
  // *************************************************

  // determine a random color that will be used for contact and normal
  Real c_x = (Real) rand() / RAND_MAX;
  Real c_y = (Real) rand() / RAND_MAX;
  Real c_z = (Real) rand() / RAND_MAX;

  // write the transform for the contact point
  out << "Transform {" << std::endl;
  out << "  translation "; 
  out << contact_point[X] << " " << contact_point[Y] << " " << contact_point[Z] << std::endl;
  out << "  children " << endl;

  // write the shape node, using default appearance
  out << "  Shape {" << std::endl;
  out << "    appearance Appearance { material Material {" << std::endl;
  out << "      transparency 0" << std::endl;
  out << "      shininess 0.2" << std::endl;
  out << "      ambientIntensity 0.2" << std::endl;
  out << "      emissiveColor 0 0 0" << std::endl;
  out << "      specularColor 0 0 0" << std::endl;
  out << "      diffuseColor " << c_x << " " << c_y << " " << c_z << std::endl;
  out << "      }}" << std::endl;

  // write the geometry (a sphere)
  out << "  geometry Sphere {" << std::endl; 
  out << "    radius " << sphere_radius << " }}} # end sphere, shape, transform " << std::endl;

  // *************************************************
  // now, write the normal
  // *************************************************

  // determine the normal edge
  Vector3 normal_start = contact_point;
  Vector3 normal_stop = normal_start + contact_normal*normal_length;

  // write the shape node, using default appearance
  out << "Shape {" << std::endl;
  out << "  appearance Appearance { material Material {" << std::endl;
  out << "    transparency 0" << std::endl;
  out << "    shininess 0.2" << std::endl;
  out << "    ambientIntensity 0.2" << std::endl;
  out << "    emissiveColor 0 0 0" << std::endl;
  out << "    specularColor 0 0 0" << std::endl;
  out << "    diffuseColor " << c_x << " " << c_y << " " << c_z << std::endl;
  out << "    }}" << std::endl;

  // write the geometry
  out << "  geometry IndexedLineSet {" << std::endl; 
  out << "    coord Coordinate { point [ ";
  out << normal_start[X] << " " << normal_start[Y] << " " << normal_start[Z] << ", ";
  out << normal_stop[X] << " " << normal_stop[Y] << " " << normal_stop[Z] << " ] } " << std::endl;
  out << "    coordIndex [ 0, 1, -1 ] }}" << std::endl;

  // **********************************************
  // determine the axis-angle rotation for the cone
  // **********************************************

  // first compose an arbitrary vector d
  Vector3 d(1,1,1);
  if (std::fabs(contact_normal[X]) > std::fabs(contact_normal[Y]))
  {
    if (std::fabs(contact_normal[X]) > std::fabs(contact_normal[Z]))
      d[X] = 0;
    else
      d[Z] = 0;
  }
  else
  {
    if (std::fabs(contact_normal[Y]) > std::fabs(contact_normal[Z]))
      d[Y] = 0;
    else
      d[Z] = 0;
  }
    
  // compute the cross product of the normal and the vector
  Vector3 x = Vector3::normalize(Vector3::cross(contact_normal, d));
  Vector3 y;
  y.copy_from(contact_normal);
  Vector3 z = Vector3::normalize(Vector3::cross(x, contact_normal));

  // compute theta and the axis of rotation
  Real theta = std::acos((x[X] + y[Y] + z[Z] - 1)/2);
  Vector3 axis(z[Y] - y[Z], x[Z] - z[X], y[X] - x[Y]);
  axis *= -(1.0/(2 * std::sin(theta)));
    
  // finally, write the cone to show the normal's direction
  out << "Transform {" << std::endl;
  out << "  rotation ";
   out  << axis[X] <<" "<< axis[1] <<" "<< axis[Z] <<" "<< theta << std::endl;
  out << "  translation ";
   out << normal_stop[X] <<" "<< normal_stop[Y] <<" "<< normal_stop[Z];
  out << std::endl;
  out << "  children [" << std::endl;
  out << "    Shape {" << std::endl;
  out << "      appearance Appearance { material Material {" << std::endl;
  out << "        transparency 0" << std::endl;
  out << "        shininess 0.2" << std::endl;
  out << "        ambientIntensity 0.2" << std::endl;
  out << "        emissiveColor 0 0 0" << std::endl;
  out << "        specularColor 0 0 0" << std::endl;
  out << "        diffuseColor " << c_x << " " << c_y << " " << c_z << std::endl;
  out << "        }}" << std::endl;
  out << "      geometry Cone {" << std::endl;
  out << "        bottomRadius " << sphere_radius << std::endl;
  out << "        height " << (normal_length * .1) << std::endl;
  out << "      } } ] }" << std::endl;
  out.close();
}

/// Determines the set of contact tangents
void Event::determine_contact_tangents()
{
  assert(event_type == Event::eContact);

  // determine an orthonormal basis using the two contact tangents
  Vector3::determine_orthonormal_basis(contact_normal, contact_tan1, contact_tan2);

/*
  // determine the projection matrix
  Matrix3 RT = Matrix3::transpose(CompGeom::calc_3D_to_2D_matrix(contact_normal));

  // setup some necessary constants
  const Real TWO_INV_NK = (Real) 1.0/contact_NK;
  const Real TWO_PI_INV_NK = M_PI * TWO_INV_NK;

  // setup the tangents
  for (unsigned i=0; i< contact_NK; i++)
  {
    Vector2 dx(std::cos(i*TWO_PI_INV_NK), std::sin(i*TWO_PI_INV_NK));
    contact_tangents[i] = CompGeom::to_3D(dx, RT);
  }
*/
}

/// Determines the type of event (impacting, resting, or separating)
Event::EventClass Event::determine_event_class(Real tol) const
{
  // get the event velocity
  Real vel = calc_event_vel();

  FILE_LOG(LOG_SIMULATOR) << "-- event type: " << event_type << " velocity: " << vel << std::endl;

  // if the velocity is less than zero, we have
  // an impacting contact
  if (vel > tol)
    return eSeparating;
  else if (vel < -tol)
    return eImpacting;
  else
    return eResting;
}

