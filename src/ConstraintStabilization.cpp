/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/**
 * Thoughts about separating bodies:
 *
 * If we really want to separate bodies with minimum work, then we'll have to
 * solve an optimal motion planning problem (fun!)- not going to be fast enough 
 * for most simulations. 
 *
 * If we apply a force to the bodies at the "wrong" contact point (like pushing
 * the box out of the plane), this will not lead to an increase in signed
 * distance.
 *
 * If we separate bodies using the vector of minimum translational distance as
 * the normal and we put a contact point halfway through the two bodies
 * (this actually will have no effect on torque), then that will locally
 * minimize the work necessary to separate the bodies. A root finding process
 * can be used to determine when the bodies will be separated.
 *
 * The other concern is what to do when the bodies are separated. If two
 * bodies are separated, then then the root finding process can determine
 * when (if ever) they come into contact.
 */

#include <map>
#include <Moby/Types.h>
#include <Moby/ConstraintSimulator.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/ConstraintStabilization.h>

using namespace Ravelin;
using namespace Moby;
using boost::weak_ptr;
using boost::shared_ptr;
using boost::shared_array;
using boost::const_pointer_cast; 
using boost::dynamic_pointer_cast; 
using std::map;
using std::vector;
using std::make_pair;
using std::endl;
using std::list;

ConstraintStabilization::ConstraintStabilization()
{
  // set tolerance to NEAR_ZERO by default
  eps = NEAR_ZERO;
}

/// Saves the velocities before constraint stabilization
void ConstraintStabilization::save_velocities(shared_ptr<ConstraintSimulator> sim, vector<VectorNd>& qd)
{
  qd.clear();
  for (unsigned i = 0; i< sim->_bodies.size(); i++)
  {
    qd.push_back(VectorNd());
    sim->_bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, qd.back());
  }
}

/// Restores the velocities after constraint stabilization
void ConstraintStabilization::restore_velocities(shared_ptr<ConstraintSimulator> sim, const vector<VectorNd>& qd)
{
  for (unsigned i = 0; i< sim->_bodies.size(); i++)
    sim->_bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, qd[i]);
}

/// Gets the minimum pairwise distance
double ConstraintStabilization::get_min_pairwise_dist(const vector<PairwiseDistInfo>& pdi)
{
  // set distance to infinite initially
  double dist = std::numeric_limits<double>::max();

  for (unsigned i=0; i< pdi.size(); i++)
    dist = std::min(dist, pdi[i].dist);

  return dist;
}

/// Stabilizes the constraints in the simulator
void ConstraintStabilization::stabilize(shared_ptr<ConstraintSimulator> sim)
{
  FILE_LOG(LOG_SIMULATOR)<< "======constraint stabilization start======"<<std::endl;
  VectorNd dq, q, v;
  std::vector<UnilateralConstraintProblemData> pd;
  std::vector<VectorNd> qd_save;

  std::map<DynamicBodyPtr, unsigned> body_index_map;

  // save the generalized velocities
  save_velocities(sim, qd_save);

  get_body_configurations(q, sim);
  generate_body_index_map(body_index_map, sim);

  // get the pairwise distances
  vector<PairwiseDistInfo>& pdi = sim->_pairwise_distances;

  // see whether any pairwise distances are below epsilon
  double min_dist = get_min_pairwise_dist(pdi);

  FILE_LOG(LOG_SIMULATOR) <<"minimum pairwise distance (before stabilization loop): "<<min_dist<<std::endl;
  unsigned iterations = 0;
  while (min_dist < eps)
  {
    // zero body velocities first (we only want to change positions based on
    // our updates)
    for (unsigned i=0; i< sim->_bodies.size(); i++)
    {
      sim->_bodies[i]->get_generalized_velocity(DynamicBody::eSpatial, v);
      v.set_zero();
      sim->_bodies[i]->set_generalized_velocity(DynamicBody::eSpatial, v);
    }

    // compute problem data (get M, N, alpha, etc.) 
    compute_problem_data(pd, sim);

    // determine dq's
    dq.set_zero(q.size());
    for (unsigned i=0; i< pd.size(); i++)
      determine_dq(pd[i], dq, body_index_map);

    // determine s and update q 
    update_q(dq, q, sim);

    // update minimum distance
    min_dist = get_min_pairwise_dist(pdi);  
    iterations++;
  }

  // restore the generalized velocities
  restore_velocities(sim, qd_save);

  FILE_LOG(LOG_SIMULATOR) << iterations << " iterations required" << std::endl;
  FILE_LOG(LOG_SIMULATOR) <<"=====constraint stabilization end ======" << std::endl;
}

/// Adds unilateral constraints for joint limits in an articulated body
void ConstraintStabilization::add_articulate_limit_constraint(std::vector<UnilateralConstraint>& constraints, ArticulatedBodyPtr ab)
{
  std::vector<UnilateralConstraint> limits;
  ab->find_limit_constraints(std::back_inserter(limits));
  constraints.insert(constraints.end(), limits.begin(), limits.end());
}

/// Adds contact constraints from a pair of bodies
void ConstraintStabilization::add_contact_constraints(std::vector<UnilateralConstraint>& constraints, RigidBodyPtr rb1, RigidBodyPtr rb2, shared_ptr<ConstraintSimulator> sim)
{
  boost::shared_ptr<const Pose3d> GLOBAL;
  Point3d p1, p2;
  std::list<CollisionGeometryPtr>& cgs1 = rb1->geometries;
  std::list<CollisionGeometryPtr>& cgs2 = rb2->geometries;
 
  BOOST_FOREACH(CollisionGeometryPtr cg1 , cgs1)
  {
    BOOST_FOREACH(CollisionGeometryPtr cg2, cgs2)
    {
      // only process each pair once
      if (cg1.get() >= cg2.get())
        continue;

      double dist = CollisionGeometry::calc_signed_dist(cg1, cg2, p1, p2);

      // case 1: bodies are "kissing" 
      if (std::fabs(dist) < NEAR_ZERO)
        sim->_coldet->find_contacts(cg1,cg2, constraints);
      // case 2: bodies are separated      
      else if (dist >= NEAR_ZERO)
      {
        boost::shared_ptr<const Pose3d> GLOBAL;
        boost::shared_ptr<const Pose3d> pose1(p1.pose);
        boost::shared_ptr<const Pose3d> pose2(p2.pose);
        Transform3d _1T2 = Ravelin::Pose3d::calc_relative_pose(pose1, pose2);
        Transform3d _1TG = Ravelin::Pose3d::calc_relative_pose(pose1, GLOBAL);
        Point3d p1_2 = _1T2.transform_point(p1);
        Point3d p1_g = _1TG.transform_point(p1);
        Ravelin::Vector3d normal = p2-p1_2;
        normal.normalize();
        UnilateralConstraint uc = CollisionDetection::create_contact(cg1, cg2, p1_g, normal, dist);
        //FILE_LOG(LOG_SIMULATOR) << "p1: " << p1_g << std::endl << "normal" << normal << std::endl << "dist" << dist<<std::endl;
        constraints.insert(constraints.end(), uc);
      }
      // case 3: bodies are properly interpenetrating
      else
      {
        // get the two rigid bodies
        RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(cg1->get_single_body());
        RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(cg2->get_single_body());

        // get the two poses of the rigid bodies
        boost::shared_ptr<const Pose3d> pose1 = rb1->get_inertial_pose();
        boost::shared_ptr<const Pose3d> pose2 = rb2->get_inertial_pose();

        // get the centers of mass in the global frame
        Vector3d rb1_c(0.0, 0.0, 0.0, pose1);
        Vector3d rb2_c(0.0, 0.0, 0.0, pose2);
        Vector3d rb1_c0 = Pose3d::transform_point(GLOBAL, rb1_c);
        Vector3d rb2_c0 = Pose3d::transform_point(GLOBAL, rb2_c);       

        // setup the contact point directly between them
        Point3d cp = rb1_c0*0.5 + rb2_c0*0.5;

        // setup the normal to point toward body 1
        Vector3d normal = Vector3d::normalize(rb1_c0 - rb2_c0);

        UnilateralConstraint uc = CollisionDetection::create_contact(cg1, cg2, cp, normal, dist);
        //FILE_LOG(LOG_SIMULATOR) << "p1: " << p1_g << std::endl << "normal" << normal << std::endl << "dist" << dist<<std::endl;
        constraints.insert(constraints.end(), uc);
      }
    }
  }
}

/// Computes the constraint data
void ConstraintStabilization::compute_problem_data(std::vector<UnilateralConstraintProblemData>& pd_vector, shared_ptr<ConstraintSimulator> sim)
{
  std::vector<UnilateralConstraint> constraints;

  // clear the problem data vector 
  pd_vector.clear();

  // get all bodies
  const std::vector<DynamicBodyPtr>& bodies = sim->_bodies;
  // TODO: Evan add sweep and prune

  // 1) for each pair of bodies in kissing contact, add as many 
  //    UnilateralConstraint objects to constraints as there are 
  //    points of contact between the bodies 
  //    (call _sim->_coldet->find_contacts(.))
  //FILE_LOG(LOG_SIMULATOR) << "*******start adding constraints*******" << std::endl;
  BOOST_FOREACH(DynamicBodyPtr D_body1, bodies)
  {
    RigidBodyPtr rb1, rb2;
    
    if((rb1 = boost::dynamic_pointer_cast<RigidBody>(D_body1)))
    {
      BOOST_FOREACH(DynamicBodyPtr D_body2, bodies)
      {
        // if the two pointer points to the same body, then no need to add contact
        if(D_body1 == D_body2)
        {
          continue;
        }

        //RigidBody
        if((rb2 = boost::dynamic_pointer_cast<RigidBody>(D_body2)))
        {
          add_contact_constraints(constraints, rb1,rb2, sim);
        }
        else
        {
          ArticulatedBodyPtr ab2 = dynamic_pointer_cast<ArticulatedBody>(D_body2);
          add_articulate_limit_constraint(constraints, ab2);
          const std::vector<RigidBodyPtr>& ls2 = ab2->get_links();
          BOOST_FOREACH(RigidBodyPtr l2, ls2)
          {
            rb2 = dynamic_pointer_cast<RigidBody> (l2);
            add_contact_constraints(constraints, rb1, rb2, sim);
          }
        }
      }
    }
    // body 1 is a articulated body
    else
    {
      ArticulatedBodyPtr ab1 = dynamic_pointer_cast<ArticulatedBody>(D_body1);
      add_articulate_limit_constraint(constraints, ab1);
      std::vector<RigidBodyPtr> ls1 = ab1->get_links();

      BOOST_FOREACH(DynamicBodyPtr D_body2, bodies)
      {
        // if the two pointer points to the same body, then no need to add contact
        if(D_body1 == D_body2)
        {
          continue;
        }

        // since the two pointer are not pointing to the same body, it is ok to start iterating through the first  
        BOOST_FOREACH(RigidBodyPtr l1, ls1)
        {
          rb1 = l1;
          //RigidBody
          if((rb2 = boost::dynamic_pointer_cast<RigidBody>(D_body2)))
          {
            add_contact_constraints(constraints, rb1,rb2, sim);
          }
          else
          {
            ArticulatedBodyPtr ab2 = dynamic_pointer_cast<ArticulatedBody>(D_body2);
            add_articulate_limit_constraint(constraints, ab2);
            const std::vector<RigidBodyPtr>& ls2 = ab2->get_links();
            BOOST_FOREACH(RigidBodyPtr l2, ls2)
            {
              rb2 = l2;
              add_contact_constraints(constraints, rb1, rb2, sim);
            }
          }
        }
      }
    }
  }
  //FILE_LOG(LOG_SIMULATOR) << "constraints added" << std::endl;
  // 2) for each articulated body, add as many UnilateralConstraint objects as
  //    there are joints at their limits


  // 3) for each pair of bodies in interpenetrating contact, add a single
  //    point of contact at the deepest interpenetrating point with normal
  //    in the direction of the signed distance function. 

  // find islands
  list<list<UnilateralConstraint*> > islands;
  UnilateralConstraint::determine_connected_constraints(constraints, islands);

  // process islands
  BOOST_FOREACH(list<UnilateralConstraint*>& island, islands)
  {
    // setup a UnilateralConstraintProblemData object
    pd_vector.push_back(UnilateralConstraintProblemData());
    UnilateralConstraintProblemData& pd = pd_vector.back();

    // put the constraint into the appropriate place
    BOOST_FOREACH(UnilateralConstraint* c, island)
    { 
      if (c->constraint_type == UnilateralConstraint::eContact)
        pd.contact_constraints.push_back(c);
      else
        pd.limit_constraints.push_back(c);
    }

    // set number of contact and limit constraints
    pd.N_CONTACTS = pd.contact_constraints.size();
    pd.N_LIMITS = pd.limit_constraints.size(); 
    pd.N_K_TOTAL = 0;
    pd.N_TRUE_CONE = 0;

    // now set the unilateral constraint data
    set_unilateral_constraint_data(pd);
    
    // set the elements of Cn_v and L_v
    // L_v is always set to zero
    // Cn_v is set to the signed distance between the two bodies
    Point3d pa,pb;
    for (unsigned i = 0; i < pd.contact_constraints.size(); ++i)
    {
      pd.Cn_v[i] = CollisionGeometry::calc_signed_dist(pd.contact_constraints[i]->contact_geom1, pd.contact_constraints[i]->contact_geom2, pa, pb) - eps;
    }
  }
}

/// Gets the super body (articulated if any)
DynamicBodyPtr ConstraintStabilization::get_super_body(SingleBodyPtr sb)
{
  ArticulatedBodyPtr ab = sb->get_articulated_body();
  if (ab)
    return ab;
  else
    return sb;
}

/// Computes the data to the LCP / QP problems
void ConstraintStabilization::set_unilateral_constraint_data(UnilateralConstraintProblemData& pd)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  MatrixNd MM;
  VectorNd v;

  // determine set of "super" bodies from contact constraints
  pd.super_bodies.clear();
  for (unsigned i=0; i< pd.contact_constraints.size(); i++)
  {
    pd.super_bodies.push_back(get_super_body(pd.contact_constraints[i]->contact_geom1->get_single_body()));
    pd.super_bodies.push_back(get_super_body(pd.contact_constraints[i]->contact_geom2->get_single_body()));
  }

  // determine set of "super" bodies from limit constraints
  for (unsigned i=0; i< pd.limit_constraints.size(); i++)
  {
    RigidBodyPtr outboard = pd.limit_constraints[i]->limit_joint->get_outboard_link();
    pd.super_bodies.push_back(get_super_body(outboard));
  }

  // make super bodies vector unique
  std::sort(pd.super_bodies.begin(), pd.super_bodies.end());
  pd.super_bodies.erase(std::unique(pd.super_bodies.begin(), pd.super_bodies.end()), pd.super_bodies.end());

  // set total number of generalized coordinates
  pd.N_GC = 0;
  for (unsigned i=0; i< pd.super_bodies.size(); i++)
    pd.N_GC += pd.super_bodies[i]->num_generalized_coordinates(DynamicBody::eSpatial);

  // initialize constants and set easy to set constants
  pd.N_CONTACTS = pd.contact_constraints.size();
  pd.N_LIMITS = pd.limit_constraints.size();

  // setup constants related to articulated bodies
  for (unsigned i=0; i< pd.super_bodies.size(); i++)
  {
    ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(pd.super_bodies[i]);
    if (abody) {
      pd.N_CONSTRAINT_EQNS_IMP += abody->num_constraint_eqns_implicit();
    }
  }

  // setup number of friction polygon edges / true cones
  pd.N_K_TOTAL = 0;
  pd.N_LIN_CONE = 0;
  pd.N_TRUE_CONE = 0; 

  // initialize the problem matrices / vectors
  pd.Cn_iM_CnT.set_zero(pd.N_CONTACTS, pd.N_CONTACTS);
  pd.Cn_iM_CsT.set_zero(pd.N_CONTACTS, pd.N_CONTACTS);
  pd.Cn_iM_CtT.set_zero(pd.N_CONTACTS, pd.N_CONTACTS);
  pd.Cn_iM_LT.set_zero(pd.N_CONTACTS, pd.N_LIMITS);
  pd.Cn_iM_JxT.set_zero(pd.N_CONTACTS, pd.N_CONSTRAINT_EQNS_IMP);
  pd.Cs_iM_CsT.set_zero(pd.N_CONTACTS, pd.N_CONTACTS);
  pd.Cs_iM_CtT.set_zero(pd.N_CONTACTS, pd.N_CONTACTS);
  pd.Cs_iM_LT.set_zero(pd.N_CONTACTS, pd.N_LIMITS);
  pd.Cs_iM_JxT.set_zero(pd.N_CONTACTS, pd.N_CONSTRAINT_EQNS_IMP);
  pd.Ct_iM_CtT.set_zero(pd.N_CONTACTS, pd.N_CONTACTS);
  pd.Ct_iM_LT.set_zero(pd.N_CONTACTS, pd.N_LIMITS);
  pd.Ct_iM_JxT.set_zero(pd.N_CONTACTS, pd.N_CONSTRAINT_EQNS_IMP);
  pd.L_iM_LT.set_zero(pd.N_LIMITS, pd.N_LIMITS);
  pd.L_iM_JxT.set_zero(pd.N_LIMITS, pd.N_CONSTRAINT_EQNS_IMP);
  pd.Jx_iM_JxT.set_zero(pd.N_CONSTRAINT_EQNS_IMP, pd.N_CONSTRAINT_EQNS_IMP);
  pd.Cn_v.set_zero(pd.N_CONTACTS);
  pd.Cs_v.set_zero(pd.N_CONTACTS);
  pd.Ct_v.set_zero(pd.N_CONTACTS);
  pd.L_v.set_zero(pd.N_LIMITS);
  pd.Jx_v.set_zero(pd.N_CONSTRAINT_EQNS_IMP);
  pd.cn.set_zero(pd.N_CONTACTS);
  pd.cs.set_zero(pd.N_CONTACTS);
  pd.ct.set_zero(pd.N_CONTACTS);
  pd.l.set_zero(pd.N_LIMITS);
  pd.alpha_x.set_zero(pd.N_CONSTRAINT_EQNS_IMP);

  // setup indices
  pd.CN_IDX = 0;
  pd.CS_IDX = pd.CN_IDX + pd.N_CONTACTS;
  pd.CT_IDX = pd.CS_IDX;
  pd.NCS_IDX = pd.CS_IDX;
  pd.NCT_IDX = pd.CS_IDX;
  pd.L_IDX = pd.CS_IDX;
  pd.ALPHA_X_IDX = pd.L_IDX + pd.N_LIMITS;
  pd.N_VARS = pd.ALPHA_X_IDX + pd.N_CONSTRAINT_EQNS_IMP;

  // get iterators to the proper matrices
  RowIteratord CnCn = pd.Cn_iM_CnT.row_iterator_begin();
  RowIteratord CnCs = pd.Cn_iM_CsT.row_iterator_begin();
  RowIteratord CnCt = pd.Cn_iM_CtT.row_iterator_begin();
  RowIteratord CsCs = pd.Cs_iM_CsT.row_iterator_begin();
  RowIteratord CsCt = pd.Cs_iM_CtT.row_iterator_begin();
  RowIteratord CtCt = pd.Ct_iM_CtT.row_iterator_begin();

  // process contact constraints, setting up matrices
  for (unsigned i=0; i< pd.contact_constraints.size(); i++)
  {
    // compute cross constraint data for contact constraints
    for (unsigned j=0; j< pd.contact_constraints.size(); j++)
    {
      // reset MM
      MM.set_zero(3, 3);

      // check whether i==j (single contact constraint)
      if (i == j)
      {
        // compute matrix / vector for contact constraint i
        v.set_zero(3);
        pd.contact_constraints[i]->compute_constraint_data(MM, v);

        // setup appropriate part of contact inertia matrices
        RowIteratord_const data = MM.row_iterator_begin();
        *CnCn = *data++;
      }

      // advance the iterators
      CnCn++;
    }

    // compute cross constraint data for contact/limit constraints
    for (unsigned j=0; j< pd.limit_constraints.size(); j++)
    {
      // reset MM
      MM.set_zero(3, 1);

      // compute matrix for cross constraint
      pd.contact_constraints[i]->compute_cross_constraint_data(*pd.limit_constraints[j], MM);

      // setup appropriate parts of contact / limit inertia matrices
      ColumnIteratord_const data = MM.column_iterator_begin();
      pd.Cn_iM_LT(i,j) = *data++;
    }
  }

  // process limit constraints, setting up matrices
  for (unsigned i=0; i< pd.limit_constraints.size(); i++)
  {
    // compute matrix / vector for contact constraint i
    pd.limit_constraints[i]->compute_constraint_data(MM, v);

    // setup appropriate entry of limit inertia matrix and limit velocity
    pd.L_iM_LT(i,i) = MM.data()[0];

    // compute cross/cross limit constraint data
    for (unsigned j=i+1; j< pd.limit_constraints.size(); j++)
    {
      // reset MM
      MM.resize(1,1);

      // compute matrix for cross constraint
      pd.limit_constraints[i]->compute_cross_constraint_data(*pd.limit_constraints[j], MM);

      // setup appropriate part of limit / limit inertia matrix
      pd.L_iM_LT(i,j) = pd.L_iM_LT(j,i) = MM.data()[0];
    }

    // NOTE: cross data has already been computed for contact/limit constraints
  }
}

/// Computes deltaq by solving a linear complementarity problem
void ConstraintStabilization::determine_dq(UnilateralConstraintProblemData& pd, VectorNd& dqm, const std::map<DynamicBodyPtr, unsigned>& body_index_map)
{
  VectorNd z, dq_sub;

  // initialize the LCP matrix and LCP vector
  MatrixNd MM(pd.N_CONTACTS + pd.N_LIMITS, pd.N_CONTACTS + pd.N_LIMITS);
  VectorNd qq(pd.N_CONTACTS + pd.N_LIMITS);

  // setup the LCP matrix and LCP vector
  MM.block(0, pd.N_CONTACTS, 0, pd.N_CONTACTS) = pd.Cn_iM_CnT;
  MM.block(0, pd.N_CONTACTS, pd.N_CONTACTS, MM.columns()) = pd.Cn_iM_LT;
  SharedMatrixNd L_iM_CnT_block = MM.block(pd.N_CONTACTS, MM.rows(), 0, pd.N_CONTACTS);
  MatrixNd::transpose(pd.Cn_iM_LT, L_iM_CnT_block);
  MM.block(pd.N_CONTACTS, MM.rows(), pd.N_CONTACTS, MM.columns()) = pd.L_iM_LT;
  qq.segment(0, pd.N_CONTACTS) = pd.Cn_v;
  qq.segment(pd.N_CONTACTS, qq.size()) = pd.L_v;

  FILE_LOG(LOG_SIMULATOR) << "# of constraints in determine_dq(): " << qq.size() << std::endl;

  // solve N*inv(M)*N'*dq = N*alpha for impulses 
  if (!_lcp.lcp_fast(MM, qq, z))
    _lcp.lcp_lemke_regularized(MM, qq, z);

  // update velocities
  update_from_stacked(z, pd);
  update_velocities(pd);

  // populate dq 
  for(unsigned i = 0; i< pd.super_bodies.size(); i++)
  {
    assert(body_index_map.find(pd.super_bodies[i]) != body_index_map.end());
    unsigned start = (body_index_map.find(pd.super_bodies[i]))->second;
    unsigned coord_num = pd.super_bodies[i]->num_generalized_coordinates(DynamicBody::eEuler);
    pd.super_bodies[i]->get_generalized_velocity(DynamicBody::eEuler, dq_sub);  
    dqm.segment(start, start+coord_num) = dq_sub;
  }
}

/// Updates determined impulses in UnilateralConstraintProblemData based on an LCP solution
void ConstraintStabilization::update_from_stacked(const VectorNd& z, UnilateralConstraintProblemData& pd)
{
  // update the problem data
  pd.update_from_stacked_qp(z);

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // save contact impulses
  for (unsigned i=0; i< pd.N_CONTACTS; i++)
  {
    // setup the contact frame
    P->q.set_identity();
    P->x = pd.contact_constraints[i]->contact_point;

    // setup the impulse in the contact frame
    Vector3d j;
    j = pd.contact_constraints[i]->contact_normal * pd.cn[i];
    j += pd.contact_constraints[i]->contact_tan1 * pd.cs[i];
    j += pd.contact_constraints[i]->contact_tan2 * pd.ct[i];

    // setup the spatial impulse
    SMomentumd jx(boost::const_pointer_cast<const Pose3d>(P));
    jx.set_linear(j);

    // transform the impulse to the global frame
    pd.contact_constraints[i]->contact_impulse += Pose3d::transform(GLOBAL, jx);
  }

  // save limit impulses
  for (unsigned i=0; i< pd.N_LIMITS; i++)
  {
    double limit_impulse = (pd.limit_constraints[i]->limit_upper) ? -pd.l[i] : pd.l[i];
    pd.limit_constraints[i]->limit_impulse += limit_impulse;
  }
}

/// Updates the body velocities based on problem data
void ConstraintStabilization::update_velocities(const UnilateralConstraintProblemData& pd)
{
  VectorNd v; 
  map<DynamicBodyPtr, VectorNd> gj;
  map<DynamicBodyPtr, VectorNd>::iterator gj_iter;

  // loop over all contact constraints first
  for (unsigned i=0; i< pd.contact_constraints.size(); i++)
  {
    // get the contact force
    const UnilateralConstraint& e = *pd.contact_constraints[i];
    SForced w(e.contact_impulse);

    // get the two single bodies of the contact
    SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
    SingleBodyPtr sb2 = e.contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr b1 = sb1->get_super_body();
    DynamicBodyPtr b2 = sb2->get_super_body();

    // convert force on first body to generalized forces
    if ((gj_iter = gj.find(b1)) == gj.end())
      b1->convert_to_generalized_force(sb1, w, gj[b1]);
    else
    {
      b1->convert_to_generalized_force(sb1, w, v);
      gj_iter->second += v;
    }

    // convert force on second body to generalized forces
    if ((gj_iter = gj.find(b2)) == gj.end())
      b2->convert_to_generalized_force(sb2, -w, gj[b2]);
    else
    {
      b2->convert_to_generalized_force(sb2, -w, v);
      gj_iter->second += v;
    }
  }

  // loop over all limit constraints next
  for (unsigned i=0; i< pd.limit_constraints.size(); i++)
  {
    const UnilateralConstraint& e = *pd.limit_constraints[i];
    ArticulatedBodyPtr ab = e.limit_joint->get_articulated_body();

    // get the iterator for the articulated body
    gj_iter = gj.find(ab);

    // apply limit impulses to bodies in independent coordinates
    if (dynamic_pointer_cast<RCArticulatedBody>(ab))
    {
      // get the index of the joint
      unsigned idx = e.limit_joint->get_coord_index() + e.limit_dof;

      // initialize the vector if necessary
      if (gj_iter == gj.end())
      {
        gj[ab].set_zero(ab->num_generalized_coordinates(DynamicBody::eSpatial));
        gj_iter = gj.find(ab);
      }

      // set the limit force
      gj_iter->second[idx] += e.limit_impulse;
    }
    else
    {
      // TODO: handle bodies in absolute coordinates here
      assert(false);
    }
  }

  // TODO: apply constraint impulses

  // apply all generalized impacts
  for (map<DynamicBodyPtr, VectorNd>::const_iterator i = gj.begin(); i != gj.end(); i++)
    i->first->apply_generalized_impulse(i->second);
}

/// Updates q doing a backtracking line search
/**
 * Two body version
 */
void ConstraintStabilization::update_q(const VectorNd& dq, VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  VectorNd qstar, grad;

  // copy the pairwise distances
  vector<PairwiseDistInfo>& pdi = sim->_pairwise_distances;
  vector<PairwiseDistInfo> pdi_old = pdi;

  // find the pairwise distances at q + dq
  qstar = dq;
  qstar += q;
  update_body_configurations(qstar, sim);
  sim->calc_pairwise_distances();

  // we may have to find roots for all pairwise distance functions
  vector<bool> bracket;
  for (unsigned i=0; i< pdi.size(); i++)
  {
    // look for penetrating and then not penetrating
    if (pdi_old[i].dist < 0.0 && pdi[i].dist > 0.0)
      bracket.push_back(true);
    else if (pdi_old[i].dist > 0.0 && pdi[i].dist < 0.0)
      bracket.push_back(true);
    else
      bracket.push_back(false); 
  }

  // NOTE: if there is no sign change at the endpoints, there may still be
  //       a sign change at some point in between, which means that after the
  //       root finding process, it is possible that two bodies could be
  //       interpenetrating that we didn't expect to be interpenetrating
  //       or that two bodies could be separated by an unexpected amount 

  // setup initial t
  double t = 1.0;

  // iterate over all brackets
  for (unsigned i=0; i< bracket.size(); i++)
  {
    // verify that we check this bracket 
    if (!bracket[i])
      continue;

    // call Ridder's method to determine new t
    double root = ridders(0, t, pdi_old[i].dist, pdi[i].dist, i, dq, q, sim);
    if (root > 0.0 && root < 1.0)
      t = std::min(root, t);
  }

  // all done? update q
  qstar = dq;
  qstar *= t;
  q += qstar;
}

/// Gets the body configurations, placing them into q 
void ConstraintStabilization::get_body_configurations(VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{  
  const std::vector<DynamicBodyPtr>& bodies = sim->_bodies;
  unsigned NGC = 0;

  // resize the vector appropriately
  BOOST_FOREACH(DynamicBodyPtr body, bodies)
    NGC += body->num_generalized_coordinates(DynamicBody::eEuler);
  q.resize(NGC);

  // set the appropriate part of the vector
  unsigned start = 0;
  BOOST_FOREACH(DynamicBodyPtr body, bodies){
    SharedVectorNd body_gcs = q.segment(start, start + body->num_generalized_coordinates(DynamicBody::eEuler));
    body->get_generalized_coordinates(DynamicBody::eEuler, body_gcs);
  }
}

/// Computes mapping from bodies to generalized coordinate indices 
void ConstraintStabilization::generate_body_index_map(std::map<DynamicBodyPtr, unsigned>& body_index_map, shared_ptr<ConstraintSimulator> sim)
{
  const std::vector<DynamicBodyPtr>& bodies = sim->_bodies;
  unsigned cur_index = 0;

  BOOST_FOREACH(DynamicBodyPtr body, bodies){
    body_index_map.insert(std::make_pair(body, cur_index));
    cur_index += body->num_generalized_coordinates(DynamicBody::eEuler);
  }
}

/// Updates the body configurations given q
void ConstraintStabilization::update_body_configurations(const VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  const std::vector<DynamicBodyPtr>& bodies = sim->_bodies;
  unsigned last = 0;
  BOOST_FOREACH(DynamicBodyPtr body, bodies){
    unsigned ngc = body->num_generalized_coordinates(DynamicBody::eEuler);
    Ravelin::SharedConstVectorNd gc_shared = q.segment(last,last+ngc);
    body->set_generalized_coordinates(DynamicBody::eEuler, gc_shared);
    last = ngc;
  }
}

/// sign function
static double sign(double x, double y)
{
  if (y > 0.0)
    return std::fabs(x);
  else 
    return -std::fabs(x);
}

/// Evaluates the function for root finding
double ConstraintStabilization::eval(double t, unsigned i, const VectorNd& dq, const VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  static VectorNd qstar;

  // setup qstar
  qstar = dq;
  qstar *= t;
  qstar += q; 

  // update body configurations
  update_body_configurations(qstar, sim);

  // compute new pairwise distance information
  sim->calc_pairwise_distances(); 

  return sim->_pairwise_distances[i].dist - eps;
}

/// Ridders method for root finding
double ConstraintStabilization::ridders(double x1, double x2, double fl, double fh, unsigned idx, const VectorNd& dq, const VectorNd& q, shared_ptr<ConstraintSimulator> sim)
{
  const unsigned MAX_ITERATIONS = 5;
  const double TOL = 1e-4;
  double ans=std::numeric_limits<double>::max(),fm,fnew,s,xh,xl,xm,xnew;

  if ((fl > 0.0 && fh < 0.0) || (fl < 0.0 && fh > 0.0)) 
  {
    xl=x1;
    xh=x2;
    ans = std::numeric_limits<double>::max();
    for (unsigned j=0;j< MAX_ITERATIONS;j++) 
    { 
      xm=0.5*(xl+xh);
      fm=eval(xm, idx, dq, q, sim); // First of two function evaluations per iteration
      s=std::sqrt(fm*fm-fl*fh); 
      if (s == 0.0) 
        return ans;
      xnew=xm+(xm-xl)*((fl >= fh ? 1.0 : -1.0)*fm/s); // Updating formula
      ans=xnew;
      fnew=eval(ans, idx, dq, q, sim);
      if (std::fabs(fnew) < TOL)
        return (fl < 0.0) ? xl : xh; 
      if (sign(fm,fnew) != fm) 
      {
        xl=xm;
        fl=fm;
        xh=ans;
        fh=fnew;
      } 
      else if (sign(fl,fnew) != fl) 
      { 
        xh=ans;
        fh=fnew;
      } 
      else if (sign(fh,fnew) != fh) 
      {
        xl=ans; // second of two function evaluations per iteration.
        fl=fnew; // Bookkeeping to keep the root bracketed on next iteration.
      } 
      else 
        assert(false);
    }

    // if here, maximum iterations have been exceeded
  }
  else 
  {
     if (fl == 0.0) 
       return x1;
     if (fh == 0.0) 
       return x2;

     throw std::runtime_error("Root was not bracketed");
  }

  return 0.0; // Never get here. 
}
