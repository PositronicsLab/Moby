#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <limits>
#include <set>
#include <cmath>
#include <numeric>

#include <Ravelin/LinAlgd.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/VectorNd.h>

#include <Moby/RCArticulatedBody.h>
#include <Moby/Constants.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/NumericalException.h>
#include <Moby/SustainedUnilateralConstraintSolveFailException.h>
#include <Moby/CompGeom.h>

#include <Moby/UnilateralConstraint.h>
#include <Moby/SustainedUnilateralConstraintHandler.h>
#include <Moby/SustainedUnilateralConstraintProblemData.h>

using namespace Ravelin;
using namespace Moby;
using std::pair;
using std::list;
using std::vector;
using std::map;
using std::multimap;
using std::set;
using std::endl;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;

/// Sets up the default parameters for the sustained unilateral handler
SustainedUnilateralConstraintHandler::SustainedUnilateralConstraintHandler(){}

// Processes sustained unilateral constraints
void SustainedUnilateralConstraintHandler::process_constraints(const vector<UnilateralConstraint>& constraints)
{
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::process_constraints() entered";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;

  // apply the method to all constraints
  if (!constraints.empty())
    apply_model(constraints);

  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::process_constraints() exited" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
}

/// Applies the model to a set of constraints
/**
 * \param constraints a set of constraints
 */
void SustainedUnilateralConstraintHandler::apply_model(const vector<UnilateralConstraint>& constraints)
{
  // **********************************************************
  // determine sets of connected constraints
  // **********************************************************
  list<list<UnilateralConstraint*> > groups;
  UnilateralConstraint::determine_connected_constraints(constraints, groups);
  UnilateralConstraint::remove_inactive_groups(groups);

  // **********************************************************
  // do method for each connected set
  // **********************************************************
  for (list<list<UnilateralConstraint*> >::iterator i = groups.begin(); i != groups.end(); i++)
  {
      // copy the list of constraints
      list<UnilateralConstraint*> rconstraints = *i;

      FILE_LOG(LOG_CONSTRAINT) << " -- pre-acceleration (all constraints: " << std::endl;
      for (list<UnilateralConstraint*>::iterator j = i->begin(); j != i->end(); j++)
        FILE_LOG(LOG_CONSTRAINT) << "    constraint: " << std::endl << **j;

      // determine a reduced set of constraints
      UnilateralConstraint::determine_minimal_set(rconstraints);

      // look to see whether all constraints have zero Coulomb friction
      bool all_frictionless = true;
      BOOST_FOREACH(UnilateralConstraint* e, rconstraints)
        if (e->constraint_type == UnilateralConstraint::eContact && e->contact_mu_coulomb > 0.0)
        {
          all_frictionless = false;
          break;
        }

      // apply model to the reduced constraints
      if (false && all_frictionless)
        apply_purely_viscous_model_to_connected_constraints(rconstraints);
      else
        apply_model_to_connected_constraints(rconstraints);
  }
}

/// Applies the Coulomb / viscous model to a set of connected constraints
/**
 * \param constraints a set of connected constraints
 */
void SustainedUnilateralConstraintHandler::apply_model_to_connected_constraints(const list<UnilateralConstraint*>& constraints)
{
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::apply_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // initialize constants and set easy to set constants
  _epd.contact_constraints.clear();
  _epd.limit_constraints.clear();
  _epd.N_CONTACTS = 0;
  _epd.N_LIMITS = 0;
  for (unsigned i=0; i< constraints.size(); i++)
  {
    if (_epd.constraints[i]->constraint_type == UnilateralConstraint::eContact)
    {
      _epd.N_CONTACTS++;
      _epd.contact_constraints.push_back(_epd.constraints[i]);
    }
    else if (_epd.constraints[i]->constraint_type == UnilateralConstraint::eLimit)
    {
      _epd.N_LIMITS++;
      _epd.limit_constraints.push_back(_epd.constraints[i]);
    }
  }

  // clear all constraint forces
  for (unsigned i=0; i< _epd.N_CONTACTS; i++)
    _epd.contact_constraints[i]->contact_impulse.set_zero(GLOBAL);
  for (unsigned i=0; i< _epd.N_LIMITS; i++)
    _epd.limit_constraints[i]->limit_impulse = 0.0;

  // compute sliding velocities
  _cs_visc.resize(_epd.N_CONTACTS);
  RowIteratord cs_visc_iter = _cs_visc.row_iterator_begin();
  for (unsigned i=0; i< _epd.N_CONTACTS; i++)
    _cs_visc[i] = _epd.contact_constraints[i]->calc_contact_vel(_epd.contact_constraints[i]->contact_tan1);

  // compute viscous friction terms
  cs_visc_iter = _cs_visc.row_iterator_begin();
  for (unsigned i=0; i< _epd.N_CONTACTS; i++, cs_visc_iter++)
    (*cs_visc_iter) *= _epd.contact_constraints[i]->contact_mu_viscous;

  // add in viscous friction forces and recompute dynamics
  _epd.cs = _cs_visc;
  bool zero_visc_force = (_epd.cs.norm_inf() < NEAR_ZERO);

  // recompute system dynamics, if necessary
  if (!zero_visc_force)
  {
    // setup a temporary frame
    shared_ptr<Pose3d> P(new Pose3d);

    // save normal contact impulses
    for (unsigned i=0; i< _epd.contact_constraints.size(); i++)
    {
      // verify that the constraint type is a contact
      assert(_epd.contact_constraints[i]->constraint_type == UnilateralConstraint::eContact);

      // setup the contact frame
      P->q.set_identity();
      P->x = _epd.contact_constraints[i]->contact_point;

      // setup the impulse in the contact frame
      Vector3d f = _epd.contact_constraints[i]->contact_tan1 * _epd.cs[i];

      // setup the spatial force
      SForced fx(boost::const_pointer_cast<const Pose3d>(P));
      fx.set_force(f);

      // transform the impulse to the global frame
      _epd.constraints[i]->contact_impulse = Pose3d::transform(GLOBAL, fx);
    }

    // apply contact forces and recompute dynamics
    apply_forces(_epd);

    BOOST_FOREACH(DynamicBodyPtr db, _epd.super_bodies)
      db->calc_fwd_dyn();
  }

  // compute all LCP problem data
  compute_problem_data(_epd);

  // solve the linear complementarity problem 
  VectorNd z;
  if (!solve_coulomb_lcp(_epd, z))
    throw SustainedUnilateralConstraintSolveFailException();

  FILE_LOG(LOG_CONSTRAINT) << "Resting constraint forces : " << z << std::endl;

  // apply FORCES
  apply_forces(_epd);

  // calculate problem data again (for debugging)
  BOOST_FOREACH(DynamicBodyPtr db, _epd.super_bodies)
    db->calc_fwd_dyn();
  compute_problem_data(_epd);
  FILE_LOG(LOG_CONSTRAINT) << "(problem data double-check for debugging) " << _epd.Cn_a << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "Cn * a: " << _epd.Cn_a << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "Cs * a: " << _epd.Cs_a << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "Ct * a: " << _epd.Ct_a << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "L * a: "  << _epd.L_a << std::endl; 

  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::apply_model_to_connected_constraints() exiting" << endl;
}

/// Applies the purely viscous model to a set of connected constraints
/**
 * \param constraints a set of connected constraints
 */
void SustainedUnilateralConstraintHandler::apply_purely_viscous_model_to_connected_constraints(const list<UnilateralConstraint*>& constraints)
{
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::apply_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // initialize constants and set easy to set constants
  _epd.contact_constraints.clear();
  _epd.limit_constraints.clear();
  _epd.N_CONTACTS = 0;
  _epd.N_LIMITS = 0;
  for (unsigned i=0; i< constraints.size(); i++)
  {
    if (_epd.constraints[i]->constraint_type == UnilateralConstraint::eContact)
    {
      _epd.N_CONTACTS++;
      _epd.contact_constraints.push_back(_epd.constraints[i]);
    }
    else if (_epd.constraints[i]->constraint_type == UnilateralConstraint::eLimit)
    {
      _epd.N_LIMITS++;
      _epd.limit_constraints.push_back(_epd.constraints[i]);
    }
  }

  // clear all constraint forces
  for (unsigned i=0; i< _epd.N_CONTACTS; i++)
    _epd.contact_constraints[i]->contact_impulse.set_zero(GLOBAL);
  for (unsigned i=0; i< _epd.N_LIMITS; i++)
    _epd.limit_constraints[i]->limit_impulse = 0.0;

  // compute sliding velocities
  _cs_visc.resize(_epd.N_CONTACTS);
  RowIteratord cs_visc_iter = _cs_visc.row_iterator_begin();
  for (unsigned i=0; i< _epd.N_CONTACTS; i++)
    _cs_visc[i] = _epd.contact_constraints[i]->calc_contact_vel(_epd.contact_constraints[i]->contact_tan1);

  // compute viscous friction terms
  cs_visc_iter = _cs_visc.row_iterator_begin();
  for (unsigned i=0; i< _epd.N_CONTACTS; i++, cs_visc_iter++)
    (*cs_visc_iter) *= _epd.contact_constraints[i]->contact_mu_viscous;

  // add in viscous friction forces and recompute dynamics
  _epd.cs = _cs_visc;
  bool zero_visc_force = (_epd.cs.norm_inf() < NEAR_ZERO);

  // recompute system dynamics, if necessary
  if (!zero_visc_force)
  {
    // setup a temporary frame
    shared_ptr<Pose3d> P(new Pose3d);

    // save normal contact impulses
    for (unsigned i=0; i< _epd.constraints.size(); i++)
    {
      // verify that the constraint type is a contact
      assert(_epd.constraints[i]->constraint_type == UnilateralConstraint::eContact);

      // setup the contact frame
      P->q.set_identity();
      P->x = _epd.constraints[i]->contact_point;

      // setup the impulse in the contact frame
      Vector3d f = _epd.constraints[i]->contact_tan1 * _epd.cs[i];

      // setup the spatial force
      SForced fx(boost::const_pointer_cast<const Pose3d>(P));
      fx.set_force(f);

      // transform the impulse to the global frame
      _epd.constraints[i]->contact_impulse = Pose3d::transform(GLOBAL, fx);
    }

    // apply contact forces and recompute dynamics
    apply_forces(_epd);

    BOOST_FOREACH(DynamicBodyPtr db, _epd.super_bodies)
      db->calc_fwd_dyn();
  }

  // compute problem data
  compute_problem_data(_epd);

  // solve the (non-frictional) linear complementarity problem to determine
  // the kappa constant
  VectorNd z;
  if (!solve_purely_viscous_lcp(_epd, z))
    throw SustainedUnilateralConstraintSolveFailException();

  FILE_LOG(LOG_CONSTRAINT) << "Resting constraint forces : " << z << std::endl;

  // apply FORCES
  apply_forces(_epd);

  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::apply_model_to_connected_constraints() exiting" << endl;
}

/// Applies resting constraint forces to bodies and saves the generalized forces
void SustainedUnilateralConstraintHandler::apply_forces(const SustainedUnilateralConstraintProblemData& q)
{
  map<DynamicBodyPtr, VectorNd> gj;
  map<DynamicBodyPtr, VectorNd>::iterator gj_iter;

  // loop over all contact contacts first
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    // get the contact force
    const UnilateralConstraint& c = *q.contact_constraints[i];
    SForced w(c.contact_impulse);

    // get the two single bodies of the contact
    SingleBodyPtr sb1 = c.contact_geom1->get_single_body();
    SingleBodyPtr sb2 = c.contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr b1 = sb1->get_super_body();
    DynamicBodyPtr b2 = sb2->get_super_body();

    // convert force on first body to generalized forces
    if ((gj_iter = gj.find(b1)) == gj.end())
      b1->convert_to_generalized_force(sb1, w, gj[b1]);
    else
    {
      b1->convert_to_generalized_force(sb1, w, _workv);
      gj_iter->second += _workv;
    }

    // convert force on second body to generalized forces
    if ((gj_iter = gj.find(b2)) == gj.end())
      b2->convert_to_generalized_force(sb2, -w, gj[b2]);
    else
    {
      b2->convert_to_generalized_force(sb2, -w, _workv);
      gj_iter->second += _workv;
    }
  }

  // loop over all limit constraints next
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    const UnilateralConstraint& e = *q.limit_constraints[i];
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

  // apply all generalized forces
  for (map<DynamicBodyPtr, VectorNd>::const_iterator i = gj.begin(); i != gj.end(); i++)
  {
    // apply the force
    i->first->add_generalized_force(i->second);
  }
}

/// Computes the data to the LCP / QP problems
void SustainedUnilateralConstraintHandler::compute_problem_data(SustainedUnilateralConstraintProblemData& q)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  SAFESTATIC MatrixNd workM;
  SAFESTATIC VectorNd workv;

// use the newer method - it works for sliding friction
compute_problem_data2(q);
return;

  // determine set of "super" bodies from the constraints
  q.super_bodies.clear();
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    q.super_bodies.push_back(get_super_body(q.contact_constraints[i]->contact_geom1->get_single_body()));
    q.super_bodies.push_back(get_super_body(q.contact_constraints[i]->contact_geom2->get_single_body()));
  }

  // make super bodies vector unique
  std::sort(q.super_bodies.begin(), q.super_bodies.end());
  q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

  // setup contact working set
  q.contact_working_set.clear();
  q.contact_working_set.resize(q.N_CONTACTS, true);

  // compute number of friction polygon edges
  q.N_STICKING = 0;
  for (unsigned i=0; i<  q.contact_constraints.size(); i++)
  {
    if (q.contact_constraints[i]->get_friction_type() == UnilateralConstraint::eSticking)
    {
      q.N_STICKING++;
      if ( q.contact_constraints[i]->contact_NK < UINF)
        q.N_K_TOTAL +=  q.contact_constraints[i]->contact_NK;
    }
    else if ( q.contact_constraints[i]->contact_NK == UINF)
      break;
  }

  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::compute_problem_data(.)" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  " << q.N_STICKING << " sticking contacts" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  " << (q.N_CONTACTS - q.N_STICKING) << " sliding contacts" << std::endl;

  // initialize the problem matrices / vectors
  q.Cn_iM_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_iM_CsT.set_zero(q.N_CONTACTS, q.N_STICKING);
  q.Cn_iM_CtT.set_zero(q.N_CONTACTS, q.N_STICKING);
  q.Cs_iM_CnT.set_zero(q.N_STICKING, q.N_CONTACTS);
  q.Cs_iM_CsT.set_zero(q.N_STICKING, q.N_STICKING);
  q.Cs_iM_CtT.set_zero(q.N_STICKING, q.N_STICKING);
  q.Ct_iM_CnT.set_zero(q.N_STICKING, q.N_CONTACTS);
  q.Ct_iM_CsT.set_zero(q.N_STICKING, q.N_STICKING);
  q.Ct_iM_CtT.set_zero(q.N_STICKING, q.N_STICKING);
  q.Cn_a.set_zero(q.N_CONTACTS);
  q.Cs_a.set_zero(q.N_STICKING);
  q.Ct_a.set_zero(q.N_STICKING);
  q.cn.set_zero(q.N_CONTACTS);
  q.cs.set_zero(q.N_CONTACTS);
  q.ct.set_zero(q.N_CONTACTS);

  // setup indices
  q.CN_IDX = 0;
  q.CS_IDX = q.CN_IDX + q.N_CONTACTS;
  q.CT_IDX = q.CS_IDX + q.N_STICKING;
  q.NCS_IDX = q.CT_IDX + q.N_STICKING;
  q.NCT_IDX = q.NCS_IDX + q.N_STICKING;
  // TODO: add constraint computation and cross computation methods to Joint

  // setup the number of variables
  q.N_VARS = q.NCT_IDX + q.N_STICKING;

  // get iterators to the proper matrices
  RowIteratord CnCn = q.Cn_iM_CnT.row_iterator_begin();
  RowIteratord CnCs = q.Cn_iM_CsT.row_iterator_begin();
  RowIteratord CnCt = q.Cn_iM_CtT.row_iterator_begin();
  RowIteratord CsCn = q.Cs_iM_CnT.row_iterator_begin();
  RowIteratord CsCs = q.Cs_iM_CsT.row_iterator_begin();
  RowIteratord CsCt = q.Cs_iM_CtT.row_iterator_begin();
  RowIteratord CtCn = q.Ct_iM_CnT.row_iterator_begin();
  RowIteratord CtCs = q.Ct_iM_CsT.row_iterator_begin();
  RowIteratord CtCt = q.Ct_iM_CtT.row_iterator_begin();

  // process contact constraints, setting up matrices
  for (unsigned i=0, k=0; i<  q.contact_constraints.size(); i++)
  {
    const UnilateralConstraint* ci =  q.contact_constraints[i];
    const unsigned ROWS = (ci->get_friction_type() == UnilateralConstraint::eSticking) ? 3 : 1;

    // this approach only works for contact constraints at the moment
    assert(ci->constraint_type == UnilateralConstraint::eContact);

    // compute cross constraint data for contact constraints
    for (unsigned j=0; j<  q.contact_constraints.size(); j++)
    {
      const UnilateralConstraint* cj =  q.contact_constraints[j];
      const unsigned COLS = (cj->get_friction_type() == UnilateralConstraint::eSticking) ? 3 : 1;

      // reset workM
      workM.set_zero(ROWS, COLS);

      // check whether i==j (single contact constraint)
      if (i == j)
      {
        // compute matrix / vector for contact constraint i
        workv.set_zero(ROWS);
         q.contact_constraints[i]->compute_constraint_data(workM, workv);

        if (ROWS == 3)
        {
          // setup appropriate parts of contact inertia matrices
          RowIteratord_const data = workM.row_iterator_begin();
          *CnCn = *data++;  *CnCs = *data++;  *CnCt = *data++;
          *CsCn = *data++;  *CsCs = *data++;  *CsCt = *data++;
          *CtCn = *data++;  *CtCs = *data++;  *CtCt = *data;

          // setup appropriate parts of contact velocities
          data = workv.row_iterator_begin();
          q.Cn_a[i] = *data++;
          q.Cs_a[k] = *data++;
          q.Ct_a[k] = *data;

          // update k (NOTE: we need k b/c some contacts may be slipping)
          k++;
        }
        else
        {
          // setup appropriate part of contact inertia matrices
          *CnCn = *workM.row_iterator_begin();

          // setup appropriate part of contact velocities
          q.Cn_a[i] = *workv.row_iterator_begin();
        }
      }
      else
      {
        // compute matrix for cross constraint
         q.contact_constraints[i]->compute_cross_constraint_data(* q.contact_constraints[j], workM);

        if (ROWS == 3)
        {
          if (COLS == 3)
          {
            // setup appropriate parts of contact inertia matrices
            RowIteratord_const data = workM.row_iterator_begin();
            *CnCn = *data++;  *CnCs = *data++;  *CnCt = *data++;
            *CsCn = *data++;  *CsCs = *data++;  *CsCt = *data++;
            *CtCn = *data++;  *CtCs = *data++;  *CtCt = *data;
          }
          else  // columns = 1
          {
            // setup appropriate parts of contact inertia matrices
            ColumnIteratord_const data = workM.column_iterator_begin();
            *CnCn = *data++;
            *CsCn = *data++;
            *CtCn = *data;
          }
        }
        else
        {
          if (COLS == 3)
          {
            // setup appropriate parts of contact inertia matrices
            RowIteratord_const data = workM.row_iterator_begin();
            *CnCn = *data++;
            *CnCs = *data++;
            *CnCt = *data;
          }
          else  // rows = 1 and columns = 1
          {
            // setup appropriate parts of contact inertia matrices
            *CnCn = *workM.column_iterator_begin();
          }
        }
      }

      // advance appropriate iterators
      if (ROWS == 3)
      {
        if (COLS == 3)
        {
          CnCn++;  CnCs++;  CnCt++;
          CsCn++;  CsCs++;  CsCt++;
          CtCn++;  CtCs++;  CtCt++;
        }
        else
        {
          CnCn++;
          CsCn++;
          CtCn++;
        }
      }
      else
      {
        if (COLS == 3)
        {
          CnCn++;  CnCs++;  CnCt++;
        }
        else
        {
          CnCn++;
        }
      }
    }
  }
}

/// Computes the data to the LCP / QP problems -- second approach
void SustainedUnilateralConstraintHandler::compute_problem_data2(SustainedUnilateralConstraintProblemData& q)
{
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  VectorNd m_l_v, m_cn_v, m_cs_v, m_ct_v;
  VectorNd p_l_v, p_cn_v, p_cs_v, p_ct_v;

  // determine set of "super" bodies from constraints
  q.super_bodies.clear();

  // loop over all contact constraints first
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    const UnilateralConstraint& e = *q.contact_constraints[i];

    // get the two single bodies of the contact
    SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
    SingleBodyPtr sb2 = e.contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr b1 = sb1->get_super_body();
    DynamicBodyPtr b2 = sb2->get_super_body();

    // add the super bodies
    q.super_bodies.push_back(b1);
    q.super_bodies.push_back(b2);
  }

  // loop over all limit constraints next
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    const UnilateralConstraint& e = *q.limit_constraints[i];
    ArticulatedBodyPtr ab = e.limit_joint->get_articulated_body();
    q.super_bodies.push_back(ab);
  }

  // make super bodies vector unique
  std::sort(q.super_bodies.begin(), q.super_bodies.end());
  q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

  // setup contact working set
  q.contact_working_set.clear();
  q.contact_working_set.resize(q.N_CONTACTS, true);

  // compute number of friction polygon edges
  q.N_STICKING = 0;
  for (unsigned i=0; i<  q.contact_constraints.size(); i++)
  {
    if (q.contact_constraints[i]->constraint_type == UnilateralConstraint::eContact)
    {
      // update the number of sticking contacts and number of friction edges
      if (q.contact_constraints[i]->get_friction_type() == UnilateralConstraint::eSticking)
      {
        q.N_STICKING++;
        q.N_K_TOTAL +=  q.contact_constraints[i]->contact_NK;
      }
    }
    else if ( q.contact_constraints[i]->contact_NK == UINF)
      throw std::runtime_error("Infinite friction cone encountered for sustained unilateral contact!");
  }

  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::compute_problem_data(.)" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  " << q.N_STICKING << " sticking contacts" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  " << (q.N_CONTACTS - q.N_STICKING) << " sliding contacts" << std::endl;

  // initialize the problem matrices / vectors
  q.Cn_iM_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
  q.Cn_iM_CsT.set_zero(q.N_CONTACTS, q.N_STICKING);
  q.Cn_iM_CtT.set_zero(q.N_CONTACTS, q.N_STICKING);
  q.Cn_iM_LT.set_zero(q.N_CONTACTS, q.N_LIMITS);
  q.Cs_iM_CnT.set_zero(q.N_STICKING, q.N_CONTACTS);
  q.Cs_iM_CsT.set_zero(q.N_STICKING, q.N_STICKING);
  q.Cs_iM_CtT.set_zero(q.N_STICKING, q.N_STICKING);
  q.Cs_iM_LT.set_zero(q.N_STICKING, q.N_LIMITS);
  q.Ct_iM_CnT.set_zero(q.N_STICKING, q.N_CONTACTS);
  q.Ct_iM_CsT.set_zero(q.N_STICKING, q.N_STICKING);
  q.Ct_iM_CtT.set_zero(q.N_STICKING, q.N_STICKING);
  q.Ct_iM_LT.set_zero(q.N_STICKING, q.N_LIMITS);
  q.L_iM_CnT.set_zero(q.N_LIMITS, q.N_CONTACTS);
  q.L_iM_CsT.set_zero(q.N_LIMITS, q.N_STICKING);
  q.L_iM_CtT.set_zero(q.N_LIMITS, q.N_STICKING);
  q.L_iM_LT.set_zero(q.N_LIMITS, q.N_LIMITS);
  q.Cn_a.set_zero(q.N_CONTACTS);
  q.Cs_a.set_zero(q.N_STICKING);
  q.Ct_a.set_zero(q.N_STICKING);
  q.L_a.set_zero(q.N_LIMITS);
  q.cn.set_zero(q.N_CONTACTS);
  q.cs.set_zero(q.N_CONTACTS);
  q.ct.set_zero(q.N_CONTACTS);
  q.l.set_zero(q.N_LIMITS);

  // setup indices
  q.CN_IDX = 0;
  q.CS_IDX = q.CN_IDX + q.N_CONTACTS;
  q.CT_IDX = q.CS_IDX + q.N_STICKING;
  q.NCS_IDX = q.CT_IDX + q.N_STICKING;
  q.NCT_IDX = q.NCS_IDX + q.N_STICKING;
  // TODO: add constraint computation and cross computation methods to Joint

  // setup the number of variables
  q.N_VARS = q.NCT_IDX + q.N_STICKING;

  // save the velocities and forces of all bodies in contacts
  std::map<DynamicBodyPtr, VectorNd> saved_velocities, saved_forces;
  for (unsigned i=0; i< q.super_bodies.size(); i++)
  {
    DynamicBodyPtr db = q.super_bodies[i];
    db->get_generalized_velocity(DynamicBody::eSpatial, saved_velocities[db]);
    db->get_generalized_forces(saved_forces[db]);
  }

  // process contact constraints, setting up vectors
  for (unsigned i=0, k=0; i<  q.contact_constraints.size(); i++)
  {
    const UnilateralConstraint* ci =  q.contact_constraints[i];
    const unsigned ROWS = (ci->get_friction_type() == UnilateralConstraint::eSticking) ? 3 : 1;

    if (ROWS == 3)
    {
      // setup appropriate parts of contact accelerations 
      q.Cn_a[i] = q.contact_constraints[i]->calc_contact_accel(q.contact_constraints[i]->contact_normal) + q.contact_constraints[i]->calc_contact_vel(q.contact_constraints[i]->contact_normal_dot);
      q.Cs_a[k] = q.contact_constraints[i]->calc_contact_accel(q.contact_constraints[i]->contact_tan1) + q.contact_constraints[i]->calc_contact_vel(q.contact_constraints[i]->contact_tan1_dot);
      q.Ct_a[k] = q.contact_constraints[i]->calc_contact_accel(q.contact_constraints[i]->contact_tan2) + q.contact_constraints[i]->calc_contact_vel(q.contact_constraints[i]->contact_tan2_dot);

      // update k (NOTE: we need k b/c some contacts may be slipping)
      k++;
    }
    else
    {
      // setup appropriate part of contact accelerations 
      q.Cn_a[i] = q.contact_constraints[i]->calc_contact_accel(q.contact_constraints[i]->contact_normal) + q.contact_constraints[i]->calc_contact_vel(q.contact_constraints[i]->contact_normal_dot);
    }
  }

  // process limit constraints, setting up vectors
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    const UnilateralConstraint* ci =  q.limit_constraints[i];
    q.L_a[i] = q.limit_constraints[i]->calc_constraint_accel();
  }

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // compute contribution in limit direction
  // NOTE: i is contact index, k is sticking friction index
  for (unsigned i=0; i<  q.limit_constraints.size(); i++)
  {
    // get the i'th limit 
    const UnilateralConstraint* ci =  q.limit_constraints[i];

    // measure the constraint velocities
    calc_constraint_velocities(q, m_cn_v, m_cs_v, m_ct_v, m_l_v);

    // get the super body for the constraint
    ArticulatedBodyPtr super = ci->limit_joint->get_articulated_body();

    // apply the impulse along the joint limit
    unsigned index = ci->limit_joint->get_coord_index() + ci->limit_dof;
    _workv.resize(super->num_generalized_coordinates(DynamicBody::eSpatial));
    _workv.set_zero();
    _workv[index] = (ci->limit_upper) ? -1.0 : 1.0; 
    super->apply_generalized_impulse(_workv);

    // measure the constraint velocities again
    calc_constraint_velocities(q, p_cn_v, p_cs_v, p_ct_v, p_l_v);

    // compute the difference
    p_cn_v -= m_cn_v;
    p_cs_v -= m_cs_v;
    p_ct_v -= m_ct_v;
    p_l_v -= m_l_v;

    // update the appropriate column of the matrix
    q.Cn_iM_LT.column(i) = p_cn_v;    
    q.Cs_iM_LT.column(i) = p_cs_v;    
    q.Ct_iM_LT.column(i) = p_ct_v;    
    q.L_iM_LT.column(i) = p_l_v;
  }  

  // compute contribution in normal direction
  // NOTE: i is contact index
  for (unsigned i=0; i<  q.contact_constraints.size(); i++)
  {
    // get the i'th contact
    const UnilateralConstraint* ci =  q.contact_constraints[i];

    // measure the constraint velocities
    calc_constraint_velocities(q, m_cn_v, m_cs_v, m_ct_v, m_l_v);

    // setup the contact frame
    P->q.set_identity();
    P->x = ci->contact_point;

    // get the two single bodies of the contact
    SingleBodyPtr sb1 = ci->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = ci->contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr b1 = sb1->get_super_body();
    DynamicBodyPtr b2 = sb2->get_super_body();

    // apply impulse(s) in the normal direction
    Vector3d f = q.contact_constraints[i]->contact_normal;
    if (ci->get_friction_type() != UnilateralConstraint::eSticking)
      f -= q.contact_constraints[i]->contact_tan1 * ci->contact_mu_coulomb;

    // setup the spatial force
    SForced fx(boost::const_pointer_cast<const Pose3d>(P));
    fx.set_force(f);

    // convert force on bodies to generalized forces
    VectorNd f1, f2;
    b1->convert_to_generalized_force(sb1, fx, f1);
    b2->convert_to_generalized_force(sb2, -fx, f2);

    // apply the impulse
    b1->apply_generalized_impulse(f1);
    b2->apply_generalized_impulse(f2);

    // measure the constraint velocities again
    calc_constraint_velocities(q, p_cn_v, p_cs_v, p_ct_v, p_l_v);

    // compute the difference
    p_cn_v -= m_cn_v;
    p_cs_v -= m_cs_v;
    p_ct_v -= m_ct_v;
    p_l_v -= m_l_v;

    // update the appropriate column of the matrix
    q.Cn_iM_CnT.column(i) = p_cn_v;    
    q.Cs_iM_CnT.column(i) = p_cs_v;    
    q.Ct_iM_CnT.column(i) = p_ct_v;    
    q.L_iM_CnT.column(i) = p_l_v;
  }

  // compute contribution in first sticking friction direction
  // NOTE: i is contact index, k is sticking friction index
  for (unsigned i=0, k=0; i<  q.contact_constraints.size(); i++)
  {
    // get the i'th contact
    const UnilateralConstraint* ci =  q.contact_constraints[i];

    // if this is not a sticking friction constraint, skip it
    if (q.contact_constraints[i]->get_friction_type() != UnilateralConstraint::eSticking)
      continue;

    // measure the constraint velocities
    calc_constraint_velocities(q, m_cn_v, m_cs_v, m_ct_v, m_l_v);

    // setup the contact frame
    P->q.set_identity();
    P->x = ci->contact_point;

    // get the two single bodies of the contact
    SingleBodyPtr sb1 = ci->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = ci->contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr b1 = sb1->get_super_body();
    DynamicBodyPtr b2 = sb2->get_super_body();

    // apply impulse(s) in the tangent direction
    Vector3d f = q.contact_constraints[i]->contact_tan1;

    // setup the spatial force
    SForced fx(boost::const_pointer_cast<const Pose3d>(P));
    fx.set_force(f);

    // convert force on bodies to generalized forces
    VectorNd f1, f2;
    b1->convert_to_generalized_force(sb1, fx, f1);
    b2->convert_to_generalized_force(sb2, -fx, f2);

    // apply the impulse
    b1->apply_generalized_impulse(f1);
    b2->apply_generalized_impulse(f2);

    // measure the constraint velocities again
    calc_constraint_velocities(q, p_cn_v, p_cs_v, p_ct_v, p_l_v);

    // compute the difference
    p_cn_v -= m_cn_v;
    p_cs_v -= m_cs_v;
    p_ct_v -= m_ct_v;
    p_l_v -= m_l_v;

    // update the appropriate column of the matrix
    q.Cn_iM_CsT.column(k) = p_cn_v;    
    q.Cs_iM_CsT.column(k) = p_cs_v;    
    q.Ct_iM_CsT.column(k) = p_ct_v;    
    q.L_iM_CsT.column(k) = p_l_v;

    // update k - the sticking friction index
    k++;
  }

  // compute contribution in second sticking friction direction
  // NOTE: i is contact index, k is sticking friction index
  for (unsigned i=0, k=0; i<  q.contact_constraints.size(); i++)
  {
    // get the i'th contact
    const UnilateralConstraint* ci =  q.contact_constraints[i];

    // if this is not a sticking friction constraint, skip it
    if (q.contact_constraints[i]->get_friction_type() != UnilateralConstraint::eSticking)
      continue;

    // measure the constraint velocities
    calc_constraint_velocities(q, m_cn_v, m_cs_v, m_ct_v, m_l_v);

    // setup the contact frame
    P->q.set_identity();
    P->x = ci->contact_point;

    // get the two single bodies of the contact
    SingleBodyPtr sb1 = ci->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = ci->contact_geom2->get_single_body();

    // get the two super bodies
    DynamicBodyPtr b1 = sb1->get_super_body();
    DynamicBodyPtr b2 = sb2->get_super_body();

    // apply impulse(s) in the tangent direction
    Vector3d f = q.contact_constraints[i]->contact_tan2;

    // setup the spatial force
    SForced fx(boost::const_pointer_cast<const Pose3d>(P));
    fx.set_force(f);

    // convert force on bodies to generalized forces
    VectorNd f1, f2;
    b1->convert_to_generalized_force(sb1, fx, f1);
    b2->convert_to_generalized_force(sb2, -fx, f2);

    // apply the impulse
    b1->apply_generalized_impulse(f1);
    b2->apply_generalized_impulse(f2);

    // measure the constraint velocities again
    calc_constraint_velocities(q, p_cn_v, p_cs_v, p_ct_v, p_l_v);

    // compute the difference
    p_cn_v -= m_cn_v;
    p_cs_v -= m_cs_v;
    p_ct_v -= m_ct_v;
    p_l_v -= m_l_v;

    // update the appropriate column of the matrix
    q.Cn_iM_CtT.column(k) = p_cn_v;    
    q.Cs_iM_CtT.column(k) = p_cs_v;    
    q.Ct_iM_CtT.column(k) = p_ct_v;    
    q.L_iM_CtT.column(k) = p_l_v;

    // update k - the sticking friction index
    k++;
  }

  // restore the velocities and apply forces to all bodies in contacts
  for (std::map<DynamicBodyPtr, VectorNd>::const_iterator i = saved_velocities.begin(); i != saved_velocities.end(); i++)
    i->first->set_generalized_velocity(DynamicBody::eSpatial, i->second);
  for (std::map<DynamicBodyPtr, VectorNd>::const_iterator i = saved_forces.begin(); i != saved_forces.end(); i++)
  {
    i->first->reset_accumulators();
    i->first->add_generalized_force(i->second);
  }

for (unsigned i=0; i< q.super_bodies.size(); i++)
  FILE_LOG(LOG_CONSTRAINT) << "generalized force on: " << q.super_bodies[i]->id << " " << saved_forces[q.super_bodies[i]] << std::endl;
}

// Calculates constraint velocities (for compute_problem_data2(.))
void SustainedUnilateralConstraintHandler::calc_constraint_velocities(SustainedUnilateralConstraintProblemData& q, VectorNd& cn_v, VectorNd& cs_v, VectorNd& ct_v, VectorNd& l_v)
{
  // setup the size of the vectors
  cn_v.resize(q.N_CONTACTS);
  cs_v.resize(q.N_STICKING);
  ct_v.resize(q.N_STICKING);
  l_v.resize(q.N_LIMITS);

  // measure the velocities along the contact directions
  for (unsigned j=0, r=0; j< q.contact_constraints.size(); j++)
  {
    cn_v[j] = q.contact_constraints[j]->calc_contact_vel(q.contact_constraints[j]->contact_normal);
    if (q.contact_constraints[j]->get_friction_type() == UnilateralConstraint::eSticking)
    {
      cs_v[r] = q.contact_constraints[j]->calc_contact_vel(q.contact_constraints[j]->contact_tan1);
      ct_v[r] = q.contact_constraints[j]->calc_contact_vel(q.contact_constraints[j]->contact_tan2);
      r++;
    }
  }

  // measure the velocities along the limit directions
  for (unsigned j=0; j< q.limit_constraints.size(); j++)
    l_v[j] = q.limit_constraints[j]->calc_constraint_vel();
}

// Solves the sustained constraint LCP w/Coulomb and viscous friction
bool SustainedUnilateralConstraintHandler::solve_coulomb_lcp(SustainedUnilateralConstraintProblemData& q, VectorNd& z)
{
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::solve_coulomb_lcp() entered" << std::endl;

  unsigned N_LIMIT = q.N_LIMITS;
  unsigned NK_DIRS = 0;
  for(unsigned i=0,j=0,r=0;i<q.N_CONTACTS;i++)
    if(q.contact_constraints[i]->get_friction_type() == UnilateralConstraint::eSticking)
    {
      if (q.contact_constraints[i]->contact_NK > 4)
        NK_DIRS+=(q.contact_constraints[i]->contact_NK+4)/4;
      else
        NK_DIRS+=1;
    }
  unsigned N_FRICT = q.N_CONTACTS+q.N_STICKING*4;
  unsigned N_CONST = N_FRICT + N_LIMIT;

  // setup sizes
  _UL.set_zero(N_CONST, N_CONST);
  _UR.set_zero(N_CONST, NK_DIRS);
  _LL.set_zero(NK_DIRS, N_CONST);
  _MM.set_zero(_UL.rows() + _LL.rows(), _UL.columns() + _UR.columns());

  // now do upper right hand block of LCP matrix
  /*     n          r          r           r           r
  n  Cn_iM_CnT  Cn_iM_CsT  -Cn_iM_CsT   Cn_iM_CtT  -Cn_iM_CtT
  r  Cs_iM_CnT  Cs_iM_CsT  -Cs_iM_CsT   Cs_iM_CtT  -Cs_iM_CtT
  r -Cs_iM_CnT -Cs_iM_CsT   Cs_iM_CsT  -Cs_iM_CtT   Cs_iM_CtT
  r  Ct_iM_CnT  Ct_iM_CsT  -Ct_iM_CsT   Ct_iM_CtT  -Ct_iM_CtT
  r -Ct_iM_CnT -Ct_iM_CsT   Ct_iM_CsT  -Ct_iM_CtT   Ct_iM_CtT
  // Set positive submatrices
         n          r          r           r           r
  n  Cn_iM_CnT  Cn_iM_CsT               Cn_iM_CtT
  r  Cs_iM_CnT  Cs_iM_CsT               Cs_iM_CtT
  r                         Cs_iM_CsT               Cs_iM_CtT
  r  Ct_iM_CnT  Ct_iM_CsT               Ct_iM_CtT
  r                         Ct_iM_CsT               Ct_iM_CtT
  */
  _UL.set_sub_mat(0,0,q.Cn_iM_CnT);
  _UL.set_sub_mat(N_FRICT,N_FRICT,q.L_iM_LT);
  _UL.set_sub_mat(N_FRICT,0,q.L_iM_CnT);
  _UL.set_sub_mat(0,N_FRICT,q.Cn_iM_LT);
  // setup the LCP matrix

  // setup the LCP vector
  _qq.set_zero(_MM.rows());
  _qq.set_sub_vec(0,q.Cn_a);
  _qq.set_sub_vec(N_FRICT,q.L_a);

  FILE_LOG(LOG_CONSTRAINT) << "Cn*inv(M)*Cn': " << std::endl << q.Cn_iM_CnT;

  FILE_LOG(LOG_CONSTRAINT) << "Cn*inv(M)*L': " << std::endl << q.Cn_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "L*inv(M)*Cn': " << std::endl << q.L_iM_CnT;


  if(q.N_STICKING > 0){

    _UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS,q.Cs_iM_CsT);
    _UL.set_sub_mat(q.N_CONTACTS,0,q.Cs_iM_CnT);
    _UL.set_sub_mat(0,q.N_CONTACTS,q.Cn_iM_CsT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,q.N_CONTACTS+q.N_STICKING,q.Cs_iM_CsT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,0,q.Ct_iM_CnT);
    _UL.set_sub_mat(0,q.N_CONTACTS+q.N_STICKING*2,q.Cn_iM_CtT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS,q.Ct_iM_CsT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,q.N_CONTACTS+q.N_STICKING,q.Ct_iM_CsT);
    _UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS+q.N_STICKING*2,q.Cs_iM_CtT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,q.N_CONTACTS+q.N_STICKING*3,q.Cs_iM_CtT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS+q.N_STICKING*2,q.Ct_iM_CtT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,q.N_CONTACTS+q.N_STICKING*3,q.Ct_iM_CtT);

  // Joint Limits
  _UL.set_sub_mat(q.N_CONTACTS,N_FRICT,q.Cs_iM_LT);
  _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,N_FRICT,q.Ct_iM_LT);
  _UL.set_sub_mat(N_FRICT,q.N_CONTACTS,q.L_iM_CsT);
  _UL.set_sub_mat(N_FRICT,q.N_CONTACTS+q.N_STICKING*2,q.L_iM_CtT);
    // Set negative submatrices
    /*     n          r          r           r           r
    n                        -Cn_iM_CsT              -Cn_iM_CtT
    r                        -Cs_iM_CsT              -Cs_iM_CtT
    r -Cs_iM_CnT -Cs_iM_CsT              -Cs_iM_CtT
    r                        -Ct_iM_CsT              -Ct_iM_CtT
    r -Ct_iM_CnT -Ct_iM_CsT              -Ct_iM_CtT
      */

  FILE_LOG(LOG_CONSTRAINT) << "Cn*inv(M)*Cs': " << std::endl << q.Cn_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "Cn*inv(M)*Ct': " << std::endl << q.Cn_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "Cs*inv(M)*Cn': " << std::endl << q.Cs_iM_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "Cs*inv(M)*Cs': " << std::endl << q.Cs_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "Cs*inv(M)*Ct': " << std::endl << q.Cs_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "Ct*inv(M)*Cn': " << std::endl << q.Ct_iM_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "Ct*inv(M)*Cs': " << std::endl << q.Ct_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "Cs*inv(M)*L': " << std::endl << q.Cs_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "Ct*inv(M)*L': " << std::endl << q.Ct_iM_LT;
  FILE_LOG(LOG_CONSTRAINT) << "L*inv(M)*Cs': " << std::endl << q.L_iM_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "L*inv(M)*Ct': " << std::endl << q.L_iM_CtT;
  FILE_LOG(LOG_CONSTRAINT) << "L*inv(M)*L': " << std::endl << q.L_iM_LT;

    q.Cn_iM_CsT.negate();
    q.Cn_iM_CtT.negate();
    q.Cs_iM_CnT.negate();
    q.Cs_iM_CsT.negate();
    q.Cs_iM_CtT.negate();
    q.Ct_iM_CnT.negate();
    q.Ct_iM_CsT.negate();
    q.Ct_iM_CtT.negate();
    q.L_iM_CsT.negate();
    q.L_iM_CtT.negate();
    q.Cs_iM_LT.negate();
    q.Ct_iM_LT.negate();

    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,0,q.Cs_iM_CnT);
    _UL.set_sub_mat(0,q.N_CONTACTS+q.N_STICKING,q.Cn_iM_CsT);

    _UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS+q.N_STICKING,q.Cs_iM_CsT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,q.N_CONTACTS,q.Cs_iM_CsT);

    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,0,q.Ct_iM_CnT);
    _UL.set_sub_mat(0,q.N_CONTACTS+q.N_STICKING*3,q.Cn_iM_CtT);

    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,q.N_CONTACTS,q.Ct_iM_CsT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS+q.N_STICKING,q.Ct_iM_CsT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,q.N_CONTACTS+q.N_STICKING*2,q.Cs_iM_CtT);
    _UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS+q.N_STICKING*3,q.Cs_iM_CtT);

    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS+q.N_STICKING*3,q.Ct_iM_CtT);
    _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,q.N_CONTACTS+q.N_STICKING*2,q.Ct_iM_CtT);
  // Joint Limits
  _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,N_FRICT,q.Cs_iM_LT);
  _UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,N_FRICT,q.Ct_iM_LT);
  _UL.set_sub_mat(N_FRICT,q.N_CONTACTS+q.N_STICKING,q.L_iM_CsT);
  _UL.set_sub_mat(N_FRICT,q.N_CONTACTS+q.N_STICKING*3,q.L_iM_CtT);

    // lower left & upper right block of matrix
    for(unsigned i=0,j=0,r=0;i<q.N_CONTACTS;i++)
    {
      const UnilateralConstraint* ci =  q.contact_constraints[i];
      if(ci->get_friction_type() == UnilateralConstraint::eSticking)
      {
        if (ci->contact_NK > 4)
        {
          int nk4 = ( ci->contact_NK+4)/4;
          for(unsigned k=0;k<nk4;k++)
          {
            FILE_LOG(LOG_CONSTRAINT) << "mu_{"<< k<< ","<< i <<"}: " << ci->contact_mu_coulomb << std::endl;
            // muK
            _LL(r+k,i) = ci->contact_mu_coulomb;
            // Xe
            _LL(r+k,q.N_CONTACTS+j)                = -cos((M_PI*k)/(2.0*nk4));
            _LL(r+k,q.N_CONTACTS+q.N_STICKING+j)   = -cos((M_PI*k)/(2.0*nk4));
            // Xf
            _LL(r+k,q.N_CONTACTS+q.N_STICKING*2+j) = -sin((M_PI*k)/(2.0*nk4));
            _LL(r+k,q.N_CONTACTS+q.N_STICKING*3+j) = -sin((M_PI*k)/(2.0*nk4));
            // XeT
            _UR(q.N_CONTACTS+j,r+k)                =  cos((M_PI*k)/(2.0*nk4));
            _UR(q.N_CONTACTS+q.N_STICKING+j,r+k)   =  cos((M_PI*k)/(2.0*nk4));
            // XfT
            _UR(q.N_CONTACTS+q.N_STICKING*2+j,r+k) =  sin((M_PI*k)/(2.0*nk4));
            _UR(q.N_CONTACTS+q.N_STICKING*3+j,r+k) =  sin((M_PI*k)/(2.0*nk4));
          }
          r+=nk4;
          j++;
        }
        else
        {
          FILE_LOG(LOG_CONSTRAINT) << "mu_{"<< i <<"}: " << ci->contact_mu_coulomb << std::endl;
          // muK
          _LL(r,i) = ci->contact_mu_coulomb;
          // Xe
          _LL(r,q.N_CONTACTS+j)                = -1.0;
          _LL(r,q.N_CONTACTS+q.N_STICKING+j)   = -1.0;
          // Xf
          _LL(r,q.N_CONTACTS+q.N_STICKING*2+j) = -1.0;
          _LL(r,q.N_CONTACTS+q.N_STICKING*3+j) = -1.0;
          // XeT
          _UR(q.N_CONTACTS+j,r)                =  1.0;
          _UR(q.N_CONTACTS+q.N_STICKING+j,r)   =  1.0;
          // XfT
          _UR(q.N_CONTACTS+q.N_STICKING*2+j,r) =  1.0;
          _UR(q.N_CONTACTS+q.N_STICKING*3+j,r) =  1.0;
          r += 1;
          j++;
        }
      }
    }

    // setup the LCP matrix
    _MM.set_sub_mat(0, _UL.columns(), _UR);
    _MM.set_sub_mat(_UL.rows(), 0, _LL);

    // setup the LCP vector
    _qq.set_sub_vec(q.N_CONTACTS,q.Cs_a);
    _qq.set_sub_vec(q.N_CONTACTS+q.N_STICKING*2,q.Ct_a);
    q.Cs_a.negate();
    q.Ct_a.negate();
    _qq.set_sub_vec(q.N_CONTACTS+q.N_STICKING,q.Cs_a);
    _qq.set_sub_vec(q.N_CONTACTS+q.N_STICKING*3,q.Ct_a);
  }

  _MM.set_sub_mat(0, 0, _UL);

  FILE_LOG(LOG_CONSTRAINT) << " LCP matrix: " << std::endl << _MM;
  FILE_LOG(LOG_CONSTRAINT) << " LCP vector: " << _qq << std::endl;

  // solve the LCP
  if (!_lcp.lcp_lemke_regularized(_MM, _qq, z, -20, 1, -8))
    return false;

  for(unsigned i=0,j=0;i<q.N_CONTACTS;i++)
  {
    const UnilateralConstraint* ci =  q.contact_constraints[i];
    q.cn[i] = z[i];
    if(ci->get_friction_type() == UnilateralConstraint::eSticking)
    {
      q.cs[i] = z[q.N_CONTACTS+j] - z[q.N_CONTACTS+q.N_STICKING+j];
      q.ct[i] = z[q.N_CONTACTS+q.N_STICKING*2+j] - z[q.N_CONTACTS+q.N_STICKING*3+j];
      j++;
    }
    else
    {
      q.cs[i] = -ci->contact_mu_coulomb*q.cn[i];
      q.ct[i] = 0.0;
    }
  }

  q.l = z.segment(N_FRICT,N_FRICT+N_LIMIT);

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // save normal contact impulses
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    // verify that the constraint type is a contact
    assert(q.contact_constraints[i]->constraint_type == UnilateralConstraint::eContact);

    // setup the contact frame
    P->q.set_identity();
    P->x = q.contact_constraints[i]->contact_point;

    // setup the impulse in the contact frame
    Vector3d f;
    f = q.contact_constraints[i]->contact_normal * q.cn[i];
    f += q.contact_constraints[i]->contact_tan1 * q.cs[i];
    f += q.contact_constraints[i]->contact_tan2 * q.ct[i];

    // setup the spatial force
    SForced fx(boost::const_pointer_cast<const Pose3d>(P));
    fx.set_force(f);

    // transform the impulse to the global frame
    q.contact_constraints[i]->contact_impulse += SMomentumd(Pose3d::transform(GLOBAL, fx));
  }

  // save limit impulses
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    double limit_force = (q.limit_constraints[i]->limit_upper) ? -q.l[i] : q.l[i];
    q.limit_constraints[i]->limit_impulse += limit_force;
  }

  if (LOGGING(LOG_CONSTRAINT))
  {
    // compute LCP 'w' vector
    VectorNd w;
    _MM.mult(z, w) += _qq;

    // output new acceleration
    FILE_LOG(LOG_CONSTRAINT) << " LCP solution (w): " << w << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "new acceleration along contact normal: " << w.segment(0, q.contact_constraints.size()) << std::endl;
  }

  FILE_LOG(LOG_CONSTRAINT) << "cn " << q.cn << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "cs " << q.cs << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ct " << q.ct << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "l " << q.l << std::endl;

  FILE_LOG(LOG_CONSTRAINT) << " LCP solution (z): " << z << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::solve_coulomb_lcp() exited" << std::endl;

  return true;
}

// Solves the sustained constraint LCP w/purely viscous friction
bool SustainedUnilateralConstraintHandler::solve_purely_viscous_lcp(SustainedUnilateralConstraintProblemData& q, VectorNd& z)
{
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::solve_purely_viscous_lcp() entered" << std::endl;

  FILE_LOG(LOG_CONSTRAINT) << " LCP matrix: " << std::endl << _MM;
  FILE_LOG(LOG_CONSTRAINT) << " LCP vector: " << _qq << std::endl;

  const unsigned NCONTACTS = q.N_CONTACTS;
  const unsigned NLIMITS = q.N_LIMITS;
  const unsigned NIMP = q.N_CONSTRAINT_EQNS_IMP;

  // fix the number of variables
  q.N_VARS = NCONTACTS + NLIMITS + NIMP;

  // we do this by solving the MLCP:
  // |  A  C  | | u | + | a | = | 0 |
  // |  D  B  | | v |   | b |   | r |

  // A is the matrix Jx*inv(M)*Jx', Jx is implicit joint constraint Jacobians
  // NOTE: we assume that Jx is of full row rank (no dependent constraints)

  // u = alphax
  // v = [ cn; l ]
  // r = [ Cn*v+; L*v+ ]
  // a = v - inv(M)*S'*muv*S*v - inv(M)*T'*muv*T*v
  // b = 0

  // Assuming that C is of full row rank (no dependent joint constraints)
  // A is invertible; then we just need to solve the LCP:

  // | B - D*inv(A)*C | | v | + | b - D*inv(A)*a | = | w |
  // and use the result to solve for u:
  // u = -inv(A)*(a + Cv)

  // compute SVD of Jx*inv(M)*Jx'
  _A = q.Jx_iM_JxT;
  _LA.svd(_A, _AU, _AS, _AV);

  // setup the B matrix
  // B = [ Cn; L ]*inv(M)*[ Cn' L' ]
  _B.resize(NCONTACTS+NLIMITS, NCONTACTS+NLIMITS);
  _B.set_sub_mat(0, 0, q.Cn_iM_CnT);
  _B.set_sub_mat(0, NCONTACTS, q.Cn_iM_LT);
  _B.set_sub_mat(NCONTACTS, 0, q.Cn_iM_LT, Ravelin::eTranspose);
  _B.set_sub_mat(NCONTACTS, NCONTACTS, q.L_iM_LT);

  // setup the C matrix and compute inv(A)*C
  // C = Jx*inv(M)*[ Cn' L' ]; note: D = C'
  _C.resize(NIMP, NCONTACTS+NLIMITS);
  _C.set_sub_mat(0,0, q.Cn_iM_JxT, Ravelin::eTranspose);
  _C.set_sub_mat(0,NCONTACTS, q.L_iM_JxT, Ravelin::eTranspose);
  MatrixNd::transpose(_C, _D);
  _LA.solve_LS_fast(_AU, _AS, _AV, _C);

  // setup the a vector and compute inv(A)*a
  // a = [ Jx*v ]
  _a = q.Jx_a;
  _LA.solve_LS_fast(_AU, _AS, _AV, _a);

  // setup the b vector
  // b = [ Cn*v; L*v ]
  _b.resize(NLIMITS+NCONTACTS);
  _b.set_sub_vec(0, q.Cn_a);
  _b.set_sub_vec(NCONTACTS, q.L_a);

  // setup the LCP matrix
  _D.mult(_C, _MM);
  _MM -= _B;
  _MM.negate();

  // setup the LCP vector
  _D.mult(_a, _qq);
  _qq -= _b;
  _qq.negate();

  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::solve_lcp() entered" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_iM_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "  Cn * a: " << q.Cn_a << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  L * a: " << q.L_a << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "  LCP matrix: " << std::endl << _MM;
  FILE_LOG(LOG_CONSTRAINT) << "  LCP vector: " << _qq << std::endl;

  q.Cn_iM_CsT.negate();
  q.Cn_iM_CtT.negate();
  q.Cs_iM_CnT.negate();
  q.Cs_iM_CsT.negate();
  q.Cs_iM_CtT.negate();
  q.Ct_iM_CnT.negate();
  q.Ct_iM_CsT.negate();
  q.Ct_iM_CtT.negate();
  q.L_iM_CsT.negate();
  q.L_iM_CtT.negate();
  q.Cs_iM_LT.negate();
  q.Ct_iM_LT.negate();
  q.Cs_a.negate();
  q.Ct_a.negate();

  // solve the LCP
  if (!_lcp.lcp_fast(_MM, _qq, _v) && !_lcp.lcp_lemke_regularized(_MM, _qq, _v))
    throw std::runtime_error("Unable to solve constraint LCP!");

  // compute alphax
  // u = -inv(A)*(a + Cv)
  _C.mult(_v, _alpha_x) += _a;
  _alpha_x.negate();

  // setup the homogeneous solution
  z.set_zero(q.N_VARS);
  z.set_sub_vec(q.CN_IDX, _v);
  z.set_sub_vec(q.ALPHA_X_IDX, _alpha_x);

  FILE_LOG(LOG_CONSTRAINT) << " LCP solution (z): " << z << std::endl;

  // get contact, joint limit, and joint constraint forces
  q.cn = z.segment(q.CN_IDX, q.N_CONTACTS);
  q.l = z.segment(q.L_IDX, q.L_IDX+q.N_LIMITS);
  q.alpha_x = z.segment(q.ALPHA_X_IDX, q.ALPHA_X_IDX + q.N_CONSTRAINT_EQNS_IMP);
  q.cs = _cs_visc;
  q.ct.set_zero(q.cs.size());
  q.cs.negate();
  q.ct.negate();

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // save normal contact impulses
  for (unsigned i=0; i< q.contact_constraints.size(); i++)
  {
    // setup the contact frame
    P->q.set_identity();
    P->x = q.contact_constraints[i]->contact_point;

    // setup the impulse in the contact frame
    Vector3d f;
    f = q.contact_constraints[i]->contact_normal * q.cn[i];
    f += q.contact_constraints[i]->contact_tan1 * q.cs[i];
    f += q.contact_constraints[i]->contact_tan2 * q.ct[i];

    // setup the spatial force
    SForced fx(boost::const_pointer_cast<const Pose3d>(P));
    fx.set_force(f);

    // transform the impulse to the global frame
    q.contact_constraints[i]->contact_impulse += SMomentumd(Pose3d::transform(GLOBAL, fx));
  }

  // save limit impulses
  for (unsigned i=0; i< q.limit_constraints.size(); i++)
  {
    double limit_force = (q.limit_constraints[i]->limit_upper) ? -q.l[i] : q.l[i];
    q.limit_constraints[i]->limit_impulse += limit_force;
  }

  if (LOGGING(LOG_CONSTRAINT))
  {
    // compute LCP 'w' vector
    VectorNd w;
    _MM.mult(z, w) += _qq;

    // output new acceleration
    FILE_LOG(LOG_CONSTRAINT) << "new normal acceleration: " << w.segment(0, q.contact_constraints.size()) << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << " LCP solution (w): " << w << std::endl;
  }

  FILE_LOG(LOG_CONSTRAINT) << "cn " << q.cn << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "cs " << q.cs << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ct " << q.ct << std::endl;

  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::solve_purely_viscous_lcp() exited" << std::endl;

  return true;
}

/// Gets the super body (articulated if any)
DynamicBodyPtr SustainedUnilateralConstraintHandler::get_super_body(SingleBodyPtr sb)
{
  ArticulatedBodyPtr ab = sb->get_articulated_body();
  if (ab)
    return ab;
  else
    return sb;
}


