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

// Processes impacts
void SustainedUnilateralConstraintHandler::process_constraints(const vector<UnilateralConstraint>& contacts)
{
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::process_constraints() entered";
  FILE_LOG(LOG_CONSTRAINT) << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************";
  FILE_LOG(LOG_CONSTRAINT) << endl;

  // verify that every constraint is a contact constraint
  for (unsigned i=0; i< contacts.size(); i++)
    assert (contacts[i].constraint_type == UnilateralConstraint::eContact);

  // apply the method to all contacts
  if (!contacts.empty())
    apply_model(contacts);

  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::process_constraints() exited" << endl;
  FILE_LOG(LOG_CONSTRAINT) << "*************************************************************" << endl;
}

/// Applies the model to a set of contacts
/**
 * \param contacts a set of contacts
 */
void SustainedUnilateralConstraintHandler::apply_model(const vector<UnilateralConstraint>& contacts)
{
  list<UnilateralConstraint*> contacting;

  // **********************************************************
  // determine sets of connected contacts
  // **********************************************************
  list<list<UnilateralConstraint*> > groups;
  UnilateralConstraint::determine_connected_constraints(contacts, groups);
  UnilateralConstraint::remove_inactive_groups(groups);

  // **********************************************************
  // do method for each connected set
  // **********************************************************
  for (list<list<UnilateralConstraint*> >::iterator i = groups.begin(); i != groups.end(); i++)
  {
    // determine contact tangents
    for (list<UnilateralConstraint*>::iterator j = i->begin(); j != i->end(); j++)
      if ((*j)->constraint_type == UnilateralConstraint::eContact)
        (*j)->determine_contact_tangents();
      // copy the list of contacts
      list<UnilateralConstraint*> rcontacts = *i;

      FILE_LOG(LOG_CONSTRAINT) << " -- pre-contact acceleration (all contacts: " << std::endl;
      for (list<UnilateralConstraint*>::iterator j = i->begin(); j != i->end(); j++)
        FILE_LOG(LOG_CONSTRAINT) << "    contact: " << std::endl << **j;

      // determine a reduced set of contacts
      UnilateralConstraint::determine_minimal_set(rcontacts);

      // apply model to the reduced contacts
      apply_model_to_connected_contacts(rcontacts);
  }
}

/**
 * \param contacts a set of connected contacts
 */
void SustainedUnilateralConstraintHandler::apply_model_to_connected_contacts(const list<UnilateralConstraint*>& contacts)
{
  SAFESTATIC SustainedUnilateralConstraintProblemData epd;
  SAFESTATIC VectorNd v,a, ke_minus, ke_plus;
  SAFESTATIC vector<VectorNd> gf;
  SAFESTATIC MatrixNd M;

  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::apply_model_to_connected_contacts() entered" << endl;

  // reset problem data
  epd.reset();

  // save the contacts
  epd.constraints = vector<UnilateralConstraint*>(contacts.begin(), contacts.end());

  // compute all contact cross-terms
  compute_problem_data(epd);

  // solve the (non-frictional) linear complementarity problem to determine
  // the kappa constant
  VectorNd z;
  if (!solve_lcp(epd, z))
    throw SustainedUnilateralConstraintSolveFailException();

  FILE_LOG(LOG_CONSTRAINT) << "Resting constraint forces : " << z << std::endl;

  // apply FORCES
  apply_forces(epd);

  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::apply_model_to_connected_contacts() exiting" << endl;
}

/// Applies resting contact forces to bodies and saves the generalized forces
void SustainedUnilateralConstraintHandler::apply_forces(const SustainedUnilateralConstraintProblemData& q) const
{
  map<DynamicBodyPtr, VectorNd> gj;
  map<DynamicBodyPtr, VectorNd>::iterator gj_iter;
  VectorNd workv;

  // loop over all contact contacts first
  for (unsigned i=0; i<  q.constraints.size(); i++)
  {
    // get the contact force
    const UnilateralConstraint& c = * q.constraints[i];
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
      b1->convert_to_generalized_force(sb1, w, workv);
      gj_iter->second += workv;
    }

    // convert force on second body to generalized forces
    if ((gj_iter = gj.find(b2)) == gj.end())
      b2->convert_to_generalized_force(sb2, -w, gj[b2]);
    else
    {
      b2->convert_to_generalized_force(sb2, -w, workv);
      gj_iter->second += workv;
    }
  }

  // TODO: this needs to be done in a different way so that it can override
  // any actuator force limits for robots
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

  // determine set of "super" bodies from contact constraints
  q.super_bodies.clear();
  for (unsigned i=0; i< q.constraints.size(); i++)
  {
    q.super_bodies.push_back(get_super_body(q.constraints[i]->contact_geom1->get_single_body()));
    q.super_bodies.push_back(get_super_body(q.constraints[i]->contact_geom2->get_single_body()));
  }

  // make super bodies vector unique
  std::sort(q.super_bodies.begin(), q.super_bodies.end());
  q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

  // initialize constants and set easy to set constants
  q.N_CONTACTS =  q.constraints.size();

  // setup contact working set
  q.contact_working_set.clear();
  q.contact_working_set.resize(q.N_CONTACTS, true);

  // compute number of friction polygon edges
  q.N_STICKING = 0;
  for (unsigned i=0; i<  q.constraints.size(); i++)
  {
    q.N_STICKING += (q.constraints[i]->get_friction_type() == UnilateralConstraint::eSticking) ? 1 : 0;
    if ( q.constraints[i]->contact_NK < UINF)
    {
        q.N_K_TOTAL +=  q.constraints[i]->contact_NK/2;
    }
    else if ( q.constraints[i]->contact_NK == UINF)
      break;
  }

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
  q.CT_IDX = q.CS_IDX + q.N_CONTACTS;
  q.NCS_IDX = q.CT_IDX + q.N_CONTACTS;
  q.NCT_IDX = q.NCS_IDX + q.N_CONTACTS;
  // TODO: add constraint computation and cross computation methods to Joint

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
  for (unsigned i=0, k=0; i<  q.constraints.size(); i++)
  {
    const UnilateralConstraint* ci =  q.constraints[i];
    const unsigned ROWS = (ci->get_friction_type() == UnilateralConstraint::eSticking) ? 3 : 1;

    // compute cross constraint data for contact constraints
    for (unsigned j=0; j<  q.constraints.size(); j++)
    {
      const UnilateralConstraint* cj =  q.constraints[j];
      const unsigned COLS = (cj->get_friction_type() == UnilateralConstraint::eSticking) ? 3 : 1;

      // reset workM
      workM.set_zero(ROWS, COLS);

      // check whether i==j (single contact constraint)
      if (i == j)
      {
        // compute matrix / vector for contact constraint i
        workv.set_zero(ROWS);
         q.constraints[i]->compute_constraint_data(workM, workv);

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
         q.constraints[i]->compute_cross_constraint_data(* q.constraints[j], workM);

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

 /// Solves the Resting constraint LCP
bool SustainedUnilateralConstraintHandler::solve_lcp(SustainedUnilateralConstraintProblemData& q, VectorNd& z)
{
  SAFESTATIC MatrixNd UL, LL, MM,UR,workM;
  SAFESTATIC VectorNd qq,workv;
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::solve_lcp() entered" << std::endl;

  unsigned NK_DIRS = 0;
  for(unsigned i=0,j=0,r=0;i<q.N_CONTACTS;i++)
    if(q.constraints[i]->get_friction_type() == UnilateralConstraint::eSticking)
      NK_DIRS+=(q.constraints[i]->contact_NK+4)/4;

  // setup sizes
  UL.set_zero(q.N_CONTACTS+q.N_STICKING*4, q.N_CONTACTS+q.N_STICKING*4);
  UR.set_zero(q.N_CONTACTS+q.N_STICKING*4, NK_DIRS);
  LL.set_zero(NK_DIRS, q.N_CONTACTS+q.N_STICKING*4);
  MM.set_zero(UL.rows() + LL.rows(), UL.columns() + UR.columns());

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
  UL.set_sub_mat(0,0,q.Cn_iM_CnT);
  // setup the LCP matrix

  // setup the LCP vector
  qq.set_zero(MM.rows());
  qq.set_sub_vec(0,q.Cn_a);

  if(q.N_STICKING > 0){

    UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS,q.Cs_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS,0,q.Cs_iM_CnT);
    UL.set_sub_mat(0,q.N_CONTACTS,q.Cn_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,q.N_CONTACTS+q.N_STICKING,q.Cs_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,0,q.Ct_iM_CnT);
    UL.set_sub_mat(0,q.N_CONTACTS+q.N_STICKING*2,q.Cn_iM_CtT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS,q.Ct_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,q.N_CONTACTS+q.N_STICKING,q.Ct_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS+q.N_STICKING*2,q.Cs_iM_CtT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,q.N_CONTACTS+q.N_STICKING*3,q.Cs_iM_CtT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS+q.N_STICKING*2,q.Ct_iM_CtT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,q.N_CONTACTS+q.N_STICKING*3,q.Ct_iM_CtT);

    // Set neagtive submatrices
    /*     n          r          r           r           r
    n                        -Cn_iM_CsT              -Cn_iM_CtT
    r                        -Cs_iM_CsT              -Cs_iM_CtT
    r -Cs_iM_CnT -Cs_iM_CsT              -Cs_iM_CtT
    r                        -Ct_iM_CsT              -Ct_iM_CtT
    r -Ct_iM_CnT -Ct_iM_CsT              -Ct_iM_CtT
      */

    q.Cn_iM_CsT.negate();
    q.Cn_iM_CtT.negate();
    q.Cs_iM_CnT.negate();
    q.Cs_iM_CsT.negate();
    q.Cs_iM_CtT.negate();
    q.Ct_iM_CnT.negate();
    q.Ct_iM_CsT.negate();
    q.Ct_iM_CtT.negate();

    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,0,q.Cs_iM_CnT);
    UL.set_sub_mat(0,q.N_CONTACTS+q.N_STICKING,q.Cn_iM_CsT);

    UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS+q.N_STICKING,q.Cs_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,q.N_CONTACTS,q.Cs_iM_CsT);

    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,0,q.Ct_iM_CnT);
    UL.set_sub_mat(0,q.N_CONTACTS+q.N_STICKING*3,q.Cn_iM_CtT);

    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,q.N_CONTACTS,q.Ct_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS+q.N_STICKING,q.Ct_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING,q.N_CONTACTS+q.N_STICKING*2,q.Cs_iM_CtT);
    UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS+q.N_STICKING*3,q.Cs_iM_CtT);

    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS+q.N_STICKING*3,q.Ct_iM_CtT);
    UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*3,q.N_CONTACTS+q.N_STICKING*2,q.Ct_iM_CtT);

    // lower left & upper right block of matrix
    for(unsigned i=0,j=0,r=0;i<q.N_CONTACTS;i++)
    {
      const UnilateralConstraint* ci =  q.constraints[i];
      if(ci->get_friction_type() == UnilateralConstraint::eSticking)
      {
        int nk4 = ( ci->contact_NK+4)/4;
        for(unsigned k=0;k<nk4;k++)
        {
          // muK
          LL(r+k,i) = ci->contact_mu_coulomb;
          // Xs
          LL(r+k,q.N_CONTACTS+j)                = -cos((M_PI*k)/(2.0*nk4));
          LL(r+k,q.N_CONTACTS+q.N_STICKING+j)   = -cos((M_PI*k)/(2.0*nk4));
          // Xt
          LL(r+k,q.N_CONTACTS+q.N_STICKING*2+j) = -sin((M_PI*k)/(2.0*nk4));
          LL(r+k,q.N_CONTACTS+q.N_STICKING*3+j) = -sin((M_PI*k)/(2.0*nk4));
          // XsT
          UR(q.N_CONTACTS+j,r+k)                =  cos((M_PI*k)/(2.0*nk4));
          UR(q.N_CONTACTS+q.N_STICKING+j,r+k)   =  cos((M_PI*k)/(2.0*nk4));
          // XtT
          UR(q.N_CONTACTS+q.N_STICKING*2+j,r+k) =  sin((M_PI*k)/(2.0*nk4));
          UR(q.N_CONTACTS+q.N_STICKING*3+j,r+k) =  sin((M_PI*k)/(2.0*nk4));
        }
        r+=nk4;
        j++;
      }
    }

    // setup the LCP matrix
    MM.set_sub_mat(0, UL.columns(), UR);
    MM.set_sub_mat(UL.rows(), 0, LL);

    // setup the LCP vector
    qq.set_sub_vec(q.N_CONTACTS,q.Cs_a);
    qq.set_sub_vec(q.N_CONTACTS+q.N_STICKING*2,q.Ct_a);
    q.Cs_a.negate();
    q.Ct_a.negate();
    qq.set_sub_vec(q.N_CONTACTS+q.N_STICKING,q.Cs_a);
    qq.set_sub_vec(q.N_CONTACTS+q.N_STICKING*3,q.Ct_a);
  }

  MM.set_sub_mat(0, 0, UL);

  FILE_LOG(LOG_CONSTRAINT) << " LCP matrix: " << std::endl << MM;
  FILE_LOG(LOG_CONSTRAINT) << " LCP vector: " << qq << std::endl;

  // solve the LCP
  //if (!_lcp.lcp_lemke_regularized(MM, qq, z))
  if (!solve_lcp(MM, qq, z))
    return false; 

  for(unsigned i=0,j=0;i<q.N_CONTACTS;i++)
  {
    const UnilateralConstraint* ci =  q.constraints[i];
    q.cn[i] = z[i];
    if(ci->get_friction_type() == UnilateralConstraint::eSticking)
    {
      q.cs[i] = z[q.N_CONTACTS+j] - z[q.N_CONTACTS+q.N_STICKING+j];
      q.ct[i] = z[q.N_CONTACTS+q.N_STICKING*2+j] - z[q.N_CONTACTS+q.N_STICKING*3+j];
      j++;
    }
    else
    {
      // Negated (against dir of sliding Q?)
      q.cs[i] = -ci->contact_mu_coulomb*q.cn[i];
      q.ct[i] = 0.0;
    }
  }

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // save normal contact impulses
  for (unsigned i=0; i< q.constraints.size(); i++)
  {
    // verify that the constraint type is a contact
    assert(q.constraints[i]->constraint_type == UnilateralConstraint::eContact);

    // setup the contact frame
    P->q.set_identity();
    P->x = q.constraints[i]->contact_point;

    // setup the impulse in the contact frame
    Vector3d f;
    f = q.constraints[i]->contact_normal * q.cn[i];
    f += q.constraints[i]->contact_tan1 * q.cs[i];
    f += q.constraints[i]->contact_tan2 * q.ct[i];

    // setup the spatial force
    SForced fx(boost::const_pointer_cast<const Pose3d>(P));
    fx.set_force(f);    

    // transform the impulse to the global frame
    q.constraints[i]->contact_impulse = Pose3d::transform(GLOBAL, fx);
  }

  if (LOGGING(LOG_CONSTRAINT))
  {
    // compute LCP 'w' vector
    VectorNd w;
    MM.mult(z, w) += qq;

    // output new acceleration
    FILE_LOG(LOG_CONSTRAINT) << "new normal acceleration: " << w.segment(0, q.constraints.size()) << std::endl;
  } 

  FILE_LOG(LOG_CONSTRAINT) << "cn " << q.cn << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "cs " << q.cs << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ct " << q.ct << std::endl;

  FILE_LOG(LOG_CONSTRAINT) << " LCP result : " << z << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "SustainedUnilateralConstraintHandler::solve_lcp() exited" << std::endl;

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

bool SustainedUnilateralConstraintHandler::solve_lcp(const MatrixNd& M, const VectorNd& q, VectorNd& z)
{
  if (!_lcp.lcp_lemke(M, q, z))
    return false; 
  else
    return true;
}

