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
#include <Moby/Contact.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/SingleBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/NumericalException.h>
#include <Moby/RestingContactHandler.h>

#include <Moby/ContactProblemData.h>

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

  /// Sets up the default parameters for the impact event handler
  RestingContactHandler::RestingContactHandler(){}

  // Processes impacts
  void RestingContactHandler::process_contacts(const vector<Contact>& contacts)
  {
    FILE_LOG(LOG_EVENT) << "*************************************************************";
    FILE_LOG(LOG_EVENT) << endl;
    FILE_LOG(LOG_EVENT) << "RestingContactHandler::process_contacts() entered";
    FILE_LOG(LOG_EVENT) << endl;
    FILE_LOG(LOG_EVENT) << "*************************************************************";
    FILE_LOG(LOG_EVENT) << endl;

    // apply the method to all contacts
    if (!contacts.empty())
      apply_model(contacts);
    else
      FILE_LOG(LOG_EVENT) << " (no contacts?!)" << endl;

    FILE_LOG(LOG_EVENT) << "*************************************************************" << endl;
    FILE_LOG(LOG_EVENT) << "RestingContactHandler::process_contacts() exited" << endl;
    FILE_LOG(LOG_EVENT) << "*************************************************************" << endl;
  }

  /// Applies the model to a set of contacts
  /**
   * \param contacts a set of contacts
   */
  void RestingContactHandler::apply_model(const vector<Contact>& contacts)
  {
    list<Contact*> contacting;

    // **********************************************************
    // determine sets of connected contacts
    // **********************************************************
    list<list<Contact*> > groups;
    Contact::determine_connected_contacts(contacts, groups);
    Contact::remove_inactive_groups(groups);

    // **********************************************************
    // do method for each connected set
    // **********************************************************
    for (list<list<Contact*> >::iterator i = groups.begin(); i != groups.end(); i++)
    {
      // determine contact tangents
      for (list<Contact*>::iterator j = i->begin(); j != i->end(); j++)
        if ((*j)->contact_type == Contact::eContact)
          (*j)->determine_contact_tangents();

        // copy the list of contacts
        list<Contact*> rcontacts = *i;

        FILE_LOG(LOG_EVENT) << " -- pre-contact acceleration (all contacts: " << std::endl;
        for (list<Contact*>::iterator j = i->begin(); j != i->end(); j++)
          FILE_LOG(LOG_EVENT) << "    contact: " << std::endl << **j;

        // determine a reduced set of contacts
        Contact::determine_minimal_set(rcontacts);

        // apply model to the reduced contacts
        apply_model_to_connected_contacts(rcontacts);

        FILE_LOG(LOG_EVENT) << " -- post-contact acceleration (all contacts): " << std::endl;
        for (list<Contact*>::iterator j = i->begin(); j != i->end(); j++)
          FILE_LOG(LOG_EVENT) << "    contact: " << std::endl << **j;
    }

//    // determine whether there are any contacting contacts remaining
//    for (list<list<Contact*> >::const_iterator i = groups.begin(); i != groups.end(); i++)
//      for (list<Contact*>::const_iterator j = i->begin(); j != i->end(); j++)
//        if ((*j)->is_contacting())
//            // if there are any contacts still contacting, throw an exception
//            if (!contacting.empty())
//              throw ImpactToleranceException(contacting);
  }

  /**
   * Applies method of Drumwright and Shell to a set of connected contacts
   * \param contacts a set of connected contacts
   */
  void RestingContactHandler::apply_model_to_connected_contacts(const list<Contact*>& contacts)
  {
    double ke_minus = 0.0, ke_plus = 0.0;
    SAFESTATIC ContactProblemData epd;

    FILE_LOG(LOG_EVENT) << "RestingContactHandler::apply_model_to_connected_contacts() entered" << endl;

    // reset problem data
    epd.reset();

    // save the contacts
    epd.contacts = vector<Contact*>(contacts.begin(), contacts.end());

    // compute all contact cross-terms
    compute_problem_data(epd);

    // compute energy
    if (LOGGING(LOG_EVENT))
    {
      for (unsigned i=0; i< epd.super_bodies.size(); i++)
      {
        double ke = epd.super_bodies[i]->calc_kinetic_energy();
        FILE_LOG(LOG_EVENT) << "  body " << epd.super_bodies[i]->id << " pre-contact handling KE: " << ke << endl;
        ke_minus += ke;
      }
    }

    // solve the (non-frictional) linear complementarity problem to determine
    // the kappa constant
    VectorNd z;
    solve_lcp(epd, z);
//    epd.kappa = 10;

    // apply FORCES
    apply_forces(epd);

    // compute energy
    if (LOGGING(LOG_EVENT))
    {
      for (unsigned i=0; i< epd.super_bodies.size(); i++)
      {
        double ke = epd.super_bodies[i]->calc_kinetic_energy();
        FILE_LOG(LOG_EVENT) << "  body " << epd.super_bodies[i]->id << " post-contact handling KE: " << ke << endl;
        ke_plus += ke;
      }
      if (ke_plus > ke_minus){
          FILE_LOG(LOG_EVENT) << "warning! KE gain detected! energy before=" << ke_minus << " energy after=" << ke_plus << endl;
          FILE_LOG(LOG_EVENT) << "Resting Contacting formulation failed: try using Impact model instead" << endl;
      }
    }

    FILE_LOG(LOG_EVENT) << "RestingContactHandler::apply_model_to_connected_contacts() exiting" << endl;
  }

  /// Applies impulses to bodies
  void RestingContactHandler::apply_forces(const ContactProblemData& q) const
  {
    map<DynamicBodyPtr, VectorNd> gj;
    map<DynamicBodyPtr, VectorNd>::iterator gj_iter;
    VectorNd workv;

    // loop over all contact contacts first
    for (unsigned i=0; i< q.contacts.size(); i++)
    {
      // get the contact force
      const Contact& c = *q.contacts[i];
      SForced w(c.contact_force);
      const Point3d& p = c.contact_point;

      // get the two single bodies of the contact
      SingleBodyPtr sb1 = e.contact_geom1->get_single_body();
      SingleBodyPtr sb2 = e.contact_geom2->get_single_body();

      // get the two super bodies
      DynamicBodyPtr b1 = sb1->get_super_body();
      DynamicBodyPtr b2 = sb2->get_super_body();

      // convert force on first body to generalized forces
      if ((gj_iter = gj.find(b1)) == gj.end())
        b1->convert_to_generalized_force(sb1, w, p, gj[b1]);
      else
      {
        b1->convert_to_generalized_force(sb1, w, p, workv);
        gj_iter->second += workv;
      }

      // convert force on second body to generalized forces
      if ((gj_iter = gj.find(b2)) == gj.end())
        b2->convert_to_generalized_force(sb2, -w, p, gj[b2]);
      else
      {
        b2->convert_to_generalized_force(sb2, -w, p, workv);
        gj_iter->second += workv;
      }
    }

     // apply all generalized forces NOTE: this is simply chnaged from apply_generalized_impulse()-- scale by dt?
    for (map<DynamicBodyPtr, VectorNd>::const_iterator i = gj.begin(); i != gj.end(); i++)
      i->first->apply_generalized_force(i->second);
  }

  /// Computes the data to the LCP / QP problems
  void RestingContactHandler::compute_problem_data(ContactProblemData& q)
  {
    const unsigned UINF = std::numeric_limits<unsigned>::max();
    MatrixNd workM;
    VectorNd workv;

    // determine set of "super" bodies from contact events
    q.super_bodies.clear();
    for (unsigned i=0; i< q.contact_events.size(); i++)
    {
      q.super_bodies.push_back(get_super_body(q.contact_events[i]->contact_geom1->get_single_body()));
      q.super_bodies.push_back(get_super_body(q.contact_events[i]->contact_geom2->get_single_body()));
    }

    // determine set of "super" bodies from limit events
    for (unsigned i=0; i< q.limit_events.size(); i++)
    {
      RigidBodyPtr outboard = q.limit_events[i]->limit_joint->get_outboard_link();
      q.super_bodies.push_back(get_super_body(outboard));
    }

    // make super bodies vector unique
    std::sort(q.super_bodies.begin(), q.super_bodies.end());
    q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

    // initialize constants and set easy to set constants
    q.N_CONTACTS = q.contact_events.size();
    q.N_LIMITS = q.limit_events.size();

    // setup contact working set
    q.contact_working_set.clear();
    q.contact_working_set.resize(q.N_CONTACTS, true);

    // setup constants related to articulated bodies
    for (unsigned i=0; i< q.super_bodies.size(); i++)
    {
      ArticulatedBodyPtr abody = dynamic_pointer_cast<ArticulatedBody>(q.super_bodies[i]);
      if (abody) {
        q.N_CONSTRAINT_EQNS_IMP += abody->num_constraint_eqns_implicit();
        if (abody->use_advanced_friction_model)
        {
          q.N_CONSTRAINT_DOF_IMP += abody->num_joint_dof_implicit();
          q.N_CONSTRAINT_DOF_EXP += abody->num_joint_dof_explicit();
        }
      }
    }

    // compute number of friction polygon edges
    for (unsigned i=0; i< q.contact_events.size(); i++)
    {
      if (q.contact_events[i]->contact_NK < UINF)
      {
          q.N_K_TOTAL += q.contact_events[i]->contact_NK/2;
      }
      else if (q.contact_events[i]->contact_NK == UINF)
        break;
    }

    // initialize the problem matrices / vectors
    q.Cn_iM_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cn_iM_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cn_iM_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cs_iM_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cs_iM_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Ct_iM_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cn_a.set_zero(q.N_CONTACTS);
    q.Cs_a.set_zero(q.N_CONTACTS);
    q.Ct_a.set_zero(q.N_CONTACTS);
    q.fn.set_zero(q.N_CONTACTS);
    q.fs.set_zero(q.N_CONTACTS);
    q.ft.set_zero(q.N_CONTACTS);

    // make muK, Xs & Xt
    // NOTE: this will work until more efficient methods are implemented
    // ---------------------------------
    unsigned nk4 = (q.N_K_TOTAL+4)/4;
    q.muK.set_zero(q.N_STICKING*nk4,q.N_CONTACTS);
    q.muS.set_zero(q.N_CONTACTS,q.N_CONTACTS);
    q.Xs.set_zero(q.N_STICKING*nk4,q.N_CONTACTS);
    q.Xt.set_zero(q.N_STICKING*nk4,q.N_CONTACTS);


    for(unsigned i=0;i<q.N_CONTACTS;i++)
    {
        for(unsigned j=0;j<q.N_STICKING;)
        {
            if(!q.contact_sticking(i))
            {
                muS(i,i) = q.contact_event[i]->contact_mu_coulomb;
                continue;
            }
            nk4 = (q.contact_events[i]->contact_NK+4)/4;
            for(unsigned k=0;k<nk4;k++)
            {
                q.muK(j*nk4+k,i) = q.contact_event[i]->contact_mu_coulomb;
                q.Xs(j*nk4+k,i) = cos((M_PI*k)/(2.0*K));
                q.Xt(j*nk4+k,i) = sin((M_PI*k)/(2.0*K));
            }
            j++;
        }
    }
    // ---------------------------------

    // setup indices
    q.CN_IDX = 0;
    q.CS_IDX = q.CN_IDX + q.N_CONTACTS;
    q.CT_IDX = q.CS_IDX + q.N_CONTACTS;
    q.NCS_IDX = q.CT_IDX + q.N_CONTACTS;
    q.NCT_IDX = q.NCS_IDX + q.N_LIN_CONE;
    q.CS_U_IDX = q.NCT_IDX + q.N_LIN_CONE;
    q.CT_U_IDX = q.CS_U_IDX + q.N_TRUE_CONE;
    q.L_IDX = q.CT_U_IDX + q.N_TRUE_CONE;
    q.BETA_T_IDX = q.L_IDX + q.N_LIMITS;
    q.ALPHA_X_IDX = q.BETA_T_IDX + q.N_CONSTRAINT_DOF_EXP;
    q.BETA_X_IDX = q.ALPHA_X_IDX + q.N_CONSTRAINT_EQNS_IMP;
    q.N_VARS = q.BETA_X_IDX + q.N_CONSTRAINT_DOF_IMP;

    // TODO: add event computation and cross computation methods to Joint

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

    // process contact events, setting up matrices
    for (unsigned i=0; i< q.contacts.size(); i++)
    {
      const Contact* ci = q.contacts[i];
      const unsigned ROWS = (ci->get_friction_type() == Contact::eSticking) ? 3 : 1;
      // compute cross event data for contact events
      for (unsigned j=0; j< q.contacts.size(); j++)
      {
        const Contact* cj = q.contacts[j];
        const unsigned COLS = (cj->get_friction_type() == Contact::eSticking) ? 3 : 1;

        // reset workM
        workM.set_zero(ROWS, COLS);

        // check whether i==j (single contact event)
        if (i == j)
        {
          // compute matrix / vector for contact event i
          workv.set_zero(ROWS);
          q.contact_events[i]->compute_event_data(workM, workv);

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
            q.Cs_a[i] = *data++;
            q.Ct_a[i] = *data;
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
          // compute matrix for cross event
          q.contact_events[i]->compute_cross_event_data(*q.contact_events[j], workM);

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
              RowIteratord_const data = workM.column_iterator_begin();
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

    // process limit events, setting up matrices
    for (unsigned i=0; i< q.limit_events.size(); i++)
    {
      // compute matrix / vector for contact event i
      q.contact_events[i]->compute_event_data(workM, workv);

      // setup appropriate entry of limit inertia matrix and limit velocity
      q.L_iM_LT(i,i) = workM(0,0);
      q.L_v[i] = workv[0];

      // NOTE: cross event data has already been computed
    }
  }

  /// Solves the Resting Contact LCP
  void RestingContactHandler::solve_lcp(ContactProblemData& q, VectorNd& z)
  {
    MatrixNd workM;
    VectorNd workv;
    SAFESTATIC MatrixNd UL, LL, MM,UR;
    SAFESTATIC VectorNd qq, Cn_aplus;

    /*
      ADD to ContactProblemData:
      Cq
      muS
      muK
      Xs, Xt
      Cn_a
      Cs_a
      Ct_a
      Cn_iM_muS_CqT = N*iM*(muS*Q)'
        for:
          Cn_iM_CnT - Cn_iM_muS_CqT = N*iM*(N-muS*Q) = N*iM*N' - N*iM*(muS*Q)'
          Cs_iM_CnT - Cs_iM_muS_CqT = S*iM*(N-muS*Q) = S*iM*N' - S*iM*(muS*Q)'
          Ct_iM_CnT - Ct_iM_muS_CqT = T*iM*(N-muS*Q) = T*iM*N' - T*iM*(muS*Q)'
          N_STICKING (or N_SLIDING)
    */

    // setup sizes
    UL.resize(q.N_CONTACTS*3, q.N_CONTACTS*3);
    UR.resize(q.N_CONTACTS*3, q.N_STICKING);
    LL.resize(q.N_STICKING*q.N_K_TOTAL, q.N_CONTACTS*3);

    // now do upper right hand block of LCP matrix
    workM = Cn_iM_CnT;
    workM -= Cn_iM_muS_CqT;
    UL.set_sub_mat(0,0,workM);

    workM = Cs_iM_CnT;
    workM -= Cs_iM_muS_CqT;
    UL.set_sub_mat(q.N_CONTACTS,0,workM);

    workM = Ct_iM_CnT;
    workM -= Ct_iM_muS_CqT;
    UL.set_sub_mat(q.N_CONTACTS*2,0,workM);

    UL.set_sub_mat(0,q.N_CONTACTS,q.Cn_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS,q.Cs_iM_CsT);
    UL.set_sub_mat(q.N_CONTACTS*2,q.N_CONTACTS,q.Ct_iM_CsT);

    UL.set_sub_mat(0,q.N_CONTACTS*2,q.Cn_iM_CtT);
    UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS*2,q.Cs_iM_CtT);
    UL.set_sub_mat(q.N_CONTACTS*2,q.N_CONTACTS*2,q.Ct_iM_CtT);

    // lower left & upper right block of matrix
    UR.set_sub_mat(q.N_CONTACTS,0,q.Xs);
    UR.set_sub_mat(q.N_CONTACTS*2,0,q.Xt);

    LL.set_sub_mat(0,0,q.muK);
    LL.set_sub_mat(0,q.N_CONTACTS,q.Xs);
    LL.set_sub_mat(0,q.N_CONTACTS*2,q.Xt);

    // setup the LCP matrix
    MM.resize(q.N_CONTACTS*3 + q.N_STICKING, q.N_CONTACTS*3 + q.N_STICKING);
    MM.set_sub_mat(0, 0, UL);
    MM.set_sub_mat(0, q.N_CONTACTS*3, UR);
    MM.set_sub_mat(q.N_CONTACTS*3, 0, LL);

    // setup the LCP vector
    qq.set_zero(MM.rows());
    N.mult(aminus,workv);
    qq.set_sub_vec(0,workv);
    S.mult(aminus,workv);
    qq.set_sub_vec(q.N_CONTACTS,workv);
    T.mult(aminus,workv);
    qq.set_sub_vec(q.N_CONTACTS*2,workv);

    FILE_LOG(LOG_EVENT) << "RestingContactHandler::solve_lcp() entered" << std::endl;
    FILE_LOG(LOG_EVENT) << "  Cn * inv(M) * Cn': " << std::endl << q.Cn_iM_CnT;
    FILE_LOG(LOG_EVENT) << "  Cn * a: " << q.Cn_a << std::endl;
    FILE_LOG(LOG_EVENT) << "  LCP matrix: " << std::endl << MM;
    FILE_LOG(LOG_EVENT) << "  LCP vector: " << qq << std::endl;

    // solve the LCP
    if (!_lcp.lcp_lemke_regularized(MM, qq, z))
      throw std::runtime_error("Unable to solve event LCP!");

    z.get_sub_vec(0,q.fn);
    z.get_sub_vec(q.N_CONTACTS,q.fs);
    z.get_sub_vec(q.N_CONTACTS*2,q.ft);

    FILE_LOG(LOG_EVENT) << "  LCP result: " << z << std::endl;
    FILE_LOG(LOG_EVENT) << "  kappa: " << q.kappa << std::endl;
    FILE_LOG(LOG_EVENT) << "  f: " << f << std::endl;
    FILE_LOG(LOG_EVENT) << "RestingContactHandler::solve_lcp() exited" << std::endl;
  }

  /// Gets the super body (articulated if any)
  DynamicBodyPtr RestingContactHandler::get_super_body(SingleBodyPtr sb)
  {
    ArticulatedBodyPtr ab = sb->get_articulated_body();
    if (ab)
      return ab;
    else
      return sb;
  }
