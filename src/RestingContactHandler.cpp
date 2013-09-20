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

#include <Moby/Event.h>
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
  void RestingContactHandler::process_contacts(const vector<Event>& contacts)
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
  void RestingContactHandler::apply_model(const vector<Event>& contacts)
  {
    list<Event*> contacting;

    // **********************************************************
    // determine sets of connected contacts
    // **********************************************************
    list<list<Event*> > groups;
    Event::determine_connected_events(contacts, groups);
    Event::remove_nonimpacting_groups(groups);

    // **********************************************************
    // do method for each connected set
    // **********************************************************
    for (list<list<Event*> >::iterator i = groups.begin(); i != groups.end(); i++)
    {
      // determine contact tangents
      for (list<Event*>::iterator j = i->begin(); j != i->end(); j++)
        if ((*j)->event_type == Event::eContact)
          (*j)->determine_contact_tangents();
        // copy the list of contacts
        list<Event*> rcontacts = *i;

        FILE_LOG(LOG_EVENT) << " -- pre-contact acceleration (all contacts: " << std::endl;
        for (list<Event*>::iterator j = i->begin(); j != i->end(); j++)
          FILE_LOG(LOG_EVENT) << "    contact: " << std::endl << **j;

        // determine a reduced set of contacts
        Event::determine_minimal_set(rcontacts);

        // apply model to the reduced contacts
        apply_model_to_connected_contacts(rcontacts);

        FILE_LOG(LOG_EVENT) << " -- post-contact acceleration (all contacts): " << std::endl;
        for (list<Event*>::iterator j = i->begin(); j != i->end(); j++)
          FILE_LOG(LOG_EVENT) << "    contact: " << std::endl << **j;
    }

//    // determine whether there are any contacting contacts remaining
//    for (list<list<Event*> >::const_iterator i = groups.begin(); i != groups.end(); i++)
//      for (list<Event*>::const_iterator j = i->begin(); j != i->end(); j++)
//        if ((*j)->is_contacting())
//            // if there are any contacts still contacting, throw an exception
//            if (!contacting.empty())
//              throw ImpactToleranceException(contacting);
  }

  /**
   * Applies method of Drumwright and Shell to a set of connected contacts
   * \param contacts a set of connected contacts
   */
  void RestingContactHandler::apply_model_to_connected_contacts(const list<Event*>& contacts)
  {
    double ke_minus = 0.0, ke_plus = 0.0;
    SAFESTATIC ContactProblemData epd;

    FILE_LOG(LOG_EVENT) << "RestingContactHandler::apply_model_to_connected_contacts() entered" << endl;

    // reset problem data
    epd.reset();

    // save the contacts
    epd.events = vector<Event*>(contacts.begin(), contacts.end());

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
    std::cout << "Resting Event forces : " << z << std::endl;
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
    for (unsigned i=0; i<  q.events.size(); i++)
    {
      // get the contact force
      const Event& c = * q.events[i];
      SForced w(c.contact_impulse);
      const Point3d& p = c.contact_point;

      // get the two single bodies of the contact
      SingleBodyPtr sb1 = c.contact_geom1->get_single_body();
      SingleBodyPtr sb2 = c.contact_geom2->get_single_body();

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
      i->first->apply_generalized_impulse(i->second);
  }

  /// Computes the data to the LCP / QP problems
  void RestingContactHandler::compute_problem_data(ContactProblemData& q)
  {
    const unsigned UINF = std::numeric_limits<unsigned>::max();
    SAFESTATIC MatrixNd workM;
    SAFESTATIC VectorNd workv;

    // determine set of "super" bodies from contact events
    q.super_bodies.clear();
    for (unsigned i=0; i< q.events.size(); i++)
    {
      q.super_bodies.push_back(get_super_body(q.events[i]->contact_geom1->get_single_body()));
      q.super_bodies.push_back(get_super_body(q.events[i]->contact_geom2->get_single_body()));
    }

    // make super bodies vector unique
    std::sort(q.super_bodies.begin(), q.super_bodies.end());
    q.super_bodies.erase(std::unique(q.super_bodies.begin(), q.super_bodies.end()), q.super_bodies.end());

    // initialize constants and set easy to set constants
    q.N_CONTACTS =  q.events.size();

    // setup contact working set
    q.contact_working_set.clear();
    q.contact_working_set.resize(q.N_CONTACTS, true);

    // compute number of friction polygon edges
    for (unsigned i=0; i<  q.events.size(); i++)
    {
      if ( q.events[i]->contact_NK < UINF)
      {
          q.N_K_TOTAL +=  q.events[i]->contact_NK/2;
      }
      else if ( q.events[i]->contact_NK == UINF)
        break;
    }

    // initialize the problem matrices / vectors
    q.Cn_iM_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cn_iM_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cn_iM_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cs_iM_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cs_iM_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cs_iM_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Ct_iM_CnT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Ct_iM_CsT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Ct_iM_CtT.set_zero(q.N_CONTACTS, q.N_CONTACTS);
    q.Cn_a.set_zero(q.N_CONTACTS);
    q.Cs_a.set_zero(q.N_CONTACTS);
    q.Ct_a.set_zero(q.N_CONTACTS);
    q.cn.set_zero(q.N_CONTACTS);
    q.cs.set_zero(q.N_CONTACTS);
    q.ct.set_zero(q.N_CONTACTS);


    // setup indices
    q.CN_IDX = 0;
    q.CS_IDX = q.CN_IDX + q.N_CONTACTS;
    q.CT_IDX = q.CS_IDX + q.N_CONTACTS;
    q.NCS_IDX = q.CT_IDX + q.N_CONTACTS;
//    q.NCT_IDX = q.NCS_IDX + q.N_LIN_CONE;
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
//    q.N_STICKING = 0;
    for (unsigned i=0; i<  q.events.size(); i++)
    {
      const Event* ci =  q.events[i];
//      q.N_STICKING += (ci->get_friction_type() == Event::eSticking) ? 1 : 0;
      const unsigned ROWS = (ci->get_friction_type() == Event::eSticking) ? 3 : 1;
      // compute cross event data for contact events
      for (unsigned j=0; j<  q.events.size(); j++)
      {
        const Event* cj =  q.events[j];
        const unsigned COLS = (cj->get_friction_type() == Event::eSticking) ? 3 : 1;

        // reset workM
        workM.set_zero(ROWS, COLS);

        // check whether i==j (single contact event)
        if (i == j)
        {
          // compute matrix / vector for contact event i
          workv.set_zero(ROWS);
           q.events[i]->compute_event_data(workM, workv);

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
           q.events[i]->compute_cross_event_data(* q.events[j], workM);

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
    q.N_STICKING = q.Cn_iM_CsT.columns()/2;

  }

  /// Solves the Resting Event LCP
 void RestingContactHandler::solve_lcp(ContactProblemData& q, VectorNd& z)
 {
  SAFESTATIC MatrixNd UL, LL, MM,UR,workM;
  SAFESTATIC VectorNd qq,workv;
  FILE_LOG(LOG_EVENT) << "RestingContactHandler::solve_lcp() entered" << std::endl;

  const unsigned NK_DIRS = q.N_STICKING*((q.N_K_TOTAL+4)/4);
  // setup sizes
  UL.resize(q.N_CONTACTS+q.N_STICKING*4, q.N_CONTACTS+q.N_STICKING*4);
  UR.resize(q.N_CONTACTS+q.N_STICKING*4, NK_DIRS);
  LL.resize(NK_DIRS, q.N_CONTACTS+q.N_STICKING*4);

  // now do upper right hand block of LCP matrix
  UL.set_sub_mat(0,0,q.Cn_iM_CnT);
  UL.set_sub_mat(q.N_CONTACTS,0,q.Cs_iM_CnT);
  UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,0,q.Ct_iM_CnT);

  UL.set_sub_mat(0,q.N_CONTACTS,q.Cn_iM_CsT);
  UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS,q.Cs_iM_CsT);
  UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS,q.Ct_iM_CsT);

  UL.set_sub_mat(0,q.N_CONTACTS+q.N_STICKING*2,q.Cn_iM_CtT);
  UL.set_sub_mat(q.N_CONTACTS,q.N_CONTACTS+q.N_STICKING*2,q.Cs_iM_CtT);
  UL.set_sub_mat(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS+q.N_STICKING*2,q.Ct_iM_CtT);

  // lower left & upper right block of matrix
  for(unsigned i=0,j=0;i<q.N_CONTACTS;i++)
  {
    const Event* ci =  q.events[i];
    if(ci->get_friction_type() == Event::eSticking)
    {
      int nk4 = ( q.events[i]->contact_NK+4)/4;
      for(unsigned k=0;k<NK_DIRS;k++)
      {
        // TODO: MIGHT NEED TO NEGATE
        // muK
        LL(j*nk4+k,i) =  q.events[i]->contact_mu_coulomb;
        // Xs
        LL(j*nk4+k,j) = -cos((M_PI*k)/(2.0*nk4));
        // XsT
        UR(j,j*nk4+k) = cos((M_PI*k)/(2.0*nk4));
        // Xt
        LL(j*nk4+k,j) = -sin((M_PI*k)/(2.0*nk4));
        // XtT
        UR(j,j*nk4+k) = sin((M_PI*k)/(2.0*nk4));
      }
      j++;
    }
  }

  // setup the LCP matrix
  MM.set_zero(UL.rows() + LL.rows(), UL.columns() + UR.columns());
  MM.set_sub_mat(0, 0, UL);
  MM.set_sub_mat(0, UL.columns(), UR);
  MM.set_sub_mat(UL.rows(), 0, LL);

  // setup the LCP vector
  qq.set_zero(MM.rows());
  qq.set_sub_vec(0,q.Cn_a);
  qq.set_sub_vec(q.Cs_a.rows(),q.Cs_a);
  qq.set_sub_vec(q.Ct_a.rows(),q.Ct_a);

  FILE_LOG(LOG_EVENT) << " LCP matrix: " << std::endl << MM;
  FILE_LOG(LOG_EVENT) << " LCP vector: " << qq << std::endl;

  // solve the LCP
  if (!_lcp.lcp_lemke_regularized(MM, qq, z))
   throw std::runtime_error("Unable to solve resting contact LCP!");

  z.get_sub_vec(0,q.N_CONTACTS,q.cn);
  z.get_sub_vec(q.N_CONTACTS,q.N_CONTACTS+q.N_STICKING,q.cs);
  z.get_sub_vec(q.N_CONTACTS+q.N_STICKING,q.N_CONTACTS+q.N_STICKING*2,workv);
  q.cs -= workv;
  z.get_sub_vec(q.N_CONTACTS+q.N_STICKING*2,q.N_CONTACTS+q.N_STICKING*3,q.ct);
  z.get_sub_vec(q.N_CONTACTS+q.N_STICKING*3,q.N_CONTACTS+q.N_STICKING*4,workv);
  q.ct -= workv;

  FILE_LOG(LOG_EVENT) << " LCP result (siz [N_CONTACTS+N_STICKING*2): " << z << std::endl;
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
