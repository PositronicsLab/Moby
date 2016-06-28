#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <limits>
#include <set>
#include <cmath>
#include <numeric>
#include <Moby/permute.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/Constants.h>
#include <Moby/UnilateralConstraint.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/ConstraintSimulator.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/ImpactConstraintHandler.h>

using namespace Ravelin;
using namespace Moby;
using std::list;
using boost::shared_ptr;
using std::vector;
using std::map;
using std::endl;
using std::cerr;
using std::pair;
using std::min_element;
using boost::dynamic_pointer_cast;

/**
 * Applies Anitescu-Potra model to connected constraints
 * \param constraints a set of connected constraints
 */
void ImpactConstraintHandler::apply_ap_model_to_connected_constraints(const std::list<UnilateralConstraint*>& constraints, const list<shared_ptr<SingleBodyd> >& single_bodies)
{
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_ap_model_to_connected_constraints() entered" << endl;

  // reset problem data
  _epd.reset();

  // set the simulator
  _epd.simulator = _simulator;

  // save the constraints
  _epd.constraints = vector<UnilateralConstraint*>(constraints.begin(), constraints.end());

  // determine sets of contact and limit constraints
  _epd.partition_constraints();

  // compute all constraint cross-terms
  compute_problem_data(_epd, single_bodies);

  // clear all impulses
  for (unsigned i=0; i< _epd.N_CONTACTS; i++)
    _epd.contact_constraints[i]->contact_impulse.set_zero(GLOBAL);
  for (unsigned i=0; i< _epd.N_LIMITS; i++)
    _epd.limit_constraints[i]->limit_impulse = 0.0;

  // solve the A-P model
  apply_ap_model(_epd);

  // determine velocities due to impulse application
  update_constraint_velocities_from_impulses(_epd);

  // get the constraint violation before applying impulses
  double minv = calc_min_constraint_velocity(_epd);

  // apply restitution
  if (apply_restitution(_epd))
  {
    // determine velocities due to impulse application
    update_constraint_velocities_from_impulses(_epd);

    // check to see whether we need to solve another impact problem
    double minv_plus = calc_min_constraint_velocity(_epd);
    FILE_LOG(LOG_CONSTRAINT) << "Applying restitution" << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "  compression v+ minimum: " << minv << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "  restitution v+ minimum: " << minv_plus << std::endl;
    if (minv_plus < 0.0 && minv_plus < minv - NEAR_ZERO)
    {
      // need to solve another impact problem
      apply_ap_model(_epd);
    }
    else
      propagate_impulse_data(_epd);
  }

  // apply impulses
  apply_impulses(_epd);

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_ap_model_to_connected_constraints() exiting" << endl;
}

/// Solves the Anitescu-Potra model LCP
void ImpactConstraintHandler::apply_ap_model(UnilateralConstraintProblemData& q)
{
  /// Matrices and vectors for solving LCP
  Ravelin::MatrixNd _UL, _LL, _MM, _UR, _workM;
  Ravelin::VectorNd _qq, _workv;

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_ap_model() entered" << std::endl;

  unsigned NC = q.N_CONTACTS;

  // Num joint limit constraints
  unsigned N_LIMIT = q.N_LIMITS;

  // Num friction directions + num normal directions
  unsigned N_FRICT = NC*4 + NC;

  // Total constraints
  unsigned N_CONST = N_FRICT + N_LIMIT;

  // Num friction constraints
  unsigned NK_DIRS = 0;
  for(unsigned i=0,j=0,r=0;i<NC;i++){
    if (q.contact_constraints[i]->contact_NK > 4)
      NK_DIRS+=(q.contact_constraints[i]->contact_NK+4)/4;
    else
      NK_DIRS+=1;
  }

  // setup sizes
  _UL.set_zero(N_CONST, N_CONST);
  _UR.set_zero(N_CONST, NK_DIRS);
  _LL.set_zero(NK_DIRS, N_CONST);
  _MM.set_zero(_UL.rows() + _LL.rows(), _UL.columns() + _UR.columns());

  MatrixNd Cn_X_CnT, Cs_X_CnT,Ct_X_CnT,Ct_X_CsT,L_X_CnT,L_X_CsT,L_X_CtT, L_X_LT;
  Cn_X_CnT = q.Cn_X_CnT;
  L_X_LT = q.L_X_LT;
  Ravelin::MatrixNd::transpose(q.Cn_X_LT,L_X_CnT);
  Ravelin::MatrixNd::transpose(q.Cs_X_LT,L_X_CsT);
  Ravelin::MatrixNd::transpose(q.Ct_X_LT,L_X_CtT);
  Ravelin::MatrixNd::transpose(q.Cn_X_CsT,Cs_X_CnT);
  Ravelin::MatrixNd::transpose(q.Cn_X_CtT,Ct_X_CnT);
  Ravelin::MatrixNd::transpose(q.Cs_X_CtT,Ct_X_CsT);
  /*     n          r          r           r           r
  n  Cn_X_CnT  Cn_X_CsT  -Cn_X_CsT   Cn_X_CtT  -Cn_X_CtT
  r  Cs_X_CnT  Cs_X_CsT  -Cs_X_CsT   Cs_X_CtT  -Cs_X_CtT
  r -Cs_X_CnT -Cs_X_CsT   Cs_X_CsT  -Cs_X_CtT   Cs_X_CtT
  r  Ct_X_CnT  Ct_X_CsT  -Ct_X_CsT   Ct_X_CtT  -Ct_X_CtT
  r -Ct_X_CnT -Ct_X_CsT   Ct_X_CsT  -Ct_X_CtT   Ct_X_CtT
  */
  FILE_LOG(LOG_CONSTRAINT) << "Cn*inv(M)*Cn': " << std::endl << q.Cn_X_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "Cn*inv(M)*Cs': " << std::endl << q.Cn_X_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "Cn*inv(M)*Ct': " << std::endl << q.Cn_X_CtT;
  
  FILE_LOG(LOG_CONSTRAINT) << "Cs*inv(M)*Cn': " << std::endl << Cs_X_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "Cs*inv(M)*Cs': " << std::endl << q.Cs_X_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "Cs*inv(M)*Ct': " << std::endl << q.Cs_X_CsT;
  
  FILE_LOG(LOG_CONSTRAINT) << "Ct*inv(M)*Cn': " << std::endl << Ct_X_CnT;
  FILE_LOG(LOG_CONSTRAINT) << "Ct*inv(M)*Cs': " << std::endl << Ct_X_CsT;

  FILE_LOG(LOG_CONSTRAINT) << "L*inv(M)*L': " << std::endl << q.L_X_LT;
  FILE_LOG(LOG_CONSTRAINT) << "Cn*inv(M)*L': " << std::endl << q.Cn_X_LT;
  FILE_LOG(LOG_CONSTRAINT) << "L*inv(M)*Cn': " << std::endl << L_X_CnT;

  FILE_LOG(LOG_CONSTRAINT) << "Cs*inv(M)*L': " << std::endl << q.Cs_X_LT;
  FILE_LOG(LOG_CONSTRAINT) << "Ct*inv(M)*L': " << std::endl << q.Ct_X_LT;
  FILE_LOG(LOG_CONSTRAINT) << "L*inv(M)*Cs': " << std::endl << L_X_CsT;
  FILE_LOG(LOG_CONSTRAINT) << "L*inv(M)*Ct': " << std::endl << L_X_CtT;

  // setup contact damping
  const unsigned N_CONTACTS = q.contact_constraints.size();
  ColumnIteratord Cn_X_CnT_iter = Cn_X_CnT.column_iterator_begin(); 
  for (unsigned i=0; i< N_CONTACTS; i++, Cn_X_CnT_iter++)
  {
    const UnilateralConstraint& c = *q.contact_constraints[i];
    double compliant_layer_depth = c.contact_geom1->compliant_layer_depth +
                                   c.contact_geom2->compliant_layer_depth;
    if (c.signed_distance + compliant_layer_depth > 0)
      *Cn_X_CnT_iter += q.contact_constraints[i]->contact_damping;
  }

  // setup limit damping
  const unsigned N_LIMITS = q.limit_constraints.size();
  ColumnIteratord L_X_LT_iter = L_X_LT.column_iterator_begin(); 
  for (unsigned i=0; i< N_LIMITS; i++, L_X_LT_iter++)
  {
    const UnilateralConstraint& c = *q.limit_constraints[i];
    double compliant_layer_depth = c.limit_joint->compliant_layer_depth;
    if (c.signed_distance + compliant_layer_depth > 0)
      *L_X_LT_iter += q.limit_constraints[i]->limit_damping;
  }

  // Set positive submatrices
  /*
          n          r          r           r           r
  n  Cn_X_CnT  Cn_X_CsT               Cn_X_CtT
  r  Cs_X_CnT  Cs_X_CsT               Cs_X_CtT
  r                         Cs_X_CsT               Cs_X_CtT
  r  Ct_X_CnT  Ct_X_CsT               Ct_X_CtT
  r                         Ct_X_CsT               Ct_X_CtT
  */
  _UL.set_sub_mat(0,0,Cn_X_CnT);

  // setup the LCP matrix

  // setup the LCP vector
  _qq.set_zero(_MM.rows());
  _qq.set_sub_vec(0,q.Cn_v);

  // setup contact stiffness
  for (unsigned i=0; i< N_CONTACTS; i++)
    _qq[i] += q.contact_constraints[i]->signed_distance * q.contact_constraints[i]->contact_stiffness;

  // setup limit stiffness
  for (unsigned i=0; i< N_LIMITS; i++)
    _qq[i+N_FRICT] += q.limit_constraints[i]->signed_distance * q.limit_constraints[i]->limit_stiffness;

  _UL.set_sub_mat(NC,NC,q.Cs_X_CsT);
  _UL.set_sub_mat(NC,0,Cs_X_CnT);
  _UL.set_sub_mat(0,NC,q.Cn_X_CsT);
  _UL.set_sub_mat(NC+NC,NC+NC,q.Cs_X_CsT);
  _UL.set_sub_mat(NC+NC*2,0,Ct_X_CnT);
  _UL.set_sub_mat(0,NC+NC*2,q.Cn_X_CtT);
  _UL.set_sub_mat(NC+NC*2,NC,Ct_X_CsT);
  _UL.set_sub_mat(NC+NC*3,NC+NC,Ct_X_CsT);
  _UL.set_sub_mat(NC,NC+NC*2,q.Cs_X_CtT);
  _UL.set_sub_mat(NC+NC,NC+NC*3,q.Cs_X_CtT);
  _UL.set_sub_mat(NC+NC*2,NC+NC*2,q.Ct_X_CtT);
  _UL.set_sub_mat(NC+NC*3,NC+NC*3,q.Ct_X_CtT);

  // Joint Limits
  _UL.set_sub_mat(N_FRICT,N_FRICT,L_X_LT);
  _UL.set_sub_mat(N_FRICT,0,L_X_CnT);
  _UL.set_sub_mat(0,N_FRICT,q.Cn_X_LT);
  _UL.set_sub_mat(NC,N_FRICT,q.Cs_X_LT);
  _UL.set_sub_mat(NC+NC*2,N_FRICT,q.Ct_X_LT);
  _UL.set_sub_mat(N_FRICT,NC,L_X_CsT);
  _UL.set_sub_mat(N_FRICT,NC+NC*2,L_X_CtT);


  // Set negative submatrices
  /*     n          r          r           r           r
  n                        -Cn_X_CsT              -Cn_X_CtT
  r                        -Cs_X_CsT              -Cs_X_CtT
  r -Cs_X_CnT -Cs_X_CsT              -Cs_X_CtT
  r                        -Ct_X_CsT              -Ct_X_CtT
  r -Ct_X_CnT -Ct_X_CsT              -Ct_X_CtT
    */

  q.Cn_X_CsT.negate();
  q.Cn_X_CtT.negate();
  Cs_X_CnT.negate();
  q.Cs_X_CsT.negate();
  q.Cs_X_CtT.negate();
  Ct_X_CnT.negate();
  Ct_X_CsT.negate();
  q.Ct_X_CtT.negate();

  q.Cs_X_LT.negate();
  q.Ct_X_LT.negate();
  L_X_CsT.negate();
  L_X_CtT.negate();

  _UL.set_sub_mat(NC+NC,0,Cs_X_CnT);
  _UL.set_sub_mat(0,NC+NC,q.Cn_X_CsT);

  _UL.set_sub_mat(NC,NC+NC,q.Cs_X_CsT);
  _UL.set_sub_mat(NC+NC,NC,q.Cs_X_CsT);

  _UL.set_sub_mat(NC+NC*3,0,Ct_X_CnT);
  _UL.set_sub_mat(0,NC+NC*3,q.Cn_X_CtT);

  _UL.set_sub_mat(NC+NC*3,NC,Ct_X_CsT);
  _UL.set_sub_mat(NC+NC*2,NC+NC,Ct_X_CsT);
  _UL.set_sub_mat(NC+NC,NC+NC*2,q.Cs_X_CtT);
  _UL.set_sub_mat(NC,NC+NC*3,q.Cs_X_CtT);

  _UL.set_sub_mat(NC+NC*2,NC+NC*3,q.Ct_X_CtT);
  _UL.set_sub_mat(NC+NC*3,NC+NC*2,q.Ct_X_CtT);

  // Joint limits
  _UL.set_sub_mat(NC+NC,N_FRICT,q.Cs_X_LT);
  _UL.set_sub_mat(NC+NC*3,N_FRICT,q.Ct_X_LT);
  _UL.set_sub_mat(N_FRICT,NC+NC,L_X_CsT);
  _UL.set_sub_mat(N_FRICT,NC+NC*3,L_X_CtT);

  // lower left & upper right block of matrix
  for(unsigned i=0,j=0,r=0;i<NC;i++)
  {
    const UnilateralConstraint* ci =  q.contact_constraints[i];
    if (ci->contact_NK > 4)
    {
      int nk4 = ( ci->contact_NK+4)/4;
      for(unsigned k=0;k<nk4;k++)
      {
        FILE_LOG(LOG_CONSTRAINT) << "mu_{"<< k<< ","<< i <<"}: " << ci->contact_mu_coulomb << std::endl;

        // muK
        _LL(r+k,i)         = ci->contact_mu_coulomb;
        // Xe
        _LL(r+k,NC+j)      = -cos((M_PI*k)/(2.0*nk4));
        _LL(r+k,NC+NC+j)   = -cos((M_PI*k)/(2.0*nk4));
        // Xf
        _LL(r+k,NC+NC*2+j) = -sin((M_PI*k)/(2.0*nk4));
        _LL(r+k,NC+NC*3+j) = -sin((M_PI*k)/(2.0*nk4));
        // XeT
        _UR(NC+j,r+k)      =  cos((M_PI*k)/(2.0*nk4));
        _UR(NC+NC+j,r+k)   =  cos((M_PI*k)/(2.0*nk4));
        // XfT
        _UR(NC+NC*2+j,r+k) =  sin((M_PI*k)/(2.0*nk4));
        _UR(NC+NC*3+j,r+k) =  sin((M_PI*k)/(2.0*nk4));
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
      _LL(r,NC+j)      = -1.0;
      _LL(r,NC+NC+j)   = -1.0;
      // Xf
      _LL(r,NC+NC*2+j) = -1.0;
      _LL(r,NC+NC*3+j) = -1.0;
      // XeT
      _UR(NC+j,r)      =  1.0;
      _UR(NC+NC+j,r)   =  1.0;
      // XfT
      _UR(NC+NC*2+j,r) =  1.0;
      _UR(NC+NC*3+j,r) =  1.0;
      r += 1;
      j++;
    }
  }

  // setup the LCP matrix
  _MM.set_sub_mat(0, _UL.columns(), _UR);
  _MM.set_sub_mat(_UL.rows(), 0, _LL);

  // setup the LCP vector
  _qq.set_sub_vec(NC,q.Cs_v);
  _qq.set_sub_vec(NC+NC*2,q.Ct_v);
  q.Cs_v.negate();
  q.Ct_v.negate();
  _qq.set_sub_vec(NC+NC,q.Cs_v);
  _qq.set_sub_vec(NC+NC*3,q.Ct_v);
  _qq.set_sub_vec(N_FRICT,q.L_v);

  _MM.set_sub_mat(0, 0, _UL);

  FILE_LOG(LOG_CONSTRAINT) << " LCP matrix: " << std::endl << _MM;
  FILE_LOG(LOG_CONSTRAINT) << " LCP vector: " << _qq << std::endl;

  // Fix Negations
  q.Cn_X_CsT.negate();
  q.Cn_X_CtT.negate();
//  Cs_X_CnT.negate();
  q.Cs_X_CsT.negate();
  q.Cs_X_CtT.negate();
//  Ct_X_CnT.negate();
//  Ct_X_CsT.negate();
  q.Ct_X_CtT.negate();
  q.Cs_X_LT.negate();
  q.Ct_X_LT.negate();
  //L_X_CsT.negate();
  //L_X_CtT.negate();
  q.Cs_v.negate();
  q.Ct_v.negate();

  // solve the LCP
  VectorNd z;
  if (!_lcp.lcp_lemke_regularized(_MM, _qq, z, -20, 1, -2))
    throw std::exception();

  for(unsigned i=0,j=0;i<NC;i++)
  {
    q.cn[i] = z[i];
    q.cs[i] = z[NC+j] - z[NC+NC+j];
    q.ct[i] = z[NC+NC*2+j] - z[NC+NC*3+j];
    j++;
  }

  q.l = z.segment(N_FRICT,N_FRICT+N_LIMIT);

  // setup a temporary frame
  shared_ptr<Pose3d> P(new Pose3d);

  // propagate the impulse data
  propagate_impulse_data(q); 

  if (LOGGING(LOG_CONSTRAINT))
  {
    // compute LCP 'w' vector
    VectorNd w;
    _MM.mult(z, w) += _qq;

    // output new acceleration
    FILE_LOG(LOG_CONSTRAINT) << "new normal v: " << w.segment(0, q.contact_constraints.size()) << std::endl;
  }

  FILE_LOG(LOG_CONSTRAINT) << "cn " << q.cn << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "cs " << q.cs << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ct " << q.ct << std::endl;

  FILE_LOG(LOG_CONSTRAINT) << "l " << q.l << std::endl;

  FILE_LOG(LOG_CONSTRAINT) << " LCP result : " << z << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::apply_ap_model() exited" << std::endl;
}

