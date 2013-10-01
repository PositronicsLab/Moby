/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iostream>
#include <boost/foreach.hpp>
#include <queue>
#include <Moby/XMLTree.h>
#include <Moby/Joint.h>
#include <Moby/Constants.h>
#include <Moby/RigidBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/NumericalException.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/URDFReader.h>

using namespace Moby;
using namespace Ravelin;
using std::set;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using std::list;
using std::vector;
using std::map;
using std::string;
using std::queue;

ArticulatedBody::ArticulatedBody()
{
  use_advanced_friction_model = false;
}

/// Gets the time-derivative of the Jacobian
/**
 * Columns correspond to joint coordinate indices.
 */
MatrixNd& ArticulatedBody::calc_jacobian_dot(boost::shared_ptr<const Pose3d> frame, DynamicBodyPtr body, MatrixNd& J)
{
  const unsigned SPATIAL_DIM = 6;

  // get the number of explicit degrees of freedom
  const unsigned NIMP_DOF = num_joint_dof_explicit();

  // get the total number of degrees of freedom
  const unsigned NDOF = (is_floating_base()) ? NIMP_DOF + SPATIAL_DIM : NIMP_DOF;

  // setup the Jacobian
  J.set_zero(NDOF, SPATIAL_DIM); 

  // get the current link
  RigidBodyPtr link = dynamic_pointer_cast<RigidBody>(body);
  
  // get the base link
  RigidBodyPtr base = get_base_link();

  // loop backward through (at most one) joint for each child until we reach 
  // the parent
  while (link != base)
  {
    // get the explicit inner joint for this link
    JointPtr joint = link->get_inner_joint_explicit();

    // get the parent link
    RigidBodyPtr parent = joint->get_inboard_link(); 

    // get the coordinate index
    const unsigned CIDX = joint->get_coord_index();

    // get the spatial axes
    const vector<SVelocityd>& s = joint->get_spatial_axes_dot();

    // update J
    for (unsigned i=0; i< s.size(); i++)
    {
      SharedVectorNd v = J.row(CIDX+i);
      Pose3d::transform(frame, s[i]).transpose_to_vector(v);
    }

    // set the link to the parent link
    link = parent;
  }

  return J;
}

/// Gets the Jacobian
/**
 * Columns correspond to joint coordinate indices.
 */
MatrixNd& ArticulatedBody::calc_jacobian(boost::shared_ptr<const Pose3d> frame, DynamicBodyPtr body, MatrixNd& J)
{
  const unsigned SPATIAL_DIM = 6;

  // get the number of explicit degrees of freedom
  const unsigned NIMP_DOF = num_joint_dof_explicit();

  // get the total number of degrees of freedom
  const unsigned NDOF = (is_floating_base()) ? NIMP_DOF + SPATIAL_DIM : NIMP_DOF;

  // setup the Jacobian
  J.set_zero(NDOF,SPATIAL_DIM); 

  // get the current link
  RigidBodyPtr link = dynamic_pointer_cast<RigidBody>(body);
  
  // get the base link
  RigidBodyPtr base = get_base_link();

  // loop backward through (at most one) joint for each child until we reach 
  // the parent
  while (link != base)
  {
    // get the explicit inner joint for this link
    JointPtr joint = link->get_inner_joint_explicit();

    // get the parent link
    RigidBodyPtr parent = joint->get_inboard_link(); 

    // get the coordinate index
    const unsigned CIDX = joint->get_coord_index();

    // get the spatial axes
    const vector<SVelocityd>& s = joint->get_spatial_axes();

    // update J
    for (unsigned i=0; i< s.size(); i++)
    {
      SharedVectorNd v = J.row(CIDX+i);
      Pose3d::transform(frame, s[i]).transpose_to_vector(v);
    }

    // set the link to the parent link
    link = parent;
  }

  // if base is floating, setup Jacobian columns at the end
  if (is_floating_base())
  {
    shared_ptr<const Pose3d> bpose = base->get_mixed_pose();
    SharedMatrixNd Jbase = J.block(NIMP_DOF, NIMP_DOF+SPATIAL_DIM, 0, SPATIAL_DIM);
    Pose3d::spatial_transform_to_matrix(bpose, frame, Jbase);
  }

  return J;
}

/// Gets the maximum angular speed of the links of this articulated body
double ArticulatedBody::get_aspeed()
{
  double max_aspeed = (double) 0.0;
  for (unsigned i=0; i< _links.size(); i++)
  {
    double aspeed = _links[i]->get_aspeed();
    if (aspeed > max_aspeed)
      max_aspeed = aspeed;
  }

  return max_aspeed;
}

/// Computes the Z matrices
/*
void ArticulatedBody::compute_Z_matrices(const vector<unsigned>& loop_indices, const vector<vector<unsigned> >& loop_links, vector<MatrixNd>& Zd, vector<MatrixNd>& Z1d, vector<MatrixNd>& Z) const
{
  SAFESTATIC std::queue<RigidBodyPtr> q;
  SAFESTATIC MatrixNd Fsum, F;
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  const unsigned SPATIAL_DIM = 6;

  // determine number of joint dof to use
  shared_ptr<const RCArticulatedBody> rcab = dynamic_pointer_cast<const RCArticulatedBody>(get_this());
  const unsigned N_JOINT_DOF = (rcab) ? num_joint_dof_explicit() : num_joint_dof();

  // resize vectors
  Zd.resize(_joints.size());
  Z1d.resize(_joints.size());
  Z.resize(_joints.size());

  // loop over all joint constraints
  for (unsigned i=0; i< _joints.size(); i++)
  {
    // get the spatial axes complement for this joint
    const vector<SVelocityd>& si_bar = _joints[i]->get_spatial_axes_complement();

    // determine the target transform
    Point3d xi = _joints[i]->get_position_global();
    Pose3d Tf(&IDENTITY_3x3, &xi);

    // clear Fsum for this joint
    Fsum.set_zero(SPATIAL_DIM, N_JOINT_DOF);

    // we'll need the outboard link for this joint later
    RigidBodyPtr outboard = _joints[i]->get_outboard_link();

    // if this joint is not part of a loop, sum force transforms on all 
    // outboard links
    if (loop_indices[i] == UINF)
    {
      // process all outboard links
      q.push(outboard);
      while (!q.empty())
      {
        // get the link off of the front of the queue
        RigidBodyPtr link = q.front();
        q.pop();

        // add to sum
        Fsum += determine_F(link->get_index(), Tf, loop_indices, F);

        // add all child links to the queue; don't process any twice (so skip
        // any links with indices lower than the one being processed)
        const unsigned N_CHILD_LINKS = link->num_child_links();
        for (unsigned j=0; j< N_CHILD_LINKS; j++)
        {
          RigidBodyPtr child = link->get_child_link(j);
          if (child->get_index() > link->get_index())
            q.push(child);
        } 
      }

      // now, compute Z for this joint
      si_bar.transpose_mult(Fsum, Z[_joints[i]->get_index()]);
    }
    else
    {
      // get the index of the loop that this joint belongs in
      unsigned loop = loop_indices[i];

      // get the links in this loop
      const vector<unsigned>& llinks = loop_links[loop];

      // get the index of the outboard link for this joint in the loop
      vector<unsigned>::const_iterator link_loop_iter = std::find(llinks.begin(), llinks.end(), outboard->get_index());
      assert(link_loop_iter != llinks.end());
      unsigned link_loop_idx = link_loop_iter - llinks.begin();

      // **************************************************************
      // determine sum of F moving forward in the loop
      // **************************************************************
      for (unsigned j=link_loop_idx; j< llinks.size(); j++)
        Fsum += determine_F(llinks[j], Tf, loop_indices, F);

      // now compute forward Z for this joint
      si_bar.transpose_mult(Fsum, Zd[i]);
 
      // **************************************************************
      // clear Fsum and determine sum of F moving backward in the loop
      // **************************************************************
      Fsum.set_zero();
      for (unsigned j=0; j<link_loop_idx; j++)
        Fsum += determine_F(llinks[j], Tf, loop_indices, F);

      // now compute backward Z for this joint
      si_bar.transpose_mult(Fsum, Z1d[i]);
 
      // **************************************************************
      // clear Fsum and determine sum of F after the loop
      // **************************************************************
      Fsum.set_zero();
      for (unsigned j=0; j< llinks.size(); j++)
      {
        // get index for link j
        unsigned lj = llinks[j];

        // get all children of link lj
        const unsigned NUM_CHILDREN = _links[lj]->num_child_links();
        for (unsigned k=0; k< NUM_CHILDREN; k++)
        {
          RigidBodyPtr child = _links[lj]->get_child_link(k);
          if (!std::binary_search(llinks.begin(), llinks.end(), child->get_index()))
            q.push(child);  // link is not in loop; add to queue for processing
        }
      }
      while (!q.empty())
      {
        // get the link off of the front of the queue
        RigidBodyPtr link = q.front();
        q.pop();

        // don't process base link
        if (link->is_base())
          continue;

        // add to sum
        Fsum += determine_F(link->get_index(), Tf, loop_indices, F);

        // add all child links to the queue; don't process any twice (so skip
        // any links with indices lower than the one being processed)
        const unsigned N_CHILD_LINKS = link->num_child_links();
        for (unsigned j=0; j< N_CHILD_LINKS; j++)
        {
          RigidBodyPtr child = link->get_child_link(j);
          if (child->get_index() > link->get_index())
            q.push(child);
        } 
      }

      // now compute Z for this joint
      si_bar.transpose_mult(Fsum, Z[i]);
    }
  }
}
*/

/// Computes the 'F' matrix (which transforms joint torques to a spatial force on the given link) in the given frame
/*
MatrixNd& ArticulatedBody::determine_F(unsigned link_idx, const Pose3d& Tf, const vector<unsigned>& loop_indices, MatrixNd& F) const
{
  const unsigned SPATIAL_DIM = 6;
  const unsigned UINF = std::numeric_limits<unsigned>::max();
  SAFESTATIC std::queue<JointPtr> q;
  SAFESTATIC MatrixNd J;
  SAFESTATIC vector<SVelocityd> sx;
  SAFESTATIC vector<bool> processed;

  // determine whether the coordinate index must be modified
  shared_ptr<const RCArticulatedBody> rcab = dynamic_pointer_cast<const RCArticulatedBody>(get_this());
  const unsigned SUB = (rcab && rcab->is_floating_base()) ? 6 : 0; 

  // initialize F to all zeros
  if (rcab)
    F.set_zero(SPATIAL_DIM, num_joint_dof_explicit());
  else
    F.set_zero(SPATIAL_DIM, num_joint_dof());

  // if there is no inner explicit link, return zero
  JointPtr inner_explicit = _links[link_idx]->get_inner_joint_explicit();
  if (!inner_explicit)
    return F;

  // do the target transformation (multiplication ordering indicates that 
  // we want spatial axes transformed first to target link frame then to target
  // frame) 
  Pose3d TfX = Tf * _links[link_idx]->get_transform();

  // determine whether the given link is part of a loop
  if (loop_indices[inner_explicit->get_index()] == UINF)
  {  
    // not part of a loop: process all joints going backward
    // add inner joint to the queue
    q.push(inner_explicit);
  }
  else
  {
    // part of a loop: process all joints in the loop *and* all joints going
    // backward
    unsigned loop_idx = loop_indices[inner_explicit->get_index()];
    for (unsigned i=0; i< loop_indices.size(); i++)
      if (loop_indices[i] == loop_idx)
        q.push(_joints[i]);
  }

  // don't process any joints twice
  processed.resize(_joints.size());
  for (unsigned i=0; i< _joints.size(); i++)
    processed[i] = false;

  // process until queue is empty
  while (!q.empty())
  {
    // get the joint off of the front of the queue
    JointPtr joint = q.front();
    q.pop();

    // if the joint has already been processed, don't process it twice
    if (processed[joint->get_index()])
      continue;
    processed[joint->get_index()] = true;

    // get the outboard link and its transform
    RigidBodyPtr outboard = joint->get_outboard_link();
    const Pose3d& To = outboard->get_transform();
    
    // get its spatial axes and transform
    const vector<SVelocityd>& si = joint->get_spatial_axes(eLink);
    SpatialTransform(To, TfX).transform(si, sx);

    // set the appropriate column in F
    const unsigned ST_IDX = joint->get_coord_index() - SUB;
    F.set_sub_mat(0, ST_IDX, sx);
//    for (unsigned i=0; i< THREE_D; i++)
//      for (unsigned j=0; j< sx.columns(); j++)
//      {
//        F(i+THREE_D,ST_IDX+j) = sx(i,j);
//        F(i,ST_IDX+j) = sx(THREE_D+i,j);
//      }
//

    // add all inner joints to the queue *unless* it's an rcab and joint is
    // implicit
    RigidBodyPtr inboard = joint->get_inboard_link();
    const list<RigidBody::InnerJointData>& ijd_list = inboard->get_inner_joints_data();
    BOOST_FOREACH(const RigidBody::InnerJointData& ijd, ijd_list)
    {
      JointPtr ij(ijd.inner_joint);
      if (!rcab || ij->get_constraint_type() == Joint::eExplicit)
        q.push(ij);
    }
  }

  // compute transpose of F
  MatrixNd::transpose(F, J);

  // compute pseudo-inverse of J
  F.copy_from(J);
  try
  {
    LinAlg::pseudo_inverse(F, LinAlg::svd1);
  }
  catch (NumericalException e)
  {
    F.copy_from(J);
    LinAlg::pseudo_inverse(F, LinAlg::svd2);
  }

  return F;
}

// objective and inequality constraint functions for convex optimization
double ArticulatedBody::calc_fwd_dyn_f0(const VectorNd& x, void* data)
{
  SAFESTATIC VectorNd Gx;
  const double INFEAS_TOL = 1e-8;

  // setup constants
  const ABFwdDynOptData& opt_data = *(const ABFwdDynOptData*) data;

  // get necessary data
  const MatrixNd& G = opt_data.G;
  const VectorNd& c = opt_data.c;

  // objective function is quadratic
  G.mult(x, Gx) *= (double) 0.5;
  Gx += c;
  return x.dot(Gx);
}

// inequality constraint functions for convex optimization
void ArticulatedBody::calc_fwd_dyn_fx(const VectorNd& x, VectorNd& fc, void* data)
{
  SAFESTATIC VectorNd Gx, w, tmp, ff, DTbx, lambda;
  const double INFEAS_TOL = 1e-8;

  // optimization vector:
  // ff
  // alphax
  // betax

  // setup constants
  const ABFwdDynOptData& opt_data = *(const ABFwdDynOptData*) data;
  const unsigned N_IMPLICIT_DOF = opt_data.N_IMPLICIT_DOF;
  const unsigned N_EXPLICIT_CONSTRAINT_EQNS = opt_data.N_EXPLICIT_CONSTRAINT_EQNS;
  const unsigned N_JOINT_DOF = opt_data.N_JOINT_DOF;
  const unsigned N_LOOPS = opt_data.N_LOOPS;
  const unsigned FF_START = 0;
  const unsigned DELTA_START = N_JOINT_DOF + N_EXPLICIT_CONSTRAINT_EQNS;
  const unsigned BETA_START = N_IMPLICIT_DOF + N_EXPLICIT_CONSTRAINT_EQNS;
  const unsigned N_EXPLICIT_DOF = N_JOINT_DOF - N_IMPLICIT_DOF;

  // get necessary data
  const VectorNd& z = opt_data.z;
  const MatrixNd& R = opt_data.R;
  const MatrixNd& G = opt_data.G;
  const VectorNd& c = opt_data.c;
  const MatrixNd& Dx = opt_data.Dx;
  const vector<unsigned>& true_indices = opt_data.true_indices;
  const vector<unsigned>& loop_indices = opt_data.loop_indices;
  const vector<double>& mu_c = opt_data.mu_c;
  const vector<double>& visc = opt_data.visc;
  const VectorNd& fext = opt_data.fext;

  // setup constraint index
  unsigned index = 0;

  // compute w
  R.mult(x, w) += z;

  // compute delta >= 0 constraints
  for (unsigned i=0; i< N_LOOPS; i++)
    fc[index++] = -w[DELTA_START+i] - INFEAS_TOL;

  // compute delta <= 1 constraints
  for (unsigned i=0; i< N_LOOPS; i++)
    fc[index++] = w[DELTA_START+i] - (double) 1.0 - INFEAS_TOL;

  // compute joint friction constraints
  for (unsigned i=0; i< N_JOINT_DOF; i++)
  {
    // original equation is mu_c ||si'*F*(fext + ff + D'*betax)|| >= ||ff||
    //                   or mu_c ||si'*F*(fext + ff + D'*betax)|| >= ||beta_x||

    // get the frictional force
    unsigned jidx = true_indices[i];
    const unsigned FIDX = (i < N_IMPLICIT_DOF) ? i : BETA_START + i - N_IMPLICIT_DOF;
    double fx = w[FIDX];

    // determine the applied forces
    w.get_sub_vec(BETA_START, BETA_START+N_EXPLICIT_DOF, ff); // get betax 
    Dx.transpose_mult(ff, DTbx);
    w.get_sub_vec(FF_START, FF_START+N_IMPLICIT_DOF, ff); // get explicit fric
    ff += fext;
    ff += DTbx;

    // two cases: joint is part of a loop or not
    if (loop_indices[jidx] != std::numeric_limits<unsigned>::max())
    {
      // get the three Z's
      const MatrixNd& Zd = opt_data.Zd[jidx];
      const MatrixNd& Z1d = opt_data.Z1d[jidx];
      const MatrixNd& Z = opt_data.Z[jidx];

      // determine delta
      const double DELTA = w[DELTA_START+loop_indices[jidx]];
 
      // compute lambda
      Zd.mult(ff, lambda) *= DELTA;
      lambda += (Z1d.mult(ff, tmp) *= ((double) 1.0 - DELTA));
      lambda += Z.mult(ff, tmp);
    }
    else
    {
      const MatrixNd& Z = opt_data.Z[jidx];
      Z.mult(ff, lambda);
    } 

    // evaluate the equation
    fc[index++] = fx*fx - mu_c[i]*lambda.norm_sq() - visc[i] - INFEAS_TOL;
  }
}

/// Calculates the joint constraint forces based on generalized force
void ArticulatedBody::calc_joint_constraint_forces(const vector<unsigned>& loop_indices, const VectorNd& delta, const vector<MatrixNd>& Zd, const vector<MatrixNd>& Z1d, const vector<MatrixNd>& Z, const VectorNd& ff) const
{
  SAFESTATIC VectorNd tmp;

  for (unsigned i=0; i< _joints.size(); i++)
  {
    // get lambda for the joint
    VectorNd& lambda = _joints[i]->lambda;

    // two cases: joint is part of a loop or not
    if (loop_indices[i] != std::numeric_limits<unsigned>::max())
    {
      // determine delta
      const double DELTA = delta[loop_indices[i]];
 
      // compute lambda
      Zd[i].mult(ff, lambda) *= DELTA;
      lambda += (Z1d[i].mult(ff, tmp) *= ((double) 1.0 - DELTA));
      lambda += Z[i].mult(ff, tmp);
    }
    else
      Z[i].mult(ff, lambda);
  }
}

// objective and inequality constraint gradients for convex optimization
void ArticulatedBody::calc_fwd_dyn_grad0(const VectorNd& x, VectorNd& grad, void* data)
{
  // setup constants
  const ABFwdDynOptData& opt_data = *(const ABFwdDynOptData*) data;

  // get necessary data
  const MatrixNd& G = opt_data.G;
  const VectorNd& c = opt_data.c;

  // objective function is quadratic
  G.mult(x, grad);
  grad += c;
}

// objective and inequality constraint gradients for convex optimization
void ArticulatedBody::calc_fwd_dyn_cJac(const VectorNd& x, MatrixNd& J, void* data)
{
  SAFESTATIC VectorNd tmpv, tmpv2, grad;
  SAFESTATIC VectorNd ff, fff, wff, wdelta, Zf, Zdf, Z1df;
  SAFESTATIC MatrixNd Rd, dX;

  // setup constants
  const ABFwdDynOptData& opt_data = *(const ABFwdDynOptData*) data;
  const unsigned N_IMPLICIT_DOF = opt_data.N_IMPLICIT_DOF;
  const unsigned N_EXPLICIT_CONSTRAINT_EQNS = opt_data.N_EXPLICIT_CONSTRAINT_EQNS;
  const unsigned N_JOINT_DOF = opt_data.N_JOINT_DOF;
  const unsigned N_LOOPS = opt_data.N_LOOPS;
  const unsigned DELTA_START = N_JOINT_DOF + N_EXPLICIT_CONSTRAINT_EQNS;
  const unsigned BETA_START = N_IMPLICIT_DOF + N_EXPLICIT_CONSTRAINT_EQNS;

  // get necessary data
  const VectorNd& z = opt_data.z;
  const MatrixNd& R = opt_data.R;
  const MatrixNd& G = opt_data.G;
  const VectorNd& c = opt_data.c;
  const VectorNd& fext = opt_data.fext;
  const vector<unsigned>& true_indices = opt_data.true_indices;
  const vector<unsigned>& loop_indices = opt_data.loop_indices;

  // resize J
  J.resize(N_LOOPS*2+N_JOINT_DOF, x.size());

  // setup the row index
  unsigned index = 0;

  // constraint delta >= 0
  for (unsigned i=0; i< N_LOOPS; i++) 
  {
    R.get_row(DELTA_START+i, grad);
    grad.negate();
    J.set_row(index++, grad);
  }

  // constraint delta <= 1
  for (unsigned i=0; i< N_LOOPS; i++) 
  {
    R.get_row(DELTA_START+i, grad);
    J.set_row(index++, grad);
  }

  // joint friction constraints
  for (unsigned i=0; i< N_JOINT_DOF; i++)
  {
    // setup numerical gradient 
//    unsigned n = x.size();
//    const double h = NEAR_ZERO;
//    const double INV_H2 = 1.0 / (h*2);
//    VectorNd xx = x;
//    grad.resize(n);
//    for (unsigned i=0; i< n; i++)
//    {
//      xx[i] += h;
//      double v1 = calc_fwd_dyn_fx(xx, m, data);
//      xx[i] -= 2*h;
//      double v2 = calc_fwd_dyn_fx(xx, m, data);
//      xx[i] += h;
//      v1 -= v2;
//      v1 *= INV_H2;
//      grad[i] = v1;
//    }
//    return;

    // get the DOF and joint index
    unsigned idx = i;
    unsigned jidx = true_indices[idx];

    // determine the friction index
    const unsigned FIDX = (idx < N_IMPLICIT_DOF) ? idx : BETA_START + idx - N_IMPLICIT_DOF;

    // get necessary constants
    const double MUCSQ = opt_data.mu_c[idx];

    // get the three Z's
    const MatrixNd& Zd = opt_data.Zd[jidx];
    const MatrixNd& Z1d = opt_data.Z1d[jidx];
    const MatrixNd& Z = opt_data.Z[jidx];

    // get components of R
    const MatrixNd& R = opt_data.R;
    const MatrixNd& Rff = opt_data.Rff;
    const MatrixNd& DxTRbetax = opt_data.DxTRbetax;

    // get components of z
    const VectorNd& zff = opt_data.zff;
    const VectorNd& zbetax = opt_data.zbetax;

    // setup variables
    Rff.mult(x, wff) += zff;
    DxTRbetax.mult(x, fff) += zbetax;
    fff += wff;                               /// fff is force
    fff += fext;

    // compute gradient of force 
    dX.copy_from(Rff) += DxTRbetax;

    // two cases: joint is part of a loop or not
    if (loop_indices[jidx] != std::numeric_limits<unsigned>::max())
    {
      // setup delta-related variables
      const unsigned LOOP_IDX = DELTA_START+loop_indices[jidx];
      R.get_sub_mat(LOOP_IDX, LOOP_IDX+1, 0, R.columns(), Rd);
      const double ZDELTA = opt_data.z[LOOP_IDX];
      Rd.mult(x, wdelta);
      const double DELTA = wdelta[0] + ZDELTA;

      // compute first component of gradient (dX' * Z' * Z * f)
      Z.mult(fff, Zf);
      Z.transpose_mult(Zf, tmpv);
      dX.transpose_mult(tmpv, grad);

      // compute second component of gradient (Rd' * f * Z' * Zd * f)
      Zd.mult(fff, Zdf);
      Z.transpose_mult(Zdf, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (fff.dot(tmpv));
      grad += tmpv2;    

      // compute third component of gradient (d * dX' * Z' * Zd * f)
      dX.transpose_mult(tmpv, tmpv2);
      tmpv2 *= DELTA;
      grad += tmpv2;

      // compute fourth component of gradient (d * dX' * Zd' * Z * f)
      Z.mult(fff, Zf);
      Zd.transpose_mult(Zf, tmpv);
      dX.transpose_mult(tmpv, tmpv2) *= DELTA;
      grad += tmpv2;

      // compute fifth component of gradient (-Rd' * f' * Z' * Z1d * f)
      Z1d.mult(fff, Z1df);
      Z.transpose_mult(Z1df, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (-fff.dot(tmpv));
      grad += tmpv2;

      // compute six and seventh components of gradient 
      // ((1-d) * dX' * [(Z' * Z1d * f) + (Z1d' * Z * f)]
      Z1d.transpose_mult(Zf, tmpv);
      Z.transpose_mult(Z1df, tmpv2);
      tmpv += tmpv2;
      dX.transpose_mult(tmpv, tmpv2) *= ((double) 1.0 - DELTA);
      grad += tmpv2;

      // compute eight component of gradient (d * Rd' * f' * Zd' * Zd * f)
      Zd.transpose_mult(Zdf, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (fff.dot(tmpv) * DELTA);
      grad += tmpv2;

      // compute ninth component of gradient (d^2 * dX' * Zd' * Zd * f)
      dX.transpose_mult(tmpv, tmpv2);
      tmpv2 *= (DELTA * DELTA);
      grad += tmpv2;

      // compute 10th and 11th components of gradient 
      // (-d * Rd' * f' * Zd' * Z1d * f) and ((1-d) * Rd' * f' * Zd' * Z1d * f)
      Zd.transpose_mult(Z1df, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (fff.dot(tmpv) * (1 - DELTA*2));
      grad += tmpv2;

      // compute 12th and 13th components of gradient
      // ((1-d) * dX * R' * Zd' * Z1d * f) and ((1-d) * d * dX' * Z1d' * Zd * f)
      Zd.transpose_mult(Z1df, tmpv);
      Z1d.transpose_mult(Zdf, tmpv2);
      tmpv += tmpv2;
      dX.transpose_mult(tmpv, tmpv2);
      tmpv2 *= (DELTA - DELTA*DELTA); 
      grad += tmpv2;

      // compute 14th component of gradient (-(1-d) * Rd' * f' * Z1d' * Z1d * f)
      Z1d.transpose_mult(Z1df, tmpv);
      Rd.get_row(0, tmpv2);
      tmpv2 *= (fff.dot(tmpv) * (DELTA - (double) 1.0));
      grad += tmpv2;

      // compute 15th component of gradient ((1-d)^2 * dX' * Z1d' * Z1d * f)
      dX.transpose_mult(tmpv, tmpv2);
      tmpv2 *= ((double) 1.0 - (double) 2.0*DELTA + DELTA*DELTA);
      grad += tmpv2;

      // scale gradient
      grad *= (-MUCSQ);

      // add in gradient due to applied frictional force
      R.get_row(FIDX, ff);
      ff *= (ff.dot(x) + z[FIDX]);
      grad += ff;
    }
    else
    {
      // compute gradient
      Z.mult(fff, Zf);
      Z.transpose_mult(Zf, tmpv);
      dX.transpose_mult(tmpv, grad) *= (-MUCSQ);

      // add in gradient due to applied frictional force
      R.get_row(FIDX, ff);
      ff *= (ff.dot(x) + z[FIDX]);
      grad += ff;
    }

    J.set_row(index++, grad);
  }
}

// objective and inequality constraint gradients for convex optimization
void ArticulatedBody::calc_fwd_dyn_hess(const VectorNd& x, double objscal, const VectorNd& hlambda, const VectorNd& nu, MatrixNd& H, void* data)
{
  SAFESTATIC MatrixNd dX, tmpM, tmpM2, tmpM3, f, Rd;
  SAFESTATIC VectorNd wdelta, wff, fff, tmpv;

  // setup constants
  const ABFwdDynOptData& opt_data = *(const ABFwdDynOptData*) data;
  const unsigned N_IMPLICIT_DOF = opt_data.N_IMPLICIT_DOF;
  const unsigned N_EXPLICIT_CONSTRAINT_EQNS = opt_data.N_EXPLICIT_CONSTRAINT_EQNS;
  const unsigned N_JOINT_DOF = opt_data.N_JOINT_DOF;
  const unsigned N_LOOPS = opt_data.N_LOOPS;
  const unsigned DELTA_START = N_JOINT_DOF + N_EXPLICIT_CONSTRAINT_EQNS;
  const unsigned BETA_START = N_IMPLICIT_DOF + N_EXPLICIT_CONSTRAINT_EQNS;

  // get necessary data
  const MatrixNd& G = opt_data.G;
  const vector<unsigned>& true_indices = opt_data.true_indices;
  const vector<unsigned>& loop_indices = opt_data.loop_indices;
  const VectorNd& fext = opt_data.fext;

  // objective function 
  H.copy_from(G) *= objscal;

  for (unsigned i=0; i< N_JOINT_DOF; i++)
  {
    // get the DOF and joint index
    unsigned idx = i; 
    unsigned jidx = true_indices[idx];

    // get the frictional index 
    const unsigned FIDX = (idx < N_IMPLICIT_DOF) ? idx : BETA_START + idx - N_IMPLICIT_DOF;

    // get necessary constants
    const double MUCSQ = opt_data.mu_c[idx];

    // get the three Z's
    const MatrixNd& Zd = opt_data.Zd[jidx];
    const MatrixNd& Z1d = opt_data.Z1d[jidx];
    const MatrixNd& Z = opt_data.Z[jidx];

    // get components of R
    const MatrixNd& R = opt_data.R;
    const MatrixNd& Rff = opt_data.Rff;
    const MatrixNd& DxTRbetax = opt_data.DxTRbetax;

    // get components of z
    const VectorNd& zff = opt_data.zff;
    const VectorNd& zbetax = opt_data.zbetax;

    // setup variables
    Rff.mult(x, wff) += zff;
    DxTRbetax.mult(x, fff) += zbetax;
    fff += wff;                               /// fff is force
    fff += fext;
    f.set(fff);

    // compute gradient of force; note: gradient of delta is just Rd
    dX.copy_from(Rff) += DxTRbetax;

    // Hessian differs depending whether joint is part of a loop 
    if (loop_indices[jidx] != std::numeric_limits<unsigned>::max())
    {
      // compute first component of Hessian (dX' * Z' * Z * dX)
      Z.mult(dX, tmpM);
      Z.transpose_mult(tmpM, tmpM2);
      dX.transpose_mult(tmpM2, tmpM3);

      // compute components necessary for delta
      const unsigned LOOP_IDX = DELTA_START+loop_indices[jidx];
      R.get_sub_mat(LOOP_IDX, LOOP_IDX+1, 0, R.columns(), Rd);
      const double ZDELTA = opt_data.z[LOOP_IDX];
      Rd.mult(x, wdelta);
      const double DELTA = wdelta[0] + ZDELTA;

      // compute second and third components of Hessian 
      // ((Rd' * f' + 2d * dX') * Z' * Zd * dX)
      Rd.transpose_mult_transpose(f,tmpM);
      MatrixNd::transpose(dX, tmpM2) *= ((double) 2.0 * DELTA);
      tmpM += tmpM2;
      tmpM.mult_transpose(Z, tmpM2);
      tmpM2.mult(Zd, tmpM);
      tmpM.mult(dX, tmpM2);
      tmpM3 += tmpM2;

      // compute first part of third Hessian component (dX' * Z' * Zd * f * dD)
      Zd.mult(f, tmpM);
      Z.transpose_mult(tmpM, tmpM2);
      dX.transpose_mult(tmpM2, tmpM);
      tmpM.mult(Rd, tmpM2) *= (double) 2.0;
      tmpM3 += tmpM2;

      // 1st part of 4th component of gradient (dX' * Zd' * Z * f * dD)
      // is already accounted for above (in first part of 3rd component)

      // 2nd part of 4th component of Hessian (d * dX' * Zd' * Z * dX)
      // is already accounted for above (in second part of 3rd component)

      // compute 1st part of 5th component of Hessian (-Rd' * f' * Z' * Z1d * f)
      // -Rd' * f' * Z' * Z1d * dX and -Rd' * dX' * Z' * Z1d * f
      Z1d.mult(dX, tmpM);
      Z.transpose_mult(tmpM, tmpM2);
      f.transpose_mult(tmpM2, tmpM);
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM2.negate();
      tmpM3 += tmpM2;

      // compute 2nd part of 5th component of Hessian -dX' * Z1d' * Z * f * Rd
      Z.mult(f, tmpM);
      Z1d.transpose_mult(tmpM, tmpM2);
      dX.transpose_mult(tmpM2, tmpM);
      tmpM.mult(Rd, tmpM2);
      tmpM2.negate();
      tmpM3 += tmpM2;

      // compute 1st part of 6th Hessian component (1-d) * dX' * Z' * Z1d * dX
      // (also accounts for 1st part of 7th Hessian component)
      dX.transpose_mult_transpose(Z, tmpM);
      tmpM.mult(Z1d, tmpM2);
      tmpM2.mult(dX, tmpM);
      tmpM *= ((double) 2.0 - DELTA*(double) 2.0);
      tmpM3 += tmpM;

      // compute 2nd part of 6th Hessian component -dD * f' * Z1d' * Z * dX
      // (also accounts for 2nd part of 7th Hessian component)
      f.transpose_mult_transpose(Z1d, tmpM);
      tmpM.mult(Z, tmpM2);
      tmpM2.mult(dX, tmpM) *= (double) -2.0;
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM3 += tmpM2;      

      // compute 8th Hessian components
      // part 1 and 2: 2 * d * Rd' * f' * Zd' * Zd * dX
      f.transpose_mult_transpose(Zd, tmpM);
      tmpM.mult(Zd, tmpM2);
      tmpM2.mult(dX, tmpM);
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM2 *= (DELTA * (double) 2.0);
      tmpM3 += tmpM2;

      // part 3: dD' * Rd * f' * Zd' * Zd * f
      tmpM.copy_from(Rd) *= Zd.mult(fff, tmpv).norm_sq();
      tmpM.transpose_mult(Rd, tmpM2);
      tmpM3 += tmpM2;

      // compute 9th component of Hessian
      // part 1: (d^2 * dX' * Zd' * Zd * dX)
      Zd.mult(dX, tmpM);
      Zd.transpose_mult(tmpM, tmpM2);
      dX.transpose_mult(tmpM2, tmpM) *= (DELTA*DELTA);
      tmpM3 += tmpM;
      
      // part 2: (2*d * dD' * f' * Zd' * Zd * dX)
      f.transpose_mult_transpose(Zd, tmpM);
      tmpM.mult(Zd, tmpM2);
      tmpM2.mult(dX, tmpM) *= ((double) 2.0 * DELTA);
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM3 += tmpM2;

      // compute first part of 10th and 11th Hessian components
      // (1-2d) * Rd' * f' * Zd' * Z1d * dX
      f.transpose_mult_transpose(Zd, tmpM);
      tmpM.mult(Z1d, tmpM2);
      tmpM2.mult(dX, tmpM) *= ((double) 1.0 - DELTA*(double) 2.0);
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM3 += tmpM2;      

      // compute second part of 10th and 11th Hessian components
      // part 2: -2 * dD' * Rd * f' * Z1d' * Zd * f 
      f.transpose_mult_transpose(Z1d, tmpM2);
      tmpM2.mult(Zd, tmpM);
      tmpM.mult(f, tmpM2) *= (double) -2.0;
      tmpM.copy_from(Rd) *= tmpM2(0,0);
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM3 += tmpM2;

      // compute 12th components of Hessian
      // part 1: (1-d) * dX' * Zd' * Z1d * dX
      // NOTE: this is used to compute part 1 of 13th component of Hessian too
      dX.transpose_mult_transpose(Zd, tmpM);
      tmpM.mult(Z1d, tmpM2);
      tmpM2.mult(dX, tmpM) *= ((double) 1.0 - DELTA*(double) 2.0 + DELTA*DELTA);
      tmpM3 += tmpM;

      // part 2: -dD' * f' * Z1d' * Zd * dX
      f.transpose_mult_transpose(Z1d, tmpM);
      tmpM.mult(Zd, tmpM2);
      tmpM2.mult(dX, tmpM).negate();
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM3 += tmpM2;

      // compute 13th components of Hessian
      // part 1: (1-d) * d * dX' * Z1d' * Zd * dX
      // (already computed above) 

      // part 2: (1-d) * dD' * f' * Zd' * Z1d * dX
      // part 3: -d * dD' * f' * Zd' * Z1d * dX
      f.transpose_mult_transpose(Zd, tmpM);
      tmpM.mult(Z1d, tmpM2);
      tmpM2.mult(dX, tmpM) *= ((double) 1.0 - DELTA*(double) 2.0);
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM3 += tmpM2;

      // compute 14th Hessian components 
      // part 1:(dD' * Rd * f' * Z1d' * Z1d * f)
      tmpM.copy_from(Rd) *= Z1d.mult(fff, tmpv).norm_sq();
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM3 += tmpM2;

      // part 2: (-(1-d) * Rd' * f' * Z1d' * Z1d * dX)
      f.transpose_mult_transpose(Z1d, tmpM);
      tmpM.mult(Z1d, tmpM2);
      tmpM2.mult(dX, tmpM) *= ((double) -2.0 - DELTA*(double) 2.0); // double
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM3 += tmpM2;

      // compute 15th Hessian components
      // part 1: ((1-d)^2 * dX' * Z1d' * Z1d * dX)
      dX.transpose_mult_transpose(Z1d, tmpM);
      tmpM.mult(Z1d, tmpM2);
      tmpM2.mult(dX, tmpM) *= ((double) 1.0 - (double) 2.0*DELTA + DELTA*DELTA);
      tmpM3 += tmpM;

      // part 2: 2*(1-d)*dD' * f' * Z1d' * Z1d * dX
      f.transpose_mult_transpose(Z1d, tmpM);
      tmpM.mult(Z1d, tmpM2);
      tmpM2.mult(dX, tmpM) *= ((double) 2.0 - DELTA*(double) 2.0);
      Rd.transpose_mult(tmpM, tmpM2);
      tmpM3 += tmpM2;

      // scale Hessian
      tmpM3 *= (-MUCSQ);

      // add in frictional component
      R.get_row(FIDX, tmpv);
      VectorN::outer_prod(tmpv, tmpv, &tmpM);
      tmpM3 += tmpM;

      // finally, scale Hessian
      tmpM3 *= hlambda[i];
      H += tmpM3;
    }
    else
    {
      // compute Hessian
      Z.mult(dX, tmpM);
      Z.transpose_mult(tmpM, tmpM2);
      dX.transpose_mult(tmpM2, tmpM3);

      // scale
      tmpM3 *= (-MUCSQ);

      // add in frictional component
      R.get_row(FIDX, tmpv);
      VectorN::outer_prod(tmpv, tmpv, &tmpM);
      tmpM3 += tmpM;

      // finally, scale Hessian
      tmpM3 *= hlambda[i];
      H += tmpM3;
    }

    // setup numerical Hessian
//    unsigned n = x.size();
//    const double h = NEAR_ZERO;
//    const double INV_H2 = 1.0 / (h*2);
//    VectorNd xx = x;
//    VectorNd v1, v2;
//    H.resize(n);
//    for (unsigned i=0; i< n; i++)
//    {
//      xx[i] += h;
//      calc_fwd_dyn_grad(xx, m, v1, data);
//      xx[i] -= 2*h;
//      calc_fwd_dyn_grad(xx, m, v2, data);
//      xx[i] += h;
//      v1 -= v2;
//      v1 *= INV_H2;
//      H.set_column(i, v1);
//    }
    // average values of the Hessian
//    for (unsigned i=0; i< n; i++)
//      for (unsigned j=i+1; j< n; j++)
//        H(i,j) = H(j,i) = 0.5*(H(i,j) + H(j,i));
  }
}

/// Transforms a force and torque on a link to a given position
SVector6 ArticulatedBody::transform_force(RigidBodyPtr link, const Vector3& x) const
{
  // get the link position
  const Vector3& com = link->get_position();

  // determine the translation
  Vector3 xlat = x - com;

  // get the link external and Coriolis forces
  Vector3 f = link->sum_forces();
  Vector3 t = link->sum_torques() - link->calc_inertial_forces();

  return SVector6(f, t - Vector3::cross(xlat, f));
} 
*/

/// Determines the loop indices corresponding to each joint and the vector of links for each joint
void ArticulatedBody::find_loops(vector<unsigned>& loop_indices, vector<vector<unsigned> >& loop_links) const
{
  SAFESTATIC vector<JointPtr> loop_joints, implicit_joints;
  queue<RigidBodyPtr> q;

  // clear vectors
  loop_indices.resize(_joints.size());
  implicit_joints.clear();

  // get all implicit joints
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->get_constraint_type() == Joint::eImplicit)
      implicit_joints.push_back(_joints[i]);

  // set all loop indices to INF (indicates no loop) initially
  for (unsigned i=0; i< _joints.size(); i++)
    loop_indices[i] = std::numeric_limits<unsigned>::max();

  // look for early exit
  if (implicit_joints.empty())
    return;

  // two cases: 1) body uses *only* implicit joints and 2) body uses 
  // explicit and implicit joints
  if (_joints.size() == implicit_joints.size())
  {
    // we're going to reset implicit_joints to hold only the joints that
    // complete loops
    implicit_joints.clear();
    for (unsigned i=0; i< _joints.size(); i++)
    {
      RigidBodyPtr inboard = _joints[i]->get_inboard_link();
      RigidBodyPtr outboard = _joints[i]->get_outboard_link();

      // check for obvious loop closure
      if (inboard->get_index() > outboard->get_index())
      {
        implicit_joints.push_back(_joints[i]);
        continue;
      }

      // check for not-so-obvious loop closure
      if (!outboard->is_enabled())
      {
        // outboard is fixed; look back to see whether one of the predecessor
        // links is fixed as well
        q.push(inboard);
        while (!q.empty())
        {
          RigidBodyPtr link = q.front();
          q.pop();
          if (!link->is_enabled())
          {
            implicit_joints.push_back(_joints[i]);
            break;
          }
          const set<JointPtr>& ij = link->get_inner_joints();
          BOOST_FOREACH(JointPtr j, ij)
            q.push(RigidBodyPtr(j->get_inboard_link()));
         }
       }
     }
  }

  // reset loop links
  loop_links.clear();
  loop_links.resize(implicit_joints.size());

  // for every kinematic loop
  for (unsigned k=0; k< implicit_joints.size(); k++)
  {
    // get the implicit joint
    JointPtr ejoint = implicit_joints[k];
    RigidBodyPtr outboard = ejoint->get_outboard_link();
    bool ground_outboard = outboard->is_ground();

    // determine all joints and links in the loop by iterating backward until
    // we get back to the first link in the loop
    loop_joints.clear();
    loop_joints.push_back(ejoint);
    loop_links[k].push_back(outboard->get_index());
    RigidBodyPtr inboard = ejoint->get_inboard_link();
    while (true)
    {
      JointPtr jx = inboard->get_inner_joint_explicit();
      loop_joints.push_back(jx);
      loop_links[k].push_back(inboard->get_index());
      inboard = jx->get_inboard_link();
 
      // check for loop termination
      if (inboard == outboard)
        break;
      if (ground_outboard && inboard->is_ground())
        break;
    }

    // reverse the vector of links so that it is (almost) sorted (all but
    // last link)
    std::reverse(loop_links[k].begin(), loop_links[k].end());
    #ifndef NDEBUG
    for (unsigned i=1; i< loop_links[k].size()-1; i++)
      assert(loop_links[k][i] > loop_links[k][i-1]);
    #endif

    // setup loop indices for each joint in the loop
    for (unsigned i=0; i< loop_joints.size(); i++)
      loop_indices[loop_joints[i]->get_index()] = k;
  }
}

/// Sets the vectors of links and joints
void ArticulatedBody::set_links_and_joints(const vector<RigidBodyPtr>& links, const vector<JointPtr>& joints)
{
  // copy the vector
  _links = links;

  // setup the link in the map 
  for (unsigned i=0; i< _links.size(); i++)
  {
    _links[i]->set_index(i);
    _links[i]->set_articulated_body(get_this());
  }

  // set vector of joints
  _joints = joints;

  // iterate over each joint
  for (unsigned i=0; i< _joints.size(); i++)
  {
    _joints[i]->set_index(i);
    _joints[i]->set_articulated_body(get_this());
  }

  compile();
}

/// Gets the number of explicit joint constraint equations
unsigned ArticulatedBody::num_constraint_eqns_explicit() const
{
  unsigned neq = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->get_constraint_type() == Joint::eExplicit)
      neq += _joints[i]->num_constraint_eqns();

  return neq;
}

/// Gets the number of implicit joint constraint equations
unsigned ArticulatedBody::num_constraint_eqns_implicit() const
{
  unsigned neq = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->get_constraint_type() == Joint::eImplicit)
      neq += _joints[i]->num_constraint_eqns();

  return neq;
}

/// Gets the number of joint degrees of freedom permitted by both implicit and explicit joint constraints
unsigned ArticulatedBody::num_joint_dof() const
{
  unsigned ndof = 0;
  for (unsigned i=0; i< _joints.size(); i++)
    ndof += _joints[i]->num_dof();
  return ndof;
}

/// Finds the joint with the given name
/**
 * \return NULL if the joint wasshared_ptr<void> not found
 */
JointPtr ArticulatedBody::find_joint(const string& jointname) const
{
  for (unsigned i=0; i< _joints.size(); i++)
    if (_joints[i]->id == jointname)
      return _joints[i];
      
  return JointPtr();
}

/// Gets the adjacent links
void ArticulatedBody::get_adjacent_links(list<sorted_pair<RigidBodyPtr> >& links) const
{
  for (unsigned i=0; i< _joints.size(); i++)
  {
    RigidBodyPtr ib = _joints[i]->get_inboard_link();
    RigidBodyPtr ob = _joints[i]->get_outboard_link();
    if (ib && ob)
      links.push_back(make_sorted_pair(ib, ob));
  }
}

/// Transforms all links in the articulated body by the given transform
/**
 * The given transformation is cumulative; the links will not necessarily be set to T.
 */
void ArticulatedBody::translate(const Origin3d& x)
{
  // apply transform to all links
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->translate(x);
}

/// Transforms all links in the articulated body by the given transform
/**
 * The given transformation is cumulative; the links will not necessarily be set to T.
 */
void ArticulatedBody::rotate(const Quatd& q)
{
  // apply transform to all links
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->rotate(q);
}

/// Calculates the combined kinetic energy of all links in this body
double ArticulatedBody::calc_kinetic_energy() 
{
  double KE = 0;
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    KE += rb->calc_kinetic_energy();

  return KE;
}

/// Resets force and torque accumulators on the body
/**
 * Force and torque accumulators on all links are reset.
 */
void ArticulatedBody::reset_accumulators()
{
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->reset_accumulators();
}

/// Finds the link with the given name
RigidBodyPtr ArticulatedBody::find_link(const string& linkid) const
{
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    if (rb->id == linkid)
      return rb;

  return RigidBodyPtr();
}

/// Updates visualization for the body
void ArticulatedBody::update_visualization()
{
  BOOST_FOREACH(RigidBodyPtr rb, _links)
    rb->update_visualization(); 
}

/// Loads a MCArticulatedBody object from an XML node
void ArticulatedBody::load_from_xml(shared_ptr<const XMLTree> node, std::map<string, BasePtr>& id_map)
{
  map<string, BasePtr>::const_iterator id_iter;

  // call parent method
  DynamicBody::load_from_xml(node, id_map);

  // don't verify the node name -- this class has derived classes
  // assert(strcasecmp(node->name().c_str(), "MCArticulatedBody") == 0);

  // determine whether to use the advanced joint friction model
  XMLAttrib* jf_attr = node->get_attrib("use-advanced-joint-friction");
  if (jf_attr)
    use_advanced_friction_model = jf_attr->get_bool_value();

  // see whether to load the model from a URDF file
  XMLAttrib* urdf_attr = node->get_attrib("urdf");
  if (urdf_attr)
  {
    // get the URDF filename
    std::string urdf_fname = urdf_attr->get_string_value();

    // load robots from the URDF
    std::string robot_name;
    std::vector<RigidBodyPtr> links;
    std::vector<JointPtr> joints; 
    if (URDFReader::read(urdf_fname, robot_name, links, joints))
      set_links_and_joints(links, joints);
    else
      std::cerr << "ArticulatedBody::load_from_xml()- unable to process URDF " << urdf_fname << std::endl;

    // do no more processing
    return;
  }

  // setup a list of joint nodes to find
  list<string> joint_node_names;
  joint_node_names.push_back("RevoluteJoint");
  joint_node_names.push_back("PrismaticJoint");
  joint_node_names.push_back("SphericalJoint");
  joint_node_names.push_back("UniversalJoint");
  joint_node_names.push_back("FixedJoint");
  joint_node_names.push_back("JointPlugin");

  // read the set of joint nodes and concatenate them into a single list
  list<shared_ptr<const XMLTree> > joint_nodes = node->find_child_nodes(joint_node_names);
  
  // read the set of link nodes
  list<shared_ptr<const XMLTree> > link_nodes = node->find_child_nodes("RigidBody");

  // if there were links read or joints read, add them 
  if (!joint_nodes.empty() || !link_nodes.empty())
  {
    // setup new lists for joints and links
    list<JointPtr> joints;
    list<RigidBodyPtr> links;

    // process all link nodes
    for (list<shared_ptr<const XMLTree> >::const_iterator i = link_nodes.begin(); i != link_nodes.end(); i++)
    {
      // get the id from the node
      XMLAttrib* id = (*i)->get_attrib("id");
      if (!id)
        throw std::runtime_error("Articulated body links are required to have unique IDs in XML");

      // get the ID
      const string& ID = id->get_string_value();

      // verify that the link was read already (if it wasn't, the problem is
      // in XMLReader)
      if ((id_iter = id_map.find(ID)) == id_map.end())
        assert(false);

      // save the link
      links.push_back(dynamic_pointer_cast<RigidBody>(id_iter->second));
    }

    // process all joint nodes in the same manner
    for (list<shared_ptr<const XMLTree> >::const_iterator i = joint_nodes.begin(); i != joint_nodes.end(); i++)
    {
      // get the id from the node
      XMLAttrib* id = (*i)->get_attrib("id");
      if (!id)
        throw std::runtime_error("Articulated body joints are required to have unique IDs in XML");

      // get the ID
      const string& ID = id->get_string_value();

      // verify that the joint was read already (if it wasn't, the problem is
      // in XMLReader)
      if ((id_iter = id_map.find(ID)) == id_map.end())
        assert(false);

      // save the joints
      joints.push_back(dynamic_pointer_cast<Joint>(id_iter->second));
    }

    // set the joints and links
    set_links_and_joints(vector<RigidBodyPtr>(links.begin(), links.end()), vector<JointPtr>(joints.begin(), joints.end()));
  }
}

/// Saves this object to a XML tree
void ArticulatedBody::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // call parent method
  DynamicBody::save_to_xml(node, shared_objects);

  // (re)set the name of this node
  node->name = "ArticulatedBody";

  // save whether to use the advanced joint friction model
  node->attribs.insert(XMLAttrib("use-advanced-joint-friction", use_advanced_friction_model));

  // add all links
  for (unsigned i=0; i< _links.size(); i++)
  {
    // create the new node in the tree
    XMLTreePtr child_node(new XMLTree("RigidBody"));
    node->add_child(child_node);

    // write to this node
    _links[i]->save_to_xml(child_node, shared_objects);
  }

  // add all joints
  for (unsigned i=0; i< _joints.size(); i++)
  {
    // create the new node in the tree -- note that save_to_xml() should set
    // the node name correctly
    XMLTreePtr child_node(new XMLTree("Joint"));
    node->add_child(child_node);

    // write to this node
    _joints[i]->save_to_xml(child_node, shared_objects);
  }
}

