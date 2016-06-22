/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <sys/times.h>
#include <iomanip>
#include <boost/foreach.hpp>
#include <boost/algorithm/minmax_element.hpp>
#include <limits>
#include <set>
#include <cmath>
#include <numeric>
#include <Moby/LP.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/Constants.h>
#include <Moby/Constraint.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/RigidBody.h>
#include <Moby/Log.h>
#include <Moby/XMLTree.h>
#include <Moby/ImpactToleranceException.h>
#include <Moby/LCPSolverException.h>
#include <Moby/QP.h>
#include <Moby/ImpactConstraintHandler.h>

#undef USE_QPOASES
#ifdef USE_QPOASES
#include <Moby/qpOASES.h>
#endif


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

/// Get the total number of variables
unsigned ImpactConstraintHandler::num_variables(const vector<Constraint*>& constraints)
{
  unsigned n_vars = 0;
  for (unsigned m=0; m< constraints.size(); m++)
    n_vars += constraints[m]->num_variables(); 
  return n_vars;
}

/// Gets the number of inequality constraints
unsigned ImpactConstraintHandler::num_inequality_constraints(const vector<Constraint*>& constraints)
{
  unsigned n = 0; 

  for (unsigned i=0; i< constraints.size(); i++)
    for (unsigned j=0; j< constraints[i]->num_constraint_equations(); j++)
      if (constraints[i]->get_constraint_equation_type(j) == Constraint::eInequality)
        n++;

  return n;
}

/// Gets the number of equality constraints
unsigned ImpactConstraintHandler::num_equality_constraints(const vector<Constraint*>& constraints)
{
  unsigned n = 0; 

  for (unsigned i=0; i< constraints.size(); i++)
    for (unsigned j=0; j< constraints[i]->num_constraint_equations(); j++)
      if (constraints[i]->get_constraint_equation_type(j) == Constraint::eEquality)
        n++;

  return n;
}

/**
 * TODO: 
 * 4. Fix LHS
 **/

/// Computes restitution for each constraint 
bool ImpactConstraintHandler::apply_restitution(const vector<Constraint*>& constraints, VectorNd& x)
{
  bool applied = false; 

  // setup a work vector
  VectorNd workv;

  // setup a mapping
  map<shared_ptr<DynamicBodyd>, VectorNd> gj;

  // loop through all constraints
  for (unsigned i=0, j=0; i< constraints.size(); i++)
  {
    // get the total number of variables and impulsive variables
    unsigned nvars = constraints[i]->num_variables();

    // get the variables out
    x.get_sub_vec(j, j+nvars, workv);

    // apply restitution 
    constraints[i]->apply_restitution(workv);
    if (workv.norm() > NEAR_ZERO)
      applied = true;

    // update x
    x.set_sub_vec(j, workv);

    // update j 
    j += nvars;
  }

  return applied;
}


/// Applies impulses using the output of the QP 
void ImpactConstraintHandler::apply_impulses(const vector<Constraint*>& constraints, const VectorNd& x)
{
  // setup a work vector
  VectorNd workv;

  // setup a mapping
  map<shared_ptr<DynamicBodyd>, VectorNd> gj;

  // loop through all constraints
  for (unsigned i=0, j=0; i< constraints.size(); i++)
  {
    // get the total number of variables and impulsive variables
    unsigned nvars = constraints[i]->num_variables();
    unsigned nimp = constraints[i]->num_impulsive_variables();

    // get the variables out
    x.get_sub_vec(j, j+nimp, workv);

    // apply the impulses
    constraints[i]->apply_impulses(workv, gj);

    // update j 
    j += nvars;
  }

  // apply the impulses
  for (map<shared_ptr<DynamicBodyd>, VectorNd>::const_iterator i=gj.begin(); i != gj.end(); i++)
    i->first->apply_generalized_impulse(i->second);
}

/// Gets the maximum constraint violation
double ImpactConstraintHandler::calc_max_constraint_violation(const vector<Constraint*>& constraints)
{
  double violation = 0.0;

  for (unsigned i=0; i< constraints.size(); i++)
    for (unsigned j=0; j< constraints[i]->num_constraint_equations(); j++)
    {
      double rhs = constraints[i]->get_constraint_rhs(j, 0.0);
      if (constraints[i]->get_constraint_equation_type(j) == Constraint::eEquality)
        violation = std::max(violation, std::fabs(rhs));
      else
        violation = std::max(violation, -rhs);
    }

  return violation;
}

/// Computes the quadratic term of a matrix
void ImpactConstraintHandler::compute_quadratic_matrix(const vector<Constraint*>& constraints, unsigned N_VARS, MatrixNd& H)
{
  // get all super bodies
  vector<shared_ptr<DynamicBodyd> > supers;
  get_super_bodies(constraints.begin(), constraints.end(), std::back_inserter(supers));
  
  // store current velocities
  vector<VectorNd> qd(supers.size());
  for (unsigned i=0; i< supers.size(); i++)
    supers[i]->get_generalized_velocity(DynamicBodyd::eSpatial, qd[i]);

  // compute the quadratic term matrix
  H.resize(N_VARS,N_VARS);
  for (unsigned m=0, o=0; m< constraints.size(); m++)
  {
    for (unsigned n=0; n< constraints[m]->num_impulsive_variables(); n++)
    {
      // get the column of H
      SharedVectorNd Hcol = H.column(o);

      // measure constraint velocities before applying the impulse
      for (unsigned i=0, j=0; i< constraints.size(); i++)
        for (unsigned k=0; k< constraints[i]->num_variables(); k++)
          Hcol[j++] = -constraints[i]->calc_projected_vel(k);

      // apply a test impulse at the constraint
      constraints[m]->apply_test_impulse(n); 

      // measure constraint velocities after applying the impulse
      for (unsigned i=0, j=0; i< constraints.size(); i++)
        for (unsigned k=0; k< constraints[i]->num_variables(); k++)
        {
          if (o <= j)
            Hcol[j++] += constraints[i]->calc_projected_vel(k);
          else
          {
            Hcol[j] = H(o,j);
            j++;
          }
        }

      // update o
      o++;
    }

    // update o to bypass slack variables
    o += (constraints[m]->num_variables() - constraints[m]->num_impulsive_variables());
  }

  // restore velocities
  for (unsigned i=0; i< supers.size(); i++)
    supers[i]->set_generalized_velocity(DynamicBodyd::eSpatial, qd[i]);
}

/// Calculates inequality constraint terms for the QP (M*x >= q)
void ImpactConstraintHandler::compute_inequality_terms(const vector<Constraint*>& constraints, const MatrixNd& H, unsigned N_INEQ_CONSTRAINTS, double inv_dt, MatrixNd& M, VectorNd& q)
{
  // get the number of variables
  const unsigned N_VARS = H.rows();

  // resize M and q
  M.resize(N_INEQ_CONSTRAINTS, N_VARS);
  q.resize(N_INEQ_CONSTRAINTS);

  for (unsigned i=0, k=0, m=0; i< constraints.size(); i++)
  {
    for (unsigned j=0; j< constraints[i]->num_constraint_equations(); j++)
    {
      if (constraints[i]->get_constraint_equation_type(j) == Constraint::eInequality)
      {
        // get the row of M
        SharedVectorNd Mrow = M.row(k);
        SharedVectorNd Msub = Mrow.segment(m, m+constraints[i]->num_variables());

        // setup the entry of q
        q[k] = constraints[i]->get_constraint_rhs(j, inv_dt);

        // get the constraint index
        if (constraints[i]->get_constraint_coeff_type(j) == Constraint::eVelocityConstraint)
          Mrow = H.row(constraints[i]->get_velocity_constraint_index(j)+m);
        else
        {
          Mrow.set_zero();
          constraints[i]->get_impulse_constraint_coeffs(j, Msub);
        }

        // update the matrix row 
        k++;
      }
    }

    // update m
    m += constraints[i]->num_variables();
  }
}

/// Calculates equality constraint terms for the QP (A*x = b)
void ImpactConstraintHandler::compute_equality_terms(const vector<Constraint*>& constraints, const MatrixNd& H, unsigned N_EQ_CONSTRAINTS, double inv_dt, MatrixNd& A, VectorNd& b)
{
  // get the number of variables
  const unsigned N_VARS = H.rows();

  // resize the matrix and vector
  A.resize(N_EQ_CONSTRAINTS, N_VARS);
  b.resize(N_EQ_CONSTRAINTS);

  // populate
  for (unsigned i=0, k=0, m=0; i< constraints.size(); i++)
  {
    for (unsigned j=0; j< constraints[i]->num_constraint_equations(); j++)
    {
      if (constraints[i]->get_constraint_equation_type(j) == Constraint::eEquality)
      {
        // get the row of A 
        SharedVectorNd Arow = A.row(k);
        SharedVectorNd Asub = Arow.segment(m, m+constraints[i]->num_variables());

        // setup the entry of b 
        b[k] = constraints[i]->get_constraint_rhs(j, inv_dt);

        // get the constraint index
        if (constraints[i]->get_constraint_coeff_type(j) == Constraint::eVelocityConstraint)
          Arow = H.row(constraints[i]->get_velocity_constraint_index(j)+m);
        else
          constraints[i]->get_impulse_constraint_coeffs(j, Asub);

        // update the matrix row 
        k++;
      }
    }

    // update m
    m += constraints[i]->num_variables();
  }
}

/// Computes the linear term of a QP/LCP 
void ImpactConstraintHandler::compute_linear_term(const vector<Constraint*>& constraints, unsigned N_VARS, double inv_dt, VectorNd& c)
{
  c.resize(N_VARS);
  for (unsigned i=0, j=0; i< constraints.size(); i++)
    for (unsigned k=0; k< constraints[i]->num_variables(); k++)
      c[j++] = constraints[i]->calc_projected_stab_vel(k, inv_dt);
}

/// Determine the number of slackable equality constraints
unsigned ImpactConstraintHandler::num_slackable_constraints(const vector<Constraint*>& constraints)
{
  unsigned n = 0;

  for (unsigned i=0; i< constraints.size(); i++)
    for (unsigned k=0; k< constraints[i]->num_constraint_equations(); k++)
    {
      if (constraints[i]->get_constraint_equation_type(k) == Constraint::eEquality && constraints[i]->is_constraint_slackable(k))
        n++;
    }

  return n;
}        


/// Forms and solves the impact problem
void ImpactConstraintHandler::form_and_solve(const vector<Constraint*>& constraints, double inv_dt, unsigned N_VARS, MatrixNd& H, VectorNd& z)
{
  const double INF = (double) std::numeric_limits<double>::max();

  // determine number of inequality and equality constraints 
  const unsigned N_INEQ_CONSTRAINTS = num_inequality_constraints(constraints);
  const unsigned N_EQ_CONSTRAINTS = num_equality_constraints(constraints);

  // compute the quadratic matrix if necessary
  if (H.rows() == 0)
  {
    compute_quadratic_matrix(constraints, N_VARS, _H);

    // add in damping terms
    for (unsigned m=0, o=0; m< constraints.size(); m++)
    {
      for (unsigned n=0; n< constraints[m]->num_impulsive_variables(); n++, o++)
      {
        // get the damping term (if any)
        double damping = constraints[m]->get_damping(n);
        if (damping <= 0.0)
          continue;

        // update H
        H(o,o) += damping;
      }

      // update o to bypass slack variables
      o += (constraints[m]->num_variables() - constraints[m]->num_impulsive_variables());
    }
  }

  // compute the linear term vector
  compute_linear_term(constraints, N_VARS, inv_dt, _c);

  // setup the inequality constraints matrix and vector
  compute_inequality_terms(constraints, H, N_INEQ_CONSTRAINTS, inv_dt, _M, _q);

  // setup the equality constraints matrix and vector
  compute_equality_terms(constraints, H, N_EQ_CONSTRAINTS, inv_dt, _A, _b);

  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::form_and_solve()" << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "H: " << std::endl  << H;
  FILE_LOG(LOG_CONSTRAINT) << "c: " << _c << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "M: " << std::endl  << _M;
  FILE_LOG(LOG_CONSTRAINT) << "q: " << _q << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "A: " << std::endl << _A;
  FILE_LOG(LOG_CONSTRAINT) << "b: " << _b << std::endl;

  // determine the number of slackable constraints 
  unsigned N_SLACKABLE_EQ_CONSTRAINTS = (N_EQ_CONSTRAINTS == 0) ? 0 : num_slackable_constraints(constraints);

  // see whether we can solve using an LCP formulation
  if (N_EQ_CONSTRAINTS == 0)
  {
    // solve using an LCP formulation
    _MM.set_zero(N_VARS + N_INEQ_CONSTRAINTS, N_VARS + N_INEQ_CONSTRAINTS);
    _qq.resize(_MM.rows());

    // setup MM and qq
    _MM.block(0, N_VARS, 0, N_VARS) = _H;
    _MM.block(N_VARS, _MM.rows(), 0, N_VARS) = _M;
    SharedMatrixNd MT = _MM.block(0, N_VARS, N_VARS, _MM.rows());
    MatrixNd::transpose(_M, MT);
    MT.negate();
    _qq.segment(0, N_VARS) = _c;
    _qq.segment(N_VARS, _qq.size()) = _q;
    _qq.segment(N_VARS, _qq.size()).negate();

    // solve using the fast regularized solver first, then fall back to the
    // Lemke solver
    if (!_lcp.lcp_fast_regularized(_MM, _qq, z, -20, 4, -8))
    {
      // Lemke does not like warm starting
      z.set_zero();

      FILE_LOG(LOG_CONSTRAINT) << "-- fast LCP solver failed; dropping back to Lemke" << std::endl;
      if (!_lcp.lcp_lemke_regularized(_MM, _qq, z))
        throw LCPSolverException();
    }
    FILE_LOG(LOG_CONSTRAINT) << "LCP solution z: " << z << std::endl;
    FILE_LOG(LOG_CONSTRAINT) << "LCP pivots: " << _lcp.pivots << std::endl;
  }
  else
  {
    // setup the lower and upper bounds variables
    _lb.resize(N_VARS);
    _ub.resize(N_VARS);
    for (unsigned i=0, j=0; i< constraints.size(); i++)
      for (unsigned k=0; k< constraints[i]->num_variables(); k++)
      {
        _lb[j] = constraints[i]->get_lower_variable_limit(k);
        _ub[j] = constraints[i]->get_upper_variable_limit(k);
        j++;
      } 

    #ifdef USE_QPOASES
    // save some variables
    static qpOASES qp;

    // special consideration needed when slackable constraints present 
    if (N_SLACKABLE_EQ_CONSTRAINTS == 0)
    {
      // add tolerance to inequality constraints
      double TOL = (inequality_tolerance >= 0.0) ? inequality_tolerance : (_H.rows() + _q.rows()) * (_H.norm_inf() + _M.norm_inf()) * std::numeric_limits<double>::epsilon();
      for (ColumnIteratord iter = _q.begin(); iter != _q.end(); iter++)
        *iter -= TOL;

      // add tolerance to lower bounds (if necessary)
      for (ColumnIteratord iter = _lb.begin(); iter != _lb.end(); iter++)
        if (*iter > -INF)
          *iter -= TOL;

      // add tolerance to upper bounds (if necessary)
      for (ColumnIteratord iter = _ub.begin(); iter != _ub.end(); iter++)
        if (*iter < INF)
          *iter += TOL;

      // no slackable equality constraints; try solving QP first w/o any
      // tolerance in the constraints
      bool result = qp.qp_activeset(_H, _c, _lb, _ub, _M, _q, _A, _b, z);
      if (!result)
      {
        FILE_LOG(LOG_CONSTRAINT) << "QP activeset solution failed without tolerance in the constraints" << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << " -- solving LP to find feasible solution" << std::endl;

        // determine new number of variables (one for inequality constraints,
        // one for each equality constraint)
        unsigned NEW_VARS = N_VARS + 1 + N_EQ_CONSTRAINTS*2;

        // augment M
        _Maug.resize(_M.rows(), NEW_VARS);
        _Maug.block(0, _M.rows(), 0, N_VARS) = _M;
        _Maug.block(0, _M.rows(), N_VARS, NEW_VARS).set_zero();
        _Maug.column(N_VARS).set_one();

        // augment A
        _Aaug.resize(_A.rows(), NEW_VARS);
        _Aaug.block(0, _A.rows(), 0, N_VARS) = _A;
        _Aaug.block(0, _A.rows(), N_VARS, NEW_VARS).set_zero();
        for (unsigned i=0, j=N_VARS+1; i< _A.rows(); i++)
        {
          _Aaug(i,j++) = 1.0;
          _Aaug(i,j++) = -1.0;
        }

        // augment lb
        _lbaug.resize(NEW_VARS);
        _lbaug.segment(0, _lb.size()) = _lb;
        _lbaug.segment(_lb.size(), NEW_VARS).set_zero();

        // augment ub
        _ubaug.resize(NEW_VARS);
        _ubaug.segment(0, _ub.size()) = _ub;
        _ubaug.segment(_ub.size(), NEW_VARS).set_one() *= INF;

        // create a c for the QP
        _caug.resize(NEW_VARS);
        _caug.segment(0, N_VARS).set_zero();
        _caug.segment(N_VARS, NEW_VARS).set_one();

        FILE_LOG(LOG_CONSTRAINT) << "H (augmented): " << std::endl  << _Haug;
        FILE_LOG(LOG_CONSTRAINT) << "c (augmented): " << _caug << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "M (augmented): " << std::endl  << _Maug;
        FILE_LOG(LOG_CONSTRAINT) << "q (augmented): " << _q << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "A (augmented): " << std::endl << _Aaug;
        FILE_LOG(LOG_CONSTRAINT) << "b (augmented): " << _b << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "lb (augmented): " << _lbaug << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "ub (augmented): " << _ubaug << std::endl;

        // solve the LP, using zero for z
        z.set_zero(NEW_VARS);
        result = LP::lp_simplex(_caug, _Maug, _q, _Aaug, _b, _lbaug, _ubaug, z);
        assert(result);
        FILE_LOG(LOG_CONSTRAINT) << " -- LP solution: " << z << std::endl;

        // modify c for re-solving the QP
        _caug.segment(0, N_VARS) = _c;
        _caug.segment(N_VARS, NEW_VARS).set_zero();

        // modify H for re-solving the QP
        _Haug.resize(NEW_VARS, NEW_VARS);
        _Haug.block(0, N_VARS, 0, N_VARS) = _H;
        _Haug.block(N_VARS, NEW_VARS, 0, N_VARS).set_zero();
        _Haug.block(0, N_VARS, N_VARS, NEW_VARS).set_zero();
        _Haug.block(N_VARS, NEW_VARS, N_VARS, NEW_VARS).set_zero();

        // modify lb and ub for the new variables
        for (unsigned i=N_VARS; i< NEW_VARS; i++)
        {
          _lbaug[i] = z[i];
          _ubaug[i] = z[i];
        }

        FILE_LOG(LOG_CONSTRAINT) << "H (augmented): " << std::endl  << _Haug;
        FILE_LOG(LOG_CONSTRAINT) << "c (augmented): " << _caug << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "M (augmented): " << std::endl  << _Maug;
        FILE_LOG(LOG_CONSTRAINT) << "q (augmented): " << _q << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "A (augmented): " << std::endl << _Aaug;
        FILE_LOG(LOG_CONSTRAINT) << "b (augmented): " << _b << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "lb (augmented): " << _lbaug << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "ub (augmented): " << _ubaug << std::endl;

        // resolve the QP
        result = qp.qp_activeset_regularized(_Haug, _caug, _lbaug, _ubaug, _Maug, _q, _Aaug, _b, z);
        assert(result);
        FILE_LOG(LOG_CONSTRAINT) << "result of QP activeset with constraints: " << result << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "robust QP activeset solution: " << z << std::endl;
      }
      else
      {
        FILE_LOG(LOG_CONSTRAINT) << "QP activeset solution found" << std::endl;
        FILE_LOG(LOG_CONSTRAINT) << "QP activeset solution: " << z << std::endl;
      }
    } 
    else
    {
      // there are slackable equality constraints; find a feasible solution
      // to the problem without minimizing slack on the slackable equality 
      // constraints  

      // first rearrange A/b so that slackable equality constraints are on the
      // bottom; 1. copy slackable constraints first
      _workM.resize(_A.rows(), _A.columns());
      _workv.resize(_b.size());
      for (unsigned i=0, j=0, m=0; i< constraints.size(); i++)
        for (unsigned k=0; k< constraints[i]->num_constraint_equations(); k++)
      {
        if (constraints[i]->get_constraint_equation_type(k) == Constraint::eEquality)
        {
          if (constraints[i]->is_constraint_slackable(k))
          {
            _workM.row(j) = _A.row(m); 
            _workv[j] = _b[m];
            j++;
          }

          // update m
          m++;
        }
      }

      // rearrange A/b so that slackable equality constraints are on the
      // bottom; 1. copy non-slackable constraints next 
      for (unsigned i=0, j=N_SLACKABLE_EQ_CONSTRAINTS, m=0; i< constraints.size(); i++)
        for (unsigned k=0; k< constraints[i]->num_constraint_equations(); k++)
      {
        if (constraints[i]->get_constraint_equation_type(k) == Constraint::eEquality)
        {
          if (!constraints[i]->is_constraint_slackable(k))
          {
            _workM.row(j) = _A.row(m); 
            _workv[j] = _b[m];
            j++;
          }

          // update m
          m++;
        }
      }

      // determine new number of variables (one for inequality constraints,
      // two for unslackable equality constraints, two for each slackable
      // equality constraint)
      unsigned NEW_VARS = N_VARS + 1 + N_EQ_CONSTRAINTS*2;

      // augment M
      _Maug.resize(_M.rows(), NEW_VARS);
      _Maug.block(0, _M.rows(), 0, N_VARS) = _M;
      _Maug.block(0, _M.rows(), N_VARS, NEW_VARS).set_zero();
      _Maug.column(N_VARS).set_one();

      // augment A
      unsigned N_NON_SLACKABLE_EQ_CONSTRAINTS = N_EQ_CONSTRAINTS - N_SLACKABLE_EQ_CONSTRAINTS;
      _Aaug.resize(_A.rows(), NEW_VARS);
      _Aaug.block(0, _A.rows(), 0, N_VARS) = _A;
      _Aaug.block(0, _A.rows(), N_VARS, NEW_VARS).set_zero();
      for (unsigned i=0, j=N_VARS+1; i< _A.rows(); i++)
      {
        _Aaug(i,j++) = 1.0;
        _Aaug(i,j++) = -1.0;
      }

      // augment lb
      _lbaug.resize(NEW_VARS);
      _lbaug.segment(0, _lb.size()) = _lb;
      _lbaug.segment(_lb.size(), NEW_VARS).set_zero();

      // augment ub
      _ubaug.resize(NEW_VARS);
      _ubaug.segment(0, _ub.size()) = _ub;
      _ubaug.segment(_ub.size(), NEW_VARS).set_one() *= INF;

      // create a H for the QP
      _Haug.set_zero(NEW_VARS, NEW_VARS); 

      // create a c for the QP
      _caug.resize(NEW_VARS);
      _caug.segment(0, N_VARS).set_zero();
      _caug.segment(N_VARS, N_VARS+3).set_one();
      _caug.segment(N_VARS+3, NEW_VARS).set_zero(); // don't penalize slack
                                                    // on slackable equality
                                                    // constraints... yet.

      // solve the QP, using zero for z
      z.set_zero(NEW_VARS);
      bool result = qp.qp_activeset(_Haug, _caug, _lbaug, _ubaug, _Maug, _q, _Aaug, _b, z);
      assert(result);
      FILE_LOG(LOG_CONSTRAINT) << " -- LP solution (1): " << z << std::endl;

      // hold the slack variables for the inequality constraints and 
      // unslackable equality constraints at their present values and now
      // minimize slack on slackable equality constraints

      // create a H for the QP, penalizing only slackable equality constraints
      _Haug.set_zero(NEW_VARS, NEW_VARS); 
      _Haug.block(N_VARS+3, NEW_VARS, N_VARS+3, NEW_VARS).set_identity();

      // create a c for the QP, penalizing only slackable equality constraints
      _caug.segment(N_VARS, N_VARS+3).set_zero();
      _caug.segment(N_VARS+3, NEW_VARS).set_one();

      // modify lb and ub for the subset of new variables 
      for (unsigned i=N_VARS; i< N_VARS+3; i++)
        _lbaug[i] = _ubaug[i] = z[i];

      // resolve the QP
      result = qp.qp_activeset(_Haug, _caug, _lbaug, _ubaug, _Maug, _q, _Aaug, _b, z);
      assert(result);
      FILE_LOG(LOG_CONSTRAINT) << " -- LP solution (2): " << z << std::endl;

      // modify c for re-solving the QP - no longer penalize slackness
      _caug.segment(0, N_VARS) = _c;
      _caug.segment(N_VARS, NEW_VARS).set_zero();

      // modify H for re-solving the QP - no longer penalize slackness
      _Haug.block(0, N_VARS, 0, N_VARS) = _H;
      _Haug.block(N_VARS, NEW_VARS, 0, N_VARS).set_zero();
      _Haug.block(0, N_VARS, N_VARS, NEW_VARS).set_zero();
      _Haug.block(N_VARS, NEW_VARS, N_VARS, NEW_VARS).set_zero();

      // modify lb and ub for the new variables
      for (unsigned i=N_VARS; i< NEW_VARS; i++)
        _lbaug[i] = _ubaug[i] = z[i];

      // resolve the QP
      result = qp.qp_activeset(_Haug, _caug, _lbaug, _ubaug, _Maug, _q, _Aaug, _b, z);
      assert(result);
      FILE_LOG(LOG_CONSTRAINT) << "result of QP activeset with constraints: " << result << std::endl;
      FILE_LOG(LOG_CONSTRAINT) << "robust QP activeset solution: " << z << std::endl;
    }
    #else
    if (N_SLACKABLE_EQ_CONSTRAINTS == 0)
    {
      // convert the QP to an LCP
      VectorNd zlo, zhi;
      QP::convex_qp_to_glcp(_H, _c, _M, _q, _A, _b, _lb, _ub, _MM, _qq, zlo, zhi);

      // solve using the fast regularized solver first, then fall back to the
      // Lemke solver
      if (!_lcp.mlcp_fast_regularized(_MM, _qq, zlo, zhi, z, -16, 4, 0))
        throw LCPSolverException();
    }
    #endif
  }
}

/// Solves the quadratic program (does all of the work) as an LCP
/**
 * \note this is the version without joint friction forces
 * \param z the solution is returned here; zeros are returned at appropriate
 *        places for inactive contacts
 */
void ImpactConstraintHandler::apply_model_to_connected_constraints(const vector<Constraint*>& constraints, double inv_dt)
{
  // mark starting time
  tms cstart;
  clock_t start = times(&cstart);

  // determine the total number of variables
  const unsigned N_VARS = num_variables(constraints);

  // init H to indicate that H must be constructed
  _H.resize(0,0);

  // prepare to solve the problem using warmstarting if possible
  VectorNd z;
  if (N_VARS == _zlast.size())
    z = _zlast;

  // form and solve the problem
  form_and_solve(constraints, inv_dt, N_VARS, _H, z);

  // apply impulses
  apply_impulses(constraints, z);

  // get the constraint violation resulting from solving the QP 
  double max_vio_pre = calc_max_constraint_violation(constraints);
  
  // apply restitution
  VectorNd zplus = z;
  if (apply_restitution(constraints, zplus))
  {
    apply_impulses(constraints, zplus);

    // update z
    z += zplus;
    FILE_LOG(LOG_CONSTRAINT) << "z (after restitution): " << z << std::endl;

    // see whether we need to solve another impact problem
    // -- if the constraint violation is negative and the constraint violation
    //    is larger than after solving the inelastic impact QP, another impact 
    //    problem must be solved 
    double max_vio_post = calc_max_constraint_violation(constraints);
    if (max_vio_post > max_vio_pre + NEAR_ZERO)
    {
      // form and solve the problem
      form_and_solve(constraints, inv_dt, N_VARS, _H, z);

      // apply impulses
      apply_impulses(constraints, z);
      FILE_LOG(LOG_CONSTRAINT) << "z (from second impact phase): " << z << std::endl;
    }
  }

  // save the last solution
  _zlast = z;

  // get the elapsed time
  const long TPS = sysconf(_SC_CLK_TCK);
  tms cstop;
  clock_t stop = times(&cstop);
  double elapsed = (double) (stop - start)/TPS;
  FILE_LOG(LOG_CONSTRAINT) << "Elapsed time: " << elapsed << std::endl;
  FILE_LOG(LOG_CONSTRAINT) << "ImpactConstraintHandler::solve_qp_work() exited" << std::endl;
}


