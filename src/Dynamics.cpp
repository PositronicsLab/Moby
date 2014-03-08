/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iostream>
#include <limits>
#include <queue>
#include <Moby/ArticulatedBody.h>
#include <Moby/RigidBody.h>
#include <Moby/Dynamics.h>

using std::vector;
using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;

void Dynamics::calculate_accelerations(const std::vector<DynamicBodyPtr>& bodies, const std::vector<JointPtr>& world_joints)
{
  VectorNd f, iMf, lambda, JiMf, JTlambda, qdd;

  // if there are no world joints, each body's dynamics can be computed 
  // individually
  if (world_joints.empty())
  {
    BOOST_FOREACH(DynamicBodyPtr body, bodies)
      body->calc_fwd_dyn();
    return;
  }

  // setup indices for bodies

  // find *full row rank* J and \dot{J}
  find_J_and_dotJ(bodies, world_joints);

  // setup the vector of generalized forces on the system
  get_generalized_forces(bodies, f);

  // the dynamics equation:
  // M\ddot{q} + J'\lambda = f
  // J*\ddot{q} + \dot{J}*\dot{q} = 0

  // we can solve by first solving for \lambda
  // J*inv(M)*J'\lambda = J*inv(M)*f + \dot{J}*\dot{q}

  // once \lambda is known, \ddot{q} is:
  // \ddot{q} = inv(M)*(f - J'*\lambda)

  // TODO: compute J*inv(M)*f
  solve_generalized_inertia(f, iMf);

  // compute \dot{J}*\dot{q}

  // compute J*inv(M)*J'

  // solve the system for \lambda
  bool success = _LA.factor_chol(JiMJT);
  assert(success);
  (lambda = JiMf) += dotJdotq;
  _LA.solve_chol(JiMJT, lambda);

  // compute J'*\lambda

  // compute \ddot{q} = inv(M)*(f - J'*\lambda)
  f -= JTlambda;
  solve_generalized_inertia(f, qdd);

  // set the generalized accelerations
  
}

