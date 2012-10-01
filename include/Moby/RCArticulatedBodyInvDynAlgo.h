/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _RC_ARTICULATED_BODY_INV_DYN_ALGO_H
#define _RC_ARTICULATED_BODY_INV_DYN_ALGO_H

#include <map>
#include <Moby/Base.h>
#include <Moby/MatrixN.h>
#include <Moby/Vector3.h>

namespace Moby {

class RCArticulatedBody;

/// Wrapper class for passing data to and from inverse dynamics algorithms
class RCArticulatedBodyInvDynData
{
  public:

    /// External force applied to this link's center-of-mass (in global frame)
    Vector3 fext;

    /// External torque applied to this link's center-of-mass (in global frame)
    Vector3 text;

    /// Inner joint acceleration (desired)
    VectorN qdd;
};
  
/// Abstract class for performing inverse dynamics computation on a reduced-coordinate articulatedc body
class RCArticulatedBodyInvDynAlgo : public virtual Base
{
  public:
    RCArticulatedBodyInvDynAlgo() { }
    virtual ~RCArticulatedBodyInvDynAlgo() {}

    /// Classes must implement the method below
    /**
     * Given the current base position and velocity, the joint positions and 
     * joint velocities, the external forces on all links (including the 
     * base), and the desired joint and base accelerations, the inverse 
     * dynamics method must compute the forces applied to the joint actuators
     * to achieve these accelerations.
     * \return a map from joints to computed actuator forces
     */
    virtual std::map<boost::shared_ptr<Joint>, VectorN> calc_inv_dyn(boost::shared_ptr<RCArticulatedBody> body, const std::map<boost::shared_ptr<RigidBody>, RCArticulatedBodyInvDynData>& inv_dyn_data) = 0;
}; // end class
} // end namespace

#endif
