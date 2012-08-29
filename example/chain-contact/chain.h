#ifndef CHAIN_SYMBOLIC_H
#define CHAIN_SYMBOLIC_H

#include <Moby/RCArticulatedBody.h>
#include <Moby/RCArticulatedBodyInvDynAlgo.h>

namespace Moby {

class chain : public RCArticulatedBody, RCArticulatedBodyInvDynAlgo
{
  public:
    virtual ~chain() {}
    virtual RCArticulatedBody* factory() const { return new chain; }
    virtual MatrixN calc_jacobian_column(JointPtr joint, const Vector3& point);
    virtual MatrixN calc_jacobian_column(JointPtr joint, const Vector3& point, const Matrix4& base_transform, const std::map<JointPtr, VectorN>& q);
    virtual void update_link_transforms() const;
    virtual void update_link_velocities();
    virtual void apply_impulse(const Vector3& j, const Vector3& k, const Vector3& p, RigidBodyPtr link);
    virtual void calc_fwd_dyn();
    virtual std::map<boost::shared_ptr<Joint>, VectorN> calc_inv_dyn(boost::shared_ptr<RCArticulatedBody> body, const std::map<boost::shared_ptr<RigidBody>, RCArticulatedBodyInvDynData>& inv_dyn_data);
}; // end class
} // end namespace

#endif

