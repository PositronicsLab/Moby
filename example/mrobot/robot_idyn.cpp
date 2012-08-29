#include <cmath>
#include <Moby/RCArticulatedBody.h>
#include <Moby/Joint.h>
#include "robot.h"

// define some necessary macros
#define Power(x,y) (std::pow((x), (y)))
#define Sqrt(x) (std::sqrt((x)))
#define Abs(x) (std::fabs((x)))
#define Sin(x) (std::sin((x)))
#define Cos(x) (std::cos((x)))

using namespace Moby;

std::map<JointPtr, VectorN> robot::calc_inv_dyn(RCArticulatedBodyPtr body, const std::map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // verify that inverse dynamics link data is the proper size
  assert(body->get_links().size() == inv_dyn_data.size());

  // declare necessary variables
  Real f0baseL, f1baseL, f2baseL;
  Real t0baseL, t1baseL, t2baseL;
  Real f0wheelcleftL, f1wheelcleftL, f2wheelcleftL;
  Real t0wheelcleftL, t1wheelcleftL, t2wheelcleftL;
  Real f0wheelcrightL, f1wheelcrightL, f2wheelcrightL;
  Real t0wheelcrightL, t1wheelcrightL, t2wheelcrightL;

  // set link transforms
  Real rxxbaseL = body->get_links()[0]->get_transform()(X,X);
  Real rxybaseL = body->get_links()[0]->get_transform()(X,Y);
  Real rxzbaseL = body->get_links()[0]->get_transform()(X,Z);
  Real txbaseL = body->get_links()[0]->get_transform()(X,W);
  Real ryxbaseL = body->get_links()[0]->get_transform()(Y,X);
  Real ryybaseL = body->get_links()[0]->get_transform()(Y,Y);
  Real ryzbaseL = body->get_links()[0]->get_transform()(Y,Z);
  Real tybaseL = body->get_links()[0]->get_transform()(Y,W);
  Real rzxbaseL = body->get_links()[0]->get_transform()(Z,X);
  Real rzybaseL = body->get_links()[0]->get_transform()(Z,Y);
  Real rzzbaseL = get_links()[0]->get_transform()(Z,Z);
  Real tzbaseL = body->get_links()[0]->get_transform()(Z,W);
  Real rxxwheelcleftL = body->get_links()[1]->get_transform()(X,X);
  Real rxywheelcleftL = body->get_links()[1]->get_transform()(X,Y);
  Real rxzwheelcleftL = body->get_links()[1]->get_transform()(X,Z);
  Real txwheelcleftL = body->get_links()[1]->get_transform()(X,W);
  Real ryxwheelcleftL = body->get_links()[1]->get_transform()(Y,X);
  Real ryywheelcleftL = body->get_links()[1]->get_transform()(Y,Y);
  Real ryzwheelcleftL = body->get_links()[1]->get_transform()(Y,Z);
  Real tywheelcleftL = body->get_links()[1]->get_transform()(Y,W);
  Real rzxwheelcleftL = body->get_links()[1]->get_transform()(Z,X);
  Real rzywheelcleftL = body->get_links()[1]->get_transform()(Z,Y);
  Real rzzwheelcleftL = get_links()[1]->get_transform()(Z,Z);
  Real tzwheelcleftL = body->get_links()[1]->get_transform()(Z,W);
  Real rxxwheelcrightL = body->get_links()[2]->get_transform()(X,X);
  Real rxywheelcrightL = body->get_links()[2]->get_transform()(X,Y);
  Real rxzwheelcrightL = body->get_links()[2]->get_transform()(X,Z);
  Real txwheelcrightL = body->get_links()[2]->get_transform()(X,W);
  Real ryxwheelcrightL = body->get_links()[2]->get_transform()(Y,X);
  Real ryywheelcrightL = body->get_links()[2]->get_transform()(Y,Y);
  Real ryzwheelcrightL = body->get_links()[2]->get_transform()(Y,Z);
  Real tywheelcrightL = body->get_links()[2]->get_transform()(Y,W);
  Real rzxwheelcrightL = body->get_links()[2]->get_transform()(Z,X);
  Real rzywheelcrightL = body->get_links()[2]->get_transform()(Z,Y);
  Real rzzwheelcrightL = get_links()[2]->get_transform()(Z,Z);
  Real tzwheelcrightL = body->get_links()[2]->get_transform()(Z,W);

  // set joint velocities
  Real qd0wheelcleftL = JointPtr(body->get_links()[1]->get_inner_joint())->get_qd()[0];
  Real qd0wheelcrightL = JointPtr(body->get_links()[2]->get_inner_joint())->get_qd()[0];

  // set link velocities
  Real v0baseL = body->get_links()[0]->get_spatial_velocity(eLink)[0];
  Real v1baseL = body->get_links()[0]->get_spatial_velocity(eLink)[1];
  Real v2baseL = body->get_links()[0]->get_spatial_velocity(eLink)[2];
  Real v3baseL = body->get_links()[0]->get_spatial_velocity(eLink)[3];
  Real v4baseL = body->get_links()[0]->get_spatial_velocity(eLink)[4];
  Real v5baseL = body->get_links()[0]->get_spatial_velocity(eLink)[5];
  Real v0wheelcleftL = body->get_links()[1]->get_spatial_velocity(eLink)[0];
  Real v1wheelcleftL = body->get_links()[1]->get_spatial_velocity(eLink)[1];
  Real v2wheelcleftL = body->get_links()[1]->get_spatial_velocity(eLink)[2];
  Real v3wheelcleftL = body->get_links()[1]->get_spatial_velocity(eLink)[3];
  Real v4wheelcleftL = body->get_links()[1]->get_spatial_velocity(eLink)[4];
  Real v5wheelcleftL = body->get_links()[1]->get_spatial_velocity(eLink)[5];
  Real v0wheelcrightL = body->get_links()[2]->get_spatial_velocity(eLink)[0];
  Real v1wheelcrightL = body->get_links()[2]->get_spatial_velocity(eLink)[1];
  Real v2wheelcrightL = body->get_links()[2]->get_spatial_velocity(eLink)[2];
  Real v3wheelcrightL = body->get_links()[2]->get_spatial_velocity(eLink)[3];
  Real v4wheelcrightL = body->get_links()[2]->get_spatial_velocity(eLink)[4];
  Real v5wheelcrightL = body->get_links()[2]->get_spatial_velocity(eLink)[5];

  // set external link forces / torques for base
  if (body->is_floating_base())
  {
    assert(inv_dyn_data.find(body->get_links().front()) != inv_dyn_data.end());
    const RCArticulatedBodyInvDynData& iddata = inv_dyn_data.find(body->get_links().front())->second;
    f0baseL = iddata._fext[X];
    f1baseL = iddata._fext[Y];
    f2baseL = iddata._fext[Z];
    t0baseL = iddata._text[X];
    t1baseL = iddata._text[Y];
    t2baseL = iddata._text[Z];
  }

  // set non-base external link forces / torques and desired joint acelerations
  assert(inv_dyn_data.find(body->get_links()[1]) != inv_dyn_data.end());
  const RCArticulatedBodyInvDynData& iddata1 = inv_dyn_data.find(body->get_links()[1])->second;
  f0wheelcleftL = iddata1._fext[X];
  f1wheelcleftL = iddata1._fext[Y];
  f2wheelcleftL = iddata1._fext[Z];
  t0wheelcleftL = iddata1._text[X];
  t1wheelcleftL = iddata1._text[Y];
  t2wheelcleftL = iddata1._text[Z];
  Real qdd0wheelcleftL = iddata1._qdd[0];
  assert(inv_dyn_data.find(body->get_links()[2]) != inv_dyn_data.end());
  const RCArticulatedBodyInvDynData& iddata2 = inv_dyn_data.find(body->get_links()[2])->second;
  f0wheelcrightL = iddata2._fext[X];
  f1wheelcrightL = iddata2._fext[Y];
  f2wheelcrightL = iddata2._fext[Z];
  t0wheelcrightL = iddata2._text[X];
  t1wheelcrightL = iddata2._text[Y];
  t2wheelcrightL = iddata2._text[Z];
  Real qdd0wheelcrightL = iddata2._qdd[0];

  // init the joint to link map
  std::map<JointPtr, VectorN> force_map;

  // declare work variables

  // compute work variables
  
  

  // write forces for inverse dynamics
  force_map[JointPtr(body->get_links()[1]->get_inner_joint())] = VectorN::construct_variable(1,(0.004044*qdd0wheelcleftL - rxxwheelcleftL*t0wheelcleftL - rxywheelcleftL*t1wheelcleftL - rxzwheelcleftL*t2wheelcleftL + 2.e-6*qd0wheelcleftL*v1wheelcleftL - 2.e-6*v0wheelcleftL*v1wheelcleftL + v1wheelcleftL*(-0.000048*v1wheelcleftL - 0.002058*v2wheelcleftL) + (0.002093*v1wheelcleftL + 0.000048*v2wheelcleftL)*v2wheelcleftL));
  force_map[JointPtr(body->get_links()[2]->get_inner_joint())] = VectorN::construct_variable(1,(0.004044*qdd0wheelcrightL - rxxwheelcrightL*t0wheelcrightL - rxywheelcrightL*t1wheelcrightL - rxzwheelcrightL*t2wheelcrightL + 2.e-6*qd0wheelcrightL*v1wheelcrightL - 2.e-6*v0wheelcrightL*v1wheelcrightL + v1wheelcrightL*(-0.000048*v1wheelcrightL - 0.002058*v2wheelcrightL) + (0.002093*v1wheelcrightL + 0.000048*v2wheelcrightL)*v2wheelcrightL));
}

