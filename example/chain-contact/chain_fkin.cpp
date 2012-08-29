#include <cmath>
#include <Moby/RCArticulatedBody.h>
#include <Moby/Joint.h>
#include "chain.h"

// define some necessary macros
#define Power(x,y) (std::pow((x), (y)))
#define Sqrt(x) (std::sqrt((x)))
#define Abs(x) (std::fabs((x)))
#define Sin(x) (std::sin((x)))
#define Cos(x) (std::cos((x)))

using namespace Moby;

// write the factory so that we can load this class as a plugin
namespace Moby {
  extern "C" boost::shared_ptr<RCArticulatedBody> factory() { return boost::shared_ptr<RCArticulatedBody>(new chain); }}

void chain::update_link_transforms() const
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // setup base transform
  Real rxxbaseL = get_links().front()->get_transform()(X,X);
  Real rxybaseL = get_links().front()->get_transform()(X,Y);
  Real rxzbaseL = get_links().front()->get_transform()(X,Z);
  Real txbaseL = get_links().front()->get_transform()(X,W);
  Real ryxbaseL = get_links().front()->get_transform()(Y,X);
  Real ryybaseL = get_links().front()->get_transform()(Y,Y);
  Real ryzbaseL = get_links().front()->get_transform()(Y,Z);
  Real tybaseL = get_links().front()->get_transform()(Y,W);
  Real rzxbaseL = get_links().front()->get_transform()(Z,X);
  Real rzybaseL = get_links().front()->get_transform()(Z,Y);
  Real rzzbaseL = get_links().front()->get_transform()(Z,Z);
  Real tzbaseL = get_links().front()->get_transform()(Z,W);

  // setup q
  Real q0l1L = boost::shared_ptr<Joint>(get_links()[1]->get_inner_joint())->q[0];
  Real q0l2L = boost::shared_ptr<Joint>(get_links()[2]->get_inner_joint())->q[0];

  // declare work variables
  Real Compile_29466;
  Real Compile_29467;
  Real Compile_29468;
  Real Compile_29486;
  Real Compile_29487;
  Real Compile_29488;
  Real Compile_29489;
  Real Compile_29490;
  Real Compile_29491;
  Real Compile_29492;
  Real Compile_29493;
  Real Compile_29494;
  Real Compile_29495;

  // compute work variables
  Compile_29468 = Sin(q0l1L);Compile_29466 = Cos(q0l1L);Compile_29467 = Compile_29466;
  Compile_29487 = Cos(q0l2L);Compile_29489 = Sin(q0l1L);Compile_29486 = Cos(q0l1L);Compile_29490 = Sin(q0l2L);Compile_29493 = -Compile_29487*Compile_29489;Compile_29494 = -Compile_29486*Compile_29490;Compile_29495 = Compile_29493 + Compile_29494;Compile_29488 = Compile_29486*Compile_29487;Compile_29491 = -Compile_29489*Compile_29490;Compile_29492 = Compile_29488 + Compile_29491;

  // compute transforms
  get_links()[1]->set_transform(Matrix4(Compile_29467,-Compile_29468,0.,2.5*Compile_29468,Compile_29468,Compile_29467,0.,-2.5*Compile_29466,0.,0.,((1. - Compile_29466) + Compile_29466),0.));
  get_links()[2]->set_transform(Matrix4(Compile_29492,Compile_29495,0.,5.*Compile_29489 - 2.5*Compile_29495,Compile_29487*Compile_29489 + Compile_29486*Compile_29490,Compile_29492,0.,-5.*Compile_29486 - 2.5*Compile_29492,0.,0.,((1. - Compile_29486) + Compile_29486)*((1. - Compile_29487) + Compile_29487),0.));
}

MatrixN chain::calc_jacobian_column(JointPtr joint, const Vector3& point, const Matrix4& base_transform, const std::map<JointPtr, VectorN>& q)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // setup point variables
  Real pointx = point[X], pointy = point[Y], pointz = point[Z];

  // determine the index of the outboard link
  RigidBodyPtr joutboard(joint->get_outboard_link());
  unsigned lidx = std::numeric_limits<unsigned>::max();
  for (unsigned i=1; i< get_links().size(); i++)
    if (joutboard == get_links()[i])
    {
      lidx = i;
      break;
    }
  if (lidx == std::numeric_limits<unsigned>::max())
    return MatrixN(0,0);

  // determine default joint positions
  Real q0l1L = JointPtr(get_links()[1]->get_inner_joint())->q[0];
  Real q0l2L = JointPtr(get_links()[2]->get_inner_joint())->q[0];

  // overwrite default joint positions with q
  JointPtr innerj1(get_links()[1]->get_inner_joint());
  if (q.find(innerj1) != q.end())
  {
    const VectorN& qdes = q.find(innerj1)->second;
    q0l1L = qdes[0];
  }

  JointPtr innerj2(get_links()[2]->get_inner_joint());
  if (q.find(innerj2) != q.end())
  {
    const VectorN& qdes = q.find(innerj2)->second;
    q0l2L = qdes[0];
  }

  // write the base transform
  Real rxxbaseL = get_links().front()->get_transform()(X,X);
  Real rxybaseL = get_links().front()->get_transform()(X,Y);
  Real rxzbaseL = get_links().front()->get_transform()(X,Z);
  Real txbaseL = get_links().front()->get_transform()(X,W);
  Real ryxbaseL = get_links().front()->get_transform()(Y,X);
  Real ryybaseL = get_links().front()->get_transform()(Y,Y);
  Real ryzbaseL = get_links().front()->get_transform()(Y,Z);
  Real tybaseL = get_links().front()->get_transform()(Y,W);
  Real rzxbaseL = get_links().front()->get_transform()(Z,X);
  Real rzybaseL = get_links().front()->get_transform()(Z,Y);
  Real rzzbaseL = get_links().front()->get_transform()(Z,Z);
  Real tzbaseL = get_links().front()->get_transform()(Z,W);

  // form the Jacobian depending on which joint was specified
  if (lidx == 1)
  {
    Real Compile_29549;
    Real Compile_29551;
    Real Compile_29552;
    Real Compile_29553;
    Real Compile_29554;
    Real Compile_29559;
    Compile_29549 = Cos(q0l1L);Compile_29551 = -Compile_29549;Compile_29552 = 1. + Compile_29551;Compile_29553 = Compile_29552;Compile_29554 = Compile_29553 + Compile_29549;Compile_29559 = Sin(q0l1L);
    return MatrixN::construct_variable(6,1,2.5*Compile_29549 - pointy*Compile_29554 - 2.5*Compile_29549*Compile_29554,pointx*Compile_29554 + 2.5*Compile_29559 - 2.5*Compile_29554*Compile_29559,0.,0.,0.,Compile_29554);
  }
  else  if (lidx == 2)
  {
    Real Compile_29625;
    Real Compile_29626;
    Real Compile_29627;
    Real Compile_29628;
    Real Compile_29629;
    Real Compile_29630;
    Real Compile_29631;
    Real Compile_29632;
    Real Compile_29633;
    Real Compile_29634;
    Real Compile_29636;
    Real Compile_29637;
    Real Compile_29638;
    Real Compile_29639;
    Real Compile_29640;
    Compile_29625 = Cos(q0l1L);Compile_29630 = Cos(q0l2L);Compile_29626 = -Compile_29625;Compile_29627 = 1. + Compile_29626;Compile_29628 = Compile_29627;Compile_29629 = Compile_29628 + Compile_29625;Compile_29631 = -Compile_29630;Compile_29632 = 1. + Compile_29631;Compile_29633 = Compile_29632;Compile_29634 = Compile_29633 + Compile_29630;Compile_29636 = Compile_29625*Compile_29630;Compile_29637 = Sin(q0l1L);Compile_29638 = Sin(q0l2L);Compile_29639 = -Compile_29637*Compile_29638;Compile_29640 = Compile_29636 + Compile_29639;
    return MatrixN::construct_variable(6,2,-pointy*Compile_29629*Compile_29634 + 2.5*Compile_29640 + Compile_29629*Compile_29634*(-5.*Compile_29625 - 2.5*Compile_29640),2.5*Compile_29625 - pointy*Compile_29629 - 2.5*Compile_29625*Compile_29629,pointx*Compile_29629*Compile_29634 + 2.5*(Compile_29630*Compile_29637 + Compile_29625*Compile_29638) + Compile_29629*Compile_29634*(-5.*Compile_29637 + 2.5*(-Compile_29630*Compile_29637 - Compile_29625*Compile_29638)),pointx*Compile_29629 + 2.5*Compile_29637 - 2.5*Compile_29629*Compile_29637,0.,0.,0.,0.,0.,0.,Compile_29629*Compile_29634,Compile_29629);
  }
  else    // handle the impossible case to make compiler happy
    return MatrixN(0,0);
}

MatrixN chain::calc_jacobian_column(JointPtr joint, const Vector3& point)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // setup point variables
  Real pointx = point[X], pointy = point[Y], pointz = point[Z];

  // determine the index of the outboard link
  RigidBodyPtr joutboard(joint->get_outboard_link());
  unsigned lidx = std::numeric_limits<unsigned>::max();
  for (unsigned i=1; i< get_links().size(); i++)
    if (joutboard == get_links()[i])
    {
      lidx = i;
      break;
    }
  if (lidx == std::numeric_limits<unsigned>::max())
    return MatrixN(0,0);



  // setup link transform values
  Real rxxbaseL = get_links()[0]->get_transform()(X,X);
  Real rxybaseL = get_links()[0]->get_transform()(X,Y);
  Real rxzbaseL = get_links()[0]->get_transform()(X,Z);
  Real txbaseL = get_links()[0]->get_transform()(X,W);
  Real ryxbaseL = get_links()[0]->get_transform()(Y,X);
  Real ryybaseL = get_links()[0]->get_transform()(Y,Y);
  Real ryzbaseL = get_links()[0]->get_transform()(Y,Z);
  Real tybaseL = get_links()[0]->get_transform()(Y,W);
  Real rzxbaseL = get_links()[0]->get_transform()(Z,X);
  Real rzybaseL = get_links()[0]->get_transform()(Z,Y);
  Real rzzbaseL = get_links()[0]->get_transform()(Z,Z);
  Real tzbaseL = get_links()[0]->get_transform()(Z,W);
  Real rxxl1L = get_links()[1]->get_transform()(X,X);
  Real rxyl1L = get_links()[1]->get_transform()(X,Y);
  Real rxzl1L = get_links()[1]->get_transform()(X,Z);
  Real txl1L = get_links()[1]->get_transform()(X,W);
  Real ryxl1L = get_links()[1]->get_transform()(Y,X);
  Real ryyl1L = get_links()[1]->get_transform()(Y,Y);
  Real ryzl1L = get_links()[1]->get_transform()(Y,Z);
  Real tyl1L = get_links()[1]->get_transform()(Y,W);
  Real rzxl1L = get_links()[1]->get_transform()(Z,X);
  Real rzyl1L = get_links()[1]->get_transform()(Z,Y);
  Real rzzl1L = get_links()[1]->get_transform()(Z,Z);
  Real tzl1L = get_links()[1]->get_transform()(Z,W);
  Real rxxl2L = get_links()[2]->get_transform()(X,X);
  Real rxyl2L = get_links()[2]->get_transform()(X,Y);
  Real rxzl2L = get_links()[2]->get_transform()(X,Z);
  Real txl2L = get_links()[2]->get_transform()(X,W);
  Real ryxl2L = get_links()[2]->get_transform()(Y,X);
  Real ryyl2L = get_links()[2]->get_transform()(Y,Y);
  Real ryzl2L = get_links()[2]->get_transform()(Y,Z);
  Real tyl2L = get_links()[2]->get_transform()(Y,W);
  Real rzxl2L = get_links()[2]->get_transform()(Z,X);
  Real rzyl2L = get_links()[2]->get_transform()(Z,Y);
  Real rzzl2L = get_links()[2]->get_transform()(Z,Z);
  Real tzl2L = get_links()[2]->get_transform()(Z,W);

  // form the Jacobian depending on which joint was specified
  if (lidx == 1)
  {
    
    return MatrixN::construct_variable(6,1,2.5*rxxl1L + pointz*ryzl1L - pointy*rzzl1L + (rzzl1L*tyl1L - ryzl1L*tzl1L),-pointz*rxzl1L + 2.5*ryxl1L + pointx*rzzl1L + (-rzzl1L*txl1L + rxzl1L*tzl1L),pointy*rxzl1L - pointx*ryzl1L + 2.5*rzxl1L + (ryzl1L*txl1L - rxzl1L*tyl1L),rxzl1L,ryzl1L,rzzl1L);
  }
  else  if (lidx == 2)
  {
    
    return MatrixN::construct_variable(6,2,2.5*rxxl2L + pointz*ryzl2L - pointy*rzzl2L + (rzzl2L*tyl2L - ryzl2L*tzl2L),2.5*rxxl1L + pointz*ryzl1L - pointy*rzzl1L + (rzzl1L*tyl1L - ryzl1L*tzl1L),-pointz*rxzl2L + 2.5*ryxl2L + pointx*rzzl2L + (-rzzl2L*txl2L + rxzl2L*tzl2L),-pointz*rxzl1L + 2.5*ryxl1L + pointx*rzzl1L + (-rzzl1L*txl1L + rxzl1L*tzl1L),pointy*rxzl2L - pointx*ryzl2L + 2.5*rzxl2L + (ryzl2L*txl2L - rxzl2L*tyl2L),pointy*rxzl1L - pointx*ryzl1L + 2.5*rzxl1L + (ryzl1L*txl1L - rxzl1L*tyl1L),rxzl2L,rxzl1L,ryzl2L,ryzl1L,rzzl2L,rzzl1L);
  }
  else    // handle the impossible case to make compiler happy
    return MatrixN(0,0);
}

