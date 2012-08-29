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

// write the factory so that we can load this class as a plugin
namespace Moby {
  extern "C" boost::shared_ptr<RCArticulatedBodySymbolic> factory() { return boost::shared_ptr<RCArticulatedBodySymbolic>(new robot); }}

void robot::update_link_transforms() const
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
  Real q0wheelcleftL = boost::shared_ptr<Joint>(get_links()[1]->get_inner_joint())->get_q()[0];
  Real q0wheelcrightL = boost::shared_ptr<Joint>(get_links()[2]->get_inner_joint())->get_q()[0];

  // declare work variables
  Real Internal_89;
  Real Internal_90;
  Real Internal_91;
  Real Internal_92;
  Real Internal_93;
  Real Internal_96;
  Real Internal_133;
  Real Internal_134;
  Real Internal_135;
  Real Internal_136;
  Real Internal_137;
  Real Internal_140;

  // compute work variables
  Internal_89 = Cos(q0wheelcleftL);Internal_96 = Sin(q0wheelcleftL);Internal_90 = -Internal_89;Internal_91 = 1. + Internal_90;Internal_92 = Internal_91;Internal_93 = Internal_92 + Internal_89;
  Internal_133 = Cos(q0wheelcrightL);Internal_140 = Sin(q0wheelcrightL);Internal_134 = -Internal_133;Internal_135 = 1. + Internal_134;Internal_136 = Internal_135;Internal_137 = Internal_136 + Internal_133;

  // compute transforms
  get_links()[1]->set_transform(MatrixN::construct_variable(4,4,rxxbaseL*Internal_93,rxybaseL*Internal_89 + rxzbaseL*Internal_96,rxzbaseL*Internal_89 - rxybaseL*Internal_96,-0.116*rxxbaseL + txbaseL,ryxbaseL*Internal_93,ryybaseL*Internal_89 + ryzbaseL*Internal_96,ryzbaseL*Internal_89 - ryybaseL*Internal_96,-0.116*ryxbaseL + tybaseL,rzxbaseL*Internal_93,rzybaseL*Internal_89 + rzzbaseL*Internal_96,rzzbaseL*Internal_89 - rzybaseL*Internal_96,-0.116*rzxbaseL + tzbaseL,0.,0.,0.,1.));
  get_links()[2]->set_transform(MatrixN::construct_variable(4,4,rxxbaseL*Internal_137,rxybaseL*Internal_133 + rxzbaseL*Internal_140,rxzbaseL*Internal_133 - rxybaseL*Internal_140,0.116*rxxbaseL + txbaseL,ryxbaseL*Internal_137,ryybaseL*Internal_133 + ryzbaseL*Internal_140,ryzbaseL*Internal_133 - ryybaseL*Internal_140,0.116*ryxbaseL + tybaseL,rzxbaseL*Internal_137,rzybaseL*Internal_133 + rzzbaseL*Internal_140,rzzbaseL*Internal_133 - rzybaseL*Internal_140,0.116*rzxbaseL + tzbaseL,0.,0.,0.,1.));
}

MatrixN robot::calc_jacobian_column(JointPtr joint, const Vector3& point, const Matrix4& base_transform, const std::map<JointPtr, VectorN>& q)
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
  Real q0wheelcleftL = JointPtr(get_links()[1]->get_inner_joint())->get_q()[0];
  Real q0wheelcrightL = JointPtr(get_links()[2]->get_inner_joint())->get_q()[0];

  // overwrite default joint positions with q
  JointPtr innerj1(get_links()[1]->get_inner_joint());
  if (q.find(innerj1) != q.end())
  {
    const VectorN& qdes = q.find(innerj1)->second;
    q0wheelcleftL = qdes[0];
  }

  JointPtr innerj2(get_links()[2]->get_inner_joint());
  if (q.find(innerj2) != q.end())
  {
    const VectorN& qdes = q.find(innerj2)->second;
    q0wheelcrightL = qdes[0];
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
    Real Internal_203;
    Real Internal_204;
    Real Internal_205;
    Real Internal_206;
    Real Internal_207;
    Real Internal_209;
    Real Internal_210;
    Real Internal_211;
    Real Internal_212;
    Real Internal_213;
    Real Internal_214;
    Real Internal_215;
    Real Internal_216;
    Real Internal_218;
    Real Internal_219;
    Real Internal_220;
    Real Internal_221;
    Real Internal_222;
    Real Internal_223;
    Real Internal_224;
    Real Internal_225;
    Real Internal_230;
    Real Internal_231;
    Real Internal_232;
    Real Internal_233;
    Real Internal_234;
    Real Internal_235;
    Real Internal_236;
    Real Internal_237;
    Internal_203 = Cos(q0wheelcleftL);Internal_204 = -Internal_203;Internal_205 = 1. + Internal_204;Internal_206 = Internal_205;Internal_207 = Internal_206 + Internal_203;Internal_209 = -0.116*rxxbaseL;Internal_210 = Internal_209 + txbaseL;Internal_211 = ryxbaseL*Internal_210*Internal_207;Internal_212 = 0.116*ryxbaseL;Internal_213 = -tybaseL;Internal_214 = Internal_212 + Internal_213;Internal_215 = rxxbaseL*Internal_214*Internal_207;Internal_216 = Internal_211 + Internal_215;Internal_230 = -0.116*ryxbaseL;Internal_231 = Internal_230 + tybaseL;Internal_232 = rzxbaseL*Internal_231*Internal_207;Internal_233 = 0.116*rzxbaseL;Internal_234 = -tzbaseL;Internal_235 = Internal_233 + Internal_234;Internal_236 = ryxbaseL*Internal_235*Internal_207;Internal_237 = Internal_232 + Internal_236;Internal_218 = 0.116*rxxbaseL;Internal_219 = -txbaseL;Internal_220 = Internal_218 + Internal_219;Internal_221 = rzxbaseL*Internal_220*Internal_207;Internal_222 = -0.116*rzxbaseL;Internal_223 = Internal_222 + tzbaseL;Internal_224 = rxxbaseL*Internal_223*Internal_207;Internal_225 = Internal_221 + Internal_224;
    return MatrixN::construct_variable(6,1,rxxbaseL*Internal_207 - pointy*Internal_216 + pointz*Internal_225,ryxbaseL*Internal_207 + pointx*Internal_216 - pointz*Internal_237,rzxbaseL*Internal_207 + pointy*Internal_237 - pointx*Internal_225,Internal_237,Internal_225,Internal_216);
  }
  else  if (lidx == 2)
  {
    Real Internal_281;
    Real Internal_282;
    Real Internal_283;
    Real Internal_284;
    Real Internal_285;
    Real Internal_287;
    Real Internal_288;
    Real Internal_289;
    Real Internal_290;
    Real Internal_291;
    Real Internal_292;
    Real Internal_293;
    Real Internal_294;
    Real Internal_296;
    Real Internal_297;
    Real Internal_298;
    Real Internal_299;
    Real Internal_300;
    Real Internal_301;
    Real Internal_302;
    Real Internal_303;
    Real Internal_308;
    Real Internal_309;
    Real Internal_310;
    Real Internal_311;
    Real Internal_312;
    Real Internal_313;
    Real Internal_314;
    Real Internal_315;
    Internal_281 = Cos(q0wheelcrightL);Internal_282 = -Internal_281;Internal_283 = 1. + Internal_282;Internal_284 = Internal_283;Internal_285 = Internal_284 + Internal_281;Internal_287 = 0.116*rxxbaseL;Internal_288 = Internal_287 + txbaseL;Internal_289 = ryxbaseL*Internal_288*Internal_285;Internal_290 = -0.116*ryxbaseL;Internal_291 = -tybaseL;Internal_292 = Internal_290 + Internal_291;Internal_293 = rxxbaseL*Internal_292*Internal_285;Internal_294 = Internal_289 + Internal_293;Internal_308 = 0.116*ryxbaseL;Internal_309 = Internal_308 + tybaseL;Internal_310 = rzxbaseL*Internal_309*Internal_285;Internal_311 = -0.116*rzxbaseL;Internal_312 = -tzbaseL;Internal_313 = Internal_311 + Internal_312;Internal_314 = ryxbaseL*Internal_313*Internal_285;Internal_315 = Internal_310 + Internal_314;Internal_296 = -0.116*rxxbaseL;Internal_297 = -txbaseL;Internal_298 = Internal_296 + Internal_297;Internal_299 = rzxbaseL*Internal_298*Internal_285;Internal_300 = 0.116*rzxbaseL;Internal_301 = Internal_300 + tzbaseL;Internal_302 = rxxbaseL*Internal_301*Internal_285;Internal_303 = Internal_299 + Internal_302;
    return MatrixN::construct_variable(6,1,rxxbaseL*Internal_285 - pointy*Internal_294 + pointz*Internal_303,ryxbaseL*Internal_285 + pointx*Internal_294 - pointz*Internal_315,rzxbaseL*Internal_285 + pointy*Internal_315 - pointx*Internal_303,Internal_315,Internal_303,Internal_294);
  }
  else    // handle the impossible case to make compiler happy
    return MatrixN(0,0);
}

MatrixN robot::calc_jacobian_column(JointPtr joint, const Vector3& point)
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
  Real rxxwheelcleftL = get_links()[1]->get_transform()(X,X);
  Real rxywheelcleftL = get_links()[1]->get_transform()(X,Y);
  Real rxzwheelcleftL = get_links()[1]->get_transform()(X,Z);
  Real txwheelcleftL = get_links()[1]->get_transform()(X,W);
  Real ryxwheelcleftL = get_links()[1]->get_transform()(Y,X);
  Real ryywheelcleftL = get_links()[1]->get_transform()(Y,Y);
  Real ryzwheelcleftL = get_links()[1]->get_transform()(Y,Z);
  Real tywheelcleftL = get_links()[1]->get_transform()(Y,W);
  Real rzxwheelcleftL = get_links()[1]->get_transform()(Z,X);
  Real rzywheelcleftL = get_links()[1]->get_transform()(Z,Y);
  Real rzzwheelcleftL = get_links()[1]->get_transform()(Z,Z);
  Real tzwheelcleftL = get_links()[1]->get_transform()(Z,W);
  Real rxxwheelcrightL = get_links()[2]->get_transform()(X,X);
  Real rxywheelcrightL = get_links()[2]->get_transform()(X,Y);
  Real rxzwheelcrightL = get_links()[2]->get_transform()(X,Z);
  Real txwheelcrightL = get_links()[2]->get_transform()(X,W);
  Real ryxwheelcrightL = get_links()[2]->get_transform()(Y,X);
  Real ryywheelcrightL = get_links()[2]->get_transform()(Y,Y);
  Real ryzwheelcrightL = get_links()[2]->get_transform()(Y,Z);
  Real tywheelcrightL = get_links()[2]->get_transform()(Y,W);
  Real rzxwheelcrightL = get_links()[2]->get_transform()(Z,X);
  Real rzywheelcrightL = get_links()[2]->get_transform()(Z,Y);
  Real rzzwheelcrightL = get_links()[2]->get_transform()(Z,Z);
  Real tzwheelcrightL = get_links()[2]->get_transform()(Z,W);

  // form the Jacobian depending on which joint was specified
  if (lidx == 1)
  {
    Real Internal_355;
    Real Internal_356;
    Real Internal_357;
    Real Internal_359;
    Real Internal_360;
    Real Internal_361;
    Real Internal_366;
    Real Internal_367;
    Real Internal_368;
    Internal_355 = ryxwheelcleftL*txwheelcleftL;Internal_356 = -rxxwheelcleftL*tywheelcleftL;Internal_357 = Internal_355 + Internal_356;Internal_359 = -rzxwheelcleftL*txwheelcleftL;Internal_360 = rxxwheelcleftL*tzwheelcleftL;Internal_361 = Internal_359 + Internal_360;Internal_366 = rzxwheelcleftL*tywheelcleftL;Internal_367 = -ryxwheelcleftL*tzwheelcleftL;Internal_368 = Internal_366 + Internal_367;
    return MatrixN::construct_variable(6,1,rxxwheelcleftL - pointy*Internal_357 + pointz*Internal_361,ryxwheelcleftL + pointx*Internal_357 - pointz*Internal_368,rzxwheelcleftL - pointx*Internal_361 + pointy*Internal_368,Internal_368,Internal_361,Internal_357);
  }
  else  if (lidx == 2)
  {
    Real Internal_407;
    Real Internal_408;
    Real Internal_409;
    Real Internal_411;
    Real Internal_412;
    Real Internal_413;
    Real Internal_418;
    Real Internal_419;
    Real Internal_420;
    Internal_407 = ryxwheelcrightL*txwheelcrightL;Internal_408 = -rxxwheelcrightL*tywheelcrightL;Internal_409 = Internal_407 + Internal_408;Internal_411 = -rzxwheelcrightL*txwheelcrightL;Internal_412 = rxxwheelcrightL*tzwheelcrightL;Internal_413 = Internal_411 + Internal_412;Internal_418 = rzxwheelcrightL*tywheelcrightL;Internal_419 = -ryxwheelcrightL*tzwheelcrightL;Internal_420 = Internal_418 + Internal_419;
    return MatrixN::construct_variable(6,1,rxxwheelcrightL - pointy*Internal_409 + pointz*Internal_413,ryxwheelcrightL + pointx*Internal_409 - pointz*Internal_420,rzxwheelcrightL - pointx*Internal_413 + pointy*Internal_420,Internal_420,Internal_413,Internal_409);
  }
  else    // handle the impossible case to make compiler happy
    return MatrixN(0,0);
}

