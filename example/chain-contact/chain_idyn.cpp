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

std::map<JointPtr, VectorN> chain::calc_inv_dyn(RCArticulatedBodyPtr body, const std::map<RigidBodyPtr, RCArticulatedBodyInvDynData>& inv_dyn_data)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // verify that inverse dynamics link data is the proper size
  assert(body->get_links().size() == inv_dyn_data.size());

  // declare necessary variables
  Real f0baseL, f1baseL, f2baseL;
  Real t0baseL, t1baseL, t2baseL;
  Real f0l1L, f1l1L, f2l1L;
  Real t0l1L, t1l1L, t2l1L;
  Real f0l2L, f1l2L, f2l2L;
  Real t0l2L, t1l2L, t2l2L;

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
  Real rxxl1L = body->get_links()[1]->get_transform()(X,X);
  Real rxyl1L = body->get_links()[1]->get_transform()(X,Y);
  Real rxzl1L = body->get_links()[1]->get_transform()(X,Z);
  Real txl1L = body->get_links()[1]->get_transform()(X,W);
  Real ryxl1L = body->get_links()[1]->get_transform()(Y,X);
  Real ryyl1L = body->get_links()[1]->get_transform()(Y,Y);
  Real ryzl1L = body->get_links()[1]->get_transform()(Y,Z);
  Real tyl1L = body->get_links()[1]->get_transform()(Y,W);
  Real rzxl1L = body->get_links()[1]->get_transform()(Z,X);
  Real rzyl1L = body->get_links()[1]->get_transform()(Z,Y);
  Real rzzl1L = get_links()[1]->get_transform()(Z,Z);
  Real tzl1L = body->get_links()[1]->get_transform()(Z,W);
  Real rxxl2L = body->get_links()[2]->get_transform()(X,X);
  Real rxyl2L = body->get_links()[2]->get_transform()(X,Y);
  Real rxzl2L = body->get_links()[2]->get_transform()(X,Z);
  Real txl2L = body->get_links()[2]->get_transform()(X,W);
  Real ryxl2L = body->get_links()[2]->get_transform()(Y,X);
  Real ryyl2L = body->get_links()[2]->get_transform()(Y,Y);
  Real ryzl2L = body->get_links()[2]->get_transform()(Y,Z);
  Real tyl2L = body->get_links()[2]->get_transform()(Y,W);
  Real rzxl2L = body->get_links()[2]->get_transform()(Z,X);
  Real rzyl2L = body->get_links()[2]->get_transform()(Z,Y);
  Real rzzl2L = get_links()[2]->get_transform()(Z,Z);
  Real tzl2L = body->get_links()[2]->get_transform()(Z,W);

  // set joint positions
  Real q0l1L = JointPtr(body->get_links()[1]->get_inner_joint())->q[0];
  Real q0l2L = JointPtr(body->get_links()[2]->get_inner_joint())->q[0];

  // set joint velocities
  Real qd0l1L = JointPtr(body->get_links()[1]->get_inner_joint())->qd[0];
  Real qd0l2L = JointPtr(body->get_links()[2]->get_inner_joint())->qd[0];

  // set link velocities
  Real v0baseL = body->get_links()[0]->get_spatial_velocity(eLink)[0];
  Real v1baseL = body->get_links()[0]->get_spatial_velocity(eLink)[1];
  Real v2baseL = body->get_links()[0]->get_spatial_velocity(eLink)[2];
  Real v3baseL = body->get_links()[0]->get_spatial_velocity(eLink)[3];
  Real v4baseL = body->get_links()[0]->get_spatial_velocity(eLink)[4];
  Real v5baseL = body->get_links()[0]->get_spatial_velocity(eLink)[5];
  Real v0l1L = body->get_links()[1]->get_spatial_velocity(eLink)[0];
  Real v1l1L = body->get_links()[1]->get_spatial_velocity(eLink)[1];
  Real v2l1L = body->get_links()[1]->get_spatial_velocity(eLink)[2];
  Real v3l1L = body->get_links()[1]->get_spatial_velocity(eLink)[3];
  Real v4l1L = body->get_links()[1]->get_spatial_velocity(eLink)[4];
  Real v5l1L = body->get_links()[1]->get_spatial_velocity(eLink)[5];
  Real v0l2L = body->get_links()[2]->get_spatial_velocity(eLink)[0];
  Real v1l2L = body->get_links()[2]->get_spatial_velocity(eLink)[1];
  Real v2l2L = body->get_links()[2]->get_spatial_velocity(eLink)[2];
  Real v3l2L = body->get_links()[2]->get_spatial_velocity(eLink)[3];
  Real v4l2L = body->get_links()[2]->get_spatial_velocity(eLink)[4];
  Real v5l2L = body->get_links()[2]->get_spatial_velocity(eLink)[5];

  // set external link forces / torques for base
  if (body->floating_base)
  {
    assert(inv_dyn_data.find(body->get_links().front()) != inv_dyn_data.end());
    const RCArticulatedBodyInvDynData& iddata = inv_dyn_data.find(body->get_links().front())->second;
    f0baseL = iddata.fext[X];
    f1baseL = iddata.fext[Y];
    f2baseL = iddata.fext[Z];
    t0baseL = iddata.text[X];
    t1baseL = iddata.text[Y];
    t2baseL = iddata.text[Z];
  }

  // set non-base external link forces / torques and desired joint acelerations
  assert(inv_dyn_data.find(body->get_links()[1]) != inv_dyn_data.end());
  const RCArticulatedBodyInvDynData& iddata1 = inv_dyn_data.find(body->get_links()[1])->second;
  f0l1L = iddata1.fext[X];
  f1l1L = iddata1.fext[Y];
  f2l1L = iddata1.fext[Z];
  t0l1L = iddata1.text[X];
  t1l1L = iddata1.text[Y];
  t2l1L = iddata1.text[Z];
  Real qdd0l1L = iddata1.qdd[0];
  assert(inv_dyn_data.find(body->get_links()[2]) != inv_dyn_data.end());
  const RCArticulatedBodyInvDynData& iddata2 = inv_dyn_data.find(body->get_links()[2])->second;
  f0l2L = iddata2.fext[X];
  f1l2L = iddata2.fext[Y];
  f2l2L = iddata2.fext[Z];
  t0l2L = iddata2.text[X];
  t1l2L = iddata2.text[Y];
  t2l2L = iddata2.text[Z];
  Real qdd0l2L = iddata2.qdd[0];

  // init the joint to link map
  std::map<JointPtr, VectorN> force_map;

  // declare work variables
  Real Compile_35194;
  Real Compile_35195;
  Real Compile_35196;
  Real Compile_35198;
  Real Compile_35199;
  Real Compile_35200;
  Real Compile_35201;
  Real Compile_35202;
  Real Compile_35203;
  Real Compile_35204;
  Real Compile_35205;
  Real Compile_35206;
  Real Compile_35207;
  Real Compile_35208;
  Real Compile_35209;
  Real Compile_35210;
  Real Compile_35211;
  Real Compile_35212;
  Real Compile_35213;
  Real Compile_35214;
  Real Compile_35215;
  Real Compile_35216;
  Real Compile_35217;
  Real Compile_35218;
  Real Compile_35219;
  Real Compile_35220;
  Real Compile_35221;
  Real Compile_35222;
  Real Compile_35223;
  Real Compile_35224;
  Real Compile_35225;
  Real Compile_35226;
  Real Compile_35227;
  Real Compile_35228;
  Real Compile_35229;
  Real Compile_35230;
  Real Compile_35231;
  Real Compile_35232;
  Real Compile_35233;
  Real Compile_35234;
  Real Compile_35235;
  Real Compile_35236;
  Real Compile_35237;
  Real Compile_35238;
  Real Compile_35239;
  Real Compile_35240;
  Real Compile_35241;
  Real Compile_35242;
  Real Compile_35243;
  Real Compile_35244;
  Real Compile_35245;
  Real Compile_35246;
  Real Compile_35247;
  Real Compile_35248;
  Real Compile_35249;
  Real Compile_35250;
  Real Compile_35251;
  Real Compile_35252;
  Real Compile_35253;
  Real Compile_35254;
  Real Compile_35255;
  Real Compile_35256;
  Real Compile_35257;
  Real Compile_35258;
  Real Compile_35259;
  Real Compile_35260;
  Real Compile_35261;
  Real Compile_35262;
  Real Compile_35263;
  Real Compile_35264;
  Real Compile_35265;
  Real Compile_35266;
  Real Compile_35267;
  Real Compile_35268;
  Real Compile_35269;
  Real Compile_35270;
  Real Compile_35271;
  Real Compile_35272;
  Real Compile_35273;
  Real Compile_35274;
  Real Compile_35277;
  Real Compile_35278;
  Real Compile_35279;
  Real Compile_35280;
  Real Compile_35281;
  Real Compile_35282;
  Real Compile_35283;
  Real Compile_35284;
  Real Compile_35285;
  Real Compile_35286;
  Real Compile_35287;
  Real Compile_35288;
  Real Compile_35289;
  Real Compile_35290;
  Real Compile_35291;
  Real Compile_35292;
  Real Compile_35293;
  Real Compile_35294;
  Real Compile_35295;
  Real Compile_35296;
  Real Compile_35297;
  Real Compile_35298;
  Real Compile_35299;
  Real Compile_35300;
  Real Compile_35301;
  Real Compile_35302;
  Real Compile_35303;
  Real Compile_35304;
  Real Compile_35305;
  Real Compile_35306;
  Real Compile_35307;
  Real Compile_35308;
  Real Compile_35309;
  Real Compile_35310;
  Real Compile_35311;
  Real Compile_35313;
  Real Compile_35314;
  Real Compile_35315;
  Real Compile_35316;
  Real Compile_35317;
  Real Compile_35318;
  Real Compile_35319;
  Real Compile_35320;
  Real Compile_35321;
  Real Compile_35322;
  Real Compile_35323;
  Real Compile_35324;
  Real Compile_35325;
  Real Compile_35326;
  Real Compile_35327;
  Real Compile_35328;
  Real Compile_35329;
  Real Compile_35330;
  Real Compile_35331;
  Real Compile_35332;
  Real Compile_35333;
  Real Compile_35334;
  Real Compile_35335;
  Real Compile_35336;
  Real Compile_35337;
  Real Compile_35338;
  Real Compile_35339;
  Real Compile_35340;
  Real Compile_35341;
  Real Compile_35342;
  Real Compile_35343;
  Real Compile_35344;
  Real Compile_35345;
  Real Compile_35390;
  Real Compile_35391;
  Real Compile_35392;
  Real Compile_35393;
  Real Compile_35394;
  Real Compile_35395;
  Real Compile_35396;
  Real Compile_35397;
  Real Compile_35398;
  Real Compile_35399;
  Real Compile_35401;
  Real Compile_35402;
  Real Compile_35403;
  Real Compile_35404;
  Real Compile_35449;
  Real Compile_35450;
  Real Compile_35451;
  Real Compile_35452;
  Real Compile_35454;
  Real Compile_35455;
  Real Compile_35456;
  Real Compile_35457;
  Real Compile_35459;
  Real Compile_35460;
  Real Compile_35461;
  Real Compile_35462;
  Real Compile_35474;
  Real Compile_35475;
  Real Compile_35476;
  Real Compile_35477;
  Real Compile_35478;
  Real Compile_35479;
  Real Compile_35480;
  Real Compile_35481;
  Real Compile_35482;
  Real Compile_35483;
  Real Compile_35489;
  Real Compile_35490;
  Real Compile_35491;
  Real Compile_35492;

  // compute work variables
  Compile_35210 = -txl2L;Compile_35211 = txl1L + Compile_35210;Compile_35213 = -tyl2L;Compile_35214 = tyl1L + Compile_35213;Compile_35216 = -tzl2L;Compile_35217 = tzl1L + Compile_35216;Compile_35212 = rxxl2L*Compile_35211;Compile_35215 = ryxl2L*Compile_35214;Compile_35218 = rzxl2L*Compile_35217;Compile_35219 = Compile_35212 + Compile_35215 + Compile_35218;Compile_35225 = -rxyl2L*Compile_35211;Compile_35226 = -ryyl2L*Compile_35214;Compile_35227 = -rzyl2L*Compile_35217;Compile_35228 = Compile_35225 + Compile_35226 + Compile_35227;Compile_35198 = rxxl1L*rxzl2L;Compile_35199 = ryxl1L*ryzl2L;Compile_35200 = rzxl1L*rzzl2L;Compile_35201 = Compile_35198 + Compile_35199 + Compile_35200;Compile_35194 = 2.5*qdd0l1L;Compile_35195 = qd0l1L*v4l1L;Compile_35196 = Compile_35194 + Compile_35195;Compile_35249 = rxxl1L*rxyl2L;Compile_35250 = ryxl1L*ryyl2L;Compile_35251 = rzxl1L*rzyl2L;Compile_35252 = Compile_35249 + Compile_35250 + Compile_35251;Compile_35244 = rxzl1L*rxzl2L;Compile_35245 = ryzl1L*ryzl2L;Compile_35246 = rzzl1L*rzzl2L;Compile_35247 = Compile_35244 + Compile_35245 + Compile_35246;Compile_35221 = rxxl2L*rxzl1L;Compile_35222 = ryxl2L*ryzl1L;Compile_35223 = rzxl2L*rzzl1L;Compile_35224 = Compile_35221 + Compile_35222 + Compile_35223;Compile_35262 = rxyl1L*rxzl2L;Compile_35263 = ryyl1L*ryzl2L;Compile_35264 = rzyl1L*rzzl2L;Compile_35265 = Compile_35262 + Compile_35263 + Compile_35264;Compile_35281 = -rxxl2L*Compile_35211;Compile_35282 = -ryxl2L*Compile_35214;Compile_35283 = -rzxl2L*Compile_35217;Compile_35284 = Compile_35281 + Compile_35282 + Compile_35283;Compile_35237 = rxxl2L*rxyl1L;Compile_35238 = ryxl2L*ryyl1L;Compile_35239 = rzxl2L*rzyl1L;Compile_35240 = Compile_35237 + Compile_35238 + Compile_35239;Compile_35286 = rxzl2L*Compile_35211;Compile_35287 = ryzl2L*Compile_35214;Compile_35288 = rzzl2L*Compile_35217;Compile_35289 = Compile_35286 + Compile_35287 + Compile_35288;Compile_35206 = rxyl2L*rxzl1L;Compile_35207 = ryyl2L*ryzl1L;Compile_35208 = rzyl2L*rzzl1L;Compile_35209 = Compile_35206 + Compile_35207 + Compile_35208;Compile_35254 = rxxl1L*rxxl2L;Compile_35255 = ryxl1L*ryxl2L;Compile_35256 = rzxl1L*rzxl2L;Compile_35257 = Compile_35254 + Compile_35255 + Compile_35256;Compile_35232 = rxyl1L*rxyl2L;Compile_35233 = ryyl1L*ryyl2L;Compile_35234 = rzyl1L*rzyl2L;Compile_35235 = Compile_35232 + Compile_35233 + Compile_35234;Compile_35266 = 2.5*v2l1L;Compile_35267 = -v3l1L;Compile_35268 = Compile_35266 + Compile_35267;Compile_35318 = rxyl2L*Compile_35211;Compile_35319 = ryyl2L*Compile_35214;Compile_35320 = rzyl2L*Compile_35217;Compile_35321 = Compile_35318 + Compile_35319 + Compile_35320;Compile_35323 = -rxzl2L*Compile_35211;Compile_35324 = -ryzl2L*Compile_35214;Compile_35325 = -rzzl2L*Compile_35217;Compile_35326 = Compile_35323 + Compile_35324 + Compile_35325;Compile_35390 = -txl1L;Compile_35391 = Compile_35390 + txl2L;Compile_35393 = -tyl1L;Compile_35394 = Compile_35393 + tyl2L;Compile_35396 = -tzl1L;Compile_35397 = Compile_35396 + tzl2L;Compile_35202 = -f0l2L*rxzl2L;Compile_35203 = -f1l2L*ryzl2L;Compile_35204 = -f2l2L*rzzl2L;Compile_35205 = -15.708*v1l2L*v3l2L;Compile_35220 = Compile_35209*Compile_35219;Compile_35229 = Compile_35224*Compile_35228;Compile_35230 = Compile_35220 + Compile_35229;Compile_35231 = qdd0l1L*Compile_35230;Compile_35236 = Compile_35235*Compile_35219;Compile_35241 = Compile_35240*Compile_35228;Compile_35242 = Compile_35236 + Compile_35241;Compile_35243 = -qd0l1L*Compile_35242*v0l1L;Compile_35248 = -2.5*qd0l1L*Compile_35247*v1l1L;Compile_35253 = Compile_35252*Compile_35219;Compile_35258 = Compile_35257*Compile_35228;Compile_35259 = Compile_35253 + Compile_35258;Compile_35260 = qd0l1L*Compile_35259*v1l1L;Compile_35261 = -2.5*qd0l2L*v1l2L;Compile_35269 = qd0l1L*Compile_35265*Compile_35268;Compile_35270 = Compile_35201*Compile_35196;Compile_35271 = Compile_35231 + Compile_35243 + Compile_35248 + Compile_35260 + Compile_35261 + Compile_35269 + Compile_35270;Compile_35272 = 15.708*Compile_35271;Compile_35273 = 15.708*v0l2L*v4l2L;Compile_35274 = Compile_35202 + Compile_35203 + Compile_35204 + Compile_35205 + Compile_35272 + Compile_35273;Compile_35392 = rxxl1L*Compile_35391;Compile_35395 = ryxl1L*Compile_35394;Compile_35398 = rzxl1L*Compile_35397;Compile_35399 = Compile_35392 + Compile_35395 + Compile_35398;Compile_35401 = -rxyl1L*Compile_35391;Compile_35402 = -ryyl1L*Compile_35394;Compile_35403 = -rzyl1L*Compile_35397;Compile_35404 = Compile_35401 + Compile_35402 + Compile_35403;Compile_35277 = -f0l2L*rxyl2L;Compile_35278 = -f1l2L*ryyl2L;Compile_35279 = -f2l2L*rzyl2L;Compile_35280 = 15.708*v2l2L*v3l2L;Compile_35285 = Compile_35247*Compile_35284;Compile_35290 = Compile_35224*Compile_35289;Compile_35291 = Compile_35285 + Compile_35290;Compile_35292 = qdd0l1L*Compile_35291;Compile_35293 = Compile_35265*Compile_35284;Compile_35294 = Compile_35240*Compile_35289;Compile_35295 = Compile_35293 + Compile_35294;Compile_35296 = -qd0l1L*Compile_35295*v0l1L;Compile_35297 = -2.5*qd0l1L*Compile_35209*v1l1L;Compile_35298 = Compile_35201*Compile_35284;Compile_35299 = Compile_35257*Compile_35289;Compile_35300 = Compile_35298 + Compile_35299;Compile_35301 = qd0l1L*Compile_35300*v1l1L;Compile_35302 = qd0l1L*Compile_35235*Compile_35268;Compile_35303 = 2.5*v2l2L;Compile_35304 = -v3l2L;Compile_35305 = Compile_35303 + Compile_35304;Compile_35306 = qd0l2L*Compile_35305;Compile_35307 = Compile_35252*Compile_35196;Compile_35308 = Compile_35292 + Compile_35296 + Compile_35297 + Compile_35301 + Compile_35302 + Compile_35306 + Compile_35307;Compile_35309 = 15.708*Compile_35308;Compile_35310 = -15.708*v0l2L*v5l2L;Compile_35311 = Compile_35277 + Compile_35278 + Compile_35279 + Compile_35280 + Compile_35309 + Compile_35310;Compile_35313 = -f0l2L*rxxl2L;Compile_35314 = -f1l2L*ryxl2L;Compile_35315 = -f2l2L*rzxl2L;Compile_35316 = -15.708*v2l2L*v4l2L;Compile_35317 = 2.5*qdd0l2L;Compile_35322 = Compile_35247*Compile_35321;Compile_35327 = Compile_35209*Compile_35326;Compile_35328 = Compile_35322 + Compile_35327;Compile_35329 = qdd0l1L*Compile_35328;Compile_35330 = Compile_35265*Compile_35321;Compile_35331 = Compile_35235*Compile_35326;Compile_35332 = Compile_35330 + Compile_35331;Compile_35333 = -qd0l1L*Compile_35332*v0l1L;Compile_35334 = -2.5*qd0l1L*Compile_35224*v1l1L;Compile_35335 = Compile_35201*Compile_35321;Compile_35336 = Compile_35252*Compile_35326;Compile_35337 = Compile_35335 + Compile_35336;Compile_35338 = qd0l1L*Compile_35337*v1l1L;Compile_35339 = qd0l1L*Compile_35240*Compile_35268;Compile_35340 = Compile_35257*Compile_35196;Compile_35341 = qd0l2L*v4l2L;Compile_35342 = Compile_35317 + Compile_35329 + Compile_35333 + Compile_35334 + Compile_35338 + Compile_35339 + Compile_35340 + Compile_35341;Compile_35343 = 15.708*Compile_35342;Compile_35344 = 15.708*v1l2L*v5l2L;Compile_35345 = Compile_35313 + Compile_35314 + Compile_35315 + Compile_35316 + Compile_35343 + Compile_35344;
  Compile_35449 = rxzl1L*rxzl2L;Compile_35450 = ryzl1L*ryzl2L;Compile_35451 = rzzl1L*rzzl2L;Compile_35452 = Compile_35449 + Compile_35450 + Compile_35451;Compile_35474 = -txl2L;Compile_35475 = txl1L + Compile_35474;Compile_35477 = -tyl2L;Compile_35478 = tyl1L + Compile_35477;Compile_35480 = -tzl2L;Compile_35481 = tzl1L + Compile_35480;Compile_35454 = rxyl1L*rxzl2L;Compile_35455 = ryyl1L*ryzl2L;Compile_35456 = rzyl1L*rzzl2L;Compile_35457 = Compile_35454 + Compile_35455 + Compile_35456;Compile_35476 = rxyl2L*Compile_35475;Compile_35479 = ryyl2L*Compile_35478;Compile_35482 = rzyl2L*Compile_35481;Compile_35483 = Compile_35476 + Compile_35479 + Compile_35482;Compile_35489 = -rxzl2L*Compile_35475;Compile_35490 = -ryzl2L*Compile_35478;Compile_35491 = -rzzl2L*Compile_35481;Compile_35492 = Compile_35489 + Compile_35490 + Compile_35491;Compile_35459 = rxxl1L*rxzl2L;Compile_35460 = ryxl1L*ryzl2L;Compile_35461 = rzxl1L*rzzl2L;Compile_35462 = Compile_35459 + Compile_35460 + Compile_35461;

  // write forces for inverse dynamics
  force_map[JointPtr(body->get_links()[1]->get_inner_joint())] = 2.5*(-f0l1L*rxxl1L - f1l1L*ryxl1L - f2l1L*rzxl1L - 15.708*v2l1L*v4l1L + 15.708*Compile_35196 + Compile_35201*Compile_35274 + 15.708*v1l1L*v5l1L + Compile_35252*Compile_35311 + Compile_35257*Compile_35345) + (36.651914*qdd0l1L - rxzl1L*t0l1L - ryzl1L*t1l1L - rzzl1L*t2l1L - 28.797931999999996*v0l1L*v1l1L + Compile_35247*(-rxzl2L*t0l2L - ryzl2L*t1l2L - rzzl2L*t2l2L + 36.651914*(qdd0l2L + qdd0l1L*Compile_35247 - qd0l1L*Compile_35265*v0l1L + qd0l1L*Compile_35201*v1l1L) - 28.797931999999996*v0l2L*v1l2L) + Compile_35209*(-rxyl2L*t0l2L - ryyl2L*t1l2L - rzyl2L*t2l2L + 7.853982*(qdd0l1L*Compile_35209 - qd0l1L*Compile_35235*v0l1L - qd0l2L*v0l2L + qd0l1L*Compile_35252*v1l1L) + 0.*v0l2L*v2l2L) + Compile_35224*(-rxxl2L*t0l2L - ryxl2L*t1l2L - rzxl2L*t2l2L + 36.651914*(qdd0l1L*Compile_35224 - qd0l1L*Compile_35240*v0l1L + qd0l1L*Compile_35257*v1l1L + qd0l2L*v1l2L) + 28.797931999999996*v1l2L*v2l2L) + (Compile_35265*Compile_35399 + Compile_35201*Compile_35404)*Compile_35274 + (Compile_35235*Compile_35399 + Compile_35252*Compile_35404)*Compile_35311 + (Compile_35240*Compile_35399 + Compile_35257*Compile_35404)*Compile_35345);
  force_map[JointPtr(body->get_links()[2]->get_inner_joint())] = (-rxzl2L*t0l2L - ryzl2L*t1l2L - rzzl2L*t2l2L + 36.651914*(qdd0l2L + qdd0l1L*Compile_35452 - qd0l1L*Compile_35457*v0l1L + qd0l1L*Compile_35462*v1l1L) - 28.797931999999996*v0l2L*v1l2L) + 2.5*(-f0l2L*rxxl2L - f1l2L*ryxl2L - f2l2L*rzxl2L - 15.708*v2l2L*v4l2L + 15.708*(2.5*qdd0l2L + qdd0l1L*(Compile_35452*Compile_35483 + (rxyl2L*rxzl1L + ryyl2L*ryzl1L + rzyl2L*rzzl1L)*Compile_35492) - qd0l1L*(Compile_35457*Compile_35483 + (rxyl1L*rxyl2L + ryyl1L*ryyl2L + rzyl1L*rzyl2L)*Compile_35492)*v0l1L - 2.5*qd0l1L*(rxxl2L*rxzl1L + ryxl2L*ryzl1L + rzxl2L*rzzl1L)*v1l1L + qd0l1L*(Compile_35462*Compile_35483 + (rxxl1L*rxyl2L + ryxl1L*ryyl2L + rzxl1L*rzyl2L)*Compile_35492)*v1l1L + qd0l1L*(rxxl2L*rxyl1L + ryxl2L*ryyl1L + rzxl2L*rzyl1L)*(2.5*v2l1L - v3l1L) + (rxxl1L*rxxl2L + ryxl1L*ryxl2L + rzxl1L*rzxl2L)*(2.5*qdd0l1L + qd0l1L*v4l1L) + qd0l2L*v4l2L) + 15.708*v1l2L*v5l2L);
}

