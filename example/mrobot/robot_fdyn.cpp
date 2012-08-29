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

void robot::update_link_velocities()
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // setup base velocity
  Real v0baseL = get_links().front()->get_spatial_velocity(eLink)[0];
  Real v1baseL = get_links().front()->get_spatial_velocity(eLink)[1];
  Real v2baseL = get_links().front()->get_spatial_velocity(eLink)[2];
  Real v3baseL = get_links().front()->get_spatial_velocity(eLink)[3];
  Real v4baseL = get_links().front()->get_spatial_velocity(eLink)[4];
  Real v5baseL = get_links().front()->get_spatial_velocity(eLink)[5];

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

  // setup qd
  Real qd0wheelcleftL = boost::shared_ptr<Joint>(get_links()[1]->get_inner_joint())->get_qd()[0];
  Real qd0wheelcrightL = boost::shared_ptr<Joint>(get_links()[2]->get_inner_joint())->get_qd()[0];

  // declare work variables
  Real Internal_1013;
  Real Internal_1014;
  Real Internal_1015;
  Real Internal_1016;
  Real Internal_1018;
  Real Internal_1019;
  Real Internal_1020;
  Real Internal_1021;
  Real Internal_1037;
  Real Internal_1038;
  Real Internal_1039;
  Real Internal_1040;
  Real Internal_1042;
  Real Internal_1043;
  Real Internal_1044;
  Real Internal_1045;
  Real Internal_935;
  Real Internal_936;
  Real Internal_937;
  Real Internal_938;
  Real Internal_940;
  Real Internal_941;
  Real Internal_942;
  Real Internal_943;
  Real Internal_945;
  Real Internal_946;
  Real Internal_947;
  Real Internal_948;
  Real Internal_951;
  Real Internal_952;
  Real Internal_953;
  Real Internal_954;
  Real Internal_956;
  Real Internal_957;
  Real Internal_958;
  Real Internal_959;
  Real Internal_961;
  Real Internal_962;
  Real Internal_963;
  Real Internal_964;
  Real Internal_967;
  Real Internal_968;
  Real Internal_969;
  Real Internal_970;
  Real Internal_972;
  Real Internal_973;
  Real Internal_974;
  Real Internal_975;
  Real Internal_977;
  Real Internal_978;
  Real Internal_979;
  Real Internal_980;
  Real Internal_983;
  Real Internal_984;
  Real Internal_985;
  Real Internal_986;
  Real Internal_987;
  Real Internal_988;
  Real Internal_989;
  Real Internal_990;
  Real Internal_991;
  Real Internal_992;
  Real Internal_994;
  Real Internal_995;
  Real Internal_996;
  Real Internal_997;
  Real Internal_1070;
  Real Internal_1071;
  Real Internal_1072;
  Real Internal_1073;
  Real Internal_1075;
  Real Internal_1076;
  Real Internal_1077;
  Real Internal_1078;
  Real Internal_1080;
  Real Internal_1081;
  Real Internal_1082;
  Real Internal_1083;
  Real Internal_1086;
  Real Internal_1087;
  Real Internal_1088;
  Real Internal_1089;
  Real Internal_1091;
  Real Internal_1092;
  Real Internal_1093;
  Real Internal_1094;
  Real Internal_1096;
  Real Internal_1097;
  Real Internal_1098;
  Real Internal_1099;
  Real Internal_1102;
  Real Internal_1103;
  Real Internal_1104;
  Real Internal_1105;
  Real Internal_1107;
  Real Internal_1108;
  Real Internal_1109;
  Real Internal_1110;
  Real Internal_1112;
  Real Internal_1113;
  Real Internal_1114;
  Real Internal_1115;
  Real Internal_1118;
  Real Internal_1119;
  Real Internal_1120;
  Real Internal_1121;
  Real Internal_1122;
  Real Internal_1123;
  Real Internal_1124;
  Real Internal_1125;
  Real Internal_1126;
  Real Internal_1127;
  Real Internal_1129;
  Real Internal_1130;
  Real Internal_1131;
  Real Internal_1132;
  Real Internal_1148;
  Real Internal_1149;
  Real Internal_1150;
  Real Internal_1151;
  Real Internal_1153;
  Real Internal_1154;
  Real Internal_1155;
  Real Internal_1156;
  Real Internal_1172;
  Real Internal_1173;
  Real Internal_1174;
  Real Internal_1175;
  Real Internal_1177;
  Real Internal_1178;
  Real Internal_1179;
  Real Internal_1180;

  // compute work variables
  Internal_967 = rxxbaseL*rxzwheelcleftL;Internal_968 = ryxbaseL*ryzwheelcleftL;Internal_969 = rzxbaseL*rzzwheelcleftL;Internal_970 = Internal_967 + Internal_968 + Internal_969;Internal_951 = rxxbaseL*rxywheelcleftL;Internal_952 = ryxbaseL*ryywheelcleftL;Internal_953 = rzxbaseL*rzywheelcleftL;Internal_954 = Internal_951 + Internal_952 + Internal_953;Internal_983 = -txwheelcleftL;Internal_984 = txbaseL + Internal_983;Internal_986 = -tywheelcleftL;Internal_987 = tybaseL + Internal_986;Internal_989 = -tzwheelcleftL;Internal_990 = tzbaseL + Internal_989;Internal_972 = rxybaseL*rxzwheelcleftL;Internal_973 = ryybaseL*ryzwheelcleftL;Internal_974 = rzybaseL*rzzwheelcleftL;Internal_975 = Internal_972 + Internal_973 + Internal_974;Internal_985 = rxywheelcleftL*Internal_984;Internal_988 = ryywheelcleftL*Internal_987;Internal_991 = rzywheelcleftL*Internal_990;Internal_992 = Internal_985 + Internal_988 + Internal_991;Internal_956 = rxybaseL*rxywheelcleftL;Internal_957 = ryybaseL*ryywheelcleftL;Internal_958 = rzybaseL*rzywheelcleftL;Internal_959 = Internal_956 + Internal_957 + Internal_958;Internal_994 = -rxzwheelcleftL*Internal_984;Internal_995 = -ryzwheelcleftL*Internal_987;Internal_996 = -rzzwheelcleftL*Internal_990;Internal_997 = Internal_994 + Internal_995 + Internal_996;Internal_977 = rxzbaseL*rxzwheelcleftL;Internal_978 = ryzbaseL*ryzwheelcleftL;Internal_979 = rzzbaseL*rzzwheelcleftL;Internal_980 = Internal_977 + Internal_978 + Internal_979;Internal_961 = rxywheelcleftL*rxzbaseL;Internal_962 = ryywheelcleftL*ryzbaseL;Internal_963 = rzywheelcleftL*rzzbaseL;Internal_964 = Internal_961 + Internal_962 + Internal_963;Internal_935 = rxxbaseL*rxxwheelcleftL;Internal_936 = ryxbaseL*ryxwheelcleftL;Internal_937 = rzxbaseL*rzxwheelcleftL;Internal_938 = Internal_935 + Internal_936 + Internal_937;Internal_940 = rxxwheelcleftL*rxybaseL;Internal_941 = ryxwheelcleftL*ryybaseL;Internal_942 = rzxwheelcleftL*rzybaseL;Internal_943 = Internal_940 + Internal_941 + Internal_942;Internal_945 = rxxwheelcleftL*rxzbaseL;Internal_946 = ryxwheelcleftL*ryzbaseL;Internal_947 = rzxwheelcleftL*rzzbaseL;Internal_948 = Internal_945 + Internal_946 + Internal_947;Internal_1013 = -rxxwheelcleftL*Internal_984;Internal_1014 = -ryxwheelcleftL*Internal_987;Internal_1015 = -rzxwheelcleftL*Internal_990;Internal_1016 = Internal_1013 + Internal_1014 + Internal_1015;Internal_1018 = rxzwheelcleftL*Internal_984;Internal_1019 = ryzwheelcleftL*Internal_987;Internal_1020 = rzzwheelcleftL*Internal_990;Internal_1021 = Internal_1018 + Internal_1019 + Internal_1020;Internal_1037 = rxxwheelcleftL*Internal_984;Internal_1038 = ryxwheelcleftL*Internal_987;Internal_1039 = rzxwheelcleftL*Internal_990;Internal_1040 = Internal_1037 + Internal_1038 + Internal_1039;Internal_1042 = -rxywheelcleftL*Internal_984;Internal_1043 = -ryywheelcleftL*Internal_987;Internal_1044 = -rzywheelcleftL*Internal_990;Internal_1045 = Internal_1042 + Internal_1043 + Internal_1044;
  Internal_1102 = rxxbaseL*rxzwheelcrightL;Internal_1103 = ryxbaseL*ryzwheelcrightL;Internal_1104 = rzxbaseL*rzzwheelcrightL;Internal_1105 = Internal_1102 + Internal_1103 + Internal_1104;Internal_1086 = rxxbaseL*rxywheelcrightL;Internal_1087 = ryxbaseL*ryywheelcrightL;Internal_1088 = rzxbaseL*rzywheelcrightL;Internal_1089 = Internal_1086 + Internal_1087 + Internal_1088;Internal_1118 = -txwheelcrightL;Internal_1119 = txbaseL + Internal_1118;Internal_1121 = -tywheelcrightL;Internal_1122 = tybaseL + Internal_1121;Internal_1124 = -tzwheelcrightL;Internal_1125 = tzbaseL + Internal_1124;Internal_1107 = rxybaseL*rxzwheelcrightL;Internal_1108 = ryybaseL*ryzwheelcrightL;Internal_1109 = rzybaseL*rzzwheelcrightL;Internal_1110 = Internal_1107 + Internal_1108 + Internal_1109;Internal_1120 = rxywheelcrightL*Internal_1119;Internal_1123 = ryywheelcrightL*Internal_1122;Internal_1126 = rzywheelcrightL*Internal_1125;Internal_1127 = Internal_1120 + Internal_1123 + Internal_1126;Internal_1091 = rxybaseL*rxywheelcrightL;Internal_1092 = ryybaseL*ryywheelcrightL;Internal_1093 = rzybaseL*rzywheelcrightL;Internal_1094 = Internal_1091 + Internal_1092 + Internal_1093;Internal_1129 = -rxzwheelcrightL*Internal_1119;Internal_1130 = -ryzwheelcrightL*Internal_1122;Internal_1131 = -rzzwheelcrightL*Internal_1125;Internal_1132 = Internal_1129 + Internal_1130 + Internal_1131;Internal_1112 = rxzbaseL*rxzwheelcrightL;Internal_1113 = ryzbaseL*ryzwheelcrightL;Internal_1114 = rzzbaseL*rzzwheelcrightL;Internal_1115 = Internal_1112 + Internal_1113 + Internal_1114;Internal_1096 = rxywheelcrightL*rxzbaseL;Internal_1097 = ryywheelcrightL*ryzbaseL;Internal_1098 = rzywheelcrightL*rzzbaseL;Internal_1099 = Internal_1096 + Internal_1097 + Internal_1098;Internal_1070 = rxxbaseL*rxxwheelcrightL;Internal_1071 = ryxbaseL*ryxwheelcrightL;Internal_1072 = rzxbaseL*rzxwheelcrightL;Internal_1073 = Internal_1070 + Internal_1071 + Internal_1072;Internal_1075 = rxxwheelcrightL*rxybaseL;Internal_1076 = ryxwheelcrightL*ryybaseL;Internal_1077 = rzxwheelcrightL*rzybaseL;Internal_1078 = Internal_1075 + Internal_1076 + Internal_1077;Internal_1080 = rxxwheelcrightL*rxzbaseL;Internal_1081 = ryxwheelcrightL*ryzbaseL;Internal_1082 = rzxwheelcrightL*rzzbaseL;Internal_1083 = Internal_1080 + Internal_1081 + Internal_1082;Internal_1148 = -rxxwheelcrightL*Internal_1119;Internal_1149 = -ryxwheelcrightL*Internal_1122;Internal_1150 = -rzxwheelcrightL*Internal_1125;Internal_1151 = Internal_1148 + Internal_1149 + Internal_1150;Internal_1153 = rxzwheelcrightL*Internal_1119;Internal_1154 = ryzwheelcrightL*Internal_1122;Internal_1155 = rzzwheelcrightL*Internal_1125;Internal_1156 = Internal_1153 + Internal_1154 + Internal_1155;Internal_1172 = rxxwheelcrightL*Internal_1119;Internal_1173 = ryxwheelcrightL*Internal_1122;Internal_1174 = rzxwheelcrightL*Internal_1125;Internal_1175 = Internal_1172 + Internal_1173 + Internal_1174;Internal_1177 = -rxywheelcrightL*Internal_1119;Internal_1178 = -ryywheelcrightL*Internal_1122;Internal_1179 = -rzywheelcrightL*Internal_1125;Internal_1180 = Internal_1177 + Internal_1178 + Internal_1179;

  // set velocities
  get_links()[1]->set_spatial_velocity(VectorN::construct_variable(6,qd0wheelcleftL + Internal_938*v0baseL + Internal_943*v1baseL + Internal_948*v2baseL,Internal_954*v0baseL + Internal_959*v1baseL + Internal_964*v2baseL,Internal_970*v0baseL + Internal_975*v1baseL + Internal_980*v2baseL,(Internal_970*Internal_992 + Internal_954*Internal_997)*v0baseL + (Internal_975*Internal_992 + Internal_959*Internal_997)*v1baseL + (Internal_980*Internal_992 + Internal_964*Internal_997)*v2baseL + Internal_938*v3baseL + Internal_943*v4baseL + Internal_948*v5baseL,(Internal_970*Internal_1016 + Internal_938*Internal_1021)*v0baseL + (Internal_975*Internal_1016 + Internal_943*Internal_1021)*v1baseL + (Internal_980*Internal_1016 + Internal_948*Internal_1021)*v2baseL + Internal_954*v3baseL + Internal_959*v4baseL + Internal_964*v5baseL,(Internal_954*Internal_1040 + Internal_938*Internal_1045)*v0baseL + (Internal_959*Internal_1040 + Internal_943*Internal_1045)*v1baseL + (Internal_964*Internal_1040 + Internal_948*Internal_1045)*v2baseL + Internal_970*v3baseL + Internal_975*v4baseL + Internal_980*v5baseL), eLink);
  get_links()[2]->set_spatial_velocity(VectorN::construct_variable(6,qd0wheelcrightL + Internal_1073*v0baseL + Internal_1078*v1baseL + Internal_1083*v2baseL,Internal_1089*v0baseL + Internal_1094*v1baseL + Internal_1099*v2baseL,Internal_1105*v0baseL + Internal_1110*v1baseL + Internal_1115*v2baseL,(Internal_1105*Internal_1127 + Internal_1089*Internal_1132)*v0baseL + (Internal_1110*Internal_1127 + Internal_1094*Internal_1132)*v1baseL + (Internal_1115*Internal_1127 + Internal_1099*Internal_1132)*v2baseL + Internal_1073*v3baseL + Internal_1078*v4baseL + Internal_1083*v5baseL,(Internal_1105*Internal_1151 + Internal_1073*Internal_1156)*v0baseL + (Internal_1110*Internal_1151 + Internal_1078*Internal_1156)*v1baseL + (Internal_1115*Internal_1151 + Internal_1083*Internal_1156)*v2baseL + Internal_1089*v3baseL + Internal_1094*v4baseL + Internal_1099*v5baseL,(Internal_1089*Internal_1175 + Internal_1073*Internal_1180)*v0baseL + (Internal_1094*Internal_1175 + Internal_1078*Internal_1180)*v1baseL + (Internal_1099*Internal_1175 + Internal_1083*Internal_1180)*v2baseL + Internal_1105*v3baseL + Internal_1110*v4baseL + Internal_1115*v5baseL), eLink);
}

void robot::calc_fwd_dyn()
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;


  // set link transforms
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

  // set joint velocities
  Real qd0wheelcleftL = JointPtr(get_links()[1]->get_inner_joint())->get_qd()[0];
  Real qd0wheelcrightL = JointPtr(get_links()[2]->get_inner_joint())->get_qd()[0];

  // set actuator forces
  Real qf0wheelcleftL = JointPtr(get_links()[1]->get_inner_joint())->get_scaled_force()[0];
  Real qf0wheelcrightL = JointPtr(get_links()[2]->get_inner_joint())->get_scaled_force()[0];

  // set link velocities
  Real v0baseL = get_links()[0]->get_spatial_velocity(eLink)[0];
  Real v1baseL = get_links()[0]->get_spatial_velocity(eLink)[1];
  Real v2baseL = get_links()[0]->get_spatial_velocity(eLink)[2];
  Real v3baseL = get_links()[0]->get_spatial_velocity(eLink)[3];
  Real v4baseL = get_links()[0]->get_spatial_velocity(eLink)[4];
  Real v5baseL = get_links()[0]->get_spatial_velocity(eLink)[5];
  Real v0wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[0];
  Real v1wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[1];
  Real v2wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[2];
  Real v3wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[3];
  Real v4wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[4];
  Real v5wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[5];
  Real v0wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[0];
  Real v1wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[1];
  Real v2wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[2];
  Real v3wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[3];
  Real v4wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[4];
  Real v5wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[5];

  // set external link forces / torques
  Real f0baseL = get_links()[0]->sum_forces()[0];
  Real t0baseL = get_links()[0]->sum_forces()[0];
  Real f1baseL = get_links()[0]->sum_forces()[1];
  Real t1baseL = get_links()[0]->sum_forces()[1];
  Real f2baseL = get_links()[0]->sum_forces()[2];
  Real t2baseL = get_links()[0]->sum_forces()[2];
  Real f0wheelcleftL = get_links()[1]->sum_forces()[0];
  Real t0wheelcleftL = get_links()[1]->sum_forces()[0];
  Real f1wheelcleftL = get_links()[1]->sum_forces()[1];
  Real t1wheelcleftL = get_links()[1]->sum_forces()[1];
  Real f2wheelcleftL = get_links()[1]->sum_forces()[2];
  Real t2wheelcleftL = get_links()[1]->sum_forces()[2];
  Real f0wheelcrightL = get_links()[2]->sum_forces()[0];
  Real t0wheelcrightL = get_links()[2]->sum_forces()[0];
  Real f1wheelcrightL = get_links()[2]->sum_forces()[1];
  Real t1wheelcrightL = get_links()[2]->sum_forces()[1];
  Real f2wheelcrightL = get_links()[2]->sum_forces()[2];
  Real t2wheelcrightL = get_links()[2]->sum_forces()[2];

  // Compute forward dynamics

  get_links().front()->set_spatial_accel(VectorN::construct_variable(6,0.,0.,0.,0.,0.,0.), eLink);

  get_links()[1]->set_spatial_accel(VectorN::construct_variable(6,247.27992087042534*(qf0wheelcleftL - 2.e-6*qd0wheelcleftL*v1wheelcleftL - (-rxxwheelcleftL*t0wheelcleftL - rxywheelcleftL*t1wheelcleftL - rxzwheelcleftL*t2wheelcleftL - 2.e-6*v0wheelcleftL*v1wheelcleftL + v1wheelcleftL*(-0.000048*v1wheelcleftL - 0.002058*v2wheelcleftL) + (0.002093*v1wheelcleftL + 0.000048*v2wheelcleftL)*v2wheelcleftL)),qd0wheelcleftL*v2wheelcleftL,-qd0wheelcleftL*v1wheelcleftL,0.,qd0wheelcleftL*v5wheelcleftL,-qd0wheelcleftL*v4wheelcleftL), eLink);

  JointPtr(get_links()[1]->get_inner_joint())->set_qdd(VectorN::construct_variable(1,247.27992087042534*(qf0wheelcleftL - 2.e-6*qd0wheelcleftL*v1wheelcleftL - (-rxxwheelcleftL*t0wheelcleftL - rxywheelcleftL*t1wheelcleftL - rxzwheelcleftL*t2wheelcleftL - 2.e-6*v0wheelcleftL*v1wheelcleftL + v1wheelcleftL*(-0.000048*v1wheelcleftL - 0.002058*v2wheelcleftL) + (0.002093*v1wheelcleftL + 0.000048*v2wheelcleftL)*v2wheelcleftL))));

  get_links()[2]->set_spatial_accel(VectorN::construct_variable(6,247.27992087042534*(qf0wheelcrightL - 2.e-6*qd0wheelcrightL*v1wheelcrightL - (-rxxwheelcrightL*t0wheelcrightL - rxywheelcrightL*t1wheelcrightL - rxzwheelcrightL*t2wheelcrightL - 2.e-6*v0wheelcrightL*v1wheelcrightL + v1wheelcrightL*(-0.000048*v1wheelcrightL - 0.002058*v2wheelcrightL) + (0.002093*v1wheelcrightL + 0.000048*v2wheelcrightL)*v2wheelcrightL)),qd0wheelcrightL*v2wheelcrightL,-qd0wheelcrightL*v1wheelcrightL,0.,qd0wheelcrightL*v5wheelcrightL,-qd0wheelcrightL*v4wheelcrightL), eLink);

  JointPtr(get_links()[2]->get_inner_joint())->set_qdd(VectorN::construct_variable(1,247.27992087042534*(qf0wheelcrightL - 2.e-6*qd0wheelcrightL*v1wheelcrightL - (-rxxwheelcrightL*t0wheelcrightL - rxywheelcrightL*t1wheelcrightL - rxzwheelcrightL*t2wheelcrightL - 2.e-6*v0wheelcrightL*v1wheelcrightL + v1wheelcrightL*(-0.000048*v1wheelcrightL - 0.002058*v2wheelcrightL) + (0.002093*v1wheelcrightL + 0.000048*v2wheelcrightL)*v2wheelcrightL))));

  // setup state derivative
  set_state_deriv();
}

void robot::apply_impulse(const Vector3& j, const Vector3& k, const Vector3& point, RigidBodyPtr link)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  Real pointx = point[X], pointy = point[Y], pointz = point[Z];
  Real jx = j[X], jy = j[Y], jz = j[Z];
  Real kx = k[X], ky = k[Y], kz = k[Z];

  // set link transforms
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

  // set joint velocities
  Real qd0wheelcleftL = JointPtr(get_links()[1]->get_inner_joint())->get_qd()[0];
  Real qd0wheelcrightL = JointPtr(get_links()[2]->get_inner_joint())->get_qd()[0];

  // set link velocities
  Real v0baseL = get_links()[0]->get_spatial_velocity(eLink)[0];
  Real v1baseL = get_links()[0]->get_spatial_velocity(eLink)[1];
  Real v2baseL = get_links()[0]->get_spatial_velocity(eLink)[2];
  Real v3baseL = get_links()[0]->get_spatial_velocity(eLink)[3];
  Real v4baseL = get_links()[0]->get_spatial_velocity(eLink)[4];
  Real v5baseL = get_links()[0]->get_spatial_velocity(eLink)[5];
  Real v0wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[0];
  Real v1wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[1];
  Real v2wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[2];
  Real v3wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[3];
  Real v4wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[4];
  Real v5wheelcleftL = get_links()[1]->get_spatial_velocity(eLink)[5];
  Real v0wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[0];
  Real v1wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[1];
  Real v2wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[2];
  Real v3wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[3];
  Real v4wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[4];
  Real v5wheelcrightL = get_links()[2]->get_spatial_velocity(eLink)[5];

  // compute result of applying impulses
  if (link == get_links()[0])
  {
    // init work variables for updating link velocity

    // compute work variables for updating link velocity
    

    // update link velocity
    get_links()[0]->set_spatial_velocity(VectorN::construct_variable(6,v0baseL,v1baseL,v2baseL,v3baseL,v4baseL,v5baseL), eLink);

    // init work variables for updating link velocity

    // compute work variables for updating link velocity
    

    // update link velocity
    get_links()[1]->set_spatial_velocity(VectorN::construct_variable(6,v0wheelcleftL,v1wheelcleftL,v2wheelcleftL,v3wheelcleftL,v4wheelcleftL,v5wheelcleftL), eLink);

    // init work variables for updating joint velocity

    // compute work variables for updating joint velocity
    

    // update joint velocity
    JointPtr(get_links()[1]->get_inner_joint())->set_qd(VectorN::construct_variable(1,qd0wheelcleftL));

    // init work variables for updating link velocity

    // compute work variables for updating link velocity
    

    // update link velocity
    get_links()[2]->set_spatial_velocity(VectorN::construct_variable(6,v0wheelcrightL,v1wheelcrightL,v2wheelcrightL,v3wheelcrightL,v4wheelcrightL,v5wheelcrightL), eLink);

    // init work variables for updating joint velocity

    // compute work variables for updating joint velocity
    

    // update joint velocity
    JointPtr(get_links()[2]->get_inner_joint())->set_qd(VectorN::construct_variable(1,qd0wheelcrightL));
  }
  else  if (link == get_links()[1])
  {
    // init work variables for updating link velocity

    // compute work variables for updating link velocity
    

    // update link velocity
    get_links()[0]->set_spatial_velocity(VectorN::construct_variable(6,v0baseL,v1baseL,v2baseL,v3baseL,v4baseL,v5baseL), eLink);

    // init work variables for updating link velocity
    Real Internal_1548;
    Real Internal_1549;
    Real Internal_1550;
    Real Internal_1551;
    Real Internal_1552;
    Real Internal_1553;
    Real Internal_1554;
    Real Internal_1555;
    Real Internal_1556;
    Real Internal_1557;
    Real Internal_1559;
    Real Internal_1560;
    Real Internal_1561;
    Real Internal_1562;

    // compute work variables for updating link velocity
    Internal_1548 = -txwheelcleftL;Internal_1549 = pointx + Internal_1548;Internal_1551 = -tywheelcleftL;Internal_1552 = pointy + Internal_1551;Internal_1554 = -tzwheelcleftL;Internal_1555 = pointz + Internal_1554;Internal_1550 = rxywheelcleftL*Internal_1549;Internal_1553 = ryywheelcleftL*Internal_1552;Internal_1556 = rzywheelcleftL*Internal_1555;Internal_1557 = Internal_1550 + Internal_1553 + Internal_1556;Internal_1559 = -rxzwheelcleftL*Internal_1549;Internal_1560 = -ryzwheelcleftL*Internal_1552;Internal_1561 = -rzzwheelcleftL*Internal_1555;Internal_1562 = Internal_1559 + Internal_1560 + Internal_1561;

    // update link velocity
    get_links()[1]->set_spatial_velocity(VectorN::construct_variable(6,-247.27992087042534*(-jx*rxxwheelcleftL - jy*ryxwheelcleftL - jz*rzxwheelcleftL - kx*(rxzwheelcleftL*Internal_1557 + rxywheelcleftL*Internal_1562) - ky*(ryzwheelcleftL*Internal_1557 + ryywheelcleftL*Internal_1562) - kz*(rzzwheelcleftL*Internal_1557 + rzywheelcleftL*Internal_1562)) + v0wheelcleftL,v1wheelcleftL,v2wheelcleftL,v3wheelcleftL,v4wheelcleftL,v5wheelcleftL), eLink);

    // init work variables for updating joint velocity
    Real Internal_1605;
    Real Internal_1606;
    Real Internal_1607;
    Real Internal_1608;
    Real Internal_1609;
    Real Internal_1610;
    Real Internal_1611;
    Real Internal_1612;
    Real Internal_1613;
    Real Internal_1614;
    Real Internal_1616;
    Real Internal_1617;
    Real Internal_1618;
    Real Internal_1619;

    // compute work variables for updating joint velocity
    Internal_1605 = -txwheelcleftL;Internal_1606 = pointx + Internal_1605;Internal_1608 = -tywheelcleftL;Internal_1609 = pointy + Internal_1608;Internal_1611 = -tzwheelcleftL;Internal_1612 = pointz + Internal_1611;Internal_1607 = rxywheelcleftL*Internal_1606;Internal_1610 = ryywheelcleftL*Internal_1609;Internal_1613 = rzywheelcleftL*Internal_1612;Internal_1614 = Internal_1607 + Internal_1610 + Internal_1613;Internal_1616 = -rxzwheelcleftL*Internal_1606;Internal_1617 = -ryzwheelcleftL*Internal_1609;Internal_1618 = -rzzwheelcleftL*Internal_1612;Internal_1619 = Internal_1616 + Internal_1617 + Internal_1618;

    // update joint velocity
    JointPtr(get_links()[1]->get_inner_joint())->set_qd(VectorN::construct_variable(1,qd0wheelcleftL - 247.27992087042534*(-jx*rxxwheelcleftL - jy*ryxwheelcleftL - jz*rzxwheelcleftL - kx*(rxzwheelcleftL*Internal_1614 + rxywheelcleftL*Internal_1619) - ky*(ryzwheelcleftL*Internal_1614 + ryywheelcleftL*Internal_1619) - kz*(rzzwheelcleftL*Internal_1614 + rzywheelcleftL*Internal_1619))));

    // init work variables for updating link velocity

    // compute work variables for updating link velocity
    

    // update link velocity
    get_links()[2]->set_spatial_velocity(VectorN::construct_variable(6,v0wheelcrightL,v1wheelcrightL,v2wheelcrightL,v3wheelcrightL,v4wheelcrightL,v5wheelcrightL), eLink);

    // init work variables for updating joint velocity

    // compute work variables for updating joint velocity
    

    // update joint velocity
    JointPtr(get_links()[2]->get_inner_joint())->set_qd(VectorN::construct_variable(1,qd0wheelcrightL));
  }
  else  if (link == get_links()[2])
  {
    // init work variables for updating link velocity

    // compute work variables for updating link velocity
    

    // update link velocity
    get_links()[0]->set_spatial_velocity(VectorN::construct_variable(6,v0baseL,v1baseL,v2baseL,v3baseL,v4baseL,v5baseL), eLink);

    // init work variables for updating link velocity

    // compute work variables for updating link velocity
    

    // update link velocity
    get_links()[1]->set_spatial_velocity(VectorN::construct_variable(6,v0wheelcleftL,v1wheelcleftL,v2wheelcleftL,v3wheelcleftL,v4wheelcleftL,v5wheelcleftL), eLink);

    // init work variables for updating joint velocity

    // compute work variables for updating joint velocity
    

    // update joint velocity
    JointPtr(get_links()[1]->get_inner_joint())->set_qd(VectorN::construct_variable(1,qd0wheelcleftL));

    // init work variables for updating link velocity
    Real Internal_1770;
    Real Internal_1771;
    Real Internal_1772;
    Real Internal_1773;
    Real Internal_1774;
    Real Internal_1775;
    Real Internal_1776;
    Real Internal_1777;
    Real Internal_1778;
    Real Internal_1779;
    Real Internal_1781;
    Real Internal_1782;
    Real Internal_1783;
    Real Internal_1784;

    // compute work variables for updating link velocity
    Internal_1770 = -txwheelcrightL;Internal_1771 = pointx + Internal_1770;Internal_1773 = -tywheelcrightL;Internal_1774 = pointy + Internal_1773;Internal_1776 = -tzwheelcrightL;Internal_1777 = pointz + Internal_1776;Internal_1772 = rxywheelcrightL*Internal_1771;Internal_1775 = ryywheelcrightL*Internal_1774;Internal_1778 = rzywheelcrightL*Internal_1777;Internal_1779 = Internal_1772 + Internal_1775 + Internal_1778;Internal_1781 = -rxzwheelcrightL*Internal_1771;Internal_1782 = -ryzwheelcrightL*Internal_1774;Internal_1783 = -rzzwheelcrightL*Internal_1777;Internal_1784 = Internal_1781 + Internal_1782 + Internal_1783;

    // update link velocity
    get_links()[2]->set_spatial_velocity(VectorN::construct_variable(6,-247.27992087042534*(-jx*rxxwheelcrightL - jy*ryxwheelcrightL - jz*rzxwheelcrightL - kx*(rxzwheelcrightL*Internal_1779 + rxywheelcrightL*Internal_1784) - ky*(ryzwheelcrightL*Internal_1779 + ryywheelcrightL*Internal_1784) - kz*(rzzwheelcrightL*Internal_1779 + rzywheelcrightL*Internal_1784)) + v0wheelcrightL,v1wheelcrightL,v2wheelcrightL,v3wheelcrightL,v4wheelcrightL,v5wheelcrightL), eLink);

    // init work variables for updating joint velocity
    Real Internal_1827;
    Real Internal_1828;
    Real Internal_1829;
    Real Internal_1830;
    Real Internal_1831;
    Real Internal_1832;
    Real Internal_1833;
    Real Internal_1834;
    Real Internal_1835;
    Real Internal_1836;
    Real Internal_1838;
    Real Internal_1839;
    Real Internal_1840;
    Real Internal_1841;

    // compute work variables for updating joint velocity
    Internal_1827 = -txwheelcrightL;Internal_1828 = pointx + Internal_1827;Internal_1830 = -tywheelcrightL;Internal_1831 = pointy + Internal_1830;Internal_1833 = -tzwheelcrightL;Internal_1834 = pointz + Internal_1833;Internal_1829 = rxywheelcrightL*Internal_1828;Internal_1832 = ryywheelcrightL*Internal_1831;Internal_1835 = rzywheelcrightL*Internal_1834;Internal_1836 = Internal_1829 + Internal_1832 + Internal_1835;Internal_1838 = -rxzwheelcrightL*Internal_1828;Internal_1839 = -ryzwheelcrightL*Internal_1831;Internal_1840 = -rzzwheelcrightL*Internal_1834;Internal_1841 = Internal_1838 + Internal_1839 + Internal_1840;

    // update joint velocity
    JointPtr(get_links()[2]->get_inner_joint())->set_qd(VectorN::construct_variable(1,qd0wheelcrightL - 247.27992087042534*(-jx*rxxwheelcrightL - jy*ryxwheelcrightL - jz*rzxwheelcrightL - kx*(rxzwheelcrightL*Internal_1836 + rxywheelcrightL*Internal_1841) - ky*(ryzwheelcrightL*Internal_1836 + ryywheelcrightL*Internal_1841) - kz*(rzzwheelcrightL*Internal_1836 + rzywheelcrightL*Internal_1841))));
  }
  else
    assert(false);  // should never get here...
}
