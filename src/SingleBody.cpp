/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <Moby/ArticulatedBody.h>
#include <Moby/SingleBody.h>

using namespace Ravelin;
using namespace Moby;

/// Gets the "super" body
DynamicBodyPtr SingleBody::get_super_body() const
{
/*
  boost::shared_ptr<const ArticulatedBody> abc = get_articulated_body();
  if (abc)
  {
    ArticulatedBodyPtr ab = boost::const_pointer_cast<ArticulatedBody>(abc);
    return boost::dynamic_pointer_cast<DynamicBody>(ab);
  }
*/
  ArticulatedBodyPtr ab = get_articulated_body();
  if (ab)
    return boost::dynamic_pointer_cast<DynamicBody>(ab);
  else
  {
    BasePtr base = boost::const_pointer_cast<Base>(shared_from_this());
    return boost::dynamic_pointer_cast<DynamicBody>(base);
  }
}

double SingleBody::calc_point_vel(const Point3d& point, const Ravelin::Vector3d& dir)
{
  Vector3d pv = calc_point_vel(point);
  return Pose3d::transform_vector(dir.pose, pv).dot(dir);
}

