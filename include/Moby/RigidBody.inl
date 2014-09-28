/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

template <class OutputIterator>
OutputIterator RigidBody::get_parent_links(OutputIterator begin) const
{
  BOOST_FOREACH(JointPtr j, _inner_joints)
    *begin++ = get_parent_link(j);

  return begin;
}

template <class OutputIterator>
OutputIterator RigidBody::get_child_links(OutputIterator begin) const
{
  BOOST_FOREACH(JointPtr j, _outer_joints)
    *begin++ = get_child_link(j);

  return begin;
}

/// Sets generalized coordinates using a templated vector
template <class V>
void RigidBody::get_generalized_coordinates_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gc) 
{
  const unsigned N_SPATIAL = 6, N_EULER = 7;

  // special case: disabled body
  if (!_enabled)
  {
    gc.resize(0);
    return;
  }

  // resize vector
  switch (gctype)
  {
    case DynamicBody::eEuler:   gc.resize(N_EULER); break;
    case DynamicBody::eSpatial: gc.resize(N_SPATIAL); break;
  }

  // convert current pose to global frame
  Ravelin::Pose3d P = *_F;
  P.update_relative_pose(GLOBAL);

  // get linear components
  gc[0] = P.x[0];
  gc[1] = P.x[1];
  gc[2] = P.x[2];

  // get angular components 
  if (gctype == DynamicBody::eSpatial)
    P.q.to_rpy(gc[3], gc[4], gc[5]);
  else
  {
    // return the generalized position using Euler parameters
    assert(gctype == DynamicBody::eEuler);
    gc[3] = P.q.x;
    gc[4] = P.q.y;
    gc[5] = P.q.z;
    gc[6] = P.q.w;
  }
}

/// Sets the generalized coordinates of this rigid body (does not call articulated body)
template <class V>
void RigidBody::set_generalized_coordinates_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gc)
{
  // special case: disabled body
  if (!_enabled)
    return;

  // do easiest case first 
  if (gctype == DynamicBody::eSpatial)
  {
    // note: generalized coordinates in eSpatial are always set with regard to 
    // the global frame
    Ravelin::Origin3d x(gc[0], gc[1], gc[2]);
    Ravelin::Quatd q = Ravelin::Quatd::rpy(gc[3], gc[4], gc[5]); 

    // convert the pose to the correct relative frame
    Ravelin::Pose3d P(q, x);
    P.update_relative_pose(_F->rpose);

    // set the transform
    set_pose(P);
  }
  else
  {
    assert(gctype == DynamicBody::eEuler);

    // get the position
    Ravelin::Origin3d x(gc[0], gc[1], gc[2]);

    // get the unit quaternion
    Ravelin::Quatd q;
    q.x = gc[3];
    q.y = gc[4];
    q.z = gc[5];
    q.w = gc[6];

    // normalize the unit quaternion, just in case
    q.normalize();

    // coordinates are in the global frame; must convert them to the
    // relative frame
    Ravelin::Pose3d P(q, x);
    P.update_relative_pose(_F->rpose); 

    // set the transform
    set_pose(P);
  }
}

/// Sets the generalized velocity of this rigid body (does not call articulated body version)
template <class V>
void RigidBody::set_generalized_velocity_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gv)
{
  // special case: disabled body
  if (!_enabled)
    return;

  // get the velocity
  Ravelin::SVelocityd xd;
  xd.pose = _F2;
 
  // set the linear velocity first
  xd.set_linear(Ravelin::Vector3d(gv[0], gv[1], gv[2]));
 
  // simplest case: spatial coordinates
  if (gctype == DynamicBody::eSpatial)
    xd.set_angular(Ravelin::Vector3d(gv[3], gv[4], gv[5]));
  else
  {
    assert(gctype == DynamicBody::eEuler);

    // get the quaternion derivatives
    Ravelin::Quatd qd;
    qd.x = gv[3] * 2.0;
    qd.y = gv[4] * 2.0;
    qd.z = gv[5] * 2.0;
    qd.w = gv[6] * 2.0;

    // setup the pose
    Ravelin::Pose3d F = *_F;
    F.update_relative_pose(GLOBAL);

    // setup the angular component
    xd.set_angular(F.q.G_mult(qd.x, qd.y, qd.z, qd.w));
  }

  // set the velocity
  set_velocity(xd);
}

/// Gets the generalized velocity of this rigid body (does not call articulated body version)
template <class V>
void RigidBody::get_generalized_velocity_generic(DynamicBody::GeneralizedCoordinateType gctype, V& gv) 
{
  const unsigned N_SPATIAL = 6, N_EULER = 7;

  // special case: disabled body
  if (!_enabled)
  {
    gv.resize(0);
    return;
  }

  // resize the generalized velocity vector
  switch (gctype)
  {
    case DynamicBody::eEuler:   gv.resize(N_EULER); break;
    case DynamicBody::eSpatial: gv.resize(N_SPATIAL); break;
  }

  // get/set linear components of velocity
  Ravelin::Vector3d lv = _xdcom.get_linear();
  gv[0] = lv[0];
  gv[1] = lv[1];
  gv[2] = lv[2];

  // determine the proper generalized coordinate type
  if (gctype == DynamicBody::eSpatial)
  {
    // get/set angular components of velocity
    Ravelin::Vector3d av = _xdcom.get_angular();
    gv[3] = av[0];
    gv[4] = av[1];
    gv[5] = av[2];
  }
  else
  {
    assert(gctype == DynamicBody::eEuler);

    // going to need Euler coordinate derivatives
    Ravelin::Pose3d F = *_F;
    F.update_relative_pose(GLOBAL);
    Ravelin::Quatd qd = F.q.G_transpose_mult(_xdcom.get_angular()) * 0.5;

    // setup the angular components 
    gv[3] = qd.x;
    gv[4] = qd.y;
    gv[5] = qd.z;
    gv[6] = qd.w; 
  }
}

/// Gets the generalized acceleration of this body (does not call articulated body version)
template <class V>
void RigidBody::get_generalized_acceleration_generic(V& ga)
{
  const unsigned N_SPATIAL = 6;

  // special case: body is disabled
  if (!_enabled)
  {
    ga.resize(0);
    return;
  } 

  // setup the linear components
  ga.resize(N_SPATIAL);

  // get linear and angular components
  Ravelin::Vector3d la = _xddcom.get_linear();
  Ravelin::Vector3d aa = _xddcom.get_angular();

  // set linear components
  ga[0] = la[0];
  ga[1] = la[1];
  ga[2] = la[2];
  ga[3] = aa[0];
  ga[4] = aa[1];
  ga[5] = aa[2];
}

