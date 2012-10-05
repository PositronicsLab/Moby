/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <boost/foreach.hpp>
#include <queue>
#include <iostream>
#include <iomanip>
#include <limits>
#include <Moby/NullPointerException.h>
#include <Moby/AAngle.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/ArticulatedBody.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/XMLTree.h>
#include <Moby/Joint.h>
#include <Moby/LinAlg.h>
#include <Moby/Log.h>
#include <Moby/RigidBody.h>

using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using boost::static_pointer_cast;
using std::cerr;
using std::endl;
using std::map;
using std::list;
using std::queue;

/// Default constructor
/**
 * Constructs a rigid body with zero mass, zero inertia tensor, and center
 * of mass at [0,0,0] with position at [0,0,0], identity orientation, and zero 
 * linear and angular velocity.  Body is enabled by default.
 */
RigidBody::RigidBody()
{
  const unsigned SPATIAL_DIM = 6;

  _mass = 0;
  _inv_mass = std::numeric_limits<Real>::max();
  _J.set_zero();
  _invJ = IDENTITY_3x3 * std::numeric_limits<Real>::max();
  _x = ZEROS_3;
  _xd = ZEROS_3;
  _omega = ZEROS_3;
  _F = IDENTITY_4x4;
  _q = Quat(&_F);
  _forces = ZEROS_3;
  _torques = ZEROS_3;
  _enabled = true;
  _link_idx = std::numeric_limits<unsigned>::max();
  coulomb_coeff = VectorN::zero(SPATIAL_DIM);
  viscous_coeff = VectorN::zero(SPATIAL_DIM);
}

/// Integrates the body forward in time
void RigidBody::integrate(Real t, Real h, shared_ptr<Integrator<VectorN> > integrator)
{
  // don't attempt to integrate disabled bodies
  if (!is_enabled())
    return;

  // if this body is a link in an articulated body, don't attempt to integrate;
  // integration must occur at the articulated body
  if (!_abody.expired())
    return;

  // call parent method if still here
  DynamicBody::integrate(t, h, integrator);

  FILE_LOG(LOG_DYNAMICS) << "RigidBody::integrate()" << endl;
  FILE_LOG(LOG_DYNAMICS) << "  new transform: " << endl << _F;
  FILE_LOG(LOG_DYNAMICS) << "  new velocity: " << _xd << " (lin)  " << _omega << " (ang)" << endl;
}

/// Computes the forward dynamics for this body
void RigidBody::calc_fwd_dyn(Real dt)
{
  // if the body is free, just compute linear and angular acceleration via
  // Newton's and Euler's laws
  if (_abody.expired())
  {
    // determine the linear acceleration
    _xdd = _inv_mass * _forces;

    // determine the angular acceleration
    Matrix3 R(&_q);
    Matrix3 iJ = Matrix3::transpose(R) * _invJ * R;
    _alpha = iJ * (_torques - calc_inertial_forces());
  }
  else
  {
    // otherwise, need to call forward dynamics on the articulated body
    ArticulatedBodyPtr abody(_abody);
  
    // calculate forward dynamics on it
    abody->calc_fwd_dyn(dt);
  }
}

/// Gets the dynamic body
/**
 * If this is a rigid body, returns this; if this is a link in an articulated
 * body, returns the articulated body.
 */
shared_ptr<const DynamicBody> RigidBody::get_dynamic_body() const
{
  if (_abody.expired())
    return shared_ptr<const DynamicBody>(get_this());
  else
  {
    ArticulatedBodyPtr abody(get_articulated_body());
    return shared_ptr<const DynamicBody>(abody);
  }
}

/// Gets the dynamic body
/**
 * If this is a rigid body, returns this; if this is a link in an articulated
 * body, returns the articulated body.
 */
shared_ptr<DynamicBody> RigidBody::get_dynamic_body() 
{
  if (_abody.expired())
    return shared_ptr<DynamicBody>(get_this());
  else
  {
    ArticulatedBodyPtr abody(get_articulated_body());
    return shared_ptr<DynamicBody>(abody);
  }
}

/// Sets the body to enabled / disabled.
/**
 * If the body is disabled, the linear and angular velocity are set to zero,
 * and the body will not be updated if it is attempted to integrate its
 * equations of motion.
 */
void RigidBody::set_enabled(bool flag)
{
  // mark as enabled / disabled
  _enabled = flag;  

  // if disabled, then zero the velocities and momenta; also zero inverse mass and inverse inertia tensor
  if (!_enabled)
  {
    _xd = ZEROS_3;
    _omega = ZEROS_3;
    _inv_mass = 0;
    _invJ = ZEROS_3x3;
  }
  else
  {
    _inv_mass = 1.0/_mass;
    _invJ = Matrix3::inverse(_J);
  }
}

/// Sets the angular velocity (world frame), optionally updating the angular momentum
void RigidBody::set_avel(const Vector3& avel)
{
  // angular velocity should not be set for disabled bodies
  if (!_enabled)
    return;  
  
  // set the angular velocity
  assert(avel.is_finite());
  _omega = avel;

  // invalidate articulated body velocity, if necessary
  invalidate_velocity();
}

/// Sets the linear velocity of this body
void RigidBody::set_lvel(const Vector3& lvel) 
{ 
  if (!_enabled)
    return;
  
  // set the linear velocity  
  assert(lvel.is_finite());
  _xd = lvel; 

  // invalidate articulated body velocity, if necessary
  invalidate_velocity();
}

/// Sets the inertia tensor for this body
void RigidBody::set_inertia(const Matrix3& inertia)
{
  // check the validity of the inertia matrix
  if (_enabled)
  {
    // first, verify that the matrix is symmetric
    if ((inertia - Matrix3::transpose(inertia)).norm_inf() > NEAR_ZERO)
    {
      std::cerr << "RigidBody::set_inertia() warning - inertia matrix does";
      std::cerr << "not appear to be symmetric!" << std::endl;
      std::cerr << "Inertia matrix: " << std::endl << inertia;
    }
    else
    { 
      // get the eigenvalues of the matrix
      VectorN evals;
      MatrixN J(inertia);
      LinAlg::eig_symm(J, evals);

      // first, verify that all eigenvalues are non-negative
      Real min_eig = *std::min_element(evals.begin(), evals.end());
      if (min_eig < (Real) 0.0)
      {
        std::cerr << "RigidBody::set_inertia() warning - minimum eigenvalue of inertia matrix is " << std::endl;
        std::cerr << min_eig << " -- matrix is *at best* PSD" << std::endl;
      }
      else
      {
        // all checks ok to here. still, look for relative inertial problems.
        Real max_eig = *std::max_element(evals.begin(), evals.end());
        if (max_eig / min_eig > (Real) 1000.0)
        {
          std::cerr << "RigidBody::set_inertia() warning - ratio of maximim:minimum eigenvalues of " << std::endl;
          std::cerr << "inertia matrix is" << (max_eig/min_eig) << ":1; for numerical stability, " << std::endl;
          std::cerr << "< 100:1 is recommended." << std::endl;
        }
      }
    } 
  }

  // store the matrix and its inverse
  _J = inertia;
  
  // set the inverse inertia to zero if body is disabled; to inverse of inertia otherwise
  if (_enabled)
    _invJ = Matrix3::inverse(inertia);
  else
    _invJ = ZEROS_3x3;

  // invalidate position, just in case
  invalidate_position();
}

/// Gets the spatial isolated inertia in link coordinates
SpatialRBInertia RigidBody::get_spatial_iso_inertia(ReferenceFrameType rftype) const
{
  if (rftype == eLink)
    return SpatialRBInertia(_mass, ZEROS_3, _J);
  else
  {
    assert(rftype == eGlobal);
    return get_spatial_transform_link_to_global().transform(_mass, _J);
  }
}

/// Sets the mass of this body
void RigidBody::set_mass(Real mass)
{
  if (mass < 0.0)
    throw std::runtime_error("Called RigidBody::set_mass() with negative mass");
  _mass = mass;
  _inv_mass = (_enabled && mass != (Real) 0.0) ? (Real) 1.0/mass : (Real) 0.0;

  // invalidate position, just in case
  invalidate_position();
}

/// Sets the position for this rigid body
/**
 * Also updates the transforms for associated visualization and collision data.
 */
void RigidBody::set_position(const Vector3& v)
{
  // store the position
  _x = v;
  _F.set_translation(v);
  
  // synchronize the geometry
  synchronize();

  // invalidate the position
  invalidate_position();
}

/// Sets the orientation for this rigid body
/**
 * Also updates the transforms for associated visualization and collision data.
 */
void RigidBody::set_orientation(const Quat& q)
{
  // store the orientation
  _q = q;
  _F.set_rotation(&q);
 
  // synchronize the geometry
  synchronize();

  // invalidate the position
  invalidate_position();
}

/// Sets the current 4x4 homogeneous transformation for this body using quaternion and vector position
void RigidBody::set_transform(const Quat& q, const Vector3& x)
{
  // compute the transform
  _F = Matrix4(&q, &x);

  // save the position and orientation
  _x = x;
  _q = q;

  // synchronize the geometry
  synchronize();

  // invalidate the position
  invalidate_position();
}

/// Sets the current 4x4 homogeneous transformation for this rigid body
/**
 * Also updates the transforms for associated visualization and collision data.
 */
void RigidBody::set_transform(const Matrix4& T) 
{ 
  // store the transform
  _F = T; 
  
  // update position and quaternion orientation (redundant)
  _x = T.get_translation(); 
  _q = Quat(&T);

  // synchronize the geometry
  synchronize();

  // invalidate the position
  invalidate_position();
}

/// Adds a force at a particular point on the body
/**
 * The frame of the applied force is centered at p and is aligned with the 
 * global frame.
 * \param f the force
 * \param p the position on the body in world coordinates (i.e., not relative to the center-of-mass)
 */
void RigidBody::add_force(const Vector3& f, const Vector3& p)
{
  // do not add forces to disabled bodies
  if (!_enabled)
    return;
  
  // add the force directly to the c.o.m.
  _forces += f;
  
  // compute and add the torque
  _torques += Vector3::cross(p - _x, f);
}

/// Synchronizes associated collision mesh transforms with this transform
void RigidBody::synchronize()
{
  BOOST_FOREACH(CollisionGeometryPtr g, geometries)
    g->set_transform(_F, false);
}

/// Calculates the velocity of a point on this rigid body
Vector3 RigidBody::calc_point_vel(const Vector3& point) const
{
  // if the body is disabled, point velocity is zero
  if (!_enabled)
    return ZEROS_3;

  // compute the arm
  Vector3 arm = point - _x;

  // determine the velocity of the point
  return _xd + Vector3::cross(_omega, arm);
}

/// Calculates the separation acceleration of a common point on two bodies along the desired direction
Real RigidBody::calc_sep_accel(RigidBody& rb1, RigidBody& rb2, const Vector3& point, const Vector3& dir, const Vector3& dir_dot, Real dt)
{
  // compute forward dynamics of the two bodies
  rb1.calc_fwd_dyn(dt);
  rb2.calc_fwd_dyn(dt);

  // compute the arms
  Vector3 arm1 = point - rb1.get_position();
  Vector3 arm2 = point - rb2.get_position();

  // compute the point acceleration on body 1
  const Vector3& av1 = rb1.get_avel();
  Vector3 la1 = rb1.get_laccel();
  Vector3 aa1 = rb1.get_aaccel();
  Vector3 pa1 = la1 + Vector3::cross(aa1, arm1) + Vector3::cross(av1, Vector3::cross(av1, arm1));  

  // compute the point acceleration on body 2
  const Vector3& av2 = rb2.get_avel();
  Vector3 la2 = rb2.get_laccel();
  Vector3 aa2 = rb2.get_aaccel();
  Vector3 pa2 = la2 + Vector3::cross(aa2, arm2) + Vector3::cross(av2, Vector3::cross(av2, arm2));

  // compute the point velocity on body 1
  const Vector3& lv1 = rb1.get_lvel();
  Vector3 pv1 = lv1 + Vector3::cross(av1, arm1);

  // compute the point velocity on body 2
  const Vector3& lv2 = rb2.get_lvel();
  Vector3 pv2 = lv2 + Vector3::cross(av2, arm2);

  // compute the separation acceleration
  Real accel = Vector3::dot(dir, pa1 - pa2) + 2.0 * Vector3::dot(dir_dot, pv1 - pv2);

  return accel;
}

/// Calculates the acceleration at a point (typically somewhere on this body) along the desired direction
/**
 * \param point the point in the global frame
 * \param dir the direction in which acceleration is to be computed
 * \param dir_dot the derivative of the direction in which acceleration is to
 *        be computed
 * \return the acceleration along dir
 * \note <b>this</b> may modify the body, because it will trigger forward
 *         dynamics via calls to get_laccel() and get_raccel().
 */
Real RigidBody::calc_point_accel(const Vector3& point, const Vector3& dir, Real dt)
{
  // if the body is disabled, point acceleration is zero
  if (!_enabled)
    return 0.0;

  // compute forward dynamics 
  calc_fwd_dyn(dt);

  // get the accelerations for the body
  Vector3 la = get_laccel();
  Vector3 aa = get_aaccel();
      
  // determine the moment arm
  Vector3 arm = point - _x;

  // get the acceleration projected along the direction at the point
  Real accel = Vector3::dot(dir, la + Vector3::cross(aa, arm) + Vector3::cross(_omega, Vector3::cross(_omega, arm)));

  return accel;
}

/// Implements Base::load_from_xml()
void RigidBody::load_from_xml(XMLTreeConstPtr node, map<std::string, BasePtr>& id_map)
{
  const unsigned X = 0, Y = 1, Z = 2;
  map<std::string, BasePtr>::const_iterator id_iter;

  // load parent data
  SingleBody::load_from_xml(node, id_map);

  // ***********************************************************************
  // don't verify that the node is correct, b/c RigidBody can be subclassed
  // ***********************************************************************
 
  // read the Coulomb dampening coefficient, if provided
  const XMLAttrib* coulomb_coeff_attr = node->get_attrib("coulomb-dampening-coeff");
  if (coulomb_coeff_attr)
    coulomb_coeff_attr->get_vector_value(coulomb_coeff);

  // read the viscous dampening coefficient, if provided
  const XMLAttrib* viscous_coeff_attr = node->get_attrib("viscous-dampening-coeff");
  if (viscous_coeff_attr)
    viscous_coeff_attr->get_vector_value(viscous_coeff);
 
  // read the linear acceleration, if provided
  const XMLAttrib* laccel_attr = node->get_attrib("linear-accel");
  if (laccel_attr)
    laccel_attr->get_vector_value(_xdd);

  // read the angular acceleration, if provided
  const XMLAttrib* aaccel_attr = node->get_attrib("angular-accel");
  if (aaccel_attr)
    aaccel_attr->get_vector_value(_alpha);

  // read whether the body is enabled, if provided
  const XMLAttrib* enabled_attr = node->get_attrib("enabled");
  if (enabled_attr)
    _enabled = enabled_attr->get_bool_value();

  // read the mass, if provided
  const XMLAttrib* mass_attr = node->get_attrib("mass");
  if (mass_attr)
    set_mass(mass_attr->get_real_value());

  // read the inertia matrix, if provided
  const XMLAttrib* inertia_attr = node->get_attrib("inertia");
  if (inertia_attr)
  {
    Matrix3 J;
    inertia_attr->get_matrix_value(J);
    set_inertia(J);
  }

  // read the forces on the body, if provided
  const XMLAttrib* sumf_attr = node->get_attrib("sum-forces");
  if (sumf_attr)
    sumf_attr->get_vector_value(_forces);

  // read the sum of torques on the body, if provided
  const XMLAttrib* sumt_attr = node->get_attrib("sum-torques");
  if (sumt_attr)
    sumt_attr->get_vector_value(_torques);

  // set the collision geometries, if provided
  list<XMLTreeConstPtr> cg_nodes = node->find_child_nodes("CollisionGeometry");
  if (!cg_nodes.empty())
  {
    // ok to clear the set of geometries
    geometries.clear();

    // read in the collision geometries
    for (list<XMLTreeConstPtr>::const_iterator i = cg_nodes.begin(); i != cg_nodes.end(); i++)
    {
      // create a new CollisionGeometry object
      CollisionGeometryPtr cg(new CollisionGeometry());

      // set the single body for the geometry
      cg->set_single_body(get_this());

      // populate the CollisionGeometry object
      cg->load_from_xml(*i, id_map);

      // add the collision geometry
      geometries.push_back(cg);
    }
  }

  // look for a inertia from primitives nodes
  // NOTE: we must do this step *before* setting velocities b/c setting 
  // velocities updates momenta!
  list<XMLTreeConstPtr> ifp_nodes = node->find_child_nodes("InertiaFromPrimitive");
  if (!ifp_nodes.empty())
  {
    // set com to zero initially
    Real mass = 0;
    Matrix3 J = ZEROS_3x3;
    Vector3 com = ZEROS_3;

    // loop over all InertiaFromPrimitive nodes
    for (list<XMLTreeConstPtr>::const_iterator i = ifp_nodes.begin(); i != ifp_nodes.end(); i++)
    {
      // make sure the child node has the ID
      const XMLAttrib* pid_attr = (*i)->get_attrib("primitive-id");
      if (!pid_attr)
      {
        cerr << "RigidBody::load_from_xml() - InertiaFromPrimitive node "; 
        cerr << "has no" << endl << "  primitive-id attribute!";
        cerr << endl << "  offending node: " << endl << *node;
        continue;
      }

      // get the ID
      const std::string& ID = pid_attr->get_string_value();

      // attempt to find the ID
      if ((id_iter = id_map.find(ID)) == id_map.end())
      {
        cerr << "RigidBody::load_from_xml() - Primitive id: ";
        cerr << ID << " not found!" << endl << "  offending node: ";
        cerr << endl << *node;
        continue;
      }

      // set the relative transformation to identity for this primitive
      Matrix4 TR = IDENTITY_4x4; 

      // read the relative transformation, if specified
      const XMLAttrib* rel_transform_attr = (*i)->get_attrib("rel-transform");
      if (rel_transform_attr)
      {
        rel_transform_attr->get_matrix_value(TR);
        if (!Matrix4::valid_transform(TR))
        {
          cerr << "RigidBody::load_from_xml() warning: invalid transform ";
          cerr << endl << TR << " when reading node " << endl;
          cerr << *node << endl;
          cerr << "  --> possibly a floating-point error..." << endl;
        }
      }

      // get the primitive
      PrimitivePtr primitive = dynamic_pointer_cast<Primitive>(id_iter->second);

      // get the mass, inertia, and com from the primitive
      Real massx = primitive->get_mass();
      Matrix3 Jx = primitive->get_inertia();
      Vector3 comx = primitive->get_com();

      // transform the inertia
      Primitive::transform_inertia(massx, Jx, comx, TR, Jx, comx);

      // update the mass, inertia, and com
      mass += massx;
      J += Jx;
      com += comx*massx;
    }

    // update the com
    assert(!_enabled || mass != 0.0);
    com *= (1.0/mass);

    if (com.norm() > NEAR_ZERO)
    {
      cerr << "RigidBody::load_from_xml() - center-of-mass of rigid body is not at its origin," << endl;
      cerr << " according to its geometry.  Dynamics may not reflect geometry of body." << endl;
      cerr << " body: " << id << "  center-of-mass: " << std::setprecision(12) << com[X] << " " << std::setprecision(12) << com[Y] << " " << std::setprecision(12) << com[Z] << endl;
    }

    // set the mass and inertia of the RigidBody additively
    set_mass(get_mass() + mass);
    set_inertia(get_inertia() + J);
  }

  // read the position of the body, if provided
  const XMLAttrib* pos_attr = node->get_attrib("position");
  if (pos_attr)
  {
    Vector3 position;
    pos_attr->get_vector_value(position);
    set_position(position);
  }

  // read the orientation of the body as an axis-angle, if provided
  const XMLAttrib* aangle_attr = node->get_attrib("aangle");
  if (aangle_attr)
  {
    VectorN aangle;
    aangle_attr->get_vector_value(aangle);
    AAngle aa(&aangle);
    set_orientation(Quat(&aa));
  }
  
  // read the transform, if provided; note that we do this after reading the
  // collision geometries to allow synchronization to occur
  const XMLAttrib* T_attr = node->get_attrib("transform");
  if (T_attr)
  {
    Matrix4 T;
    T_attr->get_matrix_value(T);
    if (!Matrix4::valid_transform(T))
    {
      cerr << "RigidBody::load_from_xml() warning: invalid transform ";
      cerr << endl << T << " when reading node " << endl;
      cerr << *node << endl;
      cerr << "  --> possibly a floating-point error..." << endl;
    }
    set_transform(T);
  }

  // read the linear velocity of the body, if provided; note that we update
  // the linear momentum too, so the mass must have been read first
  const XMLAttrib* lvel_attr = node->get_attrib("linear-velocity");
  if (lvel_attr)
  {
    Vector3 lv;
    lvel_attr->get_vector_value(lv);
    set_lvel(lv);  
  }

  // read the angular velocity of the body, if provided; note that we update
  // the angular momentum too, so the inertia must have been read first
  const XMLAttrib* avel_attr = node->get_attrib("angular-velocity");
  if (avel_attr)
  {
    Vector3 av;
    avel_attr->get_vector_value(av);
    set_avel(av);  
  }

/*
  // read in the vector from the inner joint to the com in link coordinates
  const XMLAttrib* d_attr = node->get_attrib("inner-joint-to-com-vector-link");
  if (d_attr)
  {
    Vector3 d;
    d_attr->get_vector_value(d);
    set_inner_joint_to_com_vector_link(d);
  }
*/
  // read the articulated body, if given
  const XMLAttrib* ab_attr = node->get_attrib("articulated-body-id");
  if (ab_attr)
  {
    // get the ID
    const std::string& ID = ab_attr->get_string_value();

    // look for the ID -- only warn if it is not found
    if ((id_iter = id_map.find(ID)) == id_map.end())
    {
      FILE_LOG(LOG_DYNAMICS) << "RigidBody::load_from_xml() warning - ";
      FILE_LOG(LOG_DYNAMICS) << "articulated body" << endl << "  '" << ID << "' not ";
      FILE_LOG(LOG_DYNAMICS) << "found" << endl << "  ** This warning could result ";
      FILE_LOG(LOG_DYNAMICS) << "from links being constructed before articulated bodies ";
      FILE_LOG(LOG_DYNAMICS) << endl << "    and may not be serious..." << endl; 
      FILE_LOG(LOG_DYNAMICS) << "  offending node: " << endl << *node;
    }
    else
      set_articulated_body(dynamic_pointer_cast<ArticulatedBody>(id_iter->second));
  }
/*
  // get all child links, if specified
  list<XMLTreeConstPtr> child_nodes = node->find_child_nodes("ChildLink");
  if (!child_nodes.empty())
    _child_links.clear();

  for (list<XMLTreeConstPtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
  {
    // look for the link-id and com-to-outboard-vec-link attributes
    const XMLAttrib* link_id_attr = (*i)->get_attrib("link-id");
    const XMLAttrib* ovec_attr = (*i)->get_attrib("com-to-outer-vec-link");
    if (!link_id_attr || !ovec_attr)
    {
      FILE_LOG(LOG_DYNAMICS) << "RigidBody::load_from_xml() - ChildLink node ";
      FILE_LOG(LOG_DYNAMICS) << "missing" << endl << "  link-id and/or ";
      FILE_LOG(LOG_DYNAMICS) << "com-to-outer-vec-link attributes in offending node: ";
      FILE_LOG(LOG_DYNAMICS) << endl << *node;
      continue;
    }

    // get the link ID and the com-to-outer vector
    const std::string& ID = link_id_attr->get_string_value();
    Vector3 com_to_outer;
    ovec_attr->get_vector_value(com_to_outer);

    // make sure that we can find the link ID
    if ((id_iter = id_map.find(ID)) == id_map.end())
    {
      FILE_LOG(LOG_DYNAMICS) << "RigidBody::load_from_xml() - child link-id ";
      FILE_LOG(LOG_DYNAMICS) << ID << " not found in" << endl << "  offending node: ";
      FILE_LOG(LOG_DYNAMICS) << endl << "  ** This warning could result from links ";
      FILE_LOG(LOG_DYNAMICS) << " being constructed before articulated bodies ";
      FILE_LOG(LOG_DYNAMICS) << endl << "    and may not be serious..." << endl; 
      cerr << endl << *node;
    }
    else
    {
      // add the link and the vector
      RigidBodyPtr link = dynamic_pointer_cast<RigidBody>(id_iter->second);
      add_child_link_link(link, com_to_outer);
    }
  }
*/
}

/// Implements Base::save_to_xml()
void RigidBody::save_to_xml(XMLTreePtr node, list<BaseConstPtr>& shared_objects) const
{
  // save parent data
  SingleBody::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "RigidBody";

  // save the linear acceleration
  node->attribs.insert(XMLAttrib("linear-accel", _xdd)); 

  // save the angular acceleration
  node->attribs.insert(XMLAttrib("angular-accel", _alpha)); 

  // save whether the body is enabled
  node->attribs.insert(XMLAttrib("enabled", _enabled));

  // save the mass
  node->attribs.insert(XMLAttrib("mass", _mass));

  // save the inertia
  node->attribs.insert(XMLAttrib("inertia", _J));

  // save the current transform
  node->attribs.insert(XMLAttrib("transform", _F));

  // save the linear velocity
  node->attribs.insert(XMLAttrib("linear-velocity", _xd));

  // save the cumulative forces on the body
  node->attribs.insert(XMLAttrib("sum-forces", _forces));

  // save the cumulative torques on the body
  node->attribs.insert(XMLAttrib("sum-torques", _torques));

  // save the dampening coefficients
  node->attribs.insert(XMLAttrib("coulomb-coeff", coulomb_coeff));
  node->attribs.insert(XMLAttrib("viscous-coeff", viscous_coeff));

  // save all collision geometries
  BOOST_FOREACH(CollisionGeometryPtr g, geometries)
  {
    XMLTreePtr geom_node(new XMLTree("CollisionGeometry"));
    node->add_child(geom_node);
    g->save_to_xml(geom_node, shared_objects);
  }

  // save the ID articulated body (if any)
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody(_abody);
    node->attribs.insert(XMLAttrib("articulated-body-id", abody->id));
  }
/*
  // save the IDs of all child links and the vectors from the com to the outer
  // joints
  for (map<RigidBodyPtr, Vector3>::const_iterator i = _child_links.begin(); i != _child_links.end(); i++)
  {
    XMLTreePtr child_node(new XMLTree("ChildLink"));
    node->add_child(child_node);
    child_node->attribs.insert(XMLAttrib("link-id", i->first->id));
    child_node->attribs.insert(XMLAttrib("com-to-outer-vec-link", i->second));
  }
*/
}

/// Adds an inner joint for this link
/**
 * \param parent the outer link of the parent
 * \param j the joint connecting parent and this
 * \param joint_to_com_vec_joint the vector from j to the center-of-mass of this
 *        (specified in the joint frame) 
 * \param joint_to_com_vec_link the vector from j to the center-of-mass of this
 *        (specified in this frame) 
 */
void RigidBody::add_inner_joint(RigidBodyPtr parent, JointPtr j, const Vector3& joint_to_com_vec_joint, const Vector3& joint_to_com_vec_link) 
{
  // remove the inner joint if it already exists
  for (list<InnerJointData>::iterator i = _inner_joints.begin(); i != _inner_joints.end(); )
    if (JointPtr(i->inner_joint) == j)
      i = _inner_joints.erase(i);
    else
      i++;

  // add the inner joint
  _inner_joints.push_back(InnerJointData());
  _inner_joints.back().parent = parent;
  _inner_joints.back().inner_joint = j;
  _inner_joints.back().joint_to_com_vec_jf = joint_to_com_vec_joint;
  _inner_joints.back().joint_to_com_vec_of = joint_to_com_vec_link;
 
  // set the joint to point to this, if it does not already do so
  RigidBodyPtr outboard = j->get_outboard_link();
  if (outboard != get_this())
    j->set_outboard_link(get_this());

  // update the spatial axes
  j->update_spatial_axes();

  // set the parent link / inboard link pointers, if possible
  if (j->get_inboard_link() != parent)
    j->set_inboard_link(parent);

  // set the articulated body / inner joint articulated body pointers, if
  // possible
  if (!j->get_articulated_body() && !_abody.expired())
    j->set_articulated_body(ArticulatedBodyPtr(_abody));
  else if (j->get_articulated_body() && _abody.expired())
    set_articulated_body(j->get_articulated_body());

  // again, the articulated body pointers must now be equal; it is
  // conceivable that the user is updating the art. body pointers in an
  // unorthodox manner, but we'll look for this anwyway...
  #ifndef NDEBUG
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody1 = j->get_articulated_body();
    ArticulatedBodyPtr abody2(_abody);
    assert(abody1 == abody2);
  }
  #endif
}

/// Adds an outer joint for this link
/**
 * \param child the child link
 * \param j the joint connecting this and child
 * \param com_to_joint_vec_link the vector from the center-of-mass of this link
 *        to j in this link's coordinates 
 * \note replaces the outer joint if it is already attached to this link 
 */
void RigidBody::add_outer_joint(RigidBodyPtr child, JointPtr j, const Vector3& com_to_joint_vec_link) 
{
  // remove the outer joint if it already exists
  for (list<OuterJointData>::iterator i = _outer_joints.begin(); i != _outer_joints.end(); )
    if (JointPtr(i->outer_joint) == j)
      i = _outer_joints.erase(i);
    else
      i++;

  // add the outer joint
  _outer_joints.push_back(OuterJointData());
  _outer_joints.back().child = child;
  _outer_joints.back().outer_joint = j;
  _outer_joints.back().com_to_joint_vec = com_to_joint_vec_link;
 
  // set the joint to point to this, if it does not already do so
  RigidBodyPtr inboard = j->get_inboard_link();
  if (inboard != get_this())
    j->set_inboard_link(get_this());

  // update the spatial axes
  j->update_spatial_axes();

  // set the child link / outboard link pointers, if possible
  if (j->get_outboard_link() != child)
    j->set_outboard_link(child);

  // set the articulated body / inner joint articulated body pointers, if
  // possible
  if (!j->get_articulated_body() && !_abody.expired())
    j->set_articulated_body(ArticulatedBodyPtr(_abody));
  else if (j->get_articulated_body() && _abody.expired())
    set_articulated_body(j->get_articulated_body());

  // again, the articulated body pointers must now be equal; it is
  // conceivable that the user is updating the art. body pointers in an
  // unorthodox manner, but we'll look for this anwyway...
  #ifndef NDEBUG
  if (!_abody.expired())
  {
    ArticulatedBodyPtr abody1 = j->get_articulated_body();
    ArticulatedBodyPtr abody2(_abody);
    assert(abody1 == abody2);
  }
  #endif
}

/// Determines whether the given link is a child link of this
bool RigidBody::is_child_link(RigidBodyConstPtr query) const
{
  BOOST_FOREACH(const OuterJointData& odata, _outer_joints)
    if (RigidBodyPtr(odata.child) == query)
      return true;

  return false;
}

/// Gets the i'th child link
RigidBodyPtr RigidBody::get_child_link(unsigned i) const
{
  unsigned index = 0;
  BOOST_FOREACH(const OuterJointData& odata, _outer_joints)
    if (index++ == i)
      return RigidBodyPtr(odata.child);

  // still here? not found...
  return RigidBodyPtr();
}

/// Determines whether the given link is a descendant of this
/**
 * \note returns <b>true</b> if query == this
 */
bool RigidBody::is_descendant_link(RigidBodyConstPtr query) const
{
  queue<RigidBodyConstPtr> q;

  // check for query == this
  if (query == shared_from_this())
    return true;

  // add all children to the queue
  BOOST_FOREACH(const OuterJointData& odata, _outer_joints)
    q.push(RigidBodyConstPtr(odata.child));

  // continue processing children until no more children are able to be processed
  while (!q.empty())
  {
    RigidBodyConstPtr link = q.front();
    q.pop();
    if (link == query)
      return true;
    BOOST_FOREACH(const OuterJointData& odata, link->_outer_joints)
      q.push(RigidBodyConstPtr(odata.child));
  }    

  return false;
}

RigidBody::InnerJointData& RigidBody::get_inner_joint_data(RigidBodyPtr parent)
{
  BOOST_FOREACH(InnerJointData& i, _inner_joints)
    if (RigidBodyPtr(i.parent) == parent)
      return i;

  throw std::runtime_error("Requested inner joint data was not found!");
}

RigidBody::InnerJointData& RigidBody::get_inner_joint_data(JointPtr joint)
{
  BOOST_FOREACH(InnerJointData& i, _inner_joints)
    if (JointPtr(i.inner_joint) == joint)
      return i;

  throw std::runtime_error("Requested inner joint data was not found!");
}

RigidBody::OuterJointData& RigidBody::get_outer_joint_data(RigidBodyPtr child)
{
  BOOST_FOREACH(OuterJointData& o, _outer_joints)
    if (RigidBodyPtr(o.child) == child)
      return o;

  throw std::runtime_error("Requested outer joint data was not found!");
}

RigidBody::OuterJointData& RigidBody::get_outer_joint_data(JointPtr joint)
{
  BOOST_FOREACH(OuterJointData& o, _outer_joints)
    if (JointPtr(o.outer_joint) == joint)
      return o;

  throw std::runtime_error("Requested outer joint data was not found!");
}

/// Removes all joints that connect the specified child link to this link
/**
 * Returns true if the link was found. 
 */
bool RigidBody::remove_outer_joints(RigidBodyPtr child)
{
  bool found_one = false;

  for (std::list<OuterJointData>::iterator i = _outer_joints.begin(); i != _outer_joints.end(); )
    if (RigidBodyPtr(i->child) == child)
    {
      i = _outer_joints.erase(i);
      found_one = true;
    }
    else
      i++;

  return found_one; 
}

/// Removes the specified outer joint from this link
/**
 * Returns true if the link was found. 
 */
bool RigidBody::remove_outer_joint(JointPtr joint)
{
  for (std::list<OuterJointData>::iterator i = _outer_joints.begin(); i != _outer_joints.end(); i++)
    if (JointPtr(i->outer_joint) == joint)
    {
      _outer_joints.erase(i);
      return true;
    }

  return false; 
}

/// Removes all joints that connect the specified parent link to this link
/**
 * Returns true if the link was found. 
 */
bool RigidBody::remove_inner_joints(RigidBodyPtr parent)
{
  bool found_one = false;

  for (std::list<InnerJointData>::iterator i = _inner_joints.begin(); i != _inner_joints.end(); )
    if (RigidBodyPtr(i->parent) == parent)
    {
      i = _inner_joints.erase(i);
      found_one = true;
    }
    else
      i++;

  return found_one; 
}

/// Removes the specified outer joint from this link
/**
 * Returns true if the link was found. 
 */
bool RigidBody::remove_inner_joint(JointPtr joint)
{
  for (std::list<InnerJointData>::iterator i = _inner_joints.begin(); i != _inner_joints.end(); i++)
    if (JointPtr(i->inner_joint) == joint)
    {
      _inner_joints.erase(i);
      return true;
    }

  return false; 
}

/// Adds a force to the center-of-mass of this link
/**
 * The frame of the applied force is centered at the center-of-mass of the
 * body and is aligned with the global frame.
 */
void RigidBody::add_force(const Vector3& force) 
{
  if (_enabled)
    _forces += force;
}

/// Adds a torque to the center-of-mass of this body 
/**
 * The frame of the applied torque is centered at the center-of-mass of the
 * body and is aligned with the global frame.
 */
void RigidBody::add_torque(const Vector3& torque) 
{
  if (_enabled)
    _torques += torque;
}

/// Applies a linear impulse to this link at the specified point
/**
 * \param j the linear component of the impulse
 * \param k the angular component of the impulse
 * \param p the point (global frame) at which the impulse is to be applied
 */
void RigidBody::apply_impulse(const Vector3& j, const Vector3& p)
{  
  // if this is not an articulated body, just update linear and angular
  // momenta and velocites
  if (_abody.expired())
  {
    if (!_enabled)
      return;

    // update linear and angular velocities 
    _xd += j * _inv_mass;
    Matrix3 R(&_q);
    _omega +=  R * _invJ * R.transpose_mult(Vector3::cross(p - _x, j));
    assert(_xd.is_finite());
    assert(_omega.is_finite());

    // reset the force and torque accumulators for this body
    _forces = ZEROS_3;
    _torques = ZEROS_3; 
  }
  else
  {
    // get the articulated body
    ArticulatedBodyPtr abody(_abody);
  
    // apply the impulse to the articulated body
    abody->apply_impulse(j, ZEROS_3, p, get_this());
  }
}

/// Applies linear and angular impulses to this link at the specified point
/**
 * \param j the linear component of the impulse
 * \param k the angular component of the impulse
 * \param p the point (global frame) at which the impulse is to be applied
 */
void RigidBody::apply_impulse(const Vector3& j, const Vector3& k, const Vector3& p)
{
  // if this is not an articulated body, just update linear and angular
  // momenta and velocites
  if (_abody.expired())
  {
    if (!_enabled)
      return;

    // update linear and angular velocities 
    _xd += j * _inv_mass;
    Matrix3 R(&_q);
    _omega +=  R * _invJ * R.transpose_mult(k + Vector3::cross(p - _x, j));
    assert(_xd.is_finite());
    assert(_omega.is_finite());

    // reset the force and torque accumulators for this body
    _forces = ZEROS_3;
    _torques = ZEROS_3; 
  }
  else
  {
    // get the articulated body
    ArticulatedBodyPtr abody(_abody);
  
    // apply the impulse to the articulated body
    abody->apply_impulse(j, k, p, get_this());
  }
}

/// Gets the spatial velocity for this link
SVector6 RigidBody::get_spatial_velocity(ReferenceFrameType rftype) 
{
  // get transformed linear and angular velocity
  Vector3 omega = _F.transpose_mult_vector(get_avel());
  Vector3 xd = _F.transpose_mult_vector(get_lvel());

  // form the spatial velocity vector in link coordinates
  SVector6 v(omega, xd);

  if (rftype == eLink)
    return v;
  else
    return get_spatial_transform_link_to_global().transform(v);
}

/// Gets the spatial acceleration for this link 
SVector6 RigidBody::get_spatial_accel(ReferenceFrameType rftype) 
{
  // get the transpose of the rotation matrix
  Matrix3 RT;
  get_transform().get_rotation(&RT);
  RT.transpose();

  // form the spatial acceleration vector in link coordinates
  SVector6 a(RT * _alpha, RT * _xdd);

  if (rftype == eLink)
    return a;
  else
    return get_spatial_transform_link_to_global().transform(a);
}

/// Sets the spatial velocity 
void RigidBody::set_spatial_velocity(const SVector6& v, ReferenceFrameType rftype)
{
  // get the rotation matrix of this link
  Matrix3 R;
  get_transform().get_rotation(&R);

  // set the linear and angular velocities in the global frame
  if (rftype == eLink)
  {
    RigidBody::set_lvel(R * v.get_lower());
    RigidBody::set_avel(R * v.get_upper());
  }
  else
  {
    SVector6 vi = get_spatial_transform_global_to_link().transform(v);
    RigidBody::set_lvel(R * vi.get_lower());
    RigidBody::set_avel(R * vi.get_upper());
  }
}

/// Sets the spatial acceleration 
void RigidBody::set_spatial_accel(const SVector6& a, ReferenceFrameType rftype)
{
  // get the rotation matrix of this link
  Matrix3 R;
  get_transform().get_rotation(&R);

  // set the linear and angular acceleration in the global frame
  if (rftype == eLink)
  {
    _xdd = R * a.get_lower();
    _alpha = R * a.get_upper();
  }
  else
  {
    SVector6 ai = get_spatial_transform_global_to_link().transform(a);
    _xdd = R * ai.get_lower();
    _alpha = R * ai.get_upper();
  }
}

/// Returns the spatial transform from this link's parent to this link
SpatialTransform RigidBody::get_spatial_transform_forward() const
{
  if (_inner_joints.size() > 1)
    throw std::runtime_error("RigidBody::get_spatial_transform_forward() does not work with closed chains");
  RigidBodyPtr p = (_inner_joints.empty()) ? RigidBodyPtr() : RigidBodyPtr(_inner_joints.front().parent);
  const Matrix4& parent_T = (!p) ? IDENTITY_4x4 : p->get_transform();
  return SpatialTransform(parent_T, get_transform());
}

/// Returns the spatial transform from this link to this link's parent
SpatialTransform RigidBody::get_spatial_transform_backward() const
{
  if (_inner_joints.size() > 1)
    throw std::runtime_error("RigidBody::get_spatial_transform_forward() does not work with closed chains");
  RigidBodyPtr p = (_inner_joints.empty()) ? RigidBodyPtr() : RigidBodyPtr(_inner_joints.front().parent);
  const Matrix4& parent_T = (!p) ? IDENTITY_4x4 : p->get_transform();
  return SpatialTransform(get_transform(), parent_T);
}

/// Returns the spatial transform from this link to the global coordinate system
SpatialTransform RigidBody::get_spatial_transform_link_to_global() const
{
  return SpatialTransform::to_global(_F);
}

/// Returns the spatial transform from the global coordinates system to this link
/**
 * \note the spatial transform will be computed if necessary
 */
SpatialTransform RigidBody::get_spatial_transform_global_to_link() const
{
  return SpatialTransform::from_global(_F);
}

/// Adds a generalized force to this rigid body
void RigidBody::add_generalized_force(GeneralizedCoordinateType gctype, const VectorN& gf)
{
  // if body is not enabled, do nothing
  if (!_enabled)
    return;

  assert(gf.size() == num_generalized_coordinates(gctype));

  // update force accumulator
  _forces += Vector3(gf[0], gf[1], gf[2]);

  // determine the proper generalized coordinate type
  switch (gctype)
  {
    case eRodrigues:
      _torques += _q.G_mult(gf[3], gf[4], gf[5], gf[6]) * (Real) 0.5;
      break;

    case eAxisAngle:
      _torques += Vector3(gf[3], gf[4], gf[5]);
      break;
  }
}

/// Applies a generalized impulse to this rigid body
void RigidBody::apply_generalized_impulse(GeneralizedCoordinateType gctype, const VectorN& gj)
{
  // don't do anything if this body is disabled
  if (!_enabled)
    return;

  // simple error check...
  assert(gj.size() == num_generalized_coordinates(gctype));

  // clear the force and torque accumulators
  _forces = ZEROS_3;
  _torques = ZEROS_3;

  // look for easy case (axis-angle)
  if (gctype == DynamicBody::eAxisAngle)
  {
    // get the impulses
    Vector3 j(gj[0], gj[1], gj[2]);
    Vector3 k(gj[3], gj[4], gj[5]);

    // determine the change in linear velocity
    _xd += j * _inv_mass;

    // determine the change in angular velocity
    Matrix3 R(&_q);
    Matrix3 invJ = R * _invJ * Matrix3::transpose(R);
    _omega += invJ * k;
    invalidate_velocity();
  }
  else
  {
    assert(gctype == DynamicBody::eRodrigues);

    // get proper generalized inertia matrix
    MatrixN M;
    get_generalized_inertia(gctype, M);
    VectorN qd_delta = gj;
    LinAlg::solve_fast(M, qd_delta);
    VectorN gv;
    get_generalized_velocity(gctype, gv);
    gv += qd_delta;
    set_generalized_velocity(gctype, gv); 
  }
}

/// Gets the generalized position of this rigid body
VectorN& RigidBody::get_generalized_coordinates(GeneralizedCoordinateType gctype, VectorN& gc) 
{
  // special case: disabled body
  if (!_enabled)
    return gc.resize(0);

  // make sure we were not passed the wrong type of gctype
  if (gctype == DynamicBody::eAxisAngle)
    throw std::runtime_error("Unable to use axis-angle for positional generalized coordinates!");
  else
  {
    // return the generalized position using Rodrigues parameters
    assert(gctype == DynamicBody::eRodrigues);
    gc.resize(num_generalized_coordinates(gctype));
    gc[0] = _x[0];
    gc[1] = _x[1];
    gc[2] = _x[2];
    gc[3] = _q.w;
    gc[4] = _q.x;
    gc[5] = _q.y;
    gc[6] = _q.z;
    return gc; 
  }
}

/// Sets the generalized coordinates of this rigid body
void RigidBody::set_generalized_coordinates(GeneralizedCoordinateType gctype, const VectorN& gc)
{
  // special case: disabled body
  if (!_enabled)
    return;

  // make sure we were not passed the wrong type of gctype
  if (gctype == DynamicBody::eAxisAngle)
    throw std::runtime_error("Unable to use axis-angle for positional generalized coordinates!");
  else
  {
    assert(gctype == DynamicBody::eRodrigues);
    assert(gc.size() == num_generalized_coordinates(gctype));

    // get the position
    Vector3 x(gc[0], gc[1], gc[2]);

    // get the unit quaternion
    _q.w = gc[3];
    _q.x = gc[4];
    _q.y = gc[5];
    _q.z = gc[6];

    // normalize the unit quaternion, just in case
    _q.normalize();

    // set the transform
    set_transform(_q, x);
  }
}

/// Sets the generalized velocity of this rigid body
void RigidBody::set_generalized_velocity(GeneralizedCoordinateType gctype, const VectorN& gv)
{
  // special case: disabled body
  if (!_enabled)
    return;

  // set the linear velocity 
  _xd = Vector3(gv[0], gv[1], gv[2]);

  // determine the proper type
  if (gctype == DynamicBody::eAxisAngle)
  {
    // set omega
    _omega = Vector3(gv[3], gv[4], gv[5]);

    // invalidate the velocities
    invalidate_velocity();
  }
  else
  {
    assert(gctype == DynamicBody::eRodrigues);

    // get the quaternion derivatives
    Quat qd;
    qd.w = gv[3];
    qd.x = gv[4];
    qd.y = gv[5];
    qd.z = gv[6];

    // set omega
    _omega = Quat::to_omega(_q, qd);

    // invalidate the velocities
    invalidate_velocity();
  }
}

/// Gets the generalized velocity of this rigid body
VectorN& RigidBody::get_generalized_velocity(GeneralizedCoordinateType gctype, VectorN& gv) 
{
  // special case: disabled body
  if (!_enabled)
    return gv.resize(0);

  // resize the generalized velocity
  gv.resize(num_generalized_coordinates(gctype));

  // setup the linear components
  gv[0] = _xd[0];
  gv[1] = _xd[1];
  gv[2] = _xd[2];

  // determine the proper generalized coordinate type
  if (gctype == DynamicBody::eAxisAngle)
  {
    gv[3] = _omega[0];
    gv[4] = _omega[1];
    gv[5] = _omega[2];
    return gv;
  }
  else
  {
    assert(gctype == DynamicBody::eRodrigues);

    // compute the Quaternion derivatives
    Quat qd = Quat::deriv(_q, _omega);

    // setup the generalized velocity
    gv[3] = qd.w;
    gv[4] = qd.x;
    gv[5] = qd.y;
    gv[6] = qd.z; 
    return gv;
  }
}

/// Gets the generalized acceleration of this body
VectorN& RigidBody::get_generalized_acceleration(GeneralizedCoordinateType gctype, VectorN& ga)
{
  // special case: body is disabled
  if (!_enabled)
    return ga.resize(0);

  // setup the linear components
  ga.resize(num_generalized_coordinates(gctype));
  ga[0] = _xdd[0];
  ga[1] = _xdd[1];
  ga[2] = _xdd[2];

  // determine the proper generalized coordinate type
  if (gctype == DynamicBody::eAxisAngle)
  {
    ga[3] = _alpha[0];
    ga[4] = _alpha[1];
    ga[5] = _alpha[2];
  }
  else
  {
    assert(gctype == DynamicBody::eRodrigues);

    // convert alpha to a quaternion derivative
    Quat qdd = Quat::dderiv(_q, _omega, _alpha);

    // setup generalized acceleration
    ga[3] = qdd.w;
    ga[4] = qdd.x;
    ga[5] = qdd.y;
    ga[6] = qdd.z;
  }
 
  return ga;
}

/// Gets the generalized inertia of this rigid body
MatrixN& RigidBody::get_generalized_inertia(GeneralizedCoordinateType gctype, MatrixN& M) 
{
  const unsigned X = 0, Y = 1, Z = 2;

  // special case: disabled body
  if (!_enabled)
    return M.resize(0,0);

  // init the matrix to all zeros
  const unsigned NGC = num_generalized_coordinates(gctype);
  M.set_zero(NGC,NGC);

  // upper left diagonal is mass
  M(X,X) = M(Y,Y) = M(Z,Z) = _mass;

  // use the proper generalized coordinate type
  if (gctype == DynamicBody::eAxisAngle)
  {
    Matrix3 R(&_q);
    Matrix3 Jinv = R * _invJ * Matrix3::transpose(R);
    M.set_sub_mat(3, 3, Jinv);
  }
  else
  {
    assert(gctype == DynamicBody::eRodrigues);

    // set lower 4x4 ([Nikravesh, 1988, p. 295])
    Vector3 ix = _J.get_row(X);
    Vector3 iy = _J.get_row(Y);
    Vector3 iz = _J.get_row(Z);
    Quat qx = _q.L_transpose_mult(ix) * (Real) 2.0;
    Quat qy = _q.L_transpose_mult(iy) * (Real) 2.0;
    Quat qz = _q.L_transpose_mult(iz) * (Real) 2.0;
    M(3,3) = qx.w;  M(3,4) = qx.x;  M(3,5) = qx.y;  M(3,6) = qx.z;
    M(4,3) = qy.w;  M(4,4) = qy.x;  M(4,5) = qy.y;  M(4,6) = qy.z;
    M(5,3) = qz.w;  M(5,4) = qz.x;  M(5,5) = qz.y;  M(5,6) = qz.z;
    M(6,3) = _q.w;  M(6,4) = _q.x;  M(6,5) = _q.y;  M(6,6) = _q.z;
  }

  return M;
}

/// Gets the generalized external forces
VectorN& RigidBody::get_generalized_forces(GeneralizedCoordinateType gctype, VectorN& f) 
{
  // special case: disabled body
  if (!_enabled)
    return f.resize(0);

  // resize the generalized forces vector
  const unsigned NGC = (gctype == DynamicBody::eAxisAngle) ? 6 : 7;
  f.resize(NGC);

  // setup the linear components of f
  f[0] = _forces[0];
  f[1] = _forces[1];
  f[2] = _forces[2];

  // compute inertial forces
  Vector3 tau = _torques - calc_inertial_forces();

   // use the proper generalized coordinate type
  if (gctype == DynamicBody::eAxisAngle)
  {
    f[3] = tau[0];
    f[4] = tau[1];
    f[5] = tau[2];
  }
  else
  {
    assert(gctype == DynamicBody::eRodrigues);

    // get the torque in the body's frame
    Matrix3 R(&_q);
    Vector3 tau = R.transpose_mult(tau);

    // determine quaternion parameters
    Quat qd = Quat::deriv(_q, _omega);

    // set the generalized forces
    f[3] = tau[0];
    f[4] = tau[1];
    f[5] = tau[2];
    f[6] = -qd.norm_sq();
  }

  return f;
}

/// Calculates the inertial forces on this body
Vector3 RigidBody::calc_inertial_forces() const
{
  Matrix3 R = _F.get_rotation();
  Matrix3 I = R * _J * Matrix3::transpose(R);
  return Vector3::cross(_omega, I*_omega);
}

/// Converts a force to a generalized force
VectorN& RigidBody::convert_to_generalized_force(GeneralizedCoordinateType gctype, SingleBodyPtr body, const Vector3& f, const Vector3& t, VectorN& gf) 
{
  // verify that body == this
  assert(body.get() == this);

  // special case: disabled body
  if (!_enabled)
    return gf.resize(0);

  // setup the linear components
  gf.resize(num_generalized_coordinates(gctype));
  gf[0] = f[0];
  gf[1] = f[1];
  gf[2] = f[2];

  // use the proper generalized coordinate type
  if (gctype == DynamicBody::eAxisAngle)
  {
    gf[3] = t[0];
    gf[4] = t[1];
    gf[5] = t[2];
  }
  else
  {
    assert(gctype == DynamicBody::eRodrigues);

    // convert the torque to the body's coordinate system
    Matrix3 R(&_q);
    Vector3 tau = R.transpose_mult(t);

    // setup the generalized force
    gf[3] = tau[0];
    gf[4] = tau[1];
    gf[5] = tau[2];
    gf[6] = (Real) 0.0;
  }

  return gf; 
}

/// Verifies that a given transform is valid
bool RigidBody::valid_transform(const MatrixN& T, Real tol)
{
  const unsigned X = 0, Y = 1, Z = 2, W = 3;

  // first make sure it is the proper size
  if (T.rows() || T.columns() != 4)
    return false;

  // now get the upper 3x3 and make sure that it is valid
  Matrix3 R(T.get_sub_mat(X, W, X, W).begin());
  if (!Matrix3::valid_rotation_scale(R))
    return false;

  // finally make sure that the last row is set correctly
  if (std::fabs(T(W,X)) > tol || std::fabs(T(W,Y)) > tol || std::fabs(T(W,Z)) > tol)
    return false;

  return (std::fabs(T(W,W) - 1.0) < tol);
}

/// Calculates the kinetic energy of the body
Real RigidBody::calc_kinetic_energy() const
{
  if (!_enabled)
    return (Real) 0.0;

  // convert J to world frame
  Matrix3 R = _F.get_rotation();
  Matrix3 J0 = R * _J * Matrix3::transpose(R);

  Real le = _xd.norm_sq() * _mass;
  Real ae = _omega.dot(J0 * _omega);
  Real ke = (le + ae)*0.5;
//  assert(ke > -NEAR_ZERO);
  return ke;
}

/// Gets the number of generalized coordinates
unsigned RigidBody::num_generalized_coordinates(DynamicBody::GeneralizedCoordinateType gctype) const
{
  const unsigned NGC_ROD = 7, NGC_AA = 6;

  // no generalized coordinates if this body is disabled
  if (!_enabled)
    return 0;

  // look for other case where # of g.c.'s is zero: reduced coordinate 
  // articulated body with fixed base and rigid body is the base
  ArticulatedBodyPtr abody = get_articulated_body();
  RCArticulatedBodyPtr rcab = dynamic_pointer_cast<RCArticulatedBody>(abody);
  assert(!(rcab && rcab->get_base_link() != get_this()));
  if (rcab && !rcab->is_floating_base())
    return 0;

  // return the proper number of coordinates
  switch (gctype)
  {
    case DynamicBody::eRodrigues:
      return NGC_ROD;

    case DynamicBody::eAxisAngle:
      return NGC_AA;

    default:
      assert(false);
  }

  // make compiler happy
  assert(false);
  return 0;
}

/// Gets the first parent link of this link; returns NULL if there is no parent
RigidBodyPtr RigidBody::get_parent_link() const
{
  if (_inner_joints.size() > 1)
    throw std::runtime_error("Called RigidBody::get_parent_link() when multiple parent links present! It's not reasonable to call this method for links in maximal-coordinate articulated bodies."); 

  // special case (no parent!)
  if (_inner_joints.empty())
    return RigidBodyPtr(); 

  return RigidBodyPtr(_inner_joints.front().parent);
}

/// Gets the implicit inner joint of this link; returns NULL if there is no implicit inner joint
/**
 * Throws an exception if this link has multiple implicit inner joints
 */
JointPtr RigidBody::get_inner_joint_implicit() const
{
  JointPtr ij;
  BOOST_FOREACH(const InnerJointData& ijd, _inner_joints)
  {
    JointPtr joint(ijd.inner_joint);
    if (joint->get_constraint_type() == Joint::eImplicit)
    {
      if (ij)
        throw std::runtime_error("Multiple implicit joints detected for a single link!"); 
      else
        ij = joint;
    }
  }

  return ij; 
}

/// Updates the velocity of this body using impulses computed via event data
void RigidBody::update_velocity(const EventProblemData& q)
{
  // check for easy exit
  if (q.N_CONTACTS == 0 || !_enabled)
    return;

  // setup summed impulses
  Vector3 j = ZEROS_3, k = ZEROS_3;

  FILE_LOG(LOG_CONTACT) << "RigidBody::update_velocity() entered" << std::endl;

  // update velocity using contact impulses
  for (unsigned i=0, s=0; i< q.contact_events.size(); i++)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1.get() != this && sb2.get() != this)
    {
      s += 2;
      continue;
    }

    FILE_LOG(LOG_CONTACT) << "  contact impulse magnitudes: " << q.alpha_c[i] << " " << q.beta_c[s] << " " << q.beta_c[s+1] << std::endl;

    // check whether to negate normal
    bool negate = (sb2.get() == this);

    // prepare for computation
    Vector3 n = q.contact_events[i]->contact_normal;
    if (negate)
      n = -n;
    const Vector3& p = q.contact_events[i]->contact_point;
    Vector3 r = p - _x;
    Vector3 cross = Vector3::cross(r, n);

    // update j and k
    j += n*q.alpha_c[i];
    k += cross*q.alpha_c[i];

    // now do tangent impulses
    n = q.contact_events[i]->contact_tan1;
    if (negate)
      n = -n;
    cross = Vector3::cross(r, n);
    j += n*q.beta_c[s];
    k += cross*q.beta_c[s];
    s++;
    n = q.contact_events[i]->contact_tan2;
    if (negate)
      n = -n;
    cross = Vector3::cross(r, n);
    j += n*q.beta_c[s];
    k += cross*q.beta_c[s];
    s++;
  }

  // look for whether we update
  if (j.norm() < NEAR_ZERO && k.norm() < NEAR_ZERO)
    return;
  else
    apply_impulse(j, k, _x);

  FILE_LOG(LOG_CONTACT) << "  applying impulses: " << j << " / " << k << std::endl;
  FILE_LOG(LOG_CONTACT) << "  new linear velocity for " << id << ": " << _xd << std::endl;
  FILE_LOG(LOG_CONTACT) << "  new angular velocity for " << id << ": " << _omega << std::endl;
  for (unsigned i=0; i< q.contact_events.size(); i++)
    FILE_LOG(LOG_CONTACT) << "  contact velocity at " << q.contact_events[i]->contact_point << " along normal " << q.contact_events[i]->contact_normal << ": " << q.contact_events[i]->contact_normal.dot(calc_point_vel(q.contact_events[i]->contact_point)) << std::endl;
  FILE_LOG(LOG_CONTACT) << "RigidBody::update_velocity() exited" << std::endl;
}

/// Adds contributions to the event matrices
void RigidBody::update_event_data(EventProblemData& q) 
{
  if (q.N_CONTACTS == 0 || !_enabled)
    return;

  // get inertia matrix in global frame
  Matrix3 R(&_q);
  Matrix3 invJ = R * _invJ * Matrix3::transpose(R);

  // NOTE: b/c this is an individual rigid body, we don't touch the constraint
  // or limit matrices or Ji

  // 1. update Jc_iM_JcT and Jc_v
  for (unsigned i=0; i< q.contact_events.size(); i++)
  {
    // verify that it is the proper type
    assert(q.contact_events[i]->event_type == Event::eContact);

    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
      continue;

    // check whether to negate normal
    bool negate1 = (sb2 == get_this());

    // prepare for computation
    Vector3 n1 = q.contact_events[i]->contact_normal;
    const Vector3& p1 = q.contact_events[i]->contact_point;
    Vector3 r1 = p1 - _x;
    Vector3 cross1 = Vector3::cross(r1, n1);

    // update Jc_v
    if (!negate1)
      q.Jc_v[i] += n1.dot(_xd + Vector3::cross(_omega, r1));
    else
      q.Jc_v[i] -= n1.dot(_xd + Vector3::cross(_omega, r1));

    // scale n1 and cross1
    n1 *= _inv_mass;
    cross1 = invJ * cross1;

    // loop again, note: the matrices are symmetric
    for (unsigned j=i; j< q.contact_events.size(); j++)
    {
      // get the two bodies of the contact
      SingleBodyPtr sb1 = q.contact_events[j]->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = q.contact_events[j]->contact_geom2->get_single_body();

      // if neither sb1 nor sb2 matches this, continue looping
      if (sb1 != get_this() && sb2 != get_this())
        continue;

      // check whether to negate normal
      bool negate2 = (sb2 == get_this());

      // prepare for computation
      const Vector3& n2 = q.contact_events[j]->contact_normal;
      const Vector3& p2 = q.contact_events[j]->contact_point;
      Vector3 r2 = p2 - _x;
      Vector3 cross2 = Vector3::cross(r2, n2);

      // compute entry of Jc_iM_JcT
      Real sum = n1.dot(n2) + cross1.dot(cross2);
      if ((negate1 && !negate2) || (negate2 && !negate1))
        sum = -sum;
      q.Jc_iM_JcT(i, j) += sum; 
      q.Jc_iM_JcT(j, i) = q.Jc_iM_JcT(i, j);
    }
  }

  // 2. update Dc_v
  for (unsigned i=0, ii=0; i< q.contact_events.size(); i++)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
    {
      ii += 2;
      continue;
    }

    // check whether to negate normal
    bool negate = (sb2 == get_this());

    // prepare for computation
    const Vector3& p = q.contact_events[i]->contact_point;
    Vector3 r = p - _x;

    // get the contact tangents
    const Vector3& d1 = q.contact_events[i]->contact_tan1;
    const Vector3& d2 = q.contact_events[i]->contact_tan2;

    // determine velocity vector
    Vector3 vec =  _xd + Vector3::cross(_omega, r);

    // update Dc_v
    if (!negate)
    {
      q.Dc_v[ii++] += d1.dot(vec);
      q.Dc_v[ii++] += d2.dot(vec);
    }
    else
    {
      q.Dc_v[ii++] -= d1.dot(vec);
      q.Dc_v[ii++] -= d2.dot(vec);
    }
  } 

  // 3. update Dc_iM_JcT
  for (unsigned i=0; i< q.contact_events.size(); i++)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
      continue;

    // check whether to negate normal
    bool negate1 = (sb2 == get_this());

    // prepare for computation
    Vector3 n1 = q.contact_events[i]->contact_normal;
    const Vector3& p1 = q.contact_events[i]->contact_point;
    Vector3 r1 = p1 - _x;
    Vector3 cross1 = Vector3::cross(r1, n1);

    // scale n1 and cross1
    n1 *= _inv_mass;
    cross1 = invJ * cross1;

    // loop over all contacts 
    for (unsigned j=0, jj=0; j< q.contact_events.size(); j++)
    {
      // get the two bodies of the contact
      SingleBodyPtr sb1 = q.contact_events[j]->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = q.contact_events[j]->contact_geom2->get_single_body();

      // if neither sb1 nor sb2 matches this, continue looping
      if (sb1 != get_this() && sb2 != get_this())
      {
        jj += 2;
        continue;
      }

      // check whether to negate normal
      bool negate2 = (sb2 == get_this());
      bool negate = ((negate1 && !negate2) || (negate2 && !negate1));

      // prepare for computation
      const Vector3& p2 = q.contact_events[j]->contact_point;
      Vector3 r2 = p2 - _x;

      // compute cross products for both tangent directions
      const Vector3& d21 = q.contact_events[j]->contact_tan1;
      const Vector3& d22 = q.contact_events[j]->contact_tan2;
      Vector3 cross21 = Vector3::cross(r2, d21);
      Vector3 cross22 = Vector3::cross(r2, d22);

      // compute entries of Jc_iM_DcT
      Real sum1 = n1.dot(d21) + cross1.dot(cross21);
      Real sum2 = n1.dot(d22) + cross1.dot(cross22);
      if (!negate)
      {
        q.Jc_iM_DcT(i,jj++) += sum1;
        q.Jc_iM_DcT(i,jj++) += sum2;
      }
      else
      {
        q.Jc_iM_DcT(i,jj++) -= sum1;
        q.Jc_iM_DcT(i,jj++) -= sum2;
      }
    }
  }
 
  // 4. update Dc_iM_DcT
  for (unsigned i=0, ii=0; i< q.contact_events.size(); i++, ii+=2)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
      continue;

    // check whether to negate
    bool negate1 = (sb2 == get_this());

    // get the moment arm 
    Vector3 r1 = q.contact_events[i]->contact_point - _x;

    // get the contact tangents
    Vector3 d1a = q.contact_events[i]->contact_tan1;
    Vector3 d1b = q.contact_events[i]->contact_tan2;

    // compute the cross products
    Vector3 cross1a = Vector3::cross(r1, d1a);
    Vector3 cross1b = Vector3::cross(r1, d1b);

    // scale
    d1a *= _inv_mass;
    d1b *= _inv_mass;
    cross1a = invJ * cross1a;
    cross1b = invJ * cross1b;

    // loop over the remaining contacts 
    for (unsigned j=i, jj=ii; j< q.contact_events.size(); j++, jj+= 2)
    {
      // get the two bodies of the contact
      SingleBodyPtr sb1 = q.contact_events[j]->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = q.contact_events[j]->contact_geom2->get_single_body();

      // if neither sb1 nor sb2 matches this, continue looping
      if (sb1 != get_this() && sb2 != get_this())
        continue;

      // check whether to negate normal
      bool negate2 = (sb2 == get_this());
      bool negate = ((negate1 && !negate2) || (negate2 && !negate1));

      // compute the second moment arm
      Vector3 r2 = q.contact_events[j]->contact_point - _x;

      // get the two tangent vectors 
      const Vector3& d2a = q.contact_events[j]->contact_tan1;
      const Vector3& d2b = q.contact_events[j]->contact_tan2;

      // compute the cross products
      Vector3 cross2a = Vector3::cross(r2, d2a);
      Vector3 cross2b = Vector3::cross(r2, d2b);

      // compute entries of Dc_iM_DcT
      Real sum1 = d1a.dot(d2a) + cross1a.dot(cross2a);
      Real sum2 = d1a.dot(d2b) + cross1a.dot(cross2b);
      Real sum3 = d1b.dot(d2a) + cross1b.dot(cross2a);
      Real sum4 = d1b.dot(d2b) + cross1b.dot(cross2b);

      if (!negate)
      {
        q.Dc_iM_DcT(ii, jj)     += sum1; 
        q.Dc_iM_DcT(ii, jj+1)   += sum2; 
        q.Dc_iM_DcT(ii+1, jj)   += sum3; 
        q.Dc_iM_DcT(ii+1, jj+1) += sum4;
      }
      else
      {
        q.Dc_iM_DcT(ii, jj)     -= sum1; 
        q.Dc_iM_DcT(ii, jj+1)   -= sum2; 
        q.Dc_iM_DcT(ii+1, jj)   -= sum3; 
        q.Dc_iM_DcT(ii+1, jj+1) -= sum4;
      }
    }
  }

  // ensure symmetry on Dc_iM_DcT
  for (unsigned i=0, ii=0; i< q.contact_events.size(); i++, ii+=2)
  {
    // get the two bodies of the contact
    SingleBodyPtr sb1 = q.contact_events[i]->contact_geom1->get_single_body();
    SingleBodyPtr sb2 = q.contact_events[i]->contact_geom2->get_single_body();

    // if neither sb1 nor sb2 matches this, continue looping
    if (sb1 != get_this() && sb2 != get_this())
      continue;

    // loop over the remaining contacts 
    for (unsigned j=i, jj=ii; j< q.contact_events.size(); j++, jj+= 2)
    {
      // get the two bodies of the contact
      SingleBodyPtr sb1 = q.contact_events[j]->contact_geom1->get_single_body();
      SingleBodyPtr sb2 = q.contact_events[j]->contact_geom2->get_single_body();

      // if neither sb1 nor sb2 matches this, continue looping
      if (sb1 != get_this() && sb2 != get_this())
        continue;

      // enforce symmetry
      q.Dc_iM_DcT(jj,ii) = q.Dc_iM_DcT(ii,jj);
      q.Dc_iM_DcT(jj+1,ii) = q.Dc_iM_DcT(ii,jj+1);
      q.Dc_iM_DcT(jj,ii+1) = q.Dc_iM_DcT(ii+1,jj);
      q.Dc_iM_DcT(jj+1,ii+1) = q.Dc_iM_DcT(ii+1,jj+1);
    }
  } 
}

/// Determines whether this link is a "ground" (fixed link)
bool RigidBody::is_ground() const
{
  // clear easy cases
  if (!_enabled)
    return true;

  // can't be a ground if not disabled and not part of an articulated body
  if (_abody.expired())
    return false;

  // now, case will differ depending on what type of articulated body this is
  ArticulatedBodyPtr ab(_abody);
  RCArticulatedBodyPtr rcab = dynamic_pointer_cast<RCArticulatedBody>(ab);
  if (rcab)
  {
    // check whether inner implicit joints are present (if none are present, 
    // this is a base link)
    bool is_base = true;
    BOOST_FOREACH(const InnerJointData& ijd, _inner_joints)
    {
      JointPtr joint(ijd.inner_joint);
      if (joint->get_constraint_type() == Joint::eImplicit)
      {
        is_base = false;
        break;
      }
    }  

    // if this link is a base and the base is fixed, it is a ground
    if (is_base && !rcab->is_floating_base())
      return true;
  }

  // still here? can't be a ground link
  return false;
}

/// Determines whether this link is the base
bool RigidBody::is_base() const
{
  // clear easy cases
  if (_abody.expired())
    return true;

  // check whether no implicit joints are present
  BOOST_FOREACH(const InnerJointData& ijd, _inner_joints)
  {
    JointPtr joint(ijd.inner_joint);
    if (joint->get_constraint_type() == Joint::eImplicit)
      return false;
  }

  // no implicit joints... it's the base
  return true;
}

/// Invalidates an articulated body position state
void RigidBody::invalidate_position()
{
  // only relevant for articulated bodies...
  if (_abody.expired())
    return;

  // get the articulated body and invalidate it
  ArticulatedBodyPtr abody(_abody);
  abody->invalidate_positions();
}

/// Invalidates an articulated body velocity state
void RigidBody::invalidate_velocity()
{
  // only relevant for articulated bodies...
  if (_abody.expired())
    return;

  // get the articulated body and invalidate it
  ArticulatedBodyPtr abody(_abody);
  abody->invalidate_velocities();
}

/// Outputs the object state to the specified stream
/**
 * This method outputs all of the low-level details to the stream
 */
std::ostream& Moby::operator<<(std::ostream& out, const Moby::RigidBody& rb)
{
  // write the ID of the object
  out << "rigid body ID: " << rb.id << endl;

  // indicate whether the body is enabled
  out << "  enabled? " << rb.is_enabled() << endl;

  // write inertia info
  out << "  mass: " << rb.get_mass() << endl;
  out << "  inertia: " << endl << rb.get_inertia();

  // write positions, velocities, and accelerations
  out << "  position: " << rb.get_position() << endl;
  out << "  orientation: " << AAngle(&rb.get_transform()) << endl;
  out << "  linear velocity: " << rb.get_lvel() << endl;
  out << "  angular velocity: " << rb.get_avel() << endl;

  // write sum of forces
  out << "  sum of forces: " << rb.sum_forces() << endl;
  out << "  sum of torques: " << rb.sum_torques() << endl;

  // write the articulated body
  ArticulatedBodyPtr ab = rb.get_articulated_body();
  out << "  articulated body: " << ab << endl;

  // write all collision geometries
  out << "  collision geometries: "; 
  BOOST_FOREACH(CollisionGeometryPtr g, rb.geometries)
    out << "    " << g << endl;

  return out;
}

