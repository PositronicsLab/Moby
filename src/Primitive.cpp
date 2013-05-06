/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifdef USE_OSG
#include <osg/MatrixTransform>
#include <osg/Material>
#include <osg/Matrixd>
#endif
#include <queue>
#include <Moby/Constants.h>
#include <Moby/XMLTree.h>
#include <Moby/Primitive.h>

using boost::shared_ptr;
using namespace Ravelin;
using namespace Moby;
using std::endl;
using std::list;
using std::vector;
using std::queue;
using std::make_pair;
using std::pair;
using std::map;
using std::cerr;
using boost::dynamic_pointer_cast;

#ifdef USE_OSG
/// Copies this matrix to an OpenSceneGraph Matrixd object
static void to_osg_matrix(const Pose3d& src, osg::Matrixd& tgt)
{
  // get the rotation matrix
  Matrix3d M = src.q;

  // setup the rotation components of tgt
  const unsigned X = 0, Y = 1, Z = 2, W = 3;
  for (unsigned i=X; i<= Z; i++)
    for (unsigned j=X; j<= Z; j++)
      tgt(j,i) = M(i,j);

  // setup the translation components of tgt
  for (unsigned i=X; i<= Z; i++)
    tgt(W,i) = src.x[i];

  // set constant values of the matrix
  tgt(X,W) = tgt(Y,W) = tgt(Z,W) = (double) 0.0;
  tgt(W,W) = (double) 1.0;
}
#endif

/// Constructs a primitive under the identity transformation
Primitive::Primitive()
{
  _mass = 0.0;
  _J = ZEROS_3x3;
  _com = ZEROS_3;
  _T = IDENTITY_4x4; 
  _intersection_tolerance = 1e-5;
  _deformable = false;

  // set visualization members to NULL
  _vtransform = NULL;
}

/// Constructs a primitive with the specified transform
Primitive::Primitive(const Pose3d& T)
{
  _mass = 0.0;
  _J = ZEROS_3x3;
  _com = ZEROS_3;
  _T = T;
  _intersection_tolerance = 1e-5;
  _deformable = false;

  // set visualization members to NULL
  _vtransform = NULL;
}

Primitive::~Primitive()
{
  #ifdef USE_OSG
  if (_vtransform)
    _vtransform->unref();
  #endif
}

/// Sets the intersection tolerance for this primitive
void Primitive::set_intersection_tolerance(double tol)
{
  if (tol < (double) 0.0)
    throw std::runtime_error("Primitive::set_intersection_tolerance() - tolerance cannot be negative!");

  _intersection_tolerance = tol;
}

/// Gets the visualization for this primitive, creating it if necessary
osg::Node* Primitive::get_visualization()
{
  #ifdef USE_OSG
  // if the visualization is already created, return it
  if (_vtransform)
    return _vtransform;

  // otherwise, create the visualization objects
  _vtransform = new osg::MatrixTransform;

  // setup a random color
  const float RED = (float) rand() / RAND_MAX; 
  const float GREEN = (float) rand() / RAND_MAX; 
  const float BLUE = (float) rand() / RAND_MAX; 

  // create a new material -- random color
  _mat = new osg::Material;
  _mat->setColorMode(osg::Material::DIFFUSE);
  _mat->setDiffuse(osg::Material::FRONT, osg::Vec4(RED, GREEN, BLUE, 1.0f)); 
  _vtransform->getOrCreateStateSet()->setAttribute(_mat);

  // reference the objects
  _vtransform->ref();
  _mat->ref();

  // update the visualizaiton 
  update_visualization();  
  #endif

  return _vtransform;
}

/// Updates the visualization on the primitive 
void Primitive::update_visualization()
{
  #ifdef USE_OSG
  // if there is no visualization, exit now
  if (!_vtransform)
    return;

  // update the transform
  osg::Matrixd T;
  to_osg_matrix(_T, T);
  _vtransform->setMatrix(T);

  // remove all children of the group 
  _vtransform->removeChildren(0, _vtransform->getNumChildren());

  // add the visualization to this separator
  _vtransform->addChild(create_visualization());
  #endif
}

/// Sets the mass of this primitive
/**
 * \note sets the density too
 */
void Primitive::set_mass(double mass)
{
  _mass = mass;
  calc_mass_properties();
}

/// Sets the density of this primitive
void Primitive::set_density(double density)
{
  _density = shared_ptr<double>(new double);
  *_density = density;
  calc_mass_properties();
}

/// Implements Base::load_from_xml() for serialization
void Primitive::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // ******************************************************************
  // do *not* verify the node name, b/c Primitive is abstract 
  // ******************************************************************

  // load the parent data
  Base::load_from_xml(node, id_map);

  // read in the mass, if specified
  const XMLAttrib* mass_attr = node->get_attrib("mass");
  if (mass_attr)
  {
    _mass = mass_attr->get_real_value();
    if (_mass < 0.0)
      throw std::runtime_error("Attempting to set primitive mass to negative value");
  }

  // read in the density if specified
  const XMLAttrib* density_attr = node->get_attrib("density");
  if (density_attr)
  {
    _density = shared_ptr<double>(new double);
    *_density = density_attr->get_real_value();
    if (*_density < (double) 0.0)
      throw std::runtime_error("Attempting to set primitive density to negative value");
  }

  // read the intersection tolerance
  const XMLAttrib* itol_attr = node->get_attrib("intersection-tolerance");
  if (itol_attr)
    set_intersection_tolerance(itol_attr->get_real_value());

  // read in transformation, if specified
  Pose3d T;
  const XMLAttrib* xlat_attr = node->get_attrib("translation");
  const XMLAttrib* rpy_attr = node->get_attrib("rpy");
  const XMLAttrib* quat_attr = node->get_attrib("quat");
  if (xlat_attr && rpy_attr)
  {
    T.x = xlat_attr->get_origin_value();
    T.q = rpy_attr->get_rpy_value();
    set_pose(T);
  }
  else if (xlat_attr && quat_attr)
  {
    T.x = xlat_attr->get_origin_value();
    T.q = quat_attr->get_quat_value();
    set_pose(T);
  }
  else if (xlat_attr)
  {
    T.x = xlat_attr->get_origin_value();
    set_pose(T);
  }
  else if (rpy_attr)
  {
    T.q = rpy_attr->get_rpy_value();
    set_pose(T);
  }
  else if (quat_attr)
  {
    T.q = quat_attr->get_quat_value();
    set_pose(T);
  }

  // calculate mass properties
  calc_mass_properties();
}

/// Implements Base::save_to_xml() for serialization
void Primitive::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data
  Base::save_to_xml(node, shared_objects);

  // set the node name
  node->name = "Primitive";

  // save the density, if specified, else the mass
  if (_density)
    node->attribs.insert(XMLAttrib("density", *_density));
  else
    node->attribs.insert(XMLAttrib("mass", _mass));

  // save the transform for the primitive
  node->attribs.insert(XMLAttrib("translation", _T.x));
  node->attribs.insert(XMLAttrib("quat", _T.q));

  // add the intersection tolerance as an attribute
  node->attribs.insert(XMLAttrib("intersection-tolerance", _intersection_tolerance));
}

/// Transforms the inertia given a rotation matrix
void Primitive::transform_inertia(double mass, const Matrix3d& J_in, const Point3d& com_in, const Matrix3d& R, Matrix3d& J_out, Point3d& com_out)
{  
  const unsigned X = 0, Y = 1, Z = 2;

  // copy com_in -- this is in case com_in and com_out are equal
  Point3d com_in_copy = com_in;

  // rotate/scale the inertia
  J_out = R * J_in.mult_transpose(R);
  
  // determine the new center-of-mass
  com_out = R * com_in_copy;

  // determine r (vector from displacement to com)
  Vector3d r = com_in_copy - com_out;
  
  // displace it using the parallel axis theorem
  J_out(X,X) += mass * (r[Y] * r[Y] + r[Z] * r[Z]);
  J_out(Y,Y) += mass * (r[Z] * r[Z] + r[X] * r[X]);
  J_out(Z,Z) += mass * (r[X] * r[X] + r[Y] * r[Y]);
  J_out(X,Y) -= mass * (r[X] * r[Y]);
  J_out(Y,Z) -= mass * (r[Y] * r[Z]);
  J_out(X,Z) -= mass * (r[Z] * r[X]);
  J_out(Y,X) = J_out(X,Y);
  J_out(Z,Y) = J_out(Y,Z);
  J_out(Z,X) = J_out(X,Z);
}

/// Transforms the inertia given a transformation matrix
void Primitive::transform_inertia(double mass, const Matrix3d& J_in, const Point3d& com_in, const Pose3d& T, Matrix3d& J_out, Point3d& com_out)
{  
  const unsigned X = 0, Y = 1, Z = 2;

  // copy com_in -- this is in case com_in and com_out are equal
  Point3d com_in_copy = com_in;

  // get the rotation/scaling part of the matrix as a 3x3 and the translation
  // part of the matrix as a vector
  Matrix3d RS = T.q;

  // rotate/scale the inertia
  J_out = RS * J_in.mult_transpose(RS);
  
  // determine the new center-of-mass
  com_out = RS * com_in_copy + T.x;

  // determine r (vector from displacement to com)
  Vector3d r = com_in_copy - com_out;
  
  // displace it using the parallel axis theorem
  J_out(X,X) += mass * (r[Y] * r[Y] + r[Z] * r[Z]);
  J_out(Y,Y) += mass * (r[Z] * r[Z] + r[X] * r[X]);
  J_out(Z,Z) += mass * (r[X] * r[X] + r[Y] * r[Y]);
  J_out(X,Y) -= mass * (r[X] * r[Y]);
  J_out(Y,Z) -= mass * (r[Y] * r[Z]);
  J_out(X,Z) -= mass * (r[Z] * r[X]);
  J_out(Y,X) = J_out(X,Y);
  J_out(Z,Y) = J_out(Y,Z);
  J_out(Z,X) = J_out(X,Z);
}

/// Sets the transform for this primitive -- transforms mesh and inertial properties (if calculated)
void Primitive::set_pose(const Pose3d& T)
{
  // save the new transform
  _T = T;

  #ifdef USE_OSG
  if (_vtransform)
  {
    // update the visualization transform
    osg::Matrixd tgt;
    to_osg_matrix(T, tgt);
    _vtransform->setMatrix(tgt);
  }
  #endif
}

/// Loads the state of this primitive
void Primitive::load_state(shared_ptr<void> state)
{
  shared_ptr<PrimitiveState> ps = boost::static_pointer_cast<PrimitiveState>(state);
  _deformable = ps->deformable;
  _J = ps->J;
  _density = ps->density;
  _mass = ps->mass;
  _com = ps->com;
  _T = ps->T;
  _intersection_tolerance = ps->intersection_tolerance;
  _invalidated = true;
}

/// Saves the state of this primitive
shared_ptr<void> Primitive::save_state() const
{
  shared_ptr<PrimitiveState> ps(new PrimitiveState);
  ps->deformable = _deformable;
  ps->J = _J;
  ps->density = _density;
  ps->mass = _mass;
  ps->com = _com;
  ps->T = _T;
  ps->intersection_tolerance = _intersection_tolerance;

  return ps;
}

