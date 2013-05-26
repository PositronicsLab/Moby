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
  _F = shared_ptr<Pose3d>(new Pose3d); 
  _jF = shared_ptr<Pose3d>(new Pose3d);
  _jF->rpose = _F;
  _J.pose = _jF;
  _intersection_tolerance = 1e-5;
  _deformable = false;

  // set visualization members to NULL
  _vtransform = NULL;
}

/// Constructs a primitive with the specified transform
Primitive::Primitive(shared_ptr<const Pose3d> F)
{
  _F = shared_ptr<Pose3d>(new Pose3d);
  _jF = shared_ptr<Pose3d>(new Pose3d);
  _jF->rpose = _F;
  _J.pose = _jF;
  *_F = *F; 
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

  return (osg::Node*) _vtransform;
}

/// Updates the visualization on the primitive 
void Primitive::update_visualization()
{
  #ifdef USE_OSG
  // if there is no visualization, exit now
  if (!_vtransform)
    return;

  // convert the pose to be relative to the global frame
  Pose3d F0 = *_F;
  F0.update_relative_pose(GLOBAL);

  // update the transform
  osg::Matrixd T;
  to_osg_matrix(F0, T);
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
  _J.m = mass;
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
    if (mass_attr->get_real_value() < 0.0)
      throw std::runtime_error("Attempting to set primitive mass to negative value");
    set_mass(mass_attr->get_real_value());
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
  shared_ptr<Pose3d> F(new Pose3d);
  const XMLAttrib* xlat_attr = node->get_attrib("position");
  const XMLAttrib* rpy_attr = node->get_attrib("rpy");
  const XMLAttrib* quat_attr = node->get_attrib("quat");
  if (xlat_attr && rpy_attr)
  {
    F->x = xlat_attr->get_origin_value();
    F->q = rpy_attr->get_rpy_value();
    set_pose(F);
  }
  else if (xlat_attr && quat_attr)
  {
    F->x = xlat_attr->get_origin_value();
    F->q = quat_attr->get_quat_value();
    set_pose(F);
  }
  else if (xlat_attr)
  {
    F->x = xlat_attr->get_origin_value();
    set_pose(F);
  }
  else if (rpy_attr)
  {
    F->q = rpy_attr->get_rpy_value();
    set_pose(F);
  }
  else if (quat_attr)
  {
    F->q = quat_attr->get_quat_value();
    set_pose(F);
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
    node->attribs.insert(XMLAttrib("mass", _J.m));

  // save the transform for the primitive
  Pose3d F0 = *_F;
  F0.update_relative_pose(GLOBAL);
  node->attribs.insert(XMLAttrib("position", F0.x));
  node->attribs.insert(XMLAttrib("quat", F0.q));

  // add the intersection tolerance as an attribute
  node->attribs.insert(XMLAttrib("intersection-tolerance", _intersection_tolerance));
}

/// Sets the transform for this primitive -- transforms mesh and inertial properties (if calculated)
void Primitive::set_pose(shared_ptr<const Pose3d> F)
{
  // save the new transform
  *_F = *F;

  #ifdef USE_OSG
  if (_vtransform)
  {
    // get the pose in the global frame
    Pose3d F0 = *_F;
    F0.update_relative_pose(GLOBAL);

    // update the visualization transform
    osg::Matrixd tgt;
    to_osg_matrix(F0, tgt);
    _vtransform->setMatrix(tgt);
  }
  #endif
}

