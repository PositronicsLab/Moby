/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 * ****************************************************************************/

#ifdef USE_OSG
#include <osg/MatrixTransform>
#include <Moby/OSGGroupWrapper.h>
#endif

#include <iostream>
#include <Moby/XMLTree.h>
#include <Moby/Primitive.h>
#include <Moby/Visualizable.h>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;

Visualizable::Visualizable()
{
  // set the visualization relative pose
  _vF = shared_ptr<Pose3d>(new Pose3d);

  #ifdef USE_OSG
  // create the OSGGroupWrapper
  _vizdata = OSGGroupWrapperPtr(new OSGGroupWrapper);

  // create the transform 
  _group = new osg::MatrixTransform;

  // ref the group 
  _group->ref();

  // NOTE: the non-const cast below is potentially dangerous, but I don't
  // think it will be a problem
//  _group->addChild((osg::Group*) _vizdata->get_group());
  _group->addChild(_vizdata->get_group());
  #endif
}

Visualizable::~Visualizable()
{
  #ifdef USE_OSG
  _group->unref();
  #endif
}

/// Sets the visualization relative pose
void Visualizable::set_visualization_relative_pose(const Pose3d& P)
{
  // verify the P's relative pose is correct 
  if (P.rpose != _vF->rpose)
    throw std::runtime_error("Visualizable::set_visualization_relative_pose() - pose not relative to the proper pose");

  // set the pose
  *_vF = P;
}

/// Sets the visualization data from a OSGGroupWrapper
void Visualizable::set_visualization_data(OSGGroupWrapperPtr vdata)
{
  #ifdef USE_OSG
  // store the OSGGroupWrapper
  _vizdata = vdata;

  // get the group 
  osg::Group* vgroup = _vizdata->get_group();

  // replace and add the child
  _group->removeChildren(0, _group->getNumChildren());
  _group->addChild(vgroup);
  #endif
}

/// Sets the visualization data from a node
void Visualizable::set_visualization_data(osg::Node* vdata)
{
  #ifdef USE_OSG
  // create a new OSGGroupWrapper using the vdata
  _vizdata = OSGGroupWrapperPtr(new OSGGroupWrapper(vdata));

  // clear all children
  _group->removeChildren(0, _group->getNumChildren());

  // add the separator from the OSGGroupWrapper to this separator
  _group->addChild((osg::Node*) _vizdata->get_group());
  #endif
}

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

/// Updates the visualization using the appropriate transform
/**
 * \note in contrast to set_visualization_transform(), this method
 * actually sets the transform in the scenegraph, thereby updating
 * the display.
 * \note derived classes may override this method to set the transform
 * right before it is updated, if desired
 * \sa set_visualization_transform()
 */
void Visualizable::update_visualization()
{
  #ifdef USE_OSG
  // if there is no visualization data, quit now
  if (!_vizdata)
    return;

  // get the transform; if there is none, quit now
  shared_ptr<const Pose3d> T = get_visualization_pose();
  if (!T)
    return;

  // convert the pose to reference the global frame
  Pose3d T0 = *T;
  T0.update_relative_pose(GLOBAL);

  // update the transform
  osg::Matrixd m;
  to_osg_matrix(T0, m);
  _group->setMatrix(m);
  #endif
}

/// Gets the visualization data for this object
osg::Group* Visualizable::get_visualization_data() const
{
  #ifdef USE_OSG
  return (osg::Group*) _group;
  #else
  return NULL;
  #endif
}

/// Utility method for load_from_xml()
/**
 * This method searches for visualization-filename,
 * visualization-id attributes for
 * a given node and creates a separator based on the attribute found.
 */
osg::Group* Visualizable::construct_from_node(shared_ptr<const XMLTree> node, const std::map<std::string, BasePtr>& id_map)
{  
  std::map<std::string, BasePtr>::const_iterator id_iter; 
  osg::Group* group = NULL;

  // get the relevant attributes
  XMLAttrib* viz_id_attr = node->get_attrib("visualization-id");
  XMLAttrib* vfile_id_attr = node->get_attrib("visualization-filename");

  // check that some visualization data exists
  if (!viz_id_attr && !vfile_id_attr)
    return NULL;

  // check for mix-and-match
  if (viz_id_attr && vfile_id_attr)
  {
    std::cerr << "Visualizable::construct_from_xml() - found mix and match of ";
    std::cerr << "  visualization-id and visualization-filename ";
    std::cerr << "in offending node: " << std::endl << *node;
    return NULL;
  }

  // look for a Primitive or OSGGroup id
  if (viz_id_attr)
  {
    // get the id
    const std::string& id = viz_id_attr->get_string_value();

    // find the ID object
    if ((id_iter = id_map.find(id)) == id_map.end())
    {
      std::cerr << "Visualizable::construct_from_xml() - cannot find ";
      std::cerr << "Primitive/OSGGroup object " << std::endl;
      std::cerr << "  with id '" << id << "' in offending node: " << std::endl;
      std::cerr << *node;
      return NULL;
    }

    // look for it as an OSGGroup 
    OSGGroupWrapperPtr wrapper = boost::dynamic_pointer_cast<OSGGroupWrapper>(id_iter->second);  
    if (wrapper)
    {
      // get the group
      group = wrapper->get_group();
    }
    else
    {
      // it should be a Primitive, if it is not a OSGGroup 
      PrimitivePtr prm = boost::dynamic_pointer_cast<Primitive>(id_iter->second);  
      assert(prm);

      // create the group and add a node to it
      #ifdef USE_OSG
      group = new osg::Group;
      group->addChild(prm->get_visualization());
      #endif
    }
  }
  // visualization-filename attribute
  else
  {
    #ifdef USE_OSG
    // sanity check
    assert(vfile_id_attr);

    // get the filename
    const std::string& fname = vfile_id_attr->get_string_value();

    // create the new OSGGroup wrapper
    OSGGroupWrapperPtr wrapper(new OSGGroupWrapper(fname));

    // get the group-- and reference it 
    group = wrapper->get_group();
    group->ref();
    #endif
  }

  // return the group 
  return group;
}

/// Implements Base::load_from_xml() 
void Visualizable::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  OSGGroupWrapperPtr wrap;

  // load parent data
  Base::load_from_xml(node, id_map);

  // get the relevant attributes
  XMLAttrib* viz_id_attr = node->get_attrib("visualization-id");
  XMLAttrib* vfile_id_attr = node->get_attrib("visualization-filename");

  // get whether there are any visualization nodes
  std::list<shared_ptr<const XMLTree> > viz_nodes = node->find_child_nodes("Visualization");

  // check that some visualization data exists
  if (!viz_id_attr && !vfile_id_attr && viz_nodes.empty())
    return;

  // check for mix-and-match
  if ((viz_id_attr && vfile_id_attr) || (viz_id_attr && vfile_id_attr) ||
      (vfile_id_attr && !viz_nodes.empty()))
  {
    std::cerr << "Visualizable::load_from_xml() - found mix and match of ";
    std::cerr << "visualization-id, visualization-filename, " << std::endl;
    std::cerr << "  and/or Visualization node in offending node: " << std::endl;
    std::cerr << *node;
    return;
  }

  #ifdef USE_OSG
  // if one of the attributes has been set, get the group
  if (viz_id_attr || vfile_id_attr)
  {
    osg::Group* group = construct_from_node(node, id_map);
    if (group)
      set_visualization_data(group);
  }
  else
  {
    // if we are here, then one or more Visualization nodes encountered
    // create a new group to hold all data read..
    osg::Group* root = NULL;

    // handle all Visualization nodes
    for (std::list<shared_ptr<const XMLTree> >::const_iterator i = viz_nodes.begin(); i != viz_nodes.end(); i++)
    {
      // get the group from the child node
      osg::Group* child_group = construct_from_node(*i, id_map);

      // if the child group is NULL, continue looping
      if (!child_group)
        continue;

      // look for a visualization-rel-x attributes
      XMLAttrib* rel_origin_attr = (*i)->get_attrib("visualization-rel-origin");
      XMLAttrib* rel_rpy_attr = (*i)->get_attrib("visualization-rel-rpy");
      XMLAttrib* rel_quat_attr = (*i)->get_attrib("visualization-rel-quat");
      if (rel_origin_attr || rel_rpy_attr || rel_quat_attr)
      {
        // create a new transform group
        osg::MatrixTransform* xgroup = new osg::MatrixTransform;

        // create the transform and set it
        Pose3d Tx;
        if (rel_origin_attr)
          Tx.x = rel_origin_attr->get_origin_value();
        if (rel_rpy_attr)
          Tx.q = rel_rpy_attr->get_rpy_value();
        else if (rel_quat_attr)
          Tx.q = rel_quat_attr->get_quat_value();
        osg::Matrixd T;
        to_osg_matrix(Tx, T);
        xgroup->setMatrix(T);

        // add the child to the transform group
        xgroup->addChild(child_group); 

        // point child_sep to new_sep (its parent)
        child_group = xgroup;
      }

      // create the root, if need be
      if (!root)
        root = new osg::Group;

      // add the child group to the root group 
      root->addChild(child_group);
    }

    // set the visualization, if available
    if (root)
      set_visualization_data(root);
  }
  #endif

  // update the visualization transform
  update_visualization();
}

/// Implements Base::save_to_xml() 
void Visualizable::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the Base data
  Base::save_to_xml(node, shared_objects);

  // rename this node
  node->name = "Visualizable";

  #ifdef USE_OSG
  // if there is no visualization data, skip..
  if (!_vizdata)
    return;

  // save the OSGGroup id
  node->attribs.insert(XMLAttrib("visualization-id", _vizdata->id));

  // add the visualization data to the list of shared objects
  shared_objects.push_back(_vizdata);
  #endif
} 

