/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iostream>
#include <Inventor/nodes/SoTransform.h>
#include <Inventor/actions/SoWriteAction.h>
#include <Moby/Matrix4.h>
#include <Moby/XMLTree.h>
#include <Moby/InvalidTransformException.h>
#include <Moby/SoSeparatorWrapper.h>

using namespace Moby;

SoSeparatorWrapper::SoSeparatorWrapper()
{
  _separator = new SoSeparator;
  _separator->ref();
}

SoSeparatorWrapper::SoSeparatorWrapper(SoNode* n)
{
  _separator = new SoSeparator;
  _separator->addChild(n);
  _separator->ref();
}

/// Creates an SoSeparator wrapper given a OpenInventor or VRML filename
SoSeparatorWrapper::SoSeparatorWrapper(const std::string& fname)
{
  // open the filename and read in the file
  SoInput in;
  SbBool ok = in.openFile(fname.c_str());
  if (!ok || !(_separator = SoDB::readAll(&in)))
  {
    std::cerr << "SoSeparatorWrapper::SoSeparatorWrapper() - unable to read ";
    std::cerr << "from " << fname << "!" << std::endl;
    return;
  }

  // reference the new separator
  _separator->ref();
}

SoSeparatorWrapper::~SoSeparatorWrapper()
{
  _separator->unref();
}

/// Implements Base::load_from_xml()
void SoSeparatorWrapper::load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map)
{
  // load the Base data
  Base::load_from_xml(node, id_map);

  // verify the node name
  assert(strcasecmp(node->name.c_str(), "SoSeparator") == 0);

  // if there is no visualization data filename, return now
  const XMLAttrib* viz_fname_attr = node->get_attrib("filename");
  if (!viz_fname_attr)
    return;

  // get the filename
  const std::string& fname = viz_fname_attr->get_string_value();

  // open the filename and read in the file
  SoSeparator* sep;
  SoInput in;
  SbBool ok = in.openFile(fname.c_str());
  if (!ok || !(sep = SoDB::readAll(&in)))
  {
    std::cerr << "SoSeparatorWrapper::load_from_xml() - unable to read from ";
    std::cerr  << fname << "!" << std::endl;
    return;
  }

  // remove all children from the root separator
  _separator->removeAllChildren();

  // read in the transform, if specified
  const XMLAttrib* transform_attr = node->get_attrib("transform");
  if (transform_attr)
  {
    Matrix4 T;
    transform_attr->get_matrix_value(T);
    if (!Matrix4::valid_transform(T))
      throw InvalidTransformException(T);

    // create the SoTransform and add it to the root separator
    SoTransform* trans = new SoTransform;
    SbMatrix m;
    Matrix4::to_inventor_matrix(T, &m);
    trans->setMatrix(m);    
    _separator->addChild(trans);
  }

  // add the new separator to the root separator 
  _separator->addChild(sep); 
}

/// Implements Base::save_to_xml()
void SoSeparatorWrapper::save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const
{
  // save the Base data
  Base::save_to_xml(node, shared_objects);

  // rename this node
  node->name = "SoSeparator";

  // form the filename using the ID 
  std::string filename = "vizdata_" + id + ".iv";

  // save the visualization data 
  node->attribs.insert(XMLAttrib("filename", filename));
  SoWriteAction write_action;
  write_action.getOutput()->setBinary(false);
  write_action.getOutput()->openFile(filename.c_str());
  write_action.apply(_separator);
  write_action.getOutput()->closeFile();
}

