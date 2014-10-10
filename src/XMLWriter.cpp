/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <errno.h>
#include <iostream>
#include <fstream>
#include <Moby/Base.h>
#include <Moby/XMLTree.h>
#include <Moby/XMLWriter.h>

using boost::shared_ptr;
using namespace Moby;

/// Serializes the given objects (and all dependencies) to XML
void XMLWriter::serialize_to_xml(const std::string& fname, const std::list<shared_ptr<const Base> >& objects)
{
  // get the filename
  std::string filename = fname;

  // get the current path
  size_t BUFSIZE = 8192;
  boost::shared_array<char> cwd;
  while (true)
  {
    cwd = boost::shared_array<char>((new char[BUFSIZE]));
    if (getcwd(cwd.get(), BUFSIZE) == cwd.get())
      break;
    if (errno != ERANGE)
    {
      std::cerr << "XMLWriter::serialize_to_xml() - unable to allocate sufficient memory!" << std::endl;
      return;
    }
    BUFSIZE *= 2;
  }

  // separate the path from the filename
  size_t last_path_sep = fname.find_last_of('/');
  if (last_path_sep != std::string::npos)
  {
    // get the new working path
    std::string pathname = fname.substr(0,last_path_sep+1);

    // change to the new working path
    chdir(pathname.c_str());

    // get the new filename
    filename = fname.substr(last_path_sep+1,std::string::npos);
  }

  // create a new XMLTree
  XMLTreePtr topnode(new XMLTree("XML"));

  // create a node for Moby
  XMLTreePtr node(new XMLTree("Moby"));
  topnode->add_child(node);

  // setup a list of shared objects to be equal to the list of objects
  std::list<shared_ptr<const Base> > shared_objects = objects;

  // init a set of serialized objects
  std::set<shared_ptr<const Base> > serialized;

  // develop the XML tree until there is nothing more to serialize
  while (!shared_objects.empty())
  {
    // get the object off of the front of the queue
    shared_ptr<const Base> obj = shared_objects.front();
    assert(obj);
    shared_objects.pop_front();

    // if this object has already been serialized, skip it
    if (serialized.find(obj) != serialized.end())
      continue;

    // create a new node for this object under the parent
    XMLTreePtr new_node(new XMLTree(""));
    node->add_child(new_node);

    // serialize to this new node
    obj->save_to_xml(new_node, shared_objects);

    // indicate that the node has been serialized
    serialized.insert(obj);
  }

  // open the file for writing
  std::ofstream out(filename.c_str());

  // write the tree to the file
  out << *topnode << std::endl;

  // close the file
  out.close();

  // change back to the original directory (if any)
  chdir(cwd.get());
}

/// Serializes the given object (and all of its dependencies) to XML
void XMLWriter::serialize_to_xml(const std::string& fname, shared_ptr<const Base> object)
{
  // get the filename
  std::string filename = fname;

  // get the current path
  size_t BUFSIZE = 8192;
  boost::shared_array<char> cwd;
  while (true)
  {
    cwd = boost::shared_array<char>((new char[BUFSIZE]));
    if (getcwd(cwd.get(), BUFSIZE) == cwd.get())
      break;
    if (errno != ERANGE)
    {
      std::cerr << "XMLWriter::serialize_to_xml() - unable to allocate sufficient memory!" << std::endl;
      return;
    }
    BUFSIZE *= 2;
  }

  // separate the path from the filename
  size_t last_path_sep = fname.find_last_of('/');
  if (last_path_sep != std::string::npos)
  {
    // get the new working path
    std::string pathname = fname.substr(0,last_path_sep+1);

    // change to the new working path
    chdir(pathname.c_str());

    // get the new filename
    filename = fname.substr(last_path_sep+1,std::string::npos);
  }

  // create a new XMLTree
  XMLTreePtr topnode(new XMLTree("XML"));

  // create a node for Moby
  XMLTreePtr node(new XMLTree("Moby"));
  topnode->add_child(node);

  // setup a list of shared objects
  std::list<shared_ptr<const Base> > shared_objects;

  // add the desired object to the list of shared objects
  shared_objects.push_back(object);

  // init a set of serialized objects
  std::set<shared_ptr<const Base> > serialized;
  
  // develop the XML tree until there is nothing more to serialize
  while (!shared_objects.empty())
  {
    // get the object off of the front of the queue
    shared_ptr<const Base> obj = shared_objects.front();
    assert(obj);
    shared_objects.pop_front();

    // if this object has already been serialized, skip it
    if (serialized.find(obj) != serialized.end())
      continue;

    // create a new node for this object under the parent
    XMLTreePtr new_node(new XMLTree(""));
    node->add_child(new_node);

    // serialize to this new node
    obj->save_to_xml(new_node, shared_objects);

    // indicate that the node has been serialized
    serialized.insert(obj);
  }

  // open the file for writing
  std::ofstream out(filename.c_str());

  // write the tree to the file
  out << *topnode << std::endl;

  // close the file
  out.close();

  // change back to the initial working directory
  chdir(cwd.get());
}

