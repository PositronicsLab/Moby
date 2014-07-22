/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <string.h>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <sstream>
#include <Ravelin/MatrixNd.h>
#include <Moby/MissizeException.h>
#include <Moby/XMLTree.h>

using boost::shared_ptr;
using namespace Moby;
using namespace Ravelin;

/// Constructs an XMLAttrib object from a name and a string value
XMLAttrib::XMLAttrib(const std::string& name, const std::string& string_value)
{
  this->name = name;
  this->value = string_value;
  this->processed = false;
}

/// Constructs a real-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, double real_value)
{
  this->name = name;
  this->value = str(real_value);
  this->processed = false;
}

/// Constructs an integer-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, int int_value)
{
  this->name = name;
  std::ostringstream oss;
  oss << int_value;
  this->value = oss.str();
  this->processed = false;
}

/// Constructs an unsigned integer-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, unsigned unsigned_value)
{
  this->name = name;
  std::ostringstream oss;
  oss << unsigned_value;
  this->value = oss.str();
  this->processed = false;
}

/// Constructs a Origin3d-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const Origin3d& o)
{
  this->name = name;
  std::ostringstream oss;
  oss << o[0] << " " << o[1] << " " << o[2];
  this->value = oss.str();
  this->processed = false;
}

/// Constructs a Vector2d-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const Vector2d& v)
{
  this->name = name;
  std::ostringstream oss;
  oss << v[0] << " " << v[1];
  this->value = oss.str();
  this->processed = false;
}

/// Constructs a Vector3d-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const Vector3d& v)
{
  this->name = name;
  std::ostringstream oss;
  oss << v[0] << " " << v[1] << " " << v[2];
  this->value = oss.str();
  this->processed = false;
}

/// Constructs a Quatd-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const Quatd& q)
{
  this->name = name;
  std::ostringstream oss;
  oss << q.w << " " << q.x << " " << q.y << " " << q.z;
  this->value = oss.str();
  this->processed = false;
}

/// Constructs a vector-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const VectorNd& vector_value)
{
  this->processed = false;
  this->name = name;
  std::ostringstream oss;

  // if the vector is empty, return prematurely
  if (vector_value.size() == 0)
  {
    this->value = "";
    return;
  }

  // set the first value of the vector
  oss << str(vector_value[0]);

  // set subsequent values of the vector
  for (unsigned i=1; i< vector_value.size(); i++)
      oss << " " << str(vector_value[i]);

  // get the string value
  this->value = oss.str();
}

/// Constructs a vector-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const SVector6d& vector_value)
{
  this->processed = false;
  this->name = name;
  std::ostringstream oss;

  // set the first value of the vector
  oss << str(vector_value[0]);

  // set subsequent values of the vector
  for (unsigned i=1; i< vector_value.size(); i++)
    oss << " " << str(vector_value[i]);

  // get the string value
  this->value = oss.str();
}

/// Constructs a vector-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const MatrixNd& matrix_value)
{  
  this->processed = false;
  this->name = name;
  std::ostringstream oss;

  // if the matrix is empty, return prematurely
  if (matrix_value.rows() == 0 || matrix_value.columns() == 0)
  {
    this->value = "";
    return;
  }

  // set the first value of the matrix
  oss << str(matrix_value(0,0));

  // for each row of the matrix
  for (unsigned j=0; j< matrix_value.rows(); j++)
  {
    // determine column iteration
    unsigned i = (j == 0) ? 1 : 0;

    // for each column of the matrix
    for (; i< matrix_value.columns(); i++)
      oss << " " << str(matrix_value(j,i));

    // separate rows with a semicolon
    oss << ";";
  }

  // get the string value
  this->value = oss.str();
}

/// Constructs a vector-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const Matrix3d& matrix_value)
{  
  this->processed = false;
  this->name = name;
  std::ostringstream oss;

  // set the first value of the matrix
  oss << str(matrix_value(0,0));

  // for each row of the matrix
  for (unsigned j=0; j< matrix_value.rows(); j++)
  {
    // determine column iteration
    unsigned i = (j == 0) ? 1 : 0;

    // for each column of the matrix
    for (; i< matrix_value.columns(); i++)
      oss << " " << str(matrix_value(j,i));

    // separate rows with a semicolon
    oss << ";";
  }

  // get the string value
  this->value = oss.str();
}

/// Constructs a Boolean-valued attribute from the given value
XMLAttrib::XMLAttrib(const std::string& name, bool bool_value)
{
  this->processed = false;
  this->name = name;
  this->value = (bool_value) ? "true" : "false";  
}

/// Constructs a long-valued attribute from the given value
XMLAttrib::XMLAttrib(const std::string& name, long long_value)
{
  this->processed = false;
  this->name = name;
  std::ostringstream oss;
  oss << long_value;
  this->value = oss.str();
}

/// Gets a real value as a string
std::string XMLAttrib::str(double value)
{
  if (value == std::numeric_limits<double>::infinity())
    return std::string("inf");
  else if (value == -std::numeric_limits<double>::infinity())
    return std::string("-inf");
  else
  {
    std::ostringstream oss;
    oss << value;
    return oss.str();
  }
}

/// Gets a floating point value from the underlying string representation
double XMLAttrib::get_real_value()
{
  // indicate this attribute has been processed
  processed = true;

  if (strcasecmp(this->value.c_str(), "inf") == 0)
    return std::numeric_limits<double>::infinity();
  else if (strcasecmp(this->value.c_str(), "-inf") == 0)
    return -std::numeric_limits<double>::infinity();
  else
    return (double) std::atof(this->value.c_str());
}

/// Gets a Boolean value from the underlying string representation
/**
 * \note an error message will be output to stderr if the string is not either 
 *         "true" or "false" (case insensitive)
 */
bool XMLAttrib::get_bool_value()
{
  // indicate this attribute has been processed
  processed = true;

  // remove leading and trailing whitespace from the string
  unsigned idx_first = this->value.find_first_not_of(" \t\n");
  unsigned idx_last = this->value.find_last_not_of(" \t\n");
  std::string substring = this->value.substr(idx_first, idx_last - idx_first + 1);

  if (strcasecmp(substring.c_str(), "true") == 0)
    return true;
  else if (strcasecmp(substring.c_str(), "false") == 0)
    return false;
  else
  {
    std::cerr << "XMLAttrib::get_bool_value() - invalid Boolean string: " << this->value << std::endl;
    return false;
  }
}

shared_ptr<const XMLTree> XMLTree::read_from_xml(const std::string& fname)
{
  xmlDoc* doc;

  // *************************************************************
  // going to remove any path from the argument and change to that
  // path; this is done so that all files referenced from the
  // local path of the XML file are found
  // *************************************************************

  // initialize the library and look for potential ABI mismatches
  LIBXML_TEST_VERSION

  // open the file
  if ((doc = xmlReadFile(fname.c_str(), NULL, 0)) == NULL)
  {
    std::cerr << "XMLTree::read_from_xml() - unable to open file " << fname;
    std::cerr << " for reading" << std::endl;
    return shared_ptr<const XMLTree>();
  }
  
  // get the document for parsing
  xmlNode* root = xmlDocGetRootElement(doc);

  // construct the XML tree
  shared_ptr<const XMLTree> node = construct_xml_tree(root);

  // free the XML document
  xmlFreeDoc(doc);

  return node;
}

/// Constructs an XML tree from a xmlNode object
shared_ptr<const XMLTree> XMLTree::construct_xml_tree(xmlNode* root)
{
  // construct a new node
  XMLTreePtr node(new XMLTree(std::string((char*) root->name)));

  // read all of the attributes
  for (xmlAttr* a = root->properties; a; a = a->next)
  {
    std::string name((const char*) a->name);
    std::string value((const char*) a->children->content);
    node->attribs.insert(XMLAttrib(name, value));
  }

  // recursively process all children
  for (xmlNode* n = root->children; n; n=n->next)
    if (n->type == XML_ELEMENT_NODE)
      node->add_child(boost::const_pointer_cast<XMLTree>(construct_xml_tree(n)));

  // look for content
  for (xmlNode* n = root->children; n; n=n->next)
    if (n->type == XML_TEXT_NODE && n->content)
      node->content = std::string((char*) n->content); 

  return node;
}

/// Gets a list of space-delimited and/or comma-delimited strings from the underlying string value
std::list<std::string> XMLAttrib::get_strings_value()
{
  // indicate this attribute has been processed
  processed = true;

  std::list<std::string> values;
  std::string copy = this->value;
  while (true)
  {
    // get the portion of the string before the delimiter
    size_t space_idx = copy.find(' ');
    size_t comma_idx = copy.find(',');
        
    // case 1: no delimiter found
    if (space_idx == std::string::npos && comma_idx == std::string::npos)
    {
      values.push_back(copy);
      break;
    }
        
    // case 2: delimiter found
    if (space_idx != std::string::npos && space_idx < comma_idx)
    {
      values.push_back(copy.substr(0,space_idx));
      copy = copy.substr(space_idx);
    }
    else
    {
      values.push_back(copy.substr(0,comma_idx));
      copy = copy.substr(comma_idx);
    }
        
    // get the new string
    size_t firstidx = copy.find_first_not_of(" ,");
    if (firstidx == std::string::npos)
      break;
    else
      copy = copy.substr(firstidx);        
  }

  return values;
}

/// Gets a list of space-delimited and/or comma-delimited vectors from the underlying string value
void XMLAttrib::get_vector_value(VectorNd& v)
{
  // indicate this attribute has been processed
  processed = true;

  v = VectorNd::parse(value);
}

/// Returns an Origin3d value from the attribute
Origin3d XMLAttrib::get_origin_value() 
{
  // indicate this attribute has been processed
  processed = true;

  VectorNd v = VectorNd::parse(value);
  if (v.size() != 3)
    throw std::runtime_error("Unable to parse origin from vector!");
  Origin3d o;
  o.x() = v[0];
  o.y() = v[1];
  o.z() = v[2];
  return o;
}

/// Returns a quaternion value from the attribute
Quatd XMLAttrib::get_quat_value()
{
  // indicate this attribute has been processed
  processed = true;

  VectorNd v = VectorNd::parse(value);
  if (v.size() != 4)
    throw std::runtime_error("Unable to parse quaternion from vector!");
  Quatd q;
  q.w = v[0];
  q.x = v[1];
  q.y = v[2];
  q.z = v[3];
  return q;
}

/// Returns a quaternion value from a roll-pitch-yaw attribute
Quatd XMLAttrib::get_rpy_value()
{
  // indicate this attribute has been processed
  processed = true;

  VectorNd v = VectorNd::parse(value);
  if (v.size() != 3)
    throw std::runtime_error("Unable to parse roll-pitch-yaw from vector!");
  return Quatd::rpy(v[0], v[1], v[2]);
}

/// Gets a list of space-delimited and/or comma-delimited vectors from the underlying string value
void XMLAttrib::get_vector_value(Vector2d& v)
{
  // indicate this attribute has been processed
  processed = true;

  VectorNd w = VectorNd::parse(value);
  if (w.size() != v.size())
    throw MissizeException();
  v = Vector2d(w[0], w[1]);
}  

/// Gets a list of space-delimited and/or comma-delimited vectors from the underlying string value
void XMLAttrib::get_vector_value(Vector3d& v)
{
  // indicate this attribute has been processed
  processed = true;

  VectorNd w = VectorNd::parse(value);
  if (w.size() != v.size())
    throw MissizeException();
  v = Vector3d(w[0], w[1], w[2]);
}  

/// Gets a list of space-delimited and/or comma-delimited vectors from the underlying string value
void XMLAttrib::get_vector_value(SVector6d& v)
{
  // indicate this attribute has been processed
  processed = true;

  VectorNd w = VectorNd::parse(value);
  if (w.size() != v.size())
    throw MissizeException();
  v = SVector6d(w[0], w[1], w[2], w[3], w[4], w[5]);
}  

/// Gets a list of space-delimited and/or comma-delimited strings from the underlying string value
void XMLAttrib::get_matrix_value(Matrix3d& m)
{
  // indicate this attribute has been processed
  processed = true;

  MatrixNd n;
  get_matrix_value(n);
  if (n.rows() != m.rows() || n.columns() != m.rows())
    throw MissizeException();
  m = Matrix3d(n.data());
} 

/// Gets a list of space-delimited and/or comma-delimited strings from the underlying string value
void XMLAttrib::get_matrix_value(MatrixNd& m)
{
  // indicate this attribute has been processed
  processed = true;

  // construct the list of properties
  std::list<std::list<std::string> > plist;
  std::list<std::string> clist;

  // extract all space-delimited, comma-delimited, and semi-colon-delimited substrings
  std::string copy = this->value;
  while (true)
  {
    // get the portion of the string before the delimiter
    size_t space_idx = copy.find_first_of("  \t");
    size_t comma_idx = copy.find(',');
    size_t semicolon_idx = copy.find(';');
        
    // case 1: no delimiter found
    if (space_idx == std::string::npos && comma_idx == std::string::npos && semicolon_idx == std::string::npos)
    {
      clist.push_back(copy);
      plist.push_back(clist);
      break;
    }
        
    // convert not founds to maximums for sorting
    if (space_idx == std::string::npos)
      space_idx = std::numeric_limits<int>::max();
    if (comma_idx == std::string::npos)
      comma_idx = std::numeric_limits<int>::max();
    if (semicolon_idx == std::string::npos)
      semicolon_idx = std::numeric_limits<int>::max();
        
    // case 2: semicolon delimiter found first
    if (semicolon_idx < space_idx && semicolon_idx <  comma_idx)
    {
      clist.push_back(copy.substr(0,semicolon_idx));
      plist.push_back(clist);
      clist = std::list<std::string>();
      copy = copy.substr(semicolon_idx);
    }    
    // case 3: space delimiter found first
    else if (space_idx < semicolon_idx && space_idx < comma_idx)
    {
      clist.push_back(copy.substr(0,space_idx));
      copy = copy.substr(space_idx);
    }
    // case 4: comma delimiter found first
    else if (comma_idx < semicolon_idx && space_idx > comma_idx)
    {
      clist.push_back(copy.substr(0,comma_idx));
      copy = copy.substr(comma_idx);
    }
        
    // get the new string
    size_t firstidx = copy.find_first_not_of(" ,;");
    if (firstidx == std::string::npos)
      break;
    else
      copy = copy.substr(firstidx);        
  }

  // convert the lists to a Matrix
  // first, verify that all vectors are the same length
  unsigned len = plist.front().size();
  for (std::list<std::list<std::string> >::const_iterator i = plist.begin(); i != plist.end(); i++)
    if (i->size() != len)
    {
      std::cerr << "XMLAttrib::get_matrix_value() - rows are not of the same size!" << std::endl << "  offending string: " << value << std::endl;
      m.resize(0,0);
      return;
    }
      
  // resize the matrix
  m.resize(plist.size(), plist.front().size());
  unsigned r = 0;
  for (std::list<std::list<std::string> >::const_iterator i = plist.begin(); i != plist.end(); i++)
  {
    unsigned s = 0;
    for (std::list<std::string>::const_iterator j = i->begin(); j != i->end(); j++)
    {
      if (strcasecmp(j->c_str(), "inf") == 0)
        m(r,s) = std::numeric_limits<double>::infinity();
      else if (strcasecmp(j->c_str(), "-inf") == 0)
        m(r,s) = -std::numeric_limits<double>::infinity();
      else
        m(r,s) = (double) atof(j->c_str());
      s++;
    }
    r++;
  }
}

/// Sends the specified XMLAttrib to the given stream
std::ostream& Moby::operator<<(std::ostream& out, const XMLAttrib& attr)
{
  out << attr.name << "=\"" << attr.value << "\"";

  return out;
}

/// Constructs a XMLTree with no attributes
XMLTree::XMLTree(const std::string& name)
{
  this->processed = false;
  this->name = name;
}

/// Constructs a XMLTree with the specified list of attributes
XMLTree::XMLTree(const std::string& name, const std::list<XMLAttrib>& attributes)
{
  this->name = name;
  this->attribs = std::set<XMLAttrib>(attributes.begin(), attributes.end());
  this->processed = false;
}

/// Gets the specified attribute
/**
 * \return a pointer to the attribute with the specified name, or NULL if the
 *         requested attribute does not exist
 */
XMLAttrib* XMLTree::get_attrib(const std::string& attrib_name) const
{
  // construct an empty attribute
  XMLAttrib attr(attrib_name, "");

  // look for the attribute
  std::set<XMLAttrib>::const_iterator attrib_iter = this->attribs.find(attr);

  // return the appropriate value
  return (attrib_iter == this->attribs.end()) ? NULL : (XMLAttrib*) &(*attrib_iter);
}

/// Returns a list of all child nodes (not including further descendants) matching any of the names in the given list (case insensitive)
std::list<shared_ptr<const XMLTree> > XMLTree::find_child_nodes(const std::list<std::string>& names) const
{
  std::list<shared_ptr<const XMLTree> > matches;

  // construct a set out of the list of names -- all names are converted to
  // lowercase first
  std::set<std::string> names_set;
  for (std::list<std::string>::const_iterator i = names.begin(); i != names.end(); i++)
  {
    // copy the string
    std::string name = *i;

    // make it lowercase
    std::transform(name.begin(), name.end(), name.begin(), (int(*)(int)) std::tolower);

    // add it to the set
    names_set.insert(name);
  }

  // now, look over all nodes, transforming node names to lowercase, and 
  // checking for a match
  for (std::list<XMLTreePtr>::const_iterator i = this->children.begin(); i != this->children.end(); i++)
  {
    // get the name of the node, lowercased
    std::string node_name = (*i)->name;
    std::transform(node_name.begin(), node_name.end(), node_name.begin(), (int(*)(int)) std::tolower);
    
    // check for a match against the set of names
    if (names_set.find(node_name) != names_set.end())
      matches.push_back(*i);
  }  

  return matches;
}

/// Returns a list of all child nodes (not including further descendants) matching the given name (case insensitive)
std::list<shared_ptr<const XMLTree> > XMLTree::find_child_nodes(const std::string& name) const
{
  std::list<shared_ptr<const XMLTree> > matches;

  // convert the name to lowercase
  std::string name_lower = name;
  std::transform(name_lower.begin(), name_lower.end(), name_lower.begin(), (int(*)(int)) std::tolower);

  // search for a match
  for (std::list<XMLTreePtr>::const_iterator i = this->children.begin(); i != this->children.end(); i++)
  {
    // convert the node name to lowercase
    std::string node_name = (*i)->name;
    std::transform(node_name.begin(), node_name.end(), node_name.begin(), (int(*)(int)) std::tolower);
    
    // check for equality
    if (name_lower == node_name)
      matches.push_back(*i);
  }

  return matches;
}

/// Sends the specified XMLTree node to the given stream
std::ostream& Moby::operator<<(std::ostream& out, const XMLTree& node)
{
  // get the list of child nodes
  const std::list<XMLTreePtr>& child_nodes = node.children; 
  
  // get the set of attributes for this node
  const std::set<XMLAttrib>& attribs = node.attribs;

  // write the start of the tag, the node name, and the attributes
  out << "<" << node.name << " ";
  for (std::set<XMLAttrib>::const_iterator i = attribs.begin(); i != attribs.end(); i++)
    out << *i << " ";

  // if there are no child nodes, close this tag 
  if (child_nodes.empty())
    out << "/>" << std::endl;
  else
  {
    // write the closing tag
    out << ">" << std::endl;
    
    // write all child nodes, in forward order
    for (std::list<XMLTreePtr>::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
      out << **i;

    // close this tag
    out << "</" << node.name << ">" << std::endl;
  }

  return out;
}

