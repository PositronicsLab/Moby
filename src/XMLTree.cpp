/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <string.h>
#include <cctype>
#include <algorithm>
#include <iostream>
#include <limits>
#include <cmath>
#include <sstream>
#include <Moby/MissizeException.h>
#include <Moby/MatrixN.h>
#include <Moby/XMLTree.h>

using namespace Moby;

/// Constructs an XMLAttrib object from a name and a string value
XMLAttrib::XMLAttrib(const std::string& name, const std::string& string_value)
{
  this->name = name;
  this->value = string_value;
}

/// Constructs a real-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, Real real_value)
{
  this->name = name;
  this->value = str(real_value);
}

/// Constructs an integer-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, int int_value)
{
  this->name = name;
  std::ostringstream oss;
  oss << int_value;
  this->value = oss.str();
}

/// Constructs an unsigned integer-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, unsigned unsigned_value)
{
  this->name = name;
  std::ostringstream oss;
  oss << unsigned_value;
  this->value = oss.str();
}

/// Constructs a vector-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const VectorN& vector_value)
{
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
XMLAttrib::XMLAttrib(const std::string& name, const SVector6& vector_value)
{
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
XMLAttrib::XMLAttrib(const std::string& name, const MatrixN& matrix_value)
{  
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
XMLAttrib::XMLAttrib(const std::string& name, const Matrix3& matrix_value)
{  
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

/// Constructs a vector-valued attribute with the given name
XMLAttrib::XMLAttrib(const std::string& name, const Matrix4& matrix_value)
{  
  this->name = name;
  std::ostringstream oss;

  // set the first value of the matrix
  oss << str(matrix_value(0,0));

  // for each row of the matrix
  for (unsigned j=0; j< 3; j++)
  {
    // determine column iteration
    unsigned i = (j == 0) ? 1 : 0;

    // for each column of the matrix
    for (; i< matrix_value.columns(); i++)
      oss << " " << str(matrix_value(j,i));

    // separate rows with a semicolon
    oss << ";";
  }

  // write the fourth row
  oss << " 0 0 0 1";

  // get the string value
  this->value = oss.str();
}

/// Constructs a Boolean-valued attribute from the given value
XMLAttrib::XMLAttrib(const std::string& name, bool bool_value)
{
  this->name = name;
  this->value = (bool_value) ? "true" : "false";  
}

/// Constructs a long-valued attribute from the given value
XMLAttrib::XMLAttrib(const std::string& name, long long_value)
{
  this->name = name;
  std::ostringstream oss;
  oss << long_value;
  this->value = oss.str();
}

/// Gets a real value as a string
std::string XMLAttrib::str(Real value)
{
  if (value == std::numeric_limits<Real>::infinity())
    return std::string("inf");
  else if (value == -std::numeric_limits<Real>::infinity())
    return std::string("-inf");
  else
  {
    std::ostringstream oss;
    oss << value;
    return oss.str();
  }
}

/// Gets a floating point value from the underlying string representation
Real XMLAttrib::get_real_value() const
{
  if (strcasecmp(this->value.c_str(), "inf") == 0)
    return std::numeric_limits<Real>::infinity();
  else if (strcasecmp(this->value.c_str(), "-inf") == 0)
    return -std::numeric_limits<Real>::infinity();
  else
    return (Real) std::atof(this->value.c_str());
}

/// Gets a Boolean value from the underlying string representation
/**
 * \note an error message will be output to stderr if the string is not either 
 *         "true" or "false" (case insensitive)
 */
bool XMLAttrib::get_bool_value() const
{
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
    std::cerr << "XMLTree::get_bool_value() - invalid Boolean string: " << this->value << std::endl;
    return false;
  }
}

XMLTreeConstPtr XMLTree::read_from_xml(const std::string& fname)
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
    return XMLTreeConstPtr();
  }
  
  // get the document for parsing
  xmlNode* root = xmlDocGetRootElement(doc);

  // construct the XML tree
  XMLTreeConstPtr node = construct_xml_tree(root);

  // free the XML document
  xmlFreeDoc(doc);

  return node;
}

/// Constructs an XML tree from a xmlNode object
XMLTreeConstPtr XMLTree::construct_xml_tree(xmlNode* root)
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
    node->add_child(boost::const_pointer_cast<XMLTree>(construct_xml_tree(n)));

  return node;
}

/// Gets a list of space-delimited and/or comma-delimited strings from the underlying string value
std::list<std::string> XMLAttrib::get_strings_value() const
{
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
void XMLAttrib::get_vector_value(VectorN& v) const
{
  v = VectorN::parse(value);
}

/// Gets a list of space-delimited and/or comma-delimited vectors from the underlying string value
void XMLAttrib::get_vector_value(Vector2& v) const
{
  VectorN w = VectorN::parse(value);
  if (w.size() != v.size())
    throw MissizeException();
  v = Vector2(w.begin());
}  

/// Gets a list of space-delimited and/or comma-delimited vectors from the underlying string value
void XMLAttrib::get_vector_value(Vector3& v) const
{
  VectorN w = VectorN::parse(value);
  if (w.size() != v.size())
    throw MissizeException();
  v = Vector3(w.begin());
}  

/// Gets a list of space-delimited and/or comma-delimited vectors from the underlying string value
void XMLAttrib::get_vector_value(SVector6& v) const
{
  VectorN w = VectorN::parse(value);
  if (w.size() != v.size())
    throw MissizeException();
  v = SVector6(w.begin());
}  

/// Gets a list of space-delimited and/or comma-delimited strings from the underlying string value
void XMLAttrib::get_matrix_value(Matrix3& m) const
{
  MatrixN n;
  get_matrix_value(n);
  if (n.rows() != m.size() || n.columns() != m.size())
    throw MissizeException();
  m = Matrix3(n.begin());
} 

/// Gets a list of space-delimited and/or comma-delimited strings from the underlying string value
void XMLAttrib::get_matrix_value(Matrix4& m) const
{
  const unsigned X = 0, W = 3;
  MatrixN n;
  get_matrix_value(n);
  if (n.rows() != m.size() || n.columns() != m.size())
    throw MissizeException();
  MatrixN o = n.get_sub_mat(X,W,X,W+1);
  m = Matrix4(o.begin());
} 

/// Gets a list of space-delimited and/or comma-delimited strings from the underlying string value
void XMLAttrib::get_matrix_value(MatrixN& m) const
{
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
      std::cerr << "XMLTree::get_matrix_value() - rows are not of the same size!" << std::endl << "  offending string: " << value << std::endl;
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
        m(r,s) = std::numeric_limits<Real>::infinity();
      else if (strcasecmp(j->c_str(), "-inf") == 0)
        m(r,s) = -std::numeric_limits<Real>::infinity();
      else
        m(r,s) = (Real) atof(j->c_str());
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
  this->name = name;
}

/// Constructs a XMLTree with the specified list of attributes
XMLTree::XMLTree(const std::string& name, const std::list<XMLAttrib>& attributes)
{
  this->name = name;
  this->attribs = std::set<XMLAttrib>(attributes.begin(), attributes.end());
}

/// Gets the specified attribute
/**
 * \return a pointer to the attribute with the specified name, or NULL if the
 *         requested attribute does not exist
 */
const XMLAttrib* XMLTree::get_attrib(const std::string& attrib_name) const
{
  // construct an empty attribute
  XMLAttrib attr(attrib_name, "");

  // look for the attribute
  std::set<XMLAttrib>::const_iterator attrib_iter = this->attribs.find(attr);

  // return the appropriate value
  return (attrib_iter == this->attribs.end()) ? NULL : &(*attrib_iter);
}

/// Returns a list of all child nodes (not including further descendants) matching any of the names in the given list (case insensitive)
std::list<XMLTreeConstPtr> XMLTree::find_child_nodes(const std::list<std::string>& names) const
{
  std::list<XMLTreeConstPtr> matches;

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
std::list<XMLTreeConstPtr> XMLTree::find_child_nodes(const std::string& name) const
{
  std::list<XMLTreeConstPtr> matches;

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

