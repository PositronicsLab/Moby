/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_XML_TREE_H
#define _MOBY_XML_TREE_H

#include <boost/enable_shared_from_this.hpp>
#include <list>
#include <string>
#include <set>
#include <libxml/parser.h>
#include <libxml/tree.h>
#include <Moby/VectorN.h>
#include <Moby/SVector6.h>
#include <Moby/MatrixN.h>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/Types.h>

namespace Moby {

/// Attributes used for XML nodes
class XMLAttrib
{
  public:
    XMLAttrib(const std::string& name, const std::string& string_value);
    XMLAttrib(const std::string& name, Real real_value);
    XMLAttrib(const std::string& name, int int_value);
    XMLAttrib(const std::string& name, unsigned unsigned_value);
    XMLAttrib(const std::string& name, const VectorN& vector_value);
    XMLAttrib(const std::string& name, const SVector6& vector_value);
    XMLAttrib(const std::string& name, const MatrixN& matrix_value);
    XMLAttrib(const std::string& name, const Matrix3& matrix_value);
    XMLAttrib(const std::string& name, const Matrix4& matrix_value);
    XMLAttrib(const std::string& name, bool bool_value);
    XMLAttrib(const std::string& name, long long_value);
    const std::string& get_string_value() const { return value; }
    static std::string str(Real value);
    Real get_real_value() const;
    int get_int_value() const { return std::atoi(value.c_str()); }
    unsigned get_unsigned_value() const { return (unsigned) std::atoi(value.c_str()); }
    bool get_bool_value() const;
    long get_long_value() const { return std::atol(value.c_str()); }
    std::list<std::string> get_strings_value() const;
    void get_vector_value(VectorN& v) const;
    void get_vector_value(Vector2& v) const;
    void get_vector_value(Vector3& v) const;
    void get_vector_value(SVector6& v) const;
    void get_matrix_value(Matrix3& m) const;
    void get_matrix_value(Matrix4& m) const;
    void get_matrix_value(MatrixN& m) const;
    bool operator==(const XMLAttrib& a) const { return name == a.name; }
    bool operator<(const XMLAttrib& a) const { return name < a.name; }

    /// The name of the attribute
    std::string name;

    /// The value in string form
    std::string value;
};

/// An XML tree used for serialization 
class XMLTree : public boost::enable_shared_from_this<XMLTree>
{
  public:
    XMLTree(const std::string& name);
    XMLTree(const std::string& name, const std::list<XMLAttrib>& attributes);
    static XMLTreeConstPtr read_from_xml(const std::string& name);
    const XMLAttrib* get_attrib(const std::string& attrib_name) const;
    std::list<XMLTreeConstPtr> find_child_nodes(const std::string& name) const;
    std::list<XMLTreeConstPtr> find_child_nodes(const std::list<std::string>& name) const;
    std::list<XMLTreeConstPtr> find_descendant_nodes(const std::string& name) const;

    /// Adds a child tree to this tree; also sets the parent node
    void add_child(XMLTreePtr child) { children.push_back(child); child->set_parent(shared_from_this()); }

    /// Sets the parent of this tree (if any)
    void set_parent(XMLTreePtr parent) { _parent = parent; }

    /// Gets the parent of this tree (if any)
    boost::weak_ptr<XMLTree> get_parent() const { return _parent; }

    /// The set of attributes of this node
    std::set<XMLAttrib> attribs;  

    /// The list of children of this node
    std::list<XMLTreePtr> children;

    /// The name of this node
    std::string name;

    /// The ID of this node
    std::string id;

    /// The object (if any) represented by this node
    boost::shared_ptr<void> object;

  private:
    boost::weak_ptr<XMLTree> _parent;
    static XMLTreeConstPtr construct_xml_tree(xmlNode* root);
}; // end class

std::ostream& operator<<(std::ostream& out, const XMLTree& tree);
std::ostream& operator<<(std::ostream& out, const XMLAttrib& attr);

} // end namespace

#endif

