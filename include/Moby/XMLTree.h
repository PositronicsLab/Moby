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
#include <Ravelin/Origin3d.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/SVector6d.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/Matrix3d.h>
#include <Ravelin/Pose3d.h>
#include <Moby/Types.h>

namespace Moby {

/// Attributes used for XML nodes
class XMLAttrib
{
  public:
    XMLAttrib(const std::string& name, const std::string& string_value);
    XMLAttrib(const std::string& name, double real_value);
    XMLAttrib(const std::string& name, int int_value);
    XMLAttrib(const std::string& name, unsigned unsigned_value);
    XMLAttrib(const std::string& name, const std::vector<Ravelin::Twistd>& twist_value);
    XMLAttrib(const std::string& name, const Ravelin::Vector2d& vector_value);
    XMLAttrib(const std::string& name, const Ravelin::Point2d& vector_value);
    XMLAttrib(const std::string& name, const Ravelin::Vector3d& vector_value);
    XMLAttrib(const std::string& name, const Ravelin::Point3d& vector_value);
    XMLAttrib(const std::string& name, const Ravelin::VectorNd& vector_value);
    XMLAttrib(const std::string& name, const Ravelin::SVector6d& vector_value);
    XMLAttrib(const std::string& name, const Ravelin::MatrixNd& matrix_value);
    XMLAttrib(const std::string& name, const Ravelin::Matrix3d& matrix_value);
    XMLAttrib(const std::string& name, const Ravelin::Quatd& quat_value);
    XMLAttrib(const std::string& name, const Ravelin::Origin3d& origin_value);
    XMLAttrib(const std::string& name, bool bool_value);
    XMLAttrib(const std::string& name, long long_value);
    const std::string& get_string_value() const { return value; }
    static std::string str(double value);
    double get_real_value() const;
    Ravelin::Origin3d get_origin_value() const;
    Ravelin::Point3d get_point_value() const;
    std::vector<Ravelin::Twistd> get_twist_values() const;
    int get_int_value() const { return std::atoi(value.c_str()); }
    unsigned get_unsigned_value() const { return (unsigned) std::atoi(value.c_str()); }
    bool get_bool_value() const;
    long get_long_value() const { return std::atol(value.c_str()); }
    std::list<std::string> get_strings_value() const;
    Ravelin::Quatd get_quat_value() const;
    Ravelin::Quatd get_axis_angle_value() const;
    Ravelin::Quatd get_rpy_value() const;
    void get_vector_value(Ravelin::VectorNd& v) const;
    void get_vector_value(Ravelin::Vector2d& v) const;
    void get_vector_value(Ravelin::Vector3d& v) const;
    void get_vector_value(Ravelin::SVector6d& v) const;
    void get_matrix_value(Ravelin::Matrix3d& m) const;
    void get_matrix_value(Ravelin::MatrixNd& m) const;
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
    static boost::shared_ptr<const XMLTree> read_from_xml(const std::string& name);
    const XMLAttrib* get_attrib(const std::string& attrib_name) const;
    std::list<boost::shared_ptr<const XMLTree> > find_child_nodes(const std::string& name) const;
    std::list<boost::shared_ptr<const XMLTree> > find_child_nodes(const std::list<std::string>& name) const;
    std::list<boost::shared_ptr<const XMLTree> > find_descendant_nodes(const std::string& name) const;

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
    static boost::shared_ptr<const XMLTree> construct_xml_tree(xmlNode* root);
}; // end class

std::ostream& operator<<(std::ostream& out, const XMLTree& tree);
std::ostream& operator<<(std::ostream& out, const XMLAttrib& attr);

} // end namespace

#endif

