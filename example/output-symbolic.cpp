/* 
 * This utility outputs a given reduced-coordinate articulated body to a file
 * that can be used to process the kinematics and dynamics of the body 
 * symbolically.  Symbolic versions of the kinematics and dynamics are
 * generally considerably faster.
 */

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <Moby/Types.h>
#include <Moby/XMLReader.h>
#include <Moby/RCArticulatedBody.h>
#include <Moby/Joint.h>
#include <Moby/RevoluteJoint.h>
#include <Moby/PrismaticJoint.h>
#include <Moby/SphericalJoint.h>
#include <Moby/UniversalJoint.h>

using namespace Moby;

std::string generic_transform(const std::string& prefix)
{
  std::string s = "{{";
  s += "rxx";
  s += prefix;
  s += ", rxy";
  s += prefix;
  s += ", rxz";
  s += prefix;
  s += ", tx";
  s += prefix;
  s += "}, {ryx";
  s += prefix;
  s += ", ryy";
  s += prefix;
  s += ", ryz";
  s += prefix;
  s += ", ty";
  s += prefix;
  s += "}, {rzx";
  s += prefix;
  s += ", rzy";
  s += prefix;
  s += ", rzz";
  s += prefix;
  s += ", tz";
  s += prefix;
  s += "}, {0, 0, 0, 1}}";
  return s;
}

void replace(std::string& s, char c1, char c2)
{
  size_t n;
  while ((n = s.find_first_of(c1)) != std::string::npos)
    s[n] = c2;
}

/// Make string safe for Mathematica
std::string msafe(const std::string& s)
{
  std::string ns = s;
  replace(ns, '{', 'a');
  replace(ns, '}', 'b');
  replace(ns, '-', 'c');
  replace(ns, '+', 'd');
  replace(ns, '_', 'e');
  return ns;
}

/// Forms a matrix of n x m symbolic symbols
std::string form_symbolic(const std::string& prefix, const std::string& suffix, unsigned n, unsigned m)
{
  std::ostringstream s;
  s << "{ ";

  for (unsigned i=0; i< n; i++)
  {
    if (i > 0)
      s << ",";
    s << "{ ";
    s << prefix << i << "0" << suffix;
    for (unsigned j=1; j< m; j++)
      s << ", " << prefix << i << j << suffix;
    s << "} ";
  }

  s << "}";
  return s.str();
}

/// Forms a vector of n symbolic symbols
std::string form_symbolic(const std::string& prefix, const std::string& suffix, unsigned n)
{
  std::ostringstream s;
  s << "{ ";
  
  for (unsigned i=0; i< n; i++)
  {
    s << "{";
    s << prefix;
    s << i;
    s << suffix;
    s << "}";
    if (i != n-1)
      s << ",";
    s << " ";
  }
  s << "}";

  return s.str();
}

/// Converts a Moby vector to a Mathematica column vector
std::string to_mathematica_columnv(const VectorNd& v)
{
  std::ostringstream s; 
  s << "{ {" << (double) v[0] << "} ";
  for (unsigned i=1; i< v.size(); i++)
    s << ", {" << (double) v[i] << "}";
  s << " }";
  
  return s.str();
}

/// Converts a Moby matrix to Mathematica format
std::string to_mathematica(const MatrixNd& m)
{
  std::ostringstream s; 
  s << "{ ";
  for (unsigned i=0; i< m.rows(); i++)
  {
    s << "{" << (double) m(i,0);
    for (unsigned j=1; j< m.columns(); j++)
      s << ", " << (double) m(i,j);
    s << "}";
    if (i != m.rows() - 1)
      s << ",";
    s << " ";
  }
  s << "}";

  return s.str();
}

int main(int argc, char* argv[])
{
  // verify syntax correct
  if (argc != 4)
  {
    std::cerr << "syntax: output_symbolic <XML file> <body ID> <output filename>" << std::endl;
    return -1;
  }

  // get the arguments
  std::string xml_file(argv[1]);
  std::string body_ID(argv[2]);
  std::string output_file(argv[3]);

  // read in the XML file
  std::map<std::string, BasePtr> obj_map = XMLReader::read(xml_file);

  // look for the object
  BasePtr object = obj_map[body_ID];
  if (!object)
  {
    std::cerr << "output_symbolic error: unable to find object '" << body_ID << "'" << std::endl;
    return -1;
  }

  // convert the object to a RCArticulated body
  RCArticulatedBodyPtr body = boost::dynamic_pointer_cast<RCArticulatedBody>(object);
  if (!body)
  {
    std::cerr << "output_symbolic error: unable to cast body as RCArticulatedBody!" << std::endl;
    return -1;
  }

  // open the output file for writing
  std::ofstream out(output_file.c_str());
  if (out.fail())
  {
    std::cerr << "output_symbolic error: unable to open file for writing" << std::endl;
    return -1;
  }

  // write a comment
  out << "(* Generated automatically from output_symbolic *)" << std::endl << std::endl;

  /// write the base link, if the link is fixed
  const std::vector<RigidBodyPtr>& links = body->get_links();
  unsigned start_index;
  if (body->is_floating_base())
    start_index = 0;
  else
  {
    start_index = 1;
    std::string linkname = msafe(links.front()->id) + "L";
    out << linkname << "[\"mass\"] = 0;" << std::endl;
    out << linkname << "[\"inertia\"] = {{0,0,0}, {0,0,0}, {0,0,0}};" << std::endl;
    out << linkname << "[\"I\"] = " << form_symbolic("I", linkname, 6, 6) << ";" << std::endl; 
    out << linkname << "[\"transform\"] = " << to_mathematica(links.front()->get_transform()) << ";" << std::endl;
    out << linkname << "[\"parent\"] = Null;" << std::endl;
    out << linkname << "[\"v\"] = Array[0&, {6, 1}];" << std::endl;
    out << linkname << "[\"fext\"] = { {0}, {0}, {0} };" << std::endl;
    out << linkname << "[\"text\"] = { {0} , {0}, {0} };" << std::endl;
    out << linkname << "[\"fixed-base\"] = True;" << std::endl;
    out << std::endl;
  }      

  // write the links
  for (unsigned i=start_index; i< links.size(); i++)
  {
    // get the parent link (if we need it)
    RigidBodyPtr parent;
    if (i > 0)
      parent = RigidBodyPtr(links[i]->get_parent_link());

    // setup the link name
    std::string linkname = msafe(links[i]->id) + "L";

    // write link properties
//    out << linkname << "[\"index\"] = " << links[i]->get_link_index() << ";" << std::endl;
    out << linkname << "[\"mass\"] = " << links[i]->get_mass() << ";" << std::endl;
    out << linkname << "[\"inertia\"] = " << to_mathematica(links[i]->get_inertia()) << ";" << std::endl;
    out << linkname << "[\"I\"] = " << form_symbolic("I", linkname, 6, 6) << ";" << std::endl; 
    out << linkname << "[\"fext\"] = " << form_symbolic("f", linkname, 3) << ";" << std::endl;
    out << linkname << "[\"text\"] = " << form_symbolic("t", linkname, 3) << ";" << std::endl;
    out << linkname << "[\"transform\"] = " << generic_transform(linkname) << ";" << std::endl;
    out << linkname << "[\"v\"] = " << form_symbolic("v", linkname, 6) << ";" << std::endl;
    if (i == 0)
    {
      out << linkname << "[\"fixed-base\"] = False;" << std::endl;
      out << linkname << "[\"parent\"] = Null;" << std::endl;
      out << linkname << "[\"q\"] = {};" << std::endl;
      out << linkname << "[\"qd\"] = {};" << std::endl;
      out << linkname << "[\"qdd\"] = {};" << std::endl;
      out << linkname << "[\"qf\"] = {};" << std::endl;
      out << linkname << "[\"ijointtype\"] = \"none\";" << std::endl;
      out << linkname << "[\"ijointaxis\"] = {};" << std::endl;
      out << linkname << "[\"ijointcomvec\"] = {};" << std::endl;
      out << linkname << "[\"parentcomijointvec\"] = {};" << std::endl;
    }
    else
    {
      out << linkname << "[\"parent\"] = " << msafe(parent->id) << "L;" << std::endl;
      JointPtr joint(links[i]->get_inner_joint_implicit());
      out << linkname << "[\"sIs\"] = " << form_symbolic("sIs", linkname, joint->num_dof(), joint->num_dof()) << ";" << std::endl;
      out << linkname << "[\"q\"] = " << form_symbolic("q", linkname, joint->num_dof()) << ";" << std::endl;
      out << linkname << "[\"qd\"] = " << form_symbolic("qd", linkname, joint->num_dof()) << ";" << std::endl;
      out << linkname << "[\"qdd\"] = " << form_symbolic("qdd", linkname, joint->num_dof()) << ";" << std::endl;
      out << linkname << "[\"qf\"] = " << form_symbolic("qf", linkname, joint->num_dof()) << ";" << std::endl;
      out << linkname << "[\"ijointtype\"] = "; 
      if (boost::dynamic_pointer_cast<RevoluteJoint>(joint))
      {
        boost::shared_ptr<RevoluteJoint> rj = boost::dynamic_pointer_cast<RevoluteJoint>(joint);
        out << "\"revolute\";" << std::endl;
        out << linkname << "[\"ijointaxis1\"] = " << to_mathematica_columnv(rj->get_axis_local()) << ";" << std::endl;
      }
      else if (boost::dynamic_pointer_cast<PrismaticJoint>(joint))
      {
        boost::shared_ptr<PrismaticJoint> pj = boost::dynamic_pointer_cast<PrismaticJoint>(joint);
        out << "\"prismatic\";" << std::endl;
        out << linkname << "[\"ijointaxis1\"] = " << to_mathematica_columnv(pj->get_axis_local()) << ";" << std::endl;
      }
      else if (boost::dynamic_pointer_cast<SphericalJoint>(joint))
      {
        boost::shared_ptr<SphericalJoint> sj = boost::dynamic_pointer_cast<SphericalJoint>(joint);
        out << "\"spherical\";" << std::endl;
        out << linkname << "[\"ijointaxis1\"] = " << to_mathematica_columnv(sj->get_axis_local(SphericalJoint::eAxis1)) << ";" << std::endl;
        out << linkname << "[\"ijointaxis2\"] = " << to_mathematica_columnv(sj->get_axis_local(SphericalJoint::eAxis2)) << ";" << std::endl;
        out << linkname << "[\"ijointaxis3\"] = " << to_mathematica_columnv(sj->get_axis_local(SphericalJoint::eAxis3)) << ";" << std::endl;
      }
      else if (boost::dynamic_pointer_cast<UniversalJoint>(joint))
      {
        boost::shared_ptr<UniversalJoint> uj = boost::dynamic_pointer_cast<UniversalJoint>(joint);
        out << "\"universal\";" << std::endl;
        out << linkname << "[\"ijointaxis1\"] = " << to_mathematica_columnv(uj->get_axis_local(UniversalJoint::eAxis1)) << ";" << std::endl;
        out << linkname << "[\"ijointaxis2\"] = " << to_mathematica_columnv(uj->get_axis_local(UniversalJoint::eAxis2)) << ";" << std::endl;
      }
      const RigidBody::InnerJointData& ijd = links[i]->get_inner_joint_data(joint);
      const RigidBody::OuterJointData& ojd = parent->get_outer_joint_data(joint);
      out << linkname << "[\"ijointcomvec\"] = " << to_mathematica_columnv(ijd.joint_to_com_vec_of) << ";" << std::endl;
      out << linkname << "[\"parentcomijointvec\"] = " << to_mathematica_columnv(ojd.com_to_joint_vec) << ";" << std::endl;
    }
    out << std::endl;
  }

  // write the link children
  out << "(* link children *)" << std::endl;
  for (unsigned i=0; i< links.size(); i++)
  {
    std::string linkname = msafe(links[i]->id) + "L";
    out << linkname << "[\"children\"] = { ";
    std::list<RigidBodyPtr> child_links;
    links[i]->get_child_links(std::back_inserter(child_links));
    for (std::list<RigidBodyPtr>::const_iterator j = child_links.begin(); j != child_links.end(); j++)
    {
      std::string childname = msafe((*j)->id) + "L";
      out << childname;
      std::list<RigidBodyPtr>::const_iterator k = j;
      k++;
      if (k != child_links.end())
        out << ",";
      out << " ";
    }
    out << " };" << std::endl;
  }
  out << std::endl;

  // write the list of links
  out << "(* list of links *)" << std::endl;
  out << "links = { ";
  for (unsigned i=0; i< links.size(); i++)
  {
    out << msafe(links[i]->id) << "L";
    if (i != links.size() - 1)
      out << ",";
    out << " ";
  }
  out << "}; " << std::endl << std::endl;

  // setup pointers to list of links
  out << "(* setup pointers to list of links *)" << std::endl;
  for (unsigned i=0; i< links.size(); i++)
    out << msafe(links[i]->id) << "L[\"links\"] = links;" << std::endl;
  out << std::endl;

  // close the file
  out.close();      
}

