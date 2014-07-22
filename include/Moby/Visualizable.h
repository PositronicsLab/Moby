/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _VISUALIZABLE_H
#define _VISUALIZABLE_H

#include <vector>
#include <Ravelin/Pose3d.h>
#include <Moby/Types.h>
#include <Moby/Base.h>
#include <Moby/OSGGroupWrapper.h>

namespace osg { 
  class Node;
  class MatrixTransform; 
  class Group;
}

namespace Moby {

/// Class that allows for visualizing simulation data
/**
 * This class uses the OSGGroupWrapper class to permit sharing and
 * serialization of visualization data.
 */
class Visualizable : public virtual Base 
{
  public:
    Visualizable();
    Visualizable(const Visualizable* v) : Base(v) { }
    virtual ~Visualizable(); 
    virtual void update_visualization();
    void set_visualization_relative_pose(const Ravelin::Pose3d& P);
    virtual void set_visualization_data(osg::Node* vdata); 
    virtual void set_visualization_data(OSGGroupWrapperPtr vdata); 
    static osg::Group* construct_from_node(boost::shared_ptr<const XMLTree> node, const std::map<std::string, BasePtr>& id_map);
    osg::Group* get_visualization_data() const;
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);

    /// Gets the pose for this visualizable object 
    boost::shared_ptr<const Ravelin::Pose3d> get_visualization_pose() { return _vF; }

  protected:

    /// The relative pose
    boost::shared_ptr<Ravelin::Pose3d> _vF;

    /// The underlying visualization data
    OSGGroupWrapperPtr _vizdata;

    /// The top-level group (containing _vizdata)
    osg::MatrixTransform* _group;
}; // end class
} // end namespace

#endif // (ends #ifndef _VISUALIZABLE_H)

