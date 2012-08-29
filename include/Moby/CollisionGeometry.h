/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _COLLISION_GEOMETRY_H
#define _COLLISION_GEOMETRY_H

#include <stack>
#include <list>
#include <vector>
#include <map>
#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <boost/weak_ptr.hpp>
#include <Moby/Constants.h>
#include <Moby/Base.h>
#include <Moby/Triangle.h>
#include <Moby/Matrix4.h>
#include <Moby/Primitive.h>

namespace Moby {

class SingleBody;  

/// Defines collision geometry that may be used (in principle) many ways: for rigid bodies, non-rigid bodies, ...
/**
 *  In principle the geometry may be very complex, and may support
 *  nonstationarity (changing over time) or deformities (e.g., due to 
 *  collision).  Note that while the underlying geometry may be shared, it is
 *  not intended for CollisionGeometry objects to be shared.
 */
class CollisionGeometry : public virtual Base
{
  public:
    CollisionGeometry();
    void set_transform(const Matrix4& transform, bool rel_transform_accounted);  
    void write_vrml(const std::string& filename) const;
    PrimitivePtr set_geometry(PrimitivePtr primitive);
    void set_rel_transform(const Matrix4& transform, bool update_global_transform = false);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);

    template <class OutputIterator>
    OutputIterator get_sub_geometries(OutputIterator begin) const;

    /// Gets the parent of this CollisionGeometry (or NULL if there is no parent)
    boost::shared_ptr<CollisionGeometry> get_parent() const { return (_parent.expired()) ? CollisionGeometryPtr() : CollisionGeometryPtr(_parent); }
    
    /// Sets the parent of this CollisionGeometry (or NULL to indicate no parent)
    void set_parent(boost::weak_ptr<CollisionGeometry> parent) { _parent = parent; }

    /// Gets the transform (in global coordinates) for this CollisionGeometry
    const Matrix4& get_transform() const { return _transform; }

    /// Gets the relative transform for this CollisionGeometry
    const Matrix4& get_rel_transform() const { return _rel_transform; }
    
    /// Gets the single body associated with this CollisionGeometry (if any)
    boost::shared_ptr<SingleBody> get_single_body() const { return (_single_body.expired()) ? SingleBodyPtr() : SingleBodyPtr(_single_body); }

    /// Sets the single body associated with this CollisionGeometry (if any)
    void set_single_body(boost::shared_ptr<SingleBody> s) { _single_body = s; }
    
    /// The set of children of this geometry
    std::list<CollisionGeometryPtr> children;

    /// Gets the geometry for this primitive
    PrimitivePtr get_geometry() const { return _geometry; }

  protected:
    /// The adjusted (i.e., relative transform considered) transform of the CollisionGeometry
    Matrix4 _transform;
  
    /// The relative transform to the CollisionGeometry frame
    Matrix4 _rel_transform;

    /// The underlying geometry
    PrimitivePtr _geometry;

  private:
    bool _rel_transform_identity;
    boost::weak_ptr<SingleBody> _single_body;
    boost::weak_ptr<CollisionGeometry> _parent;
    std::vector<CollisionGeometryPtr> _children;
}; // end class

#include "CollisionGeometry.inl"

} // end namespace

#endif
