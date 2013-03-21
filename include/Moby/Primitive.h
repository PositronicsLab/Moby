/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _PRIMITIVE_H
#define _PRIMITIVE_H

#include <map>
#include <vector>
#include <boost/shared_ptr.hpp>
#include <Moby/Base.h>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/Triangle.h>
#include <Moby/ThickTriangle.h>
#include <Moby/Constants.h>
#include <Moby/IndexedTriArray.h>

namespace osg {
  class MatrixTransform;
  class Material;
  class Node;
  class Matrixd;
}

namespace Moby {

/// Defines a triangle-mesh-based primitive type used for inertial property calculation and geometry provisions
/**
 * The center-of-mass of derived types may be at the origin of the world,
 * or not.  Additionally, Primitive can take a transformation matrix in its
 * constructor, with which the primitive data (com, inertia matrix, and geometry)
 * can be transformed.
 */
class Primitive : public virtual Base
{
  friend class CSG;

  public:
    Primitive();
    Primitive(const Matrix4& T);
    virtual ~Primitive();
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    void update_visualization();
    void set_mass(Real mass);
    void set_density(Real density);
    virtual boost::shared_ptr<void> save_state() const;
    virtual void load_state(boost::shared_ptr<void> state);
    virtual void set_transform(const Matrix4& T);
    virtual void set_intersection_tolerance(Real tol);
    static void transform_inertia(Real mass, const Matrix3& J_in, const Vector3& com_in, const Matrix4& T, Matrix3& J_out, Vector3& com_out);
    static void transform_inertia(Real mass, const Matrix3& J_in, const Vector3& com_in, const Matrix3& R, Matrix3& J_out, Vector3& com_out);

    /// Gets the current intersection tolerance for this primitive
    Real get_intersection_tolerance() const { return _intersection_tolerance; }

    /// Gets the visualization for this primitive
    virtual osg::Node* get_visualization();
    virtual osg::Node* create_visualization() = 0;

    /// Sets whether this primitive is used for a deformable body
    virtual void set_deformable(bool flag) { _deformable = flag; }

    /// Gets the root bounding volume for this primitive
    virtual BVPtr get_BVH_root() = 0; 

    /// Returns whether this primitive is deformable
    bool is_deformable() const { return _deformable; }

    /// Gets vertices corresponding to the bounding volume
    virtual void get_vertices(BVPtr bv, std::vector<const Vector3*>& vertices) = 0; 

    /// Determines whether a point is inside/on the geometry
    /**
     * \param gs_BV the bounding volume in which the point resides
     * \param p the query point
     * \param normal the normal to the query point, if the point is inside/on
     *        the geometry
     * \return <b>true</b> if the point is inside or on the geometry, 
     *         <b>false</b> otherwise
     */
    virtual bool point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const = 0;

    /// Determines whether a line segment and the shape intersect
    /**
     * \param bv a bounding volume (to speed intersection testing)
     * \param seg the line segment
     * \param t the parameter of the intersection (seg.first + seg.second*t)
     *        (if intersection)
     * \param isect the point of intersection, on return (if any)
     * \param normal the normal to the shape at the point of intersection, 
     *          on return (if any)
     * \return <b>true</b> if intersection, <b>false</b> otherwise 
     */
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const = 0;

    /// Determines whether a line segment and the shape intersect
    /**
     * Special version for self-intersection testing (for deformable geoms).
     * \param u a pointer to the point used to determine the line segment
     * \param bv a bounding volume (to speed intersection testing)
     * \param seg the line segment
     * \param t the parameter of the intersection (seg.first + seg.second*t)
     *        (if intersection)
     * \param isect the point of intersection, on return (if any)
     * \param normal the normal to the shape at the point of intersection, 
     *          on return (if any)
     * \return <b>true</b> if intersection, <b>false</b> otherwise 
     */
    virtual bool intersect_seg(const Vector3* u, BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const
    {
      throw std::runtime_error("Primitive::intersect_seg() not defined!");
      return false;
    }

    /// Gets mesh data for the geometry with the specified bounding volume
    /**
     * \param bv the bounding data from which the corresponding mesh data will
     *           be taken
     * \return a pair consisting of an IndexedTriArray and a list of triangle
     *         indices encapsulated by bv
     */
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv) = 0;

    /// Gets the mass of this primitive
    Real get_mass() const { return _mass; }

    /// Gets the center-of-mass of this primitive
    const Vector3& get_com() const { return _com; }

    /// Gets the transform applied to this primitive 
    const Matrix4& get_transform() const { return _T; } 

    /// Gets the underlying triangle mesh for this primitive 
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh() = 0;

    /// Gets the inertia for this primitive 
    const Matrix3& get_inertia() const { return _J; }

  protected:
    virtual void calc_mass_properties() = 0;

    /// The intersection tolerance for this shape (default 1e-5)
    Real _intersection_tolerance;

    /// The 4x4 rotational/translational transform applied to this primitive
    Matrix4 _T;

    /// The center-of-mass of this primitive
    Vector3 _com;

    /// The mass of this primitive
    Real _mass;

    /// The density of this primitive
    boost::shared_ptr<Real> _density;

    /// The inertia of the primitive
    Matrix3 _J;

    /// Indicates whether the primitive's mesh or vertices have changed
    bool _invalidated;

  private:
    static void to_osg_matrix(const Matrix4& src, osg::Matrixd& tgt);

    struct PrimitiveState
    {
      bool deformable;     // whether primitive is deformable
      Matrix4 T;          // transform of the primitive
      Vector3 com;         // com of the primitive
      Real mass;         // mass of the primitive
      boost::shared_ptr<Real> density;      // density of the primitive
      Matrix3 J;         // inertia matrix of the primitive
      Real intersection_tolerance;  // intersection tolerance of this primitive
    };

    /// Whether the geometry is deformable or not
    bool _deformable;

    /// The visualization transform (possibly NULL)
    osg::MatrixTransform* _vtransform;

    /// The material for the visualization (possibly NULL)
    osg::Material* _mat;
}; // end class

} // end namespace

#endif
