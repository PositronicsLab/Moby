/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _INDEXED_TETRA_ARRAY_H
#define _INDEXED_TETRA_ARRAY_H

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <Ravelin/Vector3d.h>
#include <Ravelin/Pose3d.h>
#include <Moby/Base.h>
#include <Moby/Tetrahedron.h>
#include <Moby/IndexedTetra.h>

namespace Moby {

/// An array of tetrahedra indexed into shared vertices
class IndexedTetraArray : public Base
{
  public:
    IndexedTetraArray() {}
    IndexedTetraArray(boost::shared_ptr<const std::vector<Point3d> > vertices, const std::vector<IndexedTetra>& facets);
    IndexedTetraArray(boost::shared_ptr<const std::vector<Point3d> > vertices, boost::shared_ptr<const std::vector<IndexedTetra> > facets);

    template <class ForwardIterator1, class ForwardIterator2>
    IndexedTetraArray(ForwardIterator1 vertices, ForwardIterator1 verts_end, ForwardIterator2 facets_begin, ForwardIterator2 facets_end);

    void center();
    unsigned num_tetra() const { return (_tetra) ? _tetra->size() : 0; }
    IndexedTetraArray transform(const Ravelin::Transform3d& T) const;
    IndexedTetraArray compress_vertices() const;
    static IndexedTetraArray read_from_tetra(const std::string& filename);
    static void write_to_tetra(const IndexedTetraArray& mesh, const std::string& filename);
    void write_to_tetra(const std::string& filename) const { write_to_tetra(*this, filename); }
    static void write_to_obj(const IndexedTetraArray& mesh, const std::string& filename);
    void write_to_obj(const std::string& filename) const { write_to_obj(*this, filename); }
    IndexedTetraArray& operator=(const IndexedTetraArray& mesh);
    Tetrahedron get_tetrahedron(unsigned i) const;
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);  
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;

    /// Gets the pose that these vertices are defined in
    /**
     * The mesh never changes this pose.
     */
    boost::shared_ptr<const Ravelin::Pose3d> get_pose() const { return _pose; }

    /// Gets the pointer to the vector of tetrahedra 
    boost::shared_ptr<const std::vector<IndexedTetra> > get_tetra_pointer() const { return _tetra; }

    /// Gets the pointer to the vector of vertices
    boost::shared_ptr<const std::vector<Point3d> > get_vertices_pointer() const { return _vertices; }

    /// Gets the vector of facets
    const std::vector<IndexedTetra>& get_tetra() const { return *_tetra; }

    /// Gets the vector of verties
    const std::vector<Point3d>& get_vertices() const { return *_vertices; }

  private:
    void validate();

    /// The pose that these vertices are defined in
    boost::shared_ptr<const Ravelin::Pose3d> _pose;
    boost::shared_ptr<const std::vector<IndexedTetra> > _tetra;
    boost::shared_ptr<const std::vector<Point3d> > _vertices;

}; // end class

// include inline functions
#include <Moby/IndexedTetraArray.inl>

} // end namespace

#endif

