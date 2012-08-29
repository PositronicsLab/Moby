/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _INDEXED_TETRA_ARRAY_H
#define _INDEXED_TETRA_ARRAY_H

#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <Moby/Base.h>
#include <Moby/Vector3.h>
#include <Moby/Matrix3.h>
#include <Moby/Matrix4.h>
#include <Moby/Tetrahedron.h>
#include <Moby/IndexedTetra.h>

namespace Moby {

/// An array of tetrahedra indexed into shared vertices
class IndexedTetraArray : public Base
{
  public:
    IndexedTetraArray() {}
    IndexedTetraArray(boost::shared_ptr<const std::vector<Vector3> > vertices, const std::vector<IndexedTetra>& facets);
    IndexedTetraArray(boost::shared_ptr<const std::vector<Vector3> > vertices, boost::shared_ptr<const std::vector<IndexedTetra> > facets);

    template <class InputIterator1, class InputIterator2>
    IndexedTetraArray(InputIterator1 vertices, InputIterator1 verts_end, InputIterator2 facets_begin, InputIterator2 facets_end);

    void center();
    unsigned num_tetra() const { return (_tetra) ? _tetra->size() : 0; }
    IndexedTetraArray transform(const Matrix4& T) const;
    IndexedTetraArray rotate_scale(const Matrix3& T) const;
    IndexedTetraArray translate(const Vector3& v) const;
    IndexedTetraArray compress_vertices() const;
    static IndexedTetraArray read_from_tetra(const std::string& filename);
    static void write_to_tetra(const IndexedTetraArray& mesh, const std::string& filename);
    void write_to_tetra(const std::string& filename) const { write_to_tetra(*this, filename); }
    static void write_to_obj(const IndexedTetraArray& mesh, const std::string& filename);
    void write_to_obj(const std::string& filename) const { write_to_obj(*this, filename); }
    IndexedTetraArray& operator=(const IndexedTetraArray& mesh);
    Tetrahedron get_tetrahedron(unsigned i) const;
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);  
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;

    /// Gets the pointer to the vector of tetrahedra 
    boost::shared_ptr<const std::vector<IndexedTetra> > get_tetra_pointer() const { return _tetra; }

    /// Gets the pointer to the vector of vertices
    boost::shared_ptr<const std::vector<Vector3> > get_vertices_pointer() const { return _vertices; }

    /// Gets the vector of facets
    const std::vector<IndexedTetra>& get_tetra() const { return *_tetra; }

    /// Gets the vector of verties
    const std::vector<Vector3>& get_vertices() const { return *_vertices; }

  private:
    void validate();

    boost::shared_ptr<const std::vector<IndexedTetra> > _tetra;
    boost::shared_ptr<const std::vector<Vector3> > _vertices;

}; // end class

// include inline functions
#include <Moby/IndexedTetraArray.inl>

} // end namespace

#endif

