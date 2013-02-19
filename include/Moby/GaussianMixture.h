#ifndef _GAUSSMIX_H_
#define _GAUSSMIX_H_

#include <Moby/Primitive.h>

namespace Moby {

class GaussianMixture : public Primitive
{
  public:

    struct Gauss
    {
      Real A;        // height of the Gaussian
      Real sigma_x;  // variance in x direction
      Real sigma_y;  // variance in y direction
      Real x0;       // x coordinate of Gaussian center
      Real y0;       // y coordinate of Gaussian center
      Real th;       // planar rotation of the Gaussian
    };

    void rebuild(const std::vector<Gauss>& gauss);
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    virtual void set_transform(const Matrix4& T);
    virtual BVPtr get_BVH_root() { return boost::dynamic_pointer_cast<BV>(_root); }
    virtual void get_vertices(BVPtr bv, std::vector<const Vector3*>& vertices);
    virtual bool point_inside(BVPtr bv, const Vector3& p, Vector3& normal) const;
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, Real& t, Vector3& isect, Vector3& normal) const;
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh() { return _mesh; }
    virtual const std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> >& get_sub_mesh(BVPtr bv);

    // visualization
    #ifdef USE_OSG
    virtual osg::Node* create_visualization();
    #endif

    private:
      static Vector3 grad(const Gauss& g, Real x, Real y);
      static Real f(const Gauss& g, const Vector3& p, const Vector3& q, Real t);
      static Real df(const Gauss& g, const Vector3& p, const Vector3& q, Real t);
      static Real newton_raphson(const Gauss& g, const Vector3& p, const Vector3& q);
      static Gauss read_gauss_node(XMLTreeConstPtr node);
      static Real gauss(const Gauss& g, Real x, Real y);
      void construct_vertices();
      void construct_BVs();
      void create_mesh();

      /// Note: there are no mass properties, because this can have no mass!
      virtual void calc_mass_properties() { }

      /// Mesh pointer
      boost::shared_ptr<IndexedTriArray> _mesh;

      /// Submesh pair
      std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _mesh_pair; 

      /// The root bounding volume (OBB)
      OBBPtr _root;

      /// The set of Gaussians
      std::vector<Gauss> _gauss;

      /// The mapping from bounding volumes to Gaussians
      std::map<OBBPtr, unsigned> _obbs;

      /// Set of vertices sampled at regular intervals
      std::vector<std::vector<Vector3> > _vertices;

      /// Set of submeshes
      std::vector<std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > > _submesh; 
}; // end class

} // end namespace Moby

#endif // _GAUSSMIX_H
