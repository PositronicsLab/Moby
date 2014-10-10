#ifndef _GAUSSMIX_H_
#define _GAUSSMIX_H_

#include <Moby/Primitive.h>

namespace Moby {

class GaussianMixture : public Primitive
{
  public:

    struct Gauss
    {
      double A;        // height of the Gaussian
      double sigma_x;  // variance in x direction
      double sigma_y;  // variance in y direction
      double x0;       // x coordinate of Gaussian center
      double y0;       // y coordinate of Gaussian center
      double th;       // planar rotation of the Gaussian
    };

    void rebuild(const std::vector<Gauss>& gauss);
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual void set_pose(const Ravelin::Pose3d& T);
    virtual BVPtr get_BVH_root(CollisionGeometryPtr geom);
    virtual bool point_inside(BVPtr bv, const Point3d& p, Ravelin::Vector3d& normal) const;
    virtual void get_vertices(boost::shared_ptr<const Ravelin::Pose3d> P, std::vector<Point3d>& vertices);
    virtual bool intersect_seg(BVPtr bv, const LineSeg3& seg, double& t, Point3d& isect, Ravelin::Vector3d& normal) const;
    virtual boost::shared_ptr<const IndexedTriArray> get_mesh(boost::shared_ptr<const Ravelin::Pose3d> P) { return _mesh; }
    virtual osg::Node* create_visualization();
    virtual double calc_dist_and_normal(const Point3d& p, Ravelin::Vector3d& normal) const;
    virtual double calc_signed_dist(boost::shared_ptr<const Primitive> p, boost::shared_ptr<const Ravelin::Pose3d> pose_this, boost::shared_ptr<const Ravelin::Pose3d> pose_p, Point3d& pthis, Point3d& pb) const;

    private:
      static Ravelin::Vector3d grad(const Gauss& g, double x, double y);
      static double f(const Gauss& g, const Point3d& p, const Point3d& q, double t);
      static double df(const Gauss& g, const Point3d& p, const Point3d& q, double t);
      static double newton_raphson(const Gauss& g, const Point3d& p, const Point3d& q);
      static Gauss read_gauss_node(boost::shared_ptr<const XMLTree> node);
      static double gauss(const Gauss& g, double x, double y);
      void construct_vertices();
      void construct_BVs(CollisionGeometryPtr geom);
      void create_mesh();

      /// Note: there are no mass properties, because this can have no mass!
      virtual void calc_mass_properties() { }

      /// Mesh pointer
      boost::shared_ptr<IndexedTriArray> _mesh;

      /// Submesh pair
      std::pair<boost::shared_ptr<const IndexedTriArray>, std::list<unsigned> > _mesh_pair; 

      /// The collision geometry this mixture is associated with
      CollisionGeometryPtr _geom;

      /// The bounding volumes (OBB)
      OBBPtr _root;

      /// The set of Gaussians
      std::vector<Gauss> _gauss;

      /// The mapping from bounding volumes to Gaussians
      std::map<OBBPtr, unsigned> _obbs;

      /// Set of vertices sampled at regular intervals
      std::vector<std::vector<Point3d> > _vertices;
}; // end class

} // end namespace Moby

#endif // _GAUSSMIX_H
