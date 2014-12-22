#include <Moby/CollisionDetection.h>
#include <Moby/CCD.h>
#include <Moby/EventDrivenSimulator.h>
#define NDEBUG
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

const double R = 1.0;  // radius from center of mass to center of blade
const double W = 0.1;     // width of a blade
const unsigned N_SPOKES = 16;

class BladePlanePlugin : public CollisionDetection
{
  private:
    boost::shared_ptr<EventDrivenSimulator> sim;
    boost::shared_ptr<CCD> ccd;
    RigidBodyPtr wheel;
    RigidBodyPtr ground_body;
    CollisionGeometryPtr ground_cg, wheel_cg, right_foot_cg;

  public:
    BladePlanePlugin() {}

    virtual void set_simulator(boost::shared_ptr<EventDrivenSimulator> sim)
    {
      this->sim = sim;

      // find the necessary objects
      for (unsigned i=0; i< sim->get_dynamic_bodies().size(); i++)
      {
        DynamicBodyPtr body = sim->get_dynamic_bodies()[i];
        if (body->id == "WHEEL")
          wheel = boost::dynamic_pointer_cast<RigidBody>(body);
        else if (body->id == "GROUND")
          ground_body = dynamic_pointer_cast<RigidBody>(body);
      }

      // make sure bodies were found
      if (!wheel)
        throw std::runtime_error("Didn't find body with id 'wheel'");
      if (!ground_body)
        throw std::runtime_error("Didn't find body with id 'ground'");

      // get the geometries
      ground_cg = ground_body->geometries.front();
      wheel_cg = wheel->geometries.front();

      // create the continuous collision detection system
      ccd = boost::shared_ptr<CCD>(new CCD);
    }

    virtual ~BladePlanePlugin() {}

    /// Does broad phase collision detection
    virtual void broad_phase(double dt, const std::vector<DynamicBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
    {
      // clear to_check
      to_check.clear();

      // remove the wheel from bodies
      std::vector<DynamicBodyPtr> remainder = bodies;
      for (unsigned i=0; i< remainder.size(); i++)
        if (remainder[i] == wheel)
        {
          remainder[i] = remainder.back();
          remainder.pop_back();
          break;
        }

      // call CCD on the remainder
      ccd->broad_phase(dt, remainder, to_check);

      // now add in collision geometries for the plane and the walker
      to_check.push_back(std::make_pair(ground_cg, wheel_cg));
    }

    /// Computes a conservative advancement step for general integration
    virtual double calc_CA_step(const PairwiseDistInfo& pdi)
    {
      return ccd->calc_CA_step(pdi);
    }

    /// Computes a conservative advancement step for Euler integration
    virtual double calc_CA_Euler_step(const PairwiseDistInfo& pdi)
    {
      return ccd->calc_CA_Euler_step(pdi);
    }

    /// Calculates signed distance between a blade and a plane
    double calc_signed_dist_blade_plane(CollisionGeometryPtr blade_cg, CollisionGeometryPtr ground_cg)
    {
      Vector3d point, normal;
      return calc_signed_dist_blade_plane(blade_cg,ground_cg,point,normal);
    }

    /// Calculates signed distance between a blade and a plane
    double calc_signed_dist_blade_plane(CollisionGeometryPtr blade_cg, CollisionGeometryPtr ground_cg,Vector3d& point,Vector3d& normal)
    {
      const unsigned Y = 1;

      // setup the minimum dist
      double min_dist = std::numeric_limits<double>::max();

      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(ground_cg->get_geometry());

      // get the pose for the plane primitive
      shared_ptr<const Pose3d> Pplane = plane_geom->get_pose(ground_cg);

      // get the pose for the blade
      shared_ptr<const Pose3d> Pblade = blade_cg->get_pose();

      // get the pose of the center-of-mass of the wheel
      Transform3d pTb = Pose3d::calc_relative_pose(Pblade, Pplane);

      // loop over each spoke
      for (unsigned i=0; i< N_SPOKES; i++)
      {
        // get the value of theta
        double theta = M_PI * i * 2.0 / N_SPOKES;

        // setup the two points of each spoke
        Point3d p1(std::cos(theta)*R, W*.5, std::sin(theta)*R, Pblade);
        Point3d p2(std::cos(theta)*R, -W*.5, std::sin(theta)*R, Pblade);

        // transform the two points
        Point3d p1_plane = pTb.transform_point(p1);
        Point3d p2_plane = pTb.transform_point(p2);

        // get the signed distance
        if (p1_plane[Y] < min_dist)
        {
          min_dist = p1_plane[Y];
          point = p1_plane;
          normal = Vector3d(0,1,0,Pplane);
        }
        if (p2_plane[Y] < min_dist)
        {
          min_dist = p2_plane[Y];
          point = p2_plane;
          normal = Vector3d(0,1,0,Pplane);
        }
      }

      return min_dist;
    }

    /// Finds contacts between a blade and a plane
    virtual void find_contacts_blade_plane(CollisionGeometryPtr blade_cg, CollisionGeometryPtr ground_cg, std::vector<UnilateralConstraint>& contacts)
    {
      const unsigned Y = 1;

      // setup the minimum dist
      double min_dist = std::numeric_limits<double>::max();

      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(ground_cg->get_geometry());

      // get the pose for the plane primitive
      shared_ptr<const Pose3d> Pplane = plane_geom->get_pose(ground_cg);

      // get the pose for the blade
      shared_ptr<const Pose3d> Pblade = blade_cg->get_pose();

      // get the pose of the center-of-mass of the wheel
      Transform3d pTb = Pose3d::calc_relative_pose(Pblade, Pplane);

      // loop over each spoke
      for (unsigned i=0; i< N_SPOKES; i++)
      {
        // get the value of theta
        double theta = M_PI * i * 2.0 / N_SPOKES;

        // setup the two points of each spoke
        Point3d p1(std::cos(theta)*R, W*.5, std::sin(theta)*R, Pblade);
        Point3d p2(std::cos(theta)*R, -W*.5, std::sin(theta)*R, Pblade);

        // transform the two points
        Point3d p1_plane = pTb.transform_point(p1);
        Point3d p2_plane = pTb.transform_point(p2);

        // create contact for the first point on the blade (if appropriate) 
        if (p1_plane[Y] < min_dist)
        {
          min_dist = p1_plane[Y];
          Point3d p1_global = Pose3d::transform_point(GLOBAL, p1); 
          Vector3d normal(0,1,0,Pplane);
          Vector3d normal_global = Pose3d::transform_vector(GLOBAL, normal);
          contacts.clear();
          contacts.push_back(
            CollisionDetection::create_contact(blade_cg,ground_cg,p1_global,normal_global,min_dist)
              );
        }
        else if (p1_plane[Y] == min_dist)
        {
          Point3d p1_global = Pose3d::transform_point(GLOBAL, p1); 
          Vector3d normal(0,1,0,Pplane);
          Vector3d normal_global = Pose3d::transform_vector(GLOBAL, normal);
          contacts.push_back(
            CollisionDetection::create_contact(blade_cg,ground_cg,p1_global,normal_global,min_dist)
              );
        }

        // create contact for the second point on the blade (if appropriate) 
        if (p2_plane[Y] < min_dist)
        {
          min_dist = p2_plane[Y];
          Point3d p2_global = Pose3d::transform_point(GLOBAL, p2); 
          Vector3d normal(0,1,0,Pplane);
          Vector3d normal_global = Pose3d::transform_vector(GLOBAL, normal);
          contacts.clear();
          contacts.push_back(
            CollisionDetection::create_contact(blade_cg,ground_cg,p2_global,normal_global,min_dist)
              );
        }
        else if (p2_plane[Y] == min_dist)
        {
          Point3d p2_global = Pose3d::transform_point(GLOBAL, p2); 
          Vector3d normal(0,1,0,Pplane);
          Vector3d normal_global = Pose3d::transform_vector(GLOBAL, normal);
          contacts.push_back(
            CollisionDetection::create_contact(blade_cg,ground_cg,p2_global,normal_global,min_dist)
              );
        }      
      }
    }

    /// Finds contacts between two collision geometries
    virtual void  find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::vector<UnilateralConstraint>& contacts, double TOL = NEAR_ZERO)
    {
      if (cgA == wheel_cg && cgB == ground_cg)
        find_contacts_blade_plane(cgA, cgB, contacts);
      else if (cgA == ground_cg && cgB == wheel_cg)
        find_contacts_blade_plane(cgB, cgA, contacts);
      else
        ccd->find_contacts(cgA, cgB, contacts, TOL);
    }

    /// Computes signed distance between geometries
    virtual double calc_signed_dist(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, Point3d& pA, Point3d& pB)
    {
      // only handle specific case
      if (cgA == wheel_cg && cgB == ground_cg)
        return calc_signed_dist_blade_plane(cgA, cgB, pA, pB);
      else if (cgA == ground_cg && cgB == wheel_cg)
        return calc_signed_dist_blade_plane(cgB, cgA, pA, pB);
      else
        return ccd->calc_signed_dist(cgA, cgB, pA, pB);
    }

};

extern "C"
{
  boost::shared_ptr<CollisionDetection> factory()
  {
    return boost::shared_ptr<CollisionDetection>(new BladePlanePlugin);
  }
}
#undef NDEBUG

