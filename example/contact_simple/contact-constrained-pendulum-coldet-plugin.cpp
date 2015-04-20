#include <Moby/CollisionDetection.h>
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/CCD.h>
#include <cstdio>
#define NDEBUG
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

class PendulumColdetPlugin : public CollisionDetection
{
  private:
    boost::shared_ptr<ConstraintSimulator> sim;
    boost::shared_ptr<CCD> ccd;
    RigidBodyPtr l1;
    RigidBodyPtr world;
    CollisionGeometryPtr world_cg, l1_cg, right_foot_cg;
    static double sqr(double x) { return x*x; }

  public:
    PendulumColdetPlugin() { ccd = boost::shared_ptr<CCD>(new CCD); }

    virtual void set_simulator(boost::shared_ptr<ConstraintSimulator> sim)
    {
      this->sim = sim;

      // find the necessary objects
      for (unsigned i=0; i< sim->get_dynamic_bodies().size(); i++)
      {
        DynamicBodyPtr body = sim->get_dynamic_bodies()[i];
        if (body->id == "l1")
          l1 = boost::dynamic_pointer_cast<RigidBody>(body);
        else if (body->id == "world")
          world = dynamic_pointer_cast<RigidBody>(body);
      }

      // make sure bodies were found
      if (!l1)
        throw std::runtime_error("Didn't find body with id 'l1'");
      if (!world)
        throw std::runtime_error("Didn't find body with id 'world'");

      // get the geometries
      world_cg = world->geometries.front();
      l1_cg = l1->geometries.front();
    }

    virtual ~PendulumColdetPlugin() {}

    /// Does broad phase collision detection
    virtual void broad_phase(double dt, const std::vector<DynamicBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
    {
      // clear to_check
      to_check.clear();

      // now add in collision geometries for the plane and the walker
      to_check.push_back(std::make_pair(world_cg, l1_cg));
    }

    /// Computes a conservative advancement step for general integration
    virtual double calc_CA_step(const PairwiseDistInfo& pdi)
    {
      return ccd->calc_CA_step(pdi);
    }

    /// Computes a conservative advancement step for Euler integration
    virtual double calc_CA_Euler_step_cv(const PairwiseDistInfo& pdi)
    {
      return ccd->calc_CA_Euler_step_cv(pdi);
    }

    /// Computes a conservative advancement step for Euler integration
    virtual double calc_CA_Euler_step_ca(const PairwiseDistInfo& pdi)
    {
      return ccd->calc_CA_Euler_step_ca(pdi);
    }

    /// Calculates signed distance between a l1 and a plane
    double calc_signed_dist_l1_world(CollisionGeometryPtr l1_cg, CollisionGeometryPtr world_cg)
    {
      Vector3d point, normal;
      return calc_signed_dist_l1_world(l1_cg,world_cg,point,normal);
    }

    /// Calculates signed distance between joint constraint and world 
    double calc_signed_dist_l1_world(CollisionGeometryPtr l1_cg, CollisionGeometryPtr world_cg,Vector3d& pl1,Vector3d& pworld)
    {
      const unsigned X = 0, Y = 1, Z = 2;

      // get the pose for the world 
      shared_ptr<const Pose3d> Pworld = world_cg->get_pose(); 

      // get the pose for the l1 
      shared_ptr<const Pose3d> Pl1 = l1_cg->get_pose();

      // setup the point for the world
      pworld = Vector3d(0,0,0, Pworld);

      // get the point on l1
      pl1 = Vector3d(0,1,0, Pl1);

      // compute the difference
      return -(Pl1->transform_point(Pworld, pl1) - pworld).norm();
    }

    /// Finds contacts between the l1 and a plane
    virtual void find_contacts_l1_world(CollisionGeometryPtr l1_cg, CollisionGeometryPtr world_cg, std::vector<UnilateralConstraint>& contacts)
    {
      const unsigned X = 0, Y = 1, Z = 2;

      contacts.clear();

      // get the point on l1
      shared_ptr<const Pose3d> Pl1 = l1_cg->get_pose();
      Point3d pl1 = Vector3d(0,1,0, Pl1);
      Point3d p = Pl1->transform_point(GLOBAL, pl1);

      // setup contact point
      Point3d point(0,0,0,GLOBAL);

      // setup six normals
      Vector3d n1(0,+1,0,GLOBAL);
      Vector3d n2(0,-1,0,GLOBAL);
      Vector3d n3(0,0,+1,GLOBAL);
      Vector3d n4(0,0,-1,GLOBAL);
      Vector3d n5(+1,0,0,GLOBAL);
      Vector3d n6(-1,0,0,GLOBAL);

      // setup four contacts
      contacts.push_back(CollisionDetection::create_contact(l1_cg,world_cg,point, n1,std::min(0.0, -p[1])));
      contacts.push_back(CollisionDetection::create_contact(l1_cg,world_cg,point, n2,std::min(0.0, -p[1])));
      contacts.push_back(CollisionDetection::create_contact(l1_cg,world_cg,point, n3,std::min(0.0, -p[2])));
      contacts.push_back(CollisionDetection::create_contact(l1_cg,world_cg,point, n4,std::min(0.0, -p[2])));
      contacts.push_back(CollisionDetection::create_contact(l1_cg,world_cg,point, n5,std::min(0.0, -p[0])));
      contacts.push_back(CollisionDetection::create_contact(l1_cg,world_cg,point, n6,std::min(0.0, -p[0])));
    }

    /// Finds contacts between two collision geometries
    virtual void  find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::vector<UnilateralConstraint>& contacts, double TOL = NEAR_ZERO)
    {
      if (cgA == l1_cg && cgB == world_cg)
        find_contacts_l1_world(cgA, cgB, contacts);
      else if (cgA == world_cg && cgB == l1_cg)
        find_contacts_l1_world(cgB, cgA, contacts);
      else
        ccd->find_contacts(cgA, cgB, contacts, TOL);
    }

    /// Computes signed distance between geometries
    virtual double calc_signed_dist(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, Point3d& pA, Point3d& pB)
    {
      // only handle specific case
      if (cgA == l1_cg && cgB == world_cg)
        return calc_signed_dist_l1_world(cgA, cgB, pA, pB);
      else if (cgA == world_cg && cgB == l1_cg)
        return calc_signed_dist_l1_world(cgB, cgA, pA, pB);
    }

};

extern "C"
{
  boost::shared_ptr<CollisionDetection> factory()
  {
    return boost::shared_ptr<CollisionDetection>(new PendulumColdetPlugin);
  }
}
#undef NDEBUG

