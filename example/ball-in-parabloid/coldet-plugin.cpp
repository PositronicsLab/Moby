#include <Moby/CollisionDetection.h>
#include <Moby/CCD.h>
#include <Moby/EventDrivenSimulator.h>
#include <cstdio>
#define NDEBUG
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

// set radius for the ball here
const double R = 1.0;
      
// curvature constants for the parabloid
const double A = 10.0, B = 10.0;

class BallParabloidPlanePlugin : public CollisionDetection
{
  private:
    boost::shared_ptr<EventDrivenSimulator> sim;
    boost::shared_ptr<CCD> ccd;
    RigidBodyPtr ball;
    RigidBodyPtr parabloid;
    CollisionGeometryPtr parabloid_cg, ball_cg, right_foot_cg;
    static double sqr(double x) { return x*x; }

  public:
    BallParabloidPlanePlugin() {}

    virtual void set_simulator(boost::shared_ptr<EventDrivenSimulator> sim)
    {
      this->sim = sim;

      // find the necessary objects
      for (unsigned i=0; i< sim->get_dynamic_bodies().size(); i++)
      {
        DynamicBodyPtr body = sim->get_dynamic_bodies()[i];
        if (body->id == "ball")
          ball = boost::dynamic_pointer_cast<RigidBody>(body);
        else if (body->id == "parabloid")
          parabloid = dynamic_pointer_cast<RigidBody>(body);
      }

      // make sure bodies were found
      if (!ball)
        throw std::runtime_error("Didn't find body with id 'ball'");
      if (!parabloid)
        throw std::runtime_error("Didn't find body with id 'parabloid'");

      // get the geometries
      parabloid_cg = parabloid->geometries.front();
      ball_cg = ball->geometries.front();

      // create the continuous collision detection system
      ccd = boost::shared_ptr<CCD>(new CCD);
    }

    virtual ~BallParabloidPlanePlugin() {}

    /// Does broad phase collision detection
    virtual void broad_phase(double dt, const std::vector<DynamicBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
    {
      // clear to_check
      to_check.clear();

      // remove the ball and the parabloid from bodies
      std::vector<DynamicBodyPtr> remainder = bodies;
      for (unsigned i=0; i< remainder.size(); i++)
        if (remainder[i] == ball || remainder[i] == parabloid)
        {
          remainder[i] = remainder.back();
          remainder.pop_back();
          i--;
        }

      // call CCD on the remainder
      ccd->broad_phase(dt, remainder, to_check);

      // now add in collision geometries for the plane and the walker
      to_check.push_back(std::make_pair(parabloid_cg, ball_cg));
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

    /// Calculates signed distance between a ball and a plane
    double calc_signed_dist_ball_parabloid(CollisionGeometryPtr ball_cg, CollisionGeometryPtr parabloid_cg)
    {
      Vector3d point, normal;
      return calc_signed_dist_ball_parabloid(ball_cg,parabloid_cg,point,normal);
    }

    /// Calculates signed distance between a ball and a plane
    double calc_signed_dist_ball_parabloid(CollisionGeometryPtr ball_cg, CollisionGeometryPtr parabloid_cg,Vector3d& pball,Vector3d& pparabloid)
    {
      const unsigned X = 0, Y = 1, Z = 2;

      // setup the minimum dist
      double min_dist = std::numeric_limits<double>::max();

      // get the pose for the parabloid 
      shared_ptr<const Pose3d> Pparabloid = parabloid_cg->get_pose(); 

      // get the pose for the ball 
      shared_ptr<const Pose3d> Pball = ball_cg->get_pose();

      // get the pose of the center-of-mass of the ball
      Transform3d pTb = Pose3d::calc_relative_pose(Pball, Pparabloid);

      // compute the height of the ball above the parabloid; note: this assumes
      // that the parabloid is relatively flat
      double dist = -sqr(pTb.x[X])/(A*A) + sqr(pTb.x[Y])/(B*B) + (pTb.x[Z] - R); 

      // compute the closest point - we'll go half way between the closest
      // point on the ball and the closest point on the surface 
      double z = sqr(pTb.x[X])/sqr(A) + sqr(pTb.x[Y])/sqr(B);
      pparabloid = Point3d(pTb.x[X], pTb.x[Y], z, Pparabloid);
      pball = Point3d(0.0, 0.0, -R, Pball);

      return dist;
    }

    /// Finds contacts between the ball and a plane
    virtual void find_contacts_ball_parabloid(CollisionGeometryPtr ball_cg, CollisionGeometryPtr parabloid_cg, std::vector<UnilateralConstraint>& contacts)
    {
      const unsigned X = 0, Y = 1, Z = 2;

      // get the pose for the parabloid 
      shared_ptr<const Pose3d> Pparabloid = parabloid_cg->get_pose(); 

      // get the pose for the ball 
      shared_ptr<const Pose3d> Pball = ball_cg->get_pose();

      // get the pose of the center-of-mass of the ball
      Transform3d pTb = Pose3d::calc_relative_pose(Pball, Pparabloid);

      // compute the height of the ball above the parabloid; note: this assumes
      // that the parabloid is relatively flat
      double dist = -sqr(pTb.x[X])/(A*A) + sqr(pTb.x[Y])/(B*B) + (pTb.x[Z] - R); 

      // compute the closest point - we'll go half way between the closest
      // point on the ball and the closest point on the surface 
      double z = sqr(pTb.x[X])/sqr(A) + sqr(pTb.x[Y])/sqr(B);
      Point3d closest_parabloid(pTb.x[X], pTb.x[Y], z, Pparabloid);
      Point3d closest_ball(0.0, 0.0, -R, Pball);
      Point3d point = (Pose3d::transform_point(GLOBAL, closest_ball) +
                       Pose3d::transform_point(GLOBAL, closest_parabloid))*0.5; 

      // compute the normal to the surface
      Vector3d normal(2.0*-pTb.x[X]/(A*A), 2.0*-pTb.x[Y]/(B*B), 1.0, Pparabloid); 
      normal.normalize();
      Vector3d normal_global = Pose3d::transform_vector(GLOBAL, normal);

      contacts.clear();
      contacts.push_back(
        CollisionDetection::create_contact(ball_cg,parabloid_cg,point, normal_global,dist)
              );
    }

    /// Finds contacts between two collision geometries
    virtual void  find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::vector<UnilateralConstraint>& contacts, double TOL = NEAR_ZERO)
    {
      if (cgA == ball_cg && cgB == parabloid_cg)
        find_contacts_ball_parabloid(cgA, cgB, contacts);
      else if (cgA == parabloid_cg && cgB == ball_cg)
        find_contacts_ball_parabloid(cgB, cgA, contacts);
      else
        ccd->find_contacts(cgA, cgB, contacts, TOL);
    }

    /// Computes signed distance between geometries
    virtual double calc_signed_dist(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, Point3d& pA, Point3d& pB)
    {
      // only handle specific case
      if (cgA == ball_cg && cgB == parabloid_cg)
        return calc_signed_dist_ball_parabloid(cgA, cgB, pA, pB);
      else if (cgA == parabloid_cg && cgB == ball_cg)
        return calc_signed_dist_ball_parabloid(cgB, cgA, pA, pB);
      else
        return ccd->calc_signed_dist(cgA, cgB, pA, pB);
    }

};

extern "C"
{
  boost::shared_ptr<CollisionDetection> factory()
  {
    return boost::shared_ptr<CollisionDetection>(new BallParabloidPlanePlugin);
  }
}
#undef NDEBUG

