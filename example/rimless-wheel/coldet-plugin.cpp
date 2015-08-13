#include <Moby/CollisionDetection.h>
#include <Moby/CCD.h>
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/Log.h>
#include "params.h"
#include <cstdio>
#define NDEBUG
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Moby;

class BladePlanePlugin : public CollisionDetection
{
  private:
    boost::shared_ptr<TimeSteppingSimulator> sim;
    boost::shared_ptr<CCD> ccd;
    RigidBodyPtr wheel;
    RigidBodyPtr ground_body;
    CollisionGeometryPtr ground_cg, wheel_cg, right_foot_cg;

  public:
    BladePlanePlugin() {}

    virtual void set_simulator(boost::shared_ptr<ConstraintSimulator> sim)
    {
      this->sim = boost::dynamic_pointer_cast<TimeSteppingSimulator>(sim);

      // find the necessary objects
      for (unsigned i=0; i< sim->get_dynamic_bodies().size(); i++)
      {
        ControlledBodyPtr body = sim->get_dynamic_bodies()[i];
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
    virtual void broad_phase(double dt, const std::vector<ControlledBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
    {
      // clear to_check
      to_check.clear();

      // remove the wheel from bodies
      std::vector<ControlledBodyPtr> remainder = bodies;
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

    /// Computes a conservative advancement step for Euler integration
    virtual double calc_CA_Euler_step(const PairwiseDistInfo& pdi)
    {
      return ccd->calc_CA_Euler_step(pdi);
    }

    /// Calculates signed distance between a wheel and a plane
    double calc_signed_dist_wheel_plane(CollisionGeometryPtr wheel_cg, CollisionGeometryPtr ground_cg)
    {
      Ravelin::Vector3d point, normal;
      return calc_signed_dist_wheel_plane(wheel_cg,ground_cg,point,normal);
    }

    /// Calculates signed distance between a wheel and a plane
    double calc_signed_dist_wheel_plane(CollisionGeometryPtr wheel_cg, CollisionGeometryPtr ground_cg,Ravelin::Vector3d& pwheel,Ravelin::Vector3d& pground)
    {
      const unsigned Y = 1;

      // setup the minimum dist
      double min_dist = std::numeric_limits<double>::max();

      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(ground_cg->get_geometry());

      // get the pose for the plane primitive
      shared_ptr<const Ravelin::Pose3d> Pplane = plane_geom->get_pose(ground_cg);

      // get the pose for the wheel 
      shared_ptr<const Ravelin::Pose3d> Pwheel = wheel_cg->get_pose();
      FILE_LOG(LOG_COLDET) << "Ravelin::Pose of the wheel: " << Ravelin::Pose3d::calc_relative_pose(Pwheel, GLOBAL) << std::endl;

      // get the pose of the center-of-mass of the wheel
      Ravelin::Transform3d pTb = Ravelin::Pose3d::calc_relative_pose(Pwheel, Pplane);

      // loop over each spoke
      for (unsigned i=0; i< N_SPOKES; i++)
      {
        // get the value of theta
        double theta = M_PI * i * 2.0 / N_SPOKES;

        // setup the two points of each spoke
        Point3d p1(std::cos(theta)*R, W*.5, std::sin(theta)*R, Pwheel);
        Point3d p2(std::cos(theta)*R, -W*.5, std::sin(theta)*R, Pwheel);

        // transform the two points
        Point3d p1_plane = pTb.transform_point(p1);
        Point3d p2_plane = pTb.transform_point(p2);

        FILE_LOG(LOG_COLDET) << "distance between blade " << i << " and ground: " << std::min(p1_plane[Y], p2_plane[Y]) << std::endl;

        // get the signed distance
        if (p1_plane[Y] < min_dist)
        {
          min_dist = p1_plane[Y];
          pground = p1_plane;
          pground[Y] = 0.0;
          pwheel = p1;
        }
        if (p2_plane[Y] < min_dist)
        {
          min_dist = p2_plane[Y];
          pground = p2_plane;
          pground[Y] = 0.0;
          pwheel = p2;
        }
      }

      return min_dist;
    }

    /// Calculates signed distance between the wheel and a curved ground
/*
    double calc_signed_dist_wheel_cos(CollisionGeometryPtr wheel_cg, CollisionGeometryPtr ground_cg,Ravelin::Vector3d& point,Ravelin::Vector3d& normal)
    {
      const unsigned X = 0, Y = 1;

      // setup the minimum dist
      double min_dist = std::numeric_limits<double>::max();

      // get the pose for the ground 
      shared_ptr<const Ravelin::Pose3d> Pground = ground_cg->get_pose();

      // get the pose for the wheel 
      shared_ptr<const Ravelin::Pose3d> Pwheel = wheel_cg->get_pose();

      // get the pose of the center-of-mass of the wheel
      Ravelin::Transform3d gTw = Ravelin::Pose3d::calc_relative_pose(Pwheel, Pground);

      // loop over each spoke
      for (unsigned i=0; i< N_SPOKES; i++)
      {
        // get the value of theta
        double theta = M_PI * i * 2.0 / N_SPOKES;

        // setup the two points of each spoke
        Point3d p1(std::cos(theta)*R, W*.5, std::sin(theta)*R, Pwheel);
        Point3d p2(std::cos(theta)*R, -W*.5, std::sin(theta)*R, Pwheel);

        // transform the two points
        Point3d p1_ground = gTw.transform_point(p1);
        Point3d p2_ground = gTw.transform_point(p2);

        // get the signed distances
        double d1 = p1_ground[Y] - std::cos(p1_ground[X]*0.5);
        double d2 = p1_ground[Y] - std::cos(p2_ground[X]*0.5);

        // get the signed distance
        if (d1 < min_dist)
        {
          min_dist = d1;
          point = p1_ground;

          // the derivative at theta gives the tangent vector to the curve
          // we need a normal to this tangent, which we'll get by rotating
          // 0 1 0 by a rotation matrix around z by theta1. If the vertical
          // (y) component is negative, we rotate around z by theta1+pi 
          normal = Ravelin::Vector3d(std::cos(theta1),-std::sin(theta1),0,Pground);
        }
        if (W > 0.0 && d2 < min_dist)
        {
          min_dist = d2;
          point = p2_ground;
          normal = Ravelin::Vector3d(0,1,0,Pground);
        }
      }

      return min_dist;
    }
*/

    /// Finds contacts between the wheel and a plane
    virtual void find_contacts_wheel_plane(CollisionGeometryPtr wheel_cg, CollisionGeometryPtr ground_cg, std::vector<UnilateralConstraint>& contacts)
    {
      const unsigned Y = 1;
      const double TOL = sim->contact_dist_thresh;

      // delete the IPC token
      if (FIND_MAP)
        std::remove("IPC.token");

      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(ground_cg->get_geometry());

      // get the pose for the plane primitive
      shared_ptr<const Ravelin::Pose3d> Pplane = plane_geom->get_pose(ground_cg);

      // get the pose for the wheel 
      shared_ptr<const Ravelin::Pose3d> Pwheel = wheel_cg->get_pose();

      // get the pose of the center-of-mass of the wheel
      Ravelin::Transform3d pTb = Ravelin::Pose3d::calc_relative_pose(Pwheel, Pplane);

      // say whether spoke = 5 or spoke = 0 found
      bool zero_spoke = false;
      bool five_spoke = false;

      // clear the set of contacts
      contacts.clear();

      // loop over each spoke
      for (unsigned i=0; i< N_SPOKES; i++)
      {
        // get the value of theta
        double theta = M_PI * i * 2.0 / N_SPOKES;

        // setup the two points of each spoke
        Point3d p1(std::cos(theta)*R, W*.5, std::sin(theta)*R, Pwheel);
        Point3d p2(std::cos(theta)*R, -W*.5, std::sin(theta)*R, Pwheel);

        // transform the two points
        Point3d p1_plane = pTb.transform_point(p1);
        Point3d p2_plane = pTb.transform_point(p2);

        // setup the versions right on the plane
        Point3d p1_plane_prime = p1_plane;
        Point3d p2_plane_prime = p2_plane;
        p1_plane_prime[Y] = 0.0;
        p2_plane_prime[Y] = 0.0;

        FILE_LOG(LOG_COLDET) << "distance between blade " << i << " and ground: " << std::min(p1_plane[Y], p2_plane[Y]) << std::endl;

        // check to see whether sub-4 or four-plus spokes were created 
        if (p1_plane[Y] < TOL || p2_plane[Y] < TOL)
        {
          if (i == 0)
            zero_spoke = true;
          else if (i == 5)
            five_spoke = true;
        }

        // create contact for the first point on the blade (if appropriate) 
        if (p1_plane[Y] < TOL)
        {
          Point3d p1x_global = Ravelin::Pose3d::transform_point(GLOBAL, p1);
          Point3d p1y_global = Ravelin::Pose3d::transform_point(GLOBAL, p1_plane_prime);
          Ravelin::Vector3d normal(0,1,0,Pplane);
          Ravelin::Vector3d normal_global = Ravelin::Pose3d::transform_vector(GLOBAL, normal);
          FILE_LOG(LOG_COLDET) << "found contact between blade " << i << " and ground" << std::endl;
          contacts.push_back(
            CollisionDetection::create_contact(wheel_cg,ground_cg,(p1x_global+p1y_global)*0.5,normal_global,p1_plane[Y])
              );
        }

        // create contact for the second point on the blade (if appropriate) 
        if (W > 0.0 && p2_plane[Y] < TOL)
        {
          FILE_LOG(LOG_COLDET) << "found contact between blade " << i << " and ground" << std::endl;
          Point3d p2x_global = Ravelin::Pose3d::transform_point(GLOBAL, p2);
          Point3d p2y_global = Ravelin::Pose3d::transform_point(GLOBAL, p2_plane_prime);
          Ravelin::Vector3d normal(0,1,0,Pplane);
          Ravelin::Vector3d normal_global = Ravelin::Pose3d::transform_vector(GLOBAL, normal);
          contacts.push_back(
            CollisionDetection::create_contact(wheel_cg,ground_cg,(p2x_global+p2y_global)*0.5,normal_global,p2_plane[Y])
              );
        }
      }

      if (FIND_MAP && zero_spoke)
      {
        if (!five_spoke)
          std::cerr << "zero spoke encountered but not five spoke!" << std::endl;
        else
          std::cerr << "zero spoke encountered with five spoke! (Good)" << std::endl;
        std::ofstream out("IPC.token");
        out.close();
      }

      // output constraint violation
      double cvio = std::numeric_limits<double>::max();
      for (unsigned i=0; i< contacts.size(); i++)
        cvio = std::min(cvio, contacts[i].signed_violation);
      std::ofstream out("cvio.dat", std::ostream::app);
      out << cvio << std::endl;
      out.close();
    }

    /// Finds contacts between two collision geometries
    virtual void  find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::vector<UnilateralConstraint>& contacts, double TOL = NEAR_ZERO)
    {
      if (cgA == wheel_cg && cgB == ground_cg)
        find_contacts_wheel_plane(cgA, cgB, contacts);
      else if (cgA == ground_cg && cgB == wheel_cg)
        find_contacts_wheel_plane(cgB, cgA, contacts);
      else
        ccd->find_contacts(cgA, cgB, contacts, TOL);
    }

    /// Computes signed distance between geometries
    virtual double calc_signed_dist(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, Point3d& pA, Point3d& pB)
    {
      // only handle specific case
      if (cgA == wheel_cg && cgB == ground_cg)
        return calc_signed_dist_wheel_plane(cgA, cgB, pA, pB);
      else if (cgA == ground_cg && cgB == wheel_cg)
        return calc_signed_dist_wheel_plane(cgB, cgA, pA, pB);
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

