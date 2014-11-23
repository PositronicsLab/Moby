#include <Moby/CollisionDetection.h>
#include <Moby/CCD.h>
#include <Moby/EventDrivenSimulator.h>

using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

class TorusPlanePlugin : public CollisionDetection
{
  private:
    boost::shared_ptr<EventDrivenSimulator> sim;
    boost::shared_ptr<CCD> ccd;
    DynamicBodyPtr walker;
    RigidBodyPtr ground_body, left_foot_body, right_foot_body;
    CollisionGeometryPtr ground_cg, left_foot_cg, right_foot_cg;

  public:
    TorusPlanePlugin() {}

    virtual void set_simulator(boost::shared_ptr<EventDrivenSimulator> sim)
    {
      this->sim = sim;

      // find the necessary objects
      for (unsigned i=0; i< sim->get_dynamic_bodies().size(); i++)
      {
        DynamicBodyPtr body = sim->get_dynamic_bodies()[i];
        if (body->id == "WALKER")
        {
          walker = body;

          // get the articulated body
          ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(body);

          // look for left foot and right foot
          for (unsigned j=0; j< ab->get_links().size(); j++)
          {
            if (ab->get_links()[j]->id == "LLEG")
              left_foot_body = ab->get_links()[j];
            else if (ab->get_links()[j]->id == "RLEG")
              right_foot_body = ab->get_links()[j];
          }
        }
        else if (body->id == "GROUND")
          ground_body = dynamic_pointer_cast<RigidBody>(body);
      }

      // make sure bodies were found
      if (!walker)
        throw std::runtime_error("Didn't find body with id 'walker'");
      if (!ground_body)
        throw std::runtime_error("Didn't find body with id 'ground'");
      if (!left_foot_body)
        throw std::runtime_error("Didn't find body with id 'left_foot'");
      if (!right_foot_body)
        throw std::runtime_error("Didn't find body with id 'right_foot'");

      // get the geometries
      ground_cg = ground_body->geometries.front();
      left_foot_cg = left_foot_body->geometries.front();
      right_foot_cg = right_foot_body->geometries.front();

      // create the continuous collision detection system
      ccd = boost::shared_ptr<CCD>(new CCD);
    }

    virtual ~TorusPlanePlugin() {}

    /// Does broad phase collision detection
    virtual void broad_phase(double dt, const std::vector<DynamicBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
    {
      // clear to_check
      to_check.clear();

      // remove the walker from bodies
      std::vector<DynamicBodyPtr> remainder = bodies;
      for (unsigned i=0; i< remainder.size(); i++)
        if (remainder[i] == walker)
        {
          remainder[i] = remainder.back();
          break;
        }

      // call CCD on the remainder
      ccd->broad_phase(dt, remainder, to_check);

      // now add in collision geometries for the plane and the walker
      to_check.push_back(std::make_pair(ground_cg, left_foot_cg));
      to_check.push_back(std::make_pair(ground_cg, right_foot_cg));
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

    double fRand(double fMin, double fMax)
    {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    }

    /// Calculates signed distance between a torus and a plane
    double calc_signed_dist_torus_plane(CollisionGeometryPtr torus_cg, CollisionGeometryPtr ground_cg,Vector3d point = Vector3d(),Vector3d normal = Vector3d())
    {
      const double R = 0.1236;  // radius from center of tube to center of torus
      const double r = 0.0   ;  // radius of the tube

      shared_ptr<const Pose3d>
          Ptorus2 = shared_ptr<const Pose3d>(
                     new Pose3d(
                       Quatd::identity(),
                       Origin3d(0,0,-1),
                       torus_cg->get_single_body()->get_pose()
                       )
                     );

      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(ground_cg->get_geometry());

      // get the pose for the plane primitive
      shared_ptr<const Pose3d> Pplane = plane_geom->get_pose(ground_cg);
      // get the pose for the torus
      shared_ptr<const Pose3d> Ptorus = torus_cg->get_pose();

      Ptorus = shared_ptr<const Pose3d>(
                 new Pose3d(
                   Quatd::identity(),
                   Origin3d(0,0,-1),
                   Ptorus
                   )
                 );
      // get the transformation from the torus's space to the plane's space
      // (with y-axis up)
      Transform3d tPp = Pose3d::calc_relative_pose(Pplane, Ptorus);
      Transform3d pPt = Pose3d::calc_relative_pose(Ptorus, Pplane);

      // equation of a z-axis torus
      // f(x,y,z) = (R - sqrt(x^2 + y^2))^2 + z^2 - r^2

      // Equation for an arbitrary plane
      // g(x,y,z) = a (x - x0) + b (y - y0) + c (z - z0)
      // interpenetration:

      // Normal of plane (in torus frame)
      // plane rotation matrix
      Matrix3d tRp(tPp.q);

      // Y column of rotaton matrix is plane normal in torus frame
      Point3d n_plane(tRp.get_column(1),Ptorus);
      normal = Pose3d::transform_vector(Moby::GLOBAL,n_plane);

      // Torus axis is k_hat
      Vector3d k(0,0,1);

      // if Torus is aligned with plane:
      // Return distance torus to plane less pipe r
      double n_dot_t = n_plane.dot(k);
      if(fabs(n_dot_t) > 1.0-Moby::NEAR_ZERO){
        double t = fRand(-M_PI_2,M_PI_2);
        Point3d p_torus(R*cos(t),R*sin(t),r,Ptorus);
        point = Pose3d::transform_point(Moby::GLOBAL,p_torus);
        return (tPp.x.norm() - r);
      }

      //((n_plane x axis_torus) x axis_torus)
      Vector3d d_ring = Vector3d::cross(
                          Vector3d::cross(n_plane,k),
                          k
                          );
      d_ring.normalize();

      // if Torus is _|_ with plane:
      // Return distance torus to plane less pipe r and ring R
      if(fabs(n_dot_t) < Moby::NEAR_ZERO){
        Point3d p_torus = (R+r) * d_ring;
        p_torus.pose = Ptorus;
        point = Pose3d::transform_point(Moby::GLOBAL,p_torus);
        return (tPp.x.norm() - r - R);
      }


      //   ((d_ring x axis_torus) x n_plane ) x (d_ring x axis_torus)
      //    ^ tangent to pipe   ^
      //   ^ _|_ to plane normal            ^
      //   ^ toward plane on torus pipe                             ^
      Vector3d d_pipe = Vector3d::cross(
                          Vector3d::cross(Vector3d::cross(d_ring,k),n_plane),
                          Vector3d::cross(d_ring,k)
                          );
      d_pipe.normalize();

      //point on torus closest to plane;
      Point3d p_torus = R * d_ring + r * d_pipe;
      p_torus.pose = Ptorus;
      // TODO: find the point in the torus's space such that
      //       tPp.transform_point(.) results in the value of y closest to
      //       negative infinity
      std::cout << p_torus << std::endl;
      point = Pose3d::transform_point(Moby::GLOBAL,p_torus);

      // y element of penetration
      return pPt.transform_point(p_torus)[1];
    }

    /// Finds contacts between a torus and a plane
    virtual void find_contacts_torus_plane(CollisionGeometryPtr torus_cg, CollisionGeometryPtr ground_cg, std::vector<UnilateralConstraint>& contacts)
    {
      contacts.clear();
      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(ground_cg->get_geometry());

      // TODO: find the point in the torus's space such that
      //       tPp.transform_point(.) results in the value of y closest to
      //       negative infinity. If that value > 0.0, there is no contact.

      //

      // contact normal is always going to be [0, 0, 1] in the plane's frame,
      // assuming that z is up
//      Vector3d normal(0.0, 0.0, 1.0, plane_geom->get_pose(ground_cg));
      Vector3d point,normal;
      double violation = calc_signed_dist_torus_plane(torus_cg,ground_cg,point,normal);

      // TODO: call CollisionDetection::create_contact(.) to create the actual contact
      contacts.push_back(
            CollisionDetection::create_contact(torus_cg,ground_cg,point,normal,violation)
            );
    }

    /// Finds contacts between two collision geometries
    virtual void  find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::vector<UnilateralConstraint>& contacts, double TOL = NEAR_ZERO)
    {
      if (cgA == left_foot_cg && cgB == ground_cg)
        find_contacts_torus_plane(cgA, cgB, contacts);
      else if (cgA == ground_cg && cgB == left_foot_cg)
        find_contacts_torus_plane(cgB, cgA, contacts);
      else if (cgA == right_foot_cg && cgB == ground_cg)
        find_contacts_torus_plane(cgA, cgB, contacts);
      else if (cgA == ground_cg && cgB == right_foot_cg)
        find_contacts_torus_plane(cgB, cgA, contacts);
      else
        ccd->find_contacts(cgA, cgB, contacts, TOL);
    }

    /// Computes signed distance between geometries
    virtual double calc_signed_dist(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB)
    {
      // only handle specific case
      if (cgA == left_foot_cg && cgB == ground_cg)
        calc_signed_dist_torus_plane(cgA, cgB);
      else if (cgA == ground_cg && cgB == left_foot_cg)
        calc_signed_dist_torus_plane(cgB, cgA);
      else if (cgA == right_foot_cg && cgB == ground_cg)
        calc_signed_dist_torus_plane(cgA, cgB);
      else if (cgA == ground_cg && cgB == right_foot_cg)
        calc_signed_dist_torus_plane(cgB, cgA);
      else
        ccd->calc_signed_dist(cgA, cgB);
    }

};

extern "C" boost::shared_ptr<TorusPlanePlugin> factory()
{
  return boost::shared_ptr<TorusPlanePlugin>(new TorusPlanePlugin);
}

