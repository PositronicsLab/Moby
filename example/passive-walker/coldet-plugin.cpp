#include <Moby/CollisionDetection.h>
#include <Moby/CCD.h>
#include <Moby/TimeSteppingSimulator.h>
#define NDEBUG
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Ravelin;
using namespace Moby;

class TorusPlanePlugin : public CCD
{
  private:
    boost::shared_ptr<TimeSteppingSimulator> sim;
    ControlledBodyPtr walker;
    RigidBodyPtr ground_body, left_foot_body, right_foot_body;
    CollisionGeometryPtr ground_cg, left_foot_cg, right_foot_cg;

    protected:
      virtual double calc_next_CA_Euler_step(const PairwiseDistInfo& pdi)
      {
        return std::numeric_limits<double>::max();
      }

  public:
    TorusPlanePlugin() {}

    virtual void set_simulator(boost::shared_ptr<TimeSteppingSimulator> sim)
    {
      this->sim = sim;

      // find the necessary objects
      for (unsigned i=0; i< sim->get_dynamic_bodies().size(); i++)
      {
        ControlledBodyPtr body = sim->get_dynamic_bodies()[i];
        if (body->id == "WALKER")
        {
          walker = body;

          // get the articulated body
          ArticulatedBodyPtr ab = dynamic_pointer_cast<ArticulatedBody>(body);

          // look for left foot and right foot
          for (unsigned j=0; j< ab->get_links().size(); j++)
          {
            if (ab->get_links()[j]->body_id == "LLEG")
              left_foot_body = dynamic_pointer_cast<RigidBody>(ab->get_links()[j]);
            else if (ab->get_links()[j]->body_id == "RLEG")
              right_foot_body = dynamic_pointer_cast<RigidBody>(ab->get_links()[j]);
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
    }

    virtual ~TorusPlanePlugin() {}

    /// Does broad phase collision detection
    virtual void broad_phase(double dt, const std::vector<ControlledBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
    {
      // clear to_check
      to_check.clear();

      // remove the walker from bodies
      std::vector<ControlledBodyPtr> remainder = bodies;
      for (unsigned i=0; i< remainder.size(); i++)
        if (remainder[i] == walker)
        {
          remainder[i] = remainder.back();
          remainder.pop_back();
          break;
        }

      // call CCD on the remainder
      CCD::broad_phase(dt, remainder, to_check);

      // now add in collision geometries for the plane and the walker
      to_check.push_back(std::make_pair(ground_cg, left_foot_cg));
      to_check.push_back(std::make_pair(ground_cg, right_foot_cg));
    }

    /// Computes a conservative advancement step for Euler integration
    virtual double calc_CA_Euler_step(const PairwiseDistInfo& pdi)
    {
      return CCD::calc_CA_Euler_step(pdi);
    }

    double fRand(double fMin, double fMax)
    {
        double f = (double)rand() / RAND_MAX;
        return fMin + f * (fMax - fMin);
    }

    /// Calculates signed distance between a torus and a plane
    double calc_signed_dist_torus_plane(CollisionGeometryPtr torus_cg, CollisionGeometryPtr ground_cg)
    {
      Vector3d point, normal;
      return calc_signed_dist_torus_plane(torus_cg,ground_cg,point,normal);
    }
    /// Calculates signed distance between a torus and a plane
    double calc_signed_dist_torus_plane(CollisionGeometryPtr torus_cg, CollisionGeometryPtr ground_cg,Vector3d& point,Vector3d& normal)
    {
#ifndef NDEBUG
      std::cout << ">> start calc_signed_dist_torus_plane(.)" << std::endl;
      std::cout << "Body: "<<  torus_cg->get_single_body()->id << std::endl;
#endif
      const double R = 0.1236;  // radius from center of tube to center of torus
      const double r = 0.0   ;  // radius of the tube


      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(ground_cg->get_geometry());

      // get the pose for the plane primitive
      shared_ptr<const Pose3d> Pplane = plane_geom->get_pose(ground_cg);
      // get the pose for the torus
      shared_ptr<const Pose3d> Ptorus = torus_cg->get_pose();

      // translate foot into position on leg link
      //   and make z-axis torus orient along y axis
      //   with 90deg x-axis rotation
      Ptorus = shared_ptr<const Pose3d>(
                 new Pose3d(
                   Quatd(Matrix3d(1,0,0,0,0,-1,0,1,0)),
                   Origin3d(0,0,0),
                   Ptorus
                   )
                 );

      // get the transformation from the torus's space to the plane's space
      // (with y-axis up)
      Transform3d tPp = Pose3d::calc_relative_pose(Pplane, Ptorus);

      // Normal of plane (in torus frame)
      // plane rotation matrix
      Matrix3d tRp(tPp.q);
      // Y column of rotaton matrix (plane to torus)
      // is plane normal in torus frame
      Point3d n_plane(tRp.get_column(1),Ptorus);
      n_plane.normalize();

      // Convert to global coords for output
      normal = Pose3d::transform_vector(Moby::GLOBAL,n_plane);

      // Torus axis is z-axis in torus frame
      Vector3d k(0,0,1,Ptorus);

      // Set intitial value of distance to contact
      double d = std::numeric_limits<double>::infinity();

      // if Torus is aligned with plane:
      // Return distance torus origin to
      // closest point on plane less pipe r
      double n_dot_k = n_plane.dot(k);
      if(fabs(n_dot_k) > 1.0-Moby::NEAR_ZERO){
        // d = depth
        // p0 = plane origin, p = plane normal
        // l0 = line origin, l = line direction

        // plane origin: plane origin in torus frame
        // line origin: torus origin in torus frame
        Point3d p0(tPp.x,Ptorus), l0(0,0,0,Ptorus);

        // plane normal: plane normal in torus frame
        // line direction: torus k axis
        Vector3d n = n_plane,l = k;

        // distance torus to closest point on plane is:
        // distance torus origin to closest point on plane
        // - distance torus edge to torus origin
        d = (p0 - l0).dot(n)/(l.dot(n)) - r;
        // Contact point is a random point on the
        // circular manifold of contact
        double t = fRand(-M_PI_2,M_PI_2);
        Point3d p_torus(R*cos(t),R*sin(t),r,Ptorus);
        point = Pose3d::transform_point(Moby::GLOBAL,p_torus);
#ifndef NDEBUG
        std::cout << " -- Torus is parallel to plane"<< std::endl;
        std::cout << "Point: "<<  point << std::endl;
        std::cout << "Normal: "<<  normal << std::endl;
        std::cout << "distance: "<<  d << std::endl;
        std::cout << "<< end calc_signed_dist_torus_plane(.)" << std::endl;
#endif
        return d;
      }

      //((n_plane x axis_torus) x axis_torus)
      Vector3d d_ring = Vector3d::cross(
                          Vector3d::cross(n_plane,k),
                          k
                          );
      d_ring.normalize();

      // if Torus is _|_ with plane:
      // Return distance torus to plane less pipe r and ring R
      if(fabs(n_dot_k) < Moby::NEAR_ZERO){
        // d = depth
        // p0 = plane origin, p = plane normal
        // l0 = line origin, l = line direction

        // plane origin: plane origin in torus frame
        // line origin: torus origin in torus frame
        Point3d p0(tPp.x,Ptorus), l0(0,0,0,Ptorus);

        // plane normal: plane normal in torus frame
        // line direction: on xy-plane of torus
        //   parallel to plane normal in torus frame
        Vector3d n = n_plane,l = d_ring;
        d = (p0 - l0).dot(n)/(l.dot(n)) - (r+R);
        Point3d p_torus = d*l + l0;
        point = Pose3d::transform_point(Moby::GLOBAL,p_torus);
#ifndef NDEBUG
        std::cout << " -- Torus is perpendicular to plane"<< std::endl;
        std::cout << "Point: "<<  point << std::endl;
        std::cout << "Normal: "<<  normal << std::endl;
        std::cout << "distance: "<<  d << std::endl;
        std::cout << "<< end calc_signed_dist_torus_plane(.)" << std::endl;
#endif
        return d;
      }

      //   ((d_ring x axis_torus) x n_plane ) x (d_ring x axis_torus)
      //    ^ tangent to pipe   ^
      //   ^ _|_ to plane normal            ^
      //   ^ toward plane on torus pipe                             ^
      Vector3d d_pipe = Vector3d::cross(
                          Vector3d::cross(Vector3d::cross(d_ring,k),n_plane),
                          Vector3d::cross(-d_ring,k)
                          );
      d_pipe.normalize();


      // d = depth
      // p0 = plane origin, p = plane normal
      // l0 = line origin, l = line direction

      // plane origin: plane origin in torus frame
      // line origin: torus origin in torus frame
      Point3d p0(tPp.x,Ptorus), l0 = R * d_ring;

      // plane normal: plane normal in torus frame
      // line direction: on xy-plane of torus
      //   parallel to plane normal in torus frame
      Vector3d n = n_plane,l = d_pipe;
      d = (p0 - l0).dot(n)/(l.dot(n)) - r;

      //point on torus closest to plane;
      Point3d p_torus = R * d_ring + r * d_pipe;
      p_torus.pose = Ptorus;
      // TODO: find the point in the torus's space such that
      //       tPp.transform_point(.) results in the value of y closest to
      //       negative infinity
      point = Pose3d::transform_point(Moby::GLOBAL,p_torus);
#ifndef NDEBUG
      std::cout << "Point: "<<  point << std::endl;
      std::cout << "Normal: "<<  normal << std::endl;
      std::cout << "distance: "<<  d << std::endl;
      std::cout << "<< end calc_signed_dist_torus_plane(.)" << std::endl;
#endif
      return d;
#ifndef NDEBUG
      std::cout << "<< end calc_signed_dist_torus_plane(.)" << std::endl;
#endif

    }

    /// Finds contacts between a torus and a plane
    virtual void find_contacts_torus_plane(CollisionGeometryPtr torus_cg, CollisionGeometryPtr ground_cg, std::vector<UnilateralConstraint>& contacts)
    {
      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(ground_cg->get_geometry());

      // TODO: find the point in the torus's space such that
      //       tPp.transform_point(.) results in the value of y closest to
      //       negative infinity. If that value > 0.0, there is no contact.

      //

      // contact normal is always going to be [0, 0, 1] in the plane's frame,
      // assuming that z is up
      Vector3d point,normal;
      double violation = calc_signed_dist_torus_plane(torus_cg,ground_cg,point,normal);

      if(violation <= 0.0)
        contacts.push_back(
              CollisionDetection::create_contact(torus_cg,ground_cg,point,normal,violation)
              );

      // TODO: call CollisionDetection::create_contact(.) to create the actual contact
//      if(contacts.size() == 0)
//        contacts.push_back(
//              CollisionDetection::create_contact(torus_cg,ground_cg,point,normal,violation)
//              );
//      else if(contacts.size() >= 2)
//        return;
//      else if( (contacts[0].contact_geom1->get_single_body()->id.compare("LLEG") == 0
//            && torus_cg->get_single_body()->id.compare("RLEG") == 0)
//          || (contacts[0].contact_geom1->get_single_body()->id.compare("RLEG") == 0
//              && torus_cg->get_single_body()->id.compare("LLEG") == 0))
//        contacts.push_back(
//              CollisionDetection::create_contact(torus_cg,ground_cg,point,normal,violation)
//              );
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
        CCD::find_contacts(cgA, cgB, contacts, TOL);
    }

    /// Computes signed distance between geometries
    virtual double calc_signed_dist(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, Point3d& pA, Point3d& pB)
    {
      // only handle specific case
      if (cgA == left_foot_cg && cgB == ground_cg)
        return calc_signed_dist_torus_plane(cgA, cgB, pA, pB);
      else if (cgA == ground_cg && cgB == left_foot_cg)
        return calc_signed_dist_torus_plane(cgB, cgA, pA, pB);
      else if (cgA == right_foot_cg && cgB == ground_cg)
        return calc_signed_dist_torus_plane(cgA, cgB, pA, pB);
      else if (cgA == ground_cg && cgB == right_foot_cg)
        return calc_signed_dist_torus_plane(cgB, cgA, pA, pB);
      else
        return CCD::calc_signed_dist(cgA, cgB, pA, pB);
    }
};

extern "C"
{
  boost::shared_ptr<CollisionDetection> factory()
  {
    return boost::shared_ptr<CollisionDetection>(new TorusPlanePlugin);
  }
}
#undef NDEBUG
