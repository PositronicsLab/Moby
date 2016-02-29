#include <Moby/CollisionDetection.h>
#include <Moby/CCD.h>
#include <Moby/TimeSteppingSimulator.h>
#include <Moby/Log.h>
#include <cstdio>
#include <stdlib.h>
#define NDEBUG
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using namespace Moby;
using namespace Ravelin;

struct Ellipse
{
  Ellipse() { x = 1.0; y = 2.0; z = 1.0; }
  double x, y, z;
};

Ellipse e; 

class EllipsePlanePlugin : public CollisionDetection
{
  private:
    boost::shared_ptr<TimeSteppingSimulator> sim;
    boost::shared_ptr<CCD> ccd;
    RigidBodyPtr ellipse;
    RigidBodyPtr ground_body, wall_pos_body, wall_neg_body;
    CollisionGeometryPtr ground_cg, ellipse_cg, wall_pos_cg, wall_neg_cg;

  public:
    EllipsePlanePlugin() {}

    virtual double calc_CA_Euler_step(const PairwiseDistInfo& pdi)
    {
      if (pdi.dist <= 0.0)
        return std::numeric_limits<double>::max();
      else return ccd->calc_CA_Euler_step(pdi);
    }

    virtual void set_simulator(boost::shared_ptr<ConstraintSimulator> sim)
    {
      this->sim = boost::dynamic_pointer_cast<TimeSteppingSimulator>(sim);

      // find the necessary objects
      for (unsigned i=0; i< sim->get_dynamic_bodies().size(); i++)
      {
        ControlledBodyPtr body = sim->get_dynamic_bodies()[i];
        if (body->id == "ellipse")
          ellipse = boost::dynamic_pointer_cast<RigidBody>(body);
        else if (body->id == "ground")
          ground_body = dynamic_pointer_cast<RigidBody>(body);
        else if (body->id == "wall_pos")
          wall_pos_body = dynamic_pointer_cast<RigidBody>(body);
        else if (body->id == "wall_neg")
          wall_neg_body = dynamic_pointer_cast<RigidBody>(body);
      }

      // make sure bodies were found
      if (!ellipse)
        throw std::runtime_error("Didn't find body with id 'ellipse'");
      if (!ground_body)
        throw std::runtime_error("Didn't find body with id 'ground'");
      if (!wall_pos_body)
        throw std::runtime_error("Didn't find body with id 'wall_pos'");
      if (!wall_neg_body)
        throw std::runtime_error("Didn't find body with id 'wall_neg'");

      // get the geometries
      ground_cg = ground_body->geometries.front();
      ellipse_cg = ellipse->geometries.front();
      wall_pos_cg = wall_pos_body->geometries.front();
      wall_neg_cg = wall_neg_body->geometries.front();

      // create the continuous collision detection system
      ccd = boost::shared_ptr<CCD>(new CCD);
    }

    virtual ~EllipsePlanePlugin() {}

    /// Does broad phase collision detection
    virtual void broad_phase(double dt, const std::vector<ControlledBodyPtr>& bodies, std::vector<std::pair<CollisionGeometryPtr, CollisionGeometryPtr> >& to_check)
    {
      // clear to_check
      to_check.clear();

      // remove the ellipse from bodies
      std::vector<ControlledBodyPtr> remainder = bodies;
      for (unsigned i=0; i< remainder.size(); i++)
        if (remainder[i] == ellipse)
        {
          remainder[i] = remainder.back();
          remainder.pop_back();
          break;
        }

      // call CCD on the remainder
      ccd->broad_phase(dt, remainder, to_check);

      // now add in collision geometries for the planes and the ellipse 
      to_check.push_back(std::make_pair(ground_cg, ellipse_cg));
      to_check.push_back(std::make_pair(wall_pos_cg, ellipse_cg));
      to_check.push_back(std::make_pair(wall_neg_cg, ellipse_cg));
    }

    /// Calculates signed distance between a ellipse and a plane
    double calc_signed_dist_ellipse_plane(CollisionGeometryPtr ellipse_cg, CollisionGeometryPtr ground_cg)
    {
      Ravelin::Vector3d point, normal;
      return calc_signed_dist_ellipse_plane(ellipse_cg,ground_cg,point,normal);
    }

    /// Evaluates the point in direction d lying on the ellipse 
    double eval_point(const Ellipse& e, const VectorNd& x, const Quatd& q, const Vector3d& d)
    {
      return eval_point(e, x[0], x[1], q, d);
    }

    /// Evaluates the point in direction d lying on the ellipse 
    double eval_point(const Ellipse& e, double theta, double phi, const Quatd& q, const Vector3d& d)
    {
      // setup Cartesian coordinates
      double x = e.x * std::cos(theta) * std::sin(phi);
      double y = e.y * std::sin(theta) * std::sin(phi);
      double z = e.z * std::cos(phi);

      // setup the vector, rotate it, and compute the dot product
      return d.dot(Vector3d(q * Origin3d(x, y, z), GLOBAL));
    } 

    /// Does a grid search looking for the best values of theta and phi
    double grid_search(const Ellipse& e, const Quatd& q, const Vector3d& d, double& theta, double& phi)
    {
      const double PI2 = 2.0 * M_PI;
      const double GRAN = 0.1;

      // setup the minimum value
      double min_value = std::numeric_limits<double>::max();

      // setup best theta and best phi
      double best_theta = 0.0, best_phi = 0.0;

      // loop over values for theta and phi
      for (theta = 0.0; theta < PI2; theta += GRAN)
        for (phi = 0.0; phi < M_PI; phi += GRAN)
        {
          double value = eval_point(e, theta, phi, q, d);
          if (value < min_value)
          {
            min_value = value;
            best_theta = theta;
            best_phi = phi;
          }
        }

      // setup best theta / phi
      theta = best_theta;
      phi = best_phi;

      return min_value;
    }

    /// Polishes a grid search looking for the best values of theta and phi
    double polish(double f0, const Ellipse& e, const Quatd& q, const Vector3d& d, double& theta, double& phi)
    {
      VectorNd grad, dx;
      MatrixNd H, Hprime;
      const double DT = 1e-6, EPS = 1e-8;

      // set x
      VectorNd x(2);
      x[0] = theta;
      x[1] = phi;

      while (true)
      {
        // compute the current gradient and Hessian 
        ngrad(x, DT, e, q, d, grad);
        nhess(x, DT, e, q, d, H);

        // solve H*dx = -g
        dx = grad;
        dx.negate();
        double eps = std::numeric_limits<double>::epsilon();
        Hprime = H;
        while (true)
        {
          bool success = LinAlgd::factor_chol(Hprime);
          if (success)
            break;
          Hprime = H;
          for (unsigned i=0; i< H.rows(); i++)
            Hprime(i,i) += eps;
          eps *= 10.0;
        }
        LinAlgd::solve_chol_fast(Hprime, dx);

        // compute lambda^2 and see whether the stopping criterion is met
        double grad_dot_dx = grad.dot(dx);
        const double LAMBDA_SQ = -grad_dot_dx;
        if (LAMBDA_SQ/2.0 <= EPS)
          break;

        // compute fprime
        (_workv = dx) += x;
        double fprime = eval_point(e, _workv, q, d);

        // do backtracking line search
        double s = 1.0;
        const double ALPHA = 0.05, BETA = 0.5;
        while (fprime > f0 + ALPHA * s * grad_dot_dx)
        {
          // update s and x'
          s *= BETA;
          ((_workv = dx) *= s) += x;
          fprime = eval_point(e, _workv, q, d);
        }

        // verify improvement was possible
        if (f0 - fprime < EPS)
          break;

        // update x
        x = _workv; 

        // update f0
        f0 = fprime;
      }

      theta = x[0];
      phi = x[1];
      return f0;
    }

    /// Calculates signed distance between a ellipse and a plane
    double calc_signed_dist_ellipse_plane(CollisionGeometryPtr ellipse_cg, CollisionGeometryPtr plane_cg,Ravelin::Vector3d& pellipse,Ravelin::Vector3d& pplane)
    {
      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(plane_cg->get_geometry());

      // get the pose for the plane primitive
      shared_ptr<const Ravelin::Pose3d> Pplane = plane_geom->get_pose(plane_cg);

      // get the pose for the ellipse 
      shared_ptr<const Ravelin::Pose3d> Pellipse = ellipse_cg->get_pose();

      // get the normal to the plane (in the world frame)
      Vector3d normalP = Pose3d::transform_vector(GLOBAL, Vector3d(0.0, 1.0, 0.0, Pplane));
      
      // get the pose of the center-of-mass of the ellipse
      Ravelin::Transform3d wTb = Ravelin::Pose3d::calc_relative_pose(Pellipse, GLOBAL);

      // do a grid search to find the closest point on the ellipse to the
      // negation of the normal
      double phi, theta;
      double f0 = grid_search(e, wTb.q, normalP, theta, phi);

      // then do Newton's method to polish it up
      double fval = polish(f0, e, wTb.q, normalP, theta, phi);

      // construct the point on the ellipse in its local frame
      double x = e.x * std::cos(theta) * std::sin(phi);
      double y = e.y * std::sin(theta) * std::sin(phi);
      double z = e.z * std::cos(phi);
      pellipse[0] = x;
      pellipse[1] = y;
      pellipse[2] = z;
      pellipse.pose = Pellipse;

      // get this point transformed to the plane
      pplane = Pose3d::transform_point(Pplane, pellipse);

      // get the distance
      double dist = pplane[1];

      // set the Y value for the point on the plane to zero
      pplane[1] = 0.0;

      return dist;
    }

    /// Finds contacts between the ellipse and a plane
    virtual void find_contacts_ellipse_plane(CollisionGeometryPtr ellipse_cg, CollisionGeometryPtr plane_cg, std::vector<UnilateralConstraint>& contacts)
    {
      // get the plane primitive
      PrimitivePtr plane_geom = dynamic_pointer_cast<Primitive>(plane_cg->get_geometry());

      // get the pose for the plane primitive
      shared_ptr<const Ravelin::Pose3d> Pplane = plane_geom->get_pose(plane_cg);

      // get the pose for the ellipse 
      shared_ptr<const Ravelin::Pose3d> Pellipse = ellipse_cg->get_pose();

      // get the normal to the plane (in the world frame)
      Vector3d normalP = Pose3d::transform_vector(GLOBAL, Vector3d(0.0, 1.0, 0.0, Pplane));
      
      // get the pose of the center-of-mass of the ellipse
      Ravelin::Transform3d wTb = Ravelin::Pose3d::calc_relative_pose(Pellipse, GLOBAL);

      // do a grid search to find the closest point on the ellipse to the
      // negation of the normal
      double phi, theta;
      double f0 = grid_search(e, wTb.q, normalP, theta, phi);

      // then do Newton's method to polish it up
      polish(f0, e, wTb.q, normalP, theta, phi);

      // construct the point on the ellipse in its local frame
      double x = e.x * std::cos(theta) * std::sin(phi);
      double y = e.y * std::sin(theta) * std::sin(phi);
      double z = e.z * std::cos(phi);
      Point3d pellipse(Pellipse);
      pellipse[0] = x;
      pellipse[1] = y;
      pellipse[2] = z;

      // get this point transformed to the plane
      Point3d pplane = Pose3d::transform_point(Pplane, pellipse);

      // get the distance
      double dist = pplane[1];

      // look for contact
      if (dist > 0.0)
        return;

      // compute the contact point 
      Point3d cp = wTb.transform_point(pellipse); 

      // setup the contact point
      contacts.push_back(
        CollisionDetection::create_contact(ellipse_cg,plane_cg,cp,normalP,dist));
    }

    /// Finds contacts between two collision geometries
    virtual void  find_contacts(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::vector<UnilateralConstraint>& contacts, double TOL = NEAR_ZERO)
    {
      // updates the ellipse
      e.x = std::atoi(getenv("ELLIPSE_X_RAD"));
      e.y = std::atoi(getenv("ELLIPSE_Y_RAD"));

      if (cgA == ellipse_cg && cgB == ground_cg)
        find_contacts_ellipse_plane(cgA, cgB, contacts);
      else if (cgA == ellipse_cg && cgB == wall_pos_cg)
        find_contacts_ellipse_plane(cgA, cgB, contacts);
      else if (cgA == ellipse_cg && cgB == wall_neg_cg)
        find_contacts_ellipse_plane(cgA, cgB, contacts);
      else if (cgA == ground_cg && cgB == ellipse_cg)
        find_contacts_ellipse_plane(cgB, cgA, contacts);
      else if (cgA == wall_pos_cg && cgB == ellipse_cg)
        find_contacts_ellipse_plane(cgB, cgA, contacts);
      else if (cgA == wall_neg_cg && cgB == ellipse_cg)
        find_contacts_ellipse_plane(cgB, cgA, contacts);
      else
        ccd->find_contacts(cgA, cgB, contacts, TOL);
    }

    /// Computes signed distance between geometries
    virtual double calc_signed_dist(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, Point3d& pA, Point3d& pB)
    {
      // updates the ellipse
      e.x = std::atoi(getenv("ELLIPSE_X_RAD"));
      e.y = std::atoi(getenv("ELLIPSE_Y_RAD"));

      // only handle specific case
      if (cgA == ellipse_cg && cgB == ground_cg)
        return calc_signed_dist_ellipse_plane(cgA, cgB, pA, pB);
      else if (cgA == ellipse_cg && cgB == wall_pos_cg)
        return calc_signed_dist_ellipse_plane(cgA, cgB, pA, pB);
      else if (cgA == ellipse_cg && cgB == wall_neg_cg)
        return calc_signed_dist_ellipse_plane(cgA, cgB, pA, pB);
      else if (cgA == ground_cg && cgB == ellipse_cg)
        return calc_signed_dist_ellipse_plane(cgB, cgA, pA, pB);
      else if (cgA == wall_pos_cg && cgB == ellipse_cg)
        return calc_signed_dist_ellipse_plane(cgB, cgA, pA, pB);
      else if (cgA == wall_neg_cg && cgB == ellipse_cg)
        return calc_signed_dist_ellipse_plane(cgB, cgA, pA, pB);
      else
        return ccd->calc_signed_dist(cgA, cgB, pA, pB);
    }

  private:
    VectorNd _workv, _workv2, _workv3;

    /// Computes the gradient of a function numerically
    VectorNd& ngrad(const VectorNd& x, double h, const Ellipse& e, const Quatd& q, const Vector3d& d, VectorNd& grad)
    {
      unsigned n = x.size();
      const double inv_h2 = 0.5 / h;

      // make a copy of x
      _workv = x;

      // compute the numerical gradient - dumb way
      grad.resize(n);
      for (unsigned i=0; i< n; i++)
      {
        _workv[i] += h;
        double v1 = eval_point(e, _workv, q, d); 
        _workv[i] = x[i];
        _workv[i] -= h;
        double v2 = eval_point(e, _workv, q, d);
        _workv[i] = x[i];
        grad[i] = (v1 - v2)*inv_h2;
      }

      return grad;
    }

    /// Computes the hessian of a function numerically
    MatrixNd& nhess(const VectorNd& x, double h, const Ellipse& e, const Quatd& q, const Vector3d& d, MatrixNd& hess)
    {
      unsigned n = x.size();
      const double inv_h = 1.0 / h;

      // make a copy of x
      _workv = x;

      hess.resize(n,n);
      for (unsigned i=0; i< n; i++)
      {
        _workv[i] += h;
        ngrad(_workv, h, e, q, d, _workv2);
        _workv[i] = x[i];
        _workv[i] -= h;
        ngrad(_workv, h, e, q, d, _workv3);
        _workv[i] = x[i];
        hess.set_column(i, ((_workv2 -= _workv3)*=inv_h));
      }

      // average values of the Hessian
      const unsigned LD = hess.leading_dim();
      double* dd = hess.data();
      for (unsigned i=0; i< n; i++)
        for (unsigned j=i+1; j< n; j++)
          dd[j+LD*i] = dd[i*LD+j] = 0.5*(dd[j+LD*i] + dd[i+LD*j]);

      return hess;
    }
};

extern "C"
{
  boost::shared_ptr<CollisionDetection> factory()
  {
    return boost::shared_ptr<CollisionDetection>(new EllipsePlanePlugin);
  }
}
#undef NDEBUG

