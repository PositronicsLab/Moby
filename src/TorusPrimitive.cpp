/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifdef USE_OSG
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#endif
#include <Moby/XMLTree.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/BoundingSphere.h>
#include <Moby/Constants.h>
#include <Moby/CollisionGeometry.h>
#include <Moby/HeightmapPrimitive.h>
#include <Moby/TorusPrimitive.h>
#include <Moby/PlanePrimitive.h>
#include "OsgTorus.h"

using boost::shared_ptr;
using boost::const_pointer_cast;
using boost::dynamic_pointer_cast;
using std::vector;
using std::map;
using namespace Ravelin;
using namespace Moby;

TorusPrimitive::TorusPrimitive()
{
  // setup torus parameters to some defaults; this gives an aspect ratio
  // (major radius / minor radius) of ~2.5, which yields a donut
  _major_radius = 2.5;
  _minor_radius = 1.0;

  // calculate the mass properties
  calc_mass_properties();
}

TorusPrimitive::TorusPrimitive(const Pose3d& P) : Primitive(P)
{
  // setup torus parameters to some defaults; this gives an aspect ratio
  // (major radius / minor radius) of ~2.5, which yields a donut
  _major_radius = 2.5;
  _minor_radius = 1.0;

  // calculate the mass properties
  calc_mass_properties();
}

TorusPrimitive::TorusPrimitive(double major_radius, double minor_radius)
{
  // setup torus parameters to some defaults; this gives an aspect ratio
  // (major radius / minor radius) of ~2.5, which yields a donut
  _major_radius = major_radius;
  _minor_radius = minor_radius;

  // calculate the mass properties
  calc_mass_properties();
}

/// Creates the visualization for this primitive
osg::Node* TorusPrimitive::create_visualization()
{
  #ifdef USE_OSG
  const double HALF = (double) 0.5;
  osg::OsgTorus* torus = new osg::OsgTorus((float) _minor_radius, (float) _major_radius);
  osg::Geode* geode = torus->operator()();
  return geode;
  #else
  return NULL;
  #endif
}

/// Gets the root BVH for this torus
BVPtr TorusPrimitive::get_BVH_root(CollisionGeometryPtr geom)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the pointer to the obb 
  OBBPtr& obb = _obbs[geom];

  // create the obb, if necessary
  if (!obb)
  {
    // create the obb 
    obb = OBBPtr(new OBB);
    obb->geom = geom;

    // get the pose for the primitive and geometries
    shared_ptr<const Pose3d> P = _cg_poses[geom];

    // setup the obb center and orientation 
    obb->center = Point3d(0.0, 0.0, 0.0, P);
    obb->R.set_identity();

    // setup OBB half-lengths
    obb->l[X] = _major_radius + _minor_radius;
    obb->l[Y] = _major_radius + _minor_radius;
    obb->l[Z] = _minor_radius;
  }

  return obb;
}

/// Computes the distance from a point to this torus
double TorusPrimitive::calc_dist_and_normal(const Point3d& p, std::vector<Vector3d>& normals) const
{
}

/// Gets the vertices corresponding to this torus
void TorusPrimitive::get_vertices(shared_ptr<const Pose3d> P, vector<Point3d>& vertices) const
{
  // setup some reasonable defaults
  const unsigned CIRCLE_SWEEPS = 20;
  const unsigned TORUS_SWEEPS = 20;
  const double START_SWEEP = 0.0;
  const double END_SWEEP = M_PI * 2.0;

  // clear vertices vector
  vertices.clear();

  // verify that the pose is found
  assert(_poses.find(const_pointer_cast<Pose3d>(P)) != _poses.end());

  // TODO: note that this does not seem to be the tube radius!
  const double TUBE_RADIUS = (_major_radius - _minor_radius)*0.5;
  const double AVG_CENTER = (_minor_radius+_major_radius)*0.5;

  const double DSWEEP = (END_SWEEP - START_SWEEP)/(double)TORUS_SWEEPS;
  const double DPHI = (M_PI * 2.0) / (double)CIRCLE_SWEEPS;

  for(unsigned j=0; j<CIRCLE_SWEEPS; j++)
  {
    double phi = DPHI*(double)j;
    double cosPhi = std::cos(phi);
    double sinPhi = std::sin(phi);
    double next_cosPhi = std::cos(phi+DPHI);
    double next_sinPhi = std::sin(phi+DPHI);

    double z = TUBE_RADIUS*sinPhi;
    double yPrime = AVG_CENTER + TUBE_RADIUS*cosPhi;

    double next_z = TUBE_RADIUS*next_sinPhi;
    double next_yPrime = AVG_CENTER + TUBE_RADIUS*next_cosPhi;

    double old_x = yPrime*std::cos(-DSWEEP);
    double old_y = yPrime*std::sin(-DSWEEP);
    double old_z = z;

    for(unsigned i=0; i<TORUS_SWEEPS; ++i)
    {
      double sweep = START_SWEEP + DSWEEP*i;
      double cosSweep = std::cos(sweep);
      double sinSweep = std::sin(sweep);

      double x = yPrime*cosSweep;
      double y = yPrime*sinSweep;

      double next_x = next_yPrime*cosSweep;
      double next_y = next_yPrime*sinSweep;

      vertices.push_back( Point3d(next_x,next_y,next_z, P) );
      vertices.push_back( Point3d(x,y,z,P) );

      old_x = x; old_y = y; old_z = z;
    } // end torus loop

    // the last point
    double last_sweep = START_SWEEP + END_SWEEP;
    double cosLastSweep = std::cos(last_sweep);
    double sinLastSweep = std::sin(last_sweep);

    double x = yPrime*cosLastSweep;
    double y = yPrime*sinLastSweep;

    double next_x = next_yPrime*cosLastSweep;
    double next_y = next_yPrime*sinLastSweep;

    vertices.push_back( Point3d(next_x,next_y,next_z,P) );
    vertices.push_back( Point3d(x,y,z,P) );
  }  // end cirle loop
}

/// Calculates the signed distance from a point
double TorusPrimitive::calc_signed_dist(const Point3d& p) const
{
  throw std::runtime_error("TorusPrimitive::calc_signed_dist(.) not implemented");
}

/// Gets the mesh corresponding to this torus
shared_ptr<const IndexedTriArray> TorusPrimitive::get_mesh(shared_ptr<const Pose3d> P)
{
  throw std::runtime_error("TorusPrimitive::get_mesh(.) not implemented");
}

/// Computes the signed distance from the torus to a primitive
double TorusPrimitive::calc_signed_dist(shared_ptr<const Primitive> p, Point3d& pthis, Point3d& pp) const
{
  // attempt torus/plane
  shared_ptr<const PlanePrimitive> planep = dynamic_pointer_cast<const PlanePrimitive>(p);
  if (planep)
    return calc_signed_dist(planep, pthis, pp);

  throw std::runtime_error("Unsupported geometric pair"); 
}

/// Get random variable
double TorusPrimitive::urand(double a, double b)
{
  return (double) rand() / RAND_MAX * (b - a) + a;
}

/// Computes the signed distance from the torus to a plane
double TorusPrimitive::calc_signed_dist(shared_ptr<const PlanePrimitive> p, Point3d& pthis, Point3d& pp) const
{
  const unsigned Y = 1, Z = 2;
  const double EPS = NEAR_ZERO * 100.0;

  FILE_LOG(LOG_COLDET) << ">> start calc_signed_dist_torus_plane(.)" << std::endl;
  const double R = _major_radius; 
  const double r = _minor_radius; 

  // get the pose for the plane primitive
  shared_ptr<const Pose3d> Pplane = pp.pose; 

  // get the pose for the torus
  shared_ptr<const Pose3d> Ptorus = pthis.pose; 

  // get the transformation from the torus's space to the plane's space
  Transform3d tPp = Pose3d::calc_relative_pose(Pplane, Ptorus);

  // get plane normal in torus frame
  Vector3d n_plane = tPp.transform_vector(Vector3d(0.0, 1.0, 0.0,Pplane));
  n_plane.normalize();

  // Torus axis is z-axis in torus frame
  Vector3d k(0,0,1,Ptorus);

  // Set intitial value of distance to contact
  double d = std::numeric_limits<double>::infinity();

  // if Torus is aligned with plane:
  // Return distance torus origin to
  // closest point on plane less pipe r
  double n_dot_k = n_plane.dot(k);
  if (std::fabs(n_dot_k) > 1.0-EPS){
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

    // closest point is a random point on the circular bottom of the torus
    double t = urand(-M_PI_2,M_PI_2);
    pthis = Point3d(R*std::cos(t),R*std::sin(t),-r,Ptorus);

    // setup the point on the plane
    pp = Pose3d::transform_point(Pplane, pthis);
    pp[Y] = 0.0;

    FILE_LOG(LOG_COLDET) << " -- Torus is parallel to plane"<< std::endl;
    FILE_LOG(LOG_COLDET) << "distance: "<<  d << std::endl;
    FILE_LOG(LOG_COLDET) << "<< end calc_signed_dist_torus_plane(.)" << std::endl;
    return d;
  }

  //((n_plane x axis_torus) x axis_torus)
  Vector3d d_ring = Vector3d::cross(Vector3d::cross(n_plane,k), k);
  d_ring.normalize();

  // if Torus is _|_ with plane:
  // Return distance torus to plane less pipe r and ring R
  if(std::fabs(n_dot_k) < EPS){
    // d = depth
    // p0 = plane origin, p = plane normal
    // l0 = line origin, l = line direction

    // plane origin: plane origin in torus frame
    // line origin: torus origin in torus frame
    Point3d p0_plane(0.0, 0.0, 0.0,Pplane), l0(0,0,0,Ptorus);
    Point3d p0 = tPp.transform_point(p0_plane);

    // plane normal: plane normal in torus frame
    // line direction: on xy-plane of torus
    //   parallel to plane normal in torus frame
    Vector3d n = n_plane,l = d_ring;
    d = (p0 - l0).dot(n)/(l.dot(n)) - (r+R);
    pthis = d*l + l0;

    // setup the point on the plane
    pp = Pose3d::transform_point(Pplane, pthis);
    pp[Y] = 0.0;

    FILE_LOG(LOG_COLDET) << " -- Torus is perpendicular to plane"<< std::endl;
    FILE_LOG(LOG_COLDET) << "distance: "<<  d << std::endl;
    FILE_LOG(LOG_COLDET) << "<< end calc_signed_dist_torus_plane(.)" << std::endl;
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
  pthis = R * d_ring + r * d_pipe;

  // TODO: find the point in the torus's space such that
  //       tPp.transform_point(.) results in the value of y closest to
  //       negative infinity

  // setup the point on the plane
  pp = Pose3d::transform_point(Pplane, pthis);
  pp[Y] = 0.0;

  FILE_LOG(LOG_COLDET) << "distance: "<<  d << std::endl;
  FILE_LOG(LOG_COLDET) << "<< end calc_signed_dist_torus_plane(.)" << std::endl;

  return d;
}

/// Implements Base::load_from_xml() for serialization
void TorusPrimitive::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  // verify that the node type is TorusPrimitive
  assert(strcasecmp(node->name.c_str(), "Torus") == 0);

  // load the parent data
  Primitive::load_from_xml(node, id_map);

  // get the radius attributes, if specified
  XMLAttrib* major_rad_attr = node->get_attrib("major-radius");
  XMLAttrib* minor_rad_attr = node->get_attrib("minor-radius");

  // set lengths to zero initially
  double major_radius = 0.0, minor_radius = 0.0; 

  // get the lengths
  if (major_rad_attr) major_radius = major_rad_attr->get_real_value();
  if (minor_rad_attr) minor_radius = minor_rad_attr->get_real_value();

  // set the radii 
  set_radii(major_radius, minor_radius);
}

/// Sets the
void TorusPrimitive::set_radii(double major_radius, double minor_radius)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // set the radii
  _major_radius = major_radius;
  _minor_radius = minor_radius;

  // update each OBB 
  for (map<CollisionGeometryPtr, OBBPtr>::iterator i = _obbs.begin(); i != _obbs.end(); i++)
  {
    // setup OBB half-lengths
    i->second->l[X] = _major_radius + _minor_radius;
    i->second->l[Y] = _major_radius + _minor_radius;
    i->second->l[Z] = _minor_radius;
  }

  // recompute mass properties
  calc_mass_properties();
}

/// Implements Base::save_to_xml() for serialization
void TorusPrimitive::save_to_xml(XMLTreePtr node, std::list<shared_ptr<const Base> >& shared_objects) const
{
  // save the parent data
  Primitive::save_to_xml(node, shared_objects);

  // (re)set the node name
  node->name = "Torus";

  // save the lengths
  node->attribs.insert(XMLAttrib("major-radius", _major_radius));
  node->attribs.insert(XMLAttrib("minor-radius", _minor_radius));
}

/// Calculates mass properties for this primitive
void TorusPrimitive::calc_mass_properties()
{
  // get the current transform
  shared_ptr<const Pose3d> T = get_pose();

  // compute some constants
  const double a = _minor_radius;
  const double c = _major_radius; 
  const double M = _J.m;

  // compute the mass if necessary
  if (_density)
  {
    const double volume = M_PI * M_PI * 2.0 * a * a * c; 
    _J.m = *_density * volume;
  } 

  // compute the inertia matrix
  _J.J = Matrix3d(M*(5.0/8*a*a + 0.5*c*c), 0, 0, 0, M*(5.0/8*a*a + 0.5*c*c), 0, 0, 0, M*(0.75*a*a + c*c));
}

/*
class TorusPlanePlugin : public CollisionDetection
{
  private:
    ControlledBodyPtr walker;
    RigidBodyPtr ground_body, left_foot_body, right_foot_body;
    CollisionGeometryPtr ground_cg, left_foot_cg, right_foot_cg;

  public:

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

*/
