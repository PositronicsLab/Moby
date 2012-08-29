/*

Makes contact points and normals

*/

#include <cmath>
#include <fstream>
#include <Moby/CollisionDetection.h>
#include <Moby/Event.h>
#include <Moby/Constants.h>

using boost::dynamic_pointer_cast;
using boost::shared_ptr;
using std::vector;
using namespace Moby;

class ColdetPlugin : public CollisionDetection
{
  public:
    virtual void add_collision_geometry(CollisionGeometryPtr cg);
    virtual bool is_contact(Real dt, std::vector<Event>& contacts);
    virtual bool is_collision(Real epsilon = 0.0) { return false; }
    virtual Real calc_distances() { return 0.0; }

  private:
    CollisionGeometryPtr _platform, _base, _lwheel, _rwheel;
};

void ColdetPlugin::add_collision_geometry(CollisionGeometryPtr cg)
{
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(cg->get_single_body());
  if (rb->id.find("base") != std::string::npos)
    _base = cg;
  else if (rb->id.find("wheel-left") != std::string::npos)
    _lwheel = cg;
  else if (rb->id.find("wheel-right") != std::string::npos)
    _rwheel = cg;
  else if (rb->id.find("platform") != std::string::npos)
    _platform = cg;
  else
    std::cout << "Unrecognized geometry!" << std::endl;
}

Real rnd() { return (Real) rand() / RAND_MAX; }

bool ColdetPlugin::is_contact(Real dt, std::vector<Event>& contacts)
{
  const Real WHEEL_RADIUS = .1575;

  // setup the "ideal" contact points according to the base
  vector<Vector3> p;
  p.push_back(Vector3(-.16,-.114,-.055));
  p.push_back(Vector3(.16,-.114,-.055));

  // determine the normal vector corresponding to each point
  vector<Vector3> n;
  for (unsigned i=0; i< p.size(); i++)
  {
    Vector3 cleft(-.16,.0715,-.055);
    Vector3 cright(.16,.0715,-.055);
/*
    if (i < p.size()/2)
      n.push_back(Vector3::normalize(cleft - p[i]));
    else
      n.push_back(Vector3::normalize(cright - p[i]));
*/
    n.push_back(Vector3(0,1,0));
  }

  // determine the transformed contact points and transformed normalS
  for (unsigned i=0; i< p.size(); i++)
  {
    p[i] = _base->get_transform().mult_point(p[i]);
    n[i] = _base->get_transform().mult_vector(n[i]);
//    std::cout << "generated contact point: " << p[i] << std::endl;
//    std::cout << "          and normal: " << n[i] << std::endl;
  }

  // setup the contacts
  for (unsigned i=0; i< p.size(); i++)
  {
    Event e;
    e.t = (Real) 0.0;
    e.event_type = Event::eContact;
    e.contact_geom1 = (i < p.size()/2) ? _lwheel : _rwheel; 
    e.contact_geom2 = _platform; 
    e.contact_normal = n[i];
    e.contact_point = p[i];
    contacts.push_back(e);
  }

  // we always want to return contact
  return true;
}

extern "C"
{
  boost::shared_ptr<CollisionDetection> factory() { return boost::shared_ptr<CollisionDetection>(new ColdetPlugin); }
} // end extern "C"

