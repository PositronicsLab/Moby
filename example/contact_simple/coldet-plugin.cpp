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
using namespace Moby;

class ColdetPlugin : public CollisionDetection
{
  public:
    virtual void add_collision_geometry(CollisionGeometryPtr cg);
    virtual bool is_contact(Real dt, std::vector<Event>& contacts);
    virtual bool is_collision(Real epsilon = 0.0) { return false; }
    virtual Real calc_distances() { return 0.0; }

  private:
    std::vector<CollisionGeometryPtr> _boxes;
    CollisionGeometryPtr _ground;
};

void ColdetPlugin::add_collision_geometry(CollisionGeometryPtr cg)
{
  RigidBodyPtr rb = dynamic_pointer_cast<RigidBody>(cg->get_single_body());
  if (rb->id.find("box") != std::string::npos)
  {
    // add the box
    _boxes.push_back(cg);

    // sort all collision geometries
    std::map<std::string, CollisionGeometryPtr> geoms;
    for (unsigned i=0; i< _boxes.size(); i++)
      geoms[_boxes[i]->get_single_body()->id] = _boxes[i];

    // clear the boxes
    _boxes.clear();
    for (std::map<std::string, CollisionGeometryPtr>::const_iterator i = geoms.begin(); i != geoms.end(); i++)
      _boxes.push_back(i->second);
  }
  else
  {
    assert(rb->id == "ground");
    _ground = cg;
  }
}

bool ColdetPlugin::is_contact(Real dt, std::vector<Event>& contacts)
{
  const unsigned Y = 1;

  // create contacts between boxes
  for (unsigned i=0; i< _boxes.size()-1; i++)
  {
    // get the y position for the contacts
    RigidBodyPtr rb1 = dynamic_pointer_cast<RigidBody>(_boxes[i]->get_single_body());
    RigidBodyPtr rb2 = dynamic_pointer_cast<RigidBody>(_boxes[i+1]->get_single_body());
    Real y = (rb1->get_position()[Y] + rb2->get_position()[Y])*0.5;

    Event e;
    e.event_type = Event::eContact;
    e.t = (Real) 0.0;
    e.contact_geom1 = _boxes[i];
    e.contact_geom2 = _boxes[i+1];
    e.contact_normal = Vector3(0,-1,0);
    e.contact_point = Vector3(-0.2, y, -0.2);
    contacts.push_back(e);
    e.contact_point = Vector3(-0.2, y, 0.2);
    contacts.push_back(e);
    e.contact_point = Vector3(0.2, y, -0.2);
    contacts.push_back(e);
    e.contact_point = Vector3(0.2, y, 0.2);
    contacts.push_back(e);
  }

  // create contacts between box and the ground
  Event e;
  e.t = (Real) 0.0;
  e.event_type = Event::eContact;
  e.contact_geom1 = _ground;
  e.contact_geom2 = _boxes.front();
  e.contact_normal = Vector3(0,-1,0);
  e.contact_point = Vector3(-0.5, 0, -0.5);
  contacts.push_back(e);
  e.contact_point = Vector3(-0.5, 0, 0.5);
  contacts.push_back(e);
  e.contact_point = Vector3(0.5, 0, -0.5);
  contacts.push_back(e);
  e.contact_point = Vector3(0.5, 0, 0.5);
  contacts.push_back(e);

  return true;
}

extern "C"
{
  boost::shared_ptr<CollisionDetection> factory() { return boost::shared_ptr<CollisionDetection>(new ColdetPlugin); }
} // end extern "C"

