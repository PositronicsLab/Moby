#include <Moby/DynamicBody.h>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;
using boost::dynamic_pointer_cast;
using std::list;

/*
/// Integrates a dynamic body
void DynamicBody::integrate(double t, double h, shared_ptr<Integrator> integrator)
{
  FILE_LOG(LOG_DYNAMICS) << "DynamicBody::integrate() - integrating from " << t << " by " << h << std::endl;

  if (_kinematic_update)
  {
    FILE_LOG(LOG_DYNAMICS) << " -- body set to kinematic update --" << std::endl;
    if (controller)
      (*controller)(dynamic_pointer_cast<DynamicBody>(shared_from_this()), t, controller_arg);
    return;
  }

  shared_ptr<DynamicBody> shared_this = dynamic_pointer_cast<DynamicBody>(shared_from_this());

  get_generalized_coordinates(eEuler, gc);
  get_generalized_velocity(eSpatial, gv); // gv depends on gc
  gcgv.resize(gc.size()+gv.size());
  gcgv.set_sub_vec(0, gc);
  gcgv.set_sub_vec(gc.size(), gv);
  integrator->integrate(gcgv, &ode_both, t, h, (void*) &shared_this);
  gcgv.get_sub_vec(0, gc.size(), gc);
  gcgv.get_sub_vec(gc.size(), gcgv.size(), gv);
  // NOTE: velocity must be set first (it's computed w.r.t. old frame)
  set_generalized_coordinates(eEuler, gc);
  set_generalized_velocity(eSpatial, gv);
}

/// Returns the ODE's for position and velocity (concatenated into x)
VectorNd& DynamicBody::ode_both(const VectorNd& x, double t, double dt, void* data, VectorNd& dx)
{
  // get the dynamic body
  shared_ptr<DynamicBody>& db = *((shared_ptr<DynamicBody>*) data);
  const unsigned NGC_EUL = db->num_generalized_coordinates(eEuler);

  // get the necessary vectors
  VectorNd& xp = db->xp;
  VectorNd& xv = db->xv;
  VectorNd& xa = db->xa;

  // return the derivatives at state x
  xv.resize(NGC_EUL);
  x.get_sub_vec(0, NGC_EUL, xp);
  x.get_sub_vec(NGC_EUL, x.size(), xv);
  db->set_generalized_coordinates(DynamicBody::eEuler, xp);
  db->set_generalized_velocity(DynamicBody::eSpatial, xv);

  // we need the generalized velocity as Rodrigues coordinates
  db->get_generalized_velocity(DynamicBody::eEuler, xv);

  // check whether we could rotate too much
  #ifdef NDEBUG
  if (dt*db->get_aspeed() > M_PI)
    std::cerr << "DynamicBody::ode_both() warning- angular speed*dt " << (dt*db->get_aspeed()) << " sufficiently high to" << std::endl << "potentially miss events!" << std::endl;
  #endif

  // clear the force accumulators on the body
  db->reset_accumulators();

  // add all recurrent forces on the body
  const list<RecurrentForcePtr>& rfs = db->get_recurrent_forces();
  BOOST_FOREACH(RecurrentForcePtr rf, rfs)
    rf->add_force(db);

  // call the body's controller
  if (db->controller)
    (*db->controller)(db, t, db->controller_arg);

  // calculate forward dynamics at state x
  db->calc_fwd_dyn();
  db->get_generalized_acceleration(xa);

  dx.resize(x.size());
  dx.set_sub_vec(0, xv);
  dx.set_sub_vec(xv.size(), xa);
  return dx;
}
*/

/// Loads the body's state via XML
void DynamicBody::load_from_xml(shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map)
{
  std::map<std::string, BasePtr>::const_iterator id_iter;

  // load Base data
  Visualizable::load_from_xml(node, id_map);

  // get all recurrent forces used in the simulator -- note: this must be done
  // *after* all bodies have been loaded
  list<shared_ptr<const XMLTree> > child_nodes = node->find_child_nodes("RecurrentForce");
  if (!child_nodes.empty())
  {
    // safe to clear the vector of recurrent forces
    _rfs.clear();

    // process all child nodes
    for (list<shared_ptr<const XMLTree> >::const_iterator i = child_nodes.begin(); i != child_nodes.end(); i++)
    {
      // verify that the dynamic-body-id attribute exists
      XMLAttrib* id_attr = (*i)->get_attrib("recurrent-force-id");

      // make sure that the ID exists
      if (!id_attr)
      {
        std::cerr << "DynamicBody::load_from_xml() - no recurrent-force-id ";
        std::cerr << "attribute in tag: " << node << std::endl;
        continue;
      }

      // look for the recurrent force with that ID
      const std::string& id = id_attr->get_string_value(); 
      if ((id_iter = id_map.find(id)) == id_map.end())
      {
        std::cerr << "DynamicBody::load_from_xml() - could not find" << std::endl;
        std::cerr << "  recurrent force w/ID: " << id << " from offending node: " << std::endl << *node;
      }
      else
      {
        RecurrentForcePtr rf = dynamic_pointer_cast<RecurrentForce>(id_iter->second);
          _rfs.push_back(rf);
      }
    }
  }
}

/// Implements Base::save_to_xml()
void DynamicBody::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call the parent save_to_xml() method
  Visualizable::save_to_xml(node, shared_objects);

  // rename the node
  node->name = "DynamicBody";

  // save the IDs of all recurrent forces
  BOOST_FOREACH(RecurrentForcePtr rf, _rfs)
  {
    assert(rf);
    XMLTreePtr child_node(new XMLTree("RecurrentForce"));
    node->add_child(child_node);
    child_node->attribs.insert(XMLAttrib("recurrent-force-id", rf->id));
    shared_objects.push_back(rf);
  }
}

