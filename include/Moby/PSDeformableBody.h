/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _PS_DEFORMABLE_BODY_H
#define _PS_DEFORMABLE_BODY_H

#include <Moby/DeformableBody.h>

namespace Moby {

/// A spring in the mass-spring system
struct Spring
{
  Spring()
  {
    kp = kv = rest_len = 0;
    node1 = node2 = std::numeric_limits<unsigned>::max();
  }

  /// The spring stiffness
  double kp;

  /// The spring dampening
  double kv;

  /// The first node connected by the spring
  unsigned node1;

  /// The second node connected by the spring
  unsigned node2;

  /// The rest length of the spring
  double rest_len;
};

/// Class for deformable bodies simulated using systems of particles 
class PSDeformableBody : public DeformableBody 
{
  public:
    PSDeformableBody();
    virtual void integrate(double t, double h, boost::shared_ptr<Integrator> integrator);
    virtual void apply_impulse(const Ravelin::SMomentumd& j);
    virtual void calc_fwd_dyn(double dt);
    virtual double calc_mass() const { return _J.m; }
    virtual void set_mesh(boost::shared_ptr<const IndexedTetraArray> tetra_mesh, boost::shared_ptr<Primitive> tri_mesh);
    const Spring& get_spring(unsigned i) { return _springs[i]; }
    void set_spring(unsigned i, Spring &s) { _springs[i] = s; }
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    virtual double calc_potential_energy() const;
    virtual boost::shared_ptr<Ravelin::Pose3d> get_computation_frame_pose() const;
    virtual void set_computation_frame_type(ReferenceFrameType rftype);
    virtual double get_mass() const;

    /// The default spring stiffness constant
    double default_KP;

    /// The default spring dampening constant
    double default_KV;

  private:
    void determine_Dc_v(const std::vector<Event*>& contact_events, Ravelin::VectorNd& Dc_v) const;
    void determine_Jc_v(const std::vector<Event*>& contact_events, Ravelin::VectorNd& Jc_v) const;

    /// The set of springs
    std::vector<Spring> _springs;
}; // end class

} // end namespace

#endif

