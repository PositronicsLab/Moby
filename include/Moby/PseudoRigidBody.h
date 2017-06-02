#ifndef PSEUDO_RIGID_BODY_H
#define PSEUDO_RIGID_BODY_H

#include <experimental/optional>

#include <osg/Array>
#include <osg/Shape>
#include <osg/ShapeDrawable>
#include <osg/Geode>
#include <osg/Geometry>

#include <Ravelin/RigidBodyd.h>
#include <Moby/ControlledBody.h>
#include <Moby/RigidBody.h>

namespace Moby {

class PseudoRigidBody : public ControlledBody, public Ravelin::DynamicBodyd
{
  friend class GravityForce;

  public:
    PseudoRigidBody();
    virtual void update_visualization() override;
    void load_from_xml(boost::shared_ptr<const XMLTree> tree, std::map<std::string, BasePtr>& id_map) override;
    void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base>>& shared_objects) const override;

    /// Body is always enabled.
    bool is_enabled() const override { return true; }

    /// Gets the total mass of the compliant layer. The default value is one
    /// tenth of the mass of the rigid core.
    double get_compliant_layer_mass() const
    {
      if (_cm_layer_mass)
        return _cm_layer_mass.value();
      else
        return _core->get_inertia().m;
    }

    /// Gets the total number of point masses used in the compliant layer.
    unsigned num_masses() const { return _point_mass_states.size(); }

    /// Gets the number of constraint equations.
    unsigned num_constraint_eqns() const { return _constraints.size()*2; } 

    Ravelin::MatrixNd& calc_jacobian(boost::shared_ptr<const Ravelin::Pose3d> source_pose, boost::shared_ptr<const Ravelin::Pose3d> target_pose, boost::shared_ptr<Ravelin::DynamicBodyd> body, Ravelin::MatrixNd& J) override;
    Ravelin::MatrixNd& calc_jacobian_dot(boost::shared_ptr<const Ravelin::Pose3d> source_pose, boost::shared_ptr<const Ravelin::Pose3d> target_pose, boost::shared_ptr<Ravelin::DynamicBodyd> body, Ravelin::MatrixNd& Jdot) override;
    void reset_accumulators() override;
    void translate(const Ravelin::Origin3d& o) override;
    void rotate(const Ravelin::Quatd& q) override;
    double calc_kinetic_energy(boost::shared_ptr<const Ravelin::Pose3d> P) override;
    boost::shared_ptr<const Ravelin::Pose3d> get_gc_pose() const override;
    void set_computation_frame_type(Ravelin::ReferenceFrameType rftype) override;
    void set_generalized_forces(const Ravelin::SharedVectorNd& gf) override;
    void add_generalized_force(const Ravelin::SharedVectorNd& gf) override;
    void apply_generalized_impulse(const Ravelin::SharedVectorNd& gj) override;
    Ravelin::SharedVectorNd& get_generalized_acceleration(Ravelin::SharedVectorNd& ga) override;
    void set_generalized_acceleration(const Ravelin::SharedVectorNd& ga) override;
    void set_generalized_coordinates_euler(const Ravelin::SharedVectorNd& gc) override;
    Ravelin::SharedVectorNd& get_generalized_velocity(Ravelin::DynamicBodyd::GeneralizedCoordinateType gctype, Ravelin::SharedVectorNd& gv) override;
    void set_generalized_velocity(Ravelin::DynamicBodyd::GeneralizedCoordinateType gctype, const Ravelin::SharedVectorNd& gv) override;
    Ravelin::SharedMatrixNd& get_generalized_inertia(Ravelin::SharedMatrixNd& M) override;
    Ravelin::SharedMatrixNd& solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X) override;
    Ravelin::SharedMatrixNd& transpose_solve_generalized_inertia(const Ravelin::SharedMatrixNd& B, Ravelin::SharedMatrixNd& X) override;
    Ravelin::SharedVectorNd& solve_generalized_inertia(const Ravelin::SharedVectorNd& b, Ravelin::SharedVectorNd& x) override;
    Ravelin::SharedVectorNd& get_generalized_forces(Ravelin::SharedVectorNd& gf) override;
    Ravelin::SharedVectorNd& convert_to_generalized_force(boost::shared_ptr<Ravelin::SingleBodyd> body, const Ravelin::SForced& w, Ravelin::SharedVectorNd& gf) override;
    Ravelin::SharedVectorNd& get_generalized_coordinates_euler(Ravelin::SharedVectorNd& gc) override;
    void calc_fwd_dyn() override;
    unsigned num_generalized_coordinates(Ravelin::DynamicBodyd::GeneralizedCoordinateType gcypte) const override;

  protected:
    void prepare_to_calc_ode(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data) override;
    void prepare_to_calc_ode_sustained_constraints(Ravelin::SharedConstVectorNd& x, double t, double dt, void* data) override;
    void ode(double t, double dt, void* data, Ravelin::SharedVectorNd& dx) override;
 
  private:
    Ravelin::LinAlgd _LA;
    void update_poses();
    void update_generalized_forces_with_spring_damper_forces(Ravelin::VectorNd& gf);
    void convert_spring_force_to_generalized_force(unsigned i, double f, Ravelin::VectorNd& gf) const;
    double determine_spring_extension(unsigned i) const;
    double determine_spring_extension_dot(unsigned i) const;
    void calc_jacobian(Ravelin::MatrixNd& J);
    void calc_jacobian_dot_v(Ravelin::VectorNd& Jdot_v);
    void evaluate_constraints(Ravelin::VectorNd& C); 
    void evaluate_constraints_dot(Ravelin::VectorNd& Cdot); 
    void evaluate_constraints_ddot(Ravelin::VectorNd& Cddot); 
    Ravelin::VectorNd& transform_piston_force_to_generalized_force(unsigned piston_index, double f, Ravelin::VectorNd& gf) const;

    // All information necessary for the point mass dynamics.
    struct PointMass
    {
      PointMass(const Point3d& x)
      {
        assert(x.pose == GLOBAL);
        this->x = x;
        xdot = Point3d(0,0,0,GLOBAL);
        xddot = Point3d(0,0,0,GLOBAL);
        f = Point3d(0,0,0,GLOBAL);
      }

      Point3d x;                 // The position of the point mass.
      Ravelin::Vector3d xdot;    // The velocity of the point mass.
      Ravelin::Vector3d xddot;   // The acceleration of the point mass.
      Ravelin::Vector3d f;       // The force acting on the point mass.
    };

    // Structure for piston constraints.
    struct PistonConstraint
    {
      PistonConstraint(const Ravelin::Vector3d& u_, const Ravelin::Vector3d& v1_, const Ravelin::Vector3d& v2_, const Point3d& p_, unsigned point_mass_index_) :
        u(u_), v1(v1_), v2(v2_), p(p_), point_mass_index(point_mass_index_) {}

      // The piston axis emanating from the rigid body.
      Ravelin::Vector3d u;

      // Vectors normal to the piston axis emanating from the rigid body.
      Ravelin::Vector3d v1, v2;

      // The attachment point of the piston on the rigid body.
      Point3d p;

      // The point mass that this piston is attached to.
      unsigned point_mass_index;
    };

    // Structure for spring and damper data.
    struct SpringDamperData
    {
      enum class ConstraintType { ePiston, eTwoPointMasses };

      SpringDamperData(unsigned i, unsigned j, double len) : 
        pm1_index(i), pm2_index(j), rest_len(len)
      {
        constraint_type = ConstraintType::eTwoPointMasses;
        piston_index = std::numeric_limits<unsigned>::max();
      }

      explicit SpringDamperData(unsigned i, double len) :
        piston_index(i), rest_len(len)
      {
        constraint_type = ConstraintType::ePiston;
        pm1_index = pm2_index = std::numeric_limits<unsigned>::max();
      }

      // The resting length of the spring
      double rest_len{-1};

      // The type of the constraint.
      ConstraintType constraint_type;
      
      // The index of the two point masses, for types of eTwoPointMasses.
      unsigned pm1_index, pm2_index;

      // The index of the piston constraint, if relevant
      unsigned piston_index; 
    };

    void update_collision_geometry_primitives();
    RigidBodyPtr _core; 

    // Uniform stiffness and damping for all springs
    double _spring_k, _spring_c;

    // The vector of vertices.
    osg::Vec3Array* _vert_array{nullptr};

    // The vector of faces for the pseudo-rigid body visualization.
    std::vector<osg::DrawElementsUInt*> _osg_faces;

    // The states of the point masses.
    std::vector<PointMass> _point_mass_states;

    // Constraint information for the pistons.
    std::vector<PistonConstraint> _constraints;

    // The set of springs/dampers.
    std::vector<SpringDamperData> _springs;

    // The compliant layer mass.
    std::experimental::optional<double> _cm_layer_mass;
};

}

#endif

