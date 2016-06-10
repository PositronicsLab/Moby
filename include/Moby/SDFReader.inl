/// Special routine for creating a gearbox
template <class OutputIterator>
static OutputIterator read_gearbox(boost::shared_ptr<const XMLTree> node, const std::map<std::string, RigidBodyPtr>& link_map, RigidBodyPtr& base_link, OutputIterator output_begin)
{
  std::string name, type;
  boost::shared_ptr<Ravelin::Pose3d> P(new Ravelin::Pose3d);
  RigidBodyPtr refbody, parent, child;
  boost::shared_ptr<RevoluteJoint> rj1, rj2;
  boost::shared_ptr<Gears> gj;

  // get the id of the joint
  XMLAttrib* name_attr = node->get_attrib("name");
  if (name_attr)
    name = name_attr->get_string_value();

  // create three joints
  rj1 = boost::shared_ptr<RevoluteJoint>(new RevoluteJoint);
  rj2 = boost::shared_ptr<RevoluteJoint>(new RevoluteJoint);
  gj = boost::shared_ptr<Gears>(new Gears);

  // setup the generic components of the joints
  rj1->id = name;
  rj1->joint_id = rj1->id;
  rj2->id = name + "'";
  rj2->joint_id = rj2->id;
  gj->id = name + "_gears";
  gj->joint_id =  gj->id;

  // read in the gearbox ratio
  boost::shared_ptr<const XMLTree> gearbox_ratio_tag = find_one_tag("gearbox_ratio", node);
  assert(gearbox_ratio_tag);
  gj->set_ratio(std::atof(gearbox_ratio_tag->content.c_str()));

  // read in the name of the reference body
  boost::shared_ptr<const XMLTree> refbody_link_tag = find_one_tag("gearbox_reference_body", node);
  assert(refbody_link_tag);
  std::string refbody_link = refbody_link_tag->content;
  if (strcasecmp(refbody_link.c_str(), "world") == 0)
  {
    // joint is attached to the world; create the body, if not already done
    if (!base_link)
    {
      base_link = RigidBodyPtr(new RigidBody);
      base_link->set_enabled(false);
    }
    refbody = base_link;
  }
  else if (link_map.find(refbody_link) == link_map.end())
  {
    std::string except_string = "SDFReader::read_gearbox(.)- refbody link '" + refbody_link + "' not found";
    throw std::runtime_error(except_string.c_str());
  }
  else
    refbody = link_map.find(refbody_link)->second;
  
  // read in the name of the parent link
  boost::shared_ptr<const XMLTree> parent_tag = find_one_tag("parent", node);
  assert(parent_tag);
  std::string parent_link = parent_tag->content;
  if (strcasecmp(parent_link.c_str(), "world") == 0)
  {
    // joint is attached to the world; create the body, if not already done
    if (!base_link)
    {
      base_link = RigidBodyPtr(new RigidBody);
      base_link->set_enabled(false);
    }
    parent = base_link;
  }
  else if (link_map.find(parent_link) == link_map.end())
  {
    std::string except_string = "SDFReader::read_gearbox(.)- parent link '" + parent_link + "' not found";
    throw std::runtime_error(except_string.c_str());
  }
  else
    parent = link_map.find(parent_link)->second;

  // read in the name of the child link
  boost::shared_ptr<const XMLTree> child_tag = find_one_tag("child", node);
  assert(child_tag);
  std::string child_link = child_tag->content;
  if (link_map.find(child_link) == link_map.end())
  {
    std::string except_string = "SDFReader::read_gearbox(.)- child link '" + child_link + "' not found";
    throw std::runtime_error(except_string.c_str());
  }
  child = link_map.find(child_link)->second;

  // get the pose of the joints
  boost::shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);
  Ravelin::Pose3d Px;
  if (pose_node)
  {
    Px = read_pose(node);
    Px.rpose = child->get_pose();
    Px.update_relative_pose(GLOBAL);
    Ravelin::Vector3d loc(Px.x[0], Px.x[1], Px.x[2], GLOBAL);
    rj1->set_location(loc, parent, child);
    rj2->set_location(loc, parent, refbody);
    gj->set_location(loc, refbody, child);
  }
  else
  {
    Px.rpose = child->get_pose();
    Px.update_relative_pose(GLOBAL);
    Ravelin::Vector3d loc(Px.x[0], Px.x[1], Px.x[2], GLOBAL);
    rj1->set_location(loc, parent, child);
    rj2->set_location(loc, parent, refbody);
    gj->set_location(loc, refbody, child);
  }
  *(P.get()) = Px;

  // read the axis tag (contains limits, joint damping/friction)
  boost::shared_ptr<const XMLTree> axis_node = find_one_tag("axis", node);
  if (axis_node)
  {
    // read whether to use the parent model frame
    boost::shared_ptr<const XMLTree> parent_model_frame_node = find_one_tag("use_parent_model_frame", axis_node);

    // read the axis, if any
    boost::shared_ptr<const XMLTree> xyz_node = find_one_tag("xyz", axis_node);
    if (xyz_node)
    {
      // get the axis
      Ravelin::Vector3d axis = read_Vector3(xyz_node);

      // set the axis pose
      if (parent_model_frame_node && read_bool(parent_model_frame_node))
        axis.pose = parent->get_pose();
      else
        axis.pose = P;

      // set the axis
      rj1->set_axis(axis);
    }

    // read viscous and Coulomb friction
    boost::shared_ptr<const XMLTree> dynamics_node = find_one_tag("dynamics", axis_node);
    if (dynamics_node)
    {
      // attempt to read 'Coulomb' friction
      boost::shared_ptr<const XMLTree> friction_node = find_one_tag("friction", dynamics_node);
      if (friction_node)
        rj1->mu_fc = read_double(friction_node);

      // attempt to read viscous friction
      boost::shared_ptr<const XMLTree> damping_node = find_one_tag("damping", dynamics_node);
      if (damping_node)
        rj1->mu_fv = read_double(damping_node);
    }

    // read joint limits
    boost::shared_ptr<const XMLTree> limit_node = find_one_tag("limit", axis_node);
    if (limit_node)
    {
      // attempt to read the lower limit
      boost::shared_ptr<const XMLTree> llimit_node = find_one_tag("lower", limit_node);
      if (llimit_node)
        rj1->lolimit[0] = read_double(llimit_node);

      // attempt to read the upper limit
      boost::shared_ptr<const XMLTree> ulimit_node = find_one_tag("upper", limit_node);
      if (ulimit_node)
        rj1->hilimit[0] = read_double(ulimit_node);
    }
  }

  // read the axis tag (contains limits, joint damping/friction)
  boost::shared_ptr<const XMLTree> axis2_node = find_one_tag("axis2", node);
  if (axis2_node)
  {
    // read whether to use the parent model frame
    boost::shared_ptr<const XMLTree> parent_model_frame_node = find_one_tag("use_parent_model_frame", axis2_node);

    // read the axis, if any
    boost::shared_ptr<const XMLTree> xyz_node = find_one_tag("xyz", axis2_node);
    if (xyz_node)
    {
      // get the axis
      Ravelin::Vector3d axis = read_Vector3(xyz_node);

      // set the axis pose
      if (parent_model_frame_node && read_bool(parent_model_frame_node))
        axis.pose = parent->get_pose();
      else
        axis.pose = P;

      // set the axis
      rj2->set_axis(axis);
    }

    // read viscous and Coulomb friction
    boost::shared_ptr<const XMLTree> dynamics_node = find_one_tag("dynamics", axis2_node);
    if (dynamics_node)
    {
      // attempt to read 'Coulomb' friction
      boost::shared_ptr<const XMLTree> friction_node = find_one_tag("friction", dynamics_node);
      if (friction_node)
        rj2->mu_fc = read_double(friction_node);

      // attempt to read viscous friction
      boost::shared_ptr<const XMLTree> damping_node = find_one_tag("damping", dynamics_node);
      if (damping_node)
        rj2->mu_fv = read_double(damping_node);
    }

    // read joint limits
    boost::shared_ptr<const XMLTree> limit_node = find_one_tag("limits", axis2_node);
    if (limit_node)
    {
      // attempt to read the lower limit
      boost::shared_ptr<const XMLTree> llimit_node = find_one_tag("lower", limit_node);
      if (llimit_node)
        rj2->lolimit[1] = read_double(llimit_node);

      // attempt to read the upper limit
      boost::shared_ptr<const XMLTree> ulimit_node = find_one_tag("upper", limit_node);
      if (ulimit_node)
        rj2->hilimit[1] = read_double(ulimit_node);
    }
  }

  *output_begin++ = rj1;
  *output_begin++ = rj2;
  *output_begin++ = gj;
  return output_begin;
}

/// Reads a joint
template <class OutputIterator>
static OutputIterator read_joint(boost::shared_ptr<const XMLTree> node, const std::map<std::string, RigidBodyPtr>& link_map, RigidBodyPtr& base_link, OutputIterator output_begin)
{
  std::string name, type;
  JointPtr joint;
  boost::shared_ptr<Ravelin::Pose3d> P(new Ravelin::Pose3d);
  RigidBodyPtr parent, child;
  boost::shared_ptr<RevoluteJoint> rj;
  boost::shared_ptr<PrismaticJoint> pj;
  boost::shared_ptr<UniversalJoint> uj;
  boost::shared_ptr<SphericalJoint> sj;

  // get the id of the joint
  XMLAttrib* name_attr = node->get_attrib("name");
  if (name_attr)
    name = name_attr->get_string_value();

  // get the joint type
  XMLAttrib* type_attr = node->get_attrib("type");
  if (type_attr)
  {
    type = type_attr->get_string_value();
    std::transform(type.begin(), type.end(), type.begin(), ::tolower);
    if (type == "revolute")
    {
      rj = boost::shared_ptr<RevoluteJoint>(new RevoluteJoint);
      joint = rj;
    }
    else if (type == "prismatic")
    {
      pj = boost::shared_ptr<PrismaticJoint>(new PrismaticJoint);
      joint = pj;
    }
    else if (type == "universal")
    {
      uj = boost::shared_ptr<UniversalJoint>(new UniversalJoint);
      joint = uj;
    }
    else if (type == "ball")
    {
      sj = boost::shared_ptr<SphericalJoint>(new SphericalJoint);
      joint = sj;
    }
    else if (type == "gearbox")
      return read_gearbox(node, link_map, base_link, output_begin);
    else if (type == "revolute2")
      throw std::runtime_error("revolute2 joint unsupported");
    else if (type == "piston")
      throw std::runtime_error("piston joint unsupported");
    else
    {
      std::string err = type + " joint unsupported";
      throw std::runtime_error(err.c_str());
    }
  }

  // setup the generic components of the joint
  joint->id = name;
  joint->joint_id = name;

  // read in the name of the parent link
  boost::shared_ptr<const XMLTree> parent_tag = find_one_tag("parent", node);
  assert(parent_tag);
  std::string parent_link = parent_tag->content;
  if (strcasecmp(parent_link.c_str(), "world") == 0)
  {
    // joint is attached to the world; create the body, if not already done
    if (!base_link)
    {
      base_link = RigidBodyPtr(new RigidBody);
      base_link->set_enabled(false);
    }
    parent = base_link;
  }
  else if (link_map.find(parent_link) == link_map.end())
  {
    std::string except_string = "SDFReader::read_joint(.)- parent link '" + parent_link + "' not found";
    throw std::runtime_error(except_string.c_str());
  }
  else
    parent = link_map.find(parent_link)->second;

  // read in the name of the child link
  boost::shared_ptr<const XMLTree> child_tag = find_one_tag("child", node);
  assert(child_tag);
  std::string child_link = child_tag->content;
  if (link_map.find(child_link) == link_map.end())
  {
    std::string except_string = "SDFReader::read_joint(.)- child link '" + child_link + "' not found";
    throw std::runtime_error(except_string.c_str());
  }
  child = link_map.find(child_link)->second;

  // get the pose of the joint
  boost::shared_ptr<const XMLTree> pose_node = find_one_tag("pose", node);
  Ravelin::Pose3d Px;
  if (pose_node)
  {
    Px = read_pose(node);
    Px.rpose = child->get_pose();
    Px.update_relative_pose(GLOBAL);
    Ravelin::Vector3d loc(Px.x[0], Px.x[1], Px.x[2], GLOBAL);
    joint->set_location(loc, parent, child);
  }
  else
  {
    Px.rpose = child->get_pose();
    Px.update_relative_pose(GLOBAL);
    Ravelin::Vector3d loc(Px.x[0], Px.x[1], Px.x[2], GLOBAL);
    joint->set_location(loc, parent, child);
  }
  *(P.get()) = Px;

  // read the axis tag (contains limits, joint damping/friction)
  boost::shared_ptr<const XMLTree> axis_node = find_one_tag("axis", node);
  if (axis_node)
  {
    // read whether to use the parent model frame
    boost::shared_ptr<const XMLTree> parent_model_frame_node = find_one_tag("use_parent_model_frame", axis_node);

    // read the axis, if any
    boost::shared_ptr<const XMLTree> xyz_node = find_one_tag("xyz", axis_node);
    if (xyz_node)
    {
      // get the axis
      Ravelin::Vector3d axis = read_Vector3(xyz_node);

      // set the axis pose
      if (parent_model_frame_node && read_bool(parent_model_frame_node))
        axis.pose = parent->get_pose();
      else
        axis.pose = P;

      // set the axis
      if (rj)
        rj->set_axis(axis);
      else if (pj)
        pj->set_axis(axis);
      else if (uj)
        uj->set_axis(axis, UniversalJoint::eAxis1);
      else if (sj)
        sj->set_axis(axis, SphericalJoint::eAxis1);
    }

    // read viscous and Coulomb friction
    boost::shared_ptr<const XMLTree> dynamics_node = find_one_tag("dynamics", axis_node);
    if (dynamics_node)
    {
      // attempt to read 'Coulomb' friction
      boost::shared_ptr<const XMLTree> friction_node = find_one_tag("friction", dynamics_node);
      if (friction_node)
        joint->mu_fc = read_double(friction_node);

      // attempt to read viscous friction
      boost::shared_ptr<const XMLTree> damping_node = find_one_tag("damping", dynamics_node);
      if (damping_node)
        joint->mu_fv = read_double(damping_node);
    }

    // read joint limits
    boost::shared_ptr<const XMLTree> limit_node = find_one_tag("limit", axis_node);
    if (limit_node)
    {
      // attempt to read the lower limit
      boost::shared_ptr<const XMLTree> llimit_node = find_one_tag("lower", limit_node);
      if (llimit_node)
        joint->lolimit[0] = read_double(llimit_node);

      // attempt to read the upper limit
      boost::shared_ptr<const XMLTree> ulimit_node = find_one_tag("upper", limit_node);
      if (ulimit_node)
        joint->hilimit[0] = read_double(ulimit_node);
    }
  }

  // read the axis tag (contains limits, joint damping/friction)
  boost::shared_ptr<const XMLTree> axis2_node = find_one_tag("axis2", node);
  if (axis2_node)
  {
    // read whether to use the parent model frame
    boost::shared_ptr<const XMLTree> parent_model_frame_node = find_one_tag("use_parent_model_frame", axis2_node);

    // read the axis, if any
    boost::shared_ptr<const XMLTree> xyz_node = find_one_tag("xyz", axis2_node);
    if (xyz_node)
    {
      // get the axis
      Ravelin::Vector3d axis = read_Vector3(xyz_node);

      // set the axis pose
      if (parent_model_frame_node && read_bool(parent_model_frame_node))
        axis.pose = parent->get_pose();
      else
        axis.pose = P;

      // set the axis
      if (rj)
        rj->set_axis(axis);
      else if (pj)
        pj->set_axis(axis);
      else if (uj)
        uj->set_axis(axis, UniversalJoint::eAxis2);
      else if (sj)
        sj->set_axis(axis, SphericalJoint::eAxis2);
    }

    // read viscous and Coulomb friction
    boost::shared_ptr<const XMLTree> dynamics_node = find_one_tag("dynamics", axis2_node);
    if (dynamics_node)
    {
      // attempt to read 'Coulomb' friction
      boost::shared_ptr<const XMLTree> friction_node = find_one_tag("friction", dynamics_node);
      if (friction_node)
        joint->mu_fc = read_double(friction_node);

      // attempt to read viscous friction
      boost::shared_ptr<const XMLTree> damping_node = find_one_tag("damping", dynamics_node);
      if (damping_node)
        joint->mu_fv = read_double(damping_node);
    }

    // read joint limits
    boost::shared_ptr<const XMLTree> limit_node = find_one_tag("limits", axis2_node);
    if (limit_node)
    {
      // attempt to read the lower limit
      boost::shared_ptr<const XMLTree> llimit_node = find_one_tag("lower", limit_node);
      if (llimit_node)
        joint->lolimit[1] = read_double(llimit_node);

      // attempt to read the upper limit
      boost::shared_ptr<const XMLTree> ulimit_node = find_one_tag("upper", limit_node);
      if (ulimit_node)
        joint->hilimit[1] = read_double(ulimit_node);
    }
  }

  *output_begin++ = joint;
  return output_begin;
}

