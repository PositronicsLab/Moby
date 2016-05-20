template <class BidirectionalIterator>
void Constraint::insertion_sort(BidirectionalIterator first, BidirectionalIterator last)
{
  // exit if nothing to do
  if (first == last)
    return;

  BidirectionalIterator min = first;

  // loop
  BidirectionalIterator i = first;
  i++;
  for (; i != last; i++)
    if (*i < *min)
      min = i;

  // swap the iterators
  std::iter_swap(first, min);
  while (++first != last)
    for (BidirectionalIterator j = first; *j < *(j-1); --j)
      std::iter_swap((j-1), j);
}

template <class OutputIterator>
OutputIterator Constraint::get_super_bodies(OutputIterator output_begin) const
{
  boost::shared_ptr<Ravelin::DynamicBodyd> db1, db2;

  if (constraint_type == eContact)
  {
    // get the two single bodies
    boost::shared_ptr<Ravelin::SingleBodyd> s1 = contact_geom1->get_single_body();
    boost::shared_ptr<Ravelin::SingleBodyd> s2 = contact_geom2->get_single_body();

    // get the two super bodies
    boost::shared_ptr<Ravelin::DynamicBodyd> sb1 = s1->get_super_body();
    boost::shared_ptr<Ravelin::DynamicBodyd> sb2 = s2->get_super_body();
 
    // see whether to add the first super body
    if (s1->is_enabled())
      *output_begin++ = sb1;
    if (s2->is_enabled() && sb1 != sb2)
      *output_begin++ = sb2; 

    return output_begin;
  }
  else if (constraint_type == eLimit)
  {
    ArticulatedBodyPtr ab = boost::dynamic_pointer_cast<ArticulatedBody>(limit_joint->get_articulated_body());
    if (ab)
      *output_begin++ = ab;
    else
    {
      // get the inboard and outboard links
      boost::shared_ptr<Ravelin::RigidBodyd> inboard = limit_joint->get_inboard_link();
      boost::shared_ptr<Ravelin::RigidBodyd> outboard = limit_joint->get_outboard_link();

      // get the super bodies
      boost::shared_ptr<Ravelin::DynamicBodyd> sb1 = inboard->get_super_body();
      boost::shared_ptr<Ravelin::DynamicBodyd> sb2 = outboard->get_super_body();
 
      // see whether to add the super bodies
      if (inboard->is_enabled())
        *output_begin++ = sb1;
      if (outboard->is_enabled() && sb1 != sb2)
        *output_begin++ = sb2; 
    }

    return output_begin;
  }
  else if (constraint_type == eInverseDynamics)
  {
    ArticulatedBodyPtr ab = boost::dynamic_pointer_cast<ArticulatedBody>(inv_dyn_joint->get_articulated_body());
    if (ab)
      *output_begin++ = ab;
    else
    {
      // get the inboard and outboard links
      boost::shared_ptr<Ravelin::RigidBodyd> inboard = inv_dyn_joint->get_inboard_link();
      boost::shared_ptr<Ravelin::RigidBodyd> outboard = inv_dyn_joint->get_outboard_link();

      // get the super bodies
      boost::shared_ptr<Ravelin::DynamicBodyd> sb1 = inboard->get_super_body();
      boost::shared_ptr<Ravelin::DynamicBodyd> sb2 = outboard->get_super_body();
 
      // see whether to add the super bodies
      if (inboard->is_enabled())
        *output_begin++ = sb1;
      if (outboard->is_enabled() && sb1 != sb2)
        *output_begin++ = sb2; 
    }

    return output_begin;
  }
  else if (constraint_type == eSpringDamper)
  {
    ArticulatedBodyPtr ab = boost::dynamic_pointer_cast<ArticulatedBody>(spring_damper_joint->get_articulated_body());
    if (ab)
      *output_begin++ = ab;
    else
    {
      // get the inboard and outboard links
      boost::shared_ptr<Ravelin::RigidBodyd> inboard = spring_damper_joint->get_inboard_link();
      boost::shared_ptr<Ravelin::RigidBodyd> outboard = spring_damper_joint->get_outboard_link();

      // get the super bodies
      boost::shared_ptr<Ravelin::DynamicBodyd> sb1 = inboard->get_super_body();
      boost::shared_ptr<Ravelin::DynamicBodyd> sb2 = outboard->get_super_body();
 
      // see whether to add the super bodies
      if (inboard->is_enabled())
        *output_begin++ = sb1;
      if (outboard->is_enabled() && sb1 != sb2)
        *output_begin++ = sb2; 
    }

    return output_begin;
  }
  else if (constraint_type == eImplicitJoint)
  {
    // get the inboard and outboard links
    boost::shared_ptr<Ravelin::RigidBodyd> inboard = implicit_joint->get_inboard_link();
    boost::shared_ptr<Ravelin::RigidBodyd> outboard = implicit_joint->get_outboard_link();

    // get the super bodies
    boost::shared_ptr<Ravelin::DynamicBodyd> sb1 = inboard->get_super_body();
    boost::shared_ptr<Ravelin::DynamicBodyd> sb2 = outboard->get_super_body();
 
    // see whether to add the super bodies
    if (inboard->is_enabled())
      *output_begin++ = sb1;
    if (outboard->is_enabled() && sb1 != sb2)
      *output_begin++ = sb2; 

    return output_begin;
  }
  else
    assert(false);
} 


