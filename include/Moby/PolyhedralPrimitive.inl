template <class OutputIterator>
OutputIterator PolyhedralPrimitive::get_halfspaces(const Polyhedron& poly, boost::shared_ptr<const Ravelin::Pose3d> pose, const Ravelin::Transform3d& wTpose, OutputIterator output_begin)
{
  // get the normals and faces from polyhedron A, transforming to global frame
  for (unsigned i=0; i< poly.get_faces().size(); i++)
  {
    // get the face
    boost::shared_ptr<Polyhedron::Face> f = poly.get_faces()[i];

    // get the plane corresponding to the face; the plane will be in the
    // global frame, when it should really be defined in poseA
    Plane pu = f->get_plane();
    Ravelin::Vector3d pu_normal = pu.get_normal();
    assert(pu_normal.pose == GLOBAL); 
    pu_normal.pose = pose;
    pu.set_normal(pu_normal);

    // get the plane corresponding to the face, transformed to global frame
    Plane p = pu.transform(wTpose);

    // save the halfspace
    *output_begin++ = std::make_pair(p.get_normal(), p.offset);
  }

  return output_begin;
}

