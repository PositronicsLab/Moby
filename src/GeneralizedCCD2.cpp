/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Does a collision check for a pair of geometries (serial version)
/**
 * \param a the first geometry
 * \param b the second geometry
 * \param aTb the transform from b's frame to a's frame
 * \param bTa the transform from a's frame to b's frame
 * \param vels linear and angular velocities of bodies
 * \param contacts on return
 */
void GeneralizedCCD::check_geoms(double dt, CollisionGeometryPtr a, const PosePair& a_poses, CollisionGeometryPtr b, const PosePair& b_poses, vector<Event>& contacts)
{
  map<BVPtr, vector<const Point3d*> > a_to_test, b_to_test;

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::check_geoms() entered" << endl;
  FILE_LOG(LOG_COLDET) << "  checking geometry " << a->id << " for body " << a->get_single_body()->id << std::endl;
  FILE_LOG(LOG_COLDET) << "  against geometry " << b->id << " for body " << b->get_single_body()->id << std::endl;

  // init statistics for geometry pair
  unsigned n_bv_tests = 0;
  unsigned n_verts_tested = 0;

  // set the earliest TOC
  double earliest = (double) 1.0; 

  // setup the contact set for these two geometries 
  vector<Event> local_contacts;

  // get the primitives for a and b
  PrimitivePtr aprimitive = a->get_geometry();
  PrimitivePtr bprimitive = b->get_geometry();

  // get the two top-level BVs for a and b
  BVPtr bv_a = aprimitive->get_BVH_root(a);
  BVPtr bv_b = bprimitive->get_BVH_root(b); 

  // set poses for the collision geometries
  *a->_F = a_poses.t0;
  *b->_F = b_poses.t0;

  // compute transform between a and b at t0
  Transform3d aTb = Pose3d::calc_relative_pose(b->get_pose(), a->get_pose());

  // add the two top-level BVs to the queue for processing
  queue<BVProcess> q;
  q.push(BVProcess());
  q.back().bva = get_swept_BV(a, bv_a, a_poses);
  q.back().bvb = get_swept_BV(b, bv_b, b_poses);
  q.back().nexp = 0;
  q.back().ax = bv_a;
  q.back().bx = bv_b;

  // process until the queue is empty
  while (!q.empty())
  {
    // get the number of expansions  
    unsigned nexp = q.front().nexp;

    // get the BVs off of the queue
    BVPtr bva = q.front().bva;
    BVPtr bvb = q.front().bvb;
    BVPtr ax = q.front().ax;
    BVPtr bx = q.front().bx;

    // increment the # of bounding volume tests
    n_bv_tests++;

    // if the velocity-expanded BVs do not intersect, continue looping
    if (!BV::intersects(bva, bvb, aTb))
    {
      FILE_LOG(LOG_COLDET) << " -- bounding volumes do not intersect" << endl;
      q.pop();
      continue;
    }

    // calculate the volume for the OBBs
    double ax_vol = ax->calc_volume();
    double bx_vol = bx->calc_volume();

    // velocity expanded BVs do intersect; if both BVs are leafs OR maximum
    // depth has been reached, intersect
    // the vertices of triangles in one BV against the triangles of the other
    // (and vice versa)
    if ((ax->is_leaf() && bx->is_leaf()) || nexp >= _max_dexp)
    {
      // get testing vertex vectors for ax and bx
      vector<pair<const Point3d*, const Point3d*> >& ax_test = a_to_test[ax];
      vector<pair<const Point3d*, const Point3d*> >& bx_test = b_to_test[bx];

      // get the sets of edges for ax and bx
      aprimitive->get_edges(ax, bx_test);
      bprimitive->get_edges(bx, ax_test);
    }
    // descend into BV with greater volume
    else if ((ax_vol > bx_vol && !ax->is_leaf()) || bx->is_leaf())
    {
      // add all children of bva to the queue for processing
      BOOST_FOREACH(BVPtr bv, ax->children)
      {
        q.push(BVProcess());
        q.back().bva = get_swept_BV(a, bv, a_poses);
        q.back().bvb = bvb; 
        q.back().nexp = nexp+1;
        q.back().ax = bv;
        q.back().bx = bx;
      }
    }
    else
    {
      // add all children of bvb to the queue for processing
      BOOST_FOREACH(BVPtr bv, bx->children)
      {
        q.push(BVProcess());
        q.back().bva = bva; 
        q.back().bvb = get_swept_BV(b, bv, b_poses);
        q.back().nexp = nexp+1;
        q.back().ax = ax;
        q.back().bx = bv;
      }
    }

    // pop the element off of the front of the queue
    q.pop();
  }

  // make vectors of edges unique
  for (map<BVPtr, vector<pair<const Point3d*, const Point3d*> >::iterator i = a_to_test.begin(); i != a_to_test.end(); i++)
  {
    std::sort(i->second.begin(), i->second.end());
    i->second.erase(std::unique(i->second.begin(), i->second.end()), i->second.end());
  }
  for (map<BVPtr, vector<pair<const Point3d*, const Point3d*> > >::iterator i = b_to_test.begin(); i != b_to_test.end(); i++)
  {
    std::sort(i->second.begin(), i->second.end());
    i->second.erase(std::unique(i->second.begin(), i->second.end()), i->second.end());
  }

  // call check_edges on bounding volumes from geometry b
  for (map<BVPtr, vector<pair<const Point3d*, const Point3d*> > >::iterator i = b_to_test.begin(); i != b_to_test.end(); i++)
  {
    n_verts_tested += i->second.size();
    check_edges(dt, a, b, i->first, i->second, a_poses, b_poses, earliest, local_contacts);
  }

  // call check_edges on bounding volumes from geometry a
  for (map<BVPtr, vector<pair<const Point3d*, const Point3d*> > >::iterator i = a_to_test.begin(); i != a_to_test.end(); i++)
  {
    n_verts_tested += i->second.size();
    check_edges(dt, b, a, i->first, i->second, b_poses, a_poses, earliest, local_contacts);
  }

  if (LOGGING(LOG_COLDET))
  {
    SingleBodyPtr sba = a->get_single_body(); 
    SingleBodyPtr sbb = b->get_single_body(); 
    FILE_LOG(LOG_COLDET) << " -- contact points for bodies " << sba->id << " and " << sbb->id << endl;
    for (unsigned i=0; i< local_contacts.size(); i++)
      FILE_LOG(LOG_COLDET) << local_contacts[i] << std::endl;
  }

  // get the time-of-impact tolerance
  shared_ptr<EventDrivenSimulator> sim(simulator);
  const double TOI_TOLERANCE = std::numeric_limits<double>::epsilon();

  // sort the vector of contacts
  std::sort(local_contacts.begin(), local_contacts.end());

  // see what points to insert into the global set of contacts 
  // we return all contacts if using an event-driven method to deal with Zeno
  // points 
  if (!return_all_contacts)
  {
    for (unsigned i=0; i< local_contacts.size(); i++)
      if (local_contacts[i].t > earliest + TOI_TOLERANCE/dt)
      {
        contacts.insert(contacts.end(), local_contacts.begin(), local_contacts.begin()+i);
        break;
      }
  }
  else
  {
    #ifdef _OPENMP
    pthread_mutex_lock(&_contact_mutex);
    #endif

    // insert the contacts from these geometries into the contact map of all
    // geometries
    contacts.insert(contacts.end(), local_contacts.begin(), local_contacts.end());

    #ifdef _OPENMP
    pthread_mutex_unlock(&_contact_mutex);
    #endif
  }

  if (LOGGING(LOG_COLDET))
  {
    // determine how many pairs of bounding volumes
    std::vector<BVPtr> all_bvs;
    bv_a->get_all_BVs(std::back_inserter(all_bvs));
    unsigned n_a_bvs = all_bvs.size();
    all_bvs.clear();
    bv_b->get_all_BVs(std::back_inserter(all_bvs));
    unsigned n_b_bvs = all_bvs.size();
    FILE_LOG(LOG_COLDET) << "  number of BV vs. BV tests: " << n_bv_tests << "/" << (n_a_bvs * n_b_bvs) << std::endl;

    // output number of vertices tested 
    FILE_LOG(LOG_COLDET) << "  number of vertices tested: " << n_verts_tested << std::endl;
  }
  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::check_geoms() exited" << endl;
}

/// Checks a set of edges of geometry a against geometry b
void GeneralizedCCD::check_edges(double dt, CollisionGeometryPtr a, CollisionGeometryPtr b, BVPtr bvb, const std::vector<pair<const Point3d*, const Point3d*> >& a_edges, const PosePair& a_poses, const PosePair& b_poses, double& earliest, vector<Event>& local_contacts) const
{
  Vector3d normal;
  Point3d point;

  // get the poses
  const Pose3d& Pa_t0 = a_poses.t0;
  const Pose3d& Pa_tf = a_poses.tf;
  const Pose3d& Pb_t0 = b_poses.t0;
  const Pose3d& Pb_tf = b_poses.tf;

  // get the time-of-impact tolerance
  shared_ptr<EventDrivenSimulator> sim(simulator);
  const double TOI_TOLERANCE = std::numeric_limits<double>::epsilon();

  // get the two bodies
  if (LOGGING(LOG_COLDET))
  {
    RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(a->get_single_body()); 
    RigidBodyPtr rbb = dynamic_pointer_cast<RigidBody>(b->get_single_body()); 
    FILE_LOG(LOG_COLDET) << "  -- checking body " << rba->id << " against " << rbb->id << endl;
  }

  // populate the DStruct for checking vertices of a against b
  DStruct ds;
  populate_dstruct(&ds, a, Pa_t0, Pa_tf, b, Pb_t0, Pb_tf, bvb);

  // setup a "queue" for checking edges 
  vector<pair<double, Point3d> > Q;
  Q.clear();
  BOOST_FOREACH(pair<const Point3d*, const Point3d*>& e, a_edges)
  {
    // get edge in proper frame
    Point3d va1 = *e.first;
    Point3d va2 = *e.second;

    // get point in b's frame at time 0 (for distance sorting) 
    Point3d u1_b = ds.bTa_t0.transform_point(va1);
    Point3d u2_b = ds.bTa_t0.transform_point(va2);

    // we'll sort on inverse distance from the center of mass (origin of b frame) 
    double dist = 1.0/(u1_b.norm_sq() + u2_b.norm_sq());

    // push the edges onto the queue
    // NOTE: assumes that center of geometry of body a is its C.O.M.
    //       (if the assumption is wrong, this will only make things slower)
    Q.push_back(make_pair(dist, *e));
  }

  // sort the queue
  std::sort(Q.begin(), Q.end());

  // check all vertices of a against b
  while (!Q.empty())
  {
    // setup u in the DStruct
    ds.u1_a = Q.back().second.first;
    ds.u2_a = Q.back().second.second;

    if (LOGGING(LOG_COLDET))
    {
      RigidBodyPtr rba = dynamic_pointer_cast<RigidBody>(a->get_single_body()); 
      FILE_LOG(LOG_COLDET) << "    -- checking vertex u1: " << ds.u1_a << " of " << rba->id << endl; 
      FILE_LOG(LOG_COLDET) << "    -- checking vertex u2: " << ds.u2_a << " of " << rba->id << endl; 
      FILE_LOG(LOG_COLDET) << "     -- u1 global pos (t0): " << Pb_t0.transform_point(ds.u1_a) << endl;
      FILE_LOG(LOG_COLDET) << "     -- u2 global pos (t0): " << Pb_t0.transform_point(ds.u2_a) << endl;
    }

    // determine TOI, if any
    // NOTE: this has been modified to account for the method of dealing with
    // Zeno points in ImpulseContactSimulator; it is considerably slower as a
    // result
    double toi;
    if (!return_all_contacts)
      toi = determine_TOI_fast(0.0, earliest, &ds, point, normal);
    else
      toi = determine_TOI_fast(0.0, (double) 1.0, &ds, point, normal);
  
    // insert into the contacts set if the TOI is finite 
    if (toi < std::numeric_limits<double>::max())
    {
      // insert the contacts
      local_contacts.push_back(create_contact(toi, a, b, point, normal));
      if (toi < earliest)
        earliest = std::min(toi + TOI_TOLERANCE/dt, (double) 1.0);
    }

    // move onto the next vertex
    Q.pop_back();
  }
}

/****************************************************************************
 Methods for Drumwright-Shell algorithm begin 
****************************************************************************/

/// Implements DETERMINE-TOC() [very fast, fully linearized version]
/**
 * \param t0 the time before integration
 * \param tf the time after integration
 * \param dstruct structure passed to determine_TOI of precomputed data
 * \param pt the contact point (if any) in world coordinates
 * \param normal the normal (if any) in world coordinates
 * \return the time of impact [0, 1] such that t0 = 0 and t1 = 1; if the time 
 *         of impact is greater than 1.0, then the time of impact is not in the 
 *         interval [t0, tf] and should be discarded
 * \pre rigid body states are at time t
 */
double GeneralizedCCD::determine_TOI_fast(double t0, double tf, const DStruct* ds, Point3d& cp, Vector3d& normal) const
{
  const double INF = std::numeric_limits<double>::max();
  const unsigned X = 0, Y = 1, Z = 2;
  OBB O;
  Vector3d nalpha, nbeta, ngamma;

  // get the static collision geometry
  CollisionGeometryPtr gb = ds->gb;

  // get the primitive for gs (the body treated as static) 
  shared_ptr<const Primitive> gb_primitive = gb->get_geometry();

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI() entered" << endl;
  FILE_LOG(LOG_COLDET) << "  time t0: " << t0 << endl;
  FILE_LOG(LOG_COLDET) << "  time tf: " << tf << endl;

  // get edge in ga's frame
  const Point3d& u1_a = ds->u1_a;
  const Point3d& u2_a = ds->u2_a;

  // get the BV for gb
  BVPtr gb_BV = ds->b_BV; 
  assert(gb->get_pose() == gb_BV->get_relative_pose());

  // get useful poses
  const Pose3d& Pb_t0 = ds->Pb_t0;
  const Pose3d& Pb_tf = ds->Pb_tf;

  // get u- in gb's frame- at times t0 and tf
  Point3d u1_0 = ds->bTa_t0.transform_point(u1_a);
  Point3d u1_f = ds->bTa_tf.transform_point(u1_a);
  Point3d u2_0 = ds->bTa_t0.transform_point(u2_a);
  Point3d u2_f = ds->bTa_tf.transform_point(u2_a);

// TODO: here is where we check for intersection
/*
  // check for intersection 
  double t;
  if (gb_primitive->intersect_seg(gb_BV, LineSeg3(u1_0, u1_f), t, cp, normal))
  {
    FILE_LOG(LOG_COLDET) << "  intersection detected!  time of impact: " << t << " (true: " << (t0 + (tf-t0)*t) << ")" << endl;

    // transform time to account for ta and tb
    t = t0 + (tf-t0)*t;

    FILE_LOG(LOG_COLDET) << "     -- point intersection (untransformed): " << cp << endl;
    FILE_LOG(LOG_COLDET) << "     -- normal (untransformed): " << normal << endl;

    // since all calculations are in gs's frame; interpolate gs to time t
    // and transform contact point and normal to global coordinates
    Origin3d cpo(cp);
    Origin3d no(normal);
    double s = t/(tf-t0);
    cp = Pose3d::interpolate_transform_point(Pb_t0, Pb_tf, s, cpo);
    normal = Pose3d::interpolate_transform_vector(Pb_t0, Pb_tf, s, no);

    // look for degenerate normal
    if (std::fabs(normal.norm() - (double) 1.0) > NEAR_ZERO)
    {
      FILE_LOG(LOG_COLDET) << "    -- degenerate normal detected! (" << normal << "); not reporting intersection" << endl;
      return INF;
    }        
    
    FILE_LOG(LOG_COLDET) << "     -- point intersection: " << cp << endl;
    FILE_LOG(LOG_COLDET) << "     -- normal: " << normal << endl;
    FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI_fast() exited" << endl;

    return t;
  }
*/

  FILE_LOG(LOG_COLDET) << "  no impact detected" << endl;
  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI_fast() exited" << endl;

  // still here?  no contact...
  return INF; 
}

/// Implements DETERMINE-TOC()
/**
 * \param t0 the time before integration
 * \param tf the time after integration
 * \param dstruct structure passed to determine_TOI of precomputed data
 * \param cp the contact point (if any) in world coordinates
 * \param normal the normal (if any) in world coordinates
 * \return the time of impact [0, 1] such that t0 = 0 and t1 = 1; if the time 
 *         of impact is greater than 1.0, then the time of impact is not in the 
 *         interval [t0, tf] and should be discarded
 * \pre rigid body states are at time t
 *
 * The math behind this function (I use LaTeXiT to visualize) follows:
 *
 * We are interested in testing the path traced by a point on body B against
 * the geometry of another body S. We will do all calculations relative to
 * S's (moving) frame so that an observer at S can treat himself as stationary.
 *
 * Notation/defs: v relative linear velocity
 *                \theta relative angular velocity
 *                r vector from center-of-mass of B to point u
 *                superscricp (A) to right (vector defined in A's frame
 *                superscricps (A,B) to left and right (transform from B to A)
 *                subscricp (A) to right (vector for body A)
 * 
 * The velocity of point u on B (in S's frame) is:
 * \dot{u}^S(t)& = v^S + \theta^S \times ~^ST^0(t) [u^0(t)-x_B^0(t)]\\
 * & = v^S + \theta^S \times ~^ST^0(t) [x_B^0(t) +~^0R^B(t) r^B-x_B^0(t)]\\
 * & = v^S + \theta^S \times~^SR^0(t)~^0R^B(t)r^B
 *
 * We define ^SR^0(t) using a first-order approximation as follows:
 * ^SR^0(t)& = [~^0R^S(t_0) + (t_f-t_0)\cdot \omega_S^0\times ~^0R^S(t_0)]^T\\
 * & = ~^SR^0(t_0) - (t_f-t_0)\cdot~^SR^0(t_0)\omega_S^0\times
 *
 * Similarly, ^0R^B(t) is defined as:
 * ^0R^B(t)& = ~^0R^B(t_0) + (t_f-t_0)\cdot\omega_B^0\times~^0R^B(t_0)
 *
 */
double GeneralizedCCD::determine_TOI(double t0, double tf, const DStruct* ds, Point3d& cp, Vector3d& normal) const
{
  const double INF = std::numeric_limits<double>::max();
  const unsigned X = 0, Y = 1, Z = 2;
  OBB O;
  Vector3d nalpha, nbeta, ngamma;

  // init bisection statistics
  unsigned nbisects = 0;

  // get the static collision geometry
  CollisionGeometryPtr gb = ds->gb;

  // get the primitive for gb (the body treated as static) 
  shared_ptr<const Primitive> gb_primitive = gb->get_geometry();

  // get useful poses
  const Pose3d& Pb_t0 = ds->Pb_t0;
  const Pose3d& Pb_tf = ds->Pb_tf;

  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI() entered" << endl;
  FILE_LOG(LOG_COLDET) << "  time t0: " << t0 << endl;
  FILE_LOG(LOG_COLDET) << "  time tf: " << tf << endl;

  // get the BV for gb
  BVPtr gb_BV = ds->b_BV;

  // get transforms between ga's frame and gb's frame
  const Transform3d& bTa_t0 = ds->bTa_t0;
  const Transform3d& bTa_tf = ds->bTa_tf;

  // get edge in ga's frame (this vector remains constant over time) at t0
  const Vector3d& u1 = ds->u1_a;
  const Vector3d& u2 = ds->u2_a;
  Point3d u1_0_b = bTa_t0.transform_point(u1);
  Point3d u2_0_b = bTa_t0.transform_point(u2);

  // get u in gb's frame at tf 
  Point3d u1_f_b = bTa_tf.transform_point(u1);
  Point3d u2_f_b = bTa_tf.transform_point(u2);

// TODO: fix this
  // setup the axes of the bounding box
  Vector3d u0uf = uf_b - u0_b;
  double norm_u0uf = u0uf.norm();
  if (norm_u0uf < std::numeric_limits<double>::epsilon())
  {
    // arbitrary bounding box
    nalpha = Vector3d(1,0,0, gb->get_pose());
    nbeta = Vector3d(0,1,0, gb->get_pose());
    ngamma = Vector3d(0,0,1, gb->get_pose());
    O.R = Matrix3d::identity();
  }
  else
  {
    // primary axis of bounding box aligned w/u0uf
    nalpha = u0uf/norm_u0uf;
    Vector3d::determine_orthonormal_basis(nalpha, nbeta, ngamma);
    O.R.set_column(X, nalpha);
    O.R.set_column(Y, nbeta);
    O.R.set_column(Z, ngamma);
  }

  // setup quaternion endpoints of interpolation
  const Quatd& q0 = ds->q0;
  const Quatd& qf = ds->qf; 

  // NOTE: collision geometry (and therefore BV) should probably be updated to proper time

  // determine whether minimum/maximum deviation is between two planes;
  // if so, determine the interpolation value that yields the minimum
  // and maximum deviation
  Vector3d normal1(gb->get_pose()), normal2(gb->get_pose());
  double rho1_max = -1.0, rho1_min = -1.0, rho2_max = -1.0, rho2_min = -1.0, rho3_max = -1.0, rho3_min = -1.0;
  double max_d1 = -INF, max_d2 = -INF, max_d3 = -INF;
  double min_d1 = INF, min_d2 = INF, min_d3 = INF;
  if (bound_u(u, q0, qf, normal1, normal2))
  {
    if (nalpha.dot(normal1)*nalpha.dot(normal2) > 0)
    {
      max_d1 = calc_max_dev(u, nalpha, q0, qf, rho1_max);
      min_d1 = calc_min_dev(u, nalpha, q0, qf, rho1_min);

      // scale rho values to (tf-t0)
      rho1_max *= (tf-t0);
      rho1_min *= (tf-t0);
    }
    if (nbeta.dot(normal1)*nbeta.dot(normal2) > 0)
    {
      max_d2 = calc_max_dev(u, nbeta, q0, qf, rho2_max);
      min_d2 = calc_min_dev(u, nbeta, q0, qf, rho2_min);

      // scale rho values to (tf-t0)
      rho2_max *= (tf-t0);
      rho2_min *= (tf-t0);
    }
    if (ngamma.dot(normal1)*ngamma.dot(normal2) > 0)
    {
      max_d3 = calc_max_dev(u, ngamma, q0, qf, rho3_max);
      min_d3 = calc_min_dev(u, ngamma, q0, qf, rho3_min);

      // scale rho values to (tf-t0)
      rho3_max *= (tf-t0);
      rho3_min *= (tf-t0);
    }
  }

  // setup the queue
  typedef pair<double, double> RPair;
  priority_queue<RPair, vector<RPair>, std::greater<RPair> > Q;
  Q.push(make_pair(t0, tf));

  // process until error sufficiently low
  while (!Q.empty())
  {
    // get the element off of the top of the queue
    double tx = Q.top().first;
    double ty = Q.top().second;
    Q.pop();

    // setup delta t
    const double dt = ty - tx;
    assert(dt > 0.0);

    // interpolation parameter ranges from [0,1]; 0 corresponds to t0, 
    // 1 corresponds to tf.  Determine what tx, ty correspond to
    const double sx = tx/(tf-t0);
    const double sy = ty/(tf-t0);

    // determine point u at times tx and ty in frame b 
    Point3d ux(Transform3d::interpolate_transform_point(bTa_t0, bTa_tf, sx, Origin3d(ds->u_a)), gb->get_pose());
    Point3d uy(Transform3d::interpolate_transform_point(bTa_t0, bTa_tf, sy, Origin3d(ds->u_a)), gb->get_pose());

    FILE_LOG(LOG_COLDET) << " -- checking segment for time [" << tx << ", " << ty << "]" << endl;
    FILE_LOG(LOG_COLDET) << "  p(" << tx << ") = " << ux << "  p(" << ty << ") ~= " << uy << endl;

    // init deviation maxima/minima
    double dp_alpha = -INF, dp_beta = -INF, dp_gamma = -INF;
    double dn_alpha = INF, dn_beta = INF, dn_gamma = INF;

    // see whether this interval contains a minimum/maximum 
    if (rho1_max >= tx && rho1_max <= ty) dp_alpha = max_d1;
    if (rho1_min >= tx && rho1_min <= ty) dn_alpha = min_d1;
    if (rho2_max >= tx && rho2_max <= ty) dp_beta = max_d2;
    if (rho2_min >= tx && rho2_min <= ty) dn_beta = min_d2;
    if (rho3_max >= tx && rho3_max <= ty) dp_gamma = max_d3;
    if (rho3_min >= tx && rho3_min <= ty) dn_gamma = min_d3;
    
    // calculate deviation at endpoints
    pair<double, double> deva = calc_deviations(u, nalpha, q0, qf, sx, sy);
    pair<double, double> devb = calc_deviations(u, nbeta, q0, qf, sx, sy);
    pair<double, double> devg = calc_deviations(u, ngamma, q0, qf, sx, sy);
  
    // set deviation maxima/minima    
    dn_alpha = std::min(dn_alpha, deva.first);
    dn_beta = std::min(dn_beta, devb.first);
    dn_gamma = std::min(dn_gamma, devg.first);
    dp_alpha = std::max(dp_alpha, deva.second);
    dp_beta = std::max(dp_beta, devb.second);
    dp_gamma = std::max(dp_gamma, devg.second);

    // determine the half deviations in each direction
    double len_uxy = (uy - ux).norm();
    double d_alpha = std::max(std::fabs(dp_alpha - dn_alpha), len_uxy*(double) 0.5);
    double d_beta = std::fabs(dp_beta - dn_beta);
    double d_gamma = std::fabs(dp_gamma - dn_gamma);

    // setup the bounding box
    double nu = (d_beta + d_gamma + d_alpha)*2.0 - len_uxy;
    assert(nu > -NEAR_ZERO);
    FILE_LOG(LOG_COLDET) << " -- dalpha: [" << dn_alpha << ", " << dp_alpha << "]" << endl;
    FILE_LOG(LOG_COLDET) << " -- dbeta: [" << dn_beta << ", " << dp_beta << "]" << endl;
    FILE_LOG(LOG_COLDET) << " -- dgamma: [" << dn_gamma << ", " << dp_gamma << "]" << endl;
    FILE_LOG(LOG_COLDET) << " -- nu contribution along alpha/beta: " << nu << endl;
    Point3d cxy = (ux + uy)*0.5;
    Point3d mini = cxy - nalpha*d_alpha - nbeta*d_beta - ngamma*d_gamma;
    Point3d maxi = cxy + nalpha*d_alpha + nbeta*d_beta + ngamma*d_gamma;
    O.center = (maxi+mini)*0.5;
    O.l = O.R.transpose_mult(Origin3d((maxi-mini)*0.5));
    O.l[0] = std::fabs(O.l[0]);
    O.l[1] = std::fabs(O.l[1]);
    O.l[2] = std::fabs(O.l[2]);

    // determine whether there is an intersection with the bounding box for b
    FILE_LOG(LOG_COLDET) << " -- nu: " << nu << endl;
    FILE_LOG(LOG_COLDET) << " -- checking BV intersection of: " << endl;
    FILE_LOG(LOG_COLDET) << O << " and " << endl << gb_BV;
    if (LOGGING(LOG_COLDET))
    {
      OBBPtr gb_OBB = dynamic_pointer_cast<OBB>(gb_BV);
      if (gb_OBB)
        FILE_LOG(LOG_COLDET) << *gb_OBB << endl;
    }

    // O is in gb_BV's frame
    if (!BV::intersects(&O, gb_BV.get()))
    {
      FILE_LOG(LOG_COLDET) << "   -- no intersection; continuing looping..." << endl;
      continue;
    }

    // ************************************************************************
    // there is a bounding box intersection; bisect or intersect with triangles 
    // ************************************************************************

    if (nu > eps_tolerance)
    {
      FILE_LOG(LOG_COLDET) << "   -- intersection detected; trajectory segment must be bisected" << endl;

      // add two elements to the queue
      double ti = (tx+ty)*0.5;
      Q.push(make_pair(tx, ti));
      Q.push(make_pair(ti, ty));
      nbisects++;
    }
    else
    {
      FILE_LOG(LOG_COLDET) << "   -- intersection detected and nu less than tolerance" << endl;

      // intersect the line segment with the geometry
      double t;
      if (gb_primitive->intersect_seg(gb_BV, LineSeg3(ux, uy), t, cp, normal))
      {
        FILE_LOG(LOG_COLDET) << "  intersection detected!  time of impact: " << t << " (true: " << (tx + (ty-tx)*t) << ")" << endl;

        // transform time to account for tx and ty
        t = tx + (ty-tx)*t;

        FILE_LOG(LOG_COLDET) << "     -- point intersection (untransformed): " << cp << endl;
        FILE_LOG(LOG_COLDET) << "     -- normal (untransformed): " << normal << endl;

        // since all calculations are in gb's frame; interpolate gb to time t
        // and transform contact point and normal to global coordinates
        Origin3d cpo(cp);
        Origin3d no(normal);
        double s = t/(tf-t0);
        cp = Pose3d::interpolate_transform_point(Pb_t0, Pb_tf, s, cpo);
        normal = Pose3d::interpolate_transform_vector(Pb_t0, Pb_tf, s, no);

        // look for degenerate normal
        if (std::fabs(normal.norm() - (double) 1.0) > NEAR_ZERO)
        {
          FILE_LOG(LOG_COLDET) << "    -- degenerate normal detected! (" << normal << "); not reporting intersection" << endl;
          continue;
        }        
        
        FILE_LOG(LOG_COLDET) << "     -- point intersection: " << cp << endl;
        FILE_LOG(LOG_COLDET) << "     -- normal: " << normal << endl;
        FILE_LOG(LOG_COLDET) << "# of bisections: " << nbisects << endl;
        FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI() exited" << endl;

        return t;
      }
      else
      {
        FILE_LOG(LOG_COLDET) << "     -- intersection does not occur in interval [t0, tf]" << endl;
        FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI() exited" << endl;
      }
    }
  }

  FILE_LOG(LOG_COLDET) << "# of bisections: " << nbisects << endl;
  FILE_LOG(LOG_COLDET) << "  no impact detected" << endl;
  FILE_LOG(LOG_COLDET) << "GeneralizedCCD::determine_TOI() exited" << endl;

  // still here?  no contact...
  return INF; 
}

/// Determines the two planes that bound a vector between two rotations
bool GeneralizedCCD::bound_u(const Vector3d& u, const Quatd& q0, const Quatd& qf, Vector3d& normal1, Vector3d& normal2)
{
  // compute the rotated u
  Vector3d u0(q0*Origin3d(u), normal1.pose);
  Vector3d uf(qf*Origin3d(u), normal1.pose);

  // determine a vector perpendicular to both
  Vector3d perp = Vector3d::cross(u0, uf);
  double perp_norm = perp.norm();
  if (perp_norm < NEAR_ZERO)
    return false;
  else
    perp /= perp_norm;

  // determine the two planes
  normal1 = Vector3d::normalize(Vector3d::cross(u0, perp));
  normal2 = Vector3d::normalize(Vector3d::cross(uf, perp));

  // make sure that both vectors are on the same side of the plane
  double dot1 = u0.dot(normal2);
  double dot2 = uf.dot(normal1);
  if (dot1*dot2 < 0)
    normal2 = -normal2;
 
  return true;
}

/// Function for computing the maximum deviation in a given direction
double GeneralizedCCD::calc_deviation(double alpha, void* params)
{
  const DeviationCalc& data = *((const DeviationCalc*) params);

  // linearly interpolate the quaternions
  Quatd q = Quatd::lerp(data.q1, data.q2, alpha);

  // do the arithmetic
  return Vector3d::dot(data.d, Vector3d(q*Origin3d(data.u), data.d.pose));
}

/// Computes the minimum deviation over interval [0,1] using bracketing and Brent's method
/**
 * \param u vector from center-of-mass to point (in local frame)
 * \param d the direction vector
 * \param q1 the orientation at t=0
 * \param q2 the orientation at t=1
 * \param t on return, contains the t at which the deviation is minimized
 * \return the minimum deviation
 */
double GeneralizedCCD::calc_min_dev(const Vector3d& u, const Vector3d& d, const Quatd& q1, const Quatd& q2, double& t)
{
  const double TOL = std::sqrt(std::sqrt(std::numeric_limits<double>::epsilon()));

  double ft;

  // setup data for deviation calculations
  DeviationCalc dc;
  dc.u = u;
  dc.q1 = q1;
  dc.q2 = q2;
  dc.d = d;

  // find suitable starting value for brent's method
  double f0 = calc_deviation(0.0, &dc);
  double f1 = calc_deviation(1.0, &dc);
  t = (f0 < f1) ? 0.0 : 1.0;
  
  // call Brent's method -- note: tolerance is a bit high
  if (!brent(0.0, 1.0, t, ft, &calc_deviation, TOL, &dc))
    cerr << "GeneralizedCCD::calc_min_dev() - brent's method failed!" << endl;

  return ft;
}

/// Computes the maximum and minimum deviation in a given direction for two points
pair<double, double> GeneralizedCCD::calc_deviations(const Vector3d& u, const Vector3d& d, const Quatd& q1, const Quatd& q2, double t1, double t2)
{
  // compute the deviation in the two directions
  DeviationCalc dc;
  dc.u = u;
  dc.d = d;
  dc.q1 = q1;
  dc.q2 = q2;
  double dev1 = calc_deviation(t1, &dc);
  double dev2 = calc_deviation(t2, &dc);
  pair<double, double> dev(dev1, dev2);
  if (dev1 > dev2)
    std::swap(dev.first, dev.second);
  return dev;
}

/// Computes the maximum deviation over interval [0,1] using bracketing and Brent's method
/**
 * \param u vector from center-of-mass to point (in local frame)
 * \param d the direction vector
 * \param q1 the orientation at t=0
 * \param q2 the orientation at t=1
 * \param t on return, contains the t at which the deviation is maximized [0,1]
 * \return the maximum deviation
 */
double GeneralizedCCD::calc_max_dev(const Vector3d& u, const Vector3d& d, const Quatd& q1, const Quatd& q2, double& t)
{
  const double TOL = 1e-2;

  double ft;

  // setup data for deviation calculations
  DeviationCalc dc;
  dc.u = u;
  dc.q1 = q1;
  dc.q2 = q2;

  // set d in the deviation data because brent's method is a minimizer and we
  // want to maximize
  dc.d = -d;

  // find suitable starting value for brent's method
  double f0 = calc_deviation(0.0, &dc);
  double f1 = calc_deviation(1.0, &dc);
  t = (f0 < f1) ? 0.0 : 1.0;
  
  // call Brent's method -- note: tolerance is a bit high
  if (!brent(0.0, 1.0, t, ft, &calc_deviation, TOL, &dc))
    cerr << "GeneralizedCCD::calc_max_dev() - brent's method failed!" << endl;

  return -ft;
}

/****************************************************************************
 Methods for Drumwright-Shell algorithm end 
****************************************************************************/


