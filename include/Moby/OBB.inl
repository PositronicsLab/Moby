/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

// For outputting description of OBB (primarily for debugging purposes)
inline std::ostream& operator<<(std::ostream& out, const OBB& o)
{
  out << " (address): " << &o << std::endl;
  out << " center: " << o.center << std::endl;
  out << " axes: " << std::endl << o.R;
  out << "   (as axis-angle): " << Ravelin::AAngled(o.R) << std::endl;
  out << " half lengths: " << o.l << "  1/8 volume: " << o.calc_volume() << std::endl;
  out << " children:";
  BOOST_FOREACH(BVPtr child, o.children)
    out << " " << child;
  out << std::endl;

  return out;
} 

/// Computes the minimum OBB from a set of lower dimensional (< 3d) points
/**
 * \param begin an iterator to type Point3d
 * \param end an iterator to type Point3d
 */
template <class ForwardIterator>
OBB OBB::calc_low_dim_OBB(ForwardIterator begin, ForwardIterator end)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const double INF = std::numeric_limits<double>::max();
  const double TOL = NEAR_ZERO;  // tolerance to expand the OBB

  // get the pose of the points
  boost::shared_ptr<const Ravelin::Pose3d> P = begin->pose;

  // if the beginning is the end, return an empty OBB
  if (begin == end)
  {
    OBB empty;
    empty.center.pose = P;
    return empty;
  }

  // find non-coincident points, if possible to establish a line
  ForwardIterator i = begin;
  Point3d p1 = *i++, p2;
  for (; i != end; i++)
  {
    p2 = *i;
    if ((p1 - p2).norm() > NEAR_ZERO)
      break;
  }

  // check for zero dimensional OBB
  if (i == end)
  {
    OBB o;
    o.R.set_identity();
    o.l[X] = o.l[Y] = o.l[Z] = (double) 0.0;
    o.center = p1;

    #ifndef NDEBUG
    for (ForwardIterator i = begin; i != end; i++)
      assert(!OBB::outside(o, *i, NEAR_ZERO));
    #endif

    return o;
  }

  // check for one dimensional OBB
  for (i++; i != end; i++)
    if (!CompGeom::collinear(p1, p2, *i))
      break;
  if (i == end)
  {
    // project all points onto the line
    Ravelin::Vector3d d1 = Ravelin::Vector3d::normalize(p2 - p1);
    double min_proj = INF, max_proj = -INF;
    for (i = begin; i != end; i++)
    {
      double proj = d1.dot(*i - p1);
      if (proj < min_proj)
        min_proj = proj;
      if (proj > max_proj)
        max_proj = proj;
    }

    // determine "lowest" point on the line
    Point3d lowest = p1 + d1*min_proj;

    // setup the OBB
    Ravelin::Vector3d d2, d3;
    Ravelin::Vector3d::determine_orthonormal_basis(d1, d2, d3);
    OBB o;
    o.R.set_column(X, d1);
    o.R.set_column(Y, d2);
    o.R.set_column(Z, d3);
    o.l[X] = (max_proj - min_proj) * 0.5 + TOL;
    o.l[Y] = o.l[Z] = (double) 0.0;
    o.center = lowest + d1*o.l[X];

    #ifndef NDEBUG
    for (ForwardIterator i = begin; i != end; i++)
      assert(!OBB::outside(o, *i, NEAR_ZERO));
    #endif

    return o;
  }
  
  // ** we have a 2D set of points
  // init the 2D centroid
  Ravelin::Vector2d centroid;

  // use svd to fit points to a plane 
  double plane_offset;
  Ravelin::Vector3d normal;
  CompGeom::fit_plane(begin, end, normal, plane_offset);

  // transform all points to 2D
  boost::shared_ptr<const Ravelin::Pose2d> GLOBAL_2D;
  Ravelin::Matrix3d R = CompGeom::calc_3D_to_2D_matrix(normal);
  double offset = CompGeom::determine_3D_to_2D_offset(Ravelin::Origin3d(*begin), R);
  std::list<Point2d> points_2D;
  for (ForwardIterator i = begin; i != end; i++)
    points_2D.push_back(Point2d(CompGeom::to_2D(*i, R), GLOBAL_2D)); 

  // compute the convex hull of the points
  std::list<Point2d> hull_2D;
  CompGeom::calc_convex_hull(points_2D.begin(), points_2D.end(), std::back_inserter(hull_2D));

  // handle degeneracy
  if (hull_2D.empty())
  {
    // determine seg endpoints
    std::pair<Point2d, Point2d> ep;
    CompGeom::determine_seg_endpoints(points_2D.begin(), points_2D.end(), ep);
    centroid = (ep.first + ep.second) * 0.5;
    
    // project the endpoints and centroid to 3D
    Ravelin::Origin2d oc(centroid), o1(ep.first), o2(ep.second);
    Point3d center(CompGeom::to_3D(oc, R, offset), P);
    Point3d ep1(CompGeom::to_3D(o1, R, offset), P);
    Point3d ep2(CompGeom::to_3D(o2, R, offset), P);

    // see whether we have zero-dimensional or one-dimensional OBB
    Ravelin::Vector3d d1 = ep2 - ep1;
    double d1_len = d1.norm();
    if (d1_len < NEAR_ZERO)
    {
      // zero dimensional OBB
      OBB o;
      o.R.set_identity();
      o.center = center;
      o.l[X] = o.l[Y] = o.l[Z] = (double) 0.0;

      #ifndef NDEBUG
      for (ForwardIterator i = begin; i != end; i++)
        assert(!OBB::outside(o, *i, std::sqrt(NEAR_ZERO)));
      #endif

      return o;
    }

    // one dimensional OBB: setup an orthonormal basis
    Ravelin::Vector3d d2, d3;
    Ravelin::Vector3d::determine_orthonormal_basis(d1, d2, d3);

    // setup the OBB
    OBB o;
    o.R.set_column(X, d1);
    o.R.set_column(Y, d2);
    o.R.set_column(Z, d3);
    o.center = center;
    o.l[X] = (ep.second - centroid).norm() + TOL;
    o.l[Y] = o.l[Z] = (double) 0.0;

    #ifndef NDEBUG
    for (ForwardIterator i = begin; i != end; i++)
      assert(!OBB::outside(o, *i, std::sqrt(NEAR_ZERO)));
    #endif

    return o;
  }

  // determine the minimum bounding rectangle
  Point2d v1, v2, v3, v4;
  CompGeom::calc_min_area_bounding_rect(hull_2D.begin(), hull_2D.end(), v1, v2, v3, v4);

  // project direction to 3D
  Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);
  Ravelin::Vector3d d2(CompGeom::to_3D(Ravelin::Origin2d(v2 - v1), RT), P);
  d2.normalize();

  // get v1 and v3 in 3D
  Point3d v13d(CompGeom::to_3D(Ravelin::Origin2d(v1), RT, offset), P);
  Point3d v33d(CompGeom::to_3D(Ravelin::Origin2d(v3), RT, offset), P);
 
  // start to setup the OBB
  OBB o;
  o.R.set_column(X, normal);
  o.R.set_column(Y, d2);
  o.R.set_column(Z, Ravelin::Vector3d::cross(normal, d2));
  o.center = (v13d + v33d)*(double) 0.5; 
  o.l[X] = TOL;
  o.l[Y] = (v2 - v1).norm()*(double) 0.5 + TOL;
  o.l[Z] = (v3 - v2).norm()*(double) 0.5 + TOL;

  #ifndef NDEBUG
  for (ForwardIterator i = begin; i != end; i++)
    assert(!OBB::outside(o, *i, std::sqrt(NEAR_ZERO)));
  #endif

  return o;
}

/// Computes the minimum OBB from a set of points
/**
 * \param begin an iterator to type Point3d
 * \param end an iterator to type Point3d
 * Algorithm taken from http://www.geometrictools.com - thanks Dave Eberly!
 */
template <class ForwardIterator>
OBB OBB::calc_min_volume_OBB(ForwardIterator begin, ForwardIterator end)
{
  const double INF = std::numeric_limits<double>::max();
  const unsigned X = 0, Y = 1, Z = 2;
  const double TOL = NEAR_ZERO;  // tolerance to expand the OBB

  // get the pose
  boost::shared_ptr<const Ravelin::Pose3d> P = begin->pose;

  // compute the convex hull of the points
  PolyhedronPtr hull = CompGeom::calc_convex_hull(begin, end);
  if (!hull)
    return OBB::calc_low_dim_OBB(begin, end);

  // setup the box with minimum volume
  OBB min_box;
  double min_box_volume = INF;

  // points are three-dimensional; for each facet of the convex hull, 
  // (1) project the hull onto the plane of the facet, (2) project the hull
  // onto the normal line of each facet, (3) calculate the box from the minimum
  // area box in the plane and the height on the line

  // determine coplanar triangles (add all but the first of the coplanar
  // triangles to a list that indicates processing should not occur)
  const IndexedTriArray& hull_mesh = hull->get_mesh();
  std::vector<bool> process(hull_mesh.num_tris(), true);
  std::map<sorted_pair<unsigned>, std::list<unsigned> > ef_map = hull_mesh.determine_edge_facet_map();
  for (std::map<sorted_pair<unsigned>, std::list<unsigned> >::const_iterator i = ef_map.begin(); i != ef_map.end(); i++)
  {
    // if the edge is not coplanar, keep processing through the edges
    if (!hull_mesh.is_coplanar(i->first.first, i->first.second))
      continue;

    // it is coplanar, indicate all but the first triangle in the list should
    // not be processed
    for (std::list<unsigned>::const_iterator j = i->second.begin()++; j != i->second.end(); j++)
      process[*j] = false;
  }

  // process triangles
  const std::vector<Ravelin::Origin3d>& hull_verts = hull_mesh.get_vertices();
  std::vector<Point2d> proj_verts(hull_verts.size());
  for (unsigned i=0; i< hull_mesh.num_tris(); i++)
  {
    // make sure that we are to process this triangle
    if (!process[i])
      continue;

    // determine the plane equation for triangle i
    Triangle t = hull_mesh.get_triangle(i, P);
    Ravelin::Vector3d normal = t.calc_normal();
    Ravelin::Matrix3d R = CompGeom::calc_3D_to_2D_matrix(normal);
    double offset = CompGeom::determine_3D_to_2D_offset(Ravelin::Origin3d(t.a), R);

    // transpose R for projecting back to 3D
    Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);

    // project all points to the plane containing the facet and the line
    // of the normal
    double min_height = INF, max_height = -INF;
    for (unsigned j=0; j< hull_verts.size(); j++)
    {
      // get the hull vertex as a Point3d
      Point3d hv(hull_verts[j], P);

      // project the vertex onto the plane containing the triangle
      proj_verts[j] = CompGeom::to_2D(hv, R);

      // compute the projection onto the normal
      double dot = hv.dot(normal) - offset;
      if (dot < min_height)
        min_height = dot;
      if (dot > max_height)
        max_height = dot;
    }

    // calculate the minimum area bounding rectangle
    Point2d a, b, c, d;
    CompGeom::calc_min_area_bounding_rect(proj_verts.begin(), proj_verts.end(), a, b, c, d);

    // special case: the minimum area bounding rectangle is just a line
    if ((b-a).norm()*(c-b).norm() < NEAR_ZERO*NEAR_ZERO)
    {
      // special case: the minimum area bounding rectangle is just a point
      if ((b-a).norm() < NEAR_ZERO && (c-b).norm() < NEAR_ZERO && 
          (d-c).norm() < NEAR_ZERO)
      {
        // this will be the minimum volume box
        OBB box;
        box.R = IDENTITY_3x3;
        Ravelin::Origin2d oa(a); 
        box.center = Point3d(CompGeom::to_3D(oa, RT, offset), P); 
        box.l[X] = box.l[Y] = box.l[Z] = (double) 0.0; 
        return box;
      }

      // determine the line endpoints
      Point2d points[4] = { a, b, c, d };
      std::pair<Point2d, Point2d> ep;
      CompGeom::determine_seg_endpoints(points, points+4, ep);
      Ravelin::Origin2d ep1(ep.first), ep2(ep.second);

      // compute the centroid of the line
      Ravelin::Origin2d centroid((ep.first + ep.second)*(double) 0.5);

      // project the centroid to 3D
      Point3d centroid_3d(CompGeom::to_3D(centroid, RT, offset), P);

      // determine the direction of the line
      Ravelin::Vector3d dir(CompGeom::to_3D(ep1, RT, offset) -
                            CompGeom::to_3D(ep2, RT, offset), P);
      dir.normalize();

      // setup the bounding box -- it will be zero volume
      OBB box;
      box.R.set_column(X, normal);
      box.R.set_column(Y, dir);
      box.R.set_column(Z, Ravelin::Vector3d::cross(normal, dir));
      box.center = centroid_3d + normal*(min_height + max_height)*(double) 0.5;
      box.l[0] = (max_height - min_height) * (double) 0.5 + TOL; 
      box.l[1] = (ep.second - ep.first).norm() * (double) 0.5 + TOL;
      box.l[2] = (double) 0.0 + TOL;
      return box;
    }

    // compute the centroid of the bounding rectangle
    Ravelin::Origin2d centroid_2D((a + c) * (double) 0.5);

    // project the centroid to 3D
    Point3d centroid(CompGeom::to_3D(centroid_2D, RT, offset), P);

    // determine the bounding rectangle edges in 3D
    Ravelin::Vector3d d2(CompGeom::to_3D(Ravelin::Origin2d(b - a), RT), P);
    d2.normalize();

    // setup the OBB and calculate its volume
    OBB box;
    box.R.set_column(X, normal);
    box.R.set_column(Y, d2);
    box.R.set_column(Z, Ravelin::Vector3d::cross(normal, d2));
    assert(box.R.is_orthonormal());
    box.center = centroid + normal*(min_height + max_height)*0.5;
    box.l[X] = (max_height - min_height)*0.5 + TOL;
    box.l[Y] = (b - a).norm() * 0.5 + TOL;
    box.l[Z] = (c - b).norm() * 0.5 + TOL;
    #ifndef NDEBUG
    for (ForwardIterator i = begin; i != end; i++)
      assert(!OBB::outside(min_box, *i, NEAR_ZERO));
    #endif

    // if the volume is less than the minimum, keep the box
    double box_volume = box.calc_volume();
    if (box_volume < min_box_volume)
    {
      min_box_volume = box_volume;
      min_box = box;
    }
  }

  // the minimum volume box can also be supported by three mutually orthogonal
  // edges of the convex hull.  For each triple of orthogonal edges, compute
  // the minimum volume box for that coordinate frame by projecting the points
  // onto the axes of the frame
  std::vector<std::list<unsigned> > ve_map = hull_mesh.determine_vertex_edge_map();
  for (unsigned i=0; i< ve_map.size(); i++)
  {
    // get the list of incident vertices as a vector
    std::vector<unsigned> ivs(ve_map[i].begin(), ve_map[i].end());

    // if- for some reason- there are fewer than three edges in the vector,
    // move onto next vertex
    if (ivs.size() < 3)
      continue;

    // process each combination of three vertices
    for (unsigned j=0; j< ivs.size(); j++)
    {
      // determine edge 1
      Ravelin::Vector3d d1(hull_verts[ivs[j]] - hull_verts[i], P);
      double d1_len = d1.norm();
      if (d1_len < NEAR_ZERO)
        continue;
      d1 /= d1_len;

      for (unsigned k=j+1; k< ivs.size(); k++)
      {
        // determine edge 2
        Ravelin::Vector3d d2(hull_verts[ivs[k]] - hull_verts[i], P);
        double d2_len = d2.norm();
        if (d2_len < NEAR_ZERO)
          continue;
        d2 /= d2_len;
        
        // if e1 and e2 are not orthogonal, keep processing
        if (std::fabs(d1.dot(d2)) > NEAR_ZERO)
          continue;

        for (unsigned m=k+1; m< ivs.size(); m++)
        {
          // determine edge 3
          Ravelin::Vector3d d3(hull_verts[ivs[m]] - hull_verts[i], P);
          double d3_len = d3.norm();
          if (d3_len < NEAR_ZERO)
            continue;
          d3 /= d3_len;

          // if e1 and e3 are not orthogonal or e2 and e3 are not orthogonal,
          // keep processing edges
          if (std::fabs(d1.dot(d3)) > NEAR_ZERO || 
              std::fabs(d2.dot(d3)) > NEAR_ZERO)
            continue;

          // d1, d2, d3 are orthogonal; project all hull points onto these axes
          // to determine the box
          double d1_min_proj = INF, d1_max_proj = -INF;
          double d2_min_proj = INF, d2_max_proj = -INF;
          double d3_min_proj = INF, d3_max_proj = -INF;
          for (unsigned n=0; n< hull_verts.size(); n++)
          {
            // get point in frame
            Point3d pt(hull_verts[n] - hull_verts[i], P);

            // determine projections
            double d1_dot = d1.dot(pt);
            double d2_dot = d2.dot(pt);
            double d3_dot = d3.dot(pt);
            if (d1_dot < d1_min_proj) d1_min_proj = d1_dot;
            if (d1_dot > d1_max_proj) d1_max_proj = d1_dot;
            if (d2_dot < d2_min_proj) d2_min_proj = d2_dot;
            if (d2_dot > d2_max_proj) d2_max_proj = d2_dot;
            if (d3_dot < d3_min_proj) d3_min_proj = d3_dot;
            if (d3_dot > d3_max_proj) d3_max_proj = d3_dot;
          }

          // calculate the volume of the box
          double volume = (d1_max_proj - d1_min_proj) * 
                        (d2_max_proj - d2_min_proj) *
                        (d3_max_proj - d3_min_proj);

          // only if the volume is smaller than the smallest OBB volume do we
          // setup the box
          if (volume >= min_box_volume)
            continue;

          // store the new minimum box volume
          min_box_volume = volume;

          // determine minimum corner of the box
          Point3d corner = Point3d(hull_verts[i], P) + d1*d1_min_proj + 
                            d2*d2_min_proj + d3*d3_min_proj;

          // setup the new minimum box
          min_box.R.set_column(X, d1);
          min_box.R.set_column(Y, d2);
          min_box.R.set_column(Z, d3);
          min_box.l[X] = (d1_max_proj - d1_min_proj)*0.5;
          min_box.l[Y] = (d2_max_proj - d2_min_proj)*0.5;
          min_box.l[Z] = (d3_max_proj - d3_min_proj)*0.5;
          min_box.center = corner + d1*min_box.l[X] + d2*min_box.l[Y] +
                                    d3*min_box.l[Z];
          min_box.l[X] += TOL;
          min_box.l[Y] += TOL;
          min_box.l[Z] += TOL;
          
          #ifndef NDEBUG
          for (ForwardIterator i = begin; i != end; i++)
            assert(!OBB::outside(min_box, *i, NEAR_ZERO));
          #endif
        }
      }
    }
  }

  return min_box;
}

/// Computes an OBB from a set of points
/**
 * \param begin an iterator to type Point3d
 * \param end an iterator to type Point3d
 * Algorithm taken from [Ericson, 2005]
 */
template <class ForwardIterator>
OBB::OBB(ForwardIterator begin, ForwardIterator end)
{
  const unsigned X = 0, Y = 1, Z = 2, THREE_D = 3;
  Ravelin::Vector3d normal;
  std::list<Ravelin::Vector3d> test;
  std::set<Ravelin::Vector3d> tested;  
  boost::shared_ptr<const Ravelin::Pose3d> P;

  // setup the pose
  if (begin != end)
    P = begin->pose;

  // initialize the center to zero
  this->center.set_zero();
  this->center.pose = P;

  // compute the convex hull of the points
  PolyhedronPtr hull;
  try
  {
    hull = CompGeom::calc_convex_hull(begin, end);
  }
  catch (NumericalException e)
  {
  }
  if (!hull)
  {
    *this = OBB::calc_low_dim_OBB(begin, end);
    return;
  }  

  // determine the area and centroid of all triangles
  const std::vector<Ravelin::Origin3d>& verts = hull->get_vertices();
  unsigned n = hull->get_facets().size();
  double total_area = 0;
  std::vector<double> areas(n);
  std::vector<Point3d> centroids(n);
  for (unsigned i=0; i< n; i++)
  {
    const IndexedTri& itri = hull->get_facets()[i];
    Triangle tri(Point3d(verts[itri.a], P), 
                 Point3d(verts[itri.b], P), 
                 Point3d(verts[itri.c], P)); 
    areas[i] = tri.calc_area();
    centroids[i] = (tri.a + tri.b + tri.c);
    total_area += areas[i];
    this->center = this->center + centroids[i]*areas[i];
  }
  this->center /= (total_area*3.0);

  // compute the covariance matrix of the points
  // 1st: subtract the covariance components of the centroid
  Ravelin::Matrix3d C;
  for (unsigned i=0; i< THREE_D; i++)
    for (unsigned j=i; j< THREE_D; j++)
      C(i,j) = -this->center[i]*this->center[j];

  // 2nd: add in the area based components
  for (unsigned i=0; i< THREE_D; i++)
    for (unsigned j=i; j< THREE_D; j++)
    {
      for (unsigned k=0; k< n; k++)
      {
        const IndexedTri& itri = hull->get_facets()[k];
        Triangle tri(Point3d(verts[itri.a], P), Point3d(verts[itri.b], P), Point3d(verts[itri.c], P)); 
        const Point3d& a = tri.a;
        const Point3d& b = tri.b;
        const Point3d& c = tri.c;
        C(i,j) += areas[k]/12.0 * (centroids[k][i]*centroids[k][j] + a[i]*a[j] + b[i]*b[j] + c[i]*c[j]);
      }
      C(i,j) *= 1.0/total_area;
    }

  // 3rd: make covariance matrix symmetric
  for (unsigned i=0; i< THREE_D; i++)
    for (unsigned j=i+1; j< THREE_D; j++)
      C(j,i) = C(i,j);

  // determine the eigenvalues and eigenvectors of the covariance matrix
  Ravelin::Vector3d evals;
  SAFESTATIC FastThreadable<Ravelin::LinAlgd> LA;
  LA().eig_symm_plus(C, evals);
  
  // first eigenvector will be direction of minimum variance
  // but add all three eigenvectors
  for (unsigned i=0; i< 3; i++)
  {
    Ravelin::Vector3d col(this->center.pose);
    C.get_column(i, col);
    double nrm = col.norm();
    if (nrm < NEAR_ZERO)
      continue;
    col /= nrm;
    test.push_back(col);
  }

  // setup the minimum volume
  double min_vol = std::numeric_limits<double>::max();

  // keep going until test is empty
  assert(!test.empty());
  while (!test.empty())
  {
    // get the normal direction
    Ravelin::Vector3d normal = test.front();
    test.pop_front();

    // if this direction has already been tested, do not test it again
    if (tested.find(normal) != tested.end())
      continue;

    // align OBB with minimum bounding rectangle using the normal
    Ravelin::Vector3d d2(this->center.pose), d3(this->center.pose);
    align(begin, end, normal, d2);
    d3 = Ravelin::Vector3d::normalize(Ravelin::Vector3d::cross(normal, d2));

    // compute the lengths and the volume
    double lengths[3];
    calc_lengths(normal, d2, d3, this->center, begin, end, lengths);
    double vol = lengths[0]*lengths[1]*lengths[2];

    // if the new volume is the minimum, add the other directions for testing
    if (vol < min_vol)
    {
      // setup the OBB
      this->R.set_column(X, normal);
      this->R.set_column(Y, d2);
      this->R.set_column(Z, d3);
      this->l = Ravelin::Vector3d(lengths[0], lengths[1], lengths[2]);

      // store the minimum volume
      min_vol = vol;

      // add other directions for testing
      test.push_back(d2);
      test.push_back(d3);
    }

    // indicate that this direction has been tested
    tested.insert(normal);
    tested.insert(-normal);
  }
  assert(this->R.is_orthonormal());
} 

/// Expands this OBB (if necessary) to fit the given points
/**
 * \param begin an iterator to the beginning of a container of type Ravelin:Point3d:
 * \param end an iterator to the end of a container of type Point3d
 */
template <class ForwardIterator>
void OBB::expand_to_fit(ForwardIterator begin, ForwardIterator end)
{
  const unsigned THREE_D = 3; 

  // get the corners of the OBB as if it were an AABB
  Point3d lo = -this->l, hi = this->l;

  // process all points
  for (ForwardIterator i=begin; i != end; i++)
  {
    Point3d pt = R.transpose_mult(*i - this->center);
    for (unsigned i=0; i< THREE_D; i++)
      if (pt[i] < lo[i])
        lo[i] = pt[i];
      else if (pt[i] > hi[i])
        hi[i] = pt[i];
  }

  // half lo and hi
  lo *= 0.5;
  hi *= 0.5;

  // transform the center of the box
  this->center += R * Ravelin::Origin3d(lo+hi);

  // set the new lengths of the box
  this->l = hi - lo;
}

/// Aligns this OBB with the minimum area bounding rectangle projected along the first dimension of the OBB
template <class ForwardIterator>
void OBB::align(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& d1, Ravelin::Vector3d& d2)
{
  // project all points to the plane perpendicular to the first direction
  boost::shared_ptr<const Ravelin::Pose2d> GLOBAL_2D;
  Ravelin::Matrix3d R2d = CompGeom::calc_3D_to_2D_matrix(d1);
  std::list<Point2d> points_2D;
  for (ForwardIterator i = begin; i != end; i++)
    points_2D.push_back(Point2d(CompGeom::to_2D(*i, R2d), GLOBAL_2D)); 

  // determine the minimum bounding rectangle of the projected points
  Point2d p1, p2, p3, p4;
  CompGeom::calc_min_area_bounding_rect(points_2D.begin(), points_2D.end(), p1, p2, p3, p4);

  // project the vectors of the minimum bounding rectangle back to 3D
  R2d.transpose();
  Ravelin::Vector2d d2_2D = p2 - p1;
  if (d2_2D.norm() < NEAR_ZERO)
  {
    d2_2D = p3 - p2;
    if (d2_2D.norm() < NEAR_ZERO)
      d2_2D = Ravelin::Vector2d(1,0);
  }
  d2 = CompGeom::to_3D(Ravelin::Origin2d(d2_2D), R2d);
  d2.normalize();
}

/// Gets the eight vertices of the bounding box
/**
 * \note vertices are output in no particular order
 */
template <class OutputIterator>
OutputIterator OBB::get_vertices(OutputIterator begin) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the three axes of the OBB, scaled by axis lengths
  Ravelin::Vector3d axis1(center.pose), axis2(center.pose), axis3(center.pose);
  R.get_column(X, axis1); 
  R.get_column(Y, axis2); 
  R.get_column(Z, axis3); 
  axis1 *= l[X];
  axis2 *= l[Y];
  axis3 *= l[Z];

  // code below works
  /*  
  *begin++ = center - axis1 - axis2 - axis3;
  *begin++ = center - axis1 - axis2 + axis3;
  *begin++ = center - axis1 + axis2 - axis3;
  *begin++ = center - axis1 + axis2 + axis3;
  *begin++ = center + axis1 - axis2 - axis3;
  *begin++ = center + axis1 - axis2 + axis3;
  *begin++ = center + axis1 + axis2 - axis3;
  *begin++ = center + axis1 + axis2 + axis3;
  */

  // this should be 50% faster or so 
  Ravelin::Vector3d axis1_2 = axis1 + axis1;
  Ravelin::Vector3d axis2_2 = axis2 + axis2;
  Ravelin::Vector3d axis3_2 = axis3 + axis3;
  Ravelin::Vector3d corner = center - axis1 - axis2 - axis3;  // ---
  *begin++ = corner;
  corner += axis3_2;                                // --+
  *begin++ = corner;
  corner += axis2_2;                                // -++
  *begin++ = corner;
  corner -= axis3_2;                                // -+-
  *begin++ = corner;
  corner += axis1_2;                                // ++-
  *begin++ = corner;
  corner += axis3_2;                                // +++
  *begin++ = corner;
  corner -= axis2_2;                                // +-+ 
  *begin++ = corner;
  corner -= axis3_2;                                // +--
  *begin++ = corner; 

  return begin;
}

template <class ForwardIterator>
void OBB::calc_lengths(const Ravelin::Vector3d& d1, const Ravelin::Vector3d& d2, const Ravelin::Vector3d& d3, const Point3d& center, ForwardIterator begin, ForwardIterator end, double lengths[3])
{
  // compute the lengths
  lengths[0] = 0.0;
  lengths[1] = 0.0;
  lengths[2] = 0.0;

  for (; begin != end; begin++)
  {
    Ravelin::Vector3d v = *begin - center;
    double l0 = d1.dot(v);
    double l1 = d2.dot(v);
    double l2 = d3.dot(v);
    lengths[0] = std::max(lengths[0], std::fabs(l0));
    lengths[1] = std::max(lengths[1], std::fabs(l1));
    lengths[2] = std::max(lengths[2], std::fabs(l2));
  }

} // end method

