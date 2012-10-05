/****************************************************************************
 * Copyright 2008 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

// For outputting description of OBB (primarily for debugging purposes)
inline std::ostream& operator<<(std::ostream& out, const OBB& o)
{
  out << " (address): " << &o << std::endl;
  out << " center: " << o.center << std::endl;
  out << " axes: " << std::endl << o.R;
  out << "   (as axis-angle): " << AAngle(&o.R) << std::endl;
  out << " half lengths: " << o.l << "  1/8 volume: " << o.calc_volume() << std::endl;
  out << " children:";
  BOOST_FOREACH(BVPtr child, o.children)
    out << " " << child;
  out << std::endl;

  return out;
} 

/// Computes the minimum OBB from a set of lower dimensional (< 3d) points
/**
 * \param begin an iterator to type Vector3
 * \param end an iterator to type Vector3
 */
template <class ForwardIterator>
OBB OBB::calc_low_dim_OBB(ForwardIterator begin, ForwardIterator end)
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Real INF = std::numeric_limits<Real>::max();
  const Real TOL = NEAR_ZERO;  // tolerance to expand the OBB

  // find non-coincident points, if possible to establish a line
  ForwardIterator i = begin;
  Vector3 p1 = *i++, p2;
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
    o.R = IDENTITY_3x3;
    o.l[X] = o.l[Y] = o.l[Z] = (Real) 0.0;
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
    Vector3 d1 = Vector3::normalize(p2 - p1);
    Real min_proj = INF, max_proj = -INF;
    for (i = begin; i != end; i++)
    {
      Real proj = Vector3::dot(d1, *i - p1);
      if (proj < min_proj)
        min_proj = proj;
      if (proj > max_proj)
        max_proj = proj;
    }

    // determine "lowest" point on the line
    Vector3 lowest = p1 + d1*min_proj;

    // setup the OBB
    Vector3 d2, d3;
    Vector3::determine_orthonormal_basis(d1, d2, d3);
    OBB o;
    o.R.set_column(X, d1);
    o.R.set_column(Y, d2);
    o.R.set_column(Z, d3);
    o.l[X] = (max_proj - min_proj) * 0.5 + TOL;
    o.l[Y] = o.l[Z] = (Real) 0.0;
    o.center = lowest + d1*o.l[X];

    #ifndef NDEBUG
    for (ForwardIterator i = begin; i != end; i++)
      assert(!OBB::outside(o, *i, NEAR_ZERO));
    #endif

    return o;
  }
  
  // ** we have a 2D set of points
  // init the 2D centroid
  Vector2 centroid;

  // use svd to fit points to a plane 
  Real plane_offset;
  Vector3 normal;
  CompGeom::fit_plane(begin, end, normal, plane_offset);

  // transform all points to 2D
  Matrix3 R = CompGeom::calc_3D_to_2D_matrix(normal);
  Real offset = CompGeom::determine_3D_to_2D_offset(*begin, R);
  std::list<Vector2> points_2D;
  CompGeom::to_2D(begin, end, std::back_inserter(points_2D), R);
    
  // compute the convex hull of the points
  std::list<Vector2> hull_2D;
  CompGeom::calc_convex_hull_2D(points_2D.begin(), points_2D.end(), std::back_inserter(hull_2D));

  // handle degeneracy
  if (hull_2D.empty())
  {
    // determine seg endpoints
    std::pair<Vector2, Vector2> ep;
    CompGeom::determine_seg_endpoints(points_2D.begin(), points_2D.end(), ep);
    centroid = (ep.first + ep.second) * 0.5;
    
    // project the endpoints and centroid to 3D
    Vector3 center = CompGeom::to_3D(centroid, R, offset);
    Vector3 ep1 = CompGeom::to_3D(ep.first, R, offset);
    Vector3 ep2 = CompGeom::to_3D(ep.second, R, offset);

    // see whether we have zero-dimensional or one-dimensional OBB
    Vector3 d1 = ep2 - ep1;
    Real d1_len = d1.norm();
    if (d1_len < NEAR_ZERO)
    {
      // zero dimensional OBB
      OBB o;
      o.R = IDENTITY_3x3;
      o.center = center;
      o.l[X] = o.l[Y] = o.l[Z] = (Real) 0.0;

      #ifndef NDEBUG
      for (ForwardIterator i = begin; i != end; i++)
        assert(!OBB::outside(o, *i, std::sqrt(NEAR_ZERO)));
      #endif

      return o;
    }

    // one dimensional OBB: setup an orthonormal basis
    Vector3 d2, d3;
    Vector3::determine_orthonormal_basis(d1, d2, d3);

    // setup the OBB
    OBB o;
    o.R.set_column(X, d1);
    o.R.set_column(Y, d2);
    o.R.set_column(Z, d3);
    o.center = center;
    o.l[X] = (ep.second - centroid).norm() + TOL;
    o.l[Y] = o.l[Z] = (Real) 0.0;

    #ifndef NDEBUG
    for (ForwardIterator i = begin; i != end; i++)
      assert(!OBB::outside(o, *i, std::sqrt(NEAR_ZERO)));
    #endif

    return o;
  }

  // determine the minimum bounding rectangle
  Vector2 v1, v2, v3, v4;
  CompGeom::calc_min_area_bounding_rect(hull_2D.begin(), hull_2D.end(), v1, v2, v3, v4);

  // project direction to 3D
  Matrix3 RT = Matrix3::transpose(R);
  Vector3 d2 = Vector3::normalize(CompGeom::to_3D(v2 - v1, RT));

  // get v1 and v3 in 3D
  Vector3 v13d = CompGeom::to_3D(v1, RT, offset);
  Vector3 v33d = CompGeom::to_3D(v3, RT, offset);
 
  // start to setup the OBB
  OBB o;
  o.R.set_column(X, normal);
  o.R.set_column(Y, d2);
  o.R.set_column(Z, Vector3::cross(normal, d2));
  o.center = (v13d + v33d)*(Real) 0.5; 
  o.l[X] = TOL;
  o.l[Y] = (v2 - v1).norm()*(Real) 0.5 + TOL;
  o.l[Z] = (v3 - v2).norm()*(Real) 0.5 + TOL;

  #ifndef NDEBUG
  for (ForwardIterator i = begin; i != end; i++)
    assert(!OBB::outside(o, *i, std::sqrt(NEAR_ZERO)));
  #endif

  return o;
}

/// Computes the minimum OBB from a set of points
/**
 * \param begin an iterator to type Vector3
 * \param end an iterator to type Vector3
 * Algorithm taken from http://www.geometrictools.com - thanks Dave Eberly!
 */
template <class ForwardIterator>
OBB OBB::calc_min_volume_OBB(ForwardIterator begin, ForwardIterator end)
{
  const Real INF = std::numeric_limits<Real>::max();
  const unsigned X = 0, Y = 1, Z = 2;
  const Real TOL = NEAR_ZERO;  // tolerance to expand the OBB

  // compute the convex hull of the points
  PolyhedronPtr hull = CompGeom::calc_convex_hull_3D(begin, end);
  bool is_3D = hull;
  if (!is_3D)
    return OBB::calc_low_dim_OBB(begin, end);

  // setup the box with minimum volume
  OBB min_box;
  Real min_box_volume = INF;

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
  const std::vector<Vector3>& hull_verts = hull_mesh.get_vertices();
  std::vector<Vector2> proj_verts(hull_verts.size());
  for (unsigned i=0; i< hull_mesh.num_tris(); i++)
  {
    // make sure that we are to process this triangle
    if (!process[i])
      continue;

    // determine the plane equation for triangle i
    Triangle t = hull_mesh.get_triangle(i);
    Vector3 normal = t.calc_normal();
    Matrix3 R = CompGeom::calc_3D_to_2D_matrix(normal);
    Real offset = CompGeom::determine_3D_to_2D_offset(t.a, R);

    // transpose R for projecting back to 3D
    Matrix3 RT = Matrix3::transpose(R);

    // project all points to the plane containing the facet and the line
    // of the normal
    Real min_height = INF, max_height = -INF;
    for (unsigned j=0; j< hull_verts.size(); j++)
    {
      // project the vertex onto the plane containing the triangle
      proj_verts[j] = CompGeom::to_2D(hull_verts[j], R);

      // compute the projection onto the normal
      Real dot = Vector3::dot(hull_verts[j], normal) - offset;
      if (dot < min_height)
        min_height = dot;
      if (dot > max_height)
        max_height = dot;
    }

    // calculate the minimum area bounding rectangle
    Vector2 a, b, c, d;
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
        box.center = CompGeom::to_3D(a, RT, offset); 
        box.l[X] = box.l[Y] = box.l[Z] = (Real) 0.0; 
        return box;
      }

      // determine the line endpoints
      Vector2 points[4] = { a, b, c, d };
      std::pair<Vector2, Vector2> ep;
      CompGeom::determine_seg_endpoints(points, points+4, ep);

      // compute the centroid of the line
      Vector2 centroid = (ep.first + ep.second)*(Real) 0.5;

      // project the centroid to 3D
      Vector3 centroid_3d = CompGeom::to_3D(centroid, RT, offset);

      // determine the direction of the line
      Vector3 dir = Vector3::normalize(CompGeom::to_3D(ep.second, RT, offset) -
                                       CompGeom::to_3D(ep.first, RT, offset));

      // setup the bounding box -- it will be zero volume
      OBB box;
      box.R.set_column(X, normal);
      box.R.set_column(Y, dir);
      box.R.set_column(Z, Vector3::cross(normal, dir));
      box.center = centroid_3d + normal*(min_height + max_height)*(Real) 0.5;
      box.l[0] = (max_height - min_height) * (Real) 0.5 + TOL; 
      box.l[1] = (ep.second - ep.first).norm() * (Real) 0.5 + TOL;
      box.l[2] = (Real) 0.0 + TOL;
      return box;
    }

    // compute the centroid of the bounding rectangle
    Vector2 centroid_2D = (a + c) * (Real) 0.5;

    // project the centroid to 3D
    Vector3 centroid = CompGeom::to_3D(centroid_2D, RT, offset);

    // determine the bounding rectangle edges in 3D
    Vector3 d2 = Vector3::normalize(CompGeom::to_3D(b - a, RT));

    // setup the OBB and calculate its volume
    OBB box;
    box.R.set_column(X, normal);
    box.R.set_column(Y, d2);
    box.R.set_column(Z, Vector3::cross(normal, d2));
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
    Real box_volume = box.calc_volume();
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
      Vector3 d1 = hull_verts[ivs[j]] - hull_verts[i];
      Real d1_len = d1.norm();
      if (d1_len < NEAR_ZERO)
        continue;
      d1 /= d1_len;

      for (unsigned k=j+1; k< ivs.size(); k++)
      {
        // determine edge 2
        Vector3 d2 = hull_verts[ivs[k]] - hull_verts[i];
        Real d2_len = d2.norm();
        if (d2_len < NEAR_ZERO)
          continue;
        d2 /= d2_len;
        
        // if e1 and e2 are not orthogonal, keep processing
        if (std::fabs(d1.dot(d2)) > NEAR_ZERO)
          continue;

        for (unsigned m=k+1; m< ivs.size(); m++)
        {
          // determine edge 3
          Vector3 d3 = hull_verts[ivs[m]] - hull_verts[i];
          Real d3_len = d3.norm();
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
          Real d1_min_proj = INF, d1_max_proj = -INF;
          Real d2_min_proj = INF, d2_max_proj = -INF;
          Real d3_min_proj = INF, d3_max_proj = -INF;
          for (unsigned n=0; n< hull_verts.size(); n++)
          {
            // get point in frame
            Vector3 pt = hull_verts[n] - hull_verts[i];

            // determine projections
            Real d1_dot = d1.dot(pt);
            Real d2_dot = d2.dot(pt);
            Real d3_dot = d3.dot(pt);
            if (d1_dot < d1_min_proj) d1_min_proj = d1_dot;
            if (d1_dot > d1_max_proj) d1_max_proj = d1_dot;
            if (d2_dot < d2_min_proj) d2_min_proj = d2_dot;
            if (d2_dot > d2_max_proj) d2_max_proj = d2_dot;
            if (d3_dot < d3_min_proj) d3_min_proj = d3_dot;
            if (d3_dot > d3_max_proj) d3_max_proj = d3_dot;
          }

          // calculate the volume of the box
          Real volume = (d1_max_proj - d1_min_proj) * 
                        (d2_max_proj - d2_min_proj) *
                        (d3_max_proj - d3_min_proj);

          // only if the volume is smaller than the smallest OBB volume do we
          // setup the box
          if (volume >= min_box_volume)
            continue;

          // store the new minimum box volume
          min_box_volume = volume;

          // determine minimum corner of the box
          Vector3 corner = hull_verts[i] + d1*d1_min_proj + d2*d2_min_proj + 
                           d3*d3_min_proj;

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
 * \param begin an iterator to type Vector3
 * \param end an iterator to type Vector3
 * Algorithm taken from [Ericson, 2005]
 */
template <class ForwardIterator>
OBB::OBB(ForwardIterator begin, ForwardIterator end)
{
  const unsigned X = 0, Y = 1, Z = 2, THREE_D = 3;
  Vector3 normal;
  std::list<Vector3> test;
  std::set<Vector3> tested;  

  // initialize the center to zero
  this->center = ZEROS_3;

  // compute the convex hull of the points
  PolyhedronPtr hull = CompGeom::calc_convex_hull_3D(begin, end);
  bool is_3D = hull;
  if (!is_3D)
  {
    *this = OBB::calc_low_dim_OBB(begin, end);
    return;
  }  

  // determine the area and centroid of all triangles
  const std::vector<Vector3>& verts = hull->get_vertices();
  unsigned n = hull->get_facets().size();
  Real total_area = 0;
  std::vector<Real> areas(n);
  std::vector<Vector3> centroids(n);
  for (unsigned i=0; i< n; i++)
  {
    const IndexedTri& itri = hull->get_facets()[i];
    Triangle tri(verts[itri.a], verts[itri.b], verts[itri.c]); 
    areas[i] = tri.calc_area();
    centroids[i] = (tri.a + tri.b + tri.c);
    total_area += areas[i];
    this->center += centroids[i]*areas[i];
  }
  this->center /= (total_area*3.0);

  // compute the covariance matrix of the points
  // 1st: subtract the covariance components of the centroid
  SAFESTATIC FastThreadable<MatrixN> Cx;
  MatrixN& C = Cx();
  C.resize(THREE_D, THREE_D);
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
        Triangle tri(verts[itri.a], verts[itri.b], verts[itri.c]); 
        const Vector3& a = tri.a;
        const Vector3& b = tri.b;
        const Vector3& c = tri.c;
        C(i,j) += areas[k]/12.0 * (centroids[k][i]*centroids[k][j] + a[i]*a[j] + b[i]*b[j] + c[i]*c[j]);
      }
      C(i,j) *= 1.0/total_area;
    }

  // 3rd: make covariance matrix symmetric
  for (unsigned i=0; i< THREE_D; i++)
    for (unsigned j=i+1; j< THREE_D; j++)
      C(j,i) = C(i,j);

  // determine the eigenvalues and eigenvectors of the covariance matrix
  SAFESTATIC FastThreadable<VectorN> evalsx;
  VectorN& evals = evalsx();
  LinAlg::eig_symm_plus(C, evals);
  
  // first eigenvector will be direction of minimum variance
  // but add all three eigenvectors
  for (unsigned i=0; i< 3; i++)
  {
    Vector3 col;
    C.get_column(i, col.begin());
    Real nrm = col.norm();
    if (nrm < NEAR_ZERO)
      continue;
    col /= nrm;
    test.push_back(col);
  }

  // setup the minimum volume
  Real min_vol = std::numeric_limits<Real>::max();

  // keep going until test is empty
  assert(!test.empty());
  while (!test.empty())
  {
    // get the normal direction
    Vector3 normal = test.front();
    test.pop_front();

    // if this direction has already been tested, do not test it again
    if (tested.find(normal) != tested.end())
      continue;

    // align OBB with minimum bounding rectangle using the normal
    Vector3 d2, d3;
    align(begin, end, normal, d2);
    d3 = Vector3::normalize(Vector3::cross(normal, d2));

    // compute the lengths and the volume
    Real lengths[3];
    calc_lengths(normal, d2, d3, this->center, begin, end, lengths);
    Real vol = lengths[0]*lengths[1]*lengths[2];

    // if the new volume is the minimum, add the other directions for testing
    if (vol < min_vol)
    {
      // setup the OBB
      this->R.set_column(X, normal);
      this->R.set_column(Y, d2);
      this->R.set_column(Z, d3);
      this->l = Vector3(lengths[0], lengths[1], lengths[2]);

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
 * \param begin an iterator to the beginning of a container of type Vector3
 * \param end an iterator to the end of a container of type Vector3
 */
template <class ForwardIterator>
void OBB::expand_to_fit(ForwardIterator begin, ForwardIterator end)
{
  const unsigned THREE_D = 3; 

  // get the corners of the OBB as if it were an AABB
  Vector3 lo = -this->l, hi = this->l;

  // process all points
  for (ForwardIterator i=begin; i != end; i++)
  {
    Vector3 pt = R.transpose_mult(*i - this->center);
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
  this->center += R * (lo+hi);

  // set the new lengths of the box
  this->l = hi - lo;
}

/// Aligns this OBB with the minimum area bounding rectangle projected along the first dimension of the OBB
template <class ForwardIterator>
void OBB::align(ForwardIterator begin, ForwardIterator end, const Vector3& d1, Vector3& d2)
{
  // project all points to the plane perpendicular to the first direction
  Matrix3 R2d = CompGeom::calc_3D_to_2D_matrix(d1);
  std::list<Vector2> points_2D;
  CompGeom::to_2D(begin, end, std::back_inserter(points_2D), R2d);

  // determine the minimum bounding rectangle of the projected points
  Vector2 p1, p2, p3, p4;
  CompGeom::calc_min_area_bounding_rect(points_2D.begin(), points_2D.end(), p1, p2, p3, p4);

  // project the vectors of the minimum bounding rectangle back to 3D
  R2d.transpose();
  Vector2 d2_2D = p2 - p1;
  if (d2_2D.norm() < NEAR_ZERO)
  {
    d2_2D = p3 - p2;
    if (d2_2D.norm() < NEAR_ZERO)
      d2_2D = Vector2(1,0);
  }
  d2 = CompGeom::to_3D(d2_2D, R2d);
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
  Vector3 axis1, axis2, axis3;
  R.get_column(X, axis1.begin()); 
  R.get_column(Y, axis2.begin()); 
  R.get_column(Z, axis3.begin()); 
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
  Vector3 axis1_2 = axis1 + axis1;
  Vector3 axis2_2 = axis2 + axis2;
  Vector3 axis3_2 = axis3 + axis3;
  Vector3 corner = center - axis1 - axis2 - axis3;  // ---
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
void OBB::calc_lengths(const Vector3& d1, const Vector3& d2, const Vector3& d3, const Vector3& center, ForwardIterator begin, ForwardIterator end, Real lengths[3])
{
  // compute the lengths
  lengths[0] = 0.0;
  lengths[1] = 0.0;
  lengths[2] = 0.0;

  for (; begin != end; begin++)
  {
    Vector3 v = *begin - center;
    Real l0 = d1.dot(v);
    Real l1 = d2.dot(v);
    Real l2 = d3.dot(v);
    lengths[0] = std::max(lengths[0], std::fabs(l0));
    lengths[1] = std::max(lengths[1], std::fabs(l1));
    lengths[2] = std::max(lengths[2], std::fabs(l2));
  }

} // end method

