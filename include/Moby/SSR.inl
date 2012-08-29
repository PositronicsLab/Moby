/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Computes a SSR from a set of points
/**
 * \param begin an iterator to type Vector3
 * \param end an iterator to type Vector3
 * Algorithm taken from [Ericson, 2005]
 */
template <class ForwardIterator>
SSR::SSR(ForwardIterator begin, ForwardIterator end)
{
  const unsigned X = 0, Y = 1, Z = 2, THREE_D = 3;
  Vector3 normal;

  // compute the convex hull of the points
  PolyhedronPtr hull = CompGeom::calc_convex_hull_3D(begin, end);
  bool is_2D = !hull;
  if (is_2D)
  {
    // init the 2D centroid
    Vector2 centroid;

    // points are less than 3D; must use SVD
    Real plane_offset;
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
      std::pair<Vector2, Vector2> ep;
      CompGeom::determine_seg_endpoints(points_2D.begin(), points_2D.end(), ep);
      centroid = (ep.first + ep.second) * 0.5;
    }
    else
      // determine the centroid of the hull
      centroid = CompGeom::calc_centroid_2D(hull_2D.begin(), hull_2D.end());

    // project the centroid back to 3D
    R.transpose();
    this->center = CompGeom::to_3D(centroid, R, offset);
    assert(std::fabs(normal.dot(this->center) - plane_offset) < NEAR_ZERO);
  }
  else
  {
    // init the center to origin
    this->center = ZEROS_3;

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
    MatrixNN C(THREE_D);
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
    VectorN evals;
    MatrixNN evecs;
    LinAlg::eig_symm(C, evals, evecs);
    
    // first eigenvector will be direction of minimum variance; that's the
    // one that we want to align with
    Vector3 col;
    evecs.get_column(0, col.begin());
    normal = Vector3::normalize(col);
  }

  // now that we've determined the normal and center, align the rectangle
  Vector3 d2, d3;
  align(begin, end, normal, d2);
  d3 = Vector3::normalize(Vector3::cross(normal, d2));

  // setup R
  this->R.set_column(X, normal);
  this->R.set_column(Y, d2);
  this->R.set_column(Z, d3);

  // compute the lengths and the radius
  calc_lengths_and_radius(begin, end);
}

/// Calculates the lengths and radius for a SSR
template <class ForwardIterator>
void SSR::calc_lengths_and_radius(ForwardIterator begin, ForwardIterator end)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // setup a transform from the SSR frame to the world frame
  Matrix4 T;
  T.set_rotation(&R);
  T.set_translation(center);

  // setup the extents
  Vector2 extents((Real) 0.0, (Real) 0.0);

  // determine the lengths of the rectangle
  for (ForwardIterator i = begin; i != end; i++)
  {
    // transform the point to the SSR frame
    Vector3 p = T.inverse_mult_point(*i);

    // get the y and z points
    Real py = p[Y]; 
    Real pz = p[Z];

    // expand the rectangle as necessary
    if (std::fabs(py) > extents[X])
      extents[X] = std::fabs(py);
    if (std::fabs(pz) > extents[Y])
      extents[Y] = std::fabs(pz);
  }

  // setup the lengths
  this->l = extents*(Real) 2.0;

  // setup the plane of the rectangle
  Vector3 normal = R.get_column(X);
  Real d = normal.dot(center);

  // init radius to zero
  this->radius = 0.0;

  // radius will be greatest planar distance of any point to the rectangle plane
  for (ForwardIterator i = begin; i != end; i++)
  {
    Real sdist = normal.dot(*i) - d;
    if (this->radius < std::fabs(sdist))
      this->radius = std::fabs(sdist);
  }
}

/// Aligns this SSR with the minimum area bounding rectangle projected along the first dimension of the SSR
template <class ForwardIterator>
void SSR::align(ForwardIterator begin, ForwardIterator end, const Vector3& d1, Vector3& d2)
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

/// Outputs the SSR to the specified stream
inline std::ostream& operator<<(std::ostream& o, const SSR& s)
{
  o << "orientation: " << std::endl << s.R;
  o << "center: " << s.center << std::endl;
  o << "lengths: " << s.l << std::endl;
  o << "radius: " << s.radius << std::endl;
  return o;
}

