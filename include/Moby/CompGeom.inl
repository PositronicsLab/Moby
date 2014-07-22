/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

/// Utility method for triangulate_polygon_2D()
template <class BidirectionalIterator>
bool CompGeom::diagonalie(BidirectionalIterator a, BidirectionalIterator b, BidirectionalIterator begin, BidirectionalIterator end, double tol)
{
  BidirectionalIterator c, c1;
  
  c = begin;
  do
  {
    c1 = c;
    c1++;

    // skip edges incident to a or b
    if ((c != a) && (c1 != a) && (c != b) && (c1 != b) &&
        intersect(**a, **b, **c, **c1, tol))
      return false;

    c++;
  }
  while (c != end);

  return true;
}

/// Utility method for triangulate_polygon_2D()
template <class BidirectionalIterator>
bool CompGeom::diagonal(BidirectionalIterator a, BidirectionalIterator b, BidirectionalIterator begin, BidirectionalIterator end, double tol)
{
  return in_cone(a, b, begin, end, tol) && 
          in_cone(b, a, begin, end, tol) && 
          diagonalie(a, b, begin, end, tol);
}

/// Utility method for triangulate_polygon_2D()
template <class BidirectionalIterator>
bool CompGeom::in_cone(BidirectionalIterator a, BidirectionalIterator b, BidirectionalIterator begin, BidirectionalIterator end, double tol)
{
  // get the vertices before and after a
  BidirectionalIterator a0 = (a == begin) ? end : a;
  BidirectionalIterator a1 = a;

  a0--;
  a1++;
  if (a1 == end)
    a1 = begin;

  if (area_sign(**a, **a1, **a0, tol) != eRight)
    return area_sign(**a, **b, **a0, tol) == eLeft && area_sign(**b, **a, **a1, tol) == eLeft;
  else
    return !(area_sign(**a, **b, **a1, tol) != eRight && area_sign(**b, **a, **a0, tol) == eRight);
}

/// Triangulates a polygon (in 2D)
/**
 * \param begin a bidirectional iterator to a container of Point2d objects
 * \param end a bidirectional iterator to a container of Point2d objects
 * \param outbegin an iterator to a container of type 
 *         std::pair<const Point2d*, const Point2d*>; output will be written here
 * \param tol tolerance to test for zero
 * \return an iterator to the end of the output
 * \todo re-implement using faster O(n lg n) algorithm (currently runs in O(n^2) time)
 */
template <class BidirectionalIterator, class OutputIterator>
OutputIterator CompGeom::triangulate_polygon_2D(BidirectionalIterator begin, BidirectionalIterator end, OutputIterator outbegin, double tol)
{
  // make a list out of the polygon; simultaneously get # of vertices
  unsigned n = 0;
  std::list<const Point2d*> poly;
  for (BidirectionalIterator v = begin; v != end; v++, n++)
    poly.push_back(&*v);

  // initialize ear for all vertices
  std::map<const Point2d*, bool> ear;
  for (std::list<const Point2d*>::iterator v1 = poly.begin(); v1 != poly.end(); v1++)
  {
    // get the next vertex
    std::list<const Point2d*>::iterator v2 = v1;
    if (++v2 == poly.end())
      v2 = poly.begin();

    // get the previous vertex
    std::list<const Point2d*>::iterator v0 = (v1 == poly.begin()) ? poly.end() : v1;
    v0--;

    // set the ear
    ear[*v1] = diagonal(v0, v2, poly.begin(), poly.end(), tol);
  }

  // get the number of vertices
  while (n > 3)
  {
    // inner loop searches for an ear
    std::list<const Point2d*>::iterator v2 = poly.begin();
    do
    {
      if (ear[*v2])
      {
        // ear found.  get the next and previous vertices
        std::list<const Point2d*>::iterator v3 = v2;
        std::list<const Point2d*>::iterator v1 = (v2 == poly.begin()) ? poly.end() : v2;
        v1--;
        v3++;
        if (v3 == poly.end())
          v3 = poly.begin();

        // get next and previous vertices v3 and v1, respectively
        std::list<const Point2d*>::iterator v4 = v3;
        std::list<const Point2d*>::iterator v0 = (v1 == poly.begin()) ? poly.end() : v1;
        v0--;
        v4++;
        if (v4 == poly.end())
          v4 = poly.begin();

        // add the diagonal
        *outbegin++ = std::make_pair(*v1, *v3);

        // update earity of diagonal endpoints
        ear[*v1] = diagonal(v0, v3, poly.begin(), poly.end(), tol);
        ear[*v3] = diagonal(v1, v4, poly.begin(), poly.end(), tol);

        // cut off v2
        poly.erase(v2);
        n--;
        break;
      }
      v2++;

    }
    while (v2 != poly.end());
  }

  return outbegin;
}

/// Computes the centroid of a set of facets
/**
 * The facets may represent a polygon, a polyhedron, or even an open polyhedron.  However, the facets may not intersect.
 */
template <class ForwardIterator>
Point3d CompGeom::calc_centroid_3D(ForwardIterator first, ForwardIterator last)
{
  double area_sum = 0;

  // look for no point
  if (first == last)
    return Point3d::zero();
  
  // compute the area of each facet and contribution from each facet
  Point3d centroid = Point3d::zero(first->a.pose);
  for (ForwardIterator i = first; i != last; i++)
  {
    double area = i->calc_area();
    area_sum += area;
    centroid += area * (i->a + i->b + i->c);
  }
  
  // divide the centroid by 3*area_sum
  centroid /= (area_sum * 3);
  
  return centroid;
}

/*****************************************************************************
 Point3d* versions of functions BEGIN
 ****************************************************************************/

/**
 * Determines whether a polygon in 2D is counter-clockwise
 * \note Degenerate polygons (alternating representation) will fail!
 */
template <class ForwardIterator>
bool CompGeomSpecOne<ForwardIterator, Point3d*>::ccw(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double tol)
{
  assert(tol >= 0.0);

  for (ForwardIterator i = begin; i != end; i++)
  {
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;

    ForwardIterator k = j;
    k++;
    if (k == end)
      k = begin; 

    // compute ji and kj
    Ravelin::Vector3d ji = **j - **i;
    Ravelin::Vector3d kj = **k - **j;

    // take the cross product of the normal and the vector j i
    Ravelin::Vector3d c = Ravelin::Vector3d::cross(normal, ji);

    // prepare to determine orientation
    double dprod = c.dot(kj);
    const double TOL = tol * std::max((double) 1.0, std::max(c.norm_inf(), kj.norm_inf()));

    // determine whether k j is to the left or right of j i
    if (dprod > TOL)
      return true;
    else if (dprod < -TOL)
      return false;
    
    // still here? can't be sure - keep going
  }

  // still here?  polygon may be degenerate!
  return true;
}

template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point3d*>::calc_convex_hull(ForwardIterator source_begin, ForwardIterator source_end, OutputIterator target_begin)
{
  const unsigned X = 0, Y = 1, Z = 2;
  int exit_code;
  int curlong, totlong;
  char flags[] = "qhull Fx";
  FILE* outfile, * errfile;
  
  FILE_LOG(LOG_COMPGEOM) << "computing 3D convex hull of following points:" << std::endl;
  for (ForwardIterator i = source_begin; i != source_end; i++)
    FILE_LOG(LOG_COMPGEOM) << "  " << **i << std::endl;
  
  // setup constants for qhull
  const int DIM = 3;
  const int N_POINTS = (int) std::distance(source_begin, source_end);
  const boolT IS_MALLOC = false;
  if (N_POINTS <= 4)
    return target_begin;

  // setup qhull outputs
  if (LOGGING(LOG_COMPGEOM))
  {
    outfile=stdout;  
    errfile=stderr;
  }
  else
  {
    outfile=NULL;
    errfile=fopen("/dev/null", "w");
    assert(errfile);
  } 

  // setup the points
  std::map<coordT*, Point3d*> vertex_map;
  SAFESTATIC std::vector<coordT> qhull_points;
  qhull_points.resize(N_POINTS*DIM);
  coordT* points_begin = &qhull_points.front();
  unsigned j=0;
  for (ForwardIterator i = source_begin; i != source_end; i++)
  {
    qhull_points[j] = (**i)[X];
    qhull_points[j+1] = (**i)[Y];
    qhull_points[j+2] = (**i)[Z];
    vertex_map[points_begin+j] = *i;
    j += DIM;
  }

  // lock the qhull mutex -- qhull is non-reentrant
  #ifdef THREADSAFE
  pthread_mutex_lock(&CompGeom::_qhull_mutex);
  #endif  

  // execute qhull  
  exit_code = qh_new_qhull(DIM, N_POINTS, points_begin, IS_MALLOC, flags, outfile, errfile);
  if (exit_code != 0)
  {
    // points are not collinear.. unsure of the error...
    FILE_LOG(LOG_COMPGEOM) << "CompGeom::calc_convex_hull_3D() - unable to execute qhull on points:" << std::endl;
    for (ForwardIterator i = source_begin; i != source_end; i++)
      FILE_LOG(LOG_COMPGEOM) << "  " << **i << std::endl;

    // free qhull memory
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // release the mutex, since we're not using qhull anymore
    #ifdef THREADSAFE
    pthread_mutex_unlock(&CompGeom::_qhull_mutex);
    #endif

    // close the error stream, if necessary
    if (!LOGGING(LOG_COMPGEOM))
      fclose(errfile);

    throw NumericalException(); 
    return target_begin;
  }

  // iterate through all vertices
  vertexT* vertex;
  FORALLvertices
    (*target_begin++) = vertex_map[vertex->point];
 
  // free qhull memory
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);

  // release the qhull mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&CompGeom::_qhull_mutex);
  #endif 

  // close the error stream, if necessary
  if (!LOGGING(LOG_COMPGEOM))
    fclose(errfile);

  return target_begin;
}

template <class ForwardIterator>
unsigned CompGeomSpecOne<ForwardIterator, Point3d*>::calc_dimensionality(ForwardIterator first, ForwardIterator last, double tol)
{
  assert(tol >= 0.0);

  // make sure we can handle the case of no points
  if (first == last)
    return 0;

  // determine whether all of the points are equal (0 dimensionality)
  ForwardIterator j = first;
  for (ForwardIterator i = first; ; i++)
  {
    // if there are no more points left, everything up to this point
    // has been approximately equal
    j++;
    if (j == last)
      return 0;      

    // if the points are not equal, we can go ahead and break out
    if ((**i - **j).norm() > tol)
      break;
  }  

  // determine whether all of the points are colinear (1 dimensionality)
  // all points from first .. j-1 are coincident, we don't need to check those...
  ForwardIterator k = j;
  while (true)
  {
    // if there are no more points left, everything up to this point has been
    // colinear
    k++;
    if (k == last)
      return 1;

    // if the points are not collinear, we can go ahead and break out
    if (!CompGeom::collinear(**first, **j, **k, tol))
      break;
  }
  
  // determine whether all of the points are coplanar (2 dimensionality)
  // points first, j, k are not colinear, so these will be the basis for our plane
  Ravelin::Vector3d v1 = **j - **first;
  Ravelin::Vector3d v2 = **k - **j;
  Ravelin::Vector3d n = Ravelin::Vector3d::normalize(Ravelin::Vector3d::cross(v1, v2));
  double d = Ravelin::Vector3d::dot(n, v1);
  const double PLANE_TOL = tol * std::max((double) 1.0, std::fabs(d));
  ForwardIterator i = k;
  while (true)
  {
    // if there are no more points left, everything up to this point has been coplanar
    i++;
    if (i == last)
      return 2;

    // if the points are not coplanar, we can go ahead and break out
    if (std::fabs(Ravelin::Vector3d::dot(n, **i) - d) > PLANE_TOL)
      break;
  }

  // still here?  full dimensionality
  return 3;
}

/// Determines the endpoints for a container of collinear Ravelin::Vector3d objects
/**
 * \param begin iterator to beginning of container of type Ravelin::Vector3
 * \param end iterator to end of container of type Ravelin::Vector3
 * \param endpoints the two farthest points on the segment on return
 */
template <class ForwardIterator>
void CompGeom::determine_seg_endpoints(ForwardIterator begin, ForwardIterator end, std::pair<Point3d*, Point3d*>& endpoints)
{
  // make sure that we have been given valid input
  assert(begin != end);

  // setup initial endpoints
  endpoints.first = *begin;
  endpoints.second = *begin;
  double dist = 0; 

  for (ForwardIterator i = ++begin; i != end; i++)
  {
    // get distance from i to both current bounding points
    double dist_e = (**i - *endpoints.second).norm();
    double dist_s = (**i - *endpoints.first).norm();

    // see which distance would be greatest
    if (dist > dist_e)
    {
      // check for no change
      if (dist > dist_s)
        continue;
      else
      {
        dist = dist_s;
        endpoints.second = *i;
      }
    }
    else
    {
      if (dist_e > dist_s)
      {
        dist = dist_e;
        endpoints.first = *i;
      }
      else
      {
        dist = dist_s;
        endpoints.second = *i;
      }
    }
  }
}

/**
 * Converts a collection of Vector3d* objects to Vector2d objects
 * \param source_begin an iterator pointing to the beginning of the Vector3d 
 *        objects
 * \param source_end an iterator pointing to the end of the Vector3d objects
 * \param begin_target an iterator pointing to the beginning of the Vector2d
 *        objects
 * \param R the projection matrix from 3D to 2D
 * \return the end of the output range
 * \note the size of the target collection must be equal to the size of the
 *       source collection
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point3d*>::to_2D(ForwardIterator begin_source, ForwardIterator end_source, OutputIterator begin_target, const Ravelin::Matrix3d& R)
{
  // project the points to 2D
  for (ForwardIterator i = begin_source; i != end_source; i++, begin_target++)
    *begin_target = CompGeom::to_2D(**i, R);

  return begin_target;
}

/// Determines whether a polygon (in 3D) is convex
template <class ForwardIterator>
bool CompGeomSpecOne<ForwardIterator, Point3d*>::is_convex_polygon(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double tol)
{
  const unsigned X = 0, Y = 1, Z = 2;
  assert(tol >= 0.0);

  // get the 3D to 2D projection matrix
  Ravelin::Matrix3d R = CompGeom::calc_3D_to_2D_matrix(normal);

  // project the points to 2D
  std::list<Point2d> points_2D(std::distance(begin, end));
  CompGeom::to_2D(begin, end, points_2D.begin(), R);

  // if the 2D polygon is not ccw, make it so
  assert(CompGeom::ccw(points_2D.begin(), points_2D.end()));
//  if (!ccw(points_2D.begin(), points_2D.end()))
//    std::reverse(points_2D.begin(), points_2D.end());

  // check whether the 2D polygon is convex
  return CompGeom::is_convex_polygon(points_2D.begin(), points_2D.end(), tol);
}

/// Calculates the convex hull of a set of points that lie on a 2D manifold using quickhull
/**
 * \param source_begin an iterator to the beginning of a container of points
 * \param source_end an iterator pointing to the end of a container of points
 * \param normal the (optional) normal of the points; this will be computed if normal is zero vector
 * \param target_begin an iterator to the beginning of a container of points;
 *         on return, contains the convex hull (NOTE: size of this container
 *         must be as large as the source container)
 * \return the new end of the target container
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point3d*>::calc_convex_hull(ForwardIterator source_begin, ForwardIterator source_end, const Ravelin::Vector3d& normal, OutputIterator target_begin)
{  
  FILE_LOG(LOG_COMPGEOM) << "computing 2D convex hull of following points:" << std::endl;
  for (ForwardIterator i = source_begin; i != source_end; i++)
    FILE_LOG(LOG_COMPGEOM) << "  " << **i << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "using normal: " << normal << std::endl;
  
  // **************************************************************
  // first, we need to project the 3D surface to a 2D polygon
  // **************************************************************
  
  // determine the normal, if necessary 
  Ravelin::Vector3d n = normal;
  if (std::fabs(n.norm() - 1.0) > NEAR_ZERO)
  {
    double offset;
    CompGeom::fit_plane(source_begin, source_end, n, offset);
  }

  // compute the 3D to 2D projection matrix
  Ravelin::Matrix3d R = CompGeom::calc_3D_to_2D_matrix(n);

  // get the 2D to 3D offset
  Ravelin::Origin3d p1(**source_begin);
  double offset = CompGeom::determine_3D_to_2D_offset(p1, R);

  // get the transpose (i.e., inverse) of the rotation matrix
  Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);

  // project the points to 2D
  unsigned sz = std::distance(source_begin, source_end);
  std::vector<Point2d> points_2D(sz);
  CompGeom::to_2D(source_begin, source_end, points_2D.begin(), R);

  // compute correspondences
  std::vector<Point2d*> points_2D_ptr(sz);
  std::map<Point2d*, Point3d*> mapping;
  unsigned i=0;
  for (ForwardIterator j=source_begin; j != source_end; i++, j++)
  {
    points_2D_ptr[i] = &points_2D[i];
    mapping[&points_2D[i]] = *j;
  } 
  FILE_LOG(LOG_COMPGEOM) << "2D points:" << std::endl;
  for (unsigned i=0; i< points_2D.size(); i++)
    FILE_LOG(LOG_COMPGEOM) << "  " << points_2D[i] << std::endl;

  // compute the convex hull
  std::list<Point2d*> hull(sz);
  std::list<Point2d*>::iterator hull_end = CompGeom::calc_convex_hull(points_2D_ptr.begin(), points_2D_ptr.end(), hull.begin());

  // use the mapping to 3D
  std::list<Point3d*> hull3D;
  for (std::list<Point2d*>::iterator i = hull.begin(); i != hull_end; i++)
    hull3D.push_back(mapping[*i]);

  // reverse the hull if necessary
  if (!CompGeom::ccw(hull3D.begin(), hull3D.end(), normal))
    std::reverse(hull3D.begin(), hull3D.end());  

  // return the hull
  return std::copy(hull3D.begin(), hull3D.end(), target_begin);
}

/*****************************************************************************
 Point3d* versions of functions END 
 ****************************************************************************/

/*****************************************************************************
 Ravelin::Vector3d versions of functions BEGIN
 ****************************************************************************/

/**
 * Determines whether a polygon in 2D is counter-clockwise
 * \note Degenerate polygons (alternating representation) will fail!
 */
template <class ForwardIterator>
bool CompGeomSpecOne<ForwardIterator, Point3d>::ccw(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double tol)
{
  assert(tol >= 0.0);

  for (ForwardIterator i = begin; i != end; i++)
  {
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;

    ForwardIterator k = j;
    k++;
    if (k == end)
      k = begin; 

    // compute ji and kj
    Ravelin::Vector3d ji = *j - *i;
    Ravelin::Vector3d kj = *k - *j;

    // take the cross product of the normal and the vector j i
    Ravelin::Vector3d c = Ravelin::Vector3d::cross(normal, ji);

    // prepare to determine orientation
    double dprod = c.dot(kj);
    const double TOL = tol * std::max((double) 1.0, std::max(c.norm_inf(), kj.norm_inf()));

    // determine whether k j is to the left or right of j i
    if (dprod > TOL)
      return true;
    else if (dprod < -TOL)
      return false;
    
    // still here? can't be sure - keep going
  }

  // still here?  polygon may be degenerate!
  return true;
}

/// Attempts to fit a plane to a set of points 
/**
 * The singular value decomposition is used to determine the plane that fits
 * the points best in a least-squares sense.
 * \param points the set of points (in 3D)
 * \param normal contains the "best" normal, on return
 * \param offset the offset such that, for any point on the plane x, 
 *        <normal, x> = offset
 * \return the maximum deviation from the plane
 */
template <class ForwardIterator>
double CompGeomSpecOne<ForwardIterator, Point3d*>::fit_plane(ForwardIterator begin, ForwardIterator end, Ravelin::Vector3d& normal, double& offset)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  if (begin == end)
  {
    normal = Ravelin::Vector3d::zero();
    offset = 0.0;
    return 0.0;
  }

  // compute the mean of the data
  unsigned n = 0;
  Point3d mu = Point3d::zero((*begin)->pose);
  for (ForwardIterator i = begin; i != end; i++, n++)
    mu += **i;
  mu /= n;

  // create a matrix subtracting each point from the mean
  SAFESTATIC Ravelin::LinAlgd _LA;
  SAFESTATIC Ravelin::MatrixNd M, U, V;
  SAFESTATIC Ravelin::VectorNd S;
  M.resize(n, THREE_D);
  unsigned idx = 0;
  for (ForwardIterator i = begin; i != end; i++)
  {
    M.set_row(idx, **i - mu);
    idx++;
  }

  // take the svd of the matrix
  _LA.svd(M, U, S, V);

  // last column of V should have the singular value we want; normalize it just in case
  normal.pose = (*begin)->pose;
  normal[X] = V(X,Z);
  normal[Y] = V(Y,Z);
  normal[Z] = V(Z,Z);
  normal.normalize();
  
  // determine offset
  offset = Ravelin::Vector3d::dot(normal, mu);

  // compute distance from all points
  double max_dev = 0;
  for (ForwardIterator i = begin; i != end; i++)
    max_dev = std::max(max_dev, std::fabs(Ravelin::Vector3d::dot(normal, **i) - offset));

  return max_dev;
}

/// Attempts to fit a plane to a set of points 
/**
 * The singular value decomposition is used to determine the plane that fits
 * the points best in a least-squares sense.
 * \param points the set of points (in 3D)
 * \param normal contains the "best" normal, on return
 * \param offset the offset such that, for any point on the plane x, 
 *        <normal, x> = offset
 * \return the maximum deviation from the plane
 */
template <class ForwardIterator>
double CompGeomSpecOne<ForwardIterator, Point3d>::fit_plane(ForwardIterator begin, ForwardIterator end, Ravelin::Vector3d& normal, double& offset)
{
  const unsigned THREE_D = 3, X = 0, Y = 1, Z = 2;
  
  if(begin == end){
      normal = Ravelin::Vector3d::zero();
      offset = 0;
      return 0;
  }
  // compute the mean of the data
  unsigned n = 0;
  Point3d mu = Point3d::zero(begin->pose);
  
  for (ForwardIterator i = begin; i != end; i++, n++)
    mu += *i;
  mu /= n;

  // create a matrix subtracting each point from the mean
  SAFESTATIC Ravelin::LinAlgd _LA;
  SAFESTATIC Ravelin::MatrixNd M, U, V;
  SAFESTATIC Ravelin::VectorNd S;
  M.resize(n, THREE_D);
  unsigned idx = 0;
  for (ForwardIterator i = begin; i != end; i++)
  {
    M.set_row(idx, *i - mu);
    idx++;
  }

  // take the svd of the matrix
  _LA.svd(M, U, S, V);

  // last column of V should have the singular value we want; normalize it just in case
  normal.pose = begin->pose;
  normal[X] = V(X,Z);
  normal[Y] = V(Y,Z);
  normal[Z] = V(Z,Z);
  normal.normalize();
  
  // determine offset
  offset = Ravelin::Vector3d::dot(normal, mu);

  // compute distance from all points
  double max_dev = 0;
  for (ForwardIterator i = begin; i != end; i++)
    max_dev = std::max(max_dev, std::fabs(Ravelin::Vector3d::dot(normal, *i) - offset));

  return max_dev;
}

/// Computes the 3D convex hull of a set of points
/**
 * \param first a forward iterator for type Ravelin::Vector3*
 * \param last a forward iterator for type Ravelin::Vector3*
 * \return a pointer to the newly created polyhedron
 */
template <class ForwardIterator>
PolyhedronPtr CompGeomSpecOne<ForwardIterator, Point3d>::calc_convex_hull(ForwardIterator first, ForwardIterator last)
{
  const unsigned X = 0, Y = 1, Z = 2;
  int exit_code;
  int curlong, totlong;
  FILE* outfile, * errfile;

  // setup qhull outputs
  if (LOGGING(LOG_COMPGEOM))
  {
    outfile=stdout;  
    errfile=stderr;
  }
  else
  {
    outfile=NULL;
    errfile=fopen("/dev/null", "w");
    assert(errfile);
  } 

  // determine how many points we are processing
  unsigned sz = 0;
  for (ForwardIterator i = first; i != last; i++)
    sz++;

  // setup constants
  const int DIM = 3;
  const int N_POINTS = sz;
  const boolT IS_MALLOC = false;
  char flags[] = "qhull ";

  // make sure there are enough points
  if (N_POINTS < 4)
  {
    if (!LOGGING(LOG_COMPGEOM))
      fclose(errfile);

    return PolyhedronPtr();
  }
  
  // setup the points
  SAFESTATIC std::vector<coordT> qhull_points;
  qhull_points.resize(N_POINTS*DIM);
  coordT* points_begin = &qhull_points.front(); 
  unsigned j=0;
  for (ForwardIterator i = first; i != last; i++)
  {
    qhull_points[j] = (*i)[X];
    qhull_points[j+1] = (*i)[Y];
    qhull_points[j+2] = (*i)[Z];
    j += DIM;
  }
  FILE_LOG(LOG_COMPGEOM) << "computing 3D convex hull of: " << std::endl;
  for (ForwardIterator i = first; i != last; i++)
    FILE_LOG(LOG_COMPGEOM) << *i << std::endl;

  // lock the qhull mutex -- qhull is non-reentrant
  #ifdef THREADSAFE
  pthread_mutex_lock(&CompGeom::_qhull_mutex);
  #endif  

  // execute qhull  
  exit_code = qh_new_qhull(DIM, N_POINTS, points_begin, IS_MALLOC, flags, outfile, errfile);
  if (exit_code)
  {
    // free qhull memory
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // qhull failed -- perhaps the dimensionality is 2 rather than 3?
    #ifdef THREADSAFE
    pthread_mutex_unlock(&CompGeom::_qhull_mutex);
    #endif

    // close the error stream, if necessary
    if (!LOGGING(LOG_COMPGEOM))
      fclose(errfile);

    throw NumericalException(); 
    return PolyhedronPtr();
  }

  // construct a new vector of vertices
  std::vector<Point3d> vertices; 

  // get all of the vertices
  std::map<vertexT*, unsigned> vertex_map;
  for (vertexT* vertex=qh vertex_list;vertex && vertex->next;vertex= vertex->next)
  {
    Point3d v;
    for (unsigned i=0; i< (unsigned) DIM; i++)
      v[i] = vertex->point[i];
    vertices.push_back(v);
    vertex_map[vertex] = vertices.size()-1;
  }

  // triangulate
  qh_triangulate();
 
  // setup list of facets
  std::list<IndexedTri> facets;
 
  // get the facet information
  for (facetT* facet=qh facet_list;facet && facet->next;facet=facet->next)
  {
    if (!facet->vertices)
      continue;

    // setup a facet
    facets.push_back(IndexedTri());
    
    // get all (three) vertices in the facet
    vertexT** vertex_pointer = (vertexT**)& ((facet->vertices)->e[0].p); 
    vertexT* vertex = *vertex_pointer++;
    assert(vertex);
    facets.back().a = vertex_map[vertex];
    vertex = *vertex_pointer++;
    assert(vertex);
    facets.back().b = vertex_map[vertex];
    vertex = *vertex_pointer++;
    assert(vertex);
    facets.back().c = vertex_map[vertex];

    // must be no more than three vertices in the facet
    assert(!*vertex_pointer);

    // reverse order of the vertices if necessary 
    if (facet->toporient ^ qh_ORIENTclock)
      std::swap(facets.back().b, facets.back().c);
  }
  
  // free qhull memory
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
  assert(!curlong && !totlong);
  
  // release the qhull mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&CompGeom::_qhull_mutex);
  #endif  

  // if the there aren't enough triangles, can't create the polyhedron
  assert(facets.size() >= 4);

  // create the polyhedron and verify that it is consistent
  PolyhedronPtr polyhedron(new Polyhedron(vertices.begin(), vertices.end(), facets.begin(), facets.end()));
  assert(polyhedron->consistent());

  FILE_LOG(LOG_COMPGEOM) << "3D convex hull is:" << std::endl << *polyhedron;

  // close the error stream, if necessary
  if (!LOGGING(LOG_COMPGEOM))
    fclose(errfile);

  return polyhedron;  
}


/// Determines the dimensionality of a set of points
/**
 * \return the dimensionality (0 [point], 1 [line], 2 [plane], 3 [full space])
 */
template <class ForwardIterator>
unsigned CompGeomSpecOne<ForwardIterator, Point3d>::calc_dimensionality(ForwardIterator first, ForwardIterator last, double tol)
{
  assert(tol >= 0.0);

  // make sure we can handle the case of no points
  if (first == last)
    return 0;

  // determine whether all of the points are equal (0 dimensionality)
  ForwardIterator j = first;
  for (ForwardIterator i = first; ; i++)
  {
    // if there are no more points left, everything up to this point
    // has been approximately equal
    j++;
    if (j == last)
      return 0;      

    // if the points are not equal, we can go ahead and break out
    if ((*i - *j).norm() > tol)
      break;
  }  

  // determine whether all of the points are colinear (1 dimensionality)
  // all points from first .. j-1 are coincident, we don't need to check those...
  ForwardIterator k = j;
  while (true)
  {
    // if there are no more points left, everything up to this point has been
    // colinear
    k++;
    if (k == last)
      return 1;

    // if the points are not collinear, we can go ahead and break out
    if (!CompGeom::collinear(*first, *j, *k, tol))
      break;
  }
  
  // determine whether all of the points are coplanar (2 dimensionality)
  // points first, j, k are not colinear, so these will be the basis for our plane
  Ravelin::Vector3d v1 = *j - *first;
  Ravelin::Vector3d v2 = *k - *j;
  Ravelin::Vector3d n = Ravelin::Vector3d::normalize(Ravelin::Vector3d::cross(v1, v2));
  double d = Ravelin::Vector3d::dot(n, v1);  
  const double PLANE_TOL = tol * std::max((double) 1.0, std::fabs(d));
  ForwardIterator i = k;
  while (true)
  {
    // if there are no more points left, everything up to this point has been coplanar
    i++;
    if (i == last)
      return 2;

    // if the points are not coplanar, we can go ahead and break out
    if (std::fabs(Ravelin::Vector3d::dot(n, *i) - d) > PLANE_TOL)
      break;
  }

  // still here?  full dimensionality
  return 3;
}

/// Determines the endpoints for a container of collinear Ravelin::Vector3d objects
/**
 * \param begin iterator to beginning of container of type Ravelin::Vector3
 * \param end iterator to end of container of type Ravelin::Vector3
 * \param endpoints the two farthest points on the segment on return
 */
template <class ForwardIterator>
void CompGeom::determine_seg_endpoints(ForwardIterator begin, ForwardIterator end, std::pair<Point3d, Point3d>& endpoints)
{
  // make sure that we have been given valid input
  assert(begin != end);

  // setup initial endpoints
  endpoints.first = *begin;
  endpoints.second = *begin;
  double dist = 0; 

  for (ForwardIterator i = ++begin; i != end; i++)
  {
    // get distance from i to both current bounding points
    double dist_e = (*i - endpoints.second).norm();
    double dist_s = (*i - endpoints.first).norm();

    // see which distance would be greatest
    if (dist > dist_e)
    {
      // check for no change
      if (dist > dist_s)
        continue;
      else
      {
        dist = dist_s;
        endpoints.second = *i;
      }
    }
    else
    {
      if (dist_e > dist_s)
      {
        dist = dist_e;
        endpoints.first = *i;
      }
      else
      {
        dist = dist_s;
        endpoints.second = *i;
      }
    }
  }
}

/**
 * Converts a collection of Ravelin::Vector3d objects to Vector2 objects
 * \param source_begin an iterator pointing to the beginning of the Ravelin::Vector3d 
 *        objects
 * \param source_end an iterator pointing to the end of the Ravelin::Vector3d objects
 * \param begin_target an iterator pointing to the beginning of the Vector2
 *        objects
 * \param R the projection matrix from 3D to 2D (on return)
 * \return the end of the output range
 * \note the size of the target collection must be equal to the size of the
 *       source collection
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point3d>::to_2D(ForwardIterator begin_source, ForwardIterator end_source, OutputIterator begin_target, const Ravelin::Matrix3d& R)
{
  // project the points to 2D
  for (ForwardIterator i = begin_source; i != end_source; i++, begin_target++)
    *begin_target = CompGeom::to_2D(*i, R);

  return begin_target;
}

/// Determines whether a polygon (in 3D) is convex
template <class ForwardIterator>
bool CompGeomSpecOne<ForwardIterator, Point3d>::is_convex_polygon(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double tol)
{
  const unsigned X = 0, Y = 1, Z = 2;
  assert(tol >= 0.0);

  // get the 3D to 2D projection matrix
  Ravelin::Matrix3d R = CompGeom::calc_3D_to_2D_matrix(normal);

  // project the points to 2D
  std::list<Point2d> points_2D(std::distance(begin, end));
  CompGeom::to_2D(begin, end, points_2D.begin(), R);

  // if the 2D polygon is not ccw, make it so
  assert(CompGeom::ccw(points_2D.begin(), points_2D.end()));
//  if (!ccw(points_2D.begin(), points_2D.end()))
//    std::reverse(points_2D.begin(), points_2D.end());

  // check whether the 2D polygon is convex
  return CompGeom::is_convex_polygon(points_2D.begin(), points_2D.end(), tol);
}

/// Calculates the convex hull of a set of points that lie on a 2D manifold using quickhull
/**
 * \param source_begin an iterator to the beginning of a container of points
 * \param source_end an iterator pointing to the end of a container of points
 * \param normal the (optional) normal of the points; this will be computed if normal is zero vector
 * \param target_begin an iterator to the beginning of a container of points;
 *         on return, contains the convex hull (NOTE: size of this container
 *         must be as large as the source container)
 * \return the new end of the target container
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point3d>::calc_convex_hull(ForwardIterator source_begin, ForwardIterator source_end, const Ravelin::Vector3d& normal, OutputIterator target_begin)
{  
  FILE_LOG(LOG_COMPGEOM) << "computing 2D convex hull of following points:" << std::endl;
  for (ForwardIterator i = source_begin; i != source_end; i++)
    FILE_LOG(LOG_COMPGEOM) << "  " << *i << std::endl;
  
  // **************************************************************
  // first, we need to project the 3D surface to a 2D polygon
  // **************************************************************
  
  // determine the normal, if necessary 
  Ravelin::Vector3d n = normal;
  if (std::fabs(n.norm() - 1.0) > NEAR_ZERO)
  {
    double offset;
    CompGeom::fit_plane(source_begin, source_end, n, offset);
  }

  // compute the 3D to 2D projection matrix
  Ravelin::Matrix3d R = CompGeom::calc_3D_to_2D_matrix(n);

  // get the 2D to 3D offset
  double offset = CompGeom::determine_3D_to_2D_offset(*source_begin, R);

  // get the transpose (i.e., inverse) of the rotation matrix
  Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);

  // project the points to 2D
  unsigned sz = std::distance(source_begin, source_end);
  std::vector<Point2d> points_2D(sz);
  CompGeom::to_2D(source_begin, source_end, R, points_2D);

  // compute correspondences
  std::vector<Point2d*> points_2D_ptr(sz);
  std::map<Point2d*, const Point3d*> mapping;
  unsigned i=0;
  for (ForwardIterator j=source_begin; j != source_end; i++, j++)
  {
    points_2D_ptr[i] = &(*j);
    mapping[&points_2D[i]] = &(*j);
  } 
  FILE_LOG(LOG_COMPGEOM) << "2D points:" << std::endl;
  for (unsigned i=0; i< points_2D.size(); i++)
    FILE_LOG(LOG_COMPGEOM) << "  " << points_2D[i] << std::endl;

  // compute the convex hull
  std::list<Point2d*> hull(sz);
  std::list<Point2d*>::iterator hull_end = calc_convex_hull(points_2D_ptr.begin(), points_2D_ptr.end(), hull.begin());

  // use the mapping to 3D
  std::list<Point3d> hull3D;
  for (std::list<Point2d*>::iterator i = hull.begin(); i != hull_end; i++)
    hull3D.push_back(*mapping[*i]);

  // reverse the hull if necessary
  if (!CompGeom::ccw(hull3D.begin(), hull3D.end(), normal))
    std::reverse(hull3D.begin(), hull3D.end());  

  // return the hull
  return std::copy(hull3D.begin(), hull3D.end(), target_begin);
} 
 
/*****************************************************************************
 Ravelin::Vector3d versions of functions END 
 ****************************************************************************/

/*****************************************************************************
 Point2d versions of functions BEGIN
 ****************************************************************************/

/// Computes the minimum area bounding rectangle of a set of points
/**
 * Uses 2D convex hull and rotating calipers method; runs in O(N lg N) time
 * On return, x1, x2, x3, and x4 are the four vertices of the bounding
 * rectangle (ordered as edges).
 */
template <class ForwardIterator>
void CompGeomSpecOne<ForwardIterator, Point2d>::calc_min_area_bounding_rect(ForwardIterator begin, ForwardIterator end, Point2d& x1, Point2d& x2, Point2d& x3, Point2d& x4)
{
  enum Flag { F_NONE, F_LEFT, F_RIGHT, F_BOTTOM, F_TOP };
  const unsigned X = 0, Y = 1;
  std::pair<Point2d, Point2d> ep;

  // calculate the convex hull of the points in ccw order
  std::vector<Point2d> points, hull;
  for (; begin != end; begin++)
    points.push_back(*begin);
  CompGeom::calc_convex_hull(points.begin(), points.end(), std::back_inserter(hull));
  if (hull.empty())
  {
    // convex hull is degenerate; compute line endpoints and make that the 
    // "hull"
    CompGeom::determine_seg_endpoints(begin, end, ep);
    hull.push_back(ep.first);
    hull.push_back(ep.second);
  }
  // get the hull in CCW order
  else if (!CompGeom::ccw(hull.begin(), hull.end()))
    std::reverse(hull.begin(), hull.end());

  // make sure that no three points are colinear
  unsigned n = hull.size();
  for (unsigned i = 0; i < n; i++)
  {
    // get the next two points
    unsigned j = (i < n-1) ? i+1 : 0;
    unsigned k = (j < n-1) ? j+1 : 0;
    if (CompGeom::collinear(hull[i], hull[j], hull[k]))
    {
      // erase the middle point
      hull.erase(hull.begin()+j);

      // decrement both n and i
      i--;
      n--;
      // if n < 3, bounding rectangle is degenerate; however, we can still
      // output something useful
      if (n < 3)
      {
        x1 = hull[0];
        x2 = hull[1];
        x3 = hull[0];
        x4 = hull[1];
        return;
      }
    }
  }

  // setup unit-length edge directions of the convex polygon
  unsigned nm1 = n - 1;
  std::vector<Ravelin::Vector2d> edges(n);
  std::vector<bool> visited(n, false);
  for (unsigned i=0; i< nm1; i++)
  {
    edges[i] = hull[i+1] - hull[i];
    edges[i].normalize();
  }
  edges[nm1] = hull[0] - hull[nm1];
  edges[nm1].normalize();

  // find the smallest axis-aligned box containing the points.  Keep track
  // of the extremum indices, L (left), R (right), B (bottom), and T (top)
  // so that the following constraints are met:
  //   V[L].X <= V[i].X for all i and V[(L+1)%N].X > V[L].X
  //   V[L].X >= V[i].X for all i and V[(R+1)%N].X < V[R].X
  //   V[L].Y <= V[i].Y for all i and V[(B+1)%N].Y > V[B].X
  //   V[L].Y >= V[i].Y for all i and V[(T+1)%N].Y < V[T].X
  double xmin = hull[0][X], xmax = xmin;
  double ymin = hull[0][Y], ymax = ymin;
  unsigned Lindex = 0, Rindex = 0, Bindex = 0, Tindex = 0;
  for (unsigned i=1; i< n; i++)
  {
    if (hull[i][X] <= xmin)
    {
      xmin = hull[i][X];
      Lindex = i;
    }
    if (hull[i][X] >= xmax)
    {
      xmax = hull[i][X];
      Rindex = i;
    }
    if (hull[i][Y] <= ymin)
    {
      ymin = hull[i][Y];
      Bindex = i;
    }
    if (hull[i][Y] >= ymax)
    {
      ymax = hull[i][Y];
      Tindex = i;
    }
  }

  // apply wrap-around tests to ensure the constraints mentioned above 
  // are satisfied
  if (Lindex == nm1)
  {
    if (hull[0][X] <= xmin)
    {
      xmin = hull[0][X];
      Lindex = 0;
    }
  }
  if (Rindex == nm1)
  {
    if (hull[0][X] >= xmax)
    {
      xmax = hull[0][X];
      Rindex = 0;
    }
  }
  if (Bindex == nm1)
  {
    if (hull[0][Y] <= ymin)
    {
      ymin = hull[0][Y];
      Bindex = 0;
    }
  }
  if (Tindex == nm1)
  {
    if (hull[0][Y] >= ymax)
    {
      ymax = hull[0][Y];
      Tindex = 0;
    }
  }

  // the dimensions of the axis-aligned box; the extents store width and height
  // for now
  Point2d center((double) 0.5 * (xmin + xmax), (double) 0.5 * (ymin + ymax));
  Ravelin::Vector2d axis[2] = { Ravelin::Vector2d((double) 1.0, (double) 0.0), 
                      Ravelin::Vector2d((double) 0.0, (double) 1.0) };
  double extent[2] = { (double) 0.5 * (xmax - xmin), (double) 0.5 * (ymax - ymin) };
  double min_area_div4 = extent[0]*extent[1];

  // the rotating calipers algorithm follows...
  Ravelin::Vector2d U((double) 1.0, (double) 0.0), V((double) 0.0, (double) 1.0);
  bool done = false;
  while (!done)
  {
    // determine the edge that forms the smallest angle with the current box
    // edges
    Flag flag = F_NONE;
    double maxdot = (double) 0.0;
    double dot = U.dot(edges[Bindex]);
    if (dot > maxdot)
    {
      maxdot = dot;
      flag = F_BOTTOM;
    }
    dot = V.dot(edges[Rindex]);
    if (dot > maxdot)
    {
      maxdot = dot;
      flag = F_RIGHT;
    }
    dot = -U.dot(edges[Tindex]);
    if (dot > maxdot)
    {
      maxdot = dot;
      flag = F_TOP;
    }
    dot = -V.dot(edges[Lindex]);
    if (dot > maxdot)
    {
      maxdot = dot;
      flag = F_LEFT;
    }

    switch (flag)
    {
      case F_BOTTOM:
        if (visited[Bindex])
          done = true;
        else
        {
          // compute box axes with E[B] as an edge
          U = edges[Bindex];
          V = -U.perp();
          CompGeom::update_box(hull[Lindex], hull[Rindex], hull[Bindex], hull[Tindex], U, V, min_area_div4, center, axis, extent);

          // mark edge visited and rotate the calipers
          visited[Bindex] = true;
          if (++Bindex == n)
             Bindex = 0;
         }
         break;

      case F_RIGHT:
        if (visited[Rindex])
          done = true;
        else
        {
          // compute box axes with E[R] as an edge
          V = edges[Rindex];
          U = V.perp();
          CompGeom::update_box(hull[Lindex], hull[Rindex], hull[Bindex], hull[Tindex], U, V, min_area_div4, center, axis, extent);

          // mark edge visited and rotate the calipers
          visited[Rindex] = true;
          if (++Rindex == n)
            Rindex = 0;
        }
        break;

      case F_TOP:
        if (visited[Tindex])
          done = true;
        else
        {
          // compute box axes with E[T] as an edge
          U = -edges[Tindex];
          V = -U.perp();
          CompGeom::update_box(hull[Lindex], hull[Rindex], hull[Bindex], hull[Tindex], U, V, min_area_div4, center, axis, extent);

          // mark edge visited and rotate the calipers
          visited[Tindex] = true;
          if (++Tindex == n)
            Tindex = 0;
        }
        break;

      case F_LEFT:
        if (visited[Lindex])
          done = true;
        else
        {
          // compute box axes with E[L] as an edge
          V = -edges[Lindex];
          U = V.perp();
          CompGeom::update_box(hull[Lindex], hull[Rindex], hull[Bindex], hull[Tindex], U, V, min_area_div4, center, axis, extent);

          // mark edge visited and rotate the calipers
          visited[Lindex] = true;
          if (++Lindex == n)
            Lindex = 0;
        }
        break;

      case F_NONE:
        // the polygon is a rectangle
        done = true;
        break;
    }
  }

  // convert Eberly's representation to my own
  x1 = center - axis[X]*extent[X] - axis[Y]*extent[Y];
  x2 = center + axis[X]*extent[X] - axis[Y]*extent[Y];
  x3 = center + axis[X]*extent[X] + axis[Y]*extent[Y];
  x4 = center - axis[X]*extent[X] + axis[Y]*extent[Y];
}

/// Determines the endpoints for a container of collinear Point2d objects
/**
 * \param begin iterator to beginning of container of type Point2d
 * \param end iterator to end of container of type Point2d
 * \param endpoints the two farthest points on the segment on return
 */
template <class ForwardIterator>
void CompGeom::determine_seg_endpoints(ForwardIterator begin, ForwardIterator end, std::pair<Point2d, Point2d>& endpoints)
{
  // make sure that we have been given valid input
  assert(begin != end);

  // setup initial endpoints
  endpoints.first = *begin;
  endpoints.second = *begin;
  double dist = 0; 

  for (ForwardIterator i = ++begin; i != end; i++)
  {
    // get distance from i to both current bounding points
    double dist_e = (*i - endpoints.second).norm();
    double dist_s = (*i - endpoints.first).norm();

    // see which distance would be greatest
    if (dist > dist_e)
    {
      // check for no change
      if (dist > dist_s)
        continue;
      else
      {
        dist = dist_s;
        endpoints.second = *i;
      }
    }
    else
    {
      if (dist_e > dist_s)
      {
        dist = dist_e;
        endpoints.first = *i;
      }
      else
      {
        dist = dist_s;
        endpoints.second = *i;
      }
    }
  }
}

/// Computes the intersection of a polygon and a line segment
/**
 * \param begin iterator pointing to start of collection of Point2d elements
 * \param end iterator pointing to end of collection of Point2d elements
 * \param seg the line segment
 * \param outbegin an iterator to a container of type LineSeg2 that will 
 *        store line segments on/in the polygon on return 
 * \param tol the tolerance for parallel lines
 * \return ePolygonSegInside if the polygon 
 * \note polygon must be given in counter-clockwise order
 * \note first vertex should not appear twice
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point2d>::intersect_seg_polygon(ForwardIterator begin, ForwardIterator end, const LineSeg2& seg, OutputIterator outbegin)
{
  std::list<double> points;

  // determine whether one (or both) of the endpoints is within the polygon
  if (CompGeom::polygon_location(begin, end, seg.first) == CompGeom::ePolygonInside)
    points.push_back(0.0);
  if (CompGeom::polygon_location(begin, end, seg.second) == CompGeom::ePolygonInside)
    points.push_back(1.0);

  // determine the inverse of the length (squared) of the line segment
  double inv_seg_len_sq = 1.0/(seg.first - seg.second).norm_sq();

  // intersect all line segments
  Point2d isect1, isect2;
  for (ForwardIterator i = begin; i != end; i++)
  {
    // get the next vertex
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;
  
    // intersect the two segments
    switch (CompGeom::intersect_segs(seg, LineSeg2(*i, *j), isect1, isect2))
    {
      case CompGeom::eSegSegNoIntersect:
        break;

      case CompGeom::eSegSegIntersect:
      case CompGeom::eSegSegVertex: 
        points.push_back((isect1 - seg.second).norm_sq() * inv_seg_len_sq);
        break;

      case CompGeom::eSegSegEdge:
        points.push_back((isect1 - seg.second).norm_sq() * inv_seg_len_sq);
        points.push_back((isect2 - seg.second).norm_sq() * inv_seg_len_sq);
        break;
    }
  }

  // sort the points
  points.sort();

  // make segments out of the points
  for (std::list<double>::const_iterator i = points.begin(); i != points.rbegin().base(); i++)
  {
    std::list<double>::const_iterator j = i;
    j++;
    Point2d p1 = seg.first * (*i) + seg.second * (1.0 - (*i));
    Point2d p2 = seg.first * (*j) + seg.second * (1.0 - (*j));
    *outbegin++ = LineSeg2(p1, p2);
  }

  return outbegin;
}

/// Computes the intersection of a convex polygon and a line segment
/**
 * \param begin iterator pointing to start of collection of Point2d elements
 * \param end iterator pointing to end of collection of Point2d elements
 * \param seg the line segment
 * \param te the parameter of the line segment for the beginning of the
 *        intersection; (1-te)*seg.first + te*seg.second is point of 
 *        intersection
 * \param tl the parameter of the line segment for the end of the intersection;
 *        (1-tl)*seg.first + tl*seg.second is point of intersection
 * \param tol the tolerance for parallel lines
 * \return <b>true</b> if the two intersect and <b>false</b> otherwise
 * \note polygon must be given in counter-clockwise order
 * \note first vertex should not appear twice
 * \note taken from http://geometryalgorithms.com 
 */
template <class ForwardIterator>
bool CompGeomSpecOne<ForwardIterator, Point2d>::intersect_seg_convex_polygon(ForwardIterator begin, ForwardIterator end, const LineSeg2& seg, double& te, double& tl, double tol)
{
  assert(tol >= 0.0);

  // initialize te and tl
  te = 0;
  tl = 1;

  Point2d dS = seg.second - seg.first;

  // iterate over all vertices
  for (ForwardIterator i = begin; i != end; i++)
  {
    // get the next vertex
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;

    // get the edge
    Ravelin::Vector2d edge = *j - *i;

    // determine the outward normal of the edge
    Ravelin::Vector2d ni(edge[1], -edge[0]);
    double N = -Ravelin::Vector2d::dot(ni, seg.first-*i);
    double D = Ravelin::Vector2d::dot(dS, ni);
    if (std::fabs(D) < 0.0)
    {
      // segment is parallel to this edge
      if (N < tol)
        // first point is outside of the edge; segment cannot intersect poly
        return false;
      else
        // segment cannot enter or leave poly across this edge, process next
        continue;
    }

    double t = N / D;
    if (D < 0.0)
    {
      // segment enters polygon across this edge
      te = std::max(te, t);
      if (te > tl)
        // segment enters polygon after leaving
        return false;
    }
    else
    {
      assert(D > 0.0);
      tl = std::min(tl, t);
      if (tl < te)
        // segment leaves polygon before entering
        return false;
    }  
  }

  return true;
}

/**
 * Converts a collection of Vector2 objects to Ravelin::Vector3d objects
 * \param source_begin an iterator pointing to the beginning of the Ravelin::Vector3d 
 *        objects
 * \param source_end an iterator pointing to the end of the Ravelin::Vector3d objects
 * \param begin_target an iterator pointing to the beginning of the Vector2
 *        objects
 * \param R the projection matrix from 2D to 3D
 * \param offset the offset that must be added to the Z-coordinate of points
 *        projected back from 2D to 3D
 * \return the end of the output range
 * \note the size of the target collection must be equal to the size of the
 *       source collection
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point2d>::to_3D(ForwardIterator begin_source, ForwardIterator end_source, OutputIterator begin_target, const Ravelin::Matrix3d& RT, double offset)
{
  // project the points back to 3D
  for (ForwardIterator i = begin_source; i != end_source; i++, begin_target++)
    *begin_target = CompGeom::to_3D(*i, RT, offset);

  return begin_target;
}

/**
 * Determines whether a polygon in 2D is counter-clockwise
 * \note Degenerate polygons (alternating representation) will fail!
 */
template <class ForwardIterator>
bool CompGeomSpecOne<ForwardIterator, Point2d>::ccw(ForwardIterator begin, ForwardIterator end, double tol)
{
  assert(tol >= 0.0);

  for (ForwardIterator i = begin; i != end; i++)
  {
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;

    ForwardIterator k = j;
    k++;
    if (k == end)
      k = begin; 
    
    CompGeom::OrientationType ori = CompGeom::area_sign(*i, *j, *k, tol);
    if (ori == CompGeom::eRight)
      return false;
  }

  // still here?  polygon may be degenerate!
  return true;
}

/// Calculates the convex hull of a set of points in 2D using quickhull
/**
 * \param source_begin an iterator to the beginning of a container of Point2d*
 * \param source_end an iterator pointing to the end of a container of Point2d*
 * \param target_begin an iterator to the beginning of a container of indices;
 *         on return, contains the convex hull (NOTE: size of this container
 *         must be as large as the source container)
 * \return the new end of the target container
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point2d>::calc_convex_hull(ForwardIterator source_begin, ForwardIterator source_end, OutputIterator target_begin)
{
  const unsigned X = 0, Y = 1;
  int exit_code;
  int curlong, totlong;
  char flags[] = "qhull Fx";
  FILE* outfile, * errfile;
  
  FILE_LOG(LOG_COMPGEOM) << "computing 2D convex hull of following points:" << std::endl;
  for (ForwardIterator i = source_begin; i != source_end; i++)
    FILE_LOG(LOG_COMPGEOM) << "  " << *i << std::endl;
 
  // setup constants for qhull
  const int DIM = 2;
  const int N_POINTS = (int) std::distance(source_begin, source_end);
  const boolT IS_MALLOC = false;
  if (N_POINTS <= 2)
    return target_begin;
  
  // setup qhull outputs
  if (LOGGING(LOG_COMPGEOM))
  {
    outfile=stdout;  
    errfile=stderr;
  }
  else
  {
    outfile=NULL;
    errfile=fopen("/dev/null", "w");
    assert(errfile);
  } 

 
  // setup the points
  SAFESTATIC std::vector<coordT> qhull_points;
  qhull_points.resize(N_POINTS*DIM);
  coordT* points_begin = &qhull_points.front();
  std::map<coordT*, Point2d*> vertex_map;
  unsigned j=0;
  for (ForwardIterator i = source_begin; i != source_end; i++)
  {
    qhull_points[j] = (*i)[X];
    qhull_points[j+1] = (*i)[Y];
    vertex_map[points_begin+j] = &(*i);
    j += DIM;
  }

  // lock the qhull mutex -- qhull is non-reentrant
  #ifdef THREADSAFE
  pthread_mutex_lock(&CompGeom::_qhull_mutex);
  #endif  

  // execute qhull  
  exit_code = qh_new_qhull(DIM, N_POINTS, points_begin, IS_MALLOC, flags, outfile, errfile);
  if (exit_code != 0)
  {
    // points are not collinear.. unsure of the error...
    FILE_LOG(LOG_COMPGEOM) << "CompGeom::calc_convex_hull_2D() - unable to execute qhull on points:" << std::endl;
    for (ForwardIterator i = source_begin; i != source_end; i++)
      FILE_LOG(LOG_COMPGEOM) << "  " << *i << std::endl;

    // free qhull memory
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // release the mutex, since we're not using qhull anymore
    #ifdef THREADSAFE
    pthread_mutex_unlock(&CompGeom::_qhull_mutex);
    #endif

    // close the error stream, if necessary
    if (!LOGGING(LOG_COMPGEOM))
      fclose(errfile);

    return target_begin;
  }

  // ordered list of edges
  std::map<Point2d*, std::list<Point2d*> > edges;
  
  // iterate through all facets  
  for (facetT* facet=qh facet_list;facet && facet->next;facet=facet->next)
  {
    // setup a list of vertices for the facet
    std::list<Point2d*> facet_vertices;
    
    // get all vertices in the facet
    vertexT* vertex;
    for (vertexT** vertex_pointer = (vertexT**)& ((facet->vertices)->e[0].p); (vertex = (*vertex_pointer++));)
      facet_vertices.push_back(vertex_map[vertex->point]);
    
    // should be exactly two vertices in the list
    assert(facet_vertices.size() == 2);
    
    // store the edge in the list of edges
    edges[facet_vertices.front()].push_back(facet_vertices.back());
    edges[facet_vertices.back()].push_back(facet_vertices.front());
  }    
  
  // free qhull memory
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);

  // release the qhull mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&CompGeom::_qhull_mutex);
  #endif  

  // construct the set of processed vertex
  std::set<Point2d*> processed;
  
  // construct the hull; compute the area at the same time of the 2D polygon
  Point2d* current_vertex = edges.begin()->first;
  std::list<Point2d> hull;

  while (true)
  {
    // add the current vertex to the list
    hull.push_back(*current_vertex);
    
    // mark this vertex as processed
    processed.insert(current_vertex);
    
    // get adjacent vertices
    std::list<Point2d*>& adj_v = edges[current_vertex];
    
    // see which vertices have been processed
    if (processed.find(adj_v.front()) == processed.end())
      current_vertex = adj_v.front();
    else if (processed.find(adj_v.back()) == processed.end())
      current_vertex = adj_v.back();
    else
      break;    
  }

  // close the error stream, if necessary
  if (!LOGGING(LOG_COMPGEOM))
    fclose(errfile);

  // reverse the hull if necessary
  if (!CompGeom::ccw(hull.begin(), hull.end()))  
    return std::copy(hull.rbegin(), hull.rend(), target_begin);    
  else
    return std::copy(hull.begin(), hull.end(), target_begin);
}

/*****************************************************************************
 Point2d versions of functions END 
 ****************************************************************************/

/*****************************************************************************
 Point2d* versions of functions BEGIN
 ****************************************************************************/

/// Computes the minimum area bounding rectangle of a set of points
/**
 * Uses 2D convex hull and rotating calipers method; runs in O(N lg N) time
 * On return, x1, x2, x3, and x4 are the four vertices of the bounding
 * rectangle (ordered as edges).
 */
template <class ForwardIterator>
void CompGeomSpecOne<ForwardIterator, Point2d*>::calc_min_area_bounding_rect(ForwardIterator begin, ForwardIterator end, Point2d& x1, Point2d& x2, Point2d& x3, Point2d& x4)
{
  enum Flag { F_NONE, F_LEFT, F_RIGHT, F_BOTTOM, F_TOP };
  const unsigned X = 0, Y = 1;
  std::pair<Point2d*, Point2d*> ep;

  // calculate the convex hull of the points in ccw order
  std::vector<Point2d*> points, hull;
  for (; begin != end; begin++)
    points.push_back(&*begin);
  CompGeom::calc_convex_hull(points.begin(), points.end(), std::back_inserter(hull));
  if (hull.empty())
  {
    // convex hull is degenerate; compute line endpoints and make that the 
    // "hull"
    CompGeom::determine_seg_endpoints(begin, end, ep);
    hull.push_back(ep.first);
    hull.push_back(ep.second);
  }
  // get the hull in CCW order
  else if (!CompGeom::ccw(hull.begin(), hull.end()))
    std::reverse(hull.begin(), hull.end());

  // make sure that no three points are colinear
  unsigned n = hull.size();
  for (unsigned i = 0; i < n; i++)
  {
    // get the next two points
    unsigned j = (i < n-1) ? i+1 : 0;
    unsigned k = (j < n-1) ? j+1 : 0;
    if (CompGeom::collinear(*hull[i], *hull[j], *hull[k]))
    {
      // erase the middle point
      hull.erase(hull.begin()+j);

      // decrement both n and i
      i--;
      n--;
      // if n < 3, bounding rectangle is degenerate; however, we can still
      // output something useful
      if (n < 3)
      {
        x1 = *hull[0];
        x2 = *hull[1];
        x3 = *hull[0];
        x4 = *hull[1];
        return;
      }
    }
  }

  // setup unit-length edge directions of the convex polygon
  unsigned nm1 = n - 1;
  std::vector<Ravelin::Vector2d> edges(n);
  std::vector<bool> visited(n, false);
  for (unsigned i=0; i< nm1; i++)
  {
    edges[i] = *hull[i+1] - *hull[i];
    edges[i].normalize();
  }
  edges[nm1] = *hull[0] - *hull[nm1];
  edges[nm1].normalize();

  // find the smallest axis-aligned box containing the points.  Keep track
  // of the extremum indices, L (left), R (right), B (bottom), and T (top)
  // so that the following constraints are met:
  //   V[L].X <= V[i].X for all i and V[(L+1)%N].X > V[L].X
  //   V[L].X >= V[i].X for all i and V[(R+1)%N].X < V[R].X
  //   V[L].Y <= V[i].Y for all i and V[(B+1)%N].Y > V[B].X
  //   V[L].Y >= V[i].Y for all i and V[(T+1)%N].Y < V[T].X
  double xmin = (*hull[0])[X], xmax = xmin;
  double ymin = (*hull[0])[Y], ymax = ymin;
  unsigned Lindex = 0, Rindex = 0, Bindex = 0, Tindex = 0;
  for (unsigned i=1; i< n; i++)
  {
    if ((*hull[i])[X] <= xmin)
    {
      xmin = (*hull[i])[X];
      Lindex = i;
    }
    if ((*hull[i])[X] >= xmax)
    {
      xmax = (*hull[i])[X];
      Rindex = i;
    }
    if ((*hull[i])[Y] <= ymin)
    {
      ymin = (*hull[i])[Y];
      Bindex = i;
    }
    if ((*hull[i])[Y] >= ymax)
    {
      ymax = (*hull[i])[Y];
      Tindex = i;
    }
  }

  // apply wrap-around tests to ensure the constraints mentioned above 
  // are satisfied
  if (Lindex == nm1)
  {
    if ((*hull[0])[X] <= xmin)
    {
      xmin = (*hull[0])[X];
      Lindex = 0;
    }
  }
  if (Rindex == nm1)
  {
    if ((*hull[0])[X] >= xmax)
    {
      xmax = (*hull[0])[X];
      Rindex = 0;
    }
  }
  if (Bindex == nm1)
  {
    if ((*hull[0])[Y] <= ymin)
    {
      ymin = (*hull[0])[Y];
      Bindex = 0;
    }
  }
  if (Tindex == nm1)
  {
    if ((*hull[0])[Y] >= ymax)
    {
      ymax = (*hull[0])[Y];
      Tindex = 0;
    }
  }

  // the dimensions of the axis-aligned box; the extents store width and height
  // for now
  Point2d center((double) 0.5 * (xmin + xmax), (double) 0.5 * (ymin + ymax));
  Ravelin::Vector2d axis[2] = { Ravelin::Vector2d((double) 1.0, (double) 0.0), 
                      Ravelin::Vector2d((double) 0.0, (double) 1.0) };
  double extent[2] = { (double) 0.5 * (xmax - xmin), (double) 0.5 * (ymax - ymin) };
  double min_area_div4 = extent[0]*extent[1];

  // the rotating calipers algorithm follows...
  Ravelin::Vector2d U((double) 1.0, (double) 0.0), V((double) 0.0, (double) 1.0);
  bool done = false;
  while (!done)
  {
    // determine the edge that forms the smallest angle with the current box
    // edges
    Flag flag = F_NONE;
    double maxdot = (double) 0.0;
    double dot = U.dot(edges[Bindex]);
    if (dot > maxdot)
    {
      maxdot = dot;
      flag = F_BOTTOM;
    }
    dot = V.dot(edges[Rindex]);
    if (dot > maxdot)
    {
      maxdot = dot;
      flag = F_RIGHT;
    }
    dot = -U.dot(edges[Tindex]);
    if (dot > maxdot)
    {
      maxdot = dot;
      flag = F_TOP;
    }
    dot = -V.dot(edges[Lindex]);
    if (dot > maxdot)
    {
      maxdot = dot;
      flag = F_LEFT;
    }

    switch (flag)
    {
      case F_BOTTOM:
        if (visited[Bindex])
          done = true;
        else
        {
          // compute box axes with E[B] as an edge
          U = edges[Bindex];
          V = -U.perp();
          CompGeom::update_box(*hull[Lindex], *hull[Rindex], *hull[Bindex], *hull[Tindex], U, V, min_area_div4, center, axis, extent);

          // mark edge visited and rotate the calipers
          visited[Bindex] = true;
          if (++Bindex == n)
             Bindex = 0;
         }
         break;

      case F_RIGHT:
        if (visited[Rindex])
          done = true;
        else
        {
          // compute box axes with E[R] as an edge
          V = edges[Rindex];
          U = V.perp();
          CompGeom::update_box(*hull[Lindex], *hull[Rindex], *hull[Bindex], *hull[Tindex], U, V, min_area_div4, center, axis, extent);

          // mark edge visited and rotate the calipers
          visited[Rindex] = true;
          if (++Rindex == n)
            Rindex = 0;
        }
        break;

      case F_TOP:
        if (visited[Tindex])
          done = true;
        else
        {
          // compute box axes with E[T] as an edge
          U = -edges[Tindex];
          V = -U.perp();
          CompGeom::update_box(*hull[Lindex], *hull[Rindex], *hull[Bindex], *hull[Tindex], U, V, min_area_div4, center, axis, extent);

          // mark edge visited and rotate the calipers
          visited[Tindex] = true;
          if (++Tindex == n)
            Tindex = 0;
        }
        break;

      case F_LEFT:
        if (visited[Lindex])
          done = true;
        else
        {
          // compute box axes with E[L] as an edge
          V = -edges[Lindex];
          U = V.perp();
          CompGeom::update_box(*hull[Lindex], *hull[Rindex], *hull[Bindex], *hull[Tindex], U, V, min_area_div4, center, axis, extent);

          // mark edge visited and rotate the calipers
          visited[Lindex] = true;
          if (++Lindex == n)
            Lindex = 0;
        }
        break;

      case F_NONE:
        // the polygon is a rectangle
        done = true;
        break;
    }
  }

  // convert Eberly's representation to my own
  x1 =  center - axis[X]*extent[X] - axis[Y]*extent[Y];
  x2 =  center + axis[X]*extent[X] - axis[Y]*extent[Y];
  x3 =  center + axis[X]*extent[X] + axis[Y]*extent[Y];
  x4 =  center - axis[X]*extent[X] + axis[Y]*extent[Y];
}

/// Computes the intersection of a polygon and a line segment
/**
 * \param begin iterator pointing to start of collection of Point2d elements
 * \param end iterator pointing to end of collection of Point2d elements
 * \param seg the line segment
 * \param outbegin an iterator to a container of type LineSeg2 that will 
 *        store line segments on/in the polygon on return 
 * \param tol the tolerance for parallel lines
 * \return ePolygonSegInside if the polygon 
 * \note polygon must be given in counter-clockwise order
 * \note first vertex should not appear twice
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point2d*>::intersect_seg_polygon(ForwardIterator begin, ForwardIterator end, const LineSeg2& seg, OutputIterator outbegin)
{
  std::list<double> points;

  // determine whether one (or both) of the endpoints is within the polygon
  if (CompGeom::polygon_location(begin, end, seg.first) == CompGeom::ePolygonInside)
    points.push_back(0.0);
  if (CompGeom::polygon_location(begin, end, seg.second) == CompGeom::ePolygonInside)
    points.push_back(1.0);

  // determine the inverse of the length (squared) of the line segment
  double inv_seg_len_sq = 1.0/(seg.first - seg.second).norm_sq();

  // intersect all line segments
  Point2d isect1, isect2;
  for (ForwardIterator i = begin; i != end; i++)
  {
    // get the next vertex
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;
  
    // intersect the two segments
    switch (CompGeom::intersect_segs(seg, LineSeg2(**i, **j), isect1, isect2))
    {
      case CompGeom::eSegSegNoIntersect:
        break;

      case CompGeom::eSegSegIntersect:
      case CompGeom::eSegSegVertex: 
        points.push_back((isect1 - seg.second).norm_sq() * inv_seg_len_sq);
        break;

      case CompGeom::eSegSegEdge:
        points.push_back((isect1 - seg.second).norm_sq() * inv_seg_len_sq);
        points.push_back((isect2 - seg.second).norm_sq() * inv_seg_len_sq);
        break;
    }
  }

  // sort the points
  points.sort();

  // make segments out of the points
  for (std::list<double>::const_iterator i = points.begin(); i != points.rbegin().base(); i++)
  {
    std::list<double>::const_iterator j = i;
    j++;
    Point2d p1 = seg.first * (*i) + seg.second * (1.0 - (*i));
    Point2d p2 = seg.first * (*j) + seg.second * (1.0 - (*j));
    *outbegin++ = LineSeg2(p1, p2);
  }

  return outbegin;
}

/// Determines the endpoints for a container of collinear Point2d* objects
/**
 * \param begin iterator to beginning of container of type Point2d*
 * \param end iterator to end of container of type Point2d*
 * \param endpoints the two farthest points on the segment on return
 */
template <class ForwardIterator>
void CompGeom::determine_seg_endpoints(ForwardIterator begin, ForwardIterator end, std::pair<Point2d*, Point2d*>& endpoints)
{
  // make sure that we have been given valid input
  assert(begin != end);

  // setup initial endpoints
  endpoints.first = *begin;
  endpoints.second = *begin;
  double dist = 0; 

  for (ForwardIterator i = ++begin; i != end; i++)
  {
    // get distance from i to both current bounding points
    double dist_e = (**i - *endpoints.second).norm();
    double dist_s = (**i - *endpoints.first).norm();

    // see which distance would be greatest
    if (dist > dist_e)
    {
      // check for no change
      if (dist > dist_s)
        continue;
      else
      {
        dist = dist_s;
        endpoints.second = *i;
      }
    }
    else
    {
      if (dist_e > dist_s)
      {
        dist = dist_e;
        endpoints.first = *i;
      }
      else
      {
        dist = dist_s;
        endpoints.second = *i;
      }
    }
  }
}

/// Computes the intersection of a convex polygon and a line segment
/**
 * \param begin iterator pointing to start of collection of Point2d* elements
 * \param end iterator pointing to end of collection of Point2d* elements
 * \param seg the line segment
 * \param te the parameter of the line segment for the beginning of the
 *        intersection; (1-te)*seg.first + te*seg.second is point of 
 *        intersection
 * \param tl the parameter of the line segment for the end of the intersection;
 *        (1-tl)*seg.first + tl*seg.second is point of intersection
 * \param tol the tolerance for parallel lines
 * \return <b>true</b> if the two intersect and <b>false</b> otherwise
 * \note polygon must be given in counter-clockwise order
 * \note first vertex should not appear twice
 * \note taken from http://geometryalgorithms.com 
 */
template <class ForwardIterator>
bool CompGeomSpecOne<ForwardIterator, Point2d*>::intersect_seg_convex_polygon(ForwardIterator begin, ForwardIterator end, const LineSeg2& seg, double& te, double& tl, double tol)
{
  assert(tol >= 0.0);

  // initialize te and tl
  te = 0;
  tl = 1;

  Ravelin::Vector2d dS = seg.second - seg.first;

  // iterate over all vertices
  for (ForwardIterator i = begin; i != end; i++)
  {
    // get the next vertex
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;

    // get the edge
    Ravelin::Vector2d edge = **j - **i;

    // determine the outward normal of the edge
    Ravelin::Vector2d ni(edge[1], -edge[0]);
    double N = -Ravelin::Vector2d::dot(ni, seg.first-**i);
    double D = Ravelin::Vector2d::dot(dS, ni);
    if (std::fabs(D) < 0.0)
    {
      // segment is parallel to this edge
      if (N < tol)
        // first point is outside of the edge; segment cannot intersect poly
        return false;
      else
        // segment cannot enter or leave poly across this edge, process next
        continue;
    }

    double t = N / D;
    if (D < 0.0)
    {
      // segment enters polygon across this edge
      te = std::max(te, t);
      if (te > tl)
        // segment enters polygon after leaving
        return false;
    }
    else
    {
      assert(D > 0.0);
      tl = std::min(tl, t);
      if (tl < te)
        // segment leaves polygon before entering
        return false;
    }  
  }

  return true;
}

/**
 * Converts a collection of Ravelin::Vector2d* objects to Ravelin::Vector3d objects
 * \param source_begin an iterator pointing to the beginning of the Ravelin::Vector3d 
 *        objects
 * \param source_end an iterator pointing to the end of the Ravelin::Vector3d objects
 * \param begin_target an iterator pointing to the beginning of the Ravelin::Vector2d
 *        objects
 * \param R the projection matrix from 2D to 3D
 * \param offset the offset that must be added to the Z-coordinate of points
 *        projected back from 2D to 3D
 * \return the end of the output range
 * \note the size of the target collection must be equal to the size of the
 *       source collection
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point2d*>::to_3D(ForwardIterator begin_source, ForwardIterator end_source, OutputIterator begin_target, const Ravelin::Matrix3d& RT, double offset)
{
  // project the points back to 3D
  for (ForwardIterator i = begin_source; i != end_source; i++, begin_target++)
    *begin_target = CompGeom::to_3D(**i, RT, offset);

  return begin_target;
}

/**
 * Determines whether a polygon in 2D is counter-clockwise
 * \note Degenerate polygons (alternating representation) will fail!
 */
template <class ForwardIterator>
bool CompGeomSpecOne<ForwardIterator, Point2d*>::ccw(ForwardIterator begin, ForwardIterator end, double tol)
{
  assert(tol >= 0.0);

  for (ForwardIterator i = begin; i != end; i++)
  {
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;

    ForwardIterator k = j;
    k++;
    if (k == end)
      k = begin; 
    
    CompGeom::OrientationType ori = CompGeom::area_sign(**i, **j, **k, tol);
    if (ori == CompGeom::eRight)
      return false;
  }

  // still here?  polygon may be degenerate!
  return true;
}

/// Determines whether a polygon in 2D is convex
template <class ForwardIterator>
bool CompGeomSpecOne<ForwardIterator, Point2d*>::is_convex_polygon(ForwardIterator begin, ForwardIterator end, double tol)
{
  assert(tol >= 0.0);

  // check whether every third point is to the left of the line segment
  // formed by the two preceeding points
  for (ForwardIterator i = begin; i != end; i++)
  {
    // get the next point -- if we've gone to the end, recycle
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;

    // get the following point -- if we've gone past the end, recycle
    ForwardIterator k = j;  
    k++;
    if (k == end)
      k = begin;

    // verify that k is not to the right of j and k
    if (CompGeom::area_sign(**i, **j, **k, tol) == CompGeom::eRight)
      return false;
  }

  // all checks passed.. convex
  return true;
}

/// Calculates the convex hull of a set of points in 2D using quickhull
/**
 * \param source_begin an iterator to the beginning of a container of Point2d*
 * \param source_end an iterator pointing to the end of a container of Point2d*
 * \param target_begin an iterator to the beginning of a container of indices;
 *         on return, contains the convex hull (NOTE: size of this container
 *         must be as large as the source container)
 * \return the new end of the target container
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeomSpecTwo<ForwardIterator, OutputIterator, Point2d*>::calc_convex_hull(ForwardIterator source_begin, ForwardIterator source_end, OutputIterator target_begin)
{
  const unsigned X = 0, Y = 1;
  int exit_code;
  int curlong, totlong;
  char flags[] = "qhull Fx";
  FILE* outfile, * errfile;
  
  FILE_LOG(LOG_COMPGEOM) << "computing 2D convex hull of following points:" << std::endl;
  for (ForwardIterator i = source_begin; i != source_end; i++)
    FILE_LOG(LOG_COMPGEOM) << "  " << *i << std::endl;
  
  // setup qhull outputs
  if (LOGGING(LOG_COMPGEOM))
  {
    outfile=stdout;  
    errfile=stderr;
  }
  else
  {
    outfile=NULL;
    errfile=fopen("/dev/null", "w");
    assert(errfile);
  } 

  // setup constants for qhull
  const int DIM = 2;
  const int N_POINTS = (int) std::distance(source_begin, source_end);
  const boolT IS_MALLOC = false;
  assert(N_POINTS > 2);
  
  // setup the points
  std::map<coordT*, Point2d*> vertex_map;
  SAFESTATIC std::vector<coordT> qhull_points;
  qhull_points.resize(N_POINTS*DIM);
  coordT* points_begin = &qhull_points.front();
  unsigned j=0;
  for (ForwardIterator i = source_begin; i != source_end; i++)
  {
    qhull_points[j] = (**i)[X];
    qhull_points[j+1] = (**i)[Y];
    vertex_map[points_begin+j] = *i;
    j += DIM;
  }

  // lock the qhull mutex -- qhull is non-reentrant
  #ifdef THREADSAFE
  pthread_mutex_lock(&CompGeom::_qhull_mutex);
  #endif  

  // execute qhull  
  exit_code = qh_new_qhull(DIM, N_POINTS, points_begin, IS_MALLOC, flags, outfile, errfile);
  if (exit_code != 0)
  {
    // points are not collinear.. unsure of the error...
    FILE_LOG(LOG_COMPGEOM) << "CompGeom::calc_convex_hull_2D() - unable to execute qhull on points:" << std::endl;
    for (ForwardIterator i = source_begin; i != source_end; i++)
      FILE_LOG(LOG_COMPGEOM) << "  " << *i << std::endl;

    // free qhull memory
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // release the mutex, since we're not using qhull anymore
    #ifdef THREADSAFE
    pthread_mutex_unlock(&CompGeom::_qhull_mutex);
    #endif

    // close the error stream, if necessary
    if (!LOGGING(LOG_COMPGEOM))
      fclose(errfile);

    throw NumericalException(); 
    return target_begin;
  }

  // ordered list of edges
  std::map<Point2d*, std::list<Point2d*> > edges;
  
  // iterate through all facets  
  for (facetT* facet=qh facet_list;facet && facet->next;facet=facet->next)
  {
    // setup a list of vertices for the facet
    std::list<Point2d*> facet_vertices;
    
    // get all vertices in the facet
    vertexT* vertex;
    for (vertexT** vertex_pointer = (vertexT**)& ((facet->vertices)->e[0].p); (vertex = (*vertex_pointer++));)
      facet_vertices.push_back(vertex_map[vertex->point]);
    
    // should be exactly two vertices in the list
    assert(facet_vertices.size() == 2);
    
    // store the edge in the list of edges
    edges[facet_vertices.front()].push_back(facet_vertices.back());
    edges[facet_vertices.back()].push_back(facet_vertices.front());
  }    
  
  // free qhull memory
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);

  // release the qhull mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&CompGeom::_qhull_mutex);
  #endif  

  // construct the set of processed vertex
  std::set<Point2d*> processed;
  
  // construct the hull; compute the area at the same time of the 2D polygon
  Point2d* current_vertex = edges.begin()->first;
  std::list<Point2d*> hull;

  while (true)
  {
    // add the current vertex to the list
    hull.push_back(current_vertex);
    
    // mark this vertex as processed
    processed.insert(current_vertex);
    
    // get adjacent vertices
    std::list<Point2d*>& adj_v = edges[current_vertex];
    
    // see which vertices have been processed
    if (processed.find(adj_v.front()) == processed.end())
      current_vertex = adj_v.front();
    else if (processed.find(adj_v.back()) == processed.end())
      current_vertex = adj_v.back();
    else
      break;    
  }

  // close the error stream, if necessary
  if (!LOGGING(LOG_COMPGEOM))
    fclose(errfile);

  // reverse the hull if necessary
  if (!CompGeom::ccw(hull.begin(), hull.end()))  
    return std::copy(hull.rbegin(), hull.rend(), target_begin);    
  else
    return std::copy(hull.begin(), hull.end(), target_begin);
}

/*****************************************************************************
 Point2d* versions of functions END 
 ****************************************************************************/

/// Computes the intersection of a polygon and a line segment
/**
 * \param begin iterator pointing to start of collection of Point2d elements
 * \param end iterator pointing to end of collection of Point2d elements
 * \param seg the line segment
 * \param outbegin an iterator to a container of type LineSeg2 that will 
 *        store line segments on/in the polygon on return 
 * \param tol the tolerance for parallel lines
 * \return ePolygonSegInside if the polygon 
 * \note polygon must be given in counter-clockwise order
 * \note first vertex should not appear twice
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeom::intersect_seg_polygon(ForwardIterator begin, ForwardIterator end, const LineSeg2& seg, OutputIterator outbegin)
{
  return CompGeomSpecTwo<ForwardIterator, OutputIterator, typename std::iterator_traits<ForwardIterator>::value_type>::intersect_seg_polygon(begin, end, seg, outbegin);
}

/// Computes the 3D convex hull of a set of points
/**
 * \param first a forward iterator for type Ravelin::Vector3
 * \param last a forward iterator for type Ravelin::Vector3
 * \return a pointer to the newly created polyhedron
 */
template <class ForwardIterator>
PolyhedronPtr CompGeom::calc_convex_hull(ForwardIterator first, ForwardIterator last)
{
  return CompGeomSpecOne<ForwardIterator, Point3d>::calc_convex_hull(first, last);
}

/// Finds an interior point of a set of halfspaces using linear programming
/**
 * The method used to find the interior point is of order O(n), where n is the
 * number of halfspace constraints.
 * \param start an iterator to the start of a collection of halfspaces, each of type std::pair<Ravelin::Vector3, double>, where the vector [nx ny nz] is the normal to the halfspace and the scalar is the offset 'd'; each halfspace will satisfy the equation nx*x + ny*y + nz*z <= d
 * \param end an iterator to the end of the collection
 * \param point contains the interior point on return if successful
 * \return the minimum distance from a halfspace of the interior point; if
 *         the point is negative, there is no interior point
 */
/*
template <class ForwardIterator>
double CompGeom::find_hs_interior_point(ForwardIterator start, ForwardIterator end, Point3d& point)
{
  assert(std::numeric_limits<double>::has_infinity);
  const double inf = std::numeric_limits<double>::max();
  const unsigned D = 5;
  const unsigned N = distance(start, end);

  // setup the limits on the variables
  Ravelin::VectorNd l(D), u(D);
  l[0] = -1.0;        u[0] = 1.0;
  l[1] = -1.0;        u[1] = 1.0;
  l[2] = -1.0;        u[2] = 1.0;
  l[3] = 0;            u[3] = inf;
  l[4] = 0;            u[4] = 1.0;

  // setup the optimization vector
  Ravelin::VectorNd c = Ravelin::VectorNd::zero(D);
  c[D-1] = 1.0;

  // setup b
  Ravelin::VectorNd b = Ravelin::VectorNd::zero(N);

  // setup A
  Ravelin::MatrixNd A(N,D);
  unsigned i = 0;
  for (; start != end; start++, i++)
  {
    A(i,0) = start->first[0];
    A(i,1) = start->first[1];
    A(i,2) = start->first[2];
    A(i,3) = -start->second;
    A(i,4) = 1.0;
  }

  // do linear programming
  Ravelin::VectorNd x;
  if (!Optimization::lp(A, b, c, l, u, x))
    return -1.0;

  // verify that x[3] is not zero
  if (x[3] <= std::numeric_limits<double>::epsilon())
    return -1.0;

  // determine interior point
  point = Ravelin::Vector3d(x[0]/x[3], x[1]/x[3], x[2]/x[3]);

  // return the distance
  return x[4]/x[3];
}
*/

/// Computes the halfspace intersection, returning the result as a convex polyhedron
/**
 * \param start an iterator to the start of a collection of halfspaces, each of type std::pair<Ravelin::Vector3, double>, where the vector [nx ny nz] is the normal to the halfspace and the scalar is the offset 'd'; each halfspace will satisfy the equation nx*x + ny*y + nz*z <= d
 * \param end an iterator to the end of the collection
 * \param interior_point a point interior to all of the halfspaces
 * \return a pointer to the created polyhedron or a NULL pointer if unsuccessful
 */
template <class ForwardIterator>
PolyhedronPtr CompGeom::calc_hs_intersection(ForwardIterator start, ForwardIterator end, const Ravelin::VectorNd& interior_point)
{
  const int DIM = 4;
  const boolT IS_MALLOC = false;
  int curlong, totlong;
  const unsigned X = 0, Y = 1, Z = 2;
  FILE* outfile, * errfile;

  // setup qhull flags
  std::ostringstream flags;
  flags << "qhull H";
  flags << interior_point[X] << "," << interior_point[Y] << "," << interior_point[Z];

  // setup qhull outputs
  if (LOGGING(LOG_COMPGEOM))
  {
    outfile=stdout;  
    errfile=stderr;
  }
  else
  {
    outfile=NULL;
    errfile=fopen("/dev/null", "w");
    assert(errfile);
  } 

  FILE_LOG(LOG_COMPGEOM) << "computing halfspace intersection of: " << std::endl;
  for (ForwardIterator i = start; i != end; i++)
    FILE_LOG(LOG_COMPGEOM) << "  halfspace normal: " << i->first << "  d: " << i->second << std::endl;

  // allocate memory for halfspaces
  int nspaces = (int) distance(start, end);
  boost::shared_array<coordT> qhull_hs(new coordT[nspaces*DIM]);

  // setup halfspaces
  unsigned j=0;
  for (ForwardIterator i = start; i != end; i++)
  {
    qhull_hs[j] = i->first[X];
    qhull_hs[j+1] = i->first[Y];
    qhull_hs[j+2] = i->first[Z];
    qhull_hs[j+3] = -i->second;
    j += DIM;
  }

  // lock the qhull mutex -- qhull is non-reentrant
  #ifdef THREADSAFE
  pthread_mutex_lock(&CompGeom::_qhull_mutex);
  #endif

  // execute qhull
  int exit_code = qh_new_qhull(DIM, nspaces, qhull_hs.get(), IS_MALLOC, (char*) flags.str().c_str(), outfile, errfile);
  if (exit_code)
  {
    // free qhull memory
    qh_freeqhull(!qh_ALL);
    qh_memfreeshort(&curlong, &totlong);

    // qhull failed
    #ifdef THREADSAFE
    pthread_mutex_unlock(&CompGeom::_qhull_mutex);
    #endif

    // close the error stream, if necessary
    if (!LOGGING(LOG_COMPGEOM))
      fclose(errfile);

    throw NumericalException(); 
    return PolyhedronPtr();
  }

  // verify that the qhull dimension is correct
  assert(qh hull_dim == 3);

  // determine the intersection vertices; NOTE: this code was motivated by
  // qhull's qh_printafacet() and qh_printfacets() functions
  std::list<Point3d> points;
  for (facetT* facet=qh facet_list;facet && facet->next;facet=facet->next)
  {
    if (facet->offset > 0)
    {
      // facet has infinite offset
      return PolyhedronPtr();
    }
    coordT* point = (coordT*) qh_memalloc(qh normal_size);
    coordT* coordp = point;
    coordT* normp = facet->normal;
    coordT* feasiblep = qh feasible_point;
    if (facet->offset < -qh MINdenom)
      for (int k=qh hull_dim; k--; )
        *(coordp++) = (*(normp++) / -facet->offset) + *(feasiblep++);
    else
      for (int k=qh hull_dim; k--;)
      {
        boolT zerodiv;
        coordT* feasiblep1 = feasiblep+1;
        feasiblep = feasiblep1+1;
        *(coordp++) = qh_divzero(*(normp++), facet->offset, qh MINdenom_1, &zerodiv) + *feasiblep1  + *feasiblep;
        feasiblep++;
        if (zerodiv)
        {
          // facet has infinite offset
          qh_memfree(point, qh normal_size);
          return PolyhedronPtr();
        }
      }

    // add the point
    points.push_back(Ravelin::Vector3d(*point, *(point+1), *(point+2)));

    // free the temporary memory
    qh_memfree(point, qh normal_size);
  }

  // free qhull memory
  qh_freeqhull(!qh_ALL);
  qh_memfreeshort(&curlong, &totlong);
  assert(!curlong && !totlong);
  
  // release the qhull mutex
  #ifdef THREADSAFE
  pthread_mutex_unlock(&CompGeom::_qhull_mutex);
  #endif

  // now, calculate the convex hull of the intersection points  
  PolyhedronPtr p = calc_convex_hull(points.begin(), points.end());

  // close the error stream, if necessary
  if (!LOGGING(LOG_COMPGEOM))
    fclose(errfile);

  return p;
}

/// Computes the intersection of a convex polygon and a line segment
/**
 * \param begin iterator pointing to start of collection of Point2d elements
 * \param end iterator pointing to end of collection of Point2d elements
 * \param seg the line segment
 * \param te the parameter of the line segment for the beginning of the
 *        intersection; (1-te)*seg.first + te*seg.second is point of 
 *        intersection
 * \param tl the parameter of the line segment for the end of the intersection;
 *        (1-tl)*seg.first + tl*seg.second is point of intersection
 * \param tol the tolerance for parallel lines
 * \return <b>true</b> if the two intersect and <b>false</b> otherwise
 * \note polygon must be given in counter-clockwise order
 * \note first vertex should not appear twice
 * \note taken from http://geometryalgorithms.com 
 */
template <class ForwardIterator>
bool CompGeom::intersect_seg_convex_polygon(ForwardIterator begin, ForwardIterator end, const LineSeg2& seg, double& te, double& tl, double tol)
{
  return CompGeomSpecOne<ForwardIterator, typename std::iterator_traits<ForwardIterator>::value_type
  >::intersect_seg_convex_polygon(begin, end, seg, te, tl, tol);
}

/**
 * Converts a collection of Ravelin::Vector3d objects to Point2d objects
 * \param source_begin an iterator pointing to the beginning of the Ravelin::Vector3d 
 *        objects
 * \param source_end an iterator pointing to the end of the Ravelin::Vector3d objects
 * \param begin_target an iterator pointing to the beginning of the Point2d
 *        objects
 * \param R the projection matrix from 3D to 2D (on return)
 * \return the end of the output range
 * \note the size of the target collection must be equal to the size of the
 *       source collection
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeom::to_2D(ForwardIterator begin_source, ForwardIterator end_source, OutputIterator begin_target, const Ravelin::Matrix3d& R)
{
  return CompGeomSpecTwo<ForwardIterator, OutputIterator, typename std::iterator_traits<ForwardIterator>::value_type>::to_2D(begin_source, end_source, begin_target, R);
}

/**
 * Converts a collection of Point2d objects to Ravelin::Vector3d objects
 * \param source_begin an iterator pointing to the beginning of the Ravelin::Vector3d 
 *        objects
 * \param source_end an iterator pointing to the end of the Ravelin::Vector3d objects
 * \param begin_target an iterator pointing to the beginning of the Point2d
 *        objects
 * \param R the projection matrix from 2D to 3D
 * \param offset the offset that must be added to the Z-coordinate of points
 *        projected back from 2D to 3D
 * \return the end of the output range
 * \note the size of the target collection must be equal to the size of the
 *       source collection
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeom::to_3D(ForwardIterator begin_source, ForwardIterator end_source, OutputIterator begin_target, const Ravelin::Matrix3d& RT, double offset)
{
  return CompGeomSpecTwo<ForwardIterator, OutputIterator, typename std::iterator_traits<ForwardIterator>::value_type>::to_3D(begin_source, end_source, begin_target, RT, offset);
}

/**
 * Determines whether a polygon in 2D is counter-clockwise
 * \note Degenerate polygons (alternating representation) will fail!
 */
template <class ForwardIterator>
bool CompGeom::ccw(ForwardIterator begin, ForwardIterator end, double tol)
{
  return CompGeomSpecOne<ForwardIterator, typename std::iterator_traits<ForwardIterator>::value_type>::ccw(begin, end, tol);
}

/**
 * Determines whether a polygon in 2D is counter-clockwise
 * \note Degenerate polygons (alternating representation) will fail!
 */
template <class ForwardIterator>
bool CompGeom::ccw(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double tol)
{
  return CompGeomSpecOne<ForwardIterator, typename std::iterator_traits<ForwardIterator>::value_type>::ccw(begin, end, normal, tol);
}

/// Intersects two coplanar triangles
/**
 * \param normal a normal to both triangles
 * \param begin an iterator to a container that can hold the points of 
 *         intersection (i.e., a container of size 6 or greater); on return,
 *         the container will contain a ccw polygon (in 3D) [orientation with
 *         respect to the given normal]
 * \return the end of the container output
 */
template <class OutputIterator>
OutputIterator CompGeom::intersect_coplanar_tris(const Triangle& t1, const Triangle& t2, const Ravelin::Vector3d& normal, OutputIterator begin)
{
  const unsigned TRI_VERTS = 3;

  // verify all points are defined with respect to the same pose
  #ifndef NDEBUG
  for (unsigned i=0; i< TRI_VERTS; i++)
    assert(t1.get_vertex(i).pose == t2.get_vertex(i).pose);
  for (unsigned i=1; i< TRI_VERTS; i++)
    assert(t1.get_vertex(i).pose == t1.get_vertex(i-1).pose);
  #endif

  // project triangles to 2D
  Ravelin::Matrix3d R = calc_3D_to_2D_matrix(normal);
  double offset = determine_3D_to_2D_offset(Ravelin::Origin3d(t1.a), R);
  Point2d t1_2D[TRI_VERTS], t2_2D[TRI_VERTS];
  for (unsigned i=0; i< TRI_VERTS; i++)
  {
    t1_2D[i] = Point2d(to_2D(t1.get_vertex(i), R), boost::shared_ptr<Ravelin::Pose2d>());
    t2_2D[i] = Point2d(to_2D(t2.get_vertex(i), R), boost::shared_ptr<Ravelin::Pose2d>());
  } 

  // verify triangles are ccw and reverse if necessary
  if (!ccw(t1_2D, t1_2D + TRI_VERTS))
    std::swap(t1_2D[1], t1_2D[2]);
  if (!ccw(t2_2D, t2_2D + TRI_VERTS))
    std::swap(t2_2D[1], t2_2D[2]);

  // intersect triangles
  std::list<Point2d> points;
  intersect_tris(t1_2D, t2_2D, std::back_inserter(points));

  // project points back to 3D
  Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);
  BOOST_FOREACH(const Point2d& p, points)
  {
    Ravelin::Origin2d o(p);
    *begin++ = Point3d(to_3D(o, RT, offset), t1.a.pose);
  }

  return begin;
}

/// Intersects two polygons in 3D
/**
 * \param pbegin an iterator pointing to the container holding
 *        a ccw polygon
 * \param pend the end of the container holding the first polygon
 * \param qbegin an iterator pointing to the container holding
 *        a ccw polygon
 * \param qend the end of the container holding the second polygon
 * \param normal the normal to the polygons
 * \param isect_begin on return, the polygon of intersection will be placed
 *         here with ccw orientation; this container must be big enough to hold
 *         the result (i.e., it must be at least of size min(p,q)
 * \return the iterator pointing to the end of the container of intersection
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeom::intersect_polygons(ForwardIterator pbegin, ForwardIterator pend, ForwardIterator qbegin, ForwardIterator qend, const Ravelin::Vector3d& normal, OutputIterator isect_begin)
{
  // **************************************************************
  // first, we need to project the 3D triangles to 2D polygons
  // **************************************************************

  // verify all points are in the same frame
  #ifndef NDEBUG
  for (ForwardIterator i = pbegin; i != pend; i++)
  {
    ForwardIterator j = i;
    j++;
    if (j == pend)
      continue;
    assert(i->pose == j->pose);
  } 
  for (ForwardIterator i = qbegin; i != qend; i++)
  {
    ForwardIterator j = i;
    j++;
    if (j == qend)
      continue;
    assert(i->pose == j->pose);
  }
  assert(pbegin->pose == qbegin->pose);   
  #endif

  // R will project the points such that they lie in the plane z=0
  Ravelin::Matrix3d R = calc_3D_to_2D_matrix(normal);
  double offset = determine_3D_to_2D_offset(Ravelin::Origin3d(*pbegin), R);
  
  // convert the two polygons to 2D
  std::vector<Point2d> p, q;
  for (ForwardIterator i = pbegin; i != pend; i++)
    p.push_back(Point2d(to_2D(*i, R), boost::shared_ptr<Ravelin::Pose2d>()));
  for (ForwardIterator i = qbegin; i != qend; i++)
    q.push_back(Point2d(to_2D(*i, R), boost::shared_ptr<Ravelin::Pose2d>()));

  // do the intersection
  std::list<Point2d> isect_2D;
  std::insert_iterator<std::list<Point2d> > ii(isect_2D, isect_2D.begin());
  intersect_polygons(p.begin(), p.end(), q.begin(), q.end(), ii);

  // transform the polygon of intersection to 3D
  R.transpose();
  std::list<Point3d> polygon;
  for (std::list<Point2d>::const_iterator i = isect_2D.begin(); i != isect_2D.end(); i++)
  {
    Ravelin::Origin2d o(*i);
    polygon.push_back(Point3d(to_3D(o, R, offset), pbegin->pose));
  }

  // verify that the polygon is ccw; if not, make it so
  for (std::list<Point3d>::const_iterator i = polygon.begin(); i != polygon.end(); i++)
  {
    // get the next two points
    std::list<Point3d>::const_iterator j = i;
    j++;
    if (j == polygon.end())
      j = polygon.begin();
    std::list<Point3d>::const_iterator k = j;
    k++;
    if (k == polygon.end())
      k = polygon.begin();

    // check the cross product
    Ravelin::Vector3d cprod = Ravelin::Vector3d::cross((*j - *i), (*k - *j));
    
    // make sure that the cross product is not zero
    if (Ravelin::Vector3d::norm(cprod) < NEAR_ZERO)
      continue;

    // check the dot product of the cross product and the normal
    double dp = Ravelin::Vector3d::dot(cprod, normal);

    // if the dot product of the cross product and the normal is less than zero, then the polygon is
    // oriented incorrectly
    if (dp > NEAR_ZERO)
      break;
    else if (dp < -NEAR_ZERO)
    {
      polygon.reverse();
      break;
    }
    else
    {
      assert(false);
    }
  }

  // copy the polygon to the output
  return std::copy(polygon.begin(), polygon.end(), isect_begin);
}

/// Intersects two polygons in 2D
/**
 * \param pbegin a random access iterator pointing to the container holding
 *        a ccw polygon (of Point2d)
 * \param pend the end of the container holding the first polygon (of Point2d)
 * \param qbegin a random access iterator pointing to the container holding
 *        a ccw polygon (of Point2d)
 * \param qend the end of the container holding the second polygon (of Point2d)
 * \param isect_begin on return, the polygon of intersection will be placed
 *         here with ccw orientation; this container must be big enough to hold
 *         the result (i.e., it must be at least of size min(p,q)
 * \return the iterator pointing to the end of the container of intersection
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeom::intersect_polygons(ForwardIterator pbegin, ForwardIterator pend, ForwardIterator qbegin, ForwardIterator qend, OutputIterator isect_begin)
{
  enum tInFlag { Pin, Qin, Unknown };

  // verify that both polygons are ccw
  assert(ccw(pbegin, pend));
  assert(ccw(qbegin, qend));  

  // get how many points in p and q
  unsigned np = std::distance(pbegin, pend);
  unsigned nq = std::distance(qbegin, qend);

  // now compute their intersections  
  unsigned a = 0, b = 0, aa = 0, ba = 0;
  tInFlag inflag = Unknown;
  bool first_point = true;
  Point2d origin(0,0);
  OutputIterator current = isect_begin;

  do
  {
    unsigned a1 = (a + np-1) % np;
    unsigned b1 = (b + nq-1) % nq;

    Ravelin::Vector2d AX = pbegin[a] - pbegin[a1];
    Ravelin::Vector2d BX = qbegin[b] - qbegin[b1];

    // determine signs of cross-products
    OrientationType cross = area_sign(origin, AX, BX);
    OrientationType aHB = area_sign(qbegin[b1], qbegin[b], pbegin[a]);
    OrientationType bHA = area_sign(pbegin[a1], pbegin[a], qbegin[b]);
    
    // if A and B intersect, update inflag
    Point2d p, q;
    SegSegIntersectType code = intersect_segs(std::make_pair(pbegin[a1], pbegin[a]), std::make_pair(qbegin[b1], qbegin[b]), p, q);
    if (code == eSegSegVertex || code == eSegSegIntersect)
    {
      if (inflag == Unknown && first_point)
      {
        aa = ba = 0;
        first_point = false;
      }

      *current++ = p;
      if (aHB == eLeft)
        inflag = Pin;
      else if (bHA == eLeft)
        inflag = Qin;
    }

    // --------- Advance rules --------------
    // special cases: O'Rourke p. 262
    // special case: A and B overlap and oppositely oriented; intersection
    // is only a line segment
    if ((code == eSegSegEdge) && Ravelin::Vector2d::dot(AX,BX) < 0)
    {
      *current++ = p;
      *current++ = q;
      return current;
    }
    // special case: A and B are parallel and disjoint
    else if ((cross == eOn) && (aHB == eRight) && (bHA == eRight))
      return isect_begin;
    // special case: A and B are collinear
    else if ((cross == eOn) && (aHB == eOn) && (bHA == eOn))
    {
      // advance but do not add point to intersecting polygon
      if (inflag == Pin)
        b = advance(b, &ba, nq, inflag == Qin, qbegin[b], current);
      else
        a = advance(a, &aa, np, inflag == Pin, pbegin[a], current);
    }
    // generic cases (continued from p. 258)
    else if (cross == eOn || cross == eLeft)
    {
      if (bHA == eLeft)
        a = advance(a, &aa, np, inflag == Pin, pbegin[a], current);
      else
        b = advance(b, &ba, nq, inflag == Qin, qbegin[b], current);
    }
    else
    {
      if (aHB == eLeft)
        b = advance(b, &ba, nq, inflag == Qin, qbegin[b], current);
      else
        a = advance(a, &aa, np, inflag == Pin, pbegin[a], current);
    }
  }
  while (((aa < np) || (ba < nq)) && (aa < 2*np) && (ba < 2*nq));

  // deal with remaining special cases: not implemented
  if (inflag == Unknown)
    return isect_begin;
  
  return current;
}

/// Utility function for intersect_coplanar_tris()
/**
 * Taken from O'Rourke, p. 259.
 */
template <class OutputIterator>
unsigned CompGeom::advance(unsigned a, unsigned* aa, unsigned n, bool inside, const Point2d& p, OutputIterator& current)
{
  if (inside)
    *current = p;

  (*aa)++;
  return (a+1) % n;
}

/// Intersects two triangles in 3D and returns the points of intersection
/**
 * \param begin an iterator to a container that can hold the points of 
 *         intersection (i.e., a container of size 6 or greater of Ravelin::Vector3
 *         objects); on return, the container will contain a ccw polygon (in 3D) 
 *         [orientation with respect to the given normal]
 * \return the end of the container output
 */
template <class OutputIterator>
OutputIterator CompGeom::intersect_tris(const Triangle& t1, const Triangle& t2, OutputIterator begin)
{
  // determine whether the triangles are coplanar
  if (CompGeom::coplanar(t1, t2))
    return intersect_coplanar_tris(t1, t2, t1.calc_normal(), begin);

  // intersect the triangles
  Point3d p1, p2;
  if (!intersect_noncoplanar_tris(t1, t2, p1, p2))
    // no intersection, clear the list and return
    return begin;
  else
  {
    *begin++ = p1;
    *begin++ = p2;
    return begin;
  }
}

/// Calculates the convex hull of a set of points in 2D using quickhull
/**
 * \param source_begin an iterator to the beginning of a container of Point2d 
 *        or Point2d*
 * \param source_end an iterator pointing to the end of a container of Point2d
 *        or Point2d*
 * \param target_begin an iterator to the beginning of a container of indices;
 *         on return, contains the convex hull (NOTE: size of this container
 *         must be as large as the source container)
 * \return the new end of the target container
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeom::calc_convex_hull(ForwardIterator source_begin, ForwardIterator source_end, OutputIterator target_begin)
{
  return CompGeomSpecTwo<ForwardIterator, OutputIterator, typename std::iterator_traits<ForwardIterator>::value_type>::calc_convex_hull(source_begin, source_end, target_begin);
}

/// Calculates the convex hull of a set of points that lie on a 2D manifold using quickhull
/**
 * \param source_begin an iterator to the beginning of a container of points
 * \param source_end an iterator pointing to the end of a container of points
 * \param normal the (optional) normal of the points; this will be computed if normal is zero vector
 * \param target_begin an iterator to the beginning of a container of points;
 *         on return, contains the convex hull (NOTE: size of this container
 *         must be as large as the source container)
 * \return the new end of the target container
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeom::calc_convex_hull(ForwardIterator source_begin, ForwardIterator source_end, const Ravelin::Vector3d& normal, OutputIterator target_begin)
{
  return CompGeomSpecTwo<ForwardIterator, OutputIterator, typename std::iterator_traits<ForwardIterator>::value_type>::calc_convex_hull(source_begin, source_end, normal, target_begin);
}

/// Determines whether a polygon in 2D is convex
template <class ForwardIterator>
bool CompGeom::is_convex_polygon(ForwardIterator begin, ForwardIterator end, double tol)
{
  return CompGeomSpecOne<ForwardIterator, typename std::iterator_traits<ForwardIterator>::value_type>::is_convex_polygon(begin, end, tol);
}

/// Determines whether a polygon (in 3D) is convex
template <class ForwardIterator>
bool CompGeom::is_convex_polygon(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double tol)
{
  return CompGeomSpecOne<ForwardIterator, typename std::iterator_traits<ForwardIterator>::value_type>::is_convex_polygon(begin, end, normal, tol);
}

/// Triangulates a convex polygon in O(n)
/**
 * \param source_begin the starting iterator to a connected, ccw-oriented, 
 *         convex polygon; additionally, the first point is assumed to be 
 *         connected to the last
 * \param source_end the ending iterator to the polygon container
 * \param target_begin the starting iterator to the container of triangles
 * \return the ending iterator to the container of triangles
 */
template <class ForwardIterator, class OutputIterator>
OutputIterator CompGeom::triangulate_convex_polygon(ForwardIterator source_begin, ForwardIterator source_end, OutputIterator target_begin)
{
  FILE_LOG(LOG_COMPGEOM) << "computing triangulation of polygon:" << std::endl;
  for (ForwardIterator i = source_begin; i != source_end; i++)
    FILE_LOG(LOG_COMPGEOM) << "  " << *i << std::endl;

  // special case: polygon is empty (return nothing)
  if (std::distance(source_begin, source_end) == 0)
    return target_begin;

  // compute the center of the points, and create new points
  std::list<Point3d> new_points;
  unsigned sz = 0;
  Point3d center = Point3d::zero();
  for (ForwardIterator i = source_begin; i != source_end; i++, sz++)
  {
    center += *i;
    new_points.push_back(*i);
  }

  // special case: empty polygon (return nothing)
  if (sz == 0)
    return target_begin;

  // setup the output iterator pointer
  OutputIterator current = target_begin;

  // scale center by the mean
  center /= sz;

  // now, create triangles
  for (std::list<Point3d>::const_iterator i = new_points.begin(); i != new_points.end(); i++)
  {
    // get the next point
    std::list<Point3d>::const_iterator j = i;
    j++;
    if (j == new_points.end())
      j = new_points.begin();    
  
    // create a triangle
    *current++ = Triangle(*i, *j, center);
  }

  return current;
}

/// Attempts to fit a plane to a set of points 
/**
 * The singular value decomposition is used to determine the plane that fits
 * the points best in a least-squares sense.
 * \param points the set of points (in 3D)
 * \param normal contains the "best" normal, on return
 * \param offset the offset such that, for any point on the plane x, 
 *        <normal, x> = offset
 * \return the maximum deviation from the plane
 */
template <class ForwardIterator>
double CompGeom::fit_plane(ForwardIterator begin, ForwardIterator end, Ravelin::Vector3d& normal, double& offset)
{
  return CompGeomSpecOne<ForwardIterator, typename std::iterator_traits<ForwardIterator>::value_type>::fit_plane(begin, end, normal, offset);
}

/// Projects a set of points onto a plane
/**
 * This method is intended to be used with fit_plane() to project the set of points exactly onto
 * a plane; note that this method is unnecessary if the points fit a plane exactly.
 */
template <class ForwardIterator>
void CompGeom::project_plane(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal, double offset)
{
  // form the projection matrix P = I - normal*normal'
  Ravelin::Matrix3d P;
  Ravelin::Opsd::outer_prod(normal, -normal, P);
  P += Ravelin::Matrix3d::identity(); 
 
  // project each point onto the plane
  for (ForwardIterator i = begin; i != end; i++)
  {
    // compute the projection
    Ravelin::Vector3d x = P * (*i);
    
    // P projects onto a plane parallel to the one we want; project directly onto the one
    // we want
    double remainder = offset - Ravelin::Vector3d::dot(x, normal);
    
    // add the remainder times the normal to x, and store it
    *i = x + (normal * remainder);
  }
}

/// Determines whether a 2D point is inside a polygon
/**
 * Adapted from O'Rourke, p. 244.
 * \param polygon a vector of points that describe a polygon (orientation irrelevant); each successive 
 *         vertex in the vector is connected to the previous vector to make the 
 *         polygon (polygon.front() and polygon.back() are connected as well)
 * \param point the point to test
 * \return inside_poly if point is inside the polygon, on_vertex if the point 
 *         coincides with a vertex, on_edge if the point lies on an edge of
 *         the polygon (but not a vertex), or outside_poly if the point is 
 *         outside of the polygon
 */
template <class ForwardIterator>
CompGeom::PolygonLocationType CompGeom::polygon_location(ForwardIterator begin, ForwardIterator end, const Point2d& point)
{
  const unsigned X = 0, Y = 1;
  int l_cross = 0, r_cross = 0;  
  
  // copy the polygon to a vector, shifted so that the point is at the origin
  std::vector<Point2d> poly_copy;
  for (ForwardIterator i=begin; i != end; i++)
    poly_copy.push_back(*i - point);
  
  // for each edge e = (i-1,i); see if crosses ray
  for (unsigned i=0; i< poly_copy.size(); i++)
  {
    // check whether the point is equal to a vertex
    if (std::fabs(poly_copy[i][X]) < NEAR_ZERO && std::fabs(poly_copy[i][Y]) < NEAR_ZERO)
      return ePolygonOnVertex;
          
    // determine i1
    unsigned i1 = (i + poly_copy.size() - 1) % poly_copy.size();

     // check whether e "straddles" the x axis, with bias above, below
    bool r_strad = (poly_copy[i][Y] > 0) != (poly_copy[i1][Y] > 0);
    bool l_strad = (poly_copy[i][Y] < 0) != (poly_copy[i1][Y] < 0);
    
    if (r_strad || l_strad)
    {
      // compute intersection of e with x axis  
      long double x = (poly_copy[i][X] * poly_copy[i1][Y] - poly_copy[i1][X] * poly_copy[i][Y]) / (poly_copy[i1][Y] - poly_copy[i][Y]);

      // crosses ray if strictly positive intersection
       if (r_strad && x > 0)
            r_cross++;
      if (l_strad && x < 0)
        l_cross++;
    }
  }

  // point on an edge if L/R cross counts are not the same parity
  if ((r_cross % 2) != (l_cross % 2))
    return ePolygonOnEdge;

  // otherwise, point is inside iff an odd number of crossings
  if (r_cross % 2 == 1)
    return ePolygonInside;
  else
    return ePolygonOutside;
}

/// Computes the area of a polygon in 2D
/**
 * \param poly a counter-clockwise oriented polygon in 3D
 */
template <class ForwardIterator>
double CompGeom::calc_polygon_area(ForwardIterator begin, ForwardIterator end)
{
  const unsigned X = 0, Y = 1, Z = 2;

  FILE_LOG(LOG_COMPGEOM) << "CompGeom::calc_polygon_area() entered" << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "  points: " << std::endl;
  for (ForwardIterator i = begin; i != end; i++)
    FILE_LOG(LOG_COMPGEOM) << "    " << *i << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "CompGeom::calc_polygon_area() exited" << std::endl;
  
  // compute the area of the polygon
  double area = 0;
  for (ForwardIterator i = begin; i != end; i++)
  {
    // get the next element; wrap around
    ForwardIterator j = i;
    j++;
    if (j == end)
      j = begin;

    // update the area
    area += (*i)[X] * (*j)[Y]  -  (*j)[X] * (*i)[Y];
  }
  area *= 0.5;

  return area;
}

/// Computes the area of a polygon in 3D
/**
 * \param poly a counter-clockwise oriented polygon in 3D
 * \param normal a vector normal to the plane that contains the points
 */
template <class ForwardIterator>
double CompGeom::calc_polygon_area(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the 3D to 2D projection matrix
  Ravelin::Matrix3d R = calc_3D_to_2D_matrix(normal);

  // project the points to 2D
  std::list<Point2d> points_2D;
  to_2D(begin, end, std::back_inserter(points_2D), R);

  // make sure that 2D polygon is ccw
  assert(ccw(points_2D.begin(), points_2D.end()));

  FILE_LOG(LOG_COMPGEOM) << "CompGeom::calc_polygon_area() entered" << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "  points (2D): " << std::endl;
  for (std::list<Point2d>::const_iterator i = points_2D.begin(); i != points_2D.end(); i++)
    FILE_LOG(LOG_COMPGEOM) << "    " << *i << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "CompGeom::calc_polygon_area() exited" << std::endl;
  
  return calc_polygon_area(points_2D.begin(), points_2D.end());
}

/// Computes the centroid of points on a plane
/**
 * \param begin an iterator to a polygon in 2D
 * \param end an iterator to a polygon in 2D
 */
template <class ForwardIterator>
Point2d CompGeom::calc_centroid_2D(ForwardIterator begin, ForwardIterator end)
{
  const unsigned X = 0, Y = 1;

  // now, compute the area of the polygon
  double area = 0;
  std::list<double> a;
  for (std::list<Point2d>::const_iterator i = begin; i != end; i++)
  {
    // get the next element; wrap around
    std::list<Point2d>::const_iterator j = i;
    j++;
    if (j == end)
      j = begin;

    // update the area
    a.push_back((*i)[X] * (*j)[Y]  -  (*j)[X] * (*i)[Y]);
    area += a.back();
  }
  area *= 0.5;

  // if the area is negative, negate it
  if (area < 0.0)
    area = -area;  
  
  // compute the 2D centroid
  Point2d centroid(0.0,0.0);
  std::list<double>::const_iterator i;
  std::list<Point2d>::const_iterator j;
  for (i = a.begin(), j = begin; i != a.end(); i++, j++)
  {
    // get the next point
    std::list<Point2d>::const_iterator k = j;
    k++;
    if (k == end)
      k = begin;

    centroid += (*j + *k) * (*i);      
  }
  centroid /= (area*6.0);
  
  return centroid;
}

/// Computes the 3D (2D) centroid of points on a plane
/**
 * \param poly a counter-clockwise oriented polygon in 3D
 * \param normal a vector normal to the plane that contains the points
 */
template <class ForwardIterator>
Point3d CompGeom::calc_centroid_2D(ForwardIterator begin, ForwardIterator end, const Ravelin::Vector3d& normal)
{
  const unsigned X = 0, Y = 1;

  // determine the pose of the points
  boost::shared_ptr<Ravelin::Pose3d> P = begin->pose;

  // get the 3D to 2D projection matrix
  Ravelin::Matrix3d R = calc_3D_to_2D_matrix(normal);

  // get the offset
  double offset = determine_3D_to_2D_offset(*begin, R);

  // project the points to 2D
  std::list<Point2d> points_2D;
  std::insert_iterator<std::list<Point2d> > ii(points_2D, points_2D.begin());
  to_2D(begin, end, ii, R);

  // make sure that 2D polygon is ccw
  assert(ccw(points_2D.begin(), points_2D.end()));

  FILE_LOG(LOG_COMPGEOM) << "polygon: " << std::endl;
  for (ForwardIterator i = begin; i != end; i++)
    FILE_LOG(LOG_COMPGEOM) << *i << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "2D points:" << std::endl;
  for (std::list<Point2d>::const_iterator i = points_2D.begin(); i != points_2D.end(); i++)
    FILE_LOG(LOG_COMPGEOM) << "    " << *i << std::endl;

  // now, compute the area of the polygon
  double area = 0;
  std::list<double> a;
  for (std::list<Point2d>::const_iterator i = points_2D.begin(); i != points_2D.end(); i++)
  {
    // get the next element; wrap around
    std::list<Point2d>::const_iterator j = i;
    j++;
    if (j == points_2D.end())
      j = points_2D.begin();

    // update the area
    a.push_back((*i)[X] * (*j)[Y]  -  (*j)[X] * (*i)[Y]);
    area += a.back();
  }
  area *= 0.5;

  FILE_LOG(LOG_COMPGEOM) << "normal: " << normal << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "area: " << area << std::endl;
  assert(area >= 0);
  
  // compute the 2D centroid
  Point2d centroid(0,0);
  std::list<double>::const_iterator i;
  std::list<Point2d>::const_iterator j;
  for (i = a.begin(), j = points_2D.begin(); i != a.end(); i++, j++)
  {
    // get the next point
    std::list<Point2d>::const_iterator k = j;
    k++;
    if (k == points_2D.end())
      k = points_2D.begin();

    centroid += (*j + *k) * (*i);      
  }
  centroid /= (area*6.0);
  
  // get the transpose (i.e., inverse) of the rotation matrix
  Ravelin::Matrix3d RT = Ravelin::Matrix3d::transpose(R);
  
  FILE_LOG(LOG_COMPGEOM) << "2D centroid: " << centroid << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "RT: " << std::endl << RT;
  
  // project the centroid back to 3D
  return Point3d(to_3D(Ravelin::Origin2d(centroid), RT, offset), P);
}

/// Computes the minimum area bounding rectangle of a set of points
/**
 * Uses 2D convex hull and rotating calipers method; runs in O(N lg N) time
 * On return, x1, x2, x3, and x4 are the four vertices of the bounding
 * rectangle (ordered as edges).
 */
template <class ForwardIterator>
void CompGeom::calc_min_area_bounding_rect(ForwardIterator begin, ForwardIterator end, Point2d& x1, Point2d& x2, Point2d& x3, Point2d& x4)
{
  CompGeomSpecOne<ForwardIterator, typename std::iterator_traits<ForwardIterator>::value_type>::calc_min_area_bounding_rect(begin, end, x1, x2, x3, x4);
}

/**
 * Intersects two 2D triangles.
 * \note this method adapted from www.geometrictools.com
 */
template <class OutputIterator>
OutputIterator CompGeom::intersect_tris(const Point2d t1[3], const Point2d t2[3], OutputIterator output_begin)
{
  const unsigned X = 0, Y = 1;

  // verify that both triangles are ccw
  assert(ccw(t1, t1+3));
  assert(ccw(t2, t2+3));

  // init the potential intersection to t2 
  Point2d isects[6];
  isects[0] = t2[0];
  isects[1] = t2[1];
  isects[2] = t2[2];
  unsigned nisects = 3;

  // clip against edges
  for (unsigned i1=2, i0=0; i0 < 3; i1 = i0, i0++)
  {
    Ravelin::Vector2d kN(t1[i1][Y] - t1[i0][Y], t1[i0][X] - t1[i1][X]);
    double fC = kN.dot(t1[i1]);
    clip_convex_polygon_against_line(kN, fC, nisects, isects);

    // look for no intersection
    if (nisects == 0)
      return output_begin;
  }

  // copy to output
  for (unsigned i=0; i< nisects; i++)
    *output_begin++ = isects[i];

  return output_begin;
}

/// Intersects a line segment and a triangle in 2D
/*
 * \param seg the line segment
 * \param tri the triangle 
 * \param output_begin an iterator to the beginning of a container of Point2d; 
 *        points of intersection will be stored here on return
 * \return an iterator to the end of a container of Point2d; 
 *        points of intersection will be stored here on return
 * \note this code adapted from http://www.geometrictools.com
 */
template <class OutputIterator>
OutputIterator CompGeom::intersect_seg_tri(const LineSeg2& seg, const Point2d tri[3], OutputIterator output_begin)
{
  Point2d isect, isect2;
  CompGeom::SegTriIntersectType code = intersect_seg_tri(seg, tri, isect, isect2);
  
  switch (code)
  {
    case eSegTriNoIntersect:
      break;

    case eSegTriVertex:
    case eSegTriEdge:
    case eSegTriPlanarIntersect: 
      *output_begin++ = isect; 

    case eSegTriEdgeOverlap:
    case eSegTriInside:
      *output_begin++ = isect;
      *output_begin++ = isect2;
      break;

    default:
      break;
  }

  return output_begin;
}
    

/*
  double afDist[3];
  int aiSign[3], iPositive, iNegative, iZero;

  // this block of code adapted from TriangleLineRelations()
  Vector2 origin = (seg.first + seg.second)*0.5;
  Vector2 seg_dir = Vector2::normalize(seg.second - seg.first);
  iPositive = 0;
  iNegative = 0;
  iZero = 0;
  for (unsigned i=0; i< 3; i++)
  {
    Vector2 kDiff = tri[i] - origin;
    afDist[i] = kDiff.dot_perp(seg_dir);
    if (afDist[i] > NEAR_ZERO)
    {
      aiSign[i] = 1;
      iPositive++;
    }
    else if (afDist[i] < -NEAR_ZERO)
    {
      aiSign[i] = -1;
      iNegative++;
    }
    else
    {
      afDist[i] = 0.0;
      aiSign[i] = 0;
      iZero++;
    }
  }

  if (iPositive == 3 || iNegative == 3)
  {
    // empty set
    return output_begin;
  }
  else
  {
    double afParam[2];

}
*/

/// Gets the parameter of a point on a line, v = p0 + dir*t, -inf <= t <= inf
template <class T, class U>
double CompGeom::determine_line_param(const T& p0, const U& dir, const T& v)
{
  double dir_norm = dir.norm();
  if (dir_norm < NEAR_ZERO)
    throw NumericalException("Attempting to normalize zero length vector");
  return sgn((v - p0).norm()/dir_norm, (v - p0).dot(dir));
}

/// Intersects two line segments in 2D
/**
 * \param seg1 the first line segment
 * \param seg2 the second line segment
 * \param output_begin an iterator to the beginning of a container of Point2d; 
 *        points of intersection will be stored here on return
 * \return an iterator to the end of a container of Point2d; 
 *        points of intersection will be stored here on return
 */ 
template <class OutputIterator>
OutputIterator CompGeom::intersect_segs(const LineSeg2& s1, const LineSeg2& s2, OutputIterator output_begin)
{
  // do the intersection
  Point2d isect, isect2;
  SegSegIntersectType isect_type = intersect_segs(s1, s2, isect, isect2);

  // switch on the intersection type
  switch (isect_type)
  {
    case eSegSegNoIntersect:
      break;

    case eSegSegIntersect:
    case eSegSegVertex:
      *output_begin++ = isect;
      break; 

    case eSegSegEdge:
      *output_begin++ = isect;
      *output_begin++ = isect2;
      break; 
  }

  return output_begin;
}

