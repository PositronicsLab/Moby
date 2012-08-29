/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iostream>
#include <cmath>
#include <set>
#include <Moby/sorted_pair>
#include <Moby/Constants.h>
#include <Moby/Log.h>
#include <Moby/CompGeom.h>
#include <Moby/Plane.h>
#include <Moby/Polyhedron.h>

using namespace Moby;
using std::vector;
using std::make_pair;
using std::pair;

// ******************************************************************
// functions for debugging
// ******************************************************************

static void write_poly(const Polyhedron& p, std::ostream& out)
{
  for (unsigned i=0; i< p.get_mesh().num_tris(); i++)
  {
    Vector3 color((Real) rand()/RAND_MAX, (Real) rand()/RAND_MAX, (Real) rand()/RAND_MAX);
    Triangle t = p.get_mesh().get_triangle(i);
    t.to_vrml(out, color);
  }
  
  const vector<Vector3>& v = p.get_mesh().get_vertices();
  for (unsigned i=0; i< v.size(); i++)
  {
    out << "Transform {" << std::endl;
    out << "  translation " << v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
    out << "  scale 0.1 0.1 0.1" << std::endl;
    out << "  children Shape { " << std::endl;
    out << "    geometry Text { " << std::endl;
    out << "      string \"" << i << "\"" << std::endl;
    out << "      } } }" << std::endl;
  }
}

static bool consistent(const vector<IndexedTri>& f, unsigned a, unsigned b)
{
  unsigned count = 0;
  for (unsigned i=0; i< f.size(); i++)
    if ((f[i].a == a || f[i].b == a || f[i].c == a) &&
        (f[i].a == b || f[i].b == b || f[i].c == b))
      count++;

  return (count == 2);
}

static void write_poly2(const Polyhedron& p, std::ostream& out)
{
  const unsigned X = 0, Y = 1, Z = 2;
 
  const vector<Vector3>& v = p.get_mesh().get_vertices();
  for (unsigned i=0; i< v.size(); i++)
  {
    out << "Transform {" << std::endl;
    out << "  translation " << v[i][0] << " " << v[i][1] << " " << v[i][2] << std::endl;
    out << "  scale 0.1 0.1 0.1" << std::endl;
    out << "  children Shape { " << std::endl;
    out << "    geometry Text { " << std::endl;
    out << "      string \"" << i << "\"" << std::endl;
    out << "      } } }" << std::endl;
  }
 
  // get the vertices and the facets
  const std::vector<Vector3>& vertices = p.get_vertices();
  const std::vector<IndexedTri>& facets = p.get_facets();

  // first write the consistent edges
  out << "Shape {" << std::endl;
  out << "  appearance Appearance { material Material { diffuseColor 1 1 1 } }" << std::endl;
  out << "  geometry ";
  out << "IndexedLineSet {" << std::endl;
  out << "    coord Coordinate { point [ ";
  for (unsigned i=0; i< vertices.size(); i++)
    out << vertices[i][X] << " " << vertices[i][Y] << " " << vertices[i][Z] << ", ";
  out << " ] }" << std::endl;
  out << "    coordIndex [ ";
  for (unsigned i=0; i< facets.size(); i++)
  {
    if (consistent(facets, facets[i].a, facets[i].b))
      out << facets[i].a << " " << facets[i].b << " -1, ";
    if (consistent(facets, facets[i].b, facets[i].c))
      out << facets[i].b << " " << facets[i].c << " -1, ";
    if (consistent(facets, facets[i].a, facets[i].c))
      out << facets[i].a << " " << facets[i].c << " -1, ";
  }
  out << " ] } }" << std::endl;

  // now write the inconsistent edges
  std::set<pair<unsigned, unsigned> > marked;
  out << "Shape {" << std::endl;
  out << "  appearance Appearance { material Material { diffuseColor 1 0 0 } }" << std::endl;
  out << "  geometry ";
  out << "IndexedLineSet {" << std::endl;
  out << "    coord Coordinate { point [ ";
  for (unsigned i=0; i< vertices.size(); i++)
    out << vertices[i][X] << " " << vertices[i][Y] << " " << vertices[i][Z] << ", ";
  out << " ] }" << std::endl;
  out << "    coordIndex [ ";
  for (unsigned i=0; i< facets.size(); i++)
  {
    if (marked.find(make_pair(facets[i].a, facets[i].b)) == marked.end() && 
        !consistent(facets, facets[i].a, facets[i].b))
    {
      out << facets[i].a << " " << facets[i].b << " -1, ";
      marked.insert(make_pair(facets[i].a, facets[i].b));
std::cout << "edge: " << facets[i].a << " " << facets[i].b << std::endl;
    }
    if (marked.find(make_pair(facets[i].b, facets[i].c)) == marked.end() && 
        !consistent(facets, facets[i].b, facets[i].c))
    {
      out << facets[i].b << " " << facets[i].c << " -1, ";
      marked.insert(make_pair(facets[i].b, facets[i].c));
std::cout << "edge: " << facets[i].b << " " << facets[i].c << std::endl;
    }
    if (marked.find(make_pair(facets[i].a, facets[i].c)) == marked.end() && 
        !consistent(facets, facets[i].a, facets[i].c))
    {
      out << facets[i].a << " " << facets[i].c << " -1, ";
      marked.insert(make_pair(facets[i].a, facets[i].c));
std::cout << "edge: " << facets[i].a << " " << facets[i].c << std::endl;
    }
  }
  out << " ] } }" << std::endl;
}

static void write_poly(const Polyhedron& p, const char* fname)
{
  std::ofstream out(fname);
  out << "#VRML V2.0 utf8" << std::endl;
  write_poly2(p, out);
  out.close();
}

/// Gets the requested vertex
unsigned get_vertex(const vector<Vector3>& vertices, const Vector3& v)
{
  unsigned n = vertices.size();

  for (unsigned i=0; i< n; i++)
    if ((vertices[i] - v).norm() < NEAR_ZERO)
      return i;

  assert(false);
  return 0;
}

/// Verifies that the edge exists
bool find_edge(const vector<IndexedTri>& f, unsigned a, unsigned b)
{
  if (a == b)
    return true;

  for (unsigned i=0; i< f.size(); i++)
    if ((f[i].a == a || f[i].b == a || f[i].c == a) &&
        (f[i].b == b || f[i].b == b || f[i].c == b))
      return true;

  return false;
}

// ******************************************************************
// functions for debugging (END)
// ******************************************************************

/// Constructs a polyhedron from an indexed triangle array
Polyhedron::Polyhedron(const IndexedTriArray& mesh)
{
  // indicate that convexity has not been computed
  _convexity_computed = false;

  // copy the mesh
  _mesh = mesh;

  // calculate the bounding box
  calc_bounding_box();
}

/// Copies a polyhedron
void Polyhedron::operator=(const Polyhedron& p)
{
  _mesh = p._mesh;
  _bb_min = p._bb_min;
  _bb_max = p._bb_max;
  _convexity = p._convexity;
  _convexity_computed = p._convexity_computed;
}

/// Computes the Minkowski sum of two convex polyhedra
PolyhedronPtr Polyhedron::minkowski(Polyhedron& p1, const Matrix4& T1, Polyhedron& p2, const Matrix4& T2, bool reflect_p2)
{
  // verify that both polyhedra are convex
  if (!p1.is_convex() || !p2.is_convex())
    throw std::runtime_error("Polyhedron::minkowski() only operates on convex polyhedra");

  // we'll transform p2 to p1's frame
  Matrix4 T2_to_T1 = Matrix4::inverse_transform(T1) * T2;
  Polyhedron p2_copy = p2;
  p2_copy.transform(T2_to_T1);

  // compute the minkowski sum
  std::list<Vector3> points;
  for (unsigned i=0; i< p1.get_vertices().size(); i++)
    for (unsigned j=0; j< p2_copy.get_vertices().size(); j++)
    {
      // add vertex from p1 to vertex from p2
      Vector3 v = p1.get_vertices()[i];
      if (reflect_p2)
        v -= p2_copy.get_vertices()[j];
      else
        v += p2_copy.get_vertices()[j];
      points.push_back(v);
    }

  // compute the convex hull of the points
  return CompGeom::calc_convex_hull_3D(points.begin(), points.end()); 
}

/// Checks whether this polyhedron is degenerate
bool Polyhedron::degenerate() const
{
  const std::vector<Vector3>& vertices = get_vertices();
  const std::vector<IndexedTri>& facets = get_facets();

  for (unsigned i=0; i< facets.size(); i++)
  {
    if ((vertices[facets[i].a] - vertices[facets[i].b]).norm() < NEAR_ZERO)
      return true;
    if ((vertices[facets[i].b] - vertices[facets[i].c]).norm() < NEAR_ZERO)
      return true;
    if ((vertices[facets[i].c] - vertices[facets[i].a]).norm() < NEAR_ZERO)
      return true;
  }

  return false;
}

/// Checks whether this polyhedron is consistent
bool Polyhedron::consistent() const
{
  // determine the set of edges for this polyhedron
  std::map<sorted_pair<unsigned>, std::list<unsigned> > edges;
  const std::vector<IndexedTri>& facets = get_facets();
  for (unsigned i=0; i< facets.size(); i++)
  {
    sorted_pair<unsigned> s1(facets[i].a, facets[i].b);
    sorted_pair<unsigned> s2(facets[i].b, facets[i].c);
    sorted_pair<unsigned> s3(facets[i].c, facets[i].a);
    edges[s1].push_back(i);
    edges[s2].push_back(i);
    edges[s3].push_back(i);
  }
  
  // verify that each edge is coincident to exactly two facets
  for (std::map<sorted_pair<unsigned>, std::list<unsigned> >::const_iterator i = edges.begin(); i != edges.end(); i++)
    if (i->second.size() != 2)
      return false;

  return true;
}

/// Intersects a plane with a line segment
/**
 * \note utility function for split should return a single point
 */
Vector3 Polyhedron::intersect_plane(const Vector3& normal, Real d, const Vector3& p1, const Vector3& p2)
{
  // get a point on the plane
  Vector3 pplane = normal * d;

  FILE_LOG(LOG_COMPGEOM) << "equation of plane is: " << normal << "; " << d << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "line segment is: " << p1 << "; " << p2 << std::endl;
  
  // compute u, the point of intersection
  Real denom = Vector3::dot(normal, p2 - p1);
  assert(denom != 0.0);
  Real num = Vector3::dot(normal, pplane - p1);
  Real u = num / denom;

  FILE_LOG(LOG_COMPGEOM) << "point of intersection is " << (p1 + u*(p2-p1)) << std::endl;
  FILE_LOG(LOG_COMPGEOM) << "u = " << u << std::endl;
  assert(u >= 0 && u <= 1);
        
  return p1 + u*(p2-p1);
} 

/// Determines whether the specified point is strictly inside this polyhedron
bool Polyhedron::inside(const Vector3& point, Real tol)
{
  const unsigned THREE_D = 3;
  
  // if the polyhedron is not convex, using the location method
  if (!is_convex())
    return location(point, tol) == eInside;
  
  // check whether the point is inside the bounding box
  for (unsigned i=0; i< THREE_D; i++)
    if (point[i] < _bb_min[i] - tol || point[i] > _bb_max[i] + tol)
      return false;

  // check whether the point is outside or coplanar with any facets    
  const std::vector<Vector3>& vertices = get_vertices();
  const std::vector<IndexedTri>& facets = get_facets();
  for (unsigned i=0; i< facets.size(); i++)
    if (CompGeom::volume_sign(vertices[facets[i].a], vertices[facets[i].b], vertices[facets[i].c], point) <= tol)
      return false;

  return true;
}

/// Determines whether the specified point is in or on this polyhedron
bool Polyhedron::inside_or_on(const Vector3& point, Real tol) 
{
  const unsigned THREE_D = 3;
 
  // if the polyhedron is not convex, using the location method
  if (!is_convex())
    return location(point, tol) != eOutside;
  
  // check whether the point is inside the bounding box
  for (unsigned i=0; i< THREE_D; i++)
    if (point[i] < _bb_min[i] - tol || point[i] > _bb_max[i] + tol)
      return false;
    
  // check whether the point is outside any facets    
  const std::vector<Vector3>& vertices = get_vertices();
  const std::vector<IndexedTri>& facets = get_facets();
  for (unsigned i=0; i< facets.size(); i++)
    if (CompGeom::volume_sign(vertices[facets[i].a], vertices[facets[i].b], vertices[facets[i].c], point) < -tol)
      return false;

  return true;
}

/// Determines the location of the specified point with respect to this polyhedron
/**
 * This method handles queries when the polyhedron is not convex or more detail is necessary.
 * Adapted from O'Rourke, p. 247-250.  Runs in worst-case time O(f), where
 * f is the number of facets of the polyhedron.
 */
Polyhedron::LocationType Polyhedron::location(const Vector3& point, Real tol) const
{
  const unsigned THREE_D = 3;
  unsigned isects;

  // check whether the point is inside the bounding box
  for (unsigned i=0; i< THREE_D; i++)
    if (point[i] < _bb_min[i] - tol || point[i] > _bb_max[i] + tol)
      return eOutside;
  
  // get the length of the diagonal of the polyhedron bounding box
  long double D = 0;
  for (unsigned i=0; i< THREE_D; i++)
  {
    double v = _bb_max[i] - _bb_min[i];
    D += v*v;
  }
  D = std::sqrt(D);  

  // the bounding radius is ceil(D) + 1
  double bounding_radius = (double) (std::ceil(D) + 1);

  // loop while there is degeneracy
  const std::vector<Vector3>& vertices = get_vertices();
  const std::vector<IndexedTri>& facets = get_facets();
  while (true)
  {
    // generate a random ray of length bounding_radius
    Vector3 r_knot;
    for (unsigned i=0; i< THREE_D; i++)
      r_knot[i] = (double) rand() - (RAND_MAX / 2);
    r_knot.normalize();
    r_knot *= bounding_radius;

    // add r to the point
    Vector3 r = point + r_knot;

    // setup the segment
    LineSeg3 seg(point, r);

    // reset the number of face intersections to zero
    isects = 0;

    // for each face, determine whether there is an intersection
    bool degenerate_ray = false;
    for (unsigned i=0; i< facets.size(); i++)
    {
      // determine the intersection of the face and the ray
      Vector3 dummy;
      Triangle tri(vertices[facets[i].a], vertices[facets[i].b], vertices[facets[i].c]);
      CompGeom::SegTriIntersectType code = CompGeom::intersect_seg_tri(seg, tri, dummy, dummy, tol);

      // if ray is degenerate, break out to outer while to generate another
      if (code == CompGeom::eSegTriInside || code == CompGeom::eSegTriEdgeOverlap ||
        code == CompGeom::eSegTriPlanarIntersect || code == CompGeom::eSegTriInclVertex ||
        code == CompGeom::eSegTriInclEdge)
       {
        degenerate_ray = true;
        break;
      }

      // if ray hits face at interior point, increment crossings
      if (code == CompGeom::eSegTriInclFace)
      {
        isects++;
        continue;
      }

      // if query endpoints responds on a vertex/edge/face, return appropriately
      if (code == CompGeom::eSegTriVertex)
        return eOnVertex;
      if (code == CompGeom::eSegTriEdge)
        return eOnEdge;
      if (code == CompGeom::eSegTriFace)
        return eOnFace;

      // if ray misses triangle, do nothing...
      if (code == CompGeom::eSegTriNoIntersect)
        continue;

      // should not still be here
      assert(false);
    }

    // if the ray was not degenerate, break
    if (!degenerate_ray)
      break;
  }

  // point is strictly interior to polyhedron iff an odd number of crossings
  if (isects % 2 == 1)
    return eInside;
  else
    return eOutside;
}

/// Translates this polyhedron by the given vector
/**
 * \note none of the vertex or triangle pointers change; rather the data
 *       that they point to changes
 */
void Polyhedron::translate(const Vector3& x)
{
  // translate underlying mesh
  _mesh = _mesh.translate(x); 

  // translate the bounding box
  _bb_min += x;
  _bb_max += x;
}

/// Transforms this polyhedron by the given transformation matrix
/**
 * \note none of the vertex or triangle pointers change; rather the data
 *       that they point to changes
 */
void Polyhedron::transform(const Matrix4& T)
{
  // transform underlying mesh
  _mesh = _mesh.transform(T); 

  // compute the new bounding box
  calc_bounding_box();
}

/// Calculates the bounding box
void Polyhedron::calc_bounding_box()
{
  const unsigned THREE_D = 3;
  
  const std::vector<Vector3>& vertices = get_vertices();
  if (vertices.empty())
  {
    _bb_min = ZEROS_3;
    _bb_max = ZEROS_3;
  }
  else
  {
    // setup the bounding box
    _bb_min = vertices.front();
    _bb_max = vertices.front();
    for (unsigned i=1; i< vertices.size(); i++)
      for (unsigned j=0; j< THREE_D; j++)
      {
        if (vertices[i][j] < _bb_min[j])
          _bb_min[j] = vertices[i][j];
        else if (vertices[i][j] > _bb_max[j])
          _bb_max[j] = vertices[i][j];
      }
  }
}

/// Determine whether the given polyhedron is convex
void Polyhedron::determine_convexity()
{
  // set convexity to -inf to begin
  _convexity = -std::numeric_limits<Real>::max();

  // determine the set of edges for this polyhedron
  const std::vector<Vector3>& vertices = get_vertices();
  const std::vector<IndexedTri>& facets = get_facets();
  std::map<sorted_pair<unsigned>, std::list<unsigned> > edges;
  for (unsigned i=0; i< facets.size(); i++)
  {
    sorted_pair<unsigned> s1(facets[i].a, facets[i].b);
    sorted_pair<unsigned> s2(facets[i].b, facets[i].c);
    sorted_pair<unsigned> s3(facets[i].c, facets[i].a);
    edges[s1].push_back(i);
    edges[s2].push_back(i);
    edges[s3].push_back(i);
  }

  // determine whether each edge is convex
  for (std::map<sorted_pair<unsigned>, std::list<unsigned> >::const_iterator i = edges.begin(); i != edges.end(); i++)
  {
    // verify that there are two facets
    assert(i->second.size() == 2);

    // get the two facets
    const IndexedTri& f1 = facets[i->second.front()];
    const IndexedTri& f2 = facets[i->second.back()];

    // get the triangle corresponding to f1
    Triangle tri_f1(vertices[f1.a], vertices[f1.b], vertices[f1.c]);

    // get the plane going through face 1
    const Vector3& normal = tri_f1.calc_normal();
    Real d = tri_f1.calc_offset(normal);

    // get the vertex of f2 not on this edge
    unsigned vert;
    if ((f2.a == i->first.first && f2.b == i->first.second) ||
        (f2.a == i->first.second && f2.b == i->first.first))
      vert = f2.c;
    else if ((f2.a == i->first.first && f2.c == i->first.second) ||
                (f2.a == i->first.second && f2.c == i->first.first))
      vert = f2.b;  
    else
    {
      assert((f2.b == i->first.first && f2.c == i->first.second) ||
              (f2.b == i->first.second && f2.c == i->first.first));
      vert = f2.c;
    }  

    // if the vertex is on the positive side of the plane, the edge is non-convex
    _convexity = std::max(_convexity, normal.dot(vertices[vert]) - d);
  }

  // indicate that convexity has been calculated
  _convexity_computed = true;
}

/// Sends this polyhedron to the specified stream using VRML
void Polyhedron::to_vrml(std::ostream& out, const Polyhedron& p, Vector3 diffuse_color, bool wireframe)
{
  const unsigned X = 0, Y = 1, Z = 2;
  
  // get the vertices and the facets
  const std::vector<Vector3>& vertices = p.get_vertices();
  const std::vector<IndexedTri>& facets = p.get_facets();

  out << "Shape {" << std::endl;
  out << "  appearance Appearance { material Material { diffuseColor " << diffuse_color[0] << " " << diffuse_color[1] << " " << diffuse_color[2] << " } }" << std::endl;
  out << "  geometry ";
  if (!wireframe)
    out << "IndexedFaceSet {" << std::endl;
  else
    out << "IndexedLineSet {" << std::endl;
  out << "    coord Coordinate { point [ ";
  for (unsigned i=0; i< vertices.size(); i++)
    out << vertices[i][X] << " " << vertices[i][Y] << " " << vertices[i][Z] << ", ";
  out << " ] }" << std::endl;
  out << "    coordIndex [ ";
  if (!wireframe)
    for (unsigned i=0; i< facets.size(); i++)
      out << facets[i].a << " " << facets[i].b << " " << facets[i].c << " -1, ";
  else
    for (unsigned i=0; i< facets.size(); i++)
      out << facets[i].a << " " << facets[i].b << " " << facets[i].c << " " << facets[i].a << " -1, ";  
  out << " ] } }" << std::endl;
}

/// Utility method for calc_volume()
void Polyhedron::calc_subexpressions(Real w0, Real w1, Real w2, Real& f1, Real& f2, Real& f3, Real& g0, Real& g1, Real& g2)
{
  Real temp0 = w0 + w1;
  f1 = temp0 + w2;
  Real temp1 = w0 * w0;
  Real temp2 = temp1 + w1 * temp0;
  f2 = temp2 + w2 * f1;
  f3 = w0 * temp1 + w1 * temp2 + w2 * f2;
  g0 = f2 + w0 * (f1 + w0);
  g1 = f2 + w1 * (f1 + w1);
  g2 = f2 + w2 * (f1 + w2);
}

/// Calculates the volume of this polyhedron
Real Polyhedron::calc_volume() const
{
  Real f1x, f2x, f3x, g0x, g1x, g2x, f1y, f2y, f3y, g0y, g1y, g2y;
  Real f1z, f2z, f3z, g0z, g1z, g2z;
  const int X = 0;
  const int Y = 1;
  const int Z = 2;

  // zero integral
  std::vector<Real> integral(10, 0);

  // order: 1, x, y, z, x^2, y^2, xy, yz, zx
  const std::vector<Vector3>& vertices = get_vertices();
  const std::vector<IndexedTri>& facets = get_facets();
  for (unsigned t=0; t< facets.size(); t++)
  {
    // get vertices of triangle t
    const Real x0 = vertices[facets[t].a][X];
    const Real y0 = vertices[facets[t].a][Y];
    const Real z0 = vertices[facets[t].a][Z];
    const Real x1 = vertices[facets[t].b][X];
    const Real y1 = vertices[facets[t].b][Y];
    const Real z1 = vertices[facets[t].b][Z];
    const Real x2 = vertices[facets[t].c][X];
    const Real y2 = vertices[facets[t].c][Y];
    const Real z2 = vertices[facets[t].c][Z];

    // get edges and cross product of edges
    Real a1 = x1 - x0;
    Real b1 = y1 - y0;
    Real c1 = z1 - z0;
    Real a2 = x2 - x0;
    Real b2 = y2 - y0;
    Real c2 = z2 - z0;
    Real d0 = b1 * c2 - b2 * c1;
    Real d1 = a2 * c1 - a1 * c2;
    Real d2 = a1 * b2 - a2 * b1;

    // compute integral terms
    calc_subexpressions(x0,x1,x2,f1x,f2x,f3x,g0x,g1x,g2x);
    calc_subexpressions(y0,y1,y2,f1y,f2y,f3y,g0y,g1y,g2y);
    calc_subexpressions(z0,z1,z2,f1z,f2z,f3z,g0z,g1z,g2z);

    // update integrals
    integral[0] += d0 * f1x;
    integral[1] += d0 * f2x;
    integral[2] += d1 * f2y;
    integral[3] += d2 * f2z;
    integral[4] += d0 * f3x;
    integral[5] += d1 * f3y;
    integral[6] += d2 * f3z;
    integral[7] += d0 * (y0 * g0x + y1 * g1x + y2 * g2x);
    integral[8] += d1 * (z0 * g0y + z1 * g1y + z2 * g2y);
    integral[9] += d2 * (x0 * g0z + x1 * g1z + x2 * g2z);
  }

  return integral[0] / 6;
}

/// Computes the signed distance from the polyhedron to a point
Real Polyhedron::calc_signed_distance(const Vector3& point, unsigned& closest_facet) 
{
  // initialize minimum distance
  Real min_dist = std::numeric_limits<Real>::max();

  // init the closest point (unused)
  Vector3 closest_point;

  // compute the minimum distance to all triangles
  const std::vector<Vector3>& vertices = get_vertices();
  const std::vector<IndexedTri>& facets = get_facets();
  for (unsigned i=0; i< facets.size(); i++)
  {
    Triangle tri(vertices[facets[i].a], vertices[facets[i].b], vertices[facets[i].c]);
    Real dist = std::sqrt(Triangle::calc_sq_dist(tri, point, closest_point));
    if (dist < min_dist)
    {
      min_dist = dist;
      closest_facet = i;
    }
  }

  // sign the distance if the point is inside
  if (inside(point))
    min_dist = -min_dist;

  // return the distance
  return min_dist;
}

/// Removes triangles from mesh that are inside p (helper function for Boolean operations on polyhedra)
void Polyhedron::remove_inside(IndexedTriArray& mesh, Polyhedron& p, bool remove_shared)
{
  // get the vertices and facets
  const vector<Vector3>& verts = mesh.get_vertices();
  vector<IndexedTri> facets = mesh.get_facets();

  // check whether each vertex is inside or outside the polyhedron
  vector<LocationType> loc(verts.size());
  for (unsigned i=0; i< loc.size(); i++)
    loc[i] = p.location(verts[i]);

  // now determine the facets that are inside the mesh (and shared)
  for (unsigned i=0; i< facets.size(); i++)
  {
    // if we are removing shared, see whether all points are on the polyhedron
    if (remove_shared)
    {
      if (loc[facets[i].a] == eOnVertex && loc[facets[i].b] == eOnVertex &&
          loc[facets[i].c] == eOnVertex)
      {
        facets[i] = facets.back();
        facets.pop_back();
        i--;
        continue;
      }
    }
    
    // see whether a single point is inside 
    if (loc[facets[i].a] == eInside || loc[facets[i].b] == eInside ||
        loc[facets[i].c] == eInside)
    {
      facets[i] = facets.back();
      facets.pop_back();
      i--;
      continue;
    }
  }

  // now reinstantiate the mesh
  mesh = IndexedTriArray(verts.begin(), verts.end(), facets.begin(), facets.end());
}

/// Removes triangles from mesh that are outside p (helper function for Boolean operations on polyhedra)
void Polyhedron::remove_outside(IndexedTriArray& mesh, Polyhedron& p, bool remove_shared)
{
  // get the vertices and facets
  const vector<Vector3>& verts = mesh.get_vertices();
  vector<IndexedTri> facets = mesh.get_facets();

  // check whether each vertex is inside or outside the polyhedron
  vector<LocationType> loc(verts.size());
  for (unsigned i=0; i< loc.size(); i++)
    loc[i] = p.location(verts[i]);

  // now determine the facets that are inside the mesh (and shared)
  for (unsigned i=0; i< facets.size(); i++)
  {
    // if we are removing shared, see whether all points are on the polyhedron
    if (remove_shared)
    {
      if (loc[facets[i].a] == eOnVertex && loc[facets[i].b] == eOnVertex && 
          loc[facets[i].c] == eOnVertex)
      {
        facets[i] = facets.back();
        facets.pop_back();
        i--;
        continue;
      }
    }

    // see whether a single point lies outside 
    if (loc[facets[i].a] == eOutside || loc[facets[i].b] == eOutside ||
        loc[facets[i].c] == eOutside)
    {
      facets[i] = facets.back();
      facets.pop_back();
      i--;
    }
  }

  // now reinstantiate the mesh
  mesh = IndexedTriArray(verts.begin(), verts.end(), facets.begin(), facets.end());
}

/// Determines whether triangle a bisects triangle b (helper function for slice)
bool Polyhedron::bisects(const Triangle& a, const Triangle& b)
{
  // get signed distance of vertices from a
  Plane p(a);

  Real da = p.calc_signed_distance(b.a);
  Real db = p.calc_signed_distance(b.b);
  Real dc = p.calc_signed_distance(b.c);

  if (da > 0.0 && db > 0.0 && dc > 0.0)
    return false;
  if (da < -0.0 && db < -0.0 && dc < -0.0)
    return false;

  return true;
}

/// Determines whether a vertex is present in the mesh already
bool Polyhedron::find_vertex(const vector<Vector3>& vertices, const Vector3& v)
{
  unsigned n = vertices.size();

  for (unsigned i=0; i< n; i++)
    if ((vertices[i] - v).norm() < 1e-8)
      return true;

  return false;
}

/// Adds a vertex to a mesh if the vertex doesn't already exist in the mesh
/**
 * Returns the index of the mesh vertex.
 */
unsigned Polyhedron::add_vertex(vector<Vector3>& vertices, const Vector3& v)
{
  unsigned n = vertices.size();

  for (unsigned i=0; i< n; i++)
    if ((vertices[i] - v).norm() < NEAR_ZERO)
      return i;

  vertices.push_back(v);
  return n;
}

/// Bisects edge ab with ac and cb
/**
 * \pre ab is part of the winding convention (ccw or cw)
 */
void Polyhedron::replace_edge(const vector<Vector3>& v, vector<IndexedTri>& f, unsigned a, unsigned b, unsigned c, vector<unsigned>& del_list)
{
  // clear the deletion list
  del_list.clear();

  for (unsigned i=0; i< f.size(); i++)
  {
    bool reverse = false;
    unsigned e = std::numeric_limits<unsigned>::max();

    if (a == f[i].a && b == f[i].b)
      e = f[i].c;
    else if (a == f[i].b && b == f[i].c)
      e = f[i].a;
    else if (a == f[i].c && b == f[i].a)
      e = f[i].b;
    else if (a == f[i].a && b == f[i].c)
    {
      e = f[i].b;
      reverse = true;
    }
    else if (a == f[i].b && b == f[i].a)
    {
      e = f[i].c;
      reverse = true;
    }
    else if (a == f[i].c && b == f[i].b)
    {
      e = f[i].a;
      reverse = true;
    }

    if (e < std::numeric_limits<unsigned>::max())
    {
      if (!reverse)
      {
        del_list.push_back(i);
        f.push_back(IndexedTri(a, c, e));
        f.push_back(IndexedTri(c, b, e));
      }
      else
      {
        del_list.push_back(i);
        f.push_back(IndexedTri(b, c, e));
        f.push_back(IndexedTri(c, a, e));
      }
    }
  }
}

bool Polyhedron::bisect(const Triangle& tbi, vector<Vector3>& v, vector<IndexedTri>& f, unsigned i)
{
  vector<unsigned> del_list;

  // get tbi's normal and offset
  Vector3 tbi_normal = tbi.calc_normal();
  Real tbi_offset = tbi.calc_offset(tbi_normal);

  // get the three vertices of the triangle
  unsigned a = f[i].a;
  unsigned b = f[i].b;
  unsigned c = f[i].c;

  // get the triangle
  Triangle t(v[a], v[b], v[c]);

  // form three line segments corresponding to the three edges of the triangle
  LineSeg3 eab(t.a, t.b);
  LineSeg3 ebc(t.b, t.c);
  LineSeg3 eca(t.c, t.a);

  // intersect them with plane 
  Real tab = CompGeom::intersect_seg_plane(tbi_normal, tbi_offset, eab);
  Real tbc = CompGeom::intersect_seg_plane(tbi_normal, tbi_offset, ebc);
  Real tca = CompGeom::intersect_seg_plane(tbi_normal, tbi_offset, eca);

  // case 1: segments ab and bc intersect the plane
  if (tab >= (Real) 0.0 && tab <= (Real) 1.0 && 
      tbc >= (Real) 0.0 && tbc <= (Real) 1.0)
  {
    // determine the points of intersection
    Vector3 i_eab = eab.first + (eab.second-eab.first)*tab;
    Vector3 i_ebc = ebc.first + (ebc.second-ebc.first)*tbc;

    // replace the edges, if necessary 
    if (!find_vertex(v, i_eab))
    {
      unsigned vn1 = v.size();
      v.push_back(i_eab);
      replace_edge(v, f, a, b, vn1, del_list);
      std::sort(del_list.begin(), del_list.end());
      std::reverse(del_list.begin(), del_list.end());
      BOOST_FOREACH(unsigned j, del_list)
      {
        f[j] = f.back();
        f.pop_back();
      }
    }

    if (!find_vertex(v, i_ebc))
    {
      unsigned vn2 = v.size();
      v.push_back(i_ebc);
      replace_edge(v, f, b, c, vn2, del_list);
      std::sort(del_list.begin(), del_list.end());
      std::reverse(del_list.begin(), del_list.end());
      BOOST_FOREACH(unsigned j, del_list)
      {
        f[j] = f.back();
        f.pop_back();
      }
    }
  }
  // case 2: segments ab and ca intersect the plane
  if (tab >= (Real) 0.0 && tab <= (Real) 1.0 && 
           tca >= (Real) 0.0 && tca <= (Real) 1.0)
  {
    // determine the points of intersection
    Vector3 i_eab = eab.first + (eab.second-eab.first)*tab;
    Vector3 i_eca = eca.first + (eca.second-eca.first)*tca;

    // replace the edges, if necessary 
    if (!find_vertex(v, i_eab))
    {
      unsigned vn1 = v.size();
      v.push_back(i_eab);
      replace_edge(v, f, a, b, vn1, del_list);
      std::sort(del_list.begin(), del_list.end());
      std::reverse(del_list.begin(), del_list.end());
      BOOST_FOREACH(unsigned j, del_list)
      {
        f[j] = f.back();
        f.pop_back();
      }
    }

    if (!find_vertex(v, i_eca))
    {
      unsigned vn2 = v.size();
      v.push_back(i_eca);
      replace_edge(v, f, c, a, vn2, del_list);
      std::sort(del_list.begin(), del_list.end());
      std::reverse(del_list.begin(), del_list.end());
      BOOST_FOREACH(unsigned j, del_list)
      {
        f[j] = f.back();
        f.pop_back();
      }
    }
  }
  // case 3: segments bc and ca intersect the plane
  if (tbc >= (Real) 0.0 && tbc <= (Real) 1.0 && 
           tca >= (Real) 0.0 && tca <= (Real) 1.0)
  {
    // determine the points of intersection
    Vector3 i_ebc = ebc.first + (ebc.second-ebc.first)*tbc;
    Vector3 i_eca = eca.first + (eca.second-eca.first)*tca;

    // replace the edges, if necessary 
    if (!find_vertex(v, i_ebc))
    {
      unsigned vn1 = v.size();
      v.push_back(i_ebc);
      replace_edge(v, f, b, c, vn1, del_list);
      std::sort(del_list.begin(), del_list.end());
      std::reverse(del_list.begin(), del_list.end());
      BOOST_FOREACH(unsigned j, del_list)
      {
        f[j] = f.back();
        f.pop_back();
      }
    }

    if (!find_vertex(v, i_eca))
    {
      unsigned vn2 = v.size();
      v.push_back(i_eca);
      replace_edge(v, f, c, a, vn2, del_list);
      std::sort(del_list.begin(), del_list.end());
      std::reverse(del_list.begin(), del_list.end());
      BOOST_FOREACH(unsigned j, del_list)
      {
        f[j] = f.back();
        f.pop_back();
      }
    }
  }

  return true;
}

void Polyhedron::slice(const Polyhedron& p1, const Polyhedron& p2, IndexedTriArray& mesh1, IndexedTriArray& mesh2)
{
  // get the vertices of the two meshes -- we will modify these vectors
  vector<Vector3> v1 = p1._mesh.get_vertices();
  vector<Vector3> v2 = p2._mesh.get_vertices();

  // create two new vectors of facets instantiated with polyhedra facets
  vector<IndexedTri> nf1 = p1._mesh.get_facets();
  vector<IndexedTri> nf2 = p2._mesh.get_facets();

  for (unsigned i=0; i< p1._mesh.num_tris(); i++)
  {
    Triangle t1 = p1._mesh.get_triangle(i);

    for (unsigned j=0; j< nf2.size(); j++)
    {
      Triangle t2(v2[nf2[j].a], v2[nf2[j].b], v2[nf2[j].c]);

      if (bisects(t1, t2) && bisects(t2, t1))
        bisect(t1, v2, nf2, j);
    }
  }
  for (unsigned i=0; i< p2._mesh.num_tris(); i++)
  {
    Triangle t1 = p2._mesh.get_triangle(i);

    for (unsigned j=0; j< nf1.size(); j++)
    {
      Triangle t2(v1[nf1[j].a], v1[nf1[j].b], v1[nf1[j].c]);

      if (bisects(t1, t2) && bisects(t2, t1))
        bisect(t1, v1, nf1, j);
    }
  }

  // setup the two meshes
  mesh1 = IndexedTriArray(v1.begin(), v1.end(), nf1.begin(), nf1.end());
  mesh2 = IndexedTriArray(v2.begin(), v2.end(), nf2.begin(), nf2.end());
}

/// Intersects two non-convex polyhedra
/**
 * The running time for this algorithm is O(nm), where n and m are the number
 * of triangles of p1 and p2, respectively.  A faster algorithm can be
 * obtained if p1 and p2 are both convex, but it is not implemented here.
 */
IndexedTriArray Polyhedron::construct_intersection(Polyhedron& p1, Polyhedron& p2)
{
  // slice the two polyhedra to get two new meshes
  IndexedTriArray m1, m2;
  slice(p1, p2, m1, m2);

  // remove triangles from m1 that are outside mesh2 (and shared)
  remove_outside(m1, p2, true);

  // remove triangles from m2 that are outside mesh1 (no shared)
  remove_outside(m2, p1, false);

IndexedTriArray i = IndexedTriArray::merge(m1, m2, NEAR_ZERO).compress_vertices();
Polyhedron p(i);
if (!p.consistent())
  write_poly(p, "isect-fail.wrl");

  // merge m1 and m2 to form a new polyhedron
  return IndexedTriArray::merge(m1, m2, NEAR_ZERO);
}

/// Unions two non-convex polyhedra
IndexedTriArray Polyhedron::construct_union(Polyhedron& p1, Polyhedron& p2)
{
  // slice the two polyhedra to get two new meshes
  IndexedTriArray m1, m2;
  slice(p1._mesh, p2._mesh, m1, m2);

  // remove triangles from m1 that are within p2 (and shared)
  remove_inside(m1, p2, true);

  // remove triangles from m2 that are within p1 (no shared)
  remove_inside(m2, p1, false);

  // merge m1 and m2 to form a new polyhedron
  return IndexedTriArray::merge(m1, m2, NEAR_ZERO);
}

/// Constructs the Boolean difference p1 - p2 
IndexedTriArray Polyhedron::construct_difference(Polyhedron& p1, Polyhedron& p2)
{
  // slice the two polyhedra to get two new meshes
  IndexedTriArray m1, m2;
  slice(p1._mesh, p2._mesh, m1, m2);

  // remove triangles from m1 that are within p2 (and shared)
  remove_inside(m1, p2, true);

  // remove triangles from m2 that are outside p1 (and shared)
  remove_outside(m2, p1, true);

  // invert normals of remaining triangles from m2
  vector<IndexedTri> new_m2_facets = m2.get_facets();
  for (unsigned i=0; i< new_m2_facets.size(); i++)
    std::swap(new_m2_facets[i].b, new_m2_facets[i].c);
  m2 = IndexedTriArray(m2.get_vertices().begin(), m2.get_vertices().end(), 
                       new_m2_facets.begin(), new_m2_facets.end());

  // merge m1 and m2 to form a new polyhedron
  return IndexedTriArray::merge(m1, m2, NEAR_ZERO);
}

/// Sends the specified polyhedron to the desired output stream
std::ostream& Moby::operator<<(std::ostream& out, const Polyhedron& p)
{
  // get the vectors of vertices and facets
  const std::vector<Vector3>& vertices = p.get_vertices();
  const std::vector<IndexedTri>& facets = p.get_facets();
  
  // output vertices 
  out << "vertices: " << std::endl;
  for (unsigned i=0; i< vertices.size(); i++)
    out << "[" << i << "]: " << vertices[i] << std::endl;
  
  // print out the facets
  out << "facets: " << std::endl;
  for (unsigned i=0; i< facets.size(); i++)
  {
    out << "[" << i << "]: ";
    out << facets[i].a << " ";
    out << facets[i].b << " ";
    out << facets[i].c << std::endl;
  }
  
  return out;
}

