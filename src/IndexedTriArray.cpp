/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#include <queue>
#include <iostream>
#include <fstream>
#include <Moby/Constants.h>
#include <Moby/CompGeom.h>
#include <Moby/InvalidIndexException.h>
#include <Moby/DegenerateTriangleException.h>
#include <Moby/IndexedTriArray.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/ConePrimitive.h>

using namespace Ravelin;
using namespace Moby;
using std::cerr;
using std::endl;
using std::map;
using std::vector;
using std::list;
using std::string;
using boost::shared_ptr;

/// Initializes a mesh using a pointer to vertices
IndexedTriArray::IndexedTriArray(shared_ptr<const vector<Origin3d> > vertices, const vector<IndexedTri>& facets)
{
  _vertices = vertices;

  // setup facets 
  shared_ptr<vector<IndexedTri> > new_facets(new vector<IndexedTri>(facets));
  _facets = new_facets;

  // validate indices within range 
  validate();  

  // calculate incident facets and determine coplanar features
  calc_incident_facets();
  determine_coplanar_features();
}

/// Initializes a mesh using pointers to vertices and facets
IndexedTriArray::IndexedTriArray(shared_ptr<const vector<Origin3d> > vertices, shared_ptr<const vector<IndexedTri> > facets)
{
  _vertices = vertices;
  _facets = facets;

  // validate indices within range
  validate();  

  // calculate incident facets and determine coplanar features
  calc_incident_facets();
  determine_coplanar_features();
}

/// Calculates the volume integrals of this primitive as a triangle mesh
void IndexedTriArray::calc_volume_ints(double volume_ints[10]) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const double f_1_6th = 1.0/6.0;
  const double f_1_24th = 1.0/24.0;
  const double f_1_60th = 1.0/60.0;
  const double f_1_120th = 1.0/120.0;

  // order:  1, x, y, z, x^2, y^2, z^2, xy, yz, zx

  // get the facets and vertices in the mesh
  const std::vector<Origin3d>& verts = get_vertices(); 
  const std::vector<IndexedTri>& facets = get_facets(); 

  // process all triangles in the mesh
  for (unsigned i = 0; i < facets.size(); i++)
  {
    // get the three vertices of the triangle
    const Origin3d& v0 = verts[facets[i].a];
    const Origin3d& v1 = verts[facets[i].b];
    const Origin3d& v2 = verts[facets[i].c];

    // get the cross-product of the edges
    Vector3d cross = Vector3d::cross(Vector3d(v1 - v0, GLOBAL), Vector3d(v2 - v0, GLOBAL));

    // compute integral terms
    double tmp0, tmp1, tmp2;
    double F1x, F2x, F3x, G0x, G1x, G2x;
    tmp0 = v0[X] + v1[X];
    F1x = tmp0 + v2[X];
    tmp1 = v0[X]*v0[X];
    tmp2 = tmp1 + v1[X]*tmp0;
    F2x = tmp2 + v2[X]*F1x;
    F3x = v0[X]*tmp1 + v1[X]*tmp2 + v2[X]*F2x;
    G0x = F2x + v0[X]*(F1x + v0[X]);
    G1x = F2x + v1[X]*(F1x + v1[X]);
    G2x = F2x + v2[X]*(F1x + v2[X]);

    double F1y, F2y, F3y, G0y, G1y, G2y;
    tmp0 = v0[Y] + v1[Y];
    F1y = tmp0 + v2[Y];
    tmp1 = v0[Y]*v0[Y];
    tmp2 = tmp1 + v1[Y]*tmp0;
    F2y = tmp2 + v2[Y]*F1y;
    F3y = v0[Y]*tmp1 + v1[Y]*tmp2 + v2[Y]*F2y;
    G0y = F2y + v0[Y]*(F1y + v0[Y]);
    G1y = F2y + v1[Y]*(F1y + v1[Y]);
    G2y = F2y + v2[Y]*(F1y + v2[Y]);

    double F1z, F2z, F3z, G0z, G1z, G2z;
    tmp0 = v0[Z] + v1[Z];
    F1z = tmp0 + v2[Z];
    tmp1 = v0[Z]*v0[Z];
    tmp2 = tmp1 + v1[Z]*tmp0;
    F2z = tmp2 + v2[Z]*F1z;
    F3z = v0[Z]*tmp1 + v1[Z]*tmp2 + v2[Z]*F2z;
    G0z = F2z + v0[Z]*(F1z + v0[Z]);
    G1z = F2z + v1[Z]*(F1z + v1[Z]);
    G2z = F2z + v2[Z]*(F1z + v2[Z]);

    // update integrals
    volume_ints[0] += cross[X]*F1x;
    volume_ints[1] += cross[X]*F2x;
    volume_ints[2] += cross[Y]*F2y;
    volume_ints[3] += cross[Z]*F2z;
    volume_ints[4] += cross[X]*F3x;
    volume_ints[5] += cross[Y]*F3y;
    volume_ints[6] += cross[Z]*F3z;
    volume_ints[7] += cross[X]*(v0[Y]*G0x + v1[Y]*G1x + v2[Y]*G2x);
    volume_ints[8] += cross[Y]*(v0[Z]*G0y + v1[Z]*G1y + v2[Z]*G2y);
    volume_ints[9] += cross[Z]*(v0[X]*G0z + v1[X]*G1z + v2[X]*G2z);
  }

  volume_ints[0] *= f_1_6th;
  volume_ints[1] *= f_1_24th;
  volume_ints[2] *= f_1_24th;
  volume_ints[3] *= f_1_24th;
  volume_ints[4] *= f_1_60th;
  volume_ints[5] *= f_1_60th;
  volume_ints[6] *= f_1_60th;
  volume_ints[7] *= f_1_120th;
  volume_ints[8] *= f_1_120th;
  volume_ints[9] *= f_1_120th;
}

/// Copies one mesh to another
IndexedTriArray& IndexedTriArray::operator=(const IndexedTriArray& mesh)
{
  _vertices = mesh._vertices;
  _facets = mesh._facets;
  _incident_facets = mesh._incident_facets;
  _coplanar_verts = mesh._coplanar_verts;
  _coplanar_edges = mesh._coplanar_edges;

  return *this;
}

/// Determines a map from vertices to facet indices
vector<list<unsigned> > IndexedTriArray::determine_vertex_facet_map() const
{
  vector<list<unsigned> > m(_vertices->size());
  const vector<IndexedTri>& facets = *_facets;

  for (unsigned i=0; i< facets.size(); i++)
  {
    m[facets[i].a].push_back(i);
    m[facets[i].b].push_back(i);
    m[facets[i].c].push_back(i);
  }

  return m;
}

/// Determines a map from vertices to edges
/**
 * \return a map from vertices to edges; each element in the map repesents
 *         the second vertex of the edge (the key is the first)
 */
vector<list<unsigned> > IndexedTriArray::determine_vertex_edge_map() const
{
  vector<list<unsigned> > m(_vertices->size());
  const vector<IndexedTri>& facets = *_facets;
  for (unsigned i=0; i< facets.size(); i++)
  {
    m[facets[i].a].push_back(facets[i].b);
    m[facets[i].a].push_back(facets[i].c);
    m[facets[i].b].push_back(facets[i].a);
    m[facets[i].b].push_back(facets[i].c);
    m[facets[i].c].push_back(facets[i].a);
    m[facets[i].c].push_back(facets[i].b);
  }

  return m;
}

/// Determines a map from edges to facet indices
map<sorted_pair<unsigned>, list<unsigned> > IndexedTriArray::determine_edge_facet_map() const
{
  map<sorted_pair<unsigned>, list<unsigned> > m;
  const vector<IndexedTri>& facets = *_facets;
  for (unsigned i=0; i< facets.size(); i++)
  {
    m[make_sorted_pair(facets[i].a, facets[i].b)].push_back(i);
    m[make_sorted_pair(facets[i].b, facets[i].c)].push_back(i);
    m[make_sorted_pair(facets[i].c, facets[i].a)].push_back(i);
  }

  return m;
}

/// Gets the i'th triangle from this mesh
Triangle IndexedTriArray::get_triangle(unsigned i, shared_ptr<const Pose3d> P) const
{
  if (i >= _facets->size())
    throw InvalidIndexException();

  // get the facets and vertices as objects
  const vector<Origin3d>& vertices = get_vertices();
  const vector<IndexedTri>& facets = get_facets();

  return Triangle(Point3d(vertices[facets[i].a], P), Point3d(vertices[facets[i].b], P), Point3d(vertices[facets[i].c], P));
}

/// Verifies that all indices are within range
void IndexedTriArray::validate() const
{
  const vector<Origin3d>& vertices = get_vertices();
  const vector<IndexedTri>& facets = get_facets();

  unsigned nv = vertices.size();
  for (unsigned i=0; i< facets.size(); i++)
    if (facets[i].a >= nv || facets[i].b >= nv || facets[i].c >= nv)
      throw InvalidIndexException();
}

/// Merges two meshes together to create a new mesh
/**
 * \note this method runs in time O(nm)
 */
IndexedTriArray IndexedTriArray::merge(const IndexedTriArray& mesh1, const IndexedTriArray& mesh2, double equal_tol)
{
  // if one mesh is empty, just return the other
  if (!mesh1._vertices)
    return mesh2;
  else if (!mesh2._vertices)
    return mesh1;

  // get the vertices and facets from the first mesh
  vector<Origin3d> vertices = mesh1.get_vertices();
  vector<IndexedTri> facets = mesh1.get_facets();

  // get the vertices and facets from the first mesh
  const vector<Origin3d>& vertices2 = mesh2.get_vertices();
  vector<IndexedTri> facets2 = mesh2.get_facets();

  // setup a mapping for vertices in the second mesh
  map<unsigned, unsigned> mapping;

  // determine the mapping for vertices in the second mesh
  for (unsigned i=0; i< vertices2.size(); i++)
  {
    double closest_dist = (vertices.front() - vertices2[i]).norm();
    unsigned closest = 0;
    for (unsigned j=1; j< vertices.size(); j++)
    {
      double dist = (vertices[j] - vertices2[i]).norm();
      if (dist < closest_dist)
      {
        closest_dist = dist;
        closest = j;
      }
    }
  
    // see whether the closest is sufficiently close
    if (closest_dist <= equal_tol)
    {
      mapping[i] = closest;
    }
    else
    {
      mapping[i] = vertices.size();
      vertices.push_back(vertices2[i]);
    }
  }

  // alter indices in second mesh
  for (unsigned i=0; i< facets2.size(); i++)
  {
    facets2[i].a = mapping[facets2[i].a];
    facets2[i].b = mapping[facets2[i].b];
    facets2[i].c = mapping[facets2[i].c];
  }

  // combine facets together
  facets.insert(facets.end(), facets2.begin(), facets2.end());

  return IndexedTriArray(vertices.begin(), vertices.end(), facets.begin(), facets.end());
}

/// Transforms this mesh to a new mesh
IndexedTriArray IndexedTriArray::transform(const Transform3d& T) const
{
  IndexedTriArray it;

  // can just copy facets and incident tris
  it._facets = _facets; 
  it._incident_facets = _incident_facets;

  // need to transform vertices
  shared_ptr<vector<Origin3d> > new_vertices(new vector<Origin3d>(*_vertices));
  it._vertices = new_vertices;
  vector<Origin3d>& vertices = *new_vertices;
  for (unsigned i=0; i< vertices.size(); i++)
    vertices[i] = T.q * vertices[i] + T.x;

  return it;
}

/// Calculates facets incident to a vertex
void IndexedTriArray::calc_incident_facets()
{
  // create the array
  shared_ptr<vector<list<unsigned> > > new_incident_facets(new vector<list<unsigned> >(_vertices->size()));
  _incident_facets = new_incident_facets;
  vector<list<unsigned> >& incident_facets = *new_incident_facets;

  // get facets and vertices
  const vector<IndexedTri>& facets = get_facets();

  // set it up
  for (unsigned i=0; i< facets.size(); i++)
  {
    incident_facets[facets[i].a].push_back(i);
    incident_facets[facets[i].b].push_back(i);
    incident_facets[facets[i].c].push_back(i);
  }
}

/// Compresses the vertices used in an IndexedTriArray to create a new mesh 
IndexedTriArray IndexedTriArray::compress_vertices() const
{
  // get facets and vertices
  vector<Origin3d> vertices = get_vertices();
  vector<IndexedTri> facets = get_facets();

  // determine unused vertices
  vector<bool> unused(vertices.size(), true);
  for (unsigned i=0; i< facets.size(); i++)
  {
    unused[facets[i].a] = false;
    unused[facets[i].b] = false;
    unused[facets[i].c] = false;
  }

  // the tricky part: remove unused vertices
  for (unsigned i=unused.size()-1; i < unused.size(); i--)
  {
    // is vertex used?  skip it
    if (!unused[i])
      continue;

    // remove the vertex from the array of vertices
    vertices.erase(vertices.begin() + i);

    // all vertex indices after this one should be decremented
    for (unsigned j=0; j< facets.size(); j++)
    {
      if (facets[j].a > i)
        facets[j].a--;
      if (facets[j].b > i)
        facets[j].b--;
      if (facets[j].c > i)
        facets[j].c--;
    }
  }

  // create the new mesh
  return IndexedTriArray(vertices.begin(), vertices.end(), facets.begin(), facets.end());
}

/// Writes triangle mesh to a Wavefront OBJ file
void IndexedTriArray::write_to_obj(const IndexedTriArray& mesh, const string& filename)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get vertices and facets
  const vector<Origin3d>& vertices = mesh.get_vertices();
  const vector<IndexedTri>& facets = mesh.get_facets();

  // open the file for writing
  std::ofstream out(filename.c_str());

  // write the set of vertices
  for (unsigned i=0; i< vertices.size(); i++)
    out << "v " << vertices[i][X] << " " << vertices[i][Y] << " " << vertices[i][Z] << endl;

  // write all facets
  for (unsigned i=0; i< facets.size(); i++)
  {
    out << "f";
    out << " " << facets[i].a+1;
    out << " " << facets[i].b+1;
    out << " " << facets[i].c+1;
    out << endl;
  }

  // close the file
  out.close();
}

/// Reads triangle mesh from a Wavefront OBJ file
IndexedTriArray IndexedTriArray::read_from_obj(const string& filename)
{
  const unsigned BUF_SIZE = 2048;
  char buffer[BUF_SIZE];
  double v1, v2, v3;
  int i1, i2, i3;

  // create arrays for vertices and facets
  vector<Origin3d> vertices;
  vector<IndexedTri> facets;

  // open the file
  std::ifstream in(filename.c_str());
  if (in.fail())
  {
    cerr << "IndexedTriArray::read_obj_file() - unable to open ";
    cerr << filename << endl;
    return IndexedTriArray();
  }

  // read in the file
  while (true)
  {
    // read in the line identifier
    string id;
    in >> id;
    if (in.eof())
      break;
    
    // determine whether the line describes a vertex
    if (id == std::string("v") || id == std::string("V"))
    {
      // read in the vertex and the rest of the line
      in >> v1;
      in >> v2;
      in >> v3;
      if (in.eof())
      {
        std::string errstr = std::string("Unexpected EOF reached reading vertices from OBJ file ") + filename + std::string(" in IndexedTriArray::read_from_obj()"); 
        throw std::runtime_error(errstr.c_str());
      }

      // create and add the vertex
      vertices.push_back(Origin3d(v1, v2, v3));
      continue;
    }
  
    // determine whether the read line describes a face
    if (id[0] == 'f' || id[0] == 'F')
    {
      string f1, f2, f3;
      size_t idx;
      
      // read in the indices 
      in >> f1;
      in >> f2;
      in >> f3;
      if (in.eof())
      {
        std::string errstr = std::string("Unexpected EOF reached reading faces from OBJ file ") + filename + std::string(" in IndexedTriArray::read_from_obj()"); 
        throw std::runtime_error(errstr.c_str());
      }

      // strip out associated texture and normal vertices for f1
      idx = f1.find_first_of('/');
      if (idx != string::npos)
        f1.erase(idx, string::npos);
      
      // strip out associated texture and normal vertices for f2
      idx = f2.find_first_of('/');
      if (idx != string::npos)
        f2.erase(idx, string::npos);
      
      // strip out associated texture and normal vertices for f3
      idx = f3.find_first_of('/');
      if (idx != string::npos)
        f3.erase(idx, string::npos);
    
      // get the indices  
      i1 = std::atoi(f1.c_str());
      i2 = std::atoi(f2.c_str());
      i3 = std::atoi(f3.c_str());

      // if the indices are negative, translate them to positive; otherwise, decrement them
      if (i1 < 0)
        i1 = vertices.size() + i1;
      else
        i1--;
      if (i2 < 0)
        i2 = vertices.size() + i2;
      else
        i2--;
      if (i3 < 0)
        i3 = vertices.size() + i3;
      else
        i3--;

      // if triangle is degenerate, do not store it
      Triangle tri(Point3d(vertices[i1], GLOBAL), 
                   Point3d(vertices[i2], GLOBAL),
                   Point3d(vertices[i3], GLOBAL));
      if (tri.calc_area() < NEAR_ZERO)
        continue;

      // create the indexed triangle and store the face
      facets.push_back(IndexedTri((unsigned) i1, (unsigned) i2, (unsigned) i3));
      continue;
    }

    // unhandled tag, read the rest of the line -- it will be ignored
    in.getline(buffer, BUF_SIZE);
  }

  // close the file
  in.close();
  
  // vertex and face lists should have been read..
  if (vertices.empty() || facets.empty())
  {
    cerr << "IndexedTriArray::read_obj_file() - no vertices and/or facets ";
    cerr << "read from " << filename << endl;
    vertices.clear();
    facets.clear();
  }

  // create the indexed triangle array
  return IndexedTriArray(vertices.begin(), vertices.end(), facets.begin(), facets.end());
}

/// Method exists b/c CompGeom cannot be included from IndexedTriArray.h
bool IndexedTriArray::query_intersect_tri_tri(const Triangle& t1, const Triangle& t2)
{
  return CompGeom::query_intersect_tri_tri(t1, t2);
}

/// Determines coplanar vertices and edges of the mesh
void IndexedTriArray::determine_coplanar_features()
{
  // compute maps of vertices to facets and edges to facets in the mesh
  vector<list<unsigned> > v_f_map = determine_vertex_facet_map();
  map<sorted_pair<unsigned>, list<unsigned> > e_f_map = determine_edge_facet_map();

  // get the vertices and facets
  const vector<IndexedTri>& f = get_facets();
  const vector<Origin3d>& v = get_vertices();

  // init lists which will hold the temporary result
  list<sorted_pair<unsigned> > cp_edges;
  list<unsigned> cp_verts;

  // construct the list of coplanar vertices
  for (unsigned i=0; i< _vertices->size(); i++)
  {
    list<unsigned>::const_iterator liter = v_f_map[i].begin();
    if (liter == v_f_map[i].end())
      continue;

    // get the normal to the first face
    const IndexedTri& fi = f[*liter++];
    Vector3d fn = Triangle(Point3d(v[fi.a], GLOBAL), Point3d(v[fi.b], GLOBAL), Point3d(v[fi.c], GLOBAL)).calc_normal();

    // iterate
    while (liter != v_f_map[i].end()) 
    {
      // get the normal to the second face
      const IndexedTri& fli = f[*liter++];      
      Vector3d fn2 = Triangle(Point3d(v[fli.a], GLOBAL), Point3d(v[fli.b], GLOBAL), Point3d(v[fli.c], GLOBAL)).calc_normal();

      // check angle between the two normals -- note: we allow the normals to
      // be completely opposed
      if (std::fabs(std::fabs(fn.dot(fn2)) - 1.0) > NEAR_ZERO)
        goto skip1;
    }

    // if here, all faces are coplanar
    cp_verts.push_back(i);

    // this jump point skips the insertion into the list
    skip1:;
  }

  // construct the list of coplanar edges
  for (map<sorted_pair<unsigned>, list<unsigned> >::const_iterator i = e_f_map.begin(); i != e_f_map.end(); i++)
  {
    assert(!i->second.empty());
    list<unsigned>::const_iterator liter = i->second.begin();

    // get the normal to the first face
    const IndexedTri& fi = f[*liter++];
    Vector3d fn = Triangle(Point3d(v[fi.a], GLOBAL), Point3d(v[fi.b], GLOBAL), Point3d(v[fi.c], GLOBAL)).calc_normal();

    // iterate
    while (liter != i->second.end()) 
    {
      // get the normal to the second face
      const IndexedTri& fli = f[*liter++];      
      Vector3d fn2 = Triangle(Point3d(v[fli.a], GLOBAL), Point3d(v[fli.b], GLOBAL), Point3d(v[fli.c], GLOBAL)).calc_normal();

      // check angle between the two normals -- note: we allow the normals to
      // be completely opposed
      if (std::fabs(std::fabs(fn.dot(fn2)) - 1.0) > NEAR_ZERO)
        goto skip2;
    }

    // if here, all faces are coplanar
    cp_edges.push_back(i->first);

    // this jump point skips the insertion into the list
    skip2:;
  }  

  // sort the lists
  cp_verts.sort();
  cp_edges.sort();

  // setup the vectors
  _coplanar_verts = vector<unsigned>(cp_verts.begin(), cp_verts.end());
  _coplanar_edges = vector<sorted_pair<unsigned> >(cp_edges.begin(), cp_edges.end());
}

