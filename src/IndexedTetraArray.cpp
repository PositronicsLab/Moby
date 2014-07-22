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
#include <Moby/BoxPrimitive.h>
#include <Moby/CylinderPrimitive.h>
#include <Moby/SpherePrimitive.h>
#include <Moby/ConePrimitive.h>
#include <Moby/DegenerateTriangleException.h>
#include <Moby/XMLTree.h>
#include <Moby/IndexedTetraArray.h>

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
IndexedTetraArray::IndexedTetraArray(shared_ptr<const vector<Origin3d> > vertices, const vector<IndexedTetra>& tetra)
{
  _vertices = vertices;

  // setup tetra 
  shared_ptr<vector<IndexedTetra> > new_tetra(new vector<IndexedTetra>(tetra));
  _tetra = new_tetra;

  // validate indices within range and tetrahedra oriented correctly 
  validate();  
}

/// Initializes a mesh using pointers to vertices and tetra
IndexedTetraArray::IndexedTetraArray(shared_ptr<const vector<Origin3d> > vertices, shared_ptr<const vector<IndexedTetra> > tetra)
{
  _vertices = vertices;
  _tetra = tetra;

  // validate indices within range and tetrahedra oriented correctly
  validate();  
}

/// Copies one mesh to another
IndexedTetraArray& IndexedTetraArray::operator=(const IndexedTetraArray& mesh)
{
  _vertices = mesh._vertices;
  _tetra = mesh._tetra;

  return *this;
}

/// Gets the desired tetrahedron
Tetrahedron IndexedTetraArray::get_tetrahedron(unsigned i, shared_ptr<const Pose3d> P) const
{
  const vector<Origin3d>& vertices = get_vertices();
  const vector<IndexedTetra>& tetra = get_tetra();

  return Tetrahedron(Point3d(vertices[tetra[i].a], P), 
                     Point3d(vertices[tetra[i].b], P), 
                     Point3d(vertices[tetra[i].c], P), 
                     Point3d(vertices[tetra[i].d], P));
}

/// Implements Base::load_from_xml()
/**
 * \note if centering is done, it is done <i<before</i> any transform is applied
 */
void IndexedTetraArray::load_from_xml(shared_ptr<const XMLTree> node, map<string, BasePtr>& id_map)
{
  // load data from the Primitive 
  Base::load_from_xml(node, id_map);

  // check that this node name is correct
  assert(strcasecmp(node->name.c_str(), "TetraMesh") == 0);

  // make sure that this Tetra array has a filename specified
  XMLAttrib* fname_attr = node->get_attrib("filename");
  if (!fname_attr)
  {
    cerr << "IndexedTetraArray::load_from_xml() - trying to load a ";
    cerr << " tetrahedral mesh w/o a filename!" << endl;
    cerr << "  offending node: " << endl << *node << endl;
    return;
  }

  // read the mesh 
  string fname(fname_attr->get_string_value());
  *this = IndexedTetraArray::read_from_tetra(fname);
  
  // see whether to center the mesh
  XMLAttrib* center_attr = node->get_attrib("center");
  if (center_attr && center_attr->get_bool_value())
    center();

  // read in the transform, if specified
  Transform3d T;
  T.source = T.target = GLOBAL;
  XMLAttrib* xlat_attr = node->get_attrib("translation");
  XMLAttrib* rpy_attr = node->get_attrib("rpy"); 
  if (xlat_attr || rpy_attr)
  {
    if (xlat_attr)
      T.x = xlat_attr->get_origin_value();
    if (rpy_attr)
      T.q = xlat_attr->get_rpy_value(); 

    // transform the tetra array
    *this = IndexedTetraArray::transform(T);
  }
}

/// Implements Base::save_to_xml()
void IndexedTetraArray::save_to_xml(XMLTreePtr node, list<shared_ptr<const Base> >& shared_objects) const
{
  // call this parent's save_to_xml() method
  Base::save_to_xml(node, shared_objects);

  // set the name for this node
  node->name = "TetraMesh";

  // make a filename using "this"
  const unsigned MAX_DIGITS = 128;
  char buffer[MAX_DIGITS+1];
  sprintf(buffer, "%p", this);
  string filename = "mesh" + string(buffer) + ".tetra";

  // add the filename as an attribute
  node->attribs.insert(XMLAttrib("filename", filename));

  // do not save the array to the OBJ file if it already exists (which we
  // crudely check for using std::ifstream to avoid OS-specific calls -- note
  // that it is possible that opening a file may fails for other reasons than
  // the file does not exist)
  std::ifstream in(filename.c_str());
  if (in.fail())
    write_to_tetra(filename);
  else
    in.close();
}

/// Centers a tetrahedral mesh (returning a new mesh)
void IndexedTetraArray::center()
{
  // determine the centroid of this array of tetrahedra
  Point3d centroid = Point3d::zero(GLOBAL);
  double total_volume = 0;
  for (unsigned i=0; i< num_tetra(); i++)
  {
    Tetrahedron tet = get_tetrahedron(i, GLOBAL);
    double volume = tet.calc_volume();
    centroid += tet.calc_centroid()*volume;
    total_volume += volume;
  }
  centroid /= total_volume;

  // translate the mesh so that the centroid is at the origin
  for (unsigned i=0; i< _vertices->size(); i++)
    ((Origin3d&) (*_vertices)[i]) -= Origin3d(centroid);
}

/// Verifies that all indices are within range and orients tetrahedra CCW
void IndexedTetraArray::validate()
{
  const vector<Origin3d>& vertices = get_vertices();
  shared_ptr<const vector<IndexedTetra> > tetra_pointer = get_tetra_pointer();
  vector<IndexedTetra>& tetra = *boost::const_pointer_cast<vector<IndexedTetra> >(_tetra);

  unsigned nv = vertices.size();
  for (unsigned i=0; i< tetra.size(); i++)
    if (tetra[i].a >= nv || tetra[i].b >= nv || tetra[i].c >= nv || tetra[i].d >= nv)
      throw InvalidIndexException();

  for (unsigned i=0; i< tetra.size(); i++)
  {
    // get the vertices
    Point3d va(vertices[tetra[i].a], GLOBAL);
    Point3d vb(vertices[tetra[i].b], GLOBAL);
    Point3d vc(vertices[tetra[i].c], GLOBAL);
    Point3d vd(vertices[tetra[i].d], GLOBAL);

    // determine the centroid
    Point3d centroid = va + vb + vc + vd;
    centroid *= 0.25;

    try
    {
      // check triangle abc
      Triangle abc(va, vb, vc); 
      if (abc.calc_signed_dist(centroid) > 0)
      {
        std::swap(tetra[i].b, tetra[i].c);
        std::swap(vb, vc);
      }

      // check triangle bdc
      Triangle bdc(vb, vd, vc);
      if (bdc.calc_signed_dist(centroid) > 0)
      {
        std::swap(tetra[i].d, tetra[i].c);
        std::swap(vd, vc);
      }

      // check triangle dac
      Triangle dac(vd, va, vc);
      if (dac.calc_signed_dist(centroid) > 0)
      {
        std::swap(tetra[i].a, tetra[i].c);
        std::swap(va, vc);
      }

      // check triangle dba
      Triangle dba(vd, vb, va); 
      if (dba.calc_signed_dist(centroid) > 0)
      {
        std::swap(tetra[i].b, tetra[i].a);
        std::swap(vb, va);
      }

      // verify that centroid is now inside tetrahedron
      assert(!Tetrahedron(va, vb, vc, vd).outside(centroid));
    }
    catch (DegenerateTriangleException e)
    {
      // do nothing..  the tetrahedron is degenerate...
    }
  }
}

/// Transforms this mesh to a new mesh
IndexedTetraArray IndexedTetraArray::transform(const Transform3d& T) const
{
  IndexedTetraArray it;

  // can just copy tetra
  it._tetra = _tetra; 

  // need to transform vertices
  shared_ptr<vector<Origin3d> > new_vertices(new vector<Origin3d>(*_vertices));
  it._vertices = new_vertices;
  vector<Origin3d>& vertices = *new_vertices;
  for (unsigned i=0; i< vertices.size(); i++)
    vertices[i] = T.q * vertices[i] + T.x;

  return it;
}

/// Compresses the vertices used in an IndexedTetraArray to create a new mesh 
IndexedTetraArray IndexedTetraArray::compress_vertices() const
{
  // get tetra and vertices
  vector<Origin3d> vertices = get_vertices();
  vector<IndexedTetra> tetra = get_tetra();

  // determine unused vertices
  vector<bool> unused(vertices.size(), true);
  for (unsigned i=0; i< tetra.size(); i++)
  {
    unused[tetra[i].a] = false;
    unused[tetra[i].b] = false;
    unused[tetra[i].c] = false;
    unused[tetra[i].d] = false;
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
    for (unsigned j=0; j< tetra.size(); j++)
    {
      if (tetra[j].a > i)
        tetra[j].a--;
      if (tetra[j].b > i)
        tetra[j].b--;
      if (tetra[j].c > i)
        tetra[j].c--;
      if (tetra[j].d > i)
        tetra[j].d--;
    }
  }

  // create the new mesh
  return IndexedTetraArray(vertices.begin(), vertices.end(), tetra.begin(), tetra.end());
}

/// Writes tetrahedron mesh to a (Moby) .tetra file
void IndexedTetraArray::write_to_tetra(const IndexedTetraArray& mesh, const string& filename)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get vertices and tetra
  const vector<Origin3d>& vertices = mesh.get_vertices();
  const vector<IndexedTetra>& tetra = mesh.get_tetra();

  // open the file for writing
  std::ofstream out(filename.c_str());

  // write the set of vertices
  for (unsigned i=0; i< vertices.size(); i++)
    out << "v " << vertices[i][X] << " " << vertices[i][Y] << " " << vertices[i][Z] << endl;

  // write all tetrahedra
  for (unsigned i=0; i< tetra.size(); i++)
  {
    out << "t";
    out << " " << tetra[i].a;
    out << " " << tetra[i].b;
    out << " " << tetra[i].c;
    out << " " << tetra[i].d;
    out << endl;
  }

  // close the file
  out.close();
}

/// Writes tetrahedron mesh to a Wavefront OBJ file
void IndexedTetraArray::write_to_obj(const IndexedTetraArray& mesh, const string& filename)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get vertices and tetra
  const vector<Origin3d>& vertices = mesh.get_vertices();
  const vector<IndexedTetra>& tetra = mesh.get_tetra();

  // open the file for writing
  std::ofstream out(filename.c_str());

  // write the set of vertices
  for (unsigned i=0; i< vertices.size(); i++)
    out << "v " << vertices[i][X] << " " << vertices[i][Y] << " " << vertices[i][Z] << endl;

  // write all tetrahedra
  for (unsigned i=0; i< tetra.size(); i++)
  {
    out << "f";
    out << " " << tetra[i].a+1;
    out << " " << tetra[i].b+1;
    out << " " << tetra[i].c+1;
    out << endl;
    out << " " << tetra[i].b+1;
    out << " " << tetra[i].d+1;
    out << " " << tetra[i].c+1;
    out << endl;
    out << " " << tetra[i].d+1;
    out << " " << tetra[i].a+1;
    out << " " << tetra[i].c+1;
    out << endl;
    out << " " << tetra[i].d+1;
    out << " " << tetra[i].b+1;
    out << " " << tetra[i].a+1;
    out << endl;
  }

  // close the file
  out.close();
}

/// Reads tetrahedron mesh from a (Moby) .tetra file
IndexedTetraArray IndexedTetraArray::read_from_tetra(const string& filename)
{
  const unsigned BUF_SIZE = 2048;
  char buffer[BUF_SIZE];
  double v1, v2, v3;
  int i1, i2, i3, i4;

  // create arrays for vertices and tetra
  vector<Origin3d> vertices;
  vector<IndexedTetra> tetra;

  // open the file
  std::ifstream in(filename.c_str());
  if (in.fail())
  {
    cerr << "IndexedTetraArray::read_tetra_file() - unable to open ";
    cerr << filename << endl;
    return IndexedTetraArray();
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
        std::string errstr = std::string("Unexpected EOF reached reading vertices from tetra file ") + filename + std::string(" in IndexedTetraArray::read_from_tetra()"); 
        throw std::runtime_error(errstr.c_str());
      }

      // create and add the vertex
      vertices.push_back(Origin3d(v1, v2, v3));
      continue;
    }
  
    // determine whether the read line describes a tetrahedron
    if (id[0] == 't' || id[0] == 'T')
    {
      string t1, t2, t3, t4;
      
      // read in the indices 
      in >> t1;
      in >> t2;
      in >> t3;
      in >> t4;
      if (in.eof())
      {
        std::string errstr = std::string("Unexpected EOF reached reading faces from tetra file ") + filename + std::string(" in IndexedTetraArray::read_from_tetra()"); 
        throw std::runtime_error(errstr.c_str());
      }

      // get the indices  
      i1 = std::atoi(t1.c_str());
      i2 = std::atoi(t2.c_str());
      i3 = std::atoi(t3.c_str());
      i4 = std::atoi(t4.c_str());

      // create the indexed tetrahedron and store the face
      tetra.push_back(IndexedTetra((unsigned) i1, (unsigned) i2, (unsigned) i3, (unsigned) i4));
      continue;
    }

    // unhandled tag, read the rest of the line -- it will be ignored
    in.getline(buffer, BUF_SIZE);
  }

  // close the file
  in.close();
  
  // vertex and face lists should have been read..
  if (vertices.empty() || tetra.empty())
  {
    cerr << "IndexedTetraArray::read_from_tetra() - no vertices and/or tetra ";
    cerr << "read from " << filename << endl;
    vertices.clear();
    tetra.clear();
  }

  // create the indexed tetrahedron array
  return IndexedTetraArray(vertices.begin(), vertices.end(), tetra.begin(), tetra.end());
}


