/*
 * Converts the geometry specified in an obj file to VRML.
 * Materials, normals, etc. are *not* converted.
 */

#include <string.h>
#include <map>
#include <string>
#include <iostream>
#include <fstream>
#include <vector>
#include <Moby/Types.h>
#include <Moby/IndexedTriArray.h>

using namespace Moby;

void write_vrml_file(const std::string& filename, const IndexedTriArray& mesh)
{
  const unsigned X = 0, Y = 1, Z = 2;

  // get the vertices and facets
  const std::vector<Vector3>& vertices = mesh.get_vertices();
  const std::vector<IndexedTri>& facets = mesh.get_facets();

  // open the file for writing
  std::ofstream out(filename.c_str());

  // write the VRML header
  out << "#VRML V2.0 utf8" << std::endl;
  
  // write the shape
  out << "Shape {" << std::endl;
  out << "  geometry IndexedFaceSet {" << std::endl;
  out << "    coord Coordinate { point [" << std::endl;
  for (unsigned i = 0; i< vertices.size(); i++)
  {
    const Vector3& v = vertices[i];
    out << "      " << v[X] << " " << v[Y] << " " << v[Z] << std::endl;
  }
  out << "    ] }" << std::endl;
  out << "    coordIndex [" << std::endl;
  for (unsigned i=0; i< facets.size(); i++)
  {
    const IndexedTri& t = facets[i];
    out << "      " << t.a << " " << t.b << " ";
    out << t.c << " -1," << std::endl;
  }
  out << "    ] } }" << std::endl;

  out.close();
}
 
int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    std::cerr << "syntax: objwrl <inputfile> <outputfile>" << std::endl << std::endl;
    std::cerr << "Converts an Wavefront OBJ file to a VRML 97 file and vice versa" << std::endl;
    exit(1);
  }

  // determine the type of input file and read in the geometry
  IndexedTriArray mesh;
  std::string infile(argv[1]);
  if (infile.find(".wrl") == infile.size()-4)
    mesh = IndexedTriArray::read_from_wrl(infile);
  else if (infile.find(".obj") == infile.size()-4)
    mesh = IndexedTriArray::read_from_obj(infile);
  else
  {
    std::cerr << "  -- unknown input filename type (no .wrl/.obj extension)" << std::endl;
    exit(1);
  }
  
  // write the output
  std::string outfile(argv[2]);
  if (outfile.find(".wrl") == outfile.size()-4)
    write_vrml_file(outfile, mesh);
  else if (outfile.find(".obj") == outfile.size()-4)
    mesh.write_to_obj(outfile);
  else
  {
    std::cerr << "  -- unknown input filename type (no .wrl/.obj extension)" << std::endl;
    exit(1);
  }

  return 0;
}

