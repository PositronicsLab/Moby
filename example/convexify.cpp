#include <Moby/Polyhedron.h>
#include <Moby/CompGeom.h>
#include <fstream>

using namespace Ravelin;
using namespace Moby;

int main(int argc, char* argv[])
{
  if (argc != 3)
  {
    std::cerr << "syntax: convexify <input> <output>" << std::endl;
    std::cerr << std::endl;
    std::cerr << "convexify takes the description of a 3D geometry from the input file (a Wavefront OBJ" << std::endl;
    std::cerr << "file) and constructs a new Wavefront 3D file of the convex hull of that file." << std::endl;
    return -1;
  }

  // read in the file
  IndexedTriArray mesh = IndexedTriArray::read_from_obj(std::string(argv[1])); 

  // get the vertices
  const std::vector<Origin3d>& vertices = mesh.get_vertices();

  // compute the convex hull
  PolyhedronPtr p = CompGeom::calc_convex_hull(vertices.begin(), vertices.end());    

  // write the resulting mesh to the output file
  if (!p)
  {
    std::cerr << "convexify() - fatal error computing convex hull!  (sorry I can't tell you more!" << std::endl;
    return -1;
  }

  // write the file
  p->get_mesh().write_to_obj(std::string(argv[2]));
}

