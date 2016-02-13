#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/CompGeom.h>

using namespace Ravelin;
using namespace Moby;

int main(int argc, char* argv[])
{
  if (argc < 1)
  {
    std::cerr << "syntax: center <file1.obj> <file2.obj> ... <fileN.obj>" << std::endl;
    return -1;
  }

  std::vector<IndexedTriArray> itas;
  std::list<Triangle> tris;

  // load all files
  for (int i=1; i< argc; i++)
  {
    itas.push_back(IndexedTriArray::read_from_obj(std::string(argv[i])));
    itas.back().get_tris(std::back_inserter(tris), GLOBAL);
  }

  // compute the center-of-mass
  Origin3d com(CompGeom::calc_centroid_3D(tris.begin(), tris.end()));
  Transform3d T;
  T.x = -com;
  std::cout << "center of mass: " << com << std::endl;

  // write the result
  for (int i=1; i< argc; i++)
  {
    // translate the mesh
    IndexedTriArray xlat = itas[i-1].transform(T);

    // write the file to the new filename
    std::string fname = "centered." + std::string(argv[i]);
    std::cout << "writing centered file " << fname << std::endl; 
    IndexedTriArray::write_to_obj(xlat, fname);
  }
}
