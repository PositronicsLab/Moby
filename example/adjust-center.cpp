#include <Moby/TriangleMeshPrimitive.h>
#include <Moby/CompGeom.h>

using namespace Ravelin;
using namespace Moby;

int main(int argc, char* argv[])
{
  const unsigned X = 0, Y = 1, Z = 2;

  if (argc < 4)
  {
    std::cerr << "syntax: adjust-center <file1.obj> +x +y +z" << std::endl;
    return -1;
  }

  IndexedTriArray itas;
  std::list<Triangle> tris;

  // load file
  itas = IndexedTriArray::read_from_obj(std::string(argv[1]));
  itas.get_tris(std::back_inserter(tris), GLOBAL);

  // get the adjustments
  Origin3d adjust;
  adjust[X] = std::atof(argv[2]);
  adjust[Y] = std::atof(argv[3]);
  adjust[Z] = std::atof(argv[4]);
  Transform3d T;
  T.x = adjust;

  // write the result
  IndexedTriArray xlat = itas.transform(T);

  // write the file to the new filename
  std::string fname = "adjusted." + std::string(argv[1]);
  std::cout << "writing centered file " << fname << std::endl; 
  IndexedTriArray::write_to_obj(xlat, fname);
}
