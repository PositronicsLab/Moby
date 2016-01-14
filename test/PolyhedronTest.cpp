#include <boost/foreach.hpp>
#include <boost/shared_ptr.hpp>
#include <gtest/gtest.h>
#include <Moby/Polyhedron.h>
#include <Moby/BoxPrimitive.h>
#include <Moby/CompGeom.h>
#include <Moby/TessellatedPolyhedron.h>
#include <Moby/Log.h>


using namespace Ravelin;
using namespace Moby;

  int main() 
  {
    std::vector<Origin3d> v(8);
      v[0] = Origin3d(-1.0, -1.0, -1.0);
      v[1] = Origin3d(-1.0, -1.0, +1.0);
      v[2] = Origin3d(-1.0, +1.0, -1.0);
      v[3] = Origin3d(-1.0, +1.0, +1.0);
      v[4] = Origin3d(+1.0, -1.0, -1.0);
      v[5] = Origin3d(+1.0, -1.0, +1.0);
      v[6] = Origin3d(+1.0, +1.0, -1.0);
      v[7] = Origin3d(+1.0, +1.0, +1.0);


    TessellatedPolyhedronPtr p_tess = CompGeom::calc_convex_hull(v.begin(), v.end());
    Polyhedron p_test = p_tess->to_polyhedron();
    std::ofstream out("box.wrl");
    Polyhedron::to_vrml(out, p_test, Origin3d(1,1,1), true);
    out.close();
 }
