/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _GJK_H
#define _GJK_H

#include <iostream>
namespace Moby {

/// An implementation of the GJK algorithm 
class GJK
{
  public:
    static double do_gjk(boost::shared_ptr<const Primitive> A, boost::shared_ptr<const Primitive> B, boost::shared_ptr<const Ravelin::Pose3d> pA, boost::shared_ptr<const Ravelin::Pose3d> pB, Point3d& cpA, Point3d& cpB, unsigned max_iter = 1000);

  private:
    struct SVertex
    {
      Point3d v;      // vA - vB: the vertex in the simplex
      Point3d vA, vB; // vertices from geometry A and B
      std::ostream& output(std::ostream& out) const;

      SVertex() {}
      SVertex(const Point3d& a, const Point3d& b)
      {
        vA = a;
        vB = b;
        v = Ravelin::Pose3d::transform_point(GLOBAL, a) - 
            Ravelin::Pose3d::transform_point(GLOBAL, b);
      }
    };

    class Simplex
    {
      public:
        enum SimplexType { ePoint, eSegment, eTriangle, eTetra };
        SimplexType get_simplex_type() const { return _type; }
        Simplex(const SVertex& p) { _type = ePoint; _v1 = p; }
        void add(const SVertex& p);
        Point3d find_closest_and_simplify();
        Point3d find_closest();
        const SVertex& get_vertex(unsigned i) const;
        unsigned num_vertices() const;
        std::ostream& output(std::ostream& out) const;

      private:
        SimplexType _type;
        SVertex _v1, _v2, _v3, _v4;
    };

}; // end class

} // end namespace

#endif
