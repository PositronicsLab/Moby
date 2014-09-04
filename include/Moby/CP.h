/****************************************************************************
 * Copyright 2014 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _CP_H
#define _CP_H

#include <Ravelin/NonsquareMatrixException.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/MatrixNd.h>
#include <Ravelin/Pose3d.h>
#include <Moby/Types.h>
#include <Moby/Primitive.h>

namespace Moby {

class PolyhedralPrimitive;

/// Determines closest points/common points between two polytopes 
class CP 
{
  public:
    static double find_cpoint(boost::shared_ptr<const PolyhedralPrimitive> A, boost::shared_ptr<const PolyhedralPrimitive> B, boost::shared_ptr<const Ravelin::Pose3d> pA, boost::shared_ptr<const Ravelin::Pose3d> pB, Point3d& cpA, Point3d& cpB);

  private:
    static bool lp_seidel(const Ravelin::MatrixNd& A, const Ravelin::VectorNd& b, const Ravelin::VectorNd& c, const Ravelin::VectorNd& l, const Ravelin::VectorNd& u, Ravelin::VectorNd& x);
    static Ravelin::VectorNd& insert_component(const Ravelin::VectorNd& x, unsigned k, Ravelin::VectorNd& xn);
    static Ravelin::VectorNd& remove_component(const Ravelin::VectorNd& x, unsigned k, Ravelin::VectorNd& xn);
    static double finitize(double x);
    static void update_z(const Ravelin::MatrixNd& J, const Ravelin::VectorNd& d, unsigned iq, Ravelin::VectorNd& z);
    static void update_r(const Ravelin::MatrixNd& R, const Ravelin::VectorNd& d, unsigned iq, Ravelin::VectorNd& r);
    static bool add_constraint(Ravelin::MatrixNd& R, Ravelin::MatrixNd& J, Ravelin::VectorNd& d, unsigned& iq, double& rnorm);
    static void delete_constraint(Ravelin::MatrixNd& R, Ravelin::MatrixNd& J, std::vector<unsigned>& A, Ravelin::VectorNd& u, unsigned n, unsigned p, unsigned& iq, unsigned l);
    static double distance(double a, double b);
    static bool qp_strict_convex(Ravelin::MatrixNd& G, const Ravelin::VectorNd& c, const Ravelin::MatrixNd& A, const Ravelin::VectorNd& b, const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& x, bool G_factored);

    // variables for strictly convex QP solver
    static Ravelin::MatrixNd _R, _J;
    static Ravelin::VectorNd _d, _z, _s, _r, _np, _u, _xold, _uold;
    static std::vector<unsigned> _a, _aold, _iai;
    static std::vector<bool> _iaexcl;
}; // end class

} // end namespace

#endif
