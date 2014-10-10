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
}; // end class

} // end namespace

#endif
