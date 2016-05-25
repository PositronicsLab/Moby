/****************************************************************************
 * Copyright 2015 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _LP_H
#define _LP_H

#include <Ravelin/NonsquareMatrixException.h>
#include <Ravelin/VectorNd.h>
#include <Ravelin/MatrixNd.h>

namespace Moby {

/// Determines closest points/common points between two polytopes 
class LP 
{
  public:
    static bool lp_simplex(const Ravelin::VectorNd& c, const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, const Ravelin::MatrixNd& A, const Ravelin::VectorNd& b, const Ravelin::VectorNd& l, const Ravelin::VectorNd& u, Ravelin::VectorNd& x);
    static bool lp_seidel(const Ravelin::MatrixNd& A, const Ravelin::VectorNd& b, const Ravelin::VectorNd& c, const Ravelin::VectorNd& l, const Ravelin::VectorNd& u, Ravelin::VectorNd& x);

  private:
    static Ravelin::VectorNd& insert_component(const Ravelin::VectorNd& x, unsigned k, Ravelin::VectorNd& xn);
    static Ravelin::VectorNd& remove_component(const Ravelin::VectorNd& x, unsigned k, Ravelin::VectorNd& xn);
    static double finitize(double x);
}; // end class

} // end namespace

#endif
