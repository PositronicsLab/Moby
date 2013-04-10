/****************************************************************************
 * Copyright 2010 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _TETRAHEDRON_
#define _TETRAHEDRON_

#include <Moby/Constants.h>
#include <limits>

namespace Moby {

/// An indexed tetrahedron; should be oriented ccw
struct Tetrahedron
{
  Tetrahedron() {  }
  Tetrahedron(const Ravelin::Point3d& v1, const Ravelin::Point3d& v2, const Ravelin::Point3d& v3, const Ravelin::Point3d& v4) { a = v1; b = v2; c = v3; d = v4; }
  double calc_signed_dist(const Ravelin::Point3d& p) const;
  void determine_barycentric_coords(const Ravelin::Point3d& p, double& u, double& v, double& w) const;
  Ravelin::Point3d calc_point(double u, double v, double w) const;
  Ravelin::Point3d calc_centroid() const;
  double calc_volume() const;
  bool outside(const Ravelin::Point3d& p, double tol = NEAR_ZERO) const;

  /// First vertex of the tetrahedron
  Ravelin::Point3d a;

  /// Second vertex of the tetrahedron
  Ravelin::Point3d b;

  /// Third vertex of the tetrahedron
  Ravelin::Point3d c;

  /// Fourth vertex of the tetrahedron
  Ravelin::Point3d d;
}; // end struct

} // end namespace

#endif

