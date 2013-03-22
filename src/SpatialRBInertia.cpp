/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/Constants.h>
#include <Moby/SMatrix6N.h>
#include <Moby/SpatialRBInertia.h>

using namespace Moby;

/// Default constructor -- constructs a zero inertia matrix
SpatialRBInertia::SpatialRBInertia()
{
  m = (Real) 0.0;
  h = ZEROS_3;
  J = ZEROS_3x3;
}

/// Constructs the 6x6 spatial matrix from the given values 
SpatialRBInertia::SpatialRBInertia(Real m, const Vector3& h, const Matrix3& J)
{
  this->m = m;
  this->h = h;
  this->J = J;
}

/// Copies a spatial matrix to this one 
SpatialRBInertia& SpatialRBInertia::operator=(const SpatialRBInertia& m)
{
  this->m = m.m;
  this->h = m.h;
  this->J = m.J;
  return *this;
}

/// Creates a zero matrix
void SpatialRBInertia::set_zero()
{
  m = (Real) 0.0;
  h = ZEROS_3;
  J = ZEROS_3x3;
}

/// Multiplies the inverse of this spatial matrix by a given vector
SVector6 SpatialRBInertia::inverse_mult(const SVector6& v) const
{
  // compute the skew symmetric version of h
  Matrix3 hx = Matrix3::skew_symmetric(h);

  // compute inverse mass
  Real inv_m = (Real) 1.0/m;

  // compute the components of the matrix
  Matrix3 UR = Matrix3::inverse((hx * hx * inv_m) + J);
  Matrix3 UL = UR * hx * -inv_m;
  Matrix3 LR = Matrix3::transpose(UL);
  Matrix3 LL = ((hx * UL) - IDENTITY_3x3) * inv_m;

  // get the components of v
  Vector3 vtop = v.get_upper();
  Vector3 vbot = v.get_lower();

  // do the arithmetic
  return SVector6(UL*vtop + UR*vbot, LL*vtop + LR*vbot);
}

/// Converts this to a matrix
/**
 * \param output a 36-element array (or larger); on return, contains the
 *        matrix representation of this inertia in column-major format
 */
void SpatialRBInertia::to_matrix(Real output[]) const
{
  const unsigned X = 0, Y = 1, Z = 2;
  const Real HX = h[0], HY = h[1], HZ = h[2];

  // upper left 3x3
  output[0] = 0;       output[6] = HZ;       output[12] = -HY;
  output[1] = -HZ;     output[7] = 0;        output[13] = HX;
  output[2] = HY;      output[8] = -HX;      output[14] = 0;

  // lower left 3x3
  output[3] = J(X,X);  output[9] = J(X,Y);   output[15] = J(X,Z);
  output[4] = J(Y,X);  output[10] = J(Y,Y);  output[16] = J(Y,Z);
  output[5] = J(Z,X);  output[11] = J(Z,Y);  output[17] = J(Z,Z);

  // upper right 3x3
  output[18] = m;       output[24] = 0;        output[30] = 0;
  output[19] = 0;       output[25] = m;        output[31] = 0;
  output[20] = 0;       output[26] = 0;        output[32] = m;

  // lower right 3x3
  output[21] = 0;        output[27] = -HZ;      output[33] = HY;
  output[22] = HZ;       output[28] = 0;        output[34] = -HX;
  output[23] = -HY;      output[29] = HX;       output[35] = 0;
}

/// Multiplies this matrix by a vector and returns the result in a new vector
SVector6 SpatialRBInertia::operator*(const SVector6& v) const
{
  // get necessary components of v
  Vector3 vtop = v.get_upper();
  Vector3 vbot = v.get_lower();

  // do some precomputation
  Matrix3 hX = Matrix3::skew_symmetric(h);

  // compute top part of result
  Vector3 rtop = (vbot * m) - (hX * vtop);
  Vector3 rbot = (J * vtop) + (hX * vbot);

  return SVector6(rtop, rbot); 
}

/// Multiplies this matrix by a scalar in place
SpatialRBInertia& SpatialRBInertia::operator*=(Real scalar)
{
  m *= scalar;
  h *= scalar;
  J *= scalar;
  return *this;
}

/// Returns the negation of this matrix
SpatialRBInertia SpatialRBInertia::operator-() const
{
  SpatialRBInertia result;
  result.m = -this->m;
  result.h = -this->h;
  result.J = -this->J;
  return result;
}

/// Adds two spatial matrices
SpatialRBInertia SpatialRBInertia::operator+(const SpatialRBInertia& m) const
{
  SpatialRBInertia result;
  result.m = this->m + m.m;
  result.h = this->h + m.h;
  result.J = this->J + m.J;
  return result;
}

/// Subtracts two spatial matrices
SpatialRBInertia SpatialRBInertia::operator-(const SpatialRBInertia& m) const
{
  SpatialRBInertia result;
  result.m = this->m - m.m;
  result.h = this->h - m.h;
  result.J = this->J - m.J;
  return result;
}

/// Adds m to this in place
SpatialRBInertia& SpatialRBInertia::operator+=(const SpatialRBInertia& m)
{
  this->m += m.m;
  this->h += m.h;
  this->J += m.J;
  return *this;
}

/// Subtracts m from this in place
SpatialRBInertia& SpatialRBInertia::operator-=(const SpatialRBInertia& m)
{
  this->m -= m.m;
  this->h -= m.h;
  this->J -= m.J;
  return *this;
}

/// Multiplies a spatial matrix by a spatial matrix and returns the result in a spatial matrix
SMatrix6N& SpatialRBInertia::mult(const SMatrix6N& m, SMatrix6N& result) const
{
  const unsigned SPATIAL_DIM = 6;

  // get number of columns of m
  const unsigned NCOLS = m.columns();

  // resize the result
  result.resize(SPATIAL_DIM, NCOLS);

  // look for empty result
  if (NCOLS == 0)
  {
    result.set_zero();
    return result;
  }

  // compute the skew symmetric matrix corresponding to h
  Matrix3 hX = Matrix3::skew_symmetric(this->h);

  // carry out multiplication one column at a time
  for (unsigned i=0; i< NCOLS; i++)
  {
    SVector6 v = m.get_column(i);
    Vector3 vtop = v.get_upper();
    Vector3 vbot = v.get_lower();
    v = SVector6((vbot * this->m)-(hX * vtop), (this->J * vtop)+(hX * vbot));
    result.set_column(i, v);
  } 

  return result;
}

/// Outputs this matrix to the stream
std::ostream& Moby::operator<<(std::ostream& out, const SpatialRBInertia& m) 
{
  out << "spatial rigid body mass=" << m.m << " h = " << m.h << " J = ";
  out << std::endl << m.J;
   
  return out;
}

