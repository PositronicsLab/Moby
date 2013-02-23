/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <Moby/Constants.h>
#include <Moby/SMatrix6N.h>
#include <Moby/MatrixN.h>
#include <Moby/SpatialABInertia.h>

using namespace Moby;

/// Default constructor -- constructs a zero inertia matrix
SpatialABInertia::SpatialABInertia()
{
  M = ZEROS_3x3;
  H = ZEROS_3x3;
  J = ZEROS_3x3;
}

/// Constructs a spatial AB inertia from the given MatrixN object
SpatialABInertia::SpatialABInertia(const MatrixN& m)
{
  m.get_sub_mat(0, 3, 3, 6, M);
  m.get_sub_mat(3, 6, 3, 6, H);
  m.get_sub_mat(3, 6, 0, 3, J);
}

/// Constructs the spatial AB inertia from the given values 
SpatialABInertia::SpatialABInertia(const Matrix3& M, const Matrix3& H, const Matrix3& J)
{
  this->M = M;
  this->H = H;
  this->J = J;
}

/// Copies a spatial AB inertia to this one 
SpatialABInertia& SpatialABInertia::operator=(const SpatialABInertia& m)
{
  this->M = m.M;
  this->H = m.H;
  this->J = m.J;
  return *this;
}

/// Copies a spatial RB inertia to this one 
SpatialABInertia& SpatialABInertia::operator=(const SpatialRBInertia& m)
{
  this->M = IDENTITY_3x3 * m.m;
  this->H = Matrix3::skew_symmetric(m.h);
  this->J = m.J;
  return *this;
}

/// Creates a zero matrix
void SpatialABInertia::set_zero()
{
  M = H = J = ZEROS_3x3;
}

/// Multiplies this matrix by a vector and returns the result in a new vector
SVector6 SpatialABInertia::operator*(const SVector6& v) const
{
  // get necessary components of v
  Vector3 vtop = v.get_upper();
  Vector3 vbot = v.get_lower();

  // do some precomputation
  Matrix3 HT = Matrix3::transpose(H);

  // compute top part of result
  Vector3 rtop = HT * vtop + (M * vbot);
  Vector3 rbot = (J * vtop) + (H * vbot);

  return SVector6(rtop, rbot); 
}

/// Multiplies this matrix by a scalar in place
SpatialABInertia& SpatialABInertia::operator*=(Real scalar)
{
  M *= scalar;
  H *= scalar;
  J *= scalar;
  return *this;
}

/// Returns the negation of this matrix
SpatialABInertia SpatialABInertia::operator-() const
{
  SpatialABInertia result;
  result.M = -this->M;
  result.H = -this->H;
  result.J = -this->J;
  return result;
}

/// Adds a spatial articulated body inertia and a spatial rigid body inertia 
SpatialABInertia SpatialABInertia::operator+(const SpatialRBInertia& m) const
{
  const unsigned X = 0, Y = 1, Z = 2;

  // do some preliminary calculations
  SpatialABInertia result;
  result.M = M;
  result.H = H;
  result.J = m.J + J;

  // update M with mass
  result.M(X,X) += m.m;
  result.M(Y,Y) += m.m;
  result.M(Z,Z) += m.m;

  // update H
  Matrix3 hx = Matrix3::skew_symmetric(m.h);
  result.H += hx;

  return result;
}

/// Adds two spatial matrices
SpatialABInertia SpatialABInertia::operator+(const SpatialABInertia& m) const
{
  SpatialABInertia result;
  result.M = this->M + m.M;
  result.H = this->H + m.H;
  result.J = this->J + m.J;
  return result;
}

/// Subtracts two spatial matrices
SpatialABInertia SpatialABInertia::operator-(const SpatialABInertia& m) const
{
  SpatialABInertia result;
  result.M = this->M - m.M;
  result.H = this->H - m.H;
  result.J = this->J - m.J;
  return result;
}

/// Adds m to this in place
SpatialABInertia& SpatialABInertia::operator+=(const SpatialABInertia& m)
{
  this->M += m.M;
  this->H += m.H;
  this->J += m.J;
  return *this;
}

/// Subtracts m from this in place
SpatialABInertia& SpatialABInertia::operator-=(const SpatialABInertia& m)
{
  this->M -= m.M;
  this->H -= m.H;
  this->J -= m.J;
  return *this;
}

/// Multiplies the inverse of this spatial AB inertia by a vector
SVector6 SpatialABInertia::inverse_mult(const SVector6& v) const
{
  Matrix3 nMinv = -Matrix3::inverse(M);
  Matrix3 HT = Matrix3::transpose(H);
  Matrix3 UR = Matrix3::inverse((H * nMinv * HT) + J);
  Matrix3 UL = UR * H * nMinv;
  Matrix3 LR = Matrix3::transpose(UL);
  Matrix3 LL = nMinv * ((HT * UL) - IDENTITY_3x3);

  // get the components of v
  Vector3 vtop = v.get_upper();
  Vector3 vbot = v.get_lower();

  // do the arithmetic
  return SVector6(UL*vtop + UR*vbot, LL*vtop + LR*vbot);
}

/// Multiplies a spatial matrix by a spatial matrix and returns the result in a spatial matrix
SMatrix6N& SpatialABInertia::mult(const SMatrix6N& m, SMatrix6N& result) const
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

  // compute the transpose of H
  Matrix3 HT = Matrix3::transpose(this->H);

  // carry out multiplication one column at a time
  for (unsigned i=0; i< NCOLS; i++)
  {
    SVector6 v = m.get_column(i);
    Vector3 vtop = v.get_upper();
    Vector3 vbot = v.get_lower();
    v = SVector6((this->M * vbot)+(HT * vtop), (this->J * vtop)+(H * vbot));
    result.set_column(i, v);
  } 

  return result;
}

/// Multiplies a 6x6 matrix by a Spatial AB inertia matrix
SpatialABInertia SpatialABInertia::mult(const MatrixN& m, const SpatialABInertia& I)
{
  assert(m.rows() == 6 && m.columns() == 6);

  // get the components of m
  Matrix3 UL, UR, LL, LR;
  m.get_sub_mat(0,3,0,3,UL);
  m.get_sub_mat(0,3,3,6,UR);
  m.get_sub_mat(3,6,0,3,LL);
  m.get_sub_mat(3,6,3,6,LR);
 
  // multiply by components of I
  SpatialABInertia result;
  result.M = (UL*I.M) + (UR*I.H);
  result.J = (LL*Matrix3::transpose(I.H)) + (LR*I.J);
  result.H = (LL*I.M) + (LR*I.H);
  return result;
}

/// Outputs this matrix to the stream
std::ostream& Moby::operator<<(std::ostream& out, const SpatialABInertia& m) 
{
  out << "spatial AB H:" << std::endl << m.H;
  out << "spatial AB M:" << std::endl << m.M;
  out << "spatial AB J:" << std::endl << m.J;
   
  return out;
}

