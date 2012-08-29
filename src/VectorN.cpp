/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#include <iomanip>
#include <cstring>
#include <Moby/cblas.h>
#include <list>
#include <cmath>
#include <iostream>
#include <Moby/Constants.h>
#include <Moby/MatrixN.h>
#include <Moby/VectorN.h>

using namespace Moby;
using boost::shared_array;
using std::vector;

/// Default constructor - constructs an empty vector
VectorN::VectorN()
{
  _len = 0;
  _capacity = 0;
  _data = shared_array<Real>();
}

/// Constructs an uninitialized N-dimensional vector
VectorN::VectorN(unsigned N)
{
  _len = N;
  _capacity = N;
//  if (N > 0)
    _data = shared_array<Real>(new Real[N]);
}

/// Constructs a vector from a Vector2
VectorN::VectorN(const Vector2& v)
{
  const unsigned LEN = 2;
  _len = LEN;
  _capacity = LEN;
  _data = shared_array<Real>(new Real[LEN]);
  for (unsigned i=0; i< LEN; i++)
    _data[i] = v[i];
}

/// Constructs a vector from a Vector3
VectorN::VectorN(const Vector3& v)
{
  const unsigned LEN = 3;
  _len = LEN;
  _capacity = LEN;
  _data = shared_array<Real>(new Real[LEN]);
  for (unsigned i=0; i< LEN; i++)
    _data[i] = v[i];
}

/// Constructs a vector from a SVector6
VectorN::VectorN(const SVector6& v)
{
  const unsigned LEN = 6;
  _len = LEN;
  _capacity = LEN;
  _data = shared_array<Real>(new Real[LEN]);
  for (unsigned i=0; i< LEN; i++)
    _data[i] = v[i];
}

/// Constructs a N-dimensional vector from the given array
/**
 * \param array a N-dimensional (or larger) array
 */
VectorN::VectorN(unsigned N, const Real* array)
{
  _len = N;
  _capacity = N;
  _data = shared_array<Real>(new Real[N]);
  CBLAS::copy(N,array,1,_data.get(),1);
}

/// Constructs a N-dimensional vector from the given array
/**
 * \param array a N-dimensional (or larger) array
 * \note the VectorN object will modify the contents of the array until the
 *       VectorN object is destroyed
 */
VectorN::VectorN(unsigned N, shared_array<Real> array)
{
  _len = N;
  _capacity = N;
  _data = array;
}

/// Copy constructor
VectorN::VectorN(const VectorN& source)
{
  _len = source._len;
  _capacity = source._capacity;
  _data = shared_array<Real>(new Real[_len]);
  operator=(source);
}

/// Assignment constructor
/**
 * \note: this method was devised to copy temporaries without copying their
 *        array contents; this method changes its array pointer to point to
 *        that of source and changes source's point to NULL
 */
VectorN::VectorN(VectorN& source)
{
  _len = source._len;
  _capacity = source._capacity;
  _data = source._data;
  source._data = shared_array<Real>();
  source._len = 0;
  source._capacity = 0;
}

/// Calculates the relative error between two vectors
Real VectorN::calc_rel_err(const VectorN& vapprox, const VectorN& vtrue)
{
  if (vapprox.size() != vtrue.size())
    throw MissizeException();

  // init the relative error
  Real rel_err = (Real) 0.0;

  // get the size
  const unsigned SZ = vapprox.size();

  // compute the relative error
  for (unsigned i=0; i< SZ; i++)
    rel_err += std::pow((Real) 1.0 - vapprox[i]/vtrue[i], (Real) SZ);

  return std::pow(rel_err, (Real) 1.0/SZ);
}

/// Calculates the absolute error between two vectors
Real VectorN::calc_abs_err(const VectorN& vapprox, const VectorN& vtrue)
{
  if (vapprox.size() != vtrue.size())
    throw MissizeException();

  // init the absolute error
  Real abs_err = (Real) 0.0;

  // get the size
  const unsigned SZ = vapprox.size();

  // compute the absolute error
  for (unsigned i=0; i< SZ; i++)
    abs_err += std::pow(std::fabs(vapprox[i] - vtrue[i]), (Real) SZ);

  return std::pow(abs_err, (Real) 1.0/SZ);
}

/// Constructs a N-dimension vector from the list of double values
/**
 * \note There is no means in C++ to check the types of a list of variable
 * arguments.  If the variable arguments are not of type double, then
 * unexpected values will result. Constructing a vector using the
 * statement VectorN::construct_variable(3, 1.0, 1.0, 0) is incorrect because
 * the programmer has assumed that the integer 0 will be converted to a double
 * type; this is not the case.
 */
VectorN VectorN::construct_variable(unsigned N, ...)
{
  VectorN v(N);
  std::va_list lst;
  va_start(lst, N);
  for (unsigned i=0; i< N; i++)
    v[i] = (Real) va_arg(lst, double);
  va_end(lst);
  return v;
}

/// Gets the specified sub-vector
/**
 * \param start_idx the starting index (inclusive)
 * \param end_idx the ending index (exclusive)
 */
VectorN& VectorN::get_sub_vec(unsigned start_idx, unsigned end_idx, VectorN& v) const
{
  assert(start_idx <= end_idx);
  assert(end_idx <= size());

  // determine the subvector size
  unsigned N = end_idx - start_idx;

  // resize the subvector
  v.resize(N);

  // only copy if there is something to copy
  if (N > 0)
  {
    // use the BLAS routine for copying
    CBLAS::copy(N,_data.get()+start_idx,1, v._data.get(),1);
  }

  return v; 
}

/// Gets the specified sub-vector
/**
 * \param start_idx the starting index (inclusive)
 * \param end_idx the ending index (exclusive)
 */
VectorN VectorN::get_sub_vec(unsigned start_idx, unsigned end_idx) const
{
  // create the subvector
  VectorN subvec;

  return get_sub_vec(start_idx, end_idx, subvec); 
}

/// Sets the specified sub-vector
/**
 * \param start_idx the starting index (inclusive)
 * \param v a (end_idx - start_idx + 1)-dimensional vector
 */
VectorN& VectorN::set_sub_vec(unsigned start_idx, const VectorN& v)
{
  unsigned N = v.size();
  assert(start_idx + N <= size());

  // use the BLAS routine for copying
  if (N > 0)
    CBLAS::copy(N,v.data(),1, _data.get()+start_idx,1);

  return *this;
}

/// Computes the dot-product between two vectors
/**
 * Note: if the dimensions are unequal, an assertion will fail. 
 */
Real VectorN::dot(const VectorN& v1, const VectorN& v2)
{    
  assert(v1.size() == v2.size());
  if (v1.size() == 0)
    return (Real) 0.0;

  // use the BLAS routine for dot-product
  return CBLAS::dot(v1, v2);
}

/// Returns the normalized vector (l2 norm)
VectorN VectorN::normalize(const VectorN& v)
{
  // copy the vector
  VectorN v_copy = v;

  // normalize the copy and return it
  v_copy.normalize();

  return v_copy;
}

/// Computes the l-infinity norm of a vector
Real VectorN::norm_inf() const
{
  Real nrm = 0.0;

  for (unsigned i=0; i< _len; i++)
    nrm = std::max(nrm, std::fabs(_data[i]));

  return nrm;
}

/// Computes the l1-norm of a vector
Real VectorN::norm1(const VectorN& v)
{
  Real sum = 0.0;
  for (unsigned i=0; i< v.size(); i++)
    sum += std::fabs(v[i]);
  return sum;  
}

/// Computes the l2-norm of a vector
Real VectorN::norm(const VectorN& v)
{  
  if (v.size() == 0)
    return (Real) 0.0;

  // use the BLAS routine for l2-norm computation
  return CBLAS::nrm2(v);
}

/// Returns a N-dimensional one vector
VectorN VectorN::one(unsigned N)
{
  VectorN v(N);
  for (unsigned i=0; i< N; i++)
    v[i] = 1;
  return v;
}

/// Computes the outerproduct of two vectors
/**
 * If v1 is m-dimensional and v2 is n-dimensional, the outer product is m x n dimensional.
 */
MatrixN* VectorN::outer_prod(const VectorN& v1, const VectorN& v2, MatrixN* M)
{
  // get n and m
  unsigned m = v1.size();
  unsigned n = v2.size();

  // resize M
  M->resize(m,n);

  // set the matrix to zero 
  M->set_zero();

  // call BLAS outer-product routine
  if (m > 0 && n > 0)
  {
    MatrixN& MM = *M;
    CBLAS::ger(v1, v2, (Real) 1.0, MM); 
  }

  return M;
}

/// Sets all components of this vector to 0
VectorN& VectorN::set_zero()
{  
  unsigned dim = size();  
  for (unsigned i=0; i< dim; i++)
    _data[i] = (Real) 0.0;
  return *this;
}

/// Sets all components of this vector to 1
VectorN& VectorN::set_one()
{
  unsigned dim = size();  
  for (unsigned i=0; i< dim; i++)
    _data[i] = 1;
  return *this;
}

/// Returns a N-dimensional zero vector
VectorN VectorN::zero(unsigned N)
{
  VectorN v(N);
  v.set_zero();
  return v;
}

/// Resizes this vector, optionally preserving its existing elements
/**
 * \note memory is only reallocated if the vector grows in size; this keeps
 *       from having to do expensive reallocations if the vector shrinks
 */
VectorN& VectorN::resize(unsigned N, bool preserve)
{
  shared_array<Real> newdata;

  // if the array already is the proper size, exit
  if (_len == N)
    return *this;

  // see whether we can just downsize
  if (N < _capacity)
  {
    _len = N;
    return *this;
  }

  // create a new array
  newdata = shared_array<Real>(new Real[N]);

  // copy existing elements, if desired
  if (preserve)
    CBLAS::copy(_len, _data.get(),1,newdata.get(),1);

  // set the new data
  _data = newdata;
  _len = N;
  _capacity = N;
  return *this;
}

/// Copies another vector
VectorN& VectorN::operator=(const Vector3& source)
{
  // resize the vector if necessary
  if (_len != source.size())
    resize(source.size());

  // don't even worry about BLAS for copying
  _data[0] = source[0];
  _data[1] = source[1];
  _data[2] = source[2];

  return *this;
}

/// Copies another vector 
VectorN& VectorN::operator=(const VectorN& source)
{  
  // resize this vector if necessary
  if (_len != source.size())
    resize(source.size());

  // use the BLAS routine for copying
  if (_len > 0)
    CBLAS::copy(source.size(),source.data(),1,_data.get(),1);

  return *this;
}

/// Assigns this vector to the other
/**
 * \note: this method was devised to copy temporaries without copying their
 *        array contents; this method swaps array pointers with that of the
 *        source and changes the source's length to zero
 */
VectorN& VectorN::operator=(VectorN& source)
{
  _len = source._len;
  source._len = 0;
  std::swap(_capacity, source._capacity);
  source._capacity = 0;
  std::swap(_data, source._data);

  return *this;
}

/// Assigns this vector to a scalar
VectorN& VectorN::operator=(Real scalar)
{
  resize(1);
  _data[0] = scalar;
  return *this;
}

/// Adds this vector to another
VectorN VectorN::operator+(const VectorN& v) const
{
  VectorN this_copy = *this;
  this_copy += v;
  return this_copy;
}

/// Adds another vector to this one in place
VectorN& VectorN::operator+=(const VectorN& v)
{
  assert(v.size() == _len);

  // use the BLAS routine for addition
  if (_len > 0)
    CBLAS::axpy(v, *this, (Real) 1.0);

  return *this;
}

/// Returns the negation of this vector
VectorN VectorN::operator-() const
{
  VectorN v(_len);
  if (_len > 0)
    std::transform(begin(), end(), v.begin(), std::negate<Real>());
  return v;
}

/// Subtracts another vector from <b>this</b>
/**
 * Operation order: <b>this</b> - v
 */ 
VectorN VectorN::operator-(const VectorN& v) const
{  
  VectorN this_copy = *this;
  this_copy -= v;
  return this_copy;
}

/// Subtracts another vector from this one in place
VectorN& VectorN::operator-=(const VectorN& v)
{  
  assert(v.size() == _len);

  // use the BLAS routine for subtraction
  if (_len > 0)
    CBLAS::axpy(v, *this, (Real) -1.0);

  return *this;
}

/// Multiplies this vector by a scalar
VectorN VectorN::operator*(Real scalar) const
{
  VectorN this_copy = *this;
  this_copy *= scalar;
  return this_copy;
}

/// Multiplies this vector in place by a scalar
VectorN& VectorN::operator*=(Real scalar)
{
  if (_len > 0)
    CBLAS::scal(*this, scalar);  

  return *this;
}

/// Divides this vector by a scalar (same as vector * 1/scalar)
VectorN VectorN::operator/(Real scalar) const
{
  VectorN this_copy = *this;
  this_copy *= (1.0/scalar);
  return this_copy;
}

/// Writes a VectorN to the specified stream
std::ostream& Moby::operator<<(std::ostream& out, const VectorN& v)
{
  const unsigned OUTPUT_PRECISION = 8;

  if (v.size() == 0)
  {
    out << "(empty) ";
    return out;
  }
  out << "[";
  for (unsigned i=0; i< v.size()-1; i++)
    out << std::setprecision(OUTPUT_PRECISION) << v[i] << ", ";
  out << std::setprecision(OUTPUT_PRECISION) << v[v.size()-1] << "] ";
  return out;
}

// Reads a VectorN from the specified stream
std::istream& Moby::operator>>(std::istream& in, VectorN& v)
{
  unsigned n;
  in >> n;
  v.resize(n);

  for (unsigned i=0; i< n; i++)
    in >> v[i];

  return in;
}

/// Determines whether two numbers are relatively equal
bool VectorN::rel_equal(Real x, Real y)
{
  return (std::fabs(x-y) <= NEAR_ZERO * std::max(std::fabs(x), std::max(std::fabs(y), (Real) 1.0)));
}

/// Compares two vectors lexographically
bool VectorN::operator<(const VectorN& v) const
{
  // compare up to the shorter of the two
  unsigned shorter = std::min(size(), v.size());
  for (unsigned i=0; i< shorter; i++)
  {
    if (rel_equal(_data[i], v[i]))
      continue;
    return (_data[i] < v[i]);
  }
  
  // still here?  comparison was identical to this point; if this is shorter
  // than v, return true
  if (_len < v.size())
    return true;
  
  // vectors *may* be identical
  return false;
}

/// Compares two vectors
/**
 * \note this method exists solely for the convenience of iterators; generally
 *       epsilon_equals() should be used instead (because it specifies a comparison
 *       tolerance rather than using the default 1e-8)
 * \sa epsilon_equals()
 */
bool VectorN::operator==(const VectorN& v) const
{
  if (_len != v.size())
    return false;
  
  for (unsigned i=0; i< _len; i++)
    if (!rel_equal(_data[i], v[i]))
      return false;
    
  return true;
}

/// Compares two vectors using a tolerance epsilon
bool VectorN::epsilon_equals(const VectorN& v1, const VectorN& v2, Real epsilon)
{
  if (v1.size() != v2.size())
    return false;
  
  for (unsigned i=0; i< v1.size(); i++)
    if (std::fabs(v1[i] - v2[i]) > epsilon)
      return false;
    
  return true;  
}

/// Compares this with another vector using a tolerance epsilon
bool VectorN::epsilon_equals(const VectorN& v, Real epsilon) const
{
  return epsilon_equals(*this, v, epsilon);
}

/// Returns true if all components of this vector are not infinite and not NaN, and false otherwise
bool VectorN::is_finite() const
{
  for (unsigned i=0; i< _len; i++)
    if (std::isnan(_data[i]) || std::isinf(_data[i]))
      return false;

  return true;
}

/// Negates the matrix in place
VectorN& VectorN::negate()
{
  for (unsigned i=0; i< _len; i++)
    _data[i] = -_data[i];

  return *this;
}

/// Concatenates two vectors together and stores the result in a third
VectorN& VectorN::concat(const VectorN& v1, const VectorN& v2, VectorN& result)
{
  result.resize(v1.size() + v2.size());
  result.set_sub_vec(0, v1);
  result.set_sub_vec(v1.size(), v2);
  return result;
}

/// Concatenates two vectors together
VectorN VectorN::concat(const VectorN& v1, const VectorN& v2)
{
  VectorN result(v1.size() + v2.size());
  return concat(v1, v2, result);
}

/// Selects a subvector from this
VectorN VectorN::select(const vector<bool>& indices) const
{
  VectorN v;
  return select(indices, v);
}

/// Selects a subvector from this
VectorN& VectorN::select(const vector<bool>& indices, VectorN& result) const
{
  const unsigned n = size();
  assert(n == indices.size());

  // determine how many indices are selected and resize result vector
  unsigned nselect = 0;
  for (unsigned i=0; i< n; i++)
    if (indices[i])
      nselect++;
  result.resize(nselect);

  for (unsigned i=0, j=0; i< n; i++)
    if (indices[i])
      result[j++] = _data[i];

  return result;
}

/// Parses a string for a vector value
VectorN VectorN::parse(const std::string& s)
{
  std::list<std::string> plist;

  // make a copy of the string
  std::string copy = s;

  while (true)
  {
    // get the portion of the string before the delimiter
    size_t space_idx = copy.find_first_of(" \t");
    size_t comma_idx = copy.find(',');
        
    // case 1: no delimiter found
    if (space_idx == std::string::npos && comma_idx == std::string::npos)
    {
      plist.push_back(copy);
      break;
    }
        
    // case 2: delimiter found
    if ((space_idx != std::string::npos && space_idx < comma_idx) || comma_idx == std::string::npos)
    {
      plist.push_back(copy.substr(0,space_idx));
      copy = copy.substr(space_idx);
    }
    else
    {
      plist.push_back(copy.substr(0,comma_idx));
      copy = copy.substr(comma_idx);
    }
        
    // get the new string
    size_t firstidx = copy.find_first_not_of(" ,");
    if (firstidx == std::string::npos)
      break;
    else
      copy = copy.substr(firstidx);        
  }

  // convert the list to a Vector
  VectorN values(plist.size());
  unsigned idx = 0;
  for (std::list<std::string>::const_iterator i=plist.begin(); i != plist.end(); i++)
  {
    if (strcasecmp(i->c_str(), "inf") == 0)
      values[idx] = std::numeric_limits<Real>::infinity();
    else if (strcasecmp(i->c_str(), "-inf") == 0)
      values[idx] = -std::numeric_limits<Real>::infinity();
    else
      values[idx] = (Real) atof(i->c_str());
    idx++;
  }
  
  return values;  
}

