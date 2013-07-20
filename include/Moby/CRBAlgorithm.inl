template <class Mat>
Mat& transform_and_transpose_mult(const std::vector<Ravelin::SAxisd>& t, const std::vector<Ravelin::SMomentumd>& w, Mat& result)
{
  const unsigned M = t.size();
  const unsigned N = w.size();

  // resize the result
  result.resize(M, N);

  // get iterator to the result matrix
  Ravelin::dIterator data = result.begin();

  // iterate
  for (unsigned j=0; j< N; j++)
    for (unsigned i=0; i< M; i++)
      *data++ = t[i].dot(w[j]); 

  return result;
}

template <class Vec>
Vec& transform_and_transpose_mult(const std::vector<Ravelin::SAxisd>& t, const Ravelin::SMomentumd& w, Vec& result)
{
  const unsigned M = t.size();
  const unsigned SPATIAL_DIM = 6;

  // resize the result
  result.resize(M, SPATIAL_DIM, false);

  // get iterator to the result matrix
  Ravelin::dIterator data = result.begin();

  // iterate
  for (unsigned i=0; i< M; i++)
    *data++ = t[i].dot(w); 

  return result;
}

template <class Mat>
Mat& transform_and_transpose_mult(const std::vector<Ravelin::SAxisd>& t, const std::vector<Ravelin::SForced>& w, Mat& result)
{
  const unsigned M = t.size();
  const unsigned N = w.size();

  // resize the result
  result.resize(M, N);

  // get iterator to the result matrix
  Ravelin::dIterator data = result.begin();

  // iterate
  for (unsigned j=0; j< N; j++)
    for (unsigned i=0; i< M; i++)
      *data++ = t[i].dot(w[j]); 


  return result;
}

template <class Vec>
Vec& transform_and_transpose_mult(const std::vector<Ravelin::SAxisd>& t, const Ravelin::SForced& w, Vec& result)
{
  const unsigned M = t.size();
  const unsigned SPATIAL_DIM = 6;

  // resize the result
  result.resize(M, SPATIAL_DIM, false);

  // get iterator to the result matrix
  Ravelin::dIterator data = result.begin();

  // iterate
  for (unsigned i=0; i< M; i++)
    *data++ = t[i].dot(w); 

  return result;
}

