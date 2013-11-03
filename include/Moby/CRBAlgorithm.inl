template <class Mat>
Mat& transform_and_transpose_mult(const std::vector<Ravelin::SVelocityd>& t, const std::vector<Ravelin::SMomentumd>& w, Mat& result)
{
  const unsigned M = t.size();
  const unsigned N = w.size();

  // resize the result
  result.resize(M, N);

  // look for early exit
  if (M == 0 || N == 0)
  {
    result.set_zero();
    return result;
  }

  // get iterator to the result matrix
  Ravelin::ColumnIteratord data = result.column_iterator_begin();

  // do minimum number of necessary transformations
  if (M < N)
  {
    // iterate
    Ravelin::Pose3d::transform(w[0].pose, t, _tandt_tx);
    for (unsigned j=0; j< N; j++)
      for (unsigned i=0; i< M; i++)
        *data++ = _tandt_tx[i].dot(w[j]); 
  }
  else
  {
    // iterate
    Ravelin::Pose3d::transform(t[0].pose, w, _tandt_wx);
    for (unsigned j=0; j< N; j++)
      for (unsigned i=0; i< M; i++)
        *data++ = t[i].dot(_tandt_wx[j]); 
  }


  return result;
}

template <class Vec>
Vec& transform_and_transpose_mult(const std::vector<Ravelin::SVelocityd>& t, const Ravelin::SMomentumd& w, Vec& result)
{
  const unsigned M = t.size();
  const unsigned SPATIAL_DIM = 6;

  // resize the result
  result.resize(M, SPATIAL_DIM, false);

  // look for early exit
  if (M == 0)
    return result;

  // get iterator to the result matrix
  Ravelin::ColumnIteratord data = result.column_iterator_begin();

  // convert w to the same frame as t
  Ravelin::SMomentumd wx = Ravelin::Pose3d::transform(t[0].pose, w);

  // iterate
  for (unsigned i=0; i< M; i++)
    *data++ = t[i].dot(wx); 

  return result;
}

template <class Mat>
Mat& transform_and_transpose_mult(const std::vector<Ravelin::SVelocityd>& t, const std::vector<Ravelin::SForced>& w, Mat& result)
{
  const unsigned M = t.size();
  const unsigned N = w.size();

  // resize the result
  result.resize(M, N);

  // look for early exit
  if (M == 0 || N == 0)
  {
    result.set_zero();
    return result;
  }

  // get iterator to the result matrix
  Ravelin::ColumnIteratord data = result.column_iterator_begin();

  // do minimum number of necessary transformations
  if (M < N)
  {
    // iterate
    Ravelin::Pose3d::transform(w[0].pose, t, _tandt_tx);
    for (unsigned j=0; j< N; j++)
      for (unsigned i=0; i< M; i++)
        *data++ = _tandt_tx[i].dot(w[j]); 
  }
  else
  {
    // iterate
    Ravelin::Pose3d::transform(t[0].pose, w, _tandt_fx);
    for (unsigned j=0; j< N; j++)
      for (unsigned i=0; i< M; i++)
        *data++ = t[i].dot(_tandt_fx[j]); 
  }

  return result;
}

template <class Vec>
Vec& transform_and_transpose_mult(const std::vector<Ravelin::SVelocityd>& t, const Ravelin::SForced& w, Vec& result)
{
  const unsigned M = t.size();
  const unsigned SPATIAL_DIM = 6;

  // resize the result
  result.resize(M, 1, false);

  // see whether we can exit early
  if (M == 0)
    return result;

  // get iterator to the result matrix
  Ravelin::ColumnIteratord data = result.column_iterator_begin();

  // convert w to the same frame as t
  Ravelin::SForced wx = Ravelin::Pose3d::transform(t[0].pose, w);

  // iterate
  for (unsigned i=0; i< M; i++)
    *data++ = t[i].dot(wx); 

  return result;
}

