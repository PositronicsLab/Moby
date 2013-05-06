template <class M>
M& transform_and_transpose_mult(const std::vector<Ravelin::Twistd>& t, const std::vector<Ravelin::Wrenchd>& w, M& result)
{

  return result;
}

template <class V>
V& transform_and_transpose_mult(const std::vector<Ravelin::Twistd>& t, const Ravelin::Wrenchd& w, V& result)
{

  return result;
}

