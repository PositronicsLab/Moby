/// Gets all super bodies in constraints (from a container holding type Constraint*)
template <class ForwardIterator, class OutputIterator>
static OutputIterator get_super_bodies(ForwardIterator constraints_begin, ForwardIterator constraints_end, OutputIterator output_begin)
{
  // get all super bodies
  std::vector<boost::shared_ptr<Ravelin::DynamicBodyd> > supers;
  for (; constraints_begin != constraints_end; constraints_begin++)
    (*constraints_begin)->get_super_bodies(std::back_inserter(supers));

  // sort the container so we can use 'unique'
  std::sort(supers.begin(), supers.end());

  // copy the output
  supers.erase(std::unique(supers.begin(), supers.end()), supers.end());
  OutputIterator output_end = std::copy(supers.begin(), supers.end(), output_begin);
  return output_end;
} 
  
