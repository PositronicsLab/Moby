#ifndef _MOBY_PERMUTE_
#define _MOBY_PERMUTE_

#include <iterator>

/** 
 * This function, which acts like an STL algorithm, selects elements with
 * given integer indices from a container.
 *
 * Example:
 * double x[10] = { 1, 2, 3, 4, 5, 6, 7, 8, 9, 10 };
 * int i[4] = {3, 1, 5, 7};
 * double results[4];
 * 
 * select(x, i, i+4, results);
 * // results now contains { 4, 2, 6, 8 }
 */ 

namespace Moby {

/**
 * \param mapping_begin beginning iterator to mapping giving source indices 
 * \param mapping_end ending iterator to mapping giving source indices
 * \param input_begin beginning iterator to source
 * \param output_begin beginning iterator to target 
 */
template <class ForwardIterator, class RandomAccessIterator1, class RandomAccessIterator2>
void permute(ForwardIterator mapping_begin, ForwardIterator mapping_end, RandomAccessIterator1 input_begin, RandomAccessIterator2 output_begin)
{
  // setup the target index
  unsigned target_index = 0;

  for (ForwardIterator i = mapping_begin; i != mapping_end; i++)
    output_begin[target_index++] = input_begin[*i]; 
}

} // end namespace

#endif

