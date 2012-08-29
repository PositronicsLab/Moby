#ifndef _MOBY_FAST_THREADABLE_H_
#define _MOBY_FAST_THREADABLE_H_

#ifdef _OPENMP
#include <omp.h>
#include <vector>
#endif

namespace Moby {

/// This class exists solely to make static declaration of dynamically allocated variables safe
template <class T>
class FastThreadable
{
  private:
    #ifdef _OPENMP
    std::vector<T> _x;
    #else
    T _x;
    #endif

  public:
    FastThreadable()
    {
      #ifdef _OPENMP
      _x.resize(omp_get_max_threads());
      #endif
    }

    T& operator()()
    {
      #ifdef _OPENMP
      return _x[omp_get_thread_num()];
      #else
      return _x;
      #endif
    }
};

} // end namespace

#endif

