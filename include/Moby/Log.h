/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _MOBY_LOG_H_
#define _MOBY_LOG_H_

#include <iostream>
#include <ctime>
#include <limits>
#include <sstream>
#include <fstream>
#include <boost/shared_ptr.hpp>

namespace Moby {

#ifdef NDEBUG
#define FILE_LOG(level) if (false) Log<OutputToFile>().get(level)
#define LOGGING(level) false
#else
#define FILE_LOG(level) if ((level & Log<OutputToFile>::reporting_level) > 0) Log<OutputToFile>().get(level)
#define LOGGING(level) ((level & Log<OutputToFile>::reporting_level) > 0)
#endif

struct OutputToFile
{
  static std::ofstream stream;

  static void output(const std::string& msg);
};

template <typename OutputPolicy>
class Log
{
  public:
    Log()
    {
      #ifdef THREADSAFE
      std::cerr << "Log() warning: this class is not thread-safe!" << std::endl;
      #endif
    }

    std::ostringstream& get(unsigned level = 0)
    {
      time_t rawtime;
      std::time(&rawtime);
      tm* ptm = std::gmtime(&rawtime);
      os << "- " << ptm->tm_hour << ":" << ptm->tm_min << ":" << ptm->tm_sec;
      os << " " << level << ": ";
      message_level = level;
      return os;
    }

    ~Log()
    {
      if ((message_level & reporting_level) > 0)
        OutputPolicy::output(os.str()); 
    }

    static unsigned reporting_level;

  private:
    std::ostringstream os;
    unsigned message_level;
}; // end class

// make the default warning level zero (no reporting)
template <class T>
unsigned Log<T>::reporting_level = 0;

} // end namespace

#endif

