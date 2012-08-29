#include <Moby/Log.h>

using namespace Moby;

std::ofstream OutputToFile::stream;

void OutputToFile::output(const std::string& msg)
{
  if (!stream.is_open())
  {
    std::ofstream stderr_stream("/dev/stderr", std::ofstream::app);
    stderr_stream << msg << std::flush;
    stderr_stream.close();
  }
  else
    stream << msg << std::flush;
}

