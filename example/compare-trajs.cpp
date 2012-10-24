#include <assert.h>
#include <iostream>
#include <fstream>
#include <cstring>
#include <vector>
#include <list>
#include <cmath>
#include <limits>
#include <cstdlib>

static const unsigned BUF_SIZE = 100000;
static char buffer[BUF_SIZE];

unsigned count_lines(char* fname)
{
  std::ifstream in(fname);
  unsigned lines = 0;
  do
  {
    in.getline(buffer, BUF_SIZE);
  }
  while (!in.eof() && ++lines);

  return lines;
}

// parses a vector from a string
std::vector<double> parse(const std::string& s)
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
  std::vector<double> values(plist.size());
  unsigned idx = 0;
  for (std::list<std::string>::const_iterator i=plist.begin(); i != plist.end(); i++)
  {
    if (strcasecmp(i->c_str(), "inf") == 0)
      values[idx] = std::numeric_limits<double>::infinity();
    else if (strcasecmp(i->c_str(), "-inf") == 0)
      values[idx] = -std::numeric_limits<double>::infinity();
    else
      values[idx] = (double) atof(i->c_str());
    idx++;
  }
 
  return values;
}

// NOTE: compares using l-inf norm
double comp(const std::vector<double>& v1, const std::vector<double>& v2)
{
  assert(v1.size() == v2.size());
  double diff = 0.0;

  for (unsigned i=0; i< v1.size(); i++)
    diff = std::max(diff, std::fabs(v1[i] - v2[i]));

  return diff;
}

int main(int argc, char* argv[])
{
  if (argc < 3)
    return -1;

  // read the two files
  std::ifstream in1(argv[1]);
  std::ifstream in2(argv[2]);

  // verify we could open both of them
  if (in1.fail() || in2.fail())
  {
    std::cerr << "compare-trajs: unable to open one or both files" << std::endl;
    in1.close();
    in2.close();
    return -1;
  }

  // determine # of lines in the first file
  unsigned n = count_lines(argv[1]);
  if (n != count_lines(argv[2]))
  {
    std::cerr << "compare-trajs: unequal numbers of lines" << std::endl;
    in1.close();
    in2.close();
    return -1;
  }

  // read (n-1) lines (line n is timings)
  double max_diff = 0.0;
  for (unsigned i=1; i< n; i++)
  {
    in1.getline(buffer, BUF_SIZE);
    std::string str1(buffer);
    std::vector<double> v1 = parse(str1);    
    in2.getline(buffer, BUF_SIZE);
    std::string str2(buffer);
    std::vector<double> v2 = parse(str2);
    double diff = comp(v1, v2);
    max_diff = std::max(max_diff, diff);
  }

  // get the timing
  in1.getline(buffer, BUF_SIZE);
  std::string str1(buffer);
  std::vector<double> v1 = parse(str1);    
  in2.getline(buffer, BUF_SIZE);
  std::string str2(buffer);
  std::vector<double> v2 = parse(str2);
  assert(v1.size() == v2.size() && v1.size() == 1);

  std::cout << "maximum difference: " << max_diff << std::endl;
  std::cout << "reference timing: " << v1.front() << "  new timing: " << v2.front() << std::endl;

  return 0;
}

