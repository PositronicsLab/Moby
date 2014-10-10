/****************************************************************************
 * Copyright 2005 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _HASH_CLASSES_H
#define _HASH_CLASSES_H

#include <string.h>
#include <string>
#include <ext/hash_map>

namespace Moby {

/// Hashing function class for string objects
class StringHash
{
  public:
    size_t operator()(const std::string s) const { __gnu_cxx::hash<const char*> H; return H(s.c_str()); }
};

/// Equivalence function class for string objects (performs case-insensitive equivalence test)
class StringEqual
{
  public:
    bool operator()(const std::string s1, const std::string s2) const { return strcasecmp(s1.c_str(), s2.c_str()) == 0; }
};

} // end namespace Moby
#endif
