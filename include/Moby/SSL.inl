/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Outputs the SSL to the specified stream
inline std::ostream& operator<<(std::ostream& o, const SSL& s)
{
  o << "p1: " << std::endl << s.p1;
  o << "p2: " << std::endl << s.p2;
  o << "radius: " << s.radius << std::endl;
  return o;
}

