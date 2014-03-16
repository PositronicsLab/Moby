
void Simplex::update(const Origin3d& c)
{
  if (type == ePoint)
  {
    type = eSegment;
    _p2 = c;
  }
  else if (type == eSegment)
  {
    type = eTriangle;
    _p3 = c;
  }
}

/// Searches the simplex for the origin
bool Simplex::search_and_update(Origin3d& d)
{
  if (type == eSegment)
  {
    // setup the points
    const Origin3d& B = _p1;
    const Origin3d& A = _p2;
  }
} 

void GJK::do_gjk()
{
  // setup the initial support (a point)
  SVertex p = support(A, B);
  Simplex S = p;

  // setup the initial direction vector
  Origin3d d = -p.v;

  // GJK loop
  while (true)
  {
    // get the new point searching the current support in the direction of d
    Origin3d c = support(A, B, d);

    // compute dot product with point and direction; if dot product is negative,
    // there is no way for the origin to be in the simplex
    if (c.dot(d) < 0.0)
    {
      // no intersection, quit here
    }
    else
    {
      // add c to the simplex
      S.update(c);

      // search the simplex
      if (S.search_and_update(d))
  }
}
