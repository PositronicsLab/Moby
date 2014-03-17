/// Finds the closest point to the origin
Origin3d Simplex::find_closest(unsigned& region) const
{
  // TODO: should we return what feature?

  if (_type == ePoint)
    return _v1;
  else if (_type == eSegment)
  {
    const Origin3d& a = _v2;
    const Origin3d& b = _v1;
    Origin3d ab = b - a;
    Origin3d aO = -a;
    if (ab.dot(aO) <= 0.0)
    {
      region = 1;
      return a;
    }
    else
    {
      // TODO: ?
      region = 2;
    }
  }
  else if (_type == eTriangle)
  {
  }
  else if (_type == eTetrahedron)
  {
  }
}

/// Adds a vertex to the simplex
void Simplex::add(const SVertex& v)
{
  if (_type == ePoint)
  {
    _type = eSeg;
    _v2 = v;
  }
  else if (_type == eSeg)
  {
    _type = eTriangle;
    _v3 = v;
  }
  else if (_type == eTriangle)
  {
    _type = eTetrahedron;
    _v4 = v;
  }
  else if (_type == eTetrahedron)
    _v5 = v;
}

/// Sees whether this is the best possible simplex (if not, returns a direction)
bool Simplex::best_simplex(Origin3d& o)
{
  if (
}

/// Gets the minimum norm point in the simplex
Origin3d Simplex::get_min_norm_point() const
{
  if (type == ePoint)
    return _p1.v;
  else if (type == eSegment)
  {
    const Origin3d& b = _p1.v;
    const Origin3d& a = _p2.v;
    Origin3d ab = b - a;
    Origin3d aO = -a;
    if (ab.dot(aO) > 0.0)
  }
}

/// Updates the simplex
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

/// Refines the simplex if possible 
void Simplex::refine()
{
}

void GJK::do_gjk(Point3d& closestA, Point3d& closestB)
{
  // setup the initial support (a point)
  SVertex p = initial_support(A, B);
  Simplex S = p;

  // GJK loop
  while (true)
  {
    // find the closest point in the simplex to the origin
    Origin3d p = S.find_closest();

    // look and see whether the origin is contained in the simplex
    double pnorm = p.norm();
    if (pnorm < NEAR_ZERO)
    {
      // A and B are intersecting

      // determine the interpenetration distance

      return -dist;
    }

    // get the new supporting points and determine the new vertex
    Point3d pA = A->get_supporting_point(-p);
    Point3d pB = B->get_supporting_point(p); 
    SVertex V(pA, pB);

    // look to see whether no intersection
    if (V.v.norm() > pnorm - NEAR_ZERO)
    {
      closestA = pA;
      closestB = pB;
      return pnorm;
    }
    else
    {
      // add the new vertex to the simplex
      S.add(V);

      // simplify the simplex
      S.refine(p);
    }
  }

  assert(false);
  return 0.0;
}
