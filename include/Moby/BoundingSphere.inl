template <class ForwardIterator>
BoundingSphere::BoundingSphere(ForwardIterator begin, ForwardIterator end)
{
  const Real INF = std::numeric_limits<Real>::max();
  const unsigned THREE_D = 3;

  // determine an AABB bounding the points
  Vector3 lo(INF, INF, INF), hi(-INF, -INF, -INF);
  for (ForwardIterator i = begin; i != end; i++)
    for (unsigned j=0; j< THREE_D; j++)
    {
      if ((*i)[j] < lo[j])
        lo[j] = (*i)[j];
      if ((*i)[j] > hi[j])
        hi[j] = (*i)[j];
    }

  // center of the circle will be at the center of the AABB
  center = (lo + hi) * (Real) 0.5;

  // set the radius to zero initially
  radius = (Real) 0.0;

  // now, determine the radius
  for (ForwardIterator i = begin; i != end; i++)
  {
    Real dist = (*i - center).norm();
    if (dist > radius)
      radius = dist;
  }
}

