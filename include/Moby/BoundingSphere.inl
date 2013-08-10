template <class ForwardIterator>
BoundingSphere::BoundingSphere(ForwardIterator begin, ForwardIterator end)
{
  const double INF = std::numeric_limits<double>::max();
  const unsigned THREE_D = 3;

  // determine an AABB bounding the points
  Point3d lo(INF, INF, INF), hi(-INF, -INF, -INF);
  for (ForwardIterator i = begin; i != end; i++)
    for (unsigned j=0; j< THREE_D; j++)
    {
      if ((*i)[j] < lo[j])
        lo[j] = (*i)[j];
      if ((*i)[j] > hi[j])
        hi[j] = (*i)[j];
    }

  // center of the circle will be at the center of the AABB
  center = (lo + hi) * (double) 0.5;

  // set the radius to zero initially
  radius = (double) 0.0;

  // now, determine the radius
  for (ForwardIterator i = begin; i != end; i++)
  {
    double dist = (*i - center).norm();
    if (dist > radius)
      radius = dist;
  }
}

