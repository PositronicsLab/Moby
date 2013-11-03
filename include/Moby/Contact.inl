template <class OutputIterator>
OutputIterator Contact::get_super_bodies(OutputIterator begin) const
{
  DynamicBodyPtr db1, db2;
  unsigned nb = get_super_bodies(db1, db2);
  switch (nb)
  {
    case 0: break;
    case 1: if (db1) *begin++ = db1; break;
    case 2: if (db1) *begin++ = db1; if (db2) *begin++ = db2; break;
    default:
      assert(false);
  }
  return begin;
} 


