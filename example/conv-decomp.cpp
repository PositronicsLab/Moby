#include <getopt.h>
#include <boost/foreach.hpp>
#include <Moby/Polyhedron.h>
#include <Moby/sorted_pair>
#include <Moby/CompGeom.h>
#include <string>
#include <queue>
#include <fstream>
#include <sstream>
#include <stack>

using namespace Ravelin;
using namespace Moby;
using boost::shared_ptr;

// flags and constants
bool DEBUG = false;
bool SPLIT_DEGENERATE = false;
bool RANDOMIZE = false;
double COPLANAR_TOL = 1e-8;
double COLINEAR_TOL = 1e-8;
double SPLIT_TOL = 1e-8;
double CONVEXITY_TOL = 1e-8;
double EDGE_ON_PLANE_TOL = 1e-8;
double MAX_TRI_AREA = std::numeric_limits<double>::max();

// list of vertices -- this is necessary to keep created vertices around
// until program exit
std::list<Point3d*> new_verts;

// a triangle with pointers to vertices
struct Tri
{
  Tri(const Point3d* a, const Point3d* b, const Point3d* c)
  {
    this->a = a;
    this->b = b;
    this->c = c;
  }

  const Point3d* a;
  const Point3d* b;
  const Point3d* c;
  Vector3d calc_normal() const { return Vector3d::cross(*b - *a, *c - *b); }
  double calc_offset(const Vector3d& normal) const { return normal.dot(*a); }
};

// an edge
struct Edge
{
  Edge(const Point3d* v1, const Point3d* v2) { v = make_sorted_pair(v1, v2); }

  // the two vertices of the edge
  sorted_pair<const Point3d*> v;

  // the two facets of the edge
  shared_ptr<const Tri> f1, f2;

  bool is_convex(double concavity_tol) const
  {
    // if not both faces set, edge is certainly not convex 
    assert(f1);
    if (!f2)
      return false;

    // get the plane going through f1
    Vector3d normal = f1->calc_normal();
    double d = f1->calc_offset(normal);

    // get the vertex of f2 not on this edge
    const Point3d* vert;
    if ((f2->a == v.first && f2->b == v.second) ||
        (f2->a == v.second && f2->b == v.first))
      vert = f2->c;
    else if ((f2->a == v.first && f2->c == v.second) ||
              (f2->a == v.second && f2->c == v.first))
      vert = f2->b;  
    else
    {
      assert((f2->b == v.first && f2->c == v.second) ||
              (f2->b == v.second && f2->c == v.first));
      vert = f2->a;
    }  

    // see whether the vertex is on the positive side of the plane
    return Vector3d::dot(normal, *vert) - d <= concavity_tol;
  }

  bool operator==(const Edge& e) const
  {
    return e.v == this->v;
  }

  bool operator<(const Edge& e) const
  {
    return this->v < e.v;
  }
};

// a plane
struct Plane
{
  double a, b, c, d;
  Plane() {}
  Plane(double nx, double ny, double nz, double dd)
  {
    const double NEAR_ZERO = std::sqrt(std::numeric_limits<double>::epsilon());
    a = nx; b = ny; c = nz;
    d = dd;
    double nrm = std::sqrt(a*a + b*b + c*c);
    assert(nrm > NEAR_ZERO);
    a /= nrm;
    b /= nrm;
    c /= nrm;
  }

  void operator=(const Plane& p)
  {
    a = p.a;
    b = p.b;
    c = p.c;
    d = p.d;
  }

  double plane_eq(const Vector3d& v) const
  {
    return v[0]*a + v[1]*b + v[2]*c - d;
  }

  bool on_plane(const Point3d& v) const
  {
    const double NEAR_ZERO = std::sqrt(std::numeric_limits<double>::epsilon());
    return std::fabs(v[0]*a + v[1]*b + v[2]*c - d) < NEAR_ZERO;
  }

  bool operator==(const Plane& p) const
  {
    const double NEAR_ZERO = std::sqrt(std::numeric_limits<double>::epsilon());
    double dot = a*p.a + b*p.b + c*p.c;
    if (std::fabs(dot - 1.0) < NEAR_ZERO && std::fabs(d - p.d) < NEAR_ZERO)
      return true;
    if (std::fabs(dot + 1.0) < NEAR_ZERO && std::fabs(d + p.d) < NEAR_ZERO)
      return true;
    return false;
  }

  bool operator<(const Plane& p) const
  {
    if (a < p.a)
      return true;
    else if (a == p.a)
    {
      if (b < p.b)
        return true;
      else if (b == p.b)
      {
        if (c < p.c)
          return true;
        else if (c == p.c)
          return d < p.d;
        else
          return false;
      }
      else
        return false;
    }
    else return false;
  }
};

typedef shared_ptr<Edge> EdgePtr;

class BSP
{
  public:
    BSP() { _negative = _positive = NULL; _split = false; _parent = NULL; }
    BSP(BSP* parent) { _negative = _positive = NULL; _split = false; _parent = parent; }
    ~BSP() { if (_split) { delete _negative; delete _positive; } }
    void split(const Plane& p);
    std::list<std::list<Point3d> > get_all_points() const;
    std::list<BSP*> get_leafs() const;
    std::list<Point3d> get_points() const;
    void get_subspaces(EdgePtr e, std::list<BSP*>& spaces);
    void to_vrml(std::ostream& out) const;
    void check() const;
    bool inside_or_on(const Point3d& v) const;
    const Plane& plane() const { return _plane; }
    std::list<BSP*> find_nodes(const Point3d& v) const;

    /// Determines whether this BSP has been split
    bool is_leaf() const { return !_split; }

    /// Gets the parent subspace
    BSP* parent() const { return _parent; }

    /// Gets the positive subspace
    BSP* positive() const { return _positive; }
  
    /// Gets the negative subspace
    BSP* negative() const { return _negative; }

    /// Gets all edges of this BSP
    const std::vector<EdgePtr>& get_edges() const { return _edges; }

    /// Adds edges to a leaf BSP node
    void add_edges(const std::vector<EdgePtr>& edges) { assert(!_split); _edges.insert(_edges.end(), edges.begin(), edges.end()); }

    /// The logging stream
    static std::ofstream _log;

  private:
    void check(const std::stack<std::pair<Plane, bool> >& hs) const;
    void to_vrml(std::ostream& out, std::list<std::pair<Plane, bool> > half_spaces) const;
    std::set<EdgePtr> find_edges(const std::list<Edge>& edges) const;

    BSP* _positive;      // the positive half-space
    BSP* _negative;      // the negative half-space
    BSP* _parent;        // the parent of the BSP (if any)
    bool _split;        // true if this BSP has been split
    std::vector<EdgePtr> _edges;
    Plane _plane;        // the splitting plane

    // edges that were created from splitting
    std::map<EdgePtr, EdgePtr> _new_edges;
};

/// because _log is static, it must be declared outside of the class as well
std::ofstream BSP::_log;

/// Finds all leaf BSP nodes that v belongs to
std::list<BSP*> BSP::find_nodes(const Point3d& v) const
{
  // if we're at this point, we can assume that this node is within all
  // halfspaces above
  if (is_leaf())
  {
    std::list<BSP*> l;
    l.push_back((BSP*) this);
    return l;
  }
  else
  {
    // check where v is with regard to the splitting plane
    double val = _plane.plane_eq(v);
    if (val > SPLIT_TOL)
      return _positive->find_nodes(v);
    else if (val < -SPLIT_TOL)
      return _negative->find_nodes(v);
    else
    {
      std::list<BSP*> l1 = _positive->find_nodes(v);
      std::list<BSP*> l2 = _negative->find_nodes(v);
      l1.insert(l1.end(), l2.begin(), l2.end());
      return l1;
    }
  }
}

/// Determines whether a vertex is inside or on this BSP node
bool BSP::inside_or_on(const Point3d& point) const
{
  // get all splitting planes
  std::list<std::pair<Plane, bool> > planes;
  BSP* node = (BSP*) this;
  while (true)
  {
    BSP* new_node = node->parent();
    if (!new_node)
      break;
    bool pos = (new_node->positive() == node);
    planes.push_back(std::make_pair(new_node->plane(), pos));
    node = new_node;
  }

  // check signed distance from all splitting planes
  for (std::list<std::pair<Plane, bool> >::const_iterator j = planes.begin(); j != planes.end(); j++)
  {
    double sdist = j->first.plane_eq(point);
    if (sdist < 0.0)
    {
      if (j->second && -sdist > SPLIT_TOL)
        return false;
    }
    else if (!j->second && sdist > SPLIT_TOL)
      return false;  
  }

  // if we made it here, must be inside or on 
  return true; 
}

/// Gets all leaf nodes
std::list<BSP*> BSP::get_leafs() const
{
  // if this is a leaf node, create a list containing only this node & return it
  if (is_leaf())
  {
    std::list<BSP*> l;
    l.push_back((BSP*) this);
    return l;
  }
  else
  {
    std::list<BSP*> l1 = _positive->get_leafs();
    std::list<BSP*> l2 = _negative->get_leafs();
    l1.insert(l1.end(), l2.begin(), l2.end());
    return l1;
  }
}

/// Writes this BSP in VRML format to an output stream
void BSP::to_vrml(std::ostream& out) const
{
  const unsigned X = 0, Y = 1, Z = 2, THREE_D = 3;
  const double INF = std::numeric_limits<double>::max();

  // determine the bounding box on all points
  std::list<std::list<Point3d> > points = get_all_points();
  Point3d bb_lo(INF, INF, INF), bb_hi(-INF, -INF, -INF);
  for (std::list<std::list<Point3d> >::const_iterator i = points.begin(); i != points.end(); i++)
    for (std::list<Point3d>::const_iterator j = i->begin(); j != i->end(); j++)
      for (unsigned k=0; k< THREE_D; k++)
      {
        bb_lo[k] = std::min(bb_lo[k], (*j)[k]);
        bb_hi[k] = std::max(bb_hi[k], (*j)[k]);
      }

  // setup half-spaces accordingly
  std::list<std::pair<Plane, bool> > hs;
  hs.push_back(std::make_pair(Plane(1,0,0,bb_hi[X]), false));
  hs.push_back(std::make_pair(Plane(0,1,0,bb_hi[Y]), false));
  hs.push_back(std::make_pair(Plane(0,0,1,bb_hi[Z]), false));
  hs.push_back(std::make_pair(Plane(1,0,0,bb_lo[X]), true));
  hs.push_back(std::make_pair(Plane(0,1,0,bb_lo[Y]), true));
  hs.push_back(std::make_pair(Plane(0,0,1,bb_lo[Z]), true));

  // call private method
  to_vrml(out, hs);
}

/// Writes this BSP in VRML format to an output stream
void BSP::to_vrml(std::ostream& out, std::list<std::pair<Plane, bool> > half_spaces) const
{
  double r_color, g_color, b_color;

  // if this is a leaf, do nothing (leaves should be written by their parents)
  if (is_leaf())
    return;

  // setup the half-space for CompGeom::find_hs_interior_point()
  std::list<std::pair<Point3d, double> > hs;
  for (std::list<std::pair<Plane, bool> >::const_iterator i = half_spaces.begin(); i != half_spaces.end(); i++)
    if (i->second)
      hs.push_back(std::make_pair(Point3d(-i->first.a, -i->first.b, -i->first.c), -i->first.d));
    else
      hs.push_back(std::make_pair(Point3d(i->first.a, i->first.b, i->first.c), i->first.d));

  // add the positive half-space of this child to the set of half-spaces
  hs.push_back(std::make_pair(Point3d(-_plane.a, -_plane.b, -_plane.c), -_plane.d));

  // find the interior point to the halfspaces
  Point3d ip;
// NOTE: we've disabled the following line b/c of removal of linear programming from Moby (hence this utility does not work!) 
//  double dist = CompGeom::find_hs_interior_point(hs.begin(), hs.end(), ip);
double dist = -1.0;
  if (dist < 0.0)
  {
    std::cerr << "could not find point interior to positive half-spaces!" << std::endl;
    std::cerr << "halfspaces:" << std::endl;
    for (std::list<std::pair<Point3d, double> >::const_iterator i=hs.begin(); i != hs.end(); i++)
      std::cerr << "  normal: " << i->first << " offset: " << i->second << std::endl;

    return;
  }

  // compute the halfspace intersection
// NOTE; this function is disabled (and hence this utility will not work), unless we re-integrate linear programming with Moby
//  PolyhedronPtr hs_isect = CompGeom::calc_hs_intersection(hs.begin(), hs.end(), ip);
  PolyhedronPtr hs_isect;
  if (hs_isect && _positive->is_leaf())
  {
    // write only leaves to VRML as a wireframe
    r_color = (double) rand() / RAND_MAX;
    g_color = (double) rand() / RAND_MAX;
    b_color = (double) rand() / RAND_MAX;
    out << "# halfspace: " << _positive << std::endl;
    Polyhedron::to_vrml(out, *hs_isect, Point3d(r_color, g_color, b_color), true);
    out << std::endl;
  }

  // now, write the negative half-space
  hs.pop_back();
  hs.push_back(std::make_pair(Point3d(_plane.a, _plane.b, _plane.c), _plane.d));

  // find the interior point to the halfspaces
// NOTE: we've disabled the following line b/c of removal of linear programming from Moby (hence this utility does not work!) 
//  dist = CompGeom::find_hs_interior_point(hs.begin(), hs.end(), ip);
  dist = -1.0;
  if (dist < 0.0)
  {
    std::cerr << "could not find point interior to negative half-spaces!" << std::endl;
    std::cerr << "halfspaces:" << std::endl;
    for (std::list<std::pair<Point3d, double> >::const_iterator i=hs.begin(); i != hs.end(); i++)
      std::cerr << "  normal: " << i->first << " offset: " << i->second << std::endl;

    return;
  }

  // compute the halfspace intersection
// NOTE; this function is disabled (and hence this utility will not work), unless we re-integrate linear programming with Moby
//  hs_isect = CompGeom::calc_hs_intersection(hs.begin(), hs.end(), ip);
  if (hs_isect && _negative->is_leaf())
  {
    // write only leaves to VRML as a wireframe
    r_color = (double) rand() / RAND_MAX;
    g_color = (double) rand() / RAND_MAX;
    b_color = (double) rand() / RAND_MAX;
    out << "# halfspace: " << _negative << std::endl;
    Polyhedron::to_vrml(out, *hs_isect, Point3d(r_color, g_color, b_color), true);
    out << std::endl;
  }

  // continue recursively if positive child is internal node
  if (!_positive->is_leaf())
  {
    // add this plane to the half-spaces
    half_spaces.push_back(std::make_pair(_plane, true));

    // call recursively
    _positive->to_vrml(out, half_spaces);

    // remove the plane from the half spaces
    half_spaces.pop_back();
  }

  // continue recursively if negative child is internal node
  if (!_negative->is_leaf())
  {
    // add this plane to the half-spaces
    half_spaces.push_back(std::make_pair(_plane, false));

    // call recursively
    _negative->to_vrml(out, half_spaces);
  }
}

/// Finds all of the edges equal to the elements in the list of edges within this BSP
std::set<EdgePtr> BSP::find_edges(const std::list<Edge>& edges) const
{
  if (is_leaf())
  {
    // create the set
    std::set<EdgePtr> found_edges;

    // look through the vector of edges first
    for (unsigned i=0; i< _edges.size(); i++)
      for (std::list<Edge>::const_iterator j = edges.begin(); j != edges.end(); j++)
        if (*_edges[i] == *j)
          found_edges.insert(_edges[i]);

    // now, look through the new edges to see whether this edge was split into another
    for (std::map<EdgePtr, EdgePtr>::const_iterator i = _new_edges.begin(); i != _new_edges.end(); i++)
      for (std::list<Edge>::const_iterator j = edges.begin(); j != edges.end(); j++)
        if (*i->first == *j)
          found_edges.insert(i->first);

    return found_edges;      
  }
  else
  {
    std::set<EdgePtr> e1 = _positive->find_edges(edges);
    std::set<EdgePtr> e2 = _negative->find_edges(edges);
    e1.insert(e2.begin(), e2.end());

    // now, look through the new edges (even though this is not a leaf) to see whether 
    // this edge was split into another
    for (std::map<EdgePtr, EdgePtr>::const_iterator i = _new_edges.begin(); i != _new_edges.end(); i++)
      for (std::list<Edge>::const_iterator j = edges.begin(); j != edges.end(); j++)
        if (*i->first == *j)
          e1.insert(i->first);

    return e1;
  }
}

/// Gets all points from this (leaf) BSP only
std::list<Point3d> BSP::get_points() const
{
  assert(is_leaf());

  // determine all unique vertices
  std::set<const Point3d*> unique_points;
  for (unsigned i=0; i< _edges.size(); i++)
  {
    unique_points.insert(_edges[i]->v.first);
    unique_points.insert(_edges[i]->v.second);
  }

  // create a list and add the vertices to it
  std::list<Point3d> points;
  for (std::set<const Point3d*>::const_iterator i = unique_points.begin(); i != unique_points.end(); i++)
    points.push_back(**i);
  return points;
}

/// Gets all points from this BSP and all children
std::list<std::list<Point3d> > BSP::get_all_points() const
{
  // if this is a leaf, just add all unique points
  if (is_leaf())
  {
    // create the list
    std::list<std::list<Point3d> > points;
    points.push_back(get_points());
    return points;
  }
  else
  {
    // get the points from both children
    std::list<std::list<Point3d> > p1 = _positive->get_all_points();
    std::list<std::list<Point3d> > p2 = _negative->get_all_points();
    p1.insert(p1.end(), p2.begin(), p2.end());
    return p1;
  }
}

/// Gets all subspaces to which the given edge belongs
void BSP::get_subspaces(EdgePtr e, std::list<BSP*>& spaces)
{
  // output logging info
  _log << "BSP::get_subspaces()- checking space " << this << " for edge " << e << std::endl;

  // handle leafs and internal nodes differently...
  if (is_leaf())
  {
    _log << " -- leaf node encountered.. checking: ";

    // check for this edge in set of edges
    for (unsigned i=0; i< _edges.size(); i++)
      if (_edges[i] == e)
      {
        spaces.push_back(this);
        _log << "  !!! edge found !!!" << std::endl;
        return;
      }
  
      _log << "  edge not found" << std::endl;
  }
  else
  {
    _log << "  -- non-leaf node encountered; checking to see whether edge has been split" << std::endl;

    // see whether this edge has been split
    if (_positive->_new_edges.find(e) != _positive->_new_edges.end())
    {
      _log << "   -- edge has been split; search for edges " << _positive->_new_edges[e] << ", " << _negative->_new_edges[e] << " in child nodes" << std::endl;

      // it has been split; search for split edge in subspaces
      assert(_negative->_new_edges.find(e) != _negative->_new_edges.end());
      _positive->get_subspaces(_positive->_new_edges[e], spaces);
      _negative->get_subspaces(_negative->_new_edges[e], spaces);
    }
    else
    {
      _log << "   -- edge has not been split; searching for " << e << " in child nodes" << std::endl; 

      // look for this exact edge in the subspaces
      _positive->get_subspaces(e, spaces);
      _negative->get_subspaces(e, spaces);
    }
  }

  _log << "BSP::get_subspaces() complete" << std::endl;
}

/// Determines whether three points are collinear to some given tolerance
bool collinear(const Point3d& a, const Point3d& b, const Point3d& c, double tol)
{
  const unsigned X = 0, Y = 1, Z = 2;

  return (std::fabs((c[Z]-a[Z])*(b[Y]-a[Y]) - (b[Z]-a[Z])*(c[Y]-a[Y])) < tol &&
      std::fabs((b[Z]-a[Z])*(c[X]-a[X]) - (b[X]-a[X])*(c[Z]-a[Z])) < tol &&
      std::fabs((b[X]-a[X])*(c[Y]-a[Y]) - (b[Y]-a[Y])*(c[X]-a[X])) < tol);
}

/// Gets the distance from all vertices to the splitting halfspaces and sends output to log
void BSP::check() const
{
  check(std::stack<std::pair<Plane, bool> >());  
}

void BSP::check(const std::stack<std::pair<Plane, bool> >& hs) const
{
  const double NEAR_ZERO = std::sqrt(std::numeric_limits<double>::epsilon());

  // if this is an internal node, add splitting plane and call recursively
  // on child nodes
  if (!is_leaf())
  {
    std::stack<std::pair<Plane, bool> > hs_copy_pos = hs;
    std::stack<std::pair<Plane, bool> > hs_copy_neg = hs;
    hs_copy_pos.push(std::make_pair(_plane, true));
    hs_copy_neg.push(std::make_pair(_plane, false));
    _positive->check(hs_copy_pos);
    _negative->check(hs_copy_neg);
  }
  else
  {
    // this is an internal node; copy the stack
    std::stack<std::pair<Plane, bool> > hs_copy = hs;
    std::list<std::pair<Plane, bool> > hs_list;
    while (!hs_copy.empty())
    {
      hs_list.push_back(hs_copy.top());
      hs_copy.pop();
    }

    // determine set of vertices
    std::set<const Point3d*> verts;
    for (unsigned i=0; i< _edges.size(); i++)
    {
      verts.insert(_edges[i]->v.first);
      verts.insert(_edges[i]->v.second);
    }

    // get the signed distance of each vertex from each plane
    for (std::set<const Point3d*>::const_iterator i = verts.begin(); i != verts.end(); i++)
    {
      double max_outside_distance = 0.0;
      for (std::list<std::pair<Plane, bool> >::const_iterator j = hs_list.begin(); j != hs_list.end(); j++)
      {
        double sdist = j->first.plane_eq(**i);
        if (sdist < 0.0)
        {
          if (j->second)
            max_outside_distance = std::max(max_outside_distance, -sdist);
        }
        else if (!j->second)
          max_outside_distance = std::max(max_outside_distance, sdist);
      }
      if (max_outside_distance > NEAR_ZERO)
        _log << "  maximum outside distance from " << **i << " (" << *i << ") to halfspace " << this << " = " << max_outside_distance << std::endl;
    }
  }
}

/// Splits this BSP using the given splitting plane
void BSP::split(const Plane& plane)
{
  const double NEAR_ZERO = std::sqrt(std::numeric_limits<double>::epsilon());

  // guarantee that the BSP is not already split
  assert(!_split);

  // setup positive and negative edges
  std::vector<EdgePtr> pedges, nedges;

  // setup maps of old edges to new edges
  std::map<EdgePtr, EdgePtr> npedges, nnedges;

  // output logging info
  _log << "splitting space " << this << " using splitting plane: " << plane.a << "x + " << plane.b << "y + " << plane.c << "z = " << plane.d << std::endl;

  // get the plane normal
  Vector3d pnormal(plane.a, plane.b, plane.c);

  // determine where all edges lie with respect to the splitting plane
  for (unsigned i=0; i< _edges.size(); i++)
  {
    double val1 = plane.plane_eq(*_edges[i]->v.first);
    double val2 = plane.plane_eq(*_edges[i]->v.second);
    Vector3d nedge = (*_edges[i]->v.second - *_edges[i]->v.first);
    if (nedge.norm() > NEAR_ZERO)
      nedge.normalize();

    // output logging info
    _log << "examining edge: " << _edges[i] << std::endl;
    _log << "  -- vertices: " << _edges[i]->v.first << ", " << _edges[i]->v.second << std::endl;
    _log << "     (" << *_edges[i]->v.first << "; " << *_edges[i]->v.second << ") " << std::endl;
    _log << "  -- vertex one value: " << val1 << std::endl;
    _log << "  -- vertex two value: " << val2 << std::endl;

    if (val1 > SPLIT_TOL && val2 > SPLIT_TOL)
    {
      // output logging info
      _log << " + both vertices in positive space" << std::endl;

      // both vertices are on positive side of the plane
      pedges.push_back(_edges[i]);
    }
    else if (val1 < -SPLIT_TOL && val2 < -SPLIT_TOL)
    {
      // output logging info
      _log << " + both vertices in negative space" << std::endl;

      // both vertices are on negative side of the plane
      nedges.push_back(_edges[i]);
    }
    // old test: would not properly classify edges (nearly) on plane
    else if (std::fabs(pnormal.dot(nedge)) < EDGE_ON_PLANE_TOL)
//    else if (std::fabs(val1) <= SPLIT_TOL && std::fabs(val2) <= SPLIT_TOL)
    {
      // output logging info
      _log << " + edge is coplanar; adding to both spaces" << std::endl;

      // edges directly on the plane should be added to both sets
      pedges.push_back(_edges[i]);
      nedges.push_back(_edges[i]);
    }
    else
    {
      // edge lies across the plane; split the edge by the plane
      double t = CompGeom::intersect_seg_plane(Point3d(plane.a, plane.b, plane.c), plane.d, LineSeg3(*_edges[i]->v.first, *_edges[i]->v.second));

      // form the vertex on the plane
      Point3d* vert = new Point3d(*_edges[i]->v.first + (*_edges[i]->v.second - *_edges[i]->v.first) * t);
      assert(std::fabs(plane.plane_eq(*vert)) < NEAR_ZERO);

      // save the vertex
      new_verts.push_back(vert);

      // create two new edges from the original edge
      EdgePtr e1(new Edge(_edges[i]->v.first, vert));
      EdgePtr e2(new Edge(_edges[i]->v.second, vert));

      // determine into what halfspaces the new edges will go; b/c one of
      // either val1 or val2
      // might be near zero, we test the value with the bigger magnitude (both
      // cannot be near zero from the first 'if' test above) and then use that
      // to determine into which half-spaces the edges will be placed
      if (std::fabs(val1) > std::fabs(val2))
      {
        if (val1 > SPLIT_TOL)
        {
          pedges.push_back(e1);
          nedges.push_back(e2);
          npedges[_edges[i]] = e1;
          nnedges[_edges[i]] = e2;
        }
        else
        {
          pedges.push_back(e2);
          nedges.push_back(e1);
          npedges[_edges[i]] = e2;
          nnedges[_edges[i]] = e1;
        }
      }
      else
      {
        if (val2 > SPLIT_TOL)
        {
          pedges.push_back(e2);
          nedges.push_back(e1);
          npedges[_edges[i]] = e2;
          nnedges[_edges[i]] = e1;
        }
        else
        {
          pedges.push_back(e1);
          nedges.push_back(e2);
          npedges[_edges[i]] = e1;
          nnedges[_edges[i]] = e2;
        }
      }

      // output logging info
      _log << " + plane splits vertices" << std::endl;
      _log << "   - newly created vertex: " << vert << "   " << *vert << std::endl;
      _log << "   - newly created edges: " << e1 << ", " << e2 << std::endl;
    }

  }

  // get all vertices from positive and negative sides
  std::set<const Point3d*> pverts, nverts;
  for (unsigned i=0; i< pedges.size(); i++)
  {
    pverts.insert(pedges[i]->v.first);
    pverts.insert(pedges[i]->v.second);
  }
  for (unsigned i=0; i< nedges.size(); i++)
  {
    nverts.insert(nedges[i]->v.first);
    nverts.insert(nedges[i]->v.second);
  }

  // output logging info
  _log << "# of positive vertices: " << pverts.size() << std::endl;
  _log << "# of negative vertices: " << nverts.size() << std::endl;

  // make sure that splitting the BSP is necessary
  if (pverts.size() < 4 || nverts.size() < 4)
  {
    if (!SPLIT_DEGENERATE || pverts.empty() || nverts.empty())
    {
      _log << "too few vertices in one space; not splitting" << std::endl << std::endl << std::endl;  
      return;
    }
  }

  // determine three noncolinear points on the positive side
  if (!SPLIT_DEGENERATE)
  {
    std::set<const Point3d*>::const_iterator k = pverts.begin();
    Point3d PA = **(k++);
    Point3d PB(0,0,0), PC(0,0,0);
    while (true)
    {
      std::set<const Point3d*>::const_iterator j = k;
      j++;
      if (j == pverts.end())
      {
        k = pverts.end();
        break;
      }
      if (!collinear(PA, **k, **j, COLINEAR_TOL))
      {
        PB = **k;
        PC = **j;
        break;
      }
      k = j;
    }

    // check for all points are colinear (no split possible)
    if (k == pverts.end())
    {
      _log << "all points colinear in positive space; not splitting" << std::endl << std::endl << std::endl;  
      return;
    }

    // determine the normal of the plane through the positive points
    Vector3d Pn = Vector3d::cross(PB-PA, PC-PA);
    Pn.normalize();

    // determine three noncolinear points on the negative side
    k = nverts.begin();
    Point3d NA = **(k++);
    Point3d NB(0,0,0), NC(0,0,0);
    while (true)
    {
      std::set<const Point3d*>::const_iterator j = k;
      j++;
      if (j == nverts.end())
      {
        k = nverts.end();
        break;
      }
      if (!collinear(NA, **k, **j, COLINEAR_TOL))
      {
        NB = **k;
        NC = **j;
        break;
      }
      k = j;
    }

    // check for all points are colinear (no split possible)
    if (k == nverts.end())
    {
      _log << "all points colinear in negative space; not splitting" << std::endl << std::endl << std::endl;  
      return;
    }

    // determine the normal of the plane through the negative points
    Vector3d Nn = Vector3d::cross(NB-NA, NC-NA);
    Nn.normalize();

    // check for all points coplanar on the positive side
    double Pd = Vector3d::dot(Pn, PA);
    for (k = ++(++pverts.begin()); k != pverts.end(); k++)
      if (std::fabs(Vector3d::dot(Pn, **k) - Pd) > COPLANAR_TOL)
        break;

    // all points coplanar indicates no split possible
    if (k == pverts.end())
    {
      _log << "vertices in positive space all coplanar; not splitting" << std::endl << std::endl << std::endl;  
      return;  
    }

    // check for all points coplanar on the negative side
    double Nd = Vector3d::dot(Nn, NA);
    for (k = ++(++nverts.begin()); k != nverts.end(); k++)
      if (std::fabs(Vector3d::dot(Nn, **k) - Nd) > COPLANAR_TOL)
        break;

    // all points coplanar indicates no split possible
    if (k == nverts.end())
    {
      _log << "vertices in negative space all coplanar; not splitting" << std::endl << std::endl << std::endl;  
      return;  
    }
  }

  // indicate that the BSP has been split
  _split = true;

  // set the splitting plane
  _plane = plane;

  // create the subspaces
  _positive = new BSP(this);
  _negative = new BSP(this);

  // add the edges to the proper subspaces
  _positive->add_edges(pedges);
  _negative->add_edges(nedges);

  // set the maps
  _positive->_new_edges = npedges;
  _negative->_new_edges = nnedges;

  // clear the vector of edges
  _edges.clear();

   // output debugging info
  _log << "created new subspaces: " << _positive << ", " << _negative << std::endl;
  _log << "  vertices in subspace " << _positive << ": " << std::endl;
  std::list<Point3d> points = _positive->get_points();
  BOOST_FOREACH(const Point3d& p, points)
    _log << "   " << p << std::endl;  
  _log << "  edges in subspace " << _positive << ": " << std::endl;
  for (unsigned i=0; i< pedges.size(); i++)
    _log << "   " << pedges[i] << std::endl;
  for (std::map<EdgePtr, EdgePtr>::const_iterator i = npedges.begin(); i != npedges.end(); i++)
    _log << "   " << i->second << " [from " << i->first << "]" << std::endl;
  _log << "  vertices in subspace " << _negative << ": " << std::endl;
  points = _negative->get_points();
  BOOST_FOREACH(const Point3d& p, points)
    _log << "   " << p << std::endl;  
  _log << "  edges in subspace " << _negative << ": " << std::endl;
  for (unsigned i=0; i< nedges.size(); i++)
    _log << "   " << nedges[i] << std::endl;
  for (std::map<EdgePtr, EdgePtr>::const_iterator i = nnedges.begin(); i != nnedges.end(); i++)
    _log << "   " << i->second << " [from " << i->first << "]" << std::endl;
  _log << std::endl << std::endl;
}

// writes a collection of edges to VRML
template <class ForwardIterator>
void edges_to_vrml(std::ostream& out, ForwardIterator begin, ForwardIterator end)
{
  BSP::_log << "-- writing edges to VRML" << std::endl; 
  out << "IndexedLineSet { " << std::endl;
  out << "  coord Coordinate {" << std::endl;
  out << "    point [ " << std::endl;
  std::map<const Point3d*, unsigned> ncedgesV;
  for (ForwardIterator i = begin; i != end; i++)
  {
    BSP::_log << "  writing edge: " << *i << std::endl;
    BSP::_log << "     vertex 1: " << (*i)->v.first << "  " << *(*i)->v.first << std::endl;
    BSP::_log << "     vertex 2: " << (*i)->v.second << "  " << *(*i)->v.second << std::endl;
    out << (*(*i)->v.first)[0] << " " << (*(*i)->v.first)[1] << " " << (*(*i)->v.first)[2] << ", ";
    out << (*(*i)->v.second)[0] << " " << (*(*i)->v.second)[1] << " " << (*(*i)->v.second)[2] << ", ";
  }
  out << "] }" << std::endl << std::endl; 
  out << "  coordIndex [ ";
  unsigned idx = 0;
  for (ForwardIterator i = begin; i != end; i++)
  {
      out << idx++ << " ";
      out << idx++ << " -1, ";
  }
  out << "] }" << std::endl;
}

/// Prints out the syntax for the decomposer
void print_syntax()
{
  std::cerr << std::endl << "syntax: convex-decomp [options] <input>" << std::endl;
  std::cerr << std::endl;
  std::cerr << "convex-decomp takes the description of a 3D geometry from the input file" << std::endl;
  std::cerr << "  (a Wavefront OBJ file) and decomposes the geometry into a number of convex pieces;" << std::endl;
  std::cerr << "  each convex piece is written as a Wavefront OBJ file." << std::endl << std::endl;
  std::cerr << "options:" << std::endl;
  std::cerr << "  -d, --debug              sets decomposer to debug mode" << std::endl;
  std::cerr << "  -v, --cvxtol VAL         sets convexity tolerance to VAL (default = 1e-8)" << std::endl;
  std::cerr << "  -p, --coptol VAL         sets coplanarity tolerance to VAL (default = 1e-8)" << std::endl;
  std::cerr << "  -l, --coltol VAL         sets colinearity tolerance to VAL (default = 1e-8)" << std::endl;
  std::cerr << "  -s, --spltol VAL         sets the splitting tolerance to VAL (default = 1e-8)" << std::endl;
  std::cerr << "  -t, --strtol VAL         sets the straddling tolerance to VAL (default = 1e-8)" << std::endl;
  std::cerr << "  -r, --randomize          randomizes the order edges are split" << std::endl;
  std::cerr << "  -e, --split-degen        splits can results in degenerate halfspaces" << std::endl;
  std::cerr << "  -m, --max-tri-area VAL   maximum triangle area for generating 'extra' points" << std::endl;
}

/// Subdivides a mesh such that each triangle is no bigger than the maximum area
IndexedTriArray subdivide(const IndexedTriArray& mesh, double max_tri_area)
{
  // setup a list of triangles within the size limit and a queue of triangles
  // to process
  std::list<IndexedTri> good_tris;
  std::queue<IndexedTri> tri_q;

  // get triangles and vertices from the mesh
  std::vector<Point3d> verts = mesh.get_vertices();  
  const std::vector<IndexedTri>& facets = mesh.get_facets();  

  // copy triangles into the list or queue as appropriate
  for (unsigned i=0; i< facets.size(); i++)
  {
    Triangle t(verts[facets[i].a], verts[facets[i].b], verts[facets[i].c]);
    if (t.calc_area() > max_tri_area)
      tri_q.push(facets[i]);
    else
      good_tris.push_back(facets[i]);
  }

  // process triangles in the queue until the queue is empty
  while (!tri_q.empty())
  {
    // get the triangle off of the front of the queue
    IndexedTri t = tri_q.front();
    tri_q.pop();

    // determine the Steiner point 
    Point3d steiner = (verts[t.a] + verts[t.b] + verts[t.c])*0.333333333;

    // add it to the vertices
    verts.push_back(steiner);

    // create three new triangles
    IndexedTri t1(t.a, t.b, verts.size()-1);
    IndexedTri t2(t.b, t.c, verts.size()-1);
    IndexedTri t3(t.c, t.a, verts.size()-1);

    // if the area of the three triangles is sufficiently small, add them
    // to the good list; otherwise, add them to the queue
    Triangle tri(verts[t1.a], verts[t1.b], verts[t1.c]);
    if (tri.calc_area() < max_tri_area)
    {
      good_tris.push_back(t1);
      good_tris.push_back(t2);
      good_tris.push_back(t3);
    }
    else
    {
      tri_q.push(t1);
      tri_q.push(t2);
      tri_q.push(t3);
    }
  }

  // create new mesh with new facets
  return IndexedTriArray(verts.begin(), verts.end(), facets.begin(), facets.end());
}

/// Determines points from intersection of the BSP half-spaces
std::map<BSP*, std::list<Point3d> > determine_extra_points(const BSP& bsp, const IndexedTriArray& mesh)
{
  // setup the map
  std::map<BSP*, std::list<Point3d> > epoints;

  // subdivide the mesh
  IndexedTriArray subdiv_mesh = subdivide(mesh, MAX_TRI_AREA);

  // get vertices of the mesh
  const std::vector<Point3d>& verts = subdiv_mesh.get_vertices();

  // for every vertex, check to see whether the vertex is inside or on the
  // leaf; if it is, store it as an extra point for this leaf
  BOOST_FOREACH(const Point3d& v, verts)
  {
    // find the BSP nodes that this vertex 'belongs' to
    std::list<BSP*> nodes = bsp.find_nodes(v);
    assert(!nodes.empty());

    // add as extra points to these nodes
    for (std::list<BSP*>::const_iterator j = nodes.begin(); j != nodes.end(); j++)
      epoints[*j].push_back(v); 
  }

  return epoints;
}

/**
 * The decomposition procedure operates as follows:
 * 1. All non-convex edges, E, are identified
 * 2. Create new BSP holding all edges
 * 3. For each e in E
 * 3a. Find the BSP subnode(s) containing e; if none found, e has
 *     been removed by another splitting plane and we can continue #3
 * 3b. Split each BSP subnode using arbitrary splitting plane compatible w/e
 * 4. Intersect all half-spaces of BSP leaves (use bounding box + eps of the
 *     polyhedron as half-spaces too); result will be an intersection 
 *     polyhedron.  Test every vertex of the intersection polyhedron against
 *     all BSP leafs; if it is inside or on the leaf, add that vertex to the
 *     leaf's points
 */ 
int main(int argc, char* argv[])
{
  const unsigned X = 0, Y =1 , Z = 2;

  // setup option structure
  static struct option long_options[] = {
    {"debug", no_argument, 0, 'd'},
    {"split-degen", no_argument, 0, 'e'},
    {"ranomize", no_argument, 0, 'r'},
    {"cvxtol", required_argument, 0, 'v'},
    {"coptol", required_argument, 0, 'p'},
    {"coltol", required_argument, 0, 'l'},
    {"spltol", required_argument, 0, 's'},
    {"strtol", required_argument, 0, 't'},
    {"max-tri-area", required_argument, 0, 'm'},
    {0, 0, 0, 0}};

  // set options for getopt() (we don't want error messages printed by getopt())
  opterr = 0;
  int option_index = 0;  

  // process options
  int c;
  while ((c = getopt_long(argc, argv, "dv:p:l:s:erm:", long_options, &option_index)) != -1)
    switch (c)
    {
      case 'd':
        DEBUG = true;
        std::cout << "- debugging mode activated" << std::endl;
        break;

      case 'v':
        CONVEXITY_TOL = (double) atof(optarg);
        std::cout << "- convexity tolerance set to " << CONVEXITY_TOL << std::endl;
        break;

      case 'p':
        COPLANAR_TOL = (double) atof(optarg);
        std::cout << "- coplanar tolerance set to " << COPLANAR_TOL << std::endl;
        break;
      
      case 'l':
        COLINEAR_TOL = (double) atof(optarg);
        std::cout << "- colinear tolerance set to " << COLINEAR_TOL << std::endl;
        break;

      case 'm':
        MAX_TRI_AREA = (double) atof(optarg);
        std::cout << "- maximum triangle area set to " << MAX_TRI_AREA << std::endl;
        break;

      case 's':
        SPLIT_TOL = (double) atof(optarg);
        std::cout << "- splitting tolerance set to " << SPLIT_TOL << std::endl;
        break;

      case 't':
        EDGE_ON_PLANE_TOL = (double) atof(optarg);
        std::cout << "- straddling tolerance set to " << EDGE_ON_PLANE_TOL << std::endl;
        break;

      case 'r':
        RANDOMIZE = true;
        std::cout << "- randomizing splitting order" << std::endl;
        break;

      case 'e':
        SPLIT_DEGENERATE = true;
        std::cout << "- edges can be split into degenerate halfspaces" << std::endl;
        break;
  
      case '?':
        print_syntax();
          return -1;

      default:
        abort();
    }

  // verify that input argument is remaining
  if (optind >= argc)
  {
    print_syntax();
    return -1;
  }

  // read in the file
  IndexedTriArray t = IndexedTriArray::read_from_obj(std::string(argv[optind]));
  const std::vector<Point3d>& verts = t.get_vertices(); 
  const std::vector<IndexedTri>& tris = t.get_facets();

  // determine the root of the filename
  std::string fname_root = std::string(argv[optind]);
  fname_root = fname_root.substr(0,fname_root.find_last_of('.'));

  // determine all edges
  std::map<sorted_pair<const Point3d*>, EdgePtr> edges;
  for (unsigned i=0; i< tris.size(); i++)
  {
    // make sorted pairs from the vertices
    sorted_pair<const Point3d*> s1(&verts[tris[i].a], &verts[tris[i].b]);
    sorted_pair<const Point3d*> s2(&verts[tris[i].b], &verts[tris[i].c]);
    sorted_pair<const Point3d*> s3(&verts[tris[i].c], &verts[tris[i].a]);

    // try to find the edges in the set already
    EdgePtr e1 = edges[s1];
    EdgePtr e2 = edges[s2];
    EdgePtr e3 = edges[s3];

    // make the edges if they are not made already
    if (!e1)
    {
      e1 = EdgePtr(new Edge(&verts[tris[i].a], &verts[tris[i].b]));
      edges[s1] = e1;
    }
    if (!e2)
    {
      e2 = EdgePtr(new Edge(&verts[tris[i].b], &verts[tris[i].c]));
      edges[s2] = e2;
    }
    if (!e3)
    {
      e3 = EdgePtr(new Edge(&verts[tris[i].c], &verts[tris[i].a]));
      edges[s3] = e3;
    }

    // add the face to to the edges
    if (e1->f1)
    {
      if (e1->f2)
      {
        std::cerr << " warning: > 2 faces coincident to this edge!  skipping extra faces..." << std::endl;
        continue;
      }
      e1->f2 = shared_ptr<Tri>(new Tri(&verts[tris[i].a], &verts[tris[i].b], &verts[tris[i].c]));
    }
    else
      e1->f1 = shared_ptr<Tri>(new Tri(&verts[tris[i].a], &verts[tris[i].b], &verts[tris[i].c]));
    if (e2->f1)
    {
      if (e2->f2)
      {
        std::cerr << " warning: > 2 faces coincident to this edge!  skipping extra faces..." << std::endl;
        continue;
      }
      e2->f2 = shared_ptr<Tri>(new Tri(&verts[tris[i].a], &verts[tris[i].b], &verts[tris[i].c]));
    }
    else
      e2->f1 = shared_ptr<Tri>(new Tri(&verts[tris[i].a], &verts[tris[i].b], &verts[tris[i].c]));
    if (e3->f1)
    {
      if (e3->f2)
      {
        std::cerr << " warning: > 2 faces coincident to this edge!  skipping extra faces..." << std::endl;
        continue;
      }
      e3->f2 = shared_ptr<Tri>(new Tri(&verts[tris[i].a], &verts[tris[i].b], &verts[tris[i].c]));
    }
    else
      e3->f1 = shared_ptr<Tri>(new Tri(&verts[tris[i].a], &verts[tris[i].b], &verts[tris[i].c]));
  }

  // determine non-convex edges and setup vector of all edges
  std::set<EdgePtr> ncedges_set;
  std::vector<EdgePtr> edgesV;
  for (std::map<sorted_pair<const Point3d*>, EdgePtr>::const_iterator i = edges.begin(); i != edges.end(); i++)
  {
    edgesV.push_back(i->second);
    if (!i->second->is_convex(CONVEXITY_TOL))
      ncedges_set.insert(i->second);
  }

  // randomize non-convex edges if desired
  std::vector<EdgePtr> ncedges(ncedges_set.begin(), ncedges_set.end());
  if (RANDOMIZE)
    std::random_shuffle(ncedges.begin(), ncedges.end());

  // write non-convex edges to VRML
  if (DEBUG)
  {
    std::ofstream eout("ncedges.wrl");
    edges_to_vrml(eout, ncedges.begin(), ncedges.end());
    eout << "#VRML V2.0 utf8" << std::endl << std::endl;
    eout.close();
  }

  // output high-level info
  std::cout << ncedges.size() << " non-convex edges found" << std::endl;

  // open BSP log file
  if (DEBUG)
    BSP::_log.open("output.log");
  else
    BSP::_log.open("/dev/null");

  // add all edges to the root BSP
  BSP bsp;
  bsp.add_edges(edgesV);

  // continually split using all non-convex edges
  unsigned eidx = 0;
  for (std::vector<EdgePtr>::const_iterator i = ncedges.begin(); i != ncedges.end(); i++, eidx++)
  {
    // find BSP subspaces containing e
    std::list<BSP*> subspaces;
    bsp.get_subspaces(*i, subspaces);
    
    // output logging info
    BSP::_log << "edge " << *i << "(" << eidx << ") in " << subspaces.size() << " subspaces" << std::endl;
    for (std::list<BSP*>::const_iterator j = subspaces.begin(); j != subspaces.end(); j++)
      BSP::_log << "  subspace: " << *j << std::endl;

    // determine the plane that will effectively split the non-convex edge
    assert((*i)->f1);
    Vector3d normal = (*i)->f1->calc_normal();
    if ((*i)->f2)
    {
      normal -= (*i)->f2->calc_normal();
      normal.normalize();
    }
    double d = Vector3d::dot(normal, *(*i)->v.first);
    Plane plane(normal[X], normal[Y], normal[Z], d);

    // write non-convex edge to VRML
    if (DEBUG)
    {
      std::ostringstream str;
      str << "ncedges-" << eidx << ".wrl";
      std::ofstream eout(str.str().c_str());
      std::vector<EdgePtr>::const_iterator j = i;
      eout << "#VRML V2.0 utf8" << std::endl << std::endl;
      edges_to_vrml(eout, i, ++j);
      eout.close();
    }

    // split each subspace using the plane 
    unsigned ssidx = 0;
    for (std::list<BSP*>::const_iterator j = subspaces.begin(); j != subspaces.end(); j++, ssidx++)
    {
      // split the plane
      (*j)->split(plane);

      // verify plane split before generating debugging info
      if (!(*j)->is_leaf() && DEBUG)
      {
        // write the generated edges 
        const std::vector<EdgePtr>& pedges = (*j)->positive()->get_edges();
        const std::vector<EdgePtr>& nedges = (*j)->negative()->get_edges();
        std::ostringstream str;
        str << "edges-" << eidx << "-" << ssidx << ".wrl";
        std::ofstream eout(str.str().c_str());
        eout << "#VRML V2.0 utf8" << std::endl << std::endl;
        eout << "# edges from " << (*j)->positive() << std::endl;
        edges_to_vrml(eout, pedges.begin(), pedges.end());
        eout << std::endl;
        eout << "# edges from " << (*j)->negative() << std::endl;
        edges_to_vrml(eout, nedges.begin(), nedges.end());
        eout.close();

        // write the generated half-spaces 
        str.str("");
        str << "hs-" << eidx << "-" << ssidx << ".wrl";
        eout.open(str.str().c_str());
        eout << "#VRML V2.0 utf8" << std::endl << std::endl;
        (*j)->to_vrml(eout);
        eout.close();
      }
    }
  }

  // write out all half-spaces
  if (DEBUG && false)
  {
    std::ofstream eout("hs.wrl");
    eout << "#VRML V2.0 utf8" << std::endl << std::endl;
    bsp.to_vrml(eout);
    eout.close();
  }

  // determine all extra points
  BSP::_log << "===================================================" << std::endl;
  BSP::_log << " preparing to write half-spaces" << std::endl;
  BSP::_log << "===================================================" << std::endl;
  std::cout << "determining extra points" << std::endl;
  std::map<BSP*, std::list<Point3d> > extra_points = determine_extra_points(bsp, t);

  // get all BSP leaf nodes
  std::list<BSP*> leafs = bsp.get_leafs();
  std::cout << leafs.size() << " leaf spaces generated" << std::endl;

  // iterate over all leaf nodes
  unsigned index = 0;
  for (std::list<BSP*>::const_iterator i = leafs.begin(); i != leafs.end(); i++, index++)
  {
    if (DEBUG)
    {
      // write the edges from this BSP node
      const std::vector<EdgePtr>& edges = (*i)->get_edges();
      std::ostringstream str;
      str << "edges-" << index << ".wrl";
      std::ofstream out(str.str().c_str());
      out << "#VRML V2.0 utf8" << std::endl << std::endl;
      edges_to_vrml(out, edges.begin(), edges.end());
      out.close();
    }

    // get the points from the BSP node
    std::list<Point3d> points = (*i)->get_points();
    points.insert(points.end(), extra_points[*i].begin(), extra_points[*i].end());
    std::cout << "BSP #" << index << ": " << (*i) << " number of points: " << points.size();

    // make sure that the number of points is not too small
    if (points.size() < 4)
    {
      std::cout << " (< 4 points = degenerate; not writing)" << std::endl;
      std::cout << " points:" << std::endl;
      BOOST_FOREACH(const Point3d& pt, points)
        std::cout << "   " << pt << std::endl;
      continue;
    }

    // compute the convex hull of the points
    PolyhedronPtr p = CompGeom::calc_convex_hull(points.begin(), points.end());
    if (!p)
    {
      std::cout << " (unable to calculate hull = degenerate; not writing)" << std::endl;
      std::cout << " points:" << std::endl;
      BOOST_FOREACH(const Point3d& pt, points)
        std::cout << "   " << pt << std::endl;
      continue;
    }
    else
      std::cout << std::endl;

    // create the output filename
    std::ostringstream sstr;
    sstr << fname_root << "." << index << ".obj";

    // write the OBJ file
    p->get_mesh().write_to_obj(sstr.str());
  }

  // close the BSP log
  BSP::_log.close();

  // free memory
  BOOST_FOREACH(Point3d* v, new_verts)
    delete v;
}

