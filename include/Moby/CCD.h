/****************************************************************************
 * Copyright 2009 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _CCD_H
#define _CCD_H

#include <list>
#include <set>
#include <map>
#include <boost/shared_ptr.hpp>
#include <Moby/sorted_pair>
#include <Moby/Log.h>
#include <Moby/BV.h>

namespace Moby {

class RigidBody;
class ArticulatedBody;
class CollisionGeometry;  

/// Implements the CollisionDetection abstract class to perform exact contact finding using abstract shapes 
class CCD
{
  public:
    CCD();
    virtual ~CCD() {}
    virtual void load_from_xml(boost::shared_ptr<const XMLTree> node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<boost::shared_ptr<const Base> >& shared_objects) const;
    double calc_CA_step(DynamicBodyPtr dbA, DynamicBodyPtr dbB);
    void find_contacts(DynamicBodyPtr dbA, DynamicBodyPtr dbB, std::list<Event>& events);

  private:
    // temporary vectors
    Ravelin::VectorNd _workv, _l, _u, _c, _x;

    void find_contacts_separated(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::list<Event>& e);
    void find_contacts_not_separated(CollisionGeometryPtr cgA, CollisionGeometryPtr cgB, std::list<Event>& e);
    double calc_CA_step(RigidBodyPtr rbA, RigidBodyPtr dbB);
    double calc_max_dist_per_t(RigidBodyPtr rb, const Ravelin::Vector3d& n, const Ravelin::Vector3d& r);
    double solve_lp(const Ravelin::VectorNd& c, const Ravelin::VectorNd& l, const Ravelin::VectorNd& u, Ravelin::VectorNd& x);    
    static void to_binary(unsigned i, const Ravelin::VectorNd& l, const Ravelin::VectorNd& u, Ravelin::VectorNd& x); 
    bool intersect_BV_trees(boost::shared_ptr<BV> a, boost::shared_ptr<BV> b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b);
    static Event create_contact(CollisionGeometryPtr a, CollisionGeometryPtr b, const Point3d& point, const Ravelin::Vector3d& normal);

    template <class OutputIterator>
    OutputIterator intersect_BV_leafs(BVPtr a, BVPtr b, const Ravelin::Transform3d& aTb, CollisionGeometryPtr geom_a, CollisionGeometryPtr geom_b, OutputIterator output_begin) const;
}; // end class

} // end namespace

#endif

