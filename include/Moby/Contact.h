/****************************************************************************
 * Copyright 2006 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _MOBY_CONTACT_H
#define _MOBY_CONTACT_H

#include <iostream>
#include <list>
#include <boost/shared_ptr.hpp>
#include <Moby/Constants.h>
#include <Moby/Vector3.h>
#include <Moby/Event.h>

class SoNode;

namespace Moby {

class Polyhedron;
class CollisionGeometry;

/// Container class for describing a contact in the simulation        
class Contact : public Event
{
  public:
    enum ContactType { eSeparating, eResting, eImpacting };
    Contact();
    Contact(const Contact& c) { *this = c; }
    virtual ~Contact() {}
    void operator=(const Contact& c);    
    Real calc_normal_rvel() const;
    ContactType determine_contact_type(Real tol = (Real) 0.0) const;

    /// The contact type (none, vertex)
    ContactGeometryType contact_geometry;

    /// The point contact
    Vector3 point;
    
    /// The vector pointing outward from the contact on the first body (in world coordinates)
    Vector3 normal;  

    /// Force / impulse that has been applied due to this contact
    /**
     * The force / impulse applied to the body of the first geometry; 
     * the reverse of this force / impulse is applied to the body of the second 
     * geometry.
     */
    Vector3 force;

    /// The collision geometries involved in the contact
    CollisionGeometryPtr geom1, geom2;

    /// The coefficient of Coulomb friction
    Real mu_coulomb;

    /// The coefficient of viscous friction
    Real mu_viscous;

    /// The coefficient of restitution
    Real epsilon;

    SoNode* to_visualization_data() const;
    void write_vrml(const std::string& filename) const;
}; // end class

std::ostream& operator<<(std::ostream& o, const Contact& c);


} // end namespace Moby

#endif

