/****************************************************************************
 * Copyright 2007 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _SO_SEPARATOR_WRAPPER_H_
#define _SO_SEPARATOR_WRAPPER_H_

#include <Inventor/nodes/SoSeparator.h>
#include <Moby/Base.h>

class SoSeparator;

namespace Moby {

/// A wrapper for OpenInventor SoSeparator class, supporting serialization
class SoSeparatorWrapper : public virtual Base
{
  public:
    SoSeparatorWrapper();
    SoSeparatorWrapper(SoNode* n); 
    SoSeparatorWrapper(const std::string& filename);
    ~SoSeparatorWrapper();
    virtual void load_from_xml(XMLTreeConstPtr node, std::map<std::string, BasePtr>& id_map);
    virtual void save_to_xml(XMLTreePtr node, std::list<BaseConstPtr>& shared_objects) const;
    SoSeparator* get_separator() { return _separator; }

  private:
    SoSeparator* _separator;
}; // end class

} // end namespace

#endif

