/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public
 * License (found in COPYING).
 ****************************************************************************/

#ifndef _PENALTY_EVENT_HANDLER_H
#define _PENALTY_EVENT_HANDLER_H

#include <list>
#include <vector>
#include <map>
#include <Ravelin/LinAlgd.h>
#include <Moby/Base.h>
#include <Moby/Types.h>
#include <Moby/LCP.h>
#include <Moby/Event.h>
#include <Moby/EventProblemData.h>

namespace Moby {

/// Defines the mechanism for handling Penalty events
class PenaltyEventHandler
{
  public:
    PenaltyEventHandler();
    void process_events(const std::vector<Event>& events) const;
  private:
    void apply_model(const std::vector<Event>& events) const;
    static double sqr(double x) { return x*x; }
}; // end class
} // end namespace

#endif

