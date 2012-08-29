/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the GNU Lesser General Public 
 * License (found in COPYING).
 ****************************************************************************/

/// Gets the time for the first contact with a joint limit (assuming Euler integration, so this is only first-order accurate)
template <class OutputIterator>
OutputIterator ArticulatedBody::find_limit_events(const VectorN& q0, const VectorN& q1, Real dt, OutputIterator output_begin) 
{
  SAFESTATIC VectorN dq;

  // compute the generalized velocity that takes us from q0 to q1
  dq.copy_from(q1) -= q0;
  set_generalized_coordinates(eRodrigues, q0);
  set_generalized_velocity(eRodrigues, dq);

  for (unsigned i=0; i< _joints.size(); i++)
    for (unsigned j=0; j< _joints[i]->num_dof(); j++)
    {
      // get the current joint position and velocity
      Real q = _joints[i]->q[j];
      Real qd = _joints[i]->qd[j];

      // setup an event for this joint/dof in case we need it
      Event e;
      e.event_type = Event::eLimit;
      e.limit_joint = _joints[i];
      e.limit_dof = j;
      e.limit_epsilon = _joints[i]->limit_restitution;

      // check whether we are already at a limit
      if (q >= _joints[i]->hilimit[j])
      {
        // add event for upper limit
        e.t = (Real) 0.0;
        e.limit_upper = true;
        *output_begin++ = e;

        // check whether lower limit is also an event
        if (_joints[i]->qd[j] < (Real) 0.0)
        {
          Real toc = (_joints[i]->lolimit[j] - q)/qd;
          if (toc < (Real) 1.0)
          {
            e.t = toc;
            e.limit_upper = false;
            *output_begin++ = e; 
          }
        }
      }
      else if (q <= _joints[i]->lolimit[j])
      {
        e.t = (Real) 0.0;
        e.limit_upper = false;
        *output_begin++ = e;

        // check whether upper limit is also an event
        if (qd > (Real) 0.0)
        {
          Real toc = (_joints[i]->hilimit[j] - q)/qd;
          if (toc < (Real) 1.0)
          {
            e.t = toc;
            e.limit_upper = true;
            *output_begin++ = e; 
          }
        }
      }
      else
      {
        // only check appropriate limit
        if (qd > (Real) 0.0)
        {
          Real toc = (_joints[i]->hilimit[j] - q)/qd;
          if (toc < (Real) 1.0)
          {
            e.t = std::max((Real) 0.0, toc);
            e.limit_upper = true;
            *output_begin++ = e;
          }
        }
        else if (qd < (Real) 0.0)
        {
          Real toc = (_joints[i]->lolimit[j] - q)/qd;
          if (toc < (Real) 1.0)
          {
            e.t = std::max((Real) 0.0, toc);
            e.limit_upper = false;
            *output_begin++ = e;
          }
        }
      }
    }

  return output_begin;
}


