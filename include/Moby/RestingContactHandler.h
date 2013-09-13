#ifndef CONTACT_H
#define CONTACT_H

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

/// Defines the mechanism for handling impact events
class RestingContactHandler
{
  public:
    RestingContactHandler();
    void process_events(const std::vector<Event>& events);

  private:
    static DynamicBodyPtr get_super_body(SingleBodyPtr sb);
    void apply_model(const std::vector<Event>& events);
    void apply_model_to_connected_events(const std::list<Event*>& events);
    static void compute_problem_data(EventProblemData& epd);
    void solve_lcp(EventProblemData& epd, Ravelin::VectorNd& z);
    double calc_ke(EventProblemData& epd, const Ravelin::VectorNd& z);
    void apply_forces(const EventProblemData& epd) const;
    static void contact_select(const std::vector<int>& alpha_c_indices, const std::vector<int>& beta_nbeta_c_indices, const Ravelin::VectorNd& x, Ravelin::VectorNd& alpha_c, Ravelin::VectorNd& beta_c);
    static void contact_select(const std::vector<int>& alpha_c_indices, const std::vector<int>& beta_nbeta_c_indices, const Ravelin::MatrixNd& m, Ravelin::MatrixNd& alpha_c_rows, Ravelin::MatrixNd& beta_c_rows);
    static double sqr(double x) { return x*x; }

    Ravelin::LinAlgd _LA;
    LCP _lcp;
}; // end class
} // end namespace

#endif // CONTACT_H
