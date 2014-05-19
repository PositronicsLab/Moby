#ifndef RESTINGCONTACTHANDLER_H
#define RESTINGCONTACTHANDLER_H

#include <list>
#include <vector>
#include <map>
#include <Ravelin/LinAlgd.h>
#include <Moby/Base.h>
#include <Moby/Types.h>
#include <Moby/LCP.h>
#include <Moby/Event.h>
#include <Moby/AccelerationEventData.h>

namespace Moby {

/// Defines the mechanism for handling impact contacts
class AccelerationEventHandler
{
  public:
    AccelerationEventHandler();
    void process_events(const std::vector<Event>& contacts);

  private:
    static DynamicBodyPtr get_super_body(SingleBodyPtr sb);
    void apply_model(const std::vector<Event>& contacts);
    void apply_model_to_connected_contacts(const std::list<Event*>& contacts);
    static void compute_problem_data(AccelerationEventData& epd);
    bool solve_lcp(AccelerationEventData& epd, Ravelin::VectorNd& z);
    double calc_ke(AccelerationEventData& epd, const Ravelin::VectorNd& z);
    void apply_forces(const AccelerationEventData& epd) const;
    static void contact_select(const std::vector<int>& alpha_c_indices, const std::vector<int>& beta_nbeta_c_indices, const Ravelin::VectorNd& x, Ravelin::VectorNd& alpha_c, Ravelin::VectorNd& beta_c);
    static void contact_select(const std::vector<int>& alpha_c_indices, const std::vector<int>& beta_nbeta_c_indices, const Ravelin::MatrixNd& m, Ravelin::MatrixNd& alpha_c_rows, Ravelin::MatrixNd& beta_c_rows);
    static double sqr(double x) { return x*x; }
    bool solve_lcp(const Ravelin::MatrixNd& M, const Ravelin::VectorNd& q, Ravelin::VectorNd& z);

    Ravelin::LinAlgd _LA;
    LCP _lcp;
}; // end class

} // end namespace

#endif // RESTINGCONTACTHANDLER_H
