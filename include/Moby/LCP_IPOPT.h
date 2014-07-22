/****************************************************************************
 * Copyright 2013 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _LCP_IPOPT_H
#define _LCP_IPOPT_H

#include <list>
#include <vector>
#include <map>
#include <coin/IpTNLP.hpp>
#include <Ravelin/LinAlgd.h>
#include <Moby/Base.h>
#include <Moby/Types.h>
#include <Moby/LCP.h>
#include <Moby/Event.h>
#include <Moby/EventProblemData.h>

namespace Moby {

/// Defines the mechanism for handling impact events 
class LCP_IPOPT : public Ipopt::TNLP 
{
  public:
    LCP_IPOPT();

    /// Particular solution (before optimization), final solution (after)
    Ravelin::VectorNd z;

    /// noninterpenetration jacobian block
    Ravelin::SharedConstMatrixNd Cn_block;

    /// noninterpenetration jacobian block
    Ravelin::SharedConstMatrixNd L_block;

    /// LCP matrix
    Ravelin::MatrixNd M;

    /// LCP vector 
    Ravelin::VectorNd q;

    /// pointer to the EventProblemData
    const EventProblemData* epd;

  public:
    // virtual methods for IPOPT-based optimization
    virtual bool get_nlp_info(Ipopt::Index& n, Ipopt::Index& m, Ipopt::Index& nnz_jac_g, Ipopt::Index& nnz_h_lag, Ipopt::TNLP::IndexStyleEnum& index_style);   
    virtual bool get_bounds_info(Ipopt::Index n, Ipopt::Number* x_l, Ipopt::Number* x_u, Ipopt::Index m, Ipopt::Number* g_l, Ipopt::Number* g_u);
    virtual bool get_starting_point(Ipopt::Index n, bool init_x, Ipopt::Number* x, bool init_z, Ipopt::Number* z_L, Ipopt::Number* z_U, Ipopt::Index m, bool init_lambda, Ipopt::Number* lambda);
    virtual bool eval_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number& obj_value);
    virtual bool eval_grad_f(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number* grad_f);
    virtual bool eval_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Number* g);
    virtual bool eval_jac_g(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Index m, Ipopt::Index nele_jac, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values);
    virtual bool eval_h(Ipopt::Index n, const Ipopt::Number* x, bool new_x, Ipopt::Number obj_factor, Ipopt::Index m, const Ipopt::Number* lambda, bool new_lambda, Ipopt::Index nele_hess, Ipopt::Index* iRow, Ipopt::Index* jCol, Ipopt::Number* values);
    virtual void finalize_solution(Ipopt::SolverReturn status, Ipopt::Index n, const Ipopt::Number* x, const Ipopt::Number* z_L, const Ipopt::Number* z_U, Ipopt::Index, const Ipopt::Number* g, const Ipopt::Number* lambda, Ipopt::Number obj_value, const Ipopt::IpoptData* ip_data, Ipopt::IpoptCalculatedQuantities* ip_cq);

    // tolerance for checking for zero
    double _ZERO_EPS;

    // number of nonzeros in objective Hessian
    unsigned _nnz_h_obj;

    // the objective Hessian nonzero values
    boost::shared_array<double> _h_obj;

    // nonzero indices of the Hessian 
    boost::shared_array<unsigned> _h_iRow, _h_jCol;

    // temporaries 
    Ravelin::VectorNd _workv;
}; // end class
} // end namespace

#endif

