/****************************************************************************
 * Copyright 2011 Evan Drumwright
 * This library is distributed under the terms of the Apache V2.0 
 * License (obtainable from http://www.apache.org/licenses/LICENSE-2.0).
 ****************************************************************************/

#ifndef _NQP_IPOPT_H
#define _NQP_IPOPT_H

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
class NQP_IPOPT : public Ipopt::TNLP 
{
  public:
    NQP_IPOPT();

    /// Particular solution (before optimization), final solution (after)
    Ravelin::VectorNd z;

    /// noninterpenetration jacobian block
    Ravelin::MatrixNd Cn_block;

    /// noninterpenetration vector
    Ravelin::SharedConstVectorNd Cn_v;

    /// noninterpenetration jacobian block
    Ravelin::SharedConstMatrixNd L_block;

    /// Nullspace for optimization
    Ravelin::MatrixNd R;

    /// Quadratic optimization matrix
    Ravelin::MatrixNd H;

    /// Linear optimization vector
    Ravelin::VectorNd c;

    /// pointer to the EventProblemData
    const EventProblemData* epd;

    /// Squared Coulomb friction coefficients for contacts (size N_TRUE_CONE)
    std::vector<double> mu_c;

    /// Squared viscous terms for contacts (size N_TRUE_CONE); this is squared viscous friction coefficient times squared tangential contact velocity
    std::vector<double> mu_visc;

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

    // number of nonzeros in each nonlinear constraint Hessian
    std::vector<unsigned> _nnz_h_con;

    // the objective Hessian nonzero values
    boost::shared_array<double> _h_obj;

    // the constraint Hessian nonzero values
    std::vector<boost::shared_array<double> > _h_con;
    
    // indices (mapping from 0..nnz in constraint hessian to index in additive Hessian array
    std::vector<std::vector<unsigned> > _h_con_indices;

    // number of constant components of the constraint Jacobian
    unsigned _cJac_constant;

    // constant components of constraint Jacobian
    boost::shared_array<double> _cJac;

    // nonzero indices of the constraint Jacobian
    boost::shared_array<unsigned> _cJac_iRow, _cJac_jCol;

    // nonzero indices of the Hessian 
    boost::shared_array<unsigned> _h_iRow, _h_jCol;

    // temporaries 
    Ravelin::VectorNd _w, _workv, _workv2, _x;
    Ravelin::MatrixNd _M;

    // is the problem dense or not?
    bool _dense;

    // tolerance for the Coulomb friction and kappa constraints
    double _tol;
}; // end class
} // end namespace

#endif

