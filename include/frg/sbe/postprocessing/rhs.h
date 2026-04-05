#pragma once

#include <mymath.h>
#include <frg/rhs_base.h>
#include <frg/sbe/state.h>
#include <frg/sbe/postprocessing/state.h>

/**
 *  \brief Functor to specify the rhs calculation for both the self-energy and the vertex
 *
 *  \params const state_vec, dfdl, const t
 *  \return dfdl (state_frg_sbe_t<Model> type) after evaluation of the rhs
 */

template<typename Model, typename state_t>
class rhs_postproc_sbe_t : public rhs_base_t<Model, state_t >
{
 public:
    using rhs_base_t = ::rhs_base_t<Model, state_t >;

    void operator() ( const state_t &postproc_state, state_t &out_postproc_state, const double t ){}; 	 

    void operator() ( const state_t& sbe_state, state_postproc_t<Model> &out_postproc_state, const double t ); 	    /**< Overload call operator to calculate full rhs */
	
    rhs_postproc_sbe_t();

    ~rhs_postproc_sbe_t();

    void test( const state_t& state, state_postproc_t<Model>& postproc_state, const double t );
    /*
    void rhs_susc_sc_split_contributions( const state_t& sbe_state, state_postproc_t<Model> &out_postproc_state, const double t );

    void rhs_susc_d_split_contributions( const state_t& sbe_state, state_postproc_t<Model> &out_postproc_state, const double t );

    void rhs_susc_m_split_contributions( const state_t& sbe_state, state_postproc_t<Model> &out_postproc_state, const double t );

    void rhs_polarisation_sc_split_contributions( const state_t& sbe_state, state_postproc_t<Model> &out_postproc_state, const double t );

    void rhs_polarisation_d_split_contributions( const state_t& sbe_state, state_postproc_t<Model> &out_postproc_state, const double t );

    void rhs_polarisation_m_split_contributions( const state_t& sbe_state, state_postproc_t<Model> &out_postproc_state, const double t );
    */

 protected:

    static dcomplex eval_rhs_Sig( const idx_1p_t<Model>& se_idx, const state_t& state_vec, const gf_1p_mat_t<Model>& Gvec, double t );      /**< Return the standard SE flow equation */
    static dcomplex eval_rhs_Sig_LBG( const idx_1p_t<Model>& se_idx, const gf_lambda_t<Model>& B_sc, const gf_lambda_t<Model>& B_d, const gf_lambda_t<Model>& B_m, const gf_1p_mat_t<Model>& Gvec, double t );                     /**< Return self-energy rhs proportional to GG phi dG. */
    static dcomplex eval_rhs_Sig_LBG_sc( const idx_1p_t<Model>& se_idx, const gf_lambda_t<Model>& B_sc, const gf_1p_mat_t<Model>& Gvec, double t );
    static dcomplex eval_rhs_Sig_LBG_d( const idx_1p_t<Model>& se_idx, const gf_lambda_t<Model>& B_d, const gf_1p_mat_t<Model>& Gvec, double t );
    static dcomplex eval_rhs_Sig_LBG_m( const idx_1p_t<Model>& se_idx, const gf_lambda_t<Model>& B_m, const gf_1p_mat_t<Model>& Gvec, double t );

    static dcomplex eval_susc_verttypefunc( const state_t&, const idx_susc_t<Model>& idx, dcomplex (*verttypefunc)( const state_t&,int,int,int,int,int,int, const double), const gf_bubble_mat_t<Model>& bubble_pp, double t ); // blueprint for evaluating susc type functions

    static dcomplex eval_susc_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t ); /**< Return full postprocessing susceptibility for gf_susc_sc */        
    static dcomplex eval_suscbubble_sc( const idx_susc_t<Model>& idx,  const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_suscvert_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_bare( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_nabla_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_M_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_d_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_nabla_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_M_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_m_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_nabla_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_M_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_sc_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );
    static dcomplex eval_susc_sc_contribution_from_Virr( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );    


    static dcomplex eval_susc_d( const idx_susc_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );/**< Return 1loop flow equation for gf_susc_d */
    static dcomplex eval_suscbubble_d( const idx_susc_t<Model>& idx,  const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_suscvert_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_bare( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_nabla_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_M_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_d_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_nabla_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_M_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_m_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_nabla_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_M_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_d_contribution_from_sc_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );    
    static dcomplex eval_susc_d_contribution_from_Virr( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );    


    static dcomplex eval_susc_m( const idx_susc_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );/**< Return 1loop flow equation for gf_susc_m */
    static dcomplex eval_suscbubble_m( const idx_susc_t<Model>& idx,  const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_suscvert_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_bare( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_nabla_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_M_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_d_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_nabla_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_M_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_m_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_nabla_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_M_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_sc_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );
    static dcomplex eval_susc_m_contribution_from_Virr( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );

    // lambda eval functions
    static dcomplex eval_lambda_verttypefunc( const state_t& state_vec, const idx_lambda_t<Model>& idx, dcomplex (*verttypefunc)( const state_t&,int,int,int,int,int,int, const double), const gf_bubble_mat_t<Model>& bubble, double bubble_sign, double t);
    // same as above, but for static lambda
    static dcomplex eval_lambda_verttypefunc( const state_t& state_vec, const idx_static_lambda_t<Model>& idx, dcomplex (*verttypefunc)( const state_t&,int,int,int,int,int,int, const double), const gf_bubble_mat_t<Model>& bubble, double bubble_sign, double t);

    static dcomplex eval_lambda_contribution_from_1( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble , double t);


    // d
    static dcomplex eval_lambda_d_contribution_from_bare( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_d_contribution_from_nabla_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);
    
    static dcomplex eval_lambda_d_contribution_from_nabla_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_d_contribution_from_nabla_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);
  
    static dcomplex eval_lambda_d_contribution_from_M_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);
    
    static dcomplex eval_lambda_d_contribution_from_M_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_d_contribution_from_M_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_d_contribution_from_d_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_d_contribution_from_sc_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_d_contribution_from_m_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    // m
    static dcomplex eval_lambda_m_contribution_from_bare( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_m_contribution_from_nabla_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);
    
    static dcomplex eval_lambda_m_contribution_from_nabla_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_m_contribution_from_nabla_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_m_contribution_from_M_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);
    
    static dcomplex eval_lambda_m_contribution_from_M_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_m_contribution_from_M_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);
  
    static dcomplex eval_lambda_m_contribution_from_d_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_m_contribution_from_sc_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_m_contribution_from_m_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);


    // sc
    static dcomplex eval_lambda_sc_contribution_from_bare( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

    static dcomplex eval_lambda_sc_contribution_from_nabla_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);
    
    static dcomplex eval_lambda_sc_contribution_from_nabla_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

    static dcomplex eval_lambda_sc_contribution_from_nabla_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

  static dcomplex eval_lambda_sc_contribution_from_M_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);
    
    static dcomplex eval_lambda_sc_contribution_from_M_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

    static dcomplex eval_lambda_sc_contribution_from_M_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

  
    static dcomplex eval_lambda_sc_contribution_from_d_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

    static dcomplex eval_lambda_sc_contribution_from_sc_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

    static dcomplex eval_lambda_sc_contribution_from_m_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

    static dcomplex eval_lambda_sc_contribution_from_varphi( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

    static dcomplex eval_lambda_d_contribution_from_varphi( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_lambda_m_contribution_from_varphi( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);


    // polarisation
    static dcomplex eval_polarisation_verttypefunc( const state_t& state_vec, const idx_polarisation_t<Model>& idx, const gf_lambda_t<Model> &lambda, const gf_bubble_mat_t<Model>& bubble, double bubble_sign, double t);

    static dcomplex eval_polarisation_sc_contribution_lambda( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t, const gf_lambda_t<Model> &lambda);

    static dcomplex eval_polarisation_d_contribution_lambda( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t, const gf_lambda_t<Model> &lambda);

    static dcomplex eval_polarisation_m_contribution_lambda( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t, const gf_lambda_t<Model> &lambda);

    static dcomplex eval_polarisation_contribution_suscvert( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t, const gf_susc_t<Model> &susc_contribution);

    static dcomplex eval_polarisation_m_or_d_contribution_from_bubble( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_polarisation_sc_contribution_from_bubble( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

    static dcomplex eval_polarisation_sc_contribution_from_varphi( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);

    static dcomplex eval_polarisation_d_contribution_from_varphi( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);
    
    static dcomplex eval_polarisation_m_contribution_from_varphi( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_polarisation_sc_contribution_from_nabla_sc( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t);
    
    static dcomplex eval_polarisation_d_contribution_from_nabla_d( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex eval_polarisation_m_contribution_from_nabla_m( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t);

    static dcomplex vertex_sc( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );         /**< Return two-particle vertex-tensor element. PP notation. */

    static dcomplex vertex_d( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );        

    static dcomplex vertex_m( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t); 

    static dcomplex B_irreducible_vertex_sc( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );         /**< Return two-particle vertex-tensor element. PP notation. */

    static dcomplex B_irreducible_vertex_d( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );        

    static dcomplex B_irreducible_vertex_m( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t); 

    
    static dcomplex vertx_sc_contribution_from_bare( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );         /**< Return two-particle vertex-tensor elem    ent. PP notation. contribution from bare vertex */
    static dcomplex vertx_sc_contribution_from_nabla_sc( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_sc_contribution_from_M_sc( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_sc_contribution_from_sc_double_counting( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_sc_contribution_from_nabla_d( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );  
    static dcomplex vertx_sc_contribution_from_M_d( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_sc_contribution_from_d_double_counting( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_sc_contribution_from_nabla_m( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_sc_contribution_from_M_m( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_sc_contribution_from_m_double_counting( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );         /**< these are the double counting terms that need to be added to nabla and lambda to give parquet's Phi (GG irreducible vertex) */
    static dcomplex vertx_sc_contribution_from_Virr( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    // needed for lambda and the polarisation
    static dcomplex B_irreducible_vertex_sc_contribution_from_nabla_sc( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    
    

    static dcomplex vertx_d_contribution_from_bare( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );         /**< Return two-particle vertex-tensor elem    ent. PP notation. contribution from bare vertex */
    static dcomplex vertx_d_contribution_from_nabla_sc( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t ); 
    static dcomplex vertx_d_contribution_from_M_sc( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t ); 
    static dcomplex vertx_d_contribution_from_sc_double_counting( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_d_contribution_from_nabla_d( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_d_contribution_from_M_d( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t ); 
    static dcomplex vertx_d_contribution_from_d_double_counting( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_d_contribution_from_nabla_m( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out , const double t);
    static dcomplex vertx_d_contribution_from_M_m( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t ); 
    static dcomplex vertx_d_contribution_from_m_double_counting( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );         /**< Return two-particle vertex-tensor element. PP notation. contribution from magnetic channel */
    static dcomplex vertx_d_contribution_from_Virr( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    // needed for lambda and polarisation
    static dcomplex B_irreducible_vertex_d_contribution_from_nabla_d( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t ); 
    

    static dcomplex vertx_m_contribution_from_bare( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );         /**< Return two-particle vertex-tensor elem    ent. PP notation. contribution from bare vertex */
    static dcomplex vertx_m_contribution_from_nabla_d( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t ); 
    static dcomplex vertx_m_contribution_from_M_d( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t ); 
    static dcomplex vertx_m_contribution_from_d_double_counting( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_m_contribution_from_nabla_m( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_m_contribution_from_M_m( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t ); 
    static dcomplex vertx_m_contribution_from_m_double_counting( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );         /**< Return two-particle vertex-tensor element. PP notation. contribution from magnetic channel */
    static dcomplex vertx_m_contribution_from_nabla_sc( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_m_contribution_from_M_sc( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );  
    static dcomplex vertx_m_contribution_from_sc_double_counting( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    static dcomplex vertx_m_contribution_from_Virr( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t );
    // needed for lambda and polarisation
    static dcomplex B_irreducible_vertex_m_contribution_from_nabla_m( const state_t&, const int W, const int w_in, const int w_out, const int K, const int n_in, const int n_out, const double t ); 
    

    static dcomplex eval_rhs_Sig_Lam2PI( const idx_1p_t<Model>& se_idx, const gf_Sig_kMat_t<Model> gf_Sig_Lam2PI_kMat, double t ); /**< Return 1st MLOOP correction to std SE flow eq. */
    static dcomplex eval_B_G0_phi_sc( const idx_lambda_t<Model>& se_idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t); /**< Return B_{pp}(\Gamma_0,\phi_{pp)} */
    static dcomplex eval_B_G0_phi_d( const idx_lambda_t<Model>& se_idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t); /**< Return B_{ph}(\Gamma_0,\phi_{ph}) */
    static dcomplex eval_B_G0_phi_m( const idx_lambda_t<Model>& se_idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t); /**< Return B_{xph}(\Gamma_0,\phi_{xph}) */


  static dcomplex eval_Sig_SDE_integrand_sc( const idx_Sig_SDE_integrand_t<Model>& idx ,const state_t&, const gf_1p_mat_t<Model>& Gvec, const double t );
  static dcomplex eval_Sig_SDE_integrand_d( const idx_Sig_SDE_integrand_t<Model>& idx ,const state_t&, const gf_1p_mat_t<Model>& Gvec, const double t );
  static dcomplex eval_Sig_SDE_integrand_m( const idx_Sig_SDE_integrand_t<Model>& idx ,const state_t&, const gf_1p_mat_t<Model>& Gvec, const double t );
};

// include implementation 
#include <../src/frg/sbe/postprocessing/rhs.tpp> // the logic of postprocessing
#include <../src/frg/sbe/postprocessing/rhs_eval.tpp> // the equations for the different contributions
