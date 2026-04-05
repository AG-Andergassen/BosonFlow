#pragma once

#include <mymath.h>
#include <frg/sbe/state.h>
#include <frg/rhs_base.h>

/**
 *  \brief Functor to specify the rhs calculation for both the self-energy and the vertex
 *
 *  \params const state_vec, dfdl, const t
 *  \return dfdl (state_frg_sbe_t<Model> type) after evaluation of the rhs
 */
#include <frg/sbe/state_bosonised_M.h>


template <typename Model, typename state_t >
    class rhs_sbe_1lfrg_t : public rhs_base_t<Model, state_t >
{
 public:
    using rhs_base_t = ::rhs_base_t<Model, state_t >;

    void operator() ( const state_t& old_state, state_t &dstate_over_dt, const double t ); 	    /**< Overload call operator to calculate full rhs */
	
 rhs_sbe_1lfrg_t():
    m_lambda_sc_left_Fdot_part(),
    m_lambda_d_left_Fdot_part(),
    m_lambda_m_left_Fdot_part(),
    m_M_sc_left_Fdot_part(),
    m_M_d_left_Fdot_part(),
    m_M_m_left_Fdot_part(){}

    ~rhs_sbe_1lfrg_t();

 protected:
    void rhs_frg_1l(state_t *dstate_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const gf_1p_mat_t<Model>& Svec, const double t);

    // screened interaction w rhs
    void rhs_w_sc_1l(gf_w_t<Model> *new_w_sc_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);
    void rhs_w_d_1l(gf_w_t<Model> *new_w_d_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);
    void rhs_w_m_1l(gf_w_t<Model> *new_w_m_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);
    
    // Hedin vertex 
    void rhs_lambda_sc_1l(gf_lambda_t<Model> *new_lambda_sc_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);
  void rhs_lambda_d_1l(gf_lambda_t<Model> *new_lambda_d_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);
    void rhs_lambda_m_1l(gf_lambda_t<Model> *new_lambda_m_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    // self energy
    void rhs_Sig_1l(gf_1p_t<Model> *new_Sig_ptr, const state_t &old_state, const gf_1p_mat_t<Model>& Svec, const double t);

    void add_katanin_correction(gf_1p_mat_t<Model> *Svec_ptr, gf_1p_mat_t<Model> &Gvec, gf_1p_mat_t<Model> &Svec_init, state_t &old_state);

    // rest function
    void rhs_M_sc_1l(gf_M_t<Model> *new_M_sc_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);
    void rhs_M_d_1l(gf_M_t<Model> *new_M_d_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);
    void rhs_M_m_1l(gf_M_t<Model> *new_M_m_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void rhs_wM_sc_1l(gf_w_M_t<Model> *new_wM_sc_plus_ptr, gf_w_M_t<Model> *new_wM_sc_minus_ptr, const state_frg_sbe_bosonised_M_t<Model> &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);
    void rhs_wM_d_1l(gf_w_M_t<Model> *new_wM_d_plus_ptr, gf_w_M_t<Model> *new_wM_d_minus_ptr, const state_frg_sbe_bosonised_M_t<Model> &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);
    void rhs_wM_m_1l(gf_w_M_t<Model> *new_M_m_plus_ptr, gf_w_M_t<Model> *new_M_m_minus_ptr,  const state_frg_sbe_bosonised_M_t<Model> &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void rhs_lambdaM_sc_1l(gf_lambda_M_t<Model> *new_lambdaM_sc_plus_ptr, gf_lambda_M_t<Model> *new_lambdaM_sc_minus_ptr, const state_frg_sbe_bosonised_M_t<Model> &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);
    void rhs_lambdaM_d_1l(gf_lambda_M_t<Model> *new_lambdaM_d_plus_ptr, gf_lambda_M_t<Model> *new_lambdaM_d_minus_ptr, const state_frg_sbe_bosonised_M_t<Model> &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);
    void rhs_lambdaM_m_1l(gf_lambda_M_t<Model> *new_M_m_plus_ptr, gf_lambda_M_t<Model> *new_M_m_minus_ptr,  const state_frg_sbe_bosonised_M_t<Model> &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void rhs_frg_bosonised_M_1l(state_t *dstate_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

 private:
    
    gf_lambda_t<Model> m_lambda_sc_left_Fdot_part;
    gf_lambda_t<Model> m_lambda_d_left_Fdot_part;
    gf_lambda_t<Model> m_lambda_m_left_Fdot_part;
    gf_M_t<Model> m_M_sc_left_Fdot_part;
    gf_M_t<Model> m_M_d_left_Fdot_part;
    gf_M_t<Model> m_M_m_left_Fdot_part;


    static dcomplex eval_Sig_conv( const idx_1p_t<Model>& idx, const state_t& state_vec, const gf_1p_mat_t<Model>& Svec, double t);

    dcomplex eval_w_sc( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_pp, double t );	/**< Returns 1loop flow equation for the screened interaction in the U-pp-reducible channel w_pp */
    dcomplex eval_w_d( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph, double t );	/**< Returns 1loop flow equation for the screened interaction in the U-ph-reducible channel w_d */
    dcomplex eval_w_m( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph, double t );	/**< Returns 1loop flow equation for the screened interaction in the U-xph-reducible channel w_m */
        
    dcomplex eval_lambda_sc( const idx_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_pp, double t ); /**< Returns 1loop flow equation for the Hedin vertex in the U-pp-reducible channel lambda_sc */
    dcomplex eval_lambda_d( const idx_lambda_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph, double t );	/**< Returns 1loop flow equation for the Hedin vertex in the U-ph-reducible channel lambda_d */
    dcomplex eval_lambda_m( const idx_lambda_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph, double t );	/**< Returns 1loop flow equation for the Hedin vertex in the U-xph-reducible channel lambda_m */


    dcomplex eval_M_sc( const idx_M_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_pp, double t );  /**< Returns 1loop flow equation for the U-irreducible rest function M_sc */
    dcomplex eval_M_d( const idx_M_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph, double t );   /**< Returns 1loop flow equation for the U-irreducible rest function M_d */
    dcomplex eval_M_m( const idx_M_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph, double t );   /**< Returns 1loop flow equation for the U-irreducible rest function M_m */

  // additional contributions that come about due to a flowing interaction
    dcomplex eval_w_sc_Udot_part( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t );	/**< Returns 1loop flow equation for the screened interaction in the U-pp-reducible channel w_pp */
    dcomplex eval_w_d_Udot_part( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );	/**< Returns 1loop flow equation for the screened interaction in the U-ph-reducible channel w_d */
    dcomplex eval_w_m_Udot_part( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );	/**< Returns 1loop flow equation for the screened interaction in the U-xph-reducible channel w_m */

    dcomplex eval_lambda_sc_Udot_part( const idx_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t ); /**< Returns 1loop flow equation for the Hedin vertex in the U-pp-reducible channel lambda_sc */
    dcomplex eval_lambda_d_Udot_part( const idx_lambda_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );	/**< Returns 1loop flow equation for the Hedin vertex in the U-ph-reducible channel lambda_d */
    dcomplex eval_lambda_m_Udot_part( const idx_lambda_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );	/**< Returns 1loop flow equation for the Hedin vertex in the U-xph-reducible channel lambda_m */


    dcomplex eval_M_sc_Udot_part( const idx_M_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t );  /**< Returns 1loop flow equation for the U-irreducible rest function M_sc */
    dcomplex eval_M_d_Udot_part( const idx_M_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );   /**< Returns 1loop flow equation for the U-irreducible rest function M_d */
    dcomplex eval_M_m_Udot_part( const idx_M_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );   /**< Returns 1loop flow equation for the U-irreducible rest function M_m */

    // only precalculated left corrections
    void precalculate_lambda_sc_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);

    void precalculate_lambda_d_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void precalculate_lambda_m_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    // M (only precalculated left corrections)
    void precalculate_M_sc_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);

    void precalculate_M_d_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void precalculate_M_m_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

  
    // lambda (Hedin vertex)
    dcomplex eval_lambda_sc_left_Fdot_part( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t );

    dcomplex eval_lambda_d_left_Fdot_part( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    dcomplex eval_lambda_m_left_Fdot_part( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    // M (SBE Rest function)
    dcomplex eval_M_sc_left_Fdot_part( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t );

    dcomplex eval_M_d_left_Fdot_part( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    dcomplex eval_M_m_left_Fdot_part( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );



    // modified flow equations for the SBEb approximation (where M is neglected *b*efore taking the derivative)

    static dcomplex eval_SBEb_w_sc( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t );	/**< Returns 1loop flow equation for the screened interaction in the U-pp-reducible channel w_pp within the SBEb approximation*/
    static dcomplex eval_SBEb_w_d( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );	/**< Returns 1loop flow equation for the screened interaction in the U-ph-reducible channel w_d within the SBEb approximation */
    static dcomplex eval_SBEb_w_m( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );	/**< Returns 1loop flow equation for the screened interaction in the U-xph-reducible channel w_m within the SBEb approximation */
        
    static dcomplex eval_SBEb_lambda_sc( const idx_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t ); /**< Returns 1loop flow equation for the Hedin vertex in the U-pp-reducible channel lambda_sc within the SBEb approximation */
    static dcomplex eval_SBEb_lambda_d( const idx_lambda_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );	/**< Returns 1loop flow equation for the Hedin vertex in the U-ph-reducible channel lambda_d within the SBEb approximation */
    static dcomplex eval_SBEb_lambda_m( const idx_lambda_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t );	/**< Returns 1loop flow equation for the Hedin vertex in the U-xph-reducible channel lambda_m within the SBEb approximation */


    static dcomplex eval_wM_sc( const idx_w_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_pp, int sign, double t );  /**< Returns 1loop flow equation for the U-irreducible rest function M_sc */
    static dcomplex eval_wM_d( const idx_w_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph, int sign, double t );   /**< Returns 1loop flow equation for the U-irreducible rest function M_d */
    static dcomplex eval_wM_m( const idx_w_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph,  int sign, double t );   /**< Returns 1loop flow equation for the U-irreducible rest function M_m */

    static dcomplex eval_lambdaM_sc( const idx_lambda_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_pp, int sign, const gf_w_M_t<Model> &wM_dot, double t );  /**< Returns 1loop flow equation for the U-irreducible rest function M_sc */
    static dcomplex eval_lambdaM_d( const idx_lambda_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph, int sign, const gf_w_M_t<Model> &wM_dot, double t );   /**< Returns 1loop flow equation for the U-irreducible rest function M_d */
    static dcomplex eval_lambdaM_m( const idx_lambda_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state_vec, const gf_bubble_mat_t<Model>& bubble_GS_ph, int sign, const gf_w_M_t<Model> &wM_dot, double t );   /**< Returns 1loop flow equation for the U-irreducible rest function M_m */
 
    inline static const int nu_0(const int W){
	// todo: optimise out abs
#ifdef STATIC_CALCULATION
	return 0;
#endif 
	return div2_floor(abs(W));
    }
    inline static const int minus_nu_0(const int W){
#ifdef STATIC_CALCULATION
	return 0;
#endif
	return div2_floor(-abs(W)-1) - 1*(fold(W, 2));
    }
};


// include implementations 
#include <../src/frg/sbe/rhs_1lfrg.tpp> // for general logic of the right hand side calculation
#include <../src/frg/sbe/rhs_1lfrg_eval.tpp> // the flow equations
#include <../src/frg/sbe/rhs_1lfrg_eval_SBEb.tpp> // the flow equations for the SBEb
#include <../src/frg/sbe/rhs_1lfrg_bosonised_M_eval.tpp> // the flow equations for the bosonised M parts
#include <../src/frg/sbe/rhs_1lfrg_eval_Udot.tpp> 

