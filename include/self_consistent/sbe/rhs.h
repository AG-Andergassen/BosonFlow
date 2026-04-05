#pragma once

#include <mymath.h>
#include <frg/sbe/state.h>
#include <frg/rhs_base.h>


template <typename Model, typename state_t >
class rhs_sbe_selfconsistent_t : public rhs_base_t<Model, state_t>
{
 public:
   using rhs_base_t = ::rhs_base_t<Model, state_t >;

    virtual void operator() (const state_t &old_state, state_t &dstate_over_dt, const double t) {}; // Overload, but we leave empty. This is not the rhs of a flow equation // (TODO: maybe make it nonpure virtual)

    void vertex( const state_t& old_state, state_t &new_state, const double prev_norm ); 	    /**< Overload call operator to calculate rhs for vertex */

    void selfenergy(const state_t &state, gf_1p_t<Model> &new_Sig); 

    void update_bubbles_and_G(const state_t &state, const double t);
    	
    

    rhs_sbe_selfconsistent_t();
    ~rhs_sbe_selfconsistent_t();

 protected:
    gf_1p_mat_t<Model> m_Gvec;
    gf_bubble_mat_t<Model> m_bubble_pp, m_bubble_ph;

    // screened interaction w rhs
    void rhs_w_sc(gf_w_t<Model> *new_w_sc_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_pp, const double prev_norm);
    void rhs_w_d(gf_w_t<Model> *new_w_d_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph, const double prev_norm);
    void rhs_w_m(gf_w_t<Model> *new_w_m_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph, const double prev_norm);
    
    // Hedin vertex 
    void rhs_lambda_sc(gf_lambda_t<Model> *new_lambda_sc_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_pp);
    void rhs_lambda_d(gf_lambda_t<Model> *new_lambda_d_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph);
    void rhs_lambda_m(gf_lambda_t<Model> *new_lambda_m_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph);

    // rest function
    void rhs_M_sc(gf_M_t<Model> *new_M_sc_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_pp);
    void rhs_M_d(gf_M_t<Model> *new_M_d_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph);
    void rhs_M_m(gf_M_t<Model> *new_M_m_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph);

    // self energy
    void rhs_Sig(gf_1p_t<Model> *new_Sig_ptr, const state_t &old_state, const gf_bubble_mat_t<Model>& bubble_ph, const gf_1p_mat_t<Model>& Gvec);

 private:

    static dcomplex eval_Sig( const idx_1p_t<Model>& idx, const state_t& state_vec, const gf_1p_mat_t<Model>& Gvec);

    static dcomplex eval_Sig_conventional_SDE_magnetic( const idx_1p_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, const gf_1p_mat_t<Model>& Gvec);
    
    static dcomplex eval_Sig_sde_sbe_magnetic( const idx_1p_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec);
    static dcomplex eval_Sig_sde_sbe_superconducting( const idx_1p_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec);
    static dcomplex eval_Sig_sde_sbe_density( const idx_1p_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec);
  
    static dcomplex eval_w_sc( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp );	/**< Calculates the screened interaction in the U-pp-reducible channel w_pp at a given index */
    static dcomplex eval_w_d( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph );	/**< Calculates the screened interaction in the U-ph-reducible channel w_d at a given index */
    static dcomplex eval_w_m( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph );	/**< Calculates the screened interaction in the U-xph-reducible channel w_m at a given index */

    static dcomplex eval_w_sc_BSE( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, const gf_w_t<Model> &w_sc );	/**< Calculates the screened interaction BSE right-hand-side in the U-pp-reducible channel w_pp at a given index */
    static dcomplex eval_w_d_BSE( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, const gf_w_t<Model> &w_d );	/**< Calculates the screened interaction BSE right-hand-side in the U-ph-reducible channel w_d at a given index */
    static dcomplex eval_w_m_BSE( const idx_w_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, const gf_w_t<Model> &w_m );	/**< Calculates the screened interaction BSE right-hand-side in the U-xph-reducible channel w_m at a given index */

  
    static dcomplex eval_lambda_sc( const idx_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp ); /**< Calculates the Hedin vertex in the U-pp-reducible channel lambda_sc at a given index */
    static dcomplex eval_lambda_d( const idx_lambda_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_ph );	/**< Calculates the Hedin vertex in the U-ph-reducible channel lambda_d at a given index */
    static dcomplex eval_lambda_m( const idx_lambda_t<Model>& idx, const state_t&  state_vec, const gf_bubble_mat_t<Model>& bubble_ph );	/**< Calculates the Hedin vertex in the U-xph-reducible channel lambda_m at a given index */


    static dcomplex eval_M_sc( const idx_M_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp );  /**< Calculates the U-irreducible rest function M_sc at a given index */
    static dcomplex eval_M_d( const idx_M_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph );   /**< Calculates the U-irreducible rest function M_d at a given index */
    static dcomplex eval_M_m( const idx_M_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph );   /**< Calculates the U-irreducible rest function M_m at a given index */
};


// include implementation
#include "../src/self_consistent/sbe/rhs.tpp" // logic of the calculation of the right hand side
#include "../src/self_consistent/sbe/rhs_eval.tpp" // the self-consistent SBE equations
