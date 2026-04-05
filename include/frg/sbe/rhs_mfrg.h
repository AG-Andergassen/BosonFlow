#pragma once

#include <frg/sbe/rhs_1lfrg.h>

template <typename Model, typename state_t >
class rhs_sbe_mfrg_t : public rhs_sbe_1lfrg_t<Model, state_t>
{
public:

    rhs_sbe_mfrg_t():
    m_lambda_sc_dot_left_corrections(/*FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1*/),
	m_lambda_d_dot_left_corrections(/*FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1*/),
	m_lambda_m_dot_left_corrections(/*FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1*/),
	m_M_sc_dot_left_corrections(/*FrequencyDependenceScheme<Model>::M_BosonicGrid().positive_freqs_count, FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1*/),
	m_M_d_dot_left_corrections(/*FrequencyDependenceScheme<Model>::M_BosonicGrid().positive_freqs_count, FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1*/),
	m_M_m_dot_left_corrections(/*FrequencyDependenceScheme<Model>::M_BosonicGrid().positive_freqs_count, FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1*/){}

    using rhs_base_t = ::rhs_base_t<Model, state_t >;

    void operator() ( const state_t& old_state, state_t &dstate_over_dt, const double t );

  static int &LatestNumberOfSelfEnergyIterations()
  {
      static int last_number_of_selfenergy_iterations;
      return last_number_of_selfenergy_iterations;
  };

  static int &LatestNumberOfVertexCorrections()
  {
      static int last_number_of_vertex_corrections;
      return last_number_of_vertex_corrections;
  };

  static int &LatestTotalNumberOfVertexCorrections()
  {
      static int last_number_of_total_vertex_corrections;
      return last_number_of_total_vertex_corrections;
  };

  static int &MaxNumberOfVertexCorrections()
  {
      static int max_number_of_vertex_corrections;
      return max_number_of_vertex_corrections;
  };
  
  static std::vector<state_t> &LatestStateDotMultiloopCorrectionsOrders()
  {
#ifdef STORE_MULTILOOP_VERTEX_CORRECTIONS    
      static std::vector<state_t> latest_dstate_over_dt_multiloop_corrections_orders(1 + EXTRA_LOOP_NUM, state_t());
#else
      static std::vector<state_t> latest_dstate_over_dt_multiloop_corrections_orders(0, state_t());
#endif
      return latest_dstate_over_dt_multiloop_corrections_orders;
  };

  
protected:

    gf_lambda_t<Model> m_lambda_sc_dot_left_corrections;
    gf_lambda_t<Model> m_lambda_d_dot_left_corrections;
    gf_lambda_t<Model> m_lambda_m_dot_left_corrections;
    gf_M_t<Model> m_M_sc_dot_left_corrections;
    gf_M_t<Model> m_M_d_dot_left_corrections;
    gf_M_t<Model> m_M_m_dot_left_corrections;

    gf_lambda_t<Model> m_lambda_sc_dot_central_corrections;
    gf_lambda_t<Model> m_lambda_d_dot_central_corrections;
    gf_lambda_t<Model> m_lambda_m_dot_central_corrections;
    gf_w_t<Model> m_w_sc_dot_central_corrections;
    gf_w_t<Model> m_w_d_dot_central_corrections;
    gf_w_t<Model> m_w_m_dot_central_corrections;
    gf_M_t<Model> m_M_sc_dot_central_corrections;
    gf_M_t<Model> m_M_d_dot_central_corrections;
    gf_M_t<Model> m_M_m_dot_central_corrections;

  void add_multiloop_vertex_corrections(unsigned extra_loops_count, const double convergence_error, state_t *dstate_over_dt_multiloop_ptr, const state_t &old_state, const state_t &dstate_over_dt_1l, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

  void add_multiloop_self_energy_corrections(double *self_energy_iteration_error_ptr, state_t *dstate_over_dt_ptr, const state_t &old_state, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model>& bubble_GS_pp, const gf_bubble_mat_t<Model>& bubble_GS_ph, const gf_bubble_mat_t<Model>& bubble_GG_pp, const gf_bubble_mat_t<Model>& bubble_GG_ph, const double t);
  
    // ------- 2-loop vertex corrections -----------

    void rhs_frg_2l_corrections(state_t *dstate_over_dt_2l_corrections_ptr, const state_t &old_state, const state_t &dstate_over_dt_1l, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    // lambda (only precalculated left corrections)
    void precalculate_lambda_sc_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);

    void precalculate_lambda_d_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void precalculate_lambda_m_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    // M (only precalculated left corrections)
    void precalculate_M_sc_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);

    void precalculate_M_d_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void precalculate_M_m_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    // lambda (Hedin vertex)
    void rhs_lambda_sc_2l_corrections(gf_lambda_t<Model> *dlambda_sc_over_dt_2l_corrections_ptr, const double t);

    void rhs_lambda_d_2l_corrections(gf_lambda_t<Model> *dlambda_d_over_dt_2l_corrections_ptr, const double t);

    void rhs_lambda_m_2l_corrections(gf_lambda_t<Model> *dlambda_m_over_dt_2l_corrections_ptr, const double t);

    // M (rest function)

    void rhs_M_sc_2l_corrections(gf_M_t<Model> *dM_sc_over_dt_2l_corrections_ptr, const double t);

    void rhs_M_d_2l_corrections(gf_M_t<Model> *dM_d_over_dt_2l_corrections_ptr, const double t);

    void rhs_M_m_2l_corrections(gf_M_t<Model> *dM_m_over_dt_2l_corrections_ptr, const double t);


  
    // ------- Higher loop (3 or more) vertex corrections-----------

    void rhs_frg_higher_loop_corrections(state_t *dstate_over_dt_l_loop_corrections_ptr, const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const state_t &dstate_over_dt_l_minus_2, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    // lambda (only precalculated central corrections)
    void precalculate_lambda_sc_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);

    void precalculate_lambda_d_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void precalculate_lambda_m_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    // w (only precalculated central corrections)
    void precalculate_w_sc_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);

    void precalculate_w_d_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void precalculate_w_m_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    // M (only precalculated central corrections)
    void precalculate_M_sc_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t);

    void precalculate_M_d_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    void precalculate_M_m_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t);

    // w (screened interaction)
    void rhs_w_sc_higher_loop_corrections(gf_w_t<Model> *dw_sc_over_dt_higher_loop_corrections_ptr, const double t);

    void rhs_w_d_higher_loop_corrections(gf_w_t<Model> *dw_d_over_dt_higher_loop_corrections_ptr, const double t);

    void rhs_w_m_higher_loop_corrections(gf_w_t<Model> *dw_m_over_dt_higher_loop_corrections_ptr, const double t);

  
    // lambda (Hedin vertex)
    void rhs_lambda_sc_higher_loop_corrections(gf_lambda_t<Model> *dlambda_sc_over_dt_higher_loop_corrections_ptr, const double t);

    void rhs_lambda_d_higher_loop_corrections(gf_lambda_t<Model> *dlambda_d_over_dt_higher_loop_corrections_ptr, const double t);

    void rhs_lambda_m_higher_loop_corrections(gf_lambda_t<Model> *dlambda_m_over_dt_higher_loop_corrections_ptr, const double t);

    // M (SBE Rest function)
    void rhs_M_sc_higher_loop_corrections(gf_M_t<Model> *dM_sc_over_dt_higher_loop_corrections_ptr, const double t);

    void rhs_M_d_higher_loop_corrections(gf_M_t<Model> *dM_d_over_dt_higher_loop_corrections_ptr, const double t);

    void rhs_M_m_higher_loop_corrections(gf_M_t<Model> *dM_m_over_dt_higher_loop_corrections_ptr, const double t);

    // Self energy
    void rhs_Sig_SDE_SBE_Magnetic(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, const double t);
    void rhs_Sig_SDE_SBE_Density(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, const double t);
    void rhs_Sig_SDE_SBE_Superconducting(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, const double t);

    void rhs_Sig_SDE_Magnetic_Conventional(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model>& bubble_GS_pp, const gf_bubble_mat_t<Model>& bubble_GS_ph, const gf_bubble_mat_t<Model>& bubble_GG_pp, const gf_bubble_mat_t<Model>& bubble_GG_ph, const double t);
  
    void rhs_Sig_SDE_Density_Conventional(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model>& bubble_GS_pp, const gf_bubble_mat_t<Model>& bubble_GS_ph, const gf_bubble_mat_t<Model>& bubble_GG_pp, const gf_bubble_mat_t<Model>& bubble_GG_ph, const double t);
  

    void compute_GG_and_GS_bubbles(const gf_1p_mat_t<Model> &Gvec, const gf_1p_mat_t<Model> &Svec, gf_bubble_mat_t<Model> *bubble_GS_pp_ptr, gf_bubble_mat_t<Model> *bubble_GS_ph_ptr, gf_bubble_mat_t<Model> *bubble_GG_pp_ptr, gf_bubble_mat_t<Model> *bubble_GG_ph_ptr);

  
private:
  
    dcomplex eval_Sig_SDE_SBE_Magnetic( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, const double t );
    dcomplex eval_Sig_SDE_SBE_Density( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, const double t );
    dcomplex eval_Sig_SDE_SBE_Superconducting( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, const double t );

    dcomplex eval_Sig_SDE_Magnetic_Conventional( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model> &bubble_ph, const gf_bubble_mat_t<Model> &bubble_ph_dot,  const double t );
  
    dcomplex eval_Sig_SDE_Density_Conventional( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model> &bubble_ph, const gf_bubble_mat_t<Model> &bubble_ph_dot,  const double t );
  

    // ------- Left corrections for vertex corrections -----------

    // lambda (Hedin vertex)
    dcomplex eval_lambda_sc_left_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t );

    dcomplex eval_lambda_d_left_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    dcomplex eval_lambda_m_left_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    // M (SBE Rest function)
    dcomplex eval_M_sc_left_corrections( const idx_M_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t );

    dcomplex eval_M_d_left_corrections( const idx_M_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    dcomplex eval_M_m_left_corrections( const idx_M_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    // ------- Central corrections for vertex corrections -----------

    // lambda (Hedin vertex)
    dcomplex eval_lambda_sc_central_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t );

    dcomplex eval_lambda_d_central_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    dcomplex eval_lambda_m_central_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    // w (screened interaction)
    dcomplex eval_w_sc_central_corrections( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t );

    dcomplex eval_w_d_central_corrections( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    dcomplex eval_w_m_central_corrections( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    // M (SBE Rest function)
    dcomplex eval_M_sc_central_corrections( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t );

    dcomplex eval_M_d_central_corrections( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );

    dcomplex eval_M_m_central_corrections( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t );


    // ----- 2-loop vertex corrections eval functions (explicit evaluation of the equations on the right handside of the ODE) ----

    // lambda ( Hedin vertex)
    dcomplex eval_lambda_sc_2l_corrections( const idx_lambda_t<Model>& idx, double t );
  
    dcomplex eval_lambda_d_2l_corrections( const idx_lambda_t<Model>& idx, double t );

    dcomplex eval_lambda_m_2l_corrections( const idx_lambda_t<Model>& idx, double t );
  
    // M ( SBE Rest function )
    dcomplex eval_M_sc_2l_corrections( const idx_M_t<Model>& idx, double t );
  
    dcomplex eval_M_d_2l_corrections( const idx_M_t<Model>& idx, double t );

    dcomplex eval_M_m_2l_corrections( const idx_M_t<Model>& idx, double t );
  

    // ------ Higher loop (3 and above) vertex corrections eval functions ----
  
    // lambda (Hedin vertex) 
    dcomplex eval_lambda_sc_higher_loop_corrections( const idx_lambda_t<Model>& idx, double t );
  
    dcomplex eval_lambda_d_higher_loop_corrections( const idx_lambda_t<Model>& idx, double t );

    dcomplex eval_lambda_m_higher_loop_corrections( const idx_lambda_t<Model>& idx, double t );

    // w (screened interaction) 
    dcomplex eval_w_sc_higher_loop_corrections( const idx_w_t<Model>& idx, double t );
  
    dcomplex eval_w_d_higher_loop_corrections( const idx_w_t<Model>& idx, double t );

    dcomplex eval_w_m_higher_loop_corrections( const idx_w_t<Model>& idx, double t );
  
    // M (SBE rest function) 
    dcomplex eval_M_sc_higher_loop_corrections( const idx_M_t<Model>& idx, double t );
  
    dcomplex eval_M_d_higher_loop_corrections( const idx_M_t<Model>& idx, double t );

    dcomplex eval_M_m_higher_loop_corrections( const idx_M_t<Model>& idx, double t );

  
};

// include implementation
#include <../src/frg/sbe/rhs_mfrg.tpp> // implements logic for calculation for the right hand side
#include <../src/frg/sbe/rhs_mfrg_eval.tpp> // the flow equations


