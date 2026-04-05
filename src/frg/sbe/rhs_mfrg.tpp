#include <frg/sbe/symmetries.h>
#include <frg/sbe/symmetries_mfrg.h>
#include <frg/flows.h>
#include <iostream>

#include <params_technical.h>

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::operator() ( const state_t& old_state, state_t& dstate_over_dt, const double t ){
    rhs_base_t::print_current_scale(t);

#ifdef DEBUG_BUILD
    omp_set_num_threads(1);  	
#endif

    gf_1p_mat_t<Model> Gvec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true ); 
    gf_1p_mat_t<Model> Svec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true );
    gf_1p_mat_t<Model> Gdotvec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true );

    std::cout << " ... computing G and S" << std::endl;
    rhs_base_t::compute_Gvec(&Gvec, old_state, t);
    rhs_base_t::compute_Svec(&Svec, old_state, t);

    std::cout << " ... computing 1-loop self-energy derivative" << std::endl;

    #pragma omp parallel default(shared)
    {

    #ifdef SELFEN_FLOW

    #pragma omp single
	std::cout << " ... self-energy (standard self-energy flow equation)" << std::endl;

        rhs_sbe_1lfrg_t<Model, state_t>::rhs_Sig_1l(&dstate_over_dt.gf_Sig(), old_state, Svec, t);

    #else 
	dstate_over_dt.gf_Sig().init( [](const idx_1p_t<Model> &idx){return 0;} );
    #endif
    }

    gf_bubble_mat_t<Model> bubble_GS_pp, bubble_GS_ph;
    gf_bubble_mat_t<Model> bubble_GG_pp, bubble_GG_ph; // needed for multiloop

    state_t dstate_over_dt_1l;
    dstate_over_dt_1l.init_zero();

    LatestTotalNumberOfVertexCorrections() = 0;
    MaxNumberOfVertexCorrections() = 1;

    double self_energy_iteration_error = 1e-1;
    for (unsigned se_loop_i = 0; se_loop_i < SELFENERGY_ITERATIONS_MAX; se_loop_i++){
	std::cout << " ... computing self-energy iteration: " << se_loop_i << std::endl;

#ifndef NO_KATANIN
	std::cout << " ... adding the katanin correction to S" << std::endl;
	rhs_sbe_1lfrg_t<Model, state_t>::add_katanin_correction(&Gdotvec, Gvec, Svec, dstate_over_dt);
#else
	std::cout << " ... NO katanin correction added!!" << std::endl;
	Gdotvec = Svec;
#endif
		
	std::cout << " ... computing GS and GG bubbles" << std::endl;
	compute_GG_and_GS_bubbles(Gvec, Gdotvec, &bubble_GS_pp, &bubble_GS_ph, &bubble_GG_pp, &bubble_GG_ph);

	std::cout << " ... computing 1-loop fRG vertex derivative " << std::endl;
	rhs_sbe_1lfrg_t<Model, state_t>::rhs_frg_1l(&dstate_over_dt_1l, old_state, bubble_GS_pp, bubble_GS_ph, bubble_GG_pp, bubble_GG_ph, Svec, t);
					 
	LatestTotalNumberOfVertexCorrections() += 1;
	
	dstate_over_dt_1l.gf_Sig() = dstate_over_dt.gf_Sig();

	#ifdef STORE_MULTILOOP_VERTEX_CORRECTIONS
	    LatestStateDotMultiloopCorrectionsOrders()[0] = dstate_over_dt_1l;
	#endif
	
	std::cout << " ... adding multiloop vertex corrections " << std::endl;
#ifdef PRECOMPUTE_STATE_PROJECTIONS
	if (EXTRA_LOOP_NUM > 0)
	    dstate_over_dt_1l.precalculate_projected_Ms_and_nablas_dot(old_state, t);
#endif
	dstate_over_dt = dstate_over_dt_1l;

	//add_multiloop_vertex_corrections(EXTRA_LOOP_NUM, ABS_ERROR_VERTEX_MULTILOOP_CORRECTIONS, &dstate_over_dt, old_state, dstate_over_dt_1l, bubble_GG_pp, bubble_GG_ph, t);
	add_multiloop_vertex_corrections(EXTRA_LOOP_NUM, std::max(ABS_ERROR_VERTEX_MULTILOOP_CORRECTIONS, self_energy_iteration_error/10.0), &dstate_over_dt, old_state, dstate_over_dt_1l, bubble_GG_pp, bubble_GG_ph, t);

	std::cout << " ... adding self-energy corrections " << std::endl;
	add_multiloop_self_energy_corrections(&self_energy_iteration_error, &dstate_over_dt, old_state, Gvec, Gdotvec, bubble_GS_pp, bubble_GS_ph, bubble_GG_pp, bubble_GG_ph, t);
	std::cout << "Self-energy iteration error: " << self_energy_iteration_error << std::endl;

#ifndef FORCE_CALCULATE_ALL_SELFENERGY_ITERATIONS
	if (self_energy_iteration_error < ABS_ERROR_SELFENERGY_ITERATIONS){
	  std::cout << "Self-energy iterations converge at iteration number " << se_loop_i << std::endl;
	  LatestNumberOfSelfEnergyIterations() = se_loop_i;
	  break;
	}
#endif

    }
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::compute_GG_and_GS_bubbles(const gf_1p_mat_t<Model> &Gvec, const gf_1p_mat_t<Model> &Gdotvec, gf_bubble_mat_t<Model> *bubble_GS_pp_ptr, gf_bubble_mat_t<Model> *bubble_GS_ph_ptr, gf_bubble_mat_t<Model> *bubble_GG_pp_ptr, gf_bubble_mat_t<Model> *bubble_GG_ph_ptr)
{
    auto &bubble_GS_pp = *bubble_GS_pp_ptr;
    auto &bubble_GS_ph = *bubble_GS_ph_ptr;
    auto &bubble_GG_pp = *bubble_GG_pp_ptr;
    auto &bubble_GG_ph = *bubble_GG_ph_ptr;

    gf_1p_mat_real_t<Model> Gvec_realspace(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true), 
	Gdotvec_realspace(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true);   

#pragma omp parallel default(shared) 
    {
	rhs_base_t::compute_Gvec_real(&Gvec_realspace, Gvec);
	rhs_base_t::compute_Svec_real(&Gdotvec_realspace, Gdotvec);	
	
#pragma omp barrier //barrier to make sure G/Gdotvec_real are present
	{
	    rhs_base_t::compute_GS_bubbles(&bubble_GS_pp, &bubble_GS_ph, Gvec_realspace, Gdotvec_realspace);
	
	    // TODO: see if calculation can be optimised.. Also there will be a factor of 2 (fix later).
	    rhs_base_t::compute_GS_bubbles(&bubble_GG_pp, &bubble_GG_ph, Gvec_realspace, Gvec_realspace);
	}
    }

    bubble_GG_pp.init([&bubble_GG_pp](const idx_bubble_mat_t<Model> &idx){return 0.5 * bubble_GG_pp(idx);});
    bubble_GG_ph.init([&bubble_GG_ph](const idx_bubble_mat_t<Model> &idx){return 0.5 * bubble_GG_ph(idx);});
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::add_multiloop_self_energy_corrections(double *self_energy_iteration_error_ptr, state_t *dstate_over_dt_ptr, const state_t &old_state, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model>& bubble_GS_pp, const gf_bubble_mat_t<Model>& bubble_GS_ph, const gf_bubble_mat_t<Model>& bubble_GG_pp, const gf_bubble_mat_t<Model>& bubble_GG_ph,  const double t)
{
    double &self_energy_iteration_error = *self_energy_iteration_error_ptr;
    state_t &dstate_over_dt = *dstate_over_dt_ptr;

    state_t old_dstate_over_dt = dstate_over_dt;


#pragma omp parallel default(shared) 
    {
      
#ifdef SELFEN_SDE_SBE_MAGNETIC

    #pragma omp single
    std::cout << " ... self-energy (Schwinger-Dyson based flow in the SBE formalism (using the magnetic channel)) " << std::endl;

    rhs_Sig_SDE_SBE_Magnetic(&dstate_over_dt.gf_Sig(), old_state, dstate_over_dt, Gvec, Gdotvec, t);
#endif

#ifdef SELFEN_SDE_SBE_DENSITY

    #pragma omp single
    std::cout << " ... self-energy (Schwinger-Dyson based flow in the SBE formalism (using the density channel)) " << std::endl;

    rhs_Sig_SDE_SBE_Density(&dstate_over_dt.gf_Sig(), old_state, dstate_over_dt, Gvec, Gdotvec, t);
#endif

#ifdef SELFEN_SDE_SBE_SUPERCONDUCTING
    
    #pragma omp single
    std::cout << " ... self-energy (Schwinger-Dyson based flow in the SBE formalism (using the superconducting channel)) " << std::endl;

    rhs_Sig_SDE_SBE_Superconducting(&dstate_over_dt.gf_Sig(), old_state, dstate_over_dt, Gvec, Gdotvec, t);
#endif

#ifdef SELFEN_SDE_CONVENTIONAL_MAGNETIC
    #pragma omp single
    std::cout << " ... self-energy (conventional magnetic Schwinger-Dyson based flow) " << std::endl;
    rhs_Sig_SDE_Magnetic_Conventional(&dstate_over_dt.gf_Sig(), old_state, dstate_over_dt, Gvec, Gdotvec, bubble_GS_pp, bubble_GS_ph, bubble_GG_pp, bubble_GG_ph, t);
#endif

#ifdef SELFEN_SDE_CONVENTIONAL_DENSITY
    #pragma omp single
    std::cout << " ... self-energy (conventional density Schwinger-Dyson based flow) " << std::endl;
    rhs_Sig_SDE_Density_Conventional(&dstate_over_dt.gf_Sig(), old_state, dstate_over_dt, Gvec, Gdotvec, bubble_GS_pp, bubble_GS_ph, bubble_GG_pp, bubble_GG_ph, t);
#endif


    }

    
    self_energy_iteration_error = norm(dstate_over_dt - old_dstate_over_dt);
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::add_multiloop_vertex_corrections(unsigned extra_loops_count, const double convergence_error, state_t *dstate_over_dt_multiloop_ptr, const state_t &old_state, const state_t &dstate_over_dt_1l, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
  state_t &dstate_over_dt_multiloop = *dstate_over_dt_multiloop_ptr;
  
  dstate_over_dt_multiloop = dstate_over_dt_1l;
  
  
  if (extra_loops_count == 0){
    std::cout << " ... no multiloop vertex corrections were performed " << std::endl;
    return;
  }
    
    state_t dstate_over_dt_2l_corrections;
    dstate_over_dt_2l_corrections.m_force_w_dot_asymptotics_to_zero = true;
    dstate_over_dt_2l_corrections.init_zero();

    std::cout << " ... computing 2-loop corrections " << std::endl;
    rhs_frg_2l_corrections(&dstate_over_dt_2l_corrections, old_state, dstate_over_dt_1l, bubble_GG_pp, bubble_GG_ph, t);

    LatestTotalNumberOfVertexCorrections() += 1;
    MaxNumberOfVertexCorrections() = std::max(2, MaxNumberOfVertexCorrections());

    #ifdef STORE_MULTILOOP_VERTEX_CORRECTIONS
        LatestStateDotMultiloopCorrectionsOrders()[1] = dstate_over_dt_2l_corrections;
    #endif
    
    dstate_over_dt_multiloop += dstate_over_dt_2l_corrections;
    
    auto loop_correction_error = norm(dstate_over_dt_2l_corrections);
	
    /*#ifndef FORCE_CALCULATE_ALL_MULTILOOP_CORRECTIONS
        auto loop_correction_error = norm(dstate_over_dt_2l_corrections);
	std::cout << "Multiloop vertex correction error at loop " << 2 << " is : " << loop_correction_error << std::endl;
	if (loop_correction_error < convergence_error){
	  std::cout << "Multiloop corrections converge at loop number " << 2 << std::endl; 
	  return;
	}
    #endif
    */

    if (extra_loops_count == 1)
	   return; // only two loop
	
    #ifdef PRECOMPUTE_STATE_PROJECTIONS
        dstate_over_dt_2l_corrections.precalculate_projected_Ms_and_nablas_dot(old_state, t);
    #endif

    state_t dstate_over_dt_l_minus_2 = dstate_over_dt_1l;
    state_t dstate_over_dt_l_minus_1 = dstate_over_dt_2l_corrections;

    dstate_over_dt_l_minus_1.m_force_w_dot_asymptotics_to_zero = true;
  
    for (unsigned l = 3; l < 2+extra_loops_count; l++){
	state_t dstate_over_dt_l_loop_corrections;
	dstate_over_dt_l_loop_corrections.m_force_w_dot_asymptotics_to_zero = true;
	dstate_over_dt_l_loop_corrections.init_zero();
	
	std::cout << " ... computing " << l << "-loop vertex corrections " << std::endl;
	rhs_frg_higher_loop_corrections(&dstate_over_dt_l_loop_corrections, old_state, dstate_over_dt_l_minus_1, dstate_over_dt_l_minus_2, bubble_GG_pp, bubble_GG_ph, t);

	LatestTotalNumberOfVertexCorrections() += 1; // this is "one"
	LatestNumberOfVertexCorrections() = l; // this is "ell"
	MaxNumberOfVertexCorrections()  = std::max(MaxNumberOfVertexCorrections(), int(l));

	// TODO: figure out maybe how to make addition of state more selective, as the self-energies here for example play no role. Idea: implement "add_vertex" method to state, which only adds the vertex components
	dstate_over_dt_multiloop += dstate_over_dt_l_loop_corrections;

        #ifdef STORE_MULTILOOP_VERTEX_CORRECTIONS
	    LatestStateDotMultiloopCorrectionsOrders()[l-1] = dstate_over_dt_l_loop_corrections;
	#endif

#ifndef FORCE_CALCULATE_ALL_MULTILOOP_CORRECTIONS
	loop_correction_error = norm(dstate_over_dt_l_loop_corrections);
	std::cout << "Multiloop vertex correction error at loop " << l << " is : " << loop_correction_error << std::endl;
	
	if (loop_correction_error < convergence_error){
	    std::cout << "Multiloop corrections converge at loop number " << l << std::endl; 
	    break;
	}
#endif


	if (l+1 < 2 + extra_loops_count){
          #ifdef PRECOMPUTE_STATE_PROJECTIONS	
	      dstate_over_dt_l_loop_corrections.precalculate_projected_Ms_and_nablas_dot(old_state, t);
          #endif
	
	  dstate_over_dt_l_minus_2 = dstate_over_dt_l_minus_1;
	  dstate_over_dt_l_minus_1 = dstate_over_dt_l_loop_corrections;
	  dstate_over_dt_l_minus_1.m_force_w_dot_asymptotics_to_zero = true;
	}
	
    }
}



  // ------- 2-loop vertex corrections -----------

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_frg_2l_corrections(state_t *dstate_over_dt_2l_corrections_ptr, const state_t &old_state, const state_t &dstate_over_dt_1l, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    state_t &dstate_over_dt_2l_corrections = *dstate_over_dt_2l_corrections_ptr;

#pragma omp parallel default(shared)
    {
	// remark 1: there are no left corrections for w
	// remark 2: the left corrections are used for the 2 loop corrections at this stage

	precalculate_lambda_sc_left_corrections(old_state, dstate_over_dt_1l, bubble_GG_pp, t);
	precalculate_lambda_d_left_corrections(old_state, dstate_over_dt_1l, bubble_GG_ph, t);
	precalculate_lambda_m_left_corrections(old_state, dstate_over_dt_1l, bubble_GG_ph, t);
#if !defined(SBEa_APPROXIMATION)

	precalculate_M_sc_left_corrections(old_state, dstate_over_dt_1l, bubble_GG_pp, t);
	precalculate_M_d_left_corrections(old_state, dstate_over_dt_1l, bubble_GG_ph, t);
	precalculate_M_m_left_corrections(old_state, dstate_over_dt_1l, bubble_GG_ph, t);
#endif //!defined(SBEa_APPROXIMATION)
//    }
#pragma omp barrier
//#pragma omp parallel default(shared)
//    {
	// remark: there are no 2-loop corrections for w

#pragma omp single
	std::cout << " ... lambda (Hedin vertex) 2l corrections" << std::endl;

	rhs_lambda_sc_2l_corrections(&dstate_over_dt_2l_corrections.gf_lambda_sc(), t);
	rhs_lambda_d_2l_corrections(&dstate_over_dt_2l_corrections.gf_lambda_d(), t);
	rhs_lambda_m_2l_corrections(&dstate_over_dt_2l_corrections.gf_lambda_m(), t);

#if !defined(SBEa_APPROXIMATION)

#pragma omp single
	std::cout << " ... M (Rest function) 2l corrections" << std::endl;

	rhs_M_sc_2l_corrections(&dstate_over_dt_2l_corrections.gf_M_sc(), t);
	rhs_M_d_2l_corrections(&dstate_over_dt_2l_corrections.gf_M_d(), t);
	rhs_M_m_2l_corrections(&dstate_over_dt_2l_corrections.gf_M_m(), t);
#endif //!defined(SBEa_APPROXIMATION)
    }

}


// lambda (only precalculated left corrections)
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_lambda_sc_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_lambda_sc_LeftCorrections_Ptr()->init_batched( m_lambda_sc_dot_left_corrections, [this, &old_state, &dstate_over_dt_l_minus_1, &bubble_GG_pp, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_sc_left_corrections( idx, old_state, dstate_over_dt_l_minus_1, bubble_GG_pp, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_lambda_d_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_lambda_d_LeftCorrections_Ptr()->init_batched( m_lambda_d_dot_left_corrections, [this, &old_state, &dstate_over_dt_l_minus_1, &bubble_GG_ph, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_d_left_corrections( idx, old_state, dstate_over_dt_l_minus_1, bubble_GG_ph, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_lambda_m_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_lambda_m_LeftCorrections_Ptr()->init_batched( m_lambda_m_dot_left_corrections, [this, &old_state, &dstate_over_dt_l_minus_1, &bubble_GG_ph, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_m_left_corrections( idx, old_state, dstate_over_dt_l_minus_1, bubble_GG_ph, t ); } ); 
}

// M (only precalculated left corrections)
#if !defined(SBEa_APPROXIMATION)
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_M_sc_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_M_sc_LeftCorrections_Ptr()->init_batched( m_M_sc_dot_left_corrections, [this, &old_state, &dstate_over_dt_l_minus_1, &bubble_GG_pp, t]( const idx_M_t<Model>& idx ){ return eval_M_sc_left_corrections( idx, old_state, dstate_over_dt_l_minus_1, bubble_GG_pp, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_M_d_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_M_d_LeftCorrections_Ptr()->init_batched( m_M_d_dot_left_corrections, [this, &old_state, &dstate_over_dt_l_minus_1, &bubble_GG_ph, t]( const idx_M_t<Model>& idx ){ return eval_M_d_left_corrections( idx, old_state, dstate_over_dt_l_minus_1, bubble_GG_ph, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_M_m_left_corrections(const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_M_m_LeftCorrections_Ptr()->init_batched( m_M_m_dot_left_corrections, [this, &old_state, &dstate_over_dt_l_minus_1, &bubble_GG_ph, t]( const idx_M_t<Model>& idx ){ return eval_M_m_left_corrections( idx, old_state, dstate_over_dt_l_minus_1, bubble_GG_ph, t ); } ); 
}
#endif //!defined(SBEa_APPROXIMATION)

// lambda
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_lambda_sc_2l_corrections(gf_lambda_t<Model> *dlambda_sc_over_dt_2l_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_sc_Ptr()->init_batched( (*dlambda_sc_over_dt_2l_corrections_ptr), [this, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_sc_2l_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_lambda_d_2l_corrections(gf_lambda_t<Model> *dlambda_d_over_dt_2l_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_d_Ptr()->init_batched( (*dlambda_d_over_dt_2l_corrections_ptr), [this, t]( const idx_lambda_t<Model>& idx ){ return eval_lambda_d_2l_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_lambda_m_2l_corrections(gf_lambda_t<Model> *dlambda_m_over_dt_2l_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_m_Ptr()->init_batched( (*dlambda_m_over_dt_2l_corrections_ptr), [this, t]( const idx_lambda_t<Model>& idx ){ 
return eval_lambda_m_2l_corrections( idx, t ); 

} ); 
}

// M
#if !defined(SBEa_APPROXIMATION)
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_M_sc_2l_corrections(gf_M_t<Model> *dM_sc_over_dt_2l_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched( (*dM_sc_over_dt_2l_corrections_ptr), [this, t]( const idx_M_t<Model>& idx ){ return eval_M_sc_2l_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_M_d_2l_corrections(gf_M_t<Model> *dM_d_over_dt_2l_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched( (*dM_d_over_dt_2l_corrections_ptr), [this, t]( const idx_M_t<Model>& idx ){ return eval_M_d_2l_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_M_m_2l_corrections(gf_M_t<Model> *dM_m_over_dt_2l_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched( (*dM_m_over_dt_2l_corrections_ptr), [this, t]( const idx_M_t<Model>& idx ){ return eval_M_m_2l_corrections( idx, t ); } ); 
}
#endif //!defined(SBEa_APPROXIMATION)



  // ------- Higher loop (3 or more) vertex corrections-----------

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_frg_higher_loop_corrections(state_t *dstate_over_dt_l_loop_corrections_ptr, const state_t &old_state, const state_t &dstate_over_dt_l_minus_1, const state_t &dstate_over_dt_l_minus_2, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    state_t &dstate_over_dt_l_loop_corrections = *dstate_over_dt_l_loop_corrections_ptr;

#pragma omp parallel default(shared)
    {
	precalculate_lambda_sc_central_corrections(old_state, bubble_GG_pp, t);
	precalculate_lambda_d_central_corrections(old_state, bubble_GG_ph, t);
	precalculate_lambda_m_central_corrections(old_state, bubble_GG_ph, t);

	precalculate_w_sc_central_corrections(old_state, bubble_GG_pp, t);
	precalculate_w_d_central_corrections(old_state, bubble_GG_ph, t);
	precalculate_w_m_central_corrections(old_state, bubble_GG_ph, t);

#if !defined(SBEa_APPROXIMATION)
	precalculate_M_sc_central_corrections(old_state, bubble_GG_pp, t);
	precalculate_M_d_central_corrections(old_state, bubble_GG_ph, t);
	precalculate_M_m_central_corrections(old_state, bubble_GG_ph, t);
#endif //!defined(SBEa_APPROXIMATION)
    }

#pragma omp parallel default(shared)
    {
	precalculate_lambda_sc_left_corrections(old_state, dstate_over_dt_l_minus_1, bubble_GG_pp, t);
	precalculate_lambda_d_left_corrections(old_state, dstate_over_dt_l_minus_1, bubble_GG_ph, t);
	precalculate_lambda_m_left_corrections(old_state, dstate_over_dt_l_minus_1, bubble_GG_ph, t);

#if !defined(SBEa_APPROXIMATION)
	precalculate_M_sc_left_corrections(old_state, dstate_over_dt_l_minus_1, bubble_GG_pp, t);
	precalculate_M_d_left_corrections(old_state, dstate_over_dt_l_minus_1, bubble_GG_ph, t);
	precalculate_M_m_left_corrections(old_state, dstate_over_dt_l_minus_1, bubble_GG_ph, t);
#endif //!defined(SBEa_APPROXIMATION)
    }


#pragma omp parallel default(shared)
    {
#pragma omp single
	std::cout << " ... lambda (Hedin vertex) higher loop corrections" << std::endl;

	rhs_lambda_sc_higher_loop_corrections(&dstate_over_dt_l_loop_corrections.gf_lambda_sc(), t);
	rhs_lambda_d_higher_loop_corrections(&dstate_over_dt_l_loop_corrections.gf_lambda_d(), t);
	rhs_lambda_m_higher_loop_corrections(&dstate_over_dt_l_loop_corrections.gf_lambda_m(), t);

#pragma omp single
        std::cout << " ... w (Screened interaction) higher loop corrections" << std::endl;

	rhs_w_sc_higher_loop_corrections(&dstate_over_dt_l_loop_corrections.gf_w_sc(), t);
	rhs_w_d_higher_loop_corrections(&dstate_over_dt_l_loop_corrections.gf_w_d(), t);
	rhs_w_m_higher_loop_corrections(&dstate_over_dt_l_loop_corrections.gf_w_m(), t);

#if !defined(SBEa_APPROXIMATION)

#pragma omp single
	std::cout << " ... M (Rest function) higher loop corrections" << std::endl;

	rhs_M_sc_higher_loop_corrections(&dstate_over_dt_l_loop_corrections.gf_M_sc(), t);
	rhs_M_d_higher_loop_corrections(&dstate_over_dt_l_loop_corrections.gf_M_d(), t);
	rhs_M_m_higher_loop_corrections(&dstate_over_dt_l_loop_corrections.gf_M_m(), t);
#endif //!defined(SBEa_APPROXIMATION)
    }

}


// lambda (only precalculated central corrections)
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_lambda_sc_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_sc_Ptr()->init_batched( m_lambda_sc_dot_central_corrections, [this, &old_state, &bubble_GG_pp, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_sc_central_corrections( idx, old_state, bubble_GG_pp, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_lambda_d_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_d_Ptr()->init_batched( m_lambda_d_dot_central_corrections, [this, &old_state, &bubble_GG_ph, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_d_central_corrections( idx, old_state, bubble_GG_ph, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_lambda_m_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_m_Ptr()->init_batched( m_lambda_m_dot_central_corrections, [this, &old_state, &bubble_GG_ph, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_m_central_corrections( idx, old_state, bubble_GG_ph, t ); } ); 
}

// w (only precalculated central corrections)
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_w_sc_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_sc_Ptr()->init_batched( m_w_sc_dot_central_corrections, [this, &old_state, &bubble_GG_pp, t]( const idx_w_t<Model>& idx ){ return eval_w_sc_central_corrections( idx, old_state, bubble_GG_pp, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_w_d_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_d_Ptr()->init_batched( m_w_d_dot_central_corrections, [this, &old_state, &bubble_GG_ph, t]( const idx_w_t<Model>& idx ){ return eval_w_d_central_corrections( idx, old_state, bubble_GG_ph, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_w_m_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_m_Ptr()->init_batched( m_w_m_dot_central_corrections, [this, &old_state, &bubble_GG_ph, t]( const idx_w_t<Model>& idx ){ return eval_w_m_central_corrections( idx, old_state, bubble_GG_ph, t ); } ); 
}

// M (only precalculated central corrections)
#if !defined(SBEa_APPROXIMATION)
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_M_sc_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched( m_M_sc_dot_central_corrections, [this, &old_state, &bubble_GG_pp, t]( const idx_M_t<Model>& idx ){ return eval_M_sc_central_corrections( idx, old_state, bubble_GG_pp, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_M_d_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched( m_M_d_dot_central_corrections, [this, &old_state, &bubble_GG_ph, t]( const idx_M_t<Model>& idx ){ return eval_M_d_central_corrections( idx, old_state, bubble_GG_ph, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::precalculate_M_m_central_corrections(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched( m_M_m_dot_central_corrections, [this, &old_state, &bubble_GG_ph, t]( const idx_M_t<Model>& idx ){ return eval_M_m_central_corrections( idx, old_state, bubble_GG_ph, t ); } ); 
}
#endif //!defined(SBEa_APPROXIMATION)

// lambda
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_lambda_sc_higher_loop_corrections(gf_lambda_t<Model> *dlambda_sc_over_dt_higher_loop_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_sc_Ptr()->init_batched( (*dlambda_sc_over_dt_higher_loop_corrections_ptr), [this, t]( const idx_lambda_t<Model>& idx ){ return eval_lambda_sc_higher_loop_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_lambda_d_higher_loop_corrections(gf_lambda_t<Model> *dlambda_d_over_dt_higher_loop_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_d_Ptr()->init_batched( (*dlambda_d_over_dt_higher_loop_corrections_ptr), [this, t]( const idx_lambda_t<Model>& idx ){ return eval_lambda_d_higher_loop_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_lambda_m_higher_loop_corrections(gf_lambda_t<Model> *dlambda_m_over_dt_higher_loop_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_m_Ptr()->init_batched( (*dlambda_m_over_dt_higher_loop_corrections_ptr), [this, t]( const idx_lambda_t<Model>& idx ){ 
return eval_lambda_m_higher_loop_corrections( idx, t ); } ); 
}

// w
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_w_sc_higher_loop_corrections(gf_w_t<Model> *dw_sc_over_dt_higher_loop_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_sc_Ptr()->init_batched( (*dw_sc_over_dt_higher_loop_corrections_ptr), [this, t]( const idx_w_t<Model>& idx ){ return eval_w_sc_higher_loop_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_w_d_higher_loop_corrections(gf_w_t<Model> *dw_d_over_dt_higher_loop_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_d_Ptr()->init_batched( (*dw_d_over_dt_higher_loop_corrections_ptr), [this, t]( const idx_w_t<Model>& idx ){ return eval_w_d_higher_loop_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_w_m_higher_loop_corrections(gf_w_t<Model> *dw_m_over_dt_higher_loop_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_m_Ptr()->init_batched( (*dw_m_over_dt_higher_loop_corrections_ptr), [this, t]( const idx_w_t<Model>& idx ){ 
return eval_w_m_higher_loop_corrections( idx, t ); } ); 
}

// M
#if !defined(SBEa_APPROXIMATION)
template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_M_sc_higher_loop_corrections(gf_M_t<Model> *dM_sc_over_dt_higher_loop_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched( (*dM_sc_over_dt_higher_loop_corrections_ptr), [this, t]( const idx_M_t<Model>& idx ){ return eval_M_sc_higher_loop_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_M_d_higher_loop_corrections(gf_M_t<Model> *dM_d_over_dt_higher_loop_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched( (*dM_d_over_dt_higher_loop_corrections_ptr), [this, t]( const idx_M_t<Model>& idx ){ return eval_M_d_higher_loop_corrections( idx, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_M_m_higher_loop_corrections(gf_M_t<Model> *dM_m_over_dt_higher_loop_corrections_ptr, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched( (*dM_m_over_dt_higher_loop_corrections_ptr), [this, t]( const idx_M_t<Model>& idx ){ return eval_M_m_higher_loop_corrections( idx, t ); } ); 
}
#endif //!defined(SBEa_APPROXIMATION)

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_Sig_SDE_Magnetic_Conventional(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model>& bubble_GS_pp, const gf_bubble_mat_t<Model>& bubble_GS_ph, const gf_bubble_mat_t<Model>& bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
  Symmetries1lfRGSBE<Model>::IdxEquivClasses_sig_Ptr()->init_batched( (*dSig_over_dt_ptr), [this, &old_state, &dstate_over_dt, &Gvec, &Gdotvec, &bubble_GG_ph, &bubble_GS_ph, t]( const idx_1p_t<Model>& idx ){ return eval_Sig_SDE_Magnetic_Conventional( idx, old_state, dstate_over_dt, Gvec, Gdotvec, bubble_GG_ph, bubble_GS_ph, t ); } );
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_Sig_SDE_Density_Conventional(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model>& bubble_GS_pp, const gf_bubble_mat_t<Model>& bubble_GS_ph, const gf_bubble_mat_t<Model>& bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
  Symmetries1lfRGSBE<Model>::IdxEquivClasses_sig_Ptr()->init_batched( (*dSig_over_dt_ptr), [this, &old_state, &dstate_over_dt, &Gvec, &Gdotvec, &bubble_GG_ph, &bubble_GS_ph, t]( const idx_1p_t<Model>& idx ){ return eval_Sig_SDE_Density_Conventional( idx, old_state, dstate_over_dt, Gvec, Gdotvec, bubble_GG_ph, bubble_GS_ph, t ); } );
}



template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_Sig_SDE_SBE_Magnetic(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const double t)
{
    // todo: investigate why (parquet originaly) didn't use init_batched (otherwise, if there's no reason, it better be used)
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_sig_Ptr()->init_batched( (*dSig_over_dt_ptr), [this, &old_state, &dstate_over_dt, &Gvec, &Gdotvec, t]( const idx_1p_t<Model>& idx ){ return eval_Sig_SDE_SBE_Magnetic( idx, old_state, dstate_over_dt, Gvec, Gdotvec, t ); } );
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_Sig_SDE_SBE_Density(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const double t)
{
    // todo: investigate why (parquet originaly) didn't use init_batched (otherwise, if there's no reason, it better be used)
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_sig_Ptr()->init_batched( (*dSig_over_dt_ptr), [this, &old_state, &dstate_over_dt, &Gvec, &Gdotvec, t]( const idx_1p_t<Model>& idx ){ return eval_Sig_SDE_SBE_Density( idx, old_state, dstate_over_dt, Gvec, Gdotvec, t ); } );
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_t<Model, state_t>::rhs_Sig_SDE_SBE_Superconducting(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const state_t &dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const double t)
{
    // todo: investigate why (parquet originaly) didn't use init_batched (otherwise, if there's no reason, it better be used)
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_sig_Ptr()->init_batched( (*dSig_over_dt_ptr), [this, &old_state, &dstate_over_dt, &Gvec, &Gdotvec, t]( const idx_1p_t<Model>& idx ){ return eval_Sig_SDE_SBE_Superconducting( idx, old_state, dstate_over_dt, Gvec, Gdotvec, t ); } );
}
