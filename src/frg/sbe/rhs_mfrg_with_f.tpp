#include <frg/sbe/symmetries.h>
#include <frg/sbe/symmetries_mfrg.h>
#include <frg/flows.h>
#include <iostream>

#include <params_technical.h>

template <typename Model, typename state_t >
void rhs_sbe_mfrg_with_F_t<Model, state_t>::operator() ( const state_t& old_state, state_t& dstate_over_dt, const double t ){
    std::cout << " Preprocessing at scale " << t << std::endl;

#ifdef DEBUG_BUILD
    omp_set_num_threads(1);  	
#endif

    
    gf_1p_mat_t<Model> Gvec( FrequenciesCount::POS_1P_RANGE, Model::GetFineMomentaCount(), true ); 
    gf_1p_mat_t<Model> Svec( FrequenciesCount::POS_1P_RANGE, Model::GetFineMomentaCount(), true );
    gf_1p_mat_t<Model> Svec_init( FrequenciesCount::POS_1P_RANGE, Model::GetFineMomentaCount(), true );

    std::cout << " ... computing G and S" << std::endl;
    rhs_base_t::compute_Gvec(&Gvec, old_state, t);
    rhs_base_t::compute_Svec(&Svec_init, old_state, t);

    std::cout << " ... computing 1-loop self-energy derivative" << std::endl;

    #pragma omp parallel default(shared)
    {

    #ifdef SELFEN_FLOW

    #pragma omp single
	std::cout << " ... self-energy " << std::endl;

        rhs_sbe_1lfrg_t<Model, state_t>::rhs_Sig_1l(&dstate_over_dt.gf_Sig(), old_state, Svec_init, t);
    #else 
	dstate_over_dt.gf_Sig().init( [](const idx_1p_t<Model> &idx){return 0;} );

    #endif
    }

    gf_bubble_mat_t<Model> bubble_GS_pp, bubble_GS_ph;
    gf_bubble_mat_t<Model> bubble_GG_pp, bubble_GG_ph; // needed for multiloop

    state_t dstate_over_dt_1l;
    dstate_over_dt_1l.init_zero();

    double self_energy_iteration_error;
    for (unsigned se_loop_i = 0; se_loop_i < SELFENERGY_ITERATIONS_MAX; se_loop_i++){
	std::cout << " ... computing self-energy iteration: " << se_loop_i << std::endl;

	std::cout << " ... adding the katanin correction to S" << std::endl;
	rhs_sbe_1lfrg_t<Model, state_t>::add_katanin_correction(&Svec, Gvec, Svec_init, dstate_over_dt);

	std::cout << " ... computing GS and GG bubbles" << std::endl;
	rhs_sbe_mfrg_t<Model, state_t>::compute_GG_and_GS_bubbles(Gvec, Svec, &bubble_GS_pp, &bubble_GS_ph, &bubble_GG_pp, &bubble_GG_ph);

	std::cout << " ... computing 1-loop fRG vertex derivative " << std::endl;
	rhs_sbe_1lfrg_t<Model, state_t>::rhs_frg_1l(&dstate_over_dt_1l, old_state, bubble_GS_pp
					   , bubble_GS_ph, Svec_init, t);
	
	dstate_over_dt = dstate_over_dt_1l;
	
	std::cout << " ... adding multiloop vertex corrections " << std::endl;
	rhs_sbe_mfrg_t<Model, state_t>::add_multiloop_vertex_corrections(EXTRA_LOOP_NUM, &dstate_over_dt, old_state, dstate_over_dt_1l, bubble_GG_pp, bubble_GG_ph, t);
	
	std::cout << " ... adding self-energy corrections " << std::endl;
	add_multiloop_self_energy_corrections(&self_energy_iteration_error, &dstate_over_dt, old_state, Gvec, Svec, t);

	if (self_energy_iteration_error < ABS_ERROR_SELFENERGY_ITERATIONS)
	    break;

    }
    std::cout << " ... Correlated free energy" << std::endl;
    rhs_F_m(&dstate_over_dt.gf_f(), old_state, Gvec, Svec, t);
}

template <typename Model, typename state_t >
void rhs_sbe_mfrg_with_F_t<Model, state_t >::rhs_F_m(gf_f_t<Model> *new_F_ptr, const state_t &old_state, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, const double t)
{
  new_F_ptr->init( [this, &old_state, &Gvec, &Svec, t]( const idx_f_t<Model>& idx ){ return eval_F_m( idx, old_state, Gvec, Svec, t ); } ); 
}

