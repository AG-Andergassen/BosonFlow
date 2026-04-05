#include <frg/sbe/symmetries.h>
#include <frg/flows.h>
#include <iostream>

#include <params_technical.h>

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_with_F_t<Model, state_t>::operator() ( const state_t& old_state, state_t& dstate_over_dt, const double t ){
    std::cout << " Preprocessing at scale " << t << std::endl;

#ifdef DEBUG_BUILD
    omp_set_num_threads(1);  	
#endif

    
    gf_1p_mat_t<Model> Gvec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true ); 
    gf_1p_mat_t<Model> Svec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true );

    std::cout << " ... computing G and S" << std::endl;
    rhs_base_t::compute_Gvec(&Gvec, old_state, t);
    rhs_base_t::compute_Svec(&Svec, old_state, t);
   
    gf_bubble_mat_t<Model> bubble_GS_pp, bubble_GS_ph;

    gf_1p_mat_real_t<Model> Gvec_realspace(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true), Svec_realspace(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true);   
    
    #pragma omp parallel default(shared) 
    {
	rhs_base_t::compute_Gvec_real(&Gvec_realspace, Gvec);
	rhs_base_t::compute_Svec_real(&Svec_realspace, Svec);

#ifdef STATIC_CALCULATION
	rhs_base_t::compute_static_GS_bubbles(&bubble_GS_pp, &bubble_GS_ph, t);
#else
	if (false){
#pragma omp barrier //barrier to make sure G/Svec_real are present
	rhs_base_t::compute_GS_bubbles_trad(&bubble_GS_pp, &bubble_GS_ph, Gvec, Svec);
	}else{
#pragma omp barrier //barrier to make sure G/Svec_real are present
	rhs_base_t::compute_GS_bubbles(&bubble_GS_pp, &bubble_GS_ph, Gvec_realspace, Svec_realspace);
	}
#endif
    }

    std::cout << " ... Correlated free energy" << std::endl;
    rhs_F_1l(&dstate_over_dt.gf_f(), old_state, Gvec, Svec, t);
    
    rhs_sbe_1lfrg_t<Model, state_t>::rhs_frg_1l(&dstate_over_dt, old_state, bubble_GS_pp, bubble_GS_ph, Svec, t);

    std::cout << "gf_F: " << old_state.gf_f() << std::endl;
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_with_F_t<Model, state_t>::rhs_F_1l(gf_f_t<Model> *new_F_ptr, const state_t &old_state, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, const double t)
{
  new_F_ptr->init( [this, &old_state, &Gvec, &Svec, t]( const idx_f_t<Model>& idx ){ return eval_F_1l( idx, old_state, Gvec, Svec, t ); } ); 
}

