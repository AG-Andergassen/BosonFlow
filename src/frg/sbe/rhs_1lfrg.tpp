#include <mymath.h>
#include <cmath>
#include <complex>
#include <frg/sbe/symmetries.h>
#include <frg/sbe/symmetries_bosonised_M.h>
#include <frg/sbe/interpolators.h>
#include <frg/flows.h>

// do we need?
#include <frg/sbe/symmetries_mfrg.h>


using namespace std; 
using boost::bind;


template <typename Model, typename state_t >
rhs_sbe_1lfrg_t<Model, state_t>::~rhs_sbe_1lfrg_t()
{

}


template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::operator() ( const state_t& old_state, state_t& dstate_over_dt, const double t ){    
    rhs_base_t::print_current_scale(t);

    gf_1p_mat_t<Model> Gvec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true ); 
    gf_1p_mat_t<Model> Svec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true );
    
    rhs_base_t::compute_Gvec(&Gvec, old_state, t);
    rhs_base_t::compute_Svec(&Svec, old_state, t);
    
    std::cout << " Rhs calculation " << std::endl;

    gf_1p_mat_real_t<Model> Gvec_realspace(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true), 
	Svec_realspace(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true);   
    gf_bubble_mat_t<Model> bubble_GS_pp, bubble_GS_ph;
    gf_bubble_mat_t<Model> bubble_GG_pp, bubble_GG_ph;

    cout << " ... GS bubble " << endl;
    
#ifdef DEBUG_BUILD
    //    omp_set_num_threads(1);  	
#endif

    #pragma omp parallel default(shared) 
    {
	rhs_base_t::compute_Gvec_real(&Gvec_realspace, Gvec);
	rhs_base_t::compute_Svec_real(&Svec_realspace, Svec);

#ifdef STATIC_CALCULATION
#pragma omp barrier //barrier to make sure G/Svec_real are present
	//rhs_base_t::compute_static_GS_bubbles(&bubble_GS_pp, &bubble_GS_ph, t);
	rhs_base_t::compute_GS_bubbles(&bubble_GS_pp, &bubble_GS_ph, Gvec_realspace, Svec_realspace);
	
#else
#pragma omp barrier //barrier to make sure G/Svec_real are present
	rhs_base_t::compute_GS_bubbles(&bubble_GS_pp, &bubble_GS_ph, Gvec_realspace, Svec_realspace);
#endif
    }

    rhs_frg_1l(&dstate_over_dt, old_state, bubble_GS_pp, bubble_GS_ph, bubble_GG_pp, bubble_GG_ph, Svec, t);
}  


template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_frg_1l(state_t *dstate_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const gf_1p_mat_t<Model>& Svec, const double t)
{
    state_t &dstate_over_dt = *dstate_over_dt_ptr;

#if (defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)) && defined(F_NONZERO)
#pragma omp parallel default(shared)
    {
	precalculate_lambda_sc_left_Fdot_part(old_state, bubble_GG_pp, t);
	precalculate_lambda_d_left_Fdot_part(old_state, bubble_GG_ph, t);
	precalculate_lambda_m_left_Fdot_part(old_state, bubble_GG_ph, t);
#if !defined(SBEa_APPROXIMATION)

	precalculate_M_sc_left_Fdot_part(old_state, bubble_GG_pp, t);
	precalculate_M_d_left_Fdot_part(old_state, bubble_GG_ph, t);
	precalculate_M_m_left_Fdot_part(old_state, bubble_GG_ph, t);
#endif //!defined(SBEa_APPROXIMATION)
    }

#endif


    #pragma omp parallel default(shared)
    {

#ifdef SELFEN_FLOW 
	
#pragma omp single
	cout << " ... self-energy " << endl;
	
	rhs_Sig_1l(&dstate_over_dt.gf_Sig(), old_state, Svec, t);
#else 
	dstate_over_dt.gf_Sig().init( [](const idx_1p_t<Model> &idx){return 0;} );
#endif

#pragma omp single
	cout << " ... w (screened interaction) " << endl;
	    
	rhs_w_sc_1l(&dstate_over_dt.gf_w_sc(), old_state, bubble_GS_pp, bubble_GG_pp, t);
	rhs_w_d_1l(&dstate_over_dt.gf_w_d(), old_state, bubble_GS_ph, bubble_GG_ph, t);
	rhs_w_m_1l(&dstate_over_dt.gf_w_m(), old_state, bubble_GS_ph, bubble_GG_ph, t);

#ifndef NO_HEDIN_VERTEX_FLOW

#pragma omp single
	cout << " ... lambda (Hedin vertex)" << endl;

	rhs_lambda_sc_1l(&dstate_over_dt.gf_lambda_sc(), old_state, bubble_GS_pp, bubble_GG_pp, t); 
	rhs_lambda_d_1l(&dstate_over_dt.gf_lambda_d(), old_state, bubble_GS_ph, bubble_GG_ph, t);
	rhs_lambda_m_1l(&dstate_over_dt.gf_lambda_m(), old_state, bubble_GS_ph, bubble_GG_ph, t);

#endif
	

#if !defined(SBEa_APPROXIMATION) && !defined(SBEb_APPROXIMATION) && !defined(BOSONISE_M)// neglecting the rest fucntion flow; i.e. the SBEa approximation. In the SBEb approximation, there is too no M to flow.

#pragma omp single
	cout << " ... M (Rest functions)" << endl;

	rhs_M_sc_1l(&dstate_over_dt.gf_M_sc(), old_state, bubble_GS_pp, bubble_GG_pp, t);
	rhs_M_d_1l(&dstate_over_dt.gf_M_d(), old_state, bubble_GS_ph, bubble_GG_ph, t);
	rhs_M_m_1l(&dstate_over_dt.gf_M_m(), old_state, bubble_GS_ph, bubble_GG_ph, t);

#endif
    }

#if !defined(SBEa_APPROXIMATION) && !defined(SBEb_APPROXIMATION) && defined(BOSONISE_M)
    rhs_frg_bosonised_M_1l(dstate_over_dt_ptr, old_state, bubble_GS_pp, bubble_GS_ph, bubble_GG_pp, bubble_GG_ph, t);
#endif

    cout << " ... finished 1l " << endl;

    auto &the_state = old_state;

    coord_t<Model::dim> zero;
    unsigned p_zeroidx = std::get<0>(Model::MomentumGrid().ptr->get_idx_from_coord(zero));
    std::cout << "gf_w_sc at (W=0, K=0):" << the_state.gf_w_sc()[0][p_zeroidx] << std::endl;
    std::cout << "gf_lambda_sc at (W=0, K=0, w=0, m=0):" << the_state.gf_lambda_sc()[0][p_zeroidx][0][0] << std::endl;
//    std::cout << "gf_M_sc at (W=0, K=0, w=0, m=0, wp=0, mp=0):" << the_state.gf_M_sc()[0][p_zeroidx][0][0][0][0] << std::endl;
    std::cout << "gf_Sig at (w=0, k=0, s_in=0, s_out=0):" << the_state.gf_Sig()[0][p_zeroidx][0][0] << std::endl;

}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_frg_bosonised_M_1l(state_t *dstate_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_pp, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    state_t &dstate_over_dt_with_bosonised_M = *dstate_over_dt_ptr;
	
#pragma omp parallel default(shared)
    {
	SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_w_M_sc_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_wM_sc_plus()), [&old_state, &bubble_GS_pp, &bubble_GG_pp, t]( const idx_w_M_t<Model>& idx ){ return eval_wM_sc(idx, old_state, bubble_GS_pp, bubble_GG_pp, +1, t ); });
        SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_w_M_sc_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_wM_sc_minus()), [&old_state, &bubble_GS_pp, &bubble_GG_pp, t]( const idx_w_M_t<Model>& idx ){ return eval_wM_sc(idx, old_state, bubble_GS_pp, bubble_GG_pp, -1, t ); });
	SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_w_M_d_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_wM_d_plus()), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t]( const idx_w_M_t<Model>& idx ){ return eval_wM_d(idx, old_state, bubble_GS_ph, bubble_GG_ph, +1, t ); });
        SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_w_M_d_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_wM_d_minus()), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t]( const idx_w_M_t<Model>& idx ){ return eval_wM_d(idx, old_state, bubble_GS_ph, bubble_GG_ph, -1, t ); });
	SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_w_M_m_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_wM_m_plus()), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t]( const idx_w_M_t<Model>& idx ){ return eval_wM_m(idx, old_state, bubble_GS_ph, bubble_GG_ph, +1, t ); });
        SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_w_M_m_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_wM_m_minus()), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t]( const idx_w_M_t<Model>& idx ){ return eval_wM_m(idx, old_state, bubble_GS_ph, bubble_GG_ph, -1, t ); });

#pragma omp barrier

	SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_lambda_M_sc_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_lambdaM_sc_plus()), [&old_state, &bubble_GS_pp, &bubble_GG_pp, &dstate_over_dt_with_bosonised_M, t]( const idx_lambda_M_t<Model>& idx ){ return eval_lambdaM_sc(idx, old_state, bubble_GS_pp, bubble_GG_pp, +1, dstate_over_dt_with_bosonised_M.gf_wM_sc_plus(), t ); });
	SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_lambda_M_sc_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_lambdaM_sc_minus()), [&old_state, &bubble_GS_pp, &bubble_GG_pp, &dstate_over_dt_with_bosonised_M, t]( const idx_lambda_M_t<Model>& idx ){ return eval_lambdaM_sc(idx, old_state, bubble_GS_pp, bubble_GG_pp, -1, dstate_over_dt_with_bosonised_M.gf_wM_sc_minus(), t ); });

	SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_lambda_M_d_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_lambdaM_d_plus()), [&old_state, &bubble_GS_ph, &bubble_GG_ph, &dstate_over_dt_with_bosonised_M, t]( const idx_lambda_M_t<Model>& idx ){ return eval_lambdaM_d(idx, old_state, bubble_GS_ph, bubble_GG_ph, +1, dstate_over_dt_with_bosonised_M.gf_wM_d_plus(), t ); });
	SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_lambda_M_d_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_lambdaM_d_minus()), [&old_state, &bubble_GS_ph, &bubble_GG_ph, &dstate_over_dt_with_bosonised_M, t]( const idx_lambda_M_t<Model>& idx ){ return eval_lambdaM_d(idx, old_state, bubble_GS_ph, bubble_GG_ph, -1, dstate_over_dt_with_bosonised_M.gf_wM_d_minus(), t ); });

	SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_lambda_M_m_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_lambdaM_m_plus()), [&old_state, &bubble_GS_ph, &bubble_GG_ph, &dstate_over_dt_with_bosonised_M, t]( const idx_lambda_M_t<Model>& idx ){ return eval_lambdaM_m(idx, old_state, bubble_GS_ph, bubble_GG_ph, +1, dstate_over_dt_with_bosonised_M.gf_wM_m_plus(), t ); });
	SymmetriesfRGSBEBosonised_M<Model>::IdxEquivClasses_lambda_M_m_Ptr()->init_batched_at_sampling_indices((dstate_over_dt_with_bosonised_M.gf_lambdaM_m_minus()), [&old_state, &bubble_GS_ph, &bubble_GG_ph, &dstate_over_dt_with_bosonised_M, t]( const idx_lambda_M_t<Model>& idx ){ return eval_lambdaM_m(idx, old_state, bubble_GS_ph, bubble_GG_ph, -1, dstate_over_dt_with_bosonised_M.gf_wM_m_minus(), t ); });	
	
    }

    std::cout << "Difference Bosonised M d:" << dstate_over_dt_with_bosonised_M.wM_d_plus(0, 0, 0, 0) << "," << dstate_over_dt_with_bosonised_M.gf_M_d()[0][0][0][0][0][0] << std::endl;
    std::cout << "Difference Bosonised M m:" << dstate_over_dt_with_bosonised_M.wM_m_plus(0, 0, 0, 0) << "," << dstate_over_dt_with_bosonised_M.gf_M_m()[0][0][0][0][0][0] << std::endl;

    std::cout << "lambdaM d plus:" << dstate_over_dt_with_bosonised_M.lambdaM_d_plus(0, 0, 0, 0, 0) << std::endl;

    std::cout << "gf_Mbos_sc at (W=0, K=0, w=0, m=0, wp=0, mp=0):" << old_state.M_sc(0,0,0,0,0,0) << std::endl;
}
 
template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::add_katanin_correction(gf_1p_mat_t<Model> *Svec_ptr, gf_1p_mat_t<Model> &Gvec, gf_1p_mat_t<Model> &Svec_init, state_t &dstate_over_dt)
{
#if !defined(SELFEN_FLOW) && !defined(SELFEN_SDE_SBE_MAGNETIC) && !defined(SELFEN_SDE_SBE_DENSITY) && !defined(SELFEN_SDE_SBE_SUPERCONDUCTING)
    (*Svec_ptr) = Svec_init;
    return;
#endif

    gf_1p_mat_t<Model> SigDotvec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount() );
    // todo: this function doesn't depend on the entire state, but really just needs dstate_over_dt.gf_Sig()..
    SigDotvec.init( bind( &state_t::SigMat_big, boost::cref( dstate_over_dt ), _1) );
  
    (*Svec_ptr) = Svec_init + Gvec * SigDotvec * Gvec;  /**< Add Katanin correction */
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_w_sc_1l(gf_w_t<Model> *dw_sc_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_sc_Ptr()->init_batched_at_sampling_indices( (*dw_sc_over_dt_ptr), [&old_state, &bubble_GS_pp, &bubble_GG_pp, t, this]( const idx_w_t<Model>& idx ){ 
#ifndef SBEb_APPROXIMATION
      return eval_w_sc( idx, old_state, bubble_GS_pp, t )
	#if defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)
       	   + eval_w_sc_Udot_part( idx, old_state, bubble_GG_pp, t )
       	#endif
	; 
#else
      return eval_SBEb_w_sc( idx, old_state, bubble_GS_pp, t );
#endif
} ); 


#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_sc_Ptr()->init_batched_at_interpolating_indices( (*dw_sc_over_dt_ptr), [dw_sc_over_dt_ptr]( const idx_w_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::w_sc_Ptr()->eval_gf(*dw_sc_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_w_d_1l(gf_w_t<Model> *dw_d_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_d_Ptr()->init_batched_at_sampling_indices( (*dw_d_over_dt_ptr), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t, this]( const idx_w_t<Model>& idx ){ 
#ifndef SBEb_APPROXIMATION
      return eval_w_d( idx, old_state, bubble_GS_ph, t )
	#if defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)
	    +eval_w_d_Udot_part( idx, old_state, bubble_GG_ph, t )
	#endif
	; 
#else
    return eval_SBEb_w_d( idx, old_state, bubble_GS_ph, t ); 
#endif
} );


#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_d_Ptr()->init_batched_at_interpolating_indices( (*dw_d_over_dt_ptr), [dw_d_over_dt_ptr]( const idx_w_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::w_d_Ptr()->eval_gf(*dw_d_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_w_m_1l(gf_w_t<Model> *dw_m_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_m_Ptr()->init_batched_at_sampling_indices( (*dw_m_over_dt_ptr), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t, this]( const idx_w_t<Model>& idx ){ 
#ifndef SBEb_APPROXIMATION    
	    return eval_w_m( idx, old_state, bubble_GS_ph, t )
	      #if defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)
	         + eval_w_m_Udot_part( idx, old_state, bubble_GG_ph, t )
	      #endif
	      ;
#else
	    return eval_SBEb_w_m( idx, old_state, bubble_GS_ph, t );
#endif
    } );    


#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_m_Ptr()->init_batched_at_interpolating_indices( (*dw_m_over_dt_ptr), [dw_m_over_dt_ptr]( const idx_w_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::w_m_Ptr()->eval_gf(*dw_m_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_lambda_sc_1l(gf_lambda_t<Model> *dlambda_sc_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_sc_Ptr()->init_batched( (*dlambda_sc_over_dt_ptr), [&old_state, &bubble_GS_pp, &bubble_GG_pp, t, this]( const idx_lambda_t<Model>& idx ){ 
#ifndef SBEb_APPROXIMATION
	    return eval_lambda_sc(idx, old_state, bubble_GS_pp, t )
	      #if defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)
	         + eval_lambda_sc_Udot_part(idx, old_state, bubble_GG_pp, t )
	      #endif
	      ; 
#else
	    return eval_SBEb_lambda_sc(idx, old_state, bubble_GS_pp, t ); 
#endif	    
	} ); 
    
#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_sc_Ptr()->init_batched_at_interpolating_indices( (*dlambda_sc_over_dt_ptr), [dlambda_sc_over_dt_ptr]( const idx_lambda_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::lambda_sc_Ptr()->eval_gf(*dlambda_sc_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_lambda_d_1l(gf_lambda_t<Model> *dlambda_d_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_d_Ptr()->init_batched( (*dlambda_d_over_dt_ptr), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t, this]( const idx_lambda_t<Model>& idx ){ 
#ifndef SBEb_APPROXIMATION
	    return eval_lambda_d( idx, old_state, bubble_GS_ph, t )
	      #if defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)
	         + eval_lambda_d_Udot_part( idx, old_state, bubble_GG_ph, t )
	      #endif
	      ;
#else
	    return eval_SBEb_lambda_d( idx, old_state, bubble_GS_ph, t );
#endif
	} ); 

#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_d_Ptr()->init_batched_at_interpolating_indices( (*dlambda_d_over_dt_ptr), [dlambda_d_over_dt_ptr]( const idx_lambda_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::lambda_d_Ptr()->eval_gf(*dlambda_d_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_lambda_m_1l(gf_lambda_t<Model> *dlambda_m_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_m_Ptr()->init_batched_at_sampling_indices( (*dlambda_m_over_dt_ptr), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t, this]( const idx_lambda_t<Model>& idx ){ 
#ifndef SBEb_APPROXIMATION
	    //auto dbg_ = eval_lambda_m( idx, old_state, bubble_GS_ph, t );
	    //std::cout << dbg_ << std::endl;
	    return eval_lambda_m( idx, old_state, bubble_GS_ph, t )
	      #if defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)
		+ eval_lambda_m_Udot_part( idx, old_state,bubble_GG_ph, t )
	      #endif
	      ;
#else
	    return eval_SBEb_lambda_m( idx, old_state, bubble_GS_ph, t );
#endif
	} );


#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_m_Ptr()->init_batched_at_interpolating_indices( (*dlambda_m_over_dt_ptr), [dlambda_m_over_dt_ptr]( const idx_lambda_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::lambda_m_Ptr()->eval_gf(*dlambda_m_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
}


template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_M_sc_1l(gf_M_t<Model> *dM_sc_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched_at_sampling_indices( (*dM_sc_over_dt_ptr), [&old_state, &bubble_GS_pp, &bubble_GG_pp, t, this]( const idx_M_t<Model>& idx ){
      return eval_M_sc(idx, old_state, bubble_GS_pp, t )
	#if defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)
           + eval_M_sc_Udot_part(idx, old_state, bubble_GG_pp, t )
	#endif
	;
    } ); 

#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched_at_interpolating_indices( (*dM_sc_over_dt_ptr), [dM_sc_over_dt_ptr]( const idx_M_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::M_sc_Ptr()->eval_gf(*dM_sc_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_M_d_1l(gf_M_t<Model> *dM_d_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched_at_sampling_indices( (*dM_d_over_dt_ptr), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t, this]( const idx_M_t<Model>& idx ){
      return eval_M_d(idx, old_state, bubble_GS_ph, t )
	#if defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)
	   + eval_M_d_Udot_part(idx, old_state, bubble_GG_ph, t )
	#endif
	;
    } ); 


#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched_at_interpolating_indices( (*dM_d_over_dt_ptr), [dM_d_over_dt_ptr]( const idx_M_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::M_d_Ptr()->eval_gf(*dM_d_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_M_m_1l(gf_M_t<Model> *dM_m_over_dt_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched_at_sampling_indices( (*dM_m_over_dt_ptr), [&old_state, &bubble_GS_ph, &bubble_GG_ph, t, this]( const idx_M_t<Model>& idx ){
      return eval_M_m(idx, old_state, bubble_GS_ph, t )
	#if defined(ACTUAL_U_FLOW) || defined(FREE_U_FLOW)
	   + eval_M_m_Udot_part(idx, old_state, bubble_GG_ph, t )
	#endif
	;
    } ); 

#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched_at_interpolating_indices( (*dM_m_over_dt_ptr), [dM_m_over_dt_ptr]( const idx_M_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::M_m_Ptr()->eval_gf(*dM_m_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
    
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_Sig_1l(gf_1p_t<Model> *dSig_over_dt_ptr, const state_t &old_state, const gf_1p_mat_t<Model>& Svec, const double t)
{
	    Symmetries1lfRGSBE<Model>::IdxEquivClasses_sig_Ptr()->init_batched_at_sampling_indices( (*dSig_over_dt_ptr), [&old_state, &Svec, t, this]( const idx_1p_t<Model>& idx ){ return eval_Sig_conv( idx, old_state, Svec, t ); } );


#pragma omp barrier
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_sig_Ptr()->init_batched_at_interpolating_indices( (*dSig_over_dt_ptr), [dSig_over_dt_ptr]( const idx_1p_t<Model>& idx ){ 
	    return Interpolators1lfRGSBE<Model>::sig_Ptr()->eval_gf(*dSig_over_dt_ptr, idx, std::complex<double>(0, 0));
	} );
}


template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_wM_sc_1l(gf_w_M_t<Model> *dwM_sc_plus_over_dt_ptr, gf_w_M_t<Model> *dwM_sc_minus_over_dt_ptr, const state_frg_sbe_bosonised_M_t<Model> &old_state, const gf_bubble_mat_t<Model> &bubble_GS_pp, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    dwM_sc_plus_over_dt_ptr->init([&old_state, &bubble_GS_pp, &bubble_GG_pp, t](const idx_w_M_t<Model>& idx){ return eval_wM_sc(idx, old_state, bubble_GS_pp, +1, t); });
    
    dwM_sc_minus_over_dt_ptr->init([&old_state, &bubble_GS_pp, &bubble_GG_pp, t](const idx_w_M_t<Model>& idx){ return eval_wM_sc(idx, old_state, bubble_GS_pp, -1, t); });
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_wM_d_1l(gf_w_M_t<Model> *dwM_d_plus_over_dt_ptr, gf_w_M_t<Model> *dwM_d_minus_over_dt_ptr, const state_frg_sbe_bosonised_M_t<Model> &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    dwM_d_plus_over_dt_ptr->init([&old_state, &bubble_GS_ph, &bubble_GG_ph, t](const idx_w_M_t<Model>& idx){ return eval_wM_d(idx, old_state, bubble_GS_ph, +1, t); });
    
    dwM_d_minus_over_dt_ptr->init([&old_state, &bubble_GS_ph, &bubble_GG_ph, t](const idx_w_M_t<Model>& idx){ return eval_wM_d(idx, old_state, bubble_GS_ph, -1, t); });
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::rhs_wM_m_1l(gf_w_M_t<Model> *dwM_m_plus_over_dt_ptr, gf_w_M_t<Model> *dwM_m_minus_over_dt_ptr, const state_frg_sbe_bosonised_M_t<Model> &old_state, const gf_bubble_mat_t<Model> &bubble_GS_ph, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    dwM_m_plus_over_dt_ptr->init([&old_state, &bubble_GS_ph, &bubble_GG_ph, t](const idx_w_M_t<Model>& idx){ return eval_wM_m(idx, old_state, bubble_GS_ph, +1, t); });
    
    dwM_m_minus_over_dt_ptr->init([&old_state, &bubble_GS_ph, &bubble_GG_ph, t](const idx_w_M_t<Model>& idx){ return eval_wM_m(idx, old_state, bubble_GS_ph, -1, t); });
}

// for U-flow when the EXTENDED flag is one (\mathcal{F} \neq 0)
template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::precalculate_lambda_sc_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_lambda_sc_LeftCorrections_Ptr()->init_batched( m_lambda_sc_left_Fdot_part, [this, &old_state, &bubble_GG_pp, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_sc_left_Fdot_part( idx, old_state, bubble_GG_pp, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::precalculate_lambda_d_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_lambda_d_LeftCorrections_Ptr()->init_batched( m_lambda_d_left_Fdot_part, [this, &old_state, &bubble_GG_ph, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_d_left_Fdot_part( idx, old_state, bubble_GG_ph, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::precalculate_lambda_m_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_lambda_m_LeftCorrections_Ptr()->init_batched( m_lambda_m_left_Fdot_part, [this, &old_state, &bubble_GG_ph, t]( const idx_lambda_t<Model>& idx ){  return eval_lambda_m_left_Fdot_part( idx, old_state, bubble_GG_ph, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::precalculate_M_sc_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_pp, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_M_sc_LeftCorrections_Ptr()->init_batched( m_M_sc_left_Fdot_part, [this, &old_state, &bubble_GG_pp, t]( const idx_M_t<Model>& idx ){ return eval_M_sc_left_Fdot_part( idx, old_state, bubble_GG_pp, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::precalculate_M_d_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_M_d_LeftCorrections_Ptr()->init_batched( m_M_d_left_Fdot_part, [this, &old_state, &bubble_GG_ph, t]( const idx_M_t<Model>& idx ){ return eval_M_d_left_Fdot_part( idx, old_state, bubble_GG_ph, t ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_1lfrg_t<Model, state_t>::precalculate_M_m_left_Fdot_part(const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_GG_ph, const double t)
{
    SymmetriesmfRGSBE<Model>::IdxEquivClasses_M_m_LeftCorrections_Ptr()->init_batched( m_M_m_left_Fdot_part, [this, &old_state, &bubble_GG_ph, t]( const idx_M_t<Model>& idx ){ return eval_M_m_left_Fdot_part( idx, old_state, bubble_GG_ph, t ); } ); 
}
