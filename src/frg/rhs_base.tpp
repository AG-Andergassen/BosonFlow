#include <frg/rhs_base.h>
#include <frg/flows.h>
#include <frg/symmetries_common.h>
#include <frequencies/matsubara_space.h>
#include <frg/interpolators_common.h>
#include "cubature.h"

//#define PCUBATURE

#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif

// compute Gvec without a scale
template <typename Model, typename State>
void rhs_base_t<Model, State>::compute_Gvec(gf_1p_mat_t<Model> *Gvec_ptr, const State &old_state)
{
    auto eval_G = [&old_state](const idx_1p_mat_t<Model>& idx){
	int idx_p = idx(I1PMAT::p);
	int idx_k = Model::GetCoarseMomentumIdxFromFine(idx_p);
	int idx_w = idx(I1PMAT::w);

	MatQN G0inv;
	double w = w_val(idx_w);
	G0inv << I*w - Model::E(idx_p) + old_state.gf_delta_mu();
	return ( G0inv - old_state.SigMat(idx_w, idx_k) ).inverse();
    };

    // todo: check whether using symmetries for G is an option
    (*Gvec_ptr).init( eval_G );  
}


template <typename Model, typename State>
void rhs_base_t<Model, State>::compute_Gvec(gf_1p_mat_t<Model> *Gvec_ptr, const State &old_state, const double t)
{
    auto eval_G = [t, &old_state](const idx_1p_mat_t<Model>& idx){
	int idx_p = idx(I1PMAT::p);
	int idx_k = Model::GetCoarseMomentumIdxFromFine(idx_p);
	int idx_w = idx(I1PMAT::w);
	return fRGFlowScheme<Model>::G( idx_w, idx_p, t, old_state.SigMat( idx_w, idx_k ), old_state.gf_delta_mu() );
    };

    // todo: check whether using symmetries for G is an option
    (*Gvec_ptr).init( eval_G );  
}

template <typename Model, typename State>
void rhs_base_t<Model, State>::compute_Svec(gf_1p_mat_t<Model> *Svec_ptr, const State &old_state, const double t)
{
    auto eval_S = [t, &old_state](const idx_1p_mat_t<Model>& idx){
	int idx_p = idx(I1PMAT::p);
	int idx_k = Model::GetCoarseMomentumIdxFromFine(idx_p);
	int idx_w = idx(I1PMAT::w);
	return fRGFlowScheme<Model>::S( idx_w, idx_p, t, old_state.SigMat( idx_w, idx_k ), old_state.gf_delta_mu(), old_state.m_d_delta_mu_over_dt );
    };

    // todo: check whether using symmetries for S is an option
    (*Svec_ptr).init( eval_S ); 
}

template <typename Model, typename State>
void rhs_base_t<Model, State>::compute_Gvec_real(gf_1p_mat_real_t<Model> *Gvec_real_ptr, const gf_1p_mat_t<Model>& Gvec)
{
    SymmetriesfRGCommon<Model>::IdxEquivClasses_Gvec_real_Ptr()->init_batched_at_sampling_indices( (*Gvec_real_ptr), [&Gvec]( const idx_1p_mat_real_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_Gvec_real( idx, Gvec ); } );    /**< Fourier transform Green's function to real space */


#pragma omp barrier
    MatReal zero;
    zero.setZero(Model::GetFineMomentaCount());
    SymmetriesfRGCommon<Model>::IdxEquivClasses_Gvec_real_Ptr()->init_batched_at_interpolating_indices( (*Gvec_real_ptr), [Gvec_real_ptr, &zero]( const idx_1p_mat_real_t<Model>& idx ){ 
	    return InterpolatorsfRGCommon<Model>::Gvec_real_Ptr()->eval_gf(*Gvec_real_ptr, idx, zero);
	} ); 
}

template <typename Model, typename State>
void rhs_base_t<Model, State>::compute_Svec_real(gf_1p_mat_real_t<Model> *Svec_real_ptr, const gf_1p_mat_t<Model>& Svec)
{
    SymmetriesfRGCommon<Model>::IdxEquivClasses_Gvec_real_Ptr()->init_batched_at_sampling_indices( (*Svec_real_ptr), [&Svec]( const idx_1p_mat_real_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_Gvec_real( idx, Svec ); } );    /**< Fourier transform single-scale propagator to real space */


#pragma omp barrier
    MatReal zero;
    zero.setZero(Model::GetFineMomentaCount());
    SymmetriesfRGCommon<Model>::IdxEquivClasses_Gvec_real_Ptr()->init_batched_at_interpolating_indices( (*Svec_real_ptr), [Svec_real_ptr, &zero]( const idx_1p_mat_real_t<Model>& idx ){ 
	    return InterpolatorsfRGCommon<Model>::Gvec_real_Ptr()->eval_gf(*Svec_real_ptr, idx, zero);
	} ); 
}

template <typename Model, typename State>
void rhs_base_t<Model, State>::compute_GS_bubbles(gf_bubble_mat_t<Model> *bubble_GS_pp_ptr, gf_bubble_mat_t<Model> *bubble_GS_ph_ptr, const gf_1p_mat_real_t<Model> &Gvec_real, const gf_1p_mat_real_t<Model> &Svec_real){
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_pp_Ptr()->init_batched_at_sampling_indices( (*bubble_GS_pp_ptr), [&Svec_real, &Gvec_real]( const idx_bubble_mat_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_diag_bubble_pp( idx, Svec_real, Gvec_real ); } ); /**< Calculation of the pp excitation with scale derivative */
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_ph_Ptr()->init_batched_at_sampling_indices( (*bubble_GS_ph_ptr), [&Svec_real, &Gvec_real]( const idx_bubble_mat_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_diag_bubble_ph( idx, Svec_real, Gvec_real ); } ); /**< Calculation of the ph excitation with scale derivative */

#pragma omp barrier
    MatPatch zero;
    zero.setZero(Model::GetRefinedMomentaCount());
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_pp_Ptr()->init_batched_at_interpolating_indices( (*bubble_GS_pp_ptr), [bubble_GS_pp_ptr, &zero]( const idx_bubble_mat_t<Model>& idx ){ 
	    return InterpolatorsfRGCommon<Model>::bubblemat_pp_Ptr()->eval_gf(*bubble_GS_pp_ptr, idx, zero); 
	} ); 
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_ph_Ptr()->init_batched_at_interpolating_indices( (*bubble_GS_ph_ptr), [bubble_GS_ph_ptr, &zero]( const idx_bubble_mat_t<Model>& idx ){ 
	    return InterpolatorsfRGCommon<Model>::bubblemat_ph_Ptr()->eval_gf(*bubble_GS_ph_ptr, idx, zero); 
	} );
}

template <typename Model, typename State>
void rhs_base_t<Model, State>::compute_GG_bubbles(gf_bubble_mat_t<Model> *bubble_GG_pp_ptr, gf_bubble_mat_t<Model> *bubble_GG_ph_ptr, const gf_1p_mat_real_t<Model> &Gvec_real)
{
    //todo: really, one should be able to reuse the above function by passing two Gs. Clean it up so that it's possible.
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_pp_Ptr()->init_batched_at_sampling_indices( (*bubble_GG_pp_ptr), [&Gvec_real]( const idx_bubble_mat_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_diag_bubble_GG_pp( idx, Gvec_real ); } );  /**< Calculation of the pp excitation */
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_ph_Ptr()->init_batched_at_sampling_indices( (*bubble_GG_ph_ptr), [&Gvec_real]( const idx_bubble_mat_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_diag_bubble_GG_ph( idx, Gvec_real ); } );  /**< Calculation of the ph excitation */


#pragma omp barrier
    MatPatch zero;
    zero.setZero(Model::GetRefinedMomentaCount());
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_pp_Ptr()->init_batched_at_interpolating_indices( (*bubble_GG_pp_ptr), [bubble_GG_pp_ptr, &zero]( const idx_bubble_mat_t<Model>& idx ){ 
	    return InterpolatorsfRGCommon<Model>::bubblemat_pp_Ptr()->eval_gf(*bubble_GG_pp_ptr, idx, zero); 
	} ); 
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_ph_Ptr()->init_batched_at_interpolating_indices( (*bubble_GG_ph_ptr), [bubble_GG_ph_ptr, &zero]( const idx_bubble_mat_t<Model>& idx ){ 
	    return InterpolatorsfRGCommon<Model>::bubblemat_ph_Ptr()->eval_gf(*bubble_GG_ph_ptr, idx, zero); 
	} );

}

template <typename Model, typename State>
MatReal rhs_base_t<Model, State>::eval_Gvec_real( const idx_1p_mat_real_t<Model>& idx, const gf_1p_mat_t<Model>& Gvec) 
{
    int w  = idx( I1PMATREAL::w );
    int s1 = idx( I1PMATREAL::s_in );
    int s2 = idx( I1PMATREAL::s_out );  
    MatReal val_g = MatReal(Model::GetFineMomentaCount());
        
    for(int p_idx = 0; p_idx < Model::GetFineMomentaCount(); ++p_idx){
	Model::FFTW(omp_get_thread_num()).input[p_idx][0] = Gvec[w][p_idx](s1,s2).real(); 
	Model::FFTW(omp_get_thread_num()).input[p_idx][1] = Gvec[w][p_idx](s1,s2).imag(); 
    }

    fftw_execute_dft(Model::FFTW(omp_get_thread_num()).backward_plan, Model::FFTW(omp_get_thread_num()).input, Model::FFTW(omp_get_thread_num()).output);
       
    for(int R = 0; R < Model::GetFineMomentaCount(); ++R){
        val_g(R) =  1.0/Model::GetFineMomentaCount() * (Model::FFTW(omp_get_thread_num()).output[R][0] + I * Model::FFTW(omp_get_thread_num()).output[R][1]);
    }
    return val_g; 
}


template <typename Model, typename State>
void rhs_base_t<Model, State>::compute_GS_bubbles_trad(gf_bubble_mat_t<Model> *bubble_GS_pp_ptr, gf_bubble_mat_t<Model> *bubble_GS_ph_ptr, const gf_1p_mat_t<Model> &Gvec, const gf_1p_mat_t<Model> &Svec){
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_pp_Ptr()->init_batched( (*bubble_GS_pp_ptr), [&Svec, &Gvec]( const idx_bubble_mat_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_diag_bubble_pp_trad( idx, Svec, Gvec ); } ); 
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_ph_Ptr()->init_batched( (*bubble_GS_ph_ptr), [&Svec, &Gvec]( const idx_bubble_mat_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_diag_bubble_ph_trad( idx, Svec, Gvec ); } );
}


template <typename Model, typename State>
MatPatch rhs_base_t<Model, State>::eval_diag_bubble_pp_trad( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_t<Model>& G, const gf_1p_mat_t<Model>& S )
{
    int W   = idx( IBUBMAT::W );
    int w   = idx( IBUBMAT::w );
    int m   = idx( IBUBMAT::m );
    int n   = idx( IBUBMAT::n );
    int s4  = idx( IBUBMAT::s1 );
    int s4p = idx( IBUBMAT::s1p); 
    int s3  = idx( IBUBMAT::s2 );
    int s3p = idx( IBUBMAT::s2p );

#ifndef STATIC_CALCULATION
    int W_over_2_minus_w_minus_1 = div2_floor( W ) - w - 1;
    int w_plus_W_over_2 = w + div2_ceil( W );
#else
    int W_over_2_minus_w_minus_1 = 0;
    int w_plus_W_over_2 = 0;	    
#endif

    
    MatPatch val_patch;

    const unsigned refined_count = Model::GetRefinedMomentaCount();
    const unsigned fine_count = Model::GetFineMomentaCount();
    val_patch.setZero(refined_count);

    std::vector<unsigned> q_fine_idx_map(refined_count);
    for (unsigned q_idx = 0; q_idx < refined_count; ++q_idx) {
    q_fine_idx_map[q_idx] = Model::GetFineMomentumIdxFromCoarse(q_idx);
    }

    std::vector<unsigned> minus_p_idx_map(fine_count);
    for (unsigned p_idx = 0; p_idx < fine_count; ++p_idx) {
    minus_p_idx_map[p_idx] = Model::GetNegativeFineMomentumIdx(p_idx);
    }

    const auto& ff_m = Model::GetFormFactorInFineMomentumIdxSpace(m);
    const auto& ff_n = Model::GetFormFactorInFineMomentumIdxSpace(n);

    for (unsigned q_idx = 0; q_idx < refined_count; ++q_idx){
    const unsigned q_fine_idx = q_fine_idx_map[q_idx];
    for (unsigned p_idx = 0; p_idx < fine_count; ++p_idx){
        const unsigned q_minus_p_idx = Model::SumFineMomentaIdxes(q_fine_idx, minus_p_idx_map[p_idx]);

        val_patch(q_idx) +=  (G[ w_plus_W_over_2 ][p_idx](s4, s4p) * S[ W_over_2_minus_w_minus_1 ][q_minus_p_idx](s3, s3p) +
                  S[ w_plus_W_over_2 ][p_idx](s4, s4p) * G[ W_over_2_minus_w_minus_1 ][q_minus_p_idx](s3, s3p) )
         * std::conj(ff_m[p_idx])
         * ff_n[p_idx];
    }
    }

    return val_patch * 1.0/Model::GetFineMomentaCount();
}

template <typename Model, typename State>
MatPatch rhs_base_t<Model, State>::eval_diag_bubble_ph_trad( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_t<Model>& G, const gf_1p_mat_t<Model>& S )
{
    int W   = idx( IBUBMAT::W );
    int w   = idx( IBUBMAT::w );
    int m   = idx( IBUBMAT::m );
    int n   = idx( IBUBMAT::n );
    int s4  = idx( IBUBMAT::s1 );
    int s4p = idx( IBUBMAT::s1p); 
    int s3  = idx( IBUBMAT::s2 );
    int s3p = idx( IBUBMAT::s2p );
    
#ifndef STATIC_CALCULATION
    int w_minus_W = w - div2_floor( W );
    int w_plus_W = div2_ceil( W ) + w;
#else
    int w_minus_W = 0;
    int w_plus_W = 0;
#endif

    
    MatPatch val_patch;

    const unsigned refined_count = Model::GetRefinedMomentaCount();
    const unsigned fine_count = Model::GetFineMomentaCount();
    val_patch.setZero(refined_count);

    std::vector<unsigned> q_fine_idx_map(refined_count);
    for (unsigned q_idx = 0; q_idx < refined_count; ++q_idx) {
	q_fine_idx_map[q_idx] = Model::GetFineMomentumIdxFromCoarse(q_idx);
    }

    const auto& ff_m = Model::GetFormFactorInFineMomentumIdxSpace(m);
    const auto& ff_n = Model::GetFormFactorInFineMomentumIdxSpace(n);

    for (unsigned q_idx = 0; q_idx < refined_count; q_idx ++){
	const unsigned q_fine_idx = q_fine_idx_map[q_idx];
	for (unsigned p_idx = 0; p_idx < fine_count; p_idx ++){
	    const unsigned p_plus_q_idx = Model::SumFineMomentaIdxes(q_fine_idx, p_idx);

	    val_patch(q_idx) +=  (G[w_minus_W][p_idx](s4, s4p) * S[w_plus_W][p_plus_q_idx](s3, s3p) + S[w_minus_W][p_idx](s4, s4p) * G[w_plus_W][p_plus_q_idx](s3, s3p)) * std::conj(ff_m[p_idx]) * ff_n[p_idx];
	}
    }

    return val_patch * 1.0/Model::GetFineMomentaCount();
}


// Remark: C.Hille's thesis Eq. (2.62b)
template <typename Model, typename State>
MatPatch rhs_base_t<Model, State>::eval_diag_bubble_pp( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_real_t<Model>& Svec_real, const gf_1p_mat_real_t<Model>& Gvec_real )
{
    int W   = idx( IBUBMAT::W );
    int w   = idx( IBUBMAT::w );
    int m   = idx( IBUBMAT::m );
    int n   = idx( IBUBMAT::n );
    int s4  = idx( IBUBMAT::s1 );
    int s4p = idx( IBUBMAT::s1p); 
    int s3  = idx( IBUBMAT::s2 );
    int s3p = idx( IBUBMAT::s2p );
    

    int W_over_2_minus_w_minus_1 = div2_floor( W ) - w - 1;
    int w_plus_W_over_2 = w + div2_ceil( W );
	
    MatPatch val_patch;

    const int thread_id = omp_get_thread_num();
    auto& fftw_workspace = Model::FFTW(thread_id);
    const unsigned refined_count = Model::GetRefinedMomentaCount();
    const unsigned fine_count = Model::GetFineMomentaCount();

    std::vector<unsigned> coarse_to_fine_idx_map(refined_count);
    for (unsigned K_idx = 0; K_idx < refined_count; ++K_idx) {
	coarse_to_fine_idx_map[K_idx] = Model::GetFineMomentumIdxFromCoarse(K_idx);
    }

    val_patch.setZero(refined_count);

#ifdef SET_MIXED_BUBBLES_TO_ZERO
    if (m != n)
	return val_patch;
#endif


#ifdef STATIC_CALCULATION
    // integrate the fermionic frequency of the bubbles
    for (w = -FrequencyDependenceScheme<Model>::BubbleTU_FermionicGrid().positive_freqs_count; w < FrequencyDependenceScheme<Model>::BubbleTU_FermionicGrid().positive_freqs_count; w++)
{
	W_over_2_minus_w_minus_1 = div2_floor( W ) - w -1;
        w_plus_W_over_2 = w + div2_ceil( W );

#endif
    for (unsigned r_idx : Model::FormFactors().fourier_transform_weights_r_support_indices){
	dcomplex weight = Model::FormFactors().fourier_transform_weights[r_idx][m][n];

        if(std::abs(weight) > CHOP_ERR){     /**< 1E-16 is error on the weights exponents, cosines,... */
	    	    
            const unsigned r_fine_idx = Model::FormFactors().real_space_to_real_fourier_transform_idx_map[r_idx];
            for(int R_idx = 0; R_idx < (int)fine_count; R_idx++){
		// we can use the momenta indices sum rules, since they're arranged in a similar manner in the primitive zone
    	      unsigned r_plus_R_idx = Model::SumFineMomentaIdxes(r_fine_idx, R_idx);		

		dcomplex temp = Gvec_real[ w_plus_W_over_2 ][s4][s4p](R_idx) * Svec_real[ W_over_2_minus_w_minus_1 ][s3][s3p](r_plus_R_idx) + Svec_real[ w_plus_W_over_2 ][s4][s4p](R_idx) * Gvec_real[ W_over_2_minus_w_minus_1 ][s3][s3p](r_plus_R_idx);
                    
		fftw_workspace.input[R_idx][0] = temp.real();
	        fftw_workspace.input[R_idx][1] = temp.imag();
	    }

            fftw_execute_dft(fftw_workspace.forward_plan, fftw_workspace.input, fftw_workspace.output);

            for (unsigned K_idx = 0; K_idx < refined_count; ++K_idx){
	      const unsigned P_idx = coarse_to_fine_idx_map[K_idx];

                val_patch(K_idx) += (fftw_workspace.output[P_idx][0] + I * fftw_workspace.output[P_idx][1])  *   weight * std::conj(Model::FineMomentumGrid().exp_ikr_idx_map[P_idx][r_idx]);
            }      
        }
    }

#ifdef STATIC_CALCULATION
    }
#endif
  
    return val_patch;

} 



template <typename Model, typename State>
MatPatch rhs_base_t<Model, State>::eval_diag_bubble_ph( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_real_t<Model>& Svec_real, const gf_1p_mat_real_t<Model>& Gvec_real)
{
    int W   = idx( IBUBMAT::W );
    int w   = idx( IBUBMAT::w );
    int m   = idx( IBUBMAT::m );
    int n   = idx( IBUBMAT::n );
    int s4  = idx( IBUBMAT::s1 );
    int s4p = idx( IBUBMAT::s1p);
    int s3  = idx( IBUBMAT::s2 );
    int s3p = idx( IBUBMAT::s2p );

    int w_minus_W_over_2 = w - div2_floor( W );
    int w_plus_W_over_2 = div2_ceil( W ) + w;	
   
    MatPatch val_patch;
    const int thread_id = omp_get_thread_num();
    auto& fftw_workspace = Model::FFTW(thread_id);
    const unsigned refined_count = Model::GetRefinedMomentaCount();
    const unsigned fine_count = Model::GetFineMomentaCount();

    std::vector<unsigned> coarse_to_fine_idx_map(refined_count);
    for (unsigned K_idx = 0; K_idx < refined_count; ++K_idx) {
	coarse_to_fine_idx_map[K_idx] = Model::GetFineMomentumIdxFromCoarse(K_idx);
    }

    std::vector<unsigned> minus_fine_idx_map(fine_count);
    for (unsigned R_idx = 0; R_idx < fine_count; ++R_idx) {
	minus_fine_idx_map[R_idx] = Model::GetNegativeFineMomentumIdx(R_idx);
    }

    val_patch.setZero(refined_count);

#ifdef SET_MIXED_BUBBLES_TO_ZERO
    if (m != n)
	return val_patch;
#endif


#ifdef STATIC_CALCULATION
    for (w = -FrequencyDependenceScheme<Model>::BubbleTU_FermionicGrid().positive_freqs_count; w < FrequencyDependenceScheme<Model>::BubbleTU_FermionicGrid().positive_freqs_count; w++)
{
    w_minus_W_over_2 = w - div2_floor( W );
    w_plus_W_over_2 = div2_ceil( W ) + w;	
#endif    

    // look into Nico's stuff

    // todo: make support indices depend on m, n
    for (unsigned r_idx : Model::FormFactors().fourier_transform_weights_r_support_indices){
	/*#pragma omp single
	  std::cout << r_idx << std::endl;*/
	dcomplex weight = Model::FormFactors().fourier_transform_weights[r_idx][m][n];
        if(std::abs(weight) > CHOP_ERR){    /**< 1E-16 is error on the weights exponents, cosines,...*/
        const unsigned r_fine_idx = Model::FormFactors().real_space_to_real_fourier_transform_idx_map[r_idx];
        const unsigned minus_r_fine_idx = Model::GetNegativeFineMomentumIdx(r_fine_idx);
        for(int R_idx = 0; R_idx < (int)fine_count; R_idx++){
		// we can use the momenta indices sum/minus rules, since they're arranged in a similar manner in the primitive zone
          const unsigned R_minus_r_idx = Model::SumFineMomentaIdxes(R_idx, minus_r_fine_idx);
          const unsigned minus_R_idx = minus_fine_idx_map[R_idx];
		
		dcomplex temp = Gvec_real[ w_minus_W_over_2 ][s4][s4p](minus_R_idx) * Svec_real[ w_plus_W_over_2 ][s3][s3p](R_minus_r_idx) + Svec_real[ w_minus_W_over_2 ][s4][s4p](minus_R_idx) * Gvec_real[ w_plus_W_over_2 ][s3][s3p](R_minus_r_idx);
                    
        fftw_workspace.input[R_idx][0] = temp.real();
        fftw_workspace.input[R_idx][1] = temp.imag();
	    }
            
            fftw_execute_dft(fftw_workspace.forward_plan, fftw_workspace.input, fftw_workspace.output);

        for (unsigned K_idx = 0; K_idx < refined_count; ++K_idx){
            const unsigned P_idx = coarse_to_fine_idx_map[K_idx];

                val_patch(K_idx) += (fftw_workspace.output[P_idx][0] + I * fftw_workspace.output[P_idx][1]) * weight * Model::FineMomentumGrid().exp_ikr_idx_map[P_idx][r_idx] ;
            }
        }
    }
#ifdef STATIC_CALCULATION
}
#endif
    return val_patch;
}

template <typename Model, typename State>
MatPatch rhs_base_t<Model, State>::eval_Sig_Lam2PI_kMat( const idx_Sig_kMat_t<Model>& idx, const gf_1p_mat_real_t<Model>&  Gvec_real, const gf_1p_mat_t<Model>& Gvec, const double t )
{
    int wext   = idx( ISIGKMAT::w );

    MatPatch val_patch;
    val_patch.setZero(Model::GetRefinedMomentaCount());
    
    for(int R_idx = 0; R_idx < Model::GetFineMomentaCount(); ++R_idx){ 

	dcomplex temp(0.0,0.0); 
	for_freq( int w, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w )
	    for_freq( int wp, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), wp < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++wp ){
		temp += Model::vert_bare_real_space(R_idx, 0, 0, 0, 0) *
		    Gvec_real[w][0][0](R_idx) * 
		    Gvec_real[wp][0][0](R_idx) *
		    Gvec_real[w+wp-wext][0][0](Model::GetNegativeFineMomentumIdx(R_idx)) *
		    FrequencyDependenceScheme<Model>::IntegrationWeightVecs2D()[w][wp] *   //TODO: insert weight
		    Model::vertex_4pt_local_part_bare(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0); //todo: consider non-local
	    }
	Model::FFTW(omp_get_thread_num()).input[R_idx][0] = temp.real();
	Model::FFTW(omp_get_thread_num()).input[R_idx][1] = temp.imag();
    }
    fftw_execute_dft(Model::FFTW(omp_get_thread_num()).forward_plan, Model::FFTW(omp_get_thread_num()).input, Model::FFTW(omp_get_thread_num()).output);

    for (int K_idx = 0; K_idx < Model::GetRefinedMomentaCount(); ++K_idx){
        unsigned P_idx = Model::GetFineMomentumIdxFromCoarse(K_idx);
        
        val_patch(K_idx) = (Model::FFTW(omp_get_thread_num()).output[P_idx][0]+ I * Model::FFTW(omp_get_thread_num()).output[P_idx][1]) * -1./Beta(t)/Beta(t); 
        
        for( int p = 0; p < Model::GetFineMomentaCount(); ++p )
            for_freq( int w, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w ){
		val_patch(K_idx) += - Gvec[w][p](0,0) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w] * Model::vertex_4pt_local_part_bare(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) /(dcomplex)(Model::GetFineMomentaCount())/Beta(t);  // todo: consider non-local 
	    }
    }      
  
    return val_patch;
}



// unused. Todo: remove
template <typename Model, typename State>
MatPatch rhs_base_t<Model, State>::eval_diag_bubble_GG_pp( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_real_t<Model>& Gvec_real) 
{
    int W   = idx( IBUBMAT::W );
    int w   = idx( IBUBMAT::w );
    int m   = idx( IBUBMAT::m );
    int n   = idx( IBUBMAT::n );
    int s4  = idx( IBUBMAT::s1 );
    int s4p = idx( IBUBMAT::s1p);
    int s3  = idx( IBUBMAT::s2 );
    int s3p = idx( IBUBMAT::s2p );
    

    MatPatch val_patch;
    val_patch.setZero(Model::GetRefinedMomentaCount());
    
#ifdef SET_MIXED_BUBBLES_TO_ZERO
    if (m != n)
	return val_patch;
#endif

    for(unsigned r_idx =0; r_idx< Model::RealLattice().ptr->m_points.size(); ++r_idx){
	if(std::abs(Model::FormFactors().fourier_transform_weights[r_idx][m][n]) > CHOP_ERR){    /**< 1E-16 is error on the weights exponents, cosines,...*/
            
	    coord_units_t<Model::dim> r_units = Model::RealLattice().unit_points[r_idx];
	    const int rx = r_units(0);
	    const int ry = r_units(1);

            for(int R_idx = 0; R_idx < Model::GetFineMomentaCount(); R_idx++){        
		const int Rx = GetBubbleFFTRealCoordXFromIdx<Model>(R_idx);
		const int Ry = GetBubbleFFTRealCoordYFromIdx<Model>(R_idx);
		
		dcomplex temp = Gvec_real[ w + div2_ceil( W ) ][s4][s4p](Matidx<Model>(Rx,Ry)) * Gvec_real[ div2_floor( W ) - w - 1 ][s3][s3p](Matidx<Model>(Rx+rx,Ry+ry));
                    
		Model::FFTW(omp_get_thread_num()).input[R_idx][0] = temp.real();
		Model::FFTW(omp_get_thread_num()).input[R_idx][1] = temp.imag();
	    }
	    fftw_execute_dft(Model::FFTW(omp_get_thread_num()).forward_plan, Model::FFTW(omp_get_thread_num()).input, Model::FFTW(omp_get_thread_num()).output);

	    for (int K_idx = 0; K_idx < Model::GetRefinedMomentaCount(); ++K_idx){
		dcomplex weight = Model::FormFactors().fourier_transform_weights[r_idx][m][n];

		unsigned P_idx = Model::GetFineMomentumIdxFromCoarse(K_idx);
		
		weight *= std::conj(Model::FineMomentumGrid().exp_ikr_idx_map[P_idx][r_idx]);
		
		val_patch(K_idx) += (Model::FFTW(omp_get_thread_num()).output[P_idx][0]+ I * Model::FFTW(omp_get_thread_num()).output[P_idx][1]) * weight;  
            }  
        }
    }
  
    return val_patch;
} 

// unused. Todo: remove
template <typename Model, typename State>
MatPatch rhs_base_t<Model, State>::eval_diag_bubble_GG_ph( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_real_t<Model>& Gvec_real )
{
    int W   = idx( IBUBMAT::W );
    int w   = idx( IBUBMAT::w );
    int m   = idx( IBUBMAT::m );
    int n   = idx( IBUBMAT::n );
    int s4  = idx( IBUBMAT::s1 );
    int s4p = idx( IBUBMAT::s1p);
    int s3  = idx( IBUBMAT::s2 );
    int s3p = idx( IBUBMAT::s2p );
   

    MatPatch val_patch;
    val_patch.setZero(Model::GetRefinedMomentaCount());
   
#ifdef SET_MIXED_BUBBLES_TO_ZERO
    if (m != n)
	return val_patch;
#endif 
     
    int w_minus_W_over_2 = w - div2_floor( W );
    int w_plus_W_over_2 = div2_ceil( W ) + w;

#ifdef STATIC_CALCULATION
    for (w = -FrequencyDependenceScheme<Model>::BubbleTU_FermionicGrid().positive_freqs_count; w < FrequencyDependenceScheme<Model>::BubbleTU_FermionicGrid().positive_freqs_count; w++)
{
    w_minus_W_over_2 = w - div2_floor( W );
    w_plus_W_over_2 = div2_ceil( W ) + w;	
#endif    

    for(unsigned r_idx =0; r_idx< Model::RealLattice().ptr->m_points.size(); ++r_idx){
        if(std::abs(Model::FormFactors().fourier_transform_weights[r_idx][m][n]) > CHOP_ERR){    /**< 1E-16 is error on the weights exponents, cosines,...*/

	    coord_units_t<Model::dim> r_units = Model::RealLattice().unit_points[r_idx];
	    const int rx = r_units(0);
	    const int ry = r_units(1);
	    for(int R_idx = 0; R_idx < Model::GetFineMomentaCount(); R_idx++){
		const int Rx = GetBubbleFFTRealCoordXFromIdx<Model>(R_idx);
		const int Ry = GetBubbleFFTRealCoordYFromIdx<Model>(R_idx);
		dcomplex temp = Gvec_real[ w_minus_W_over_2 ][s4][s4p](Matidx<Model>(Model::GetFFTDim()-Rx,Model::GetFFTDim()-Ry)) * Gvec_real[w_plus_W_over_2][s3][s3p](Matidx<Model>(Rx-rx,Ry-ry));
		Model::FFTW(omp_get_thread_num()).input[R_idx][0] = temp.real();
		Model::FFTW(omp_get_thread_num()).input[R_idx][1] = temp.imag();
	    }
	    
	    fftw_execute_dft(Model::FFTW(omp_get_thread_num()).forward_plan, Model::FFTW(omp_get_thread_num()).input, Model::FFTW(omp_get_thread_num()).output);

            for (int K_idx = 0; K_idx < Model::GetRefinedMomentaCount(); ++K_idx){ 
		dcomplex weight = Model::FormFactors().fourier_transform_weights[r_idx][m][n];

		unsigned P_idx = Model::GetFineMomentumIdxFromCoarse(K_idx);

                weight *= Model::FineMomentumGrid().exp_ikr_idx_map[P_idx][r_idx];
            
                val_patch(K_idx) += (Model::FFTW(omp_get_thread_num()).output[P_idx][0]+ I * Model::FFTW(omp_get_thread_num()).output[P_idx][1]) * weight;
            }
        }
    }
#ifdef STATIC_CALCULATION
}
#endif
    return val_patch;
}

template <typename Model, typename State>
void rhs_base_t<Model, State>::print_current_scale(const double t)
{
    std::cout << " Preprocessing at scale " << t;
    
    if (fRGFlowScheme<THE_MODEL>::G_ChosenFlowName() == G_FlowSchemeName::Temperature){
	std::cout << " corresponding to BETA = " << std::exp(- LN_10 * t ); 
    }
    else if (fRGFlowScheme<Model>::G_ChosenFlowName() == G_FlowSchemeName::Interaction)
    {
	std::cout << " corresponding to U = " << Model::vertex_4pt_local_part_bare(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0) * t * t;
    }
    
    std::cout << std::endl;
}
