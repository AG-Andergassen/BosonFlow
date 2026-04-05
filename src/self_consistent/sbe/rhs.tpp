#include <mymath.h>
#include <cmath>
#include <complex>
#include <frg/sbe/symmetries.h>
#include <frg/flows.h>


using namespace std; 
using boost::bind;

template <typename Model, typename state_t >
rhs_sbe_selfconsistent_t<Model, state_t>::rhs_sbe_selfconsistent_t()
    :m_Gvec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true )  
{ 
    
}


template <typename Model, typename state_t >
rhs_sbe_selfconsistent_t<Model, state_t>::~rhs_sbe_selfconsistent_t()
{

}

template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::vertex( const state_t& old_state, state_t& new_state, const double prev_norm )
{

#ifdef DEBUG_BUILD
    omp_set_num_threads(1);  	
#endif
    
    #pragma omp parallel default(shared)
    {

#pragma omp single
	cout << " ... w (screened interaction) " << endl;
	    
	rhs_w_sc(&new_state.gf_w_sc(), old_state, m_bubble_pp, prev_norm);
	rhs_w_d(&new_state.gf_w_d(), old_state, m_bubble_ph, prev_norm);
	rhs_w_m(&new_state.gf_w_m(), old_state, m_bubble_ph, prev_norm);
	
#ifdef NO_HEDIN_VERTEX_FLOW
#pragma omp single
	cout << " ... lambda (Hedin vertex)" << endl;

    rhs_lambda_sc(&new_state.gf_lambda_sc(), old_state, m_bubble_pp); 
	rhs_lambda_d(&new_state.gf_lambda_d(), old_state, m_bubble_ph);
	rhs_lambda_m(&new_state.gf_lambda_m(), old_state, m_bubble_ph);
#endif
	
#if !defined(SBEa_APPROXIMATION) && !defined(SBEb_APPROXIMATION)

#pragma omp single
	cout << " ... M (Rest functions)" << endl;

	rhs_M_sc(&new_state.gf_M_sc(), old_state, m_bubble_pp);
	rhs_M_d(&new_state.gf_M_d(), old_state, m_bubble_ph);
	rhs_M_m(&new_state.gf_M_m(), old_state, m_bubble_ph);

#endif
	
    }

    auto &the_state = old_state;

    coord_t<Model::dim> zero;
    unsigned p_zeroidx = std::get<0>(Model::MomentumGrid().ptr->get_idx_from_coord(zero));
    std::cout << "gf_w_sc at (W=0, K=0):" << the_state.gf_w_sc()[0][p_zeroidx] << std::endl;
    std::cout << "gf_lambda_sc at (W=0, K=0, w=0, m=0):" << the_state.gf_lambda_sc()[0][p_zeroidx][0][0] << std::endl;
    std::cout << "gf_M_sc at (W=0, K=0, w=0, m=0, wp=0, mp=0):" << the_state.gf_M_sc()[0][p_zeroidx][0][0][0][0] << std::endl;
    std::cout << "gf_Sig at (w=0, k=0, s_in=0, s_out=0):" << the_state.gf_Sig()[0][p_zeroidx][0][0] << std::endl;
}


template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::update_bubbles_and_G( const state_t& state, const double t){

#ifdef DEBUG_BUILD
    omp_set_num_threads(1);	
#endif

    rhs_base_t::compute_Gvec(&m_Gvec, state, t);

    gf_1p_mat_real_t<Model> Gvec_realspace(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true);   

    cout << " ... GG bubbles " << endl;
    
    #pragma omp parallel default(shared) 
    {
	rhs_base_t::compute_Gvec_real(&Gvec_realspace, m_Gvec);
#pragma omp barrier //barrier to make sure Gvec_realspace is present
	rhs_base_t::compute_GS_bubbles(&m_bubble_pp, &m_bubble_ph, Gvec_realspace, Gvec_realspace);
	//rhs_base_t::compute_GG_bubbles(&bubble_pp, &bubble_ph, Gvec_realspace);
    }

    // TODO: after adding factor of 1/2 in the original definition, remove this
    m_bubble_pp.init([this](const idx_bubble_mat_t<Model> &idx){return 0.5 * this->m_bubble_pp(idx);});
    m_bubble_ph.init([this](const idx_bubble_mat_t<Model> &idx){return 0.5 * this->m_bubble_ph(idx);});
}  


template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::selfenergy(const state_t &old_state, gf_1p_t<Model> &new_Sig)
{   
    cout << " ... self-energy " << endl;
#ifndef STATIC_CALCULATION 
    rhs_Sig(&new_Sig, old_state, m_bubble_ph, m_Gvec);
#else
    new_Sig.init( [](const idx_1p_t<Model> &idx){return 0;} );
#endif	
}


template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_w_sc(gf_w_t<Model> *w_sc_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_pp, const double prev_norm)
{
     
#ifdef ITERATE_DYSON_EQUATION_FOR_BOSONIC_SC_PROPAGATOR
    gf_w_t<Model> w_sc_prev = old_state.gf_w_sc();
    int iterations_count = static_cast<int>(std::max(1.0, 1.0/prev_norm));
    std::cout << "... Iterations in bosonic dyson equation sc " << iterations_count << std::endl;
    for (unsigned i = 0; i < iterations_count; i++){
        Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_sc_Ptr()->init_batched( (*w_sc_ptr), [&old_state, &bubble_pp, &w_sc_prev]( const idx_w_t<Model>& idx ){ return eval_w_sc_BSE( idx, old_state, bubble_pp, w_sc_prev); } );
	#pragma omp barrier
	w_sc_prev = *w_sc_ptr;
    }
#else
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_sc_Ptr()->init_batched( (*w_sc_ptr), [&old_state, &bubble_pp]( const idx_w_t<Model>& idx ){ return eval_w_sc( idx, old_state, bubble_pp ); } );
#endif

}

template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_w_d(gf_w_t<Model> *w_d_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph, const double prev_norm)
{
#ifdef ITERATE_DYSON_EQUATION_FOR_BOSONIC_D_PROPAGATOR
    gf_w_t<Model> w_d_prev = old_state.gf_w_d();
    int iterations_count = static_cast<int>(std::max(1.0, 1.0/prev_norm));
    std::cout << "... Iterations in dyson equation sc " << iterations_count << std::endl;
    for (unsigned i = 0; i < iterations_count; i++){
        Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_d_Ptr()->init_batched( (*w_d_ptr), [&old_state, &bubble_ph, &w_d_prev]( const idx_w_t<Model>& idx ){ return eval_w_d_BSE( idx, old_state, bubble_ph, w_d_prev); } );
	#pragma omp barrier
	w_d_prev = *w_d_ptr;
    }
#else
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_d_Ptr()->init_batched( (*w_d_ptr), [&old_state, &bubble_ph]( const idx_w_t<Model>& idx ){ return eval_w_d( idx, old_state, bubble_ph ); } ); 
#endif
}

template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_w_m(gf_w_t<Model> *w_m_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph, const double prev_norm )
{
#ifdef ITERATE_DYSON_EQUATION_FOR_BOSONIC_M_PROPAGATOR
    gf_w_t<Model> w_m_prev = old_state.gf_w_m();
    int iterations_count = static_cast<int>(std::max(1.0, 1.0/prev_norm));
    std::cout << "... Iterations in dyson equation sc " << iterations_count << std::endl;
    for (unsigned i = 0; i < iterations_count; i++){
        Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_m_Ptr()->init_batched( (*w_m_ptr), [&old_state, &bubble_ph, &w_m_prev]( const idx_w_t<Model>& idx ){ return eval_w_m_BSE( idx, old_state, bubble_ph, w_m_prev); } );
	#pragma omp barrier
	w_m_prev = *w_m_ptr;
    }    
#else
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_m_Ptr()->init_batched( (*w_m_ptr), [&old_state, &bubble_ph]( const idx_w_t<Model>& idx ){ return eval_w_m( idx, old_state, bubble_ph ); } );    
#endif

}

template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_lambda_sc(gf_lambda_t<Model> *lambda_sc_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_pp)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_sc_Ptr()->init_batched( (*lambda_sc_ptr), [&old_state, &bubble_pp]( const idx_lambda_t<Model>& idx ){ return eval_lambda_sc(idx, old_state, bubble_pp); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_lambda_d(gf_lambda_t<Model> *lambda_d_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_d_Ptr()->init_batched( (*lambda_d_ptr), [&old_state, &bubble_ph]( const idx_lambda_t<Model>& idx ){ return eval_lambda_d( idx, old_state, bubble_ph); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_lambda_m(gf_lambda_t<Model> *lambda_m_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_m_Ptr()->init_batched( (*lambda_m_ptr), [&old_state, &bubble_ph]( const idx_lambda_t<Model>& idx ){ return eval_lambda_m( idx, old_state, bubble_ph ); } );
}


template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_M_sc(gf_M_t<Model> *M_sc_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_pp)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched( (*M_sc_ptr), [&old_state, &bubble_pp]( const idx_M_t<Model>& idx ){ return eval_M_sc(idx, old_state, bubble_pp ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_M_d(gf_M_t<Model> *M_d_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched( (*M_d_ptr), [&old_state, &bubble_ph]( const idx_M_t<Model>& idx ){ return eval_M_d(idx, old_state, bubble_ph ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_M_m(gf_M_t<Model> *M_m_ptr, const state_t &old_state, const gf_bubble_mat_t<Model> &bubble_ph)
{
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched( (*M_m_ptr), [&old_state, &bubble_ph]( const idx_M_t<Model>& idx ){ return eval_M_m(idx, old_state, bubble_ph ); } ); 
}

template <typename Model, typename state_t >
void rhs_sbe_selfconsistent_t<Model, state_t>::rhs_Sig(gf_1p_t<Model> *Sig_ptr, const state_t &old_state, const gf_bubble_mat_t<Model>& bubble_ph, const gf_1p_mat_t<Model>& Gvec)
{
    //todo: check self-energies are called appropriately
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_sig_Ptr()->init_batched( (*Sig_ptr), [&old_state, &Gvec, &bubble_ph]( const idx_1p_t<Model>& idx ){
#ifdef SELFEN_SDE_CONVENTIONAL_MAGNETIC
      return eval_Sig_conventional_SDE_magnetic( idx, old_state, bubble_ph, Gvec);
#elif defined(SELFEN_SDE_SBE_MAGNETIC)
      return eval_Sig_sde_sbe_magnetic( idx, old_state, Gvec);
#elif defined(SELFEN_SDE_SBE_SUPERCONDUCTING)
      return eval_Sig_sde_sbe_superconducting( idx, old_state, Gvec);
#elif defined(SELFEN_SDE_SBE_DENSITY)
      return eval_Sig_sde_sbe_density( idx, old_state, Gvec);
#else
      // Default to magnetic when no specific SDE flag is set
      return eval_Sig_sde_sbe_magnetic( idx, old_state, Gvec);
#endif
    } );
}
