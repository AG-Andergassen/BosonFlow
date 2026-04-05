#include <mymath.h>
#include <frequencies/matsubara_space.h>
#include <frg/flows.h>
#include <cmath>
#include <complex>
#include <grid.h>
#include <models/square_hubbard.h>


using namespace std; 
using boost::bind;


// screened interaction

template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_w_sc( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    dcomplex P_sc(0.0, 0.0); // polarisation

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq(int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_)
      P_sc += bubble_pp[W][w_][0][mp_][0][0][0][0](K) * state.lambda_sc( W, K, w_, mp_ );// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ). The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta. 
    P_sc /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    dcomplex U_sc = Model::B_sc(W, K, 0, 0, 0, 0);

    val = U_sc/(1.0 - U_sc * P_sc); 

    return val; 
}



template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_w_d( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    dcomplex P_d(0.0, 0.0); // polarisation

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq(int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_)
      P_d += (- bubble_ph[W][w_][0][mp_][0][0][0][0](K)) * state.lambda_d( W, K, w_, mp_ );// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ). The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta
    P_d /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    dcomplex U_d = Model::B_d(W, K, 0, 0, 0, 0);

    val = U_d/(1.0 - U_d * P_d); 

    return val; 
}


template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_w_m( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    dcomplex P_m(0.0, 0.0); // polarisation

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq(int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_)
      P_m += (- bubble_ph[W][w_][0][mp_][0][0][0][0](K)) * state.lambda_m( W, K, w_, mp_ );// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ).  The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta
    P_m /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    dcomplex U_m = Model::B_m(W, K, 0, 0, 0, 0);

    val = U_m/(1.0 - U_m * P_m); 

    return val; 
}



template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_w_sc_BSE( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp, const gf_w_t<Model> &w_sc)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    dcomplex P_sc(0.0, 0.0); // polarisation

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq(int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_)
      P_sc += bubble_pp[W][w_][0][mp_][0][0][0][0](K) * state.lambda_sc( W, K, w_, mp_ );// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ). The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta. 
    P_sc /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    dcomplex U_sc = Model::B_sc(W, K, 0, 0, 0, 0);

    val =  U_sc + w_sc[W][K] * P_sc * U_sc; 

    return val; 
}


template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_w_d_BSE( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, const gf_w_t<Model> &w_d)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    dcomplex P_d(0.0, 0.0); // polarisation

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq(int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_)
      P_d += (- bubble_ph[W][w_][0][mp_][0][0][0][0](K)) * state.lambda_d( W, K, w_, mp_ );// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ). The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta
    P_d /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    dcomplex U_d = Model::B_d(W, K, 0, 0, 0, 0);

    val = U_d + w_d[W][K] * P_d * U_d; 

    return val; 
}


template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_w_m_BSE( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, const gf_w_t<Model> &w_m)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    dcomplex P_m(0.0, 0.0); // polarisation

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq(int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_)
      P_m += (- bubble_ph[W][w_][0][mp_][0][0][0][0](K)) * state.lambda_m( W, K, w_, mp_ );// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ).  The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta
    P_m /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    dcomplex U_m = Model::B_m(W, K, 0, 0, 0, 0);

    val = U_m + w_m[W][K] * P_m * U_m; 

    return val; 
}


// Hedin vertex

template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_lambda_sc( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );
    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
      val += state.B_irreducible_vertex_sc( W, K, w, m, w_, mp_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) * bubble_pp[W][w_][mp_][0][0][0][0][0](K);// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    val /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    if (m == 0)
	val += 1.0;

    return val; 
}


template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_lambda_d( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );
    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
      val += state.B_irreducible_vertex_d( W, K, w, m, w_, mp_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) * (-bubble_ph[W][w_][mp_][0][0][0][0][0](K));// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

	// per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ). The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta
    val /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    if (m == 0)
	val += 1.0;

    return val; 
}


template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_lambda_m( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );
    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() , w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
      val += state.B_irreducible_vertex_m( W, K, w, m, w_, mp_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) * (-bubble_ph[W][w_][mp_][0][0][0][0][0](K));// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

	// per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ). The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta
    val /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    if (m == 0)
	val += 1.0;

    return val; 
}


// SBE rest function M

template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_M_sc( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp)
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ ){
      val += state.B_irreducible_vertex_sc( W, K, w_in, n_in, w_, m_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) *
	       bubble_pp[W][w_][m_][mp_][0][0][0][0](K) *
	(state.B_irreducible_vertex_sc( W, K, w_, mp_, w_out, n_out, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) - state.M_sc( W, K, w_, mp_, w_out, n_out ));//  * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }
    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ) and a matsubara sum with 1/Beta
    val /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    return val; 
}

template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_M_d( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph)
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ ){
      val += state.B_irreducible_vertex_d( W, K, w_in, n_in, w_, m_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) *
	    (-bubble_ph[W][w_][m_][mp_][0][0][0][0](K)) *
	(state.B_irreducible_vertex_d( W, K, w_, mp_, w_out, n_out, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) - state.M_d( W, K, w_, mp_, w_out, n_out ));// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }
    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ) and a matsubara sum with 1/Beta
    val /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    return val; 
}

template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_M_m( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph)
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ ){
      val += state.B_irreducible_vertex_m( W, K, w_in, n_in, w_, m_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) *
	    (-bubble_ph[W][w_][m_][mp_][0][0][0][0](K)) *
	(state.B_irreducible_vertex_m( W, K, w_, mp_, w_out, n_out, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) - state.M_m( W, K, w_, mp_, w_out, n_out ));// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }
    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ) and a matsubara sum with 1/Beta
    val /= Model::MomentumGrid().ptr->get_volume() * BETA;
    
    return val; 
}

// self-energy: SDE SBE MAGNETIC
template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_Sig_sde_sbe_magnetic( const idx_1p_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec)
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );
    int k_in = idx(I1P::k);
    int p_in = Model::GetFineMomentumIdxFromCoarse(k_in);  

    for_freq( int W_, -FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w, W_ < FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w, ++W_ ){
        for( int P_ = 0; P_ < Model::GetFineMomentaCount(); ++P_ ){
	    const unsigned K_ = Model::GetCoarseMomentumIdxFromFine(P_);
	    unsigned p_in_plus_P_ = Model::SumFineMomentaIdxes(p_in, P_);
            int w_plus_W_over_2 = w + div2_floor( W_ );
	    val += Gvec[w + W_][p_in_plus_P_](0,0) * state.mom_lambda_m(W_, K_, w_plus_W_over_2, k_in) * state.w_m(W_, K_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final);
	    
	    //val *= weight_vec[w + W_];
        }
    }
    

    val *= 1.0/BETA/Model::GetFineMomentaCount()/std::sqrt(Model::MomentumGrid().ptr->get_volume());

    return val; 
}

// self-energy: SDE SBE SUPERCONDUCTING
template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_Sig_sde_sbe_superconducting( const idx_1p_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec)
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );
    int k_in = idx(I1P::k);
    int p_in = Model::GetFineMomentumIdxFromCoarse(k_in);
    for_freq( int W_, -FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count+w + 1, W_ < FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count+w, ++W_ ){
        for( int P_ = 0; P_ < Model::GetFineMomentaCount(); ++P_ ){
	    const unsigned K_ = Model::GetCoarseMomentumIdxFromFine(P_);
	    unsigned P_minus_p_in = Model::SumFineMomentaIdxes(P_, Model::GetNegativeFineMomentumIdx(p_in));
	    unsigned K_minus_k_in = Model::GetCoarseMomentumIdxFromFine(P_minus_p_in);
            int W_over_2_minus_w = div2_floor( W_ ) - w -1;
	    val += Gvec[W_ - w - 1][P_minus_p_in](0,0) * state.mom_lambda_sc(W_, K_, W_over_2_minus_w, K_minus_k_in) * state.w_sc(W_, K_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final);
	}
    }

    val *= -1.0/BETA/Model::GetFineMomentaCount()/std::sqrt(Model::MomentumGrid().ptr->get_volume());

    return val; 
}


// self-energy: SDE SBE DENSITY
template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_Sig_sde_sbe_density( const idx_1p_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec)
{
  dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );
    int k_in = idx(I1P::k);
    int p_in = Model::GetFineMomentumIdxFromCoarse(k_in);  


    
    for_freq( int W_, -FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w, W_ < FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w, ++W_ ){
        for( int P_ = 0; P_ < Model::GetFineMomentaCount(); ++P_ ){
	    const unsigned K_ = Model::GetCoarseMomentumIdxFromFine(P_);
	    unsigned p_in_plus_P_ = Model::SumFineMomentaIdxes(p_in, P_);
            int w_plus_W_over_2 = w + div2_floor( W_ );
	    val += Gvec[w + W_][p_in_plus_P_](0,0) * state.mom_lambda_d(W_, K_, w_plus_W_over_2, k_in) * state.w_d(W_, K_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final);

	    const auto bare_U = (Model::B_d(W_, K_, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0) + BOSONIC_INTERACTION_D_SHIFT * Model::MomentumGrid().ptr->get_volume())/std::sqrt(Model::MomentumGrid().ptr->get_volume());
    
	    val -= 2. * bare_U * Gvec[w + W_][p_in_plus_P_](0,0);
        }
    }

    val *= 1.0/BETA/Model::GetFineMomentaCount()/std::sqrt(Model::MomentumGrid().ptr->get_volume());
    return val;
}


#if(false)
template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_Sig_conventional_SDE( const idx_1p_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp, const gf_1p_mat_t<Model>& Gvec)
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );
    int k_in = idx(I1P::k);
    int p_in = Model::GetFineMomentumIdxFromCoarse(k_in);  
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for( int P_ = 0; P_ < Model::GetFineMomentaCount(); ++P_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
      for_freq( int W_, 2*(-FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count+w), W_ < 2*(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count+w-1), ++W_ ){
      unsigned K_ = Model::GetCoarseMomentumIdxFromFine(P_);
      unsigned P_minus_p_in = Model::SumFineMomentaIdxes(Model::GetNegativeFineMomentumIdx(p_in), P_);
      int W_over_2_minus_w = div2_floor( W_ ) - w - 1;

      dcomplex vertex_p_in_m(0, 0);
      dcomplex bare_vertex_mp_p_in(0, 0);
      for (int n_ = 0; n_ < Model::GetMomentumFormFactorsCount(); ++n_ ){
	const dcomplex f_n_p_in = Model::GetFormFactorInMomentumIdxSpace(n_)[p_in];
	const dcomplex f_n_p_in_star = std::conj(Model::GetFormFactorInMomentumIdxSpace(n_)[p_in]);
	vertex_p_in_m += f_n_p_in_star * state.vertex_sc( W_, K_, w, n_, w_, m_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final );
	bare_vertex_mp_p_in += f_n_p_in * state.vertex_sc_bare( W_, K_, w_, mp_, w, n_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final );
      }
      //std::cout << W_ << "," << w << "," << FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count << "," << w_plus_W_over_2 << "," << 2*(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w) <<  std::endl;
	
      dcomplex bubble = 0.0;
      if (std::abs(W_)-1 < FrequencyDependenceScheme<Model>::BubbleTU_BosonicGrid().positive_freqs_count){
	bubble = (bubble_pp[W_][w_][m_][mp_][0][0][0][0](K_));
      }
      
      val += vertex_p_in_m * bubble * bare_vertex_mp_p_in * Gvec[W_over_2_minus_w ][P_minus_p_in](0,0);// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }
    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ) and a matsubara sum with 1/Beta
    val /= -Model::MomentumGrid().ptr->get_volume() * Model::MomentumGrid().ptr->get_volume() * BETA * BETA;
    
    return val; 
}
#endif

#if(false)
template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_Sig_conventional_SDE_magnetic( const idx_1p_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, const gf_1p_mat_t<Model>& Gvec)
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );
    int k_in = idx(I1P::k);
    int p_in = Model::GetFineMomentumIdxFromCoarse(k_in);  
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for( int P_ = 0; P_ < Model::GetFineMomentaCount(); ++P_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
      for_freq( int W_, 2*(-FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w), W_ < 2*(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w-1), ++W_ ){
      unsigned K_ = Model::GetCoarseMomentumIdxFromFine(P_);
      unsigned p_in_plus_P_ = Model::SumFineMomentaIdxes(p_in, P_);
      int w_plus_W_over_2 = w + div2_ceil( W_ );

      dcomplex vertex_p_in_m(0, 0);
      dcomplex bare_vertex_mp_p_in(0, 0);
      for (int n_ = 0; n_ < Model::GetMomentumFormFactorsCount(); ++n_ ){
	const dcomplex f_n_p_in = Model::GetFormFactorInMomentumIdxSpace(n_)[p_in];
	const dcomplex f_n_p_in_star = std::conj(Model::GetFormFactorInMomentumIdxSpace(n_)[p_in]);
	vertex_p_in_m += f_n_p_in_star * state.vertex_m( W_, K_, w, n_, w_, m_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final );
	bare_vertex_mp_p_in += f_n_p_in * state.vertex_m_bare( W_, K_, w_, mp_, w, n_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final );
      }
      //std::cout << W_ << "," << w << "," << FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count << "," << w_plus_W_over_2 << "," << 2*(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w) <<  std::endl;
	
      dcomplex bubble = 0.0;
      if (std::abs(W_)-1 < FrequencyDependenceScheme<Model>::BubbleTU_BosonicGrid().positive_freqs_count){
	bubble = (-bubble_ph[W_][w_][m_][mp_][0][0][0][0](K_));
      }
      
      val += vertex_p_in_m * bubble * bare_vertex_mp_p_in * Gvec[w_plus_W_over_2][p_in_plus_P_](0,0);// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }
    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ) and a matsubara sum with 1/Beta
    val /= Model::MomentumGrid().ptr->get_volume() * Model::MomentumGrid().ptr->get_volume() * BETA * BETA;
    
    return val; 
}
#endif


template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_Sig_conventional_SDE_magnetic( const idx_1p_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, const gf_1p_mat_t<Model>& Gvec)
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );
    int k_in = idx(I1P::k);
    int p_in = Model::GetFineMomentumIdxFromCoarse(k_in);  
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for( int P_ = 0; P_ < Model::GetFineMomentaCount(); ++P_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
      for_freq( int W_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), W_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++W_ ){
      int W_xph = W_ - w;
      int w_xph = w + div2_floor(W_xph);

      unsigned K_ = Model::GetCoarseMomentumIdxFromFine(P_);
      unsigned p_in_plus_P_ = Model::SumFineMomentaIdxes(p_in, P_);
      //int w_plus_W_over_2 = w + div2_ceil( W_ );

      dcomplex vertex_p_in_m(0, 0);
      dcomplex bare_vertex_mp_p_in(0, 0);
      for (int n_ = 0; n_ < Model::GetMomentumFormFactorsCount(); ++n_ ){
	const dcomplex f_n_p_in = Model::GetFormFactorInMomentumIdxSpace(n_)[p_in];
	const dcomplex f_n_p_in_star = std::conj(Model::GetFormFactorInMomentumIdxSpace(n_)[p_in]);
	vertex_p_in_m += f_n_p_in_star * state.vertex_m( W_xph, K_, w_xph, n_, w_, m_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final );
	bare_vertex_mp_p_in += f_n_p_in * (2.0*state.vertex_m_bare( W_xph, K_, w_, mp_, w_xph, n_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final) - state.vertex_m_bare( W_xph, K_, W_xph, mp_, 0, n_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final ) );
	//bare_vertex_mp_p_in += f_n_p_in * state.vertex_m_bare( W_xph, K_, w_, mp_, w_xph, n_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final );
      }
      //std::cout << W_ << "," << w << "," << FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count << "," << w_plus_W_over_2 << "," << 2*(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w) <<  std::endl;
	
      dcomplex bubble = 0.0;
      if (std::abs(W_)-1 < FrequencyDependenceScheme<Model>::BubbleTU_BosonicGrid().positive_freqs_count){
	bubble = (-bubble_ph[W_xph][w_][m_][mp_][0][0][0][0](K_));
      }
      
      val += vertex_p_in_m * bubble * bare_vertex_mp_p_in * Gvec[W_][p_in_plus_P_](0,0);// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }
    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ) and a matsubara sum with 1/Beta

    val /= Model::MomentumGrid().ptr->get_volume() * Model::MomentumGrid().ptr->get_volume() * BETA * BETA;
   
    dcomplex val_Fock( 0.0, 0.0 );  
    // caution: enable flag SET_SIG_TAIL_TO_ZERO
    for( int p_ = 0; p_ < Model::GetFineMomentaCount(); ++p_ )
    for_freq( int W_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), W_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++W_ )
    {     
	val_Fock += Gvec[W_][p_](0,0) * //FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[W_] *   /**< Minus sign from internal loop */
	    state.ferm_vertex_sc_bare( 00, 00, w, 0, W_, 0, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final );
    }      
    
    val_Fock *= 1.0/BETA/Model::GetFineMomentaCount();


 
    return val + val_Fock; 
}

/*
// pairing
template<typename Model, typename state_t>
dcomplex rhs_sbe_selfconsistent_t<Model, state_t>::eval_Sig_conventional_SDE( const idx_1p_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, const gf_1p_mat_t<Model>& Gvec)
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );
    int k_in = idx(I1P::k);
    int p_in = Model::GetFineMomentumIdxFromCoarse(k_in);  
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for( int P_ = 0; P_ < Model::GetFineMomentaCount(); ++P_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
      for_freq( int W_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), W_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++W_ ){
      int W_pp = W_ + w + 1;
      int w_pp = w + div2_ceil(W_pp);

      unsigned K_ = Model::GetCoarseMomentumIdxFromFine(P_);
      unsigned P_minus_p_in = Model::SumFineMomentaIdxes(Model::GetNegativeFineMomentumIdx(p_in), P_);
      
      dcomplex vertex_p_in_m(0, 0);
      dcomplex bare_vertex_mp_p_in(0, 0);
      for (int n_ = 0; n_ < Model::GetMomentumFormFactorsCount(); ++n_ ){
	const dcomplex f_n_p_in = Model::GetFormFactorInMomentumIdxSpace(n_)[p_in];
	const dcomplex f_n_p_in_star = std::conj(Model::GetFormFactorInMomentumIdxSpace(n_)[p_in]);
	vertex_p_in_m += f_n_p_in_star * state.vertex_m( W_pp, K_, w_pp, n_, w_, m_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final );
	bare_vertex_mp_p_in += f_n_p_in * state.vertex_m_bare( W_pp, K_, w_, mp_, w_pp, n_, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final );
      }
      //std::cout << W_ << "," << w << "," << FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count << "," << w_plus_W_over_2 << "," << 2*(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w) <<  std::endl;
	
      dcomplex bubble = 0.0;
      if (std::abs(W_)-1 < FrequencyDependenceScheme<Model>::BubbleTU_BosonicGrid().positive_freqs_count){
	bubble = bubble = (bubble_pp[W_][w_][m_][mp_][0][0][0][0](K_));
      }
      
      val += vertex_p_in_m * bubble * bare_vertex_mp_p_in * Gvec[W_][p_in_plus_P_](0,0);// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }
    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ) and a matsubara sum with 1/Beta
    val /= Model::MomentumGrid().ptr->get_volume() * Model::MomentumGrid().ptr->get_volume() * BETA * BETA;
    
    return val; 
}
*/
