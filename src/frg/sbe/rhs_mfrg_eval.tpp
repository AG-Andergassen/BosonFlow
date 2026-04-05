#include <mymath.h> //should be added from rhs_1lfrg.cpp, right?
#include <cmath> //should be added from rhs_1lfrg.cpp, right?
#include <complex> //should be added from rhs_1lfrg.cpp, right?
#include <frg/sbe/symmetries.h> //should be added from rhs_1lfrg.cpp, right?
#include <frg/flows.h> //should be added from rhs_1lfrg.cpp, right?
#include <models/square_hubbard.h>


  // ------- Left corrections for vertex corrections -----------
  //The calculation of left corrections (and their storage) enables us to avoid nested frequency loops in the eval functions for the 3-loop (or higher) vertex corrections

  // lambda (Hedin vertex)

template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_sc_left_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif

    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
	for_freq( int wp_, negative_frequency_integration_range, wp_ < positive_frequency_integration_range, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	    {
		int mpp_ = mp_;
#else
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
		temp += state.lambda_sc( W, K, wp_, mpp_ ) *
		    bubble_GG_pp[W][wp_][mpp_][mp_][0][0][0][0](K);
	    }
	    val += temp * dstate_over_dt_l_minus_1.phi_sc_bar_dot(state, W, K, wp_, mp_, w, m, t );
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_pp( W, t ) *
	    state.lambda_sc( W, K, wp_, mp_ ) *
	    dstate_over_dt_l_minus_1.phi_sc_bar_dot(state, W, K, wp_, mp_, w, m, t);
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_d_left_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif

    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
	for_freq( int wp_, negative_frequency_integration_range, wp_ < positive_frequency_integration_range, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	    {
		int mpp_ = mp_;
#else
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
		temp += state.lambda_d( W, K, wp_, mpp_ ) *
		    bubble_GG_ph[W][wp_][mpp_][mp_][0][0][0][0](K);
	    }
	    val += temp * dstate_over_dt_l_minus_1.phi_d_bar_dot(state, W, K, wp_, mp_, w, m, t );
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) *
	    state.lambda_d( W, K, wp_, mp_ ) *
	    dstate_over_dt_l_minus_1.phi_d_bar_dot(state, W, K, wp_, mp_, w, m, t );
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_m_left_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif
    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
	for_freq( int wp_, negative_frequency_integration_range, wp_ < positive_frequency_integration_range, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	    {
		int mpp_ = mp_;
#else
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
		temp += state.lambda_m( W, K, wp_, mpp_ ) *
		    bubble_GG_ph[W][wp_][mpp_][mp_][0][0][0][0](K);
	    }
	    val += temp * dstate_over_dt_l_minus_1.phi_m_bar_dot(state, W, K, wp_, mp_, w, m, t );
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) *
	    state.lambda_m( W, K, wp_, mp_ ) *
	    dstate_over_dt_l_minus_1.phi_m_bar_dot(state, W, K, wp_, mp_, w, m, t );
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


  // M (rest function)

#if !defined(SBEa_APPROXIMATION)
template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_sc_left_corrections( const idx_M_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
        for_freq( int w_, negative_frequency_integration_range, w_ < positive_frequency_integration_range, ++w_ ){
       	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
            for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
                temp += state.B_irreducible_vertex_sc( W, K, w_in, n_in, w_, mp_, t ) * bubble_GG_pp[W][w_][mp_][m_][0][0][0][0](K) ;   
            }
	    val += temp * dstate_over_dt_l_minus_1.phi_sc_bar_dot(state, W, K, w_, m_, w_out, n_out, t );
        }

#ifndef STATIC_CALCULATION         
	int w_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

        val += inv_freq_sum_normalisation(t) * state.B_irreducible_vertex_sc( W, K, w_in, n_in, w_, m_, t ) * fRGFlowScheme<Model>::asymptotic_GG_pp( W, t ) * dstate_over_dt_l_minus_1.phi_sc_bar_dot(state, W, K, w_, m_, w_out, n_out, t );
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_d_left_corrections( const idx_M_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
        for_freq( int w_, negative_frequency_integration_range, w_ < positive_frequency_integration_range, ++w_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
            for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
                temp += state.B_irreducible_vertex_d( W, K, w_in, n_in, w_, mp_, t ) * bubble_GG_ph[W][w_][mp_][m_][0][0][0][0](K) ;   
            }
	    val += temp * dstate_over_dt_l_minus_1.phi_d_bar_dot(state, W, K, w_, m_, w_out, n_out, t );
        }

#ifndef STATIC_CALCULATION         
	int w_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

        val += inv_freq_sum_normalisation(t) * state.B_irreducible_vertex_d( W, K, w_in, n_in, w_, m_, t ) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) * dstate_over_dt_l_minus_1.phi_d_bar_dot(state, W, K, w_, m_, w_out, n_out , t);
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_m_left_corrections( const idx_M_t<Model>& idx, const state_t& state, const state_t &dstate_over_dt_l_minus_1, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

    
#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
      for_freq( int w_, negative_frequency_integration_range, w_ < positive_frequency_integration_range, ++w_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
            for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
                temp += state.B_irreducible_vertex_m( W, K, w_in, n_in, w_, mp_, t ) * bubble_GG_ph[W][w_][mp_][m_][0][0][0][0](K) ;   
            }
	    val += temp * dstate_over_dt_l_minus_1.phi_m_bar_dot(state, W, K, w_, m_, w_out, n_out, t );
        }

#ifndef STATIC_CALCULATION         
	int w_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

        val += inv_freq_sum_normalisation(t) * state.B_irreducible_vertex_m( W, K, w_in, n_in, w_, m_, t ) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) * dstate_over_dt_l_minus_1.phi_m_bar_dot(state, W, K, w_, m_, w_out, n_out, t );
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}
#endif //!defined(SBEa_APPROXIMATION)



  // ------- Central corrections for vertex corrections -----------
  //The calculation of the central corrections is done below to be used subsequently for the eval functions of \dot{\lambda}^(l), \dot{w}^(l) and \dot{M}^(l) at loop orders l>2

  // lambda (Hedin vertex)

template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_sc_central_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
	        temp += bubble_GG_pp[W][wp_][mp_][mpp_][0][0][0][0](K) * state.B_irreducible_vertex_sc( W, K, wp_, mpp_, w, m, t );
	        }
	    val += m_lambda_sc_dot_left_corrections[W][K][wp_][mp_] * temp;
            }
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_d_central_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
        for_freq( int wp_,  -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
	        temp += bubble_GG_ph[W][wp_][mp_][mpp_][0][0][0][0](K) * state.B_irreducible_vertex_d( W, K, wp_, mpp_, w, m, t );
	        }
	        val += m_lambda_d_dot_left_corrections[W][K][wp_][mp_] * temp;
            }
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_m_central_corrections( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
	        temp += bubble_GG_ph[W][wp_][mp_][mpp_][0][0][0][0](K) * state.B_irreducible_vertex_m( W, K, wp_, mpp_, w, m, t );
	        }
	    val += m_lambda_m_dot_left_corrections[W][K][wp_][mp_] * temp;
            }
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


  // w (screened interaction)

template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_w_sc_central_corrections( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    // m_lambda_sc_dot_left_corrections "decays" (is zero) outside this range for W
    if (FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count <= W || W <= -FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count)
	return val;

    dcomplex w_sc = state.w_sc( W, K, t );
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
	        temp += bubble_GG_pp[W][wp_][mp_][mpp_][0][0][0][0](K) * state.lambda_sc( W, K, wp_, mpp_ );
	        }
	    val += m_lambda_sc_dot_left_corrections[W][K][wp_][mp_] * temp;
            }
    }

    val *= w_sc * w_sc/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_w_d_central_corrections( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    // m_lambda_d_dot_left_corrections "decays" (is zero) outside this range for W
    if (FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count <= W || W <= -FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count)
	return val;

    dcomplex w_d = state.w_d( W, K, t );
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
	        temp += bubble_GG_ph[W][wp_][mp_][mpp_][0][0][0][0](K) * state.lambda_d( W, K, wp_, mpp_ );
	        }
	    val += m_lambda_d_dot_left_corrections[W][K][wp_][mp_] * temp;
            }
    }

    val *= w_d * w_d/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_w_m_central_corrections( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    // m_lambda_m_dot_left_corrections "decays" (is zero) outside this range for W
    if (FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count <= W || W <= -FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count)
	return val;

    dcomplex w_m = state.w_m( W, K, t );
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
	        temp += bubble_GG_ph[W][wp_][mp_][mpp_][0][0][0][0](K) * state.lambda_m( W, K, wp_, mpp_ );
	        }
	    val += m_lambda_m_dot_left_corrections[W][K][wp_][mp_] * temp;
            }
    }

    val *= w_m * w_m/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


  // M (rest function)

#if !defined(SBEa_APPROXIMATION)
template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_sc_central_corrections( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
	        temp += bubble_GG_pp[W][wp_][mp_][mpp_][0][0][0][0](K) * state.B_irreducible_vertex_sc( W, K, wp_, mpp_, w_out, n_out, t );
	        }
	    val += m_M_sc_dot_left_corrections[W][K][w_in][n_in][wp_][mp_] * temp;
            }
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_d_central_corrections( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
	        temp += bubble_GG_ph[W][wp_][mp_][mpp_][0][0][0][0](K) * state.B_irreducible_vertex_d( W, K, wp_, mpp_, w_out, n_out, t );
	        }
	    val += m_M_d_dot_left_corrections[W][K][w_in][n_in][wp_][mp_] * temp;
            }
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_m_central_corrections( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
	        temp += bubble_GG_ph[W][wp_][mp_][mpp_][0][0][0][0](K) * state.B_irreducible_vertex_m( W, K, wp_, mpp_, w_out, n_out, t );
	    }
	    val += m_M_m_dot_left_corrections[W][K][w_in][n_in][wp_][mp_] * temp;
	}
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}
#endif //!defined(SBEa_APPROXIMATION)



  // ------- 2-loop vertex corrections -----------

  // lambda (Hedin vertex)

template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_sc_2l_corrections( const idx_lambda_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );


    val = m_lambda_sc_dot_left_corrections[W][K][w][m];

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_d_2l_corrections( const idx_lambda_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

   
    val = -m_lambda_d_dot_left_corrections[W][K][w][m];

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_m_2l_corrections( const idx_lambda_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

  
    val = -m_lambda_m_dot_left_corrections[W][K][w][m];
 
    return val;
}


  // M (rest function)

#if !defined(SBEa_APPROXIMATION)
template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_sc_2l_corrections( const idx_M_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

 
    val = m_M_sc_dot_left_corrections[W][K][w_in][n_in][w_out][n_out]+m_M_sc_dot_left_corrections[W][K][w_out][n_out][w_in][n_in];

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_d_2l_corrections( const idx_M_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );


    val = -(m_M_d_dot_left_corrections[W][K][w_in][n_in][w_out][n_out]+m_M_d_dot_left_corrections[W][K][w_out][n_out][w_in][n_in]);

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_m_2l_corrections( const idx_M_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );


    val = -(m_M_m_dot_left_corrections[W][K][w_in][n_in][w_out][n_out]+m_M_m_dot_left_corrections[W][K][w_out][n_out][w_in][n_in]);

    return val;
}
#endif //!defined(SBEa_APPROXIMATION)



  // ------- Higher loop (3 or more) vertex corrections-----------

  // lambda (Hedin vertex)

template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_sc_higher_loop_corrections( const idx_lambda_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    val = m_lambda_sc_dot_left_corrections[W][K][w][m]+m_lambda_sc_dot_central_corrections[W][K][w][m];
 
    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_d_higher_loop_corrections( const idx_lambda_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    val = -(m_lambda_d_dot_left_corrections[W][K][w][m]-m_lambda_d_dot_central_corrections[W][K][w][m]);

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_lambda_m_higher_loop_corrections( const idx_lambda_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    val = -(m_lambda_m_dot_left_corrections[W][K][w][m]-m_lambda_m_dot_central_corrections[W][K][w][m]);

    return val;
}


  // w (screened interaction)

template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_w_sc_higher_loop_corrections( const idx_w_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    val = m_w_sc_dot_central_corrections[W][K];

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_w_d_higher_loop_corrections( const idx_w_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    val = m_w_d_dot_central_corrections[W][K];
 
    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_w_m_higher_loop_corrections( const idx_w_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    val = m_w_m_dot_central_corrections[W][K];

    return val;
}


  // M (rest function)

#if !defined(SBEa_APPROXIMATION)
template <typename Model, typename state_t >
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_sc_higher_loop_corrections( const idx_M_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

    val = m_M_sc_dot_left_corrections[W][K][w_in][n_in][w_out][n_out]+m_M_sc_dot_central_corrections[W][K][w_in][n_in][w_out][n_out]+m_M_sc_dot_left_corrections[W][K][w_out][n_out][w_in][n_in];

    return val;
}


template <typename Model, typename state_t >
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_d_higher_loop_corrections( const idx_M_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

    val = -(m_M_d_dot_left_corrections[W][K][w_in][n_in][w_out][n_out]-m_M_d_dot_central_corrections[W][K][w_in][n_in][w_out][n_out]+m_M_d_dot_left_corrections[W][K][w_out][n_out][w_in][n_in]);

    return val;
}


template <typename Model, typename state_t >
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_M_m_higher_loop_corrections( const idx_M_t<Model>& idx, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );

    val = -(m_M_m_dot_left_corrections[W][K][w_in][n_in][w_out][n_out]-m_M_m_dot_central_corrections[W][K][w_in][n_in][w_out][n_out]+m_M_m_dot_left_corrections[W][K][w_out][n_out][w_in][n_in]);

    return val;
}
#endif //!defined(SBEa_APPROXIMATION)


// note: does not mention orbital indices (compared to the parquet version)...
template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_Sig_SDE_SBE_Magnetic( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const double t )
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
	    val += Gdotvec[w + W_][p_in_plus_P_](0,0) * state.mom_lambda_m(W_, K_, w_plus_W_over_2, k_in) * state.w_m(W_, K_, t);
	    val += Gvec[w + W_][p_in_plus_P_](0,0) * dstate_over_dt.mom_lambda_m_dot(W_, K_, w_plus_W_over_2, k_in) * state.w_m(W_, K_, t);
	    val += Gvec[w + W_][p_in_plus_P_](0,0) * state.mom_lambda_m(W_, K_, w_plus_W_over_2, k_in) * dstate_over_dt.w_m_dot(W_, K_, t);
	    
	    //val *= weight_vec[w + W_];
        }
    }

    return val *= 1.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount()/std::sqrt(Model::MomentumGrid().ptr->get_volume());
 
}


// note: does not mention orbital indices (compared to the parquet version)...
template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_Sig_SDE_SBE_Density( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const double t )
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
	    val += Gdotvec[w + W_][p_in_plus_P_](0,0) * state.mom_lambda_d(W_, K_, w_plus_W_over_2, k_in) * state.w_d(W_, K_, t);
	    val += Gvec[w + W_][p_in_plus_P_](0,0) * dstate_over_dt.mom_lambda_d_dot(W_, K_, w_plus_W_over_2, k_in) * state.w_d(W_, K_, t);
	    val += Gvec[w + W_][p_in_plus_P_](0,0) * state.mom_lambda_d(W_, K_, w_plus_W_over_2, k_in) * dstate_over_dt.w_d_dot(W_, K_, t);
	    const auto bare_U = (Model::B_d(W_, K_, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0) + BOSONIC_INTERACTION_D_SHIFT * Model::MomentumGrid().ptr->get_volume())/std::sqrt(Model::MomentumGrid().ptr->get_volume());
	    const auto bare_U_dot = bare_U * fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t);
    
	    val -= 2. * bare_U * Gdotvec[w + W_][p_in_plus_P_](0,0);
            val -= 2. * bare_U_dot * Gvec[w + W_][p_in_plus_P_](0,0);
	    //val *= weight_vec[w + W_];
        }
    }

    return val *= 1.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount()/std::sqrt(Model::MomentumGrid().ptr->get_volume());
 
}


// note: does not mention orbital indices (compared to the parquet version)...
template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_Sig_SDE_SBE_Superconducting( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const double t )
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
            int W_over_2_minus_w = div2_floor( W_ ) - w - 1;
	    val += Gdotvec[W_ - w - 1][P_minus_p_in](0,0) * state.mom_lambda_sc(W_, K_, W_over_2_minus_w, K_minus_k_in) * state.w_sc(W_, K_, t);
	    val += Gvec[W_ - w - 1][P_minus_p_in](0,0) * dstate_over_dt.mom_lambda_sc_dot(W_, K_, W_over_2_minus_w, K_minus_k_in) * state.w_sc(W_, K_, t);
	    val += Gvec[W_ - w - 1][P_minus_p_in](0,0) * state.mom_lambda_sc(W_, K_, W_over_2_minus_w, K_minus_k_in) * dstate_over_dt.w_sc_dot(W_, K_, t);
            //val *= weight_vec[w + W_];
        }
    }

    return val *= -1.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount()/std::sqrt(Model::MomentumGrid().ptr->get_volume());
 
}



// note: does not mention orbital indices (compared to the parquet version)...
template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_Sig_SDE_Magnetic_Conventional( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model> &bubble_ph, const gf_bubble_mat_t<Model> &bubble_ph_dot,  const double t )
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

      dcomplex vertex_p_in_m(0, 0), vertex_dot_p_in_m(0, 0);
      dcomplex bare_vertex_mp_p_in(0, 0), bare_vertex_dot_mp_p_in(0, 0);
      
      //const auto bare_Um = (Model::B_m(W_xph, K_, 0, 0, 0, 0) + BOSONIC_INTERACTION_M_SHIFT * Model::MomentumGrid().ptr->get_volume());
      //const auto bare_Um_dot = bare_Um * fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t);


      for (int n_ = 0; n_ < Model::GetMomentumFormFactorsCount(); ++n_ ){
	const dcomplex f_n_p_in = Model::GetFormFactorInMomentumIdxSpace(n_)[p_in];
	const dcomplex f_n_p_in_star = std::conj(Model::GetFormFactorInMomentumIdxSpace(n_)[p_in]);
	vertex_p_in_m += f_n_p_in_star * state.vertx_m( W_xph, K_, w_xph, n_, w_, m_, t );

	bare_vertex_mp_p_in += f_n_p_in * (2.0*state.vertex_m_bare( W_xph, K_, w_, mp_, w_xph, n_, t ) - state.vertex_m_bare( W_xph, K_, W_xph, mp_, 0, n_, t ) );
	// warning: only works for density-density interactions in frequency space; should be reviewed for other cases
	const auto B_dot1 = state.bos_vertex_m_bare_dot( W_xph  /*can be whatever*/ , K_, w_, mp_, w_xph, n_, t );
	const auto F_dot1 = state.ferm_vertex_m_bare_dot( W_xph  /*can be whatever*/ , K_, w_, mp_, w_xph, n_, t );
	const auto B_dot2 = state.bos_vertex_m_bare_dot( W_xph /*can be whatever*/, K_, W_xph, mp_, 0, n_, t );
	const auto F_dot2 = state.ferm_vertex_m_bare_dot( W_xph /*can be whatever*/, K_, W_xph, mp_, 0, n_, t );


	//B^M(q, k, k') + F(q, k, k') = U(k - k') ~  1/(\omega_0^2 + (k-k')^2)

	vertex_dot_p_in_m += f_n_p_in_star * (dstate_over_dt.phi_m_dot(state, W_xph, K_, w_xph, n_, w_, m_, t ) + dstate_over_dt.phi_m_bar_dot(state, W_xph, K_, w_xph, n_, w_, m_, t ) + state.bos_vertex_m_bare_dot( W_xph , K_, w_xph, n_, w_, m_, t ) + state.ferm_vertex_m_bare_dot( W_xph, K_, w_xph, n_, w_, m_, t ) );
	bare_vertex_dot_mp_p_in += f_n_p_in * (2.0*(B_dot1 + F_dot1) - (B_dot2 + F_dot2));
      }
      //std::cout << W_ << "," << w << "," << FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count << "," << w_plus_W_over_2 << "," << 2*(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w) <<  std::endl;
	
      dcomplex bubble = 0.0;
      dcomplex bubble_dot = 0.0;
      if (std::abs(W_)-1 < FrequencyDependenceScheme<Model>::BubbleTU_BosonicGrid().positive_freqs_count){
	bubble = (-bubble_ph[W_xph][w_][m_][mp_][0][0][0][0](K_));
	bubble_dot = (-bubble_ph_dot[W_xph][w_][m_][mp_][0][0][0][0](K_));
      }
      
      val += vertex_dot_p_in_m * bubble * bare_vertex_mp_p_in * Gvec[W_][p_in_plus_P_](0,0) +
	vertex_p_in_m * bubble_dot * bare_vertex_mp_p_in * Gvec[W_][p_in_plus_P_](0,0) +
	vertex_p_in_m * bubble * bare_vertex_dot_mp_p_in * Gvec[W_][p_in_plus_P_](0,0) +
	vertex_p_in_m * bubble * bare_vertex_mp_p_in * Gdotvec[W_][p_in_plus_P_](0,0);// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }

    //for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
    //  val += Gvec[w] * ;
    
    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ) and a matsubara sum with 1/Beta

    val /= Model::MomentumGrid().ptr->get_volume() * Model::MomentumGrid().ptr->get_volume() * inv_freq_sum_normalisation(t) * inv_freq_sum_normalisation(t); //todo: thinking about replacing BETA with inverse_frequency_factor (note, this should remain beta because regardless of cutoff it should hold. 

    dcomplex val_Fock( 0.0, 0.0 );  
    // caution: enable flag SET_SIG_TAIL_TO_ZERO
    for( int p_ = 0; p_ < Model::GetFineMomentaCount(); ++p_ )
    for_freq( int W_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), W_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++W_ )
    {     
	val_Fock += Gvec[W_][p_](0,0) *// FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[W_] *   /**< Minus sign from internal loop */
	    state.ferm_vertex_sc_bare_dot( 00, 00, w, 0, W_, 0, t );
	val_Fock += Gdotvec[W_][p_](0,0) * //FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[W_] *   /**< Minus sign from internal loop */
	    state.ferm_vertex_sc_bare( 00, 00, w, 0, W_, 0, t );
    }      
  
    val_Fock *= 1.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount();

    
    return val + val_Fock;  
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_t<Model, state_t>::eval_Sig_SDE_Density_Conventional( const idx_1p_t<Model>& idx, const state_t& state, const state_t& dstate_over_dt, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Gdotvec, const gf_bubble_mat_t<Model> &bubble_ph, const gf_bubble_mat_t<Model> &bubble_ph_dot,  const double t )
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

      dcomplex vertex_p_in_m(0, 0), vertex_dot_p_in_m(0, 0);
      dcomplex bare_vertex_mp_p_in(0, 0), bare_vertex_dot_mp_p_in(0, 0);
      
      //const auto bare_Um = (Model::B_m(W_xph, K_, 0, 0, 0, 0) + BOSONIC_INTERACTION_M_SHIFT * Model::MomentumGrid().ptr->get_volume());
      //const auto bare_Um_dot = bare_Um * fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t);


      for (int n_ = 0; n_ < Model::GetMomentumFormFactorsCount(); ++n_ ){
	const dcomplex f_n_p_in = Model::GetFormFactorInMomentumIdxSpace(n_)[p_in];
	const dcomplex f_n_p_in_star = std::conj(Model::GetFormFactorInMomentumIdxSpace(n_)[p_in]);
	vertex_p_in_m += f_n_p_in_star * state.vertx_d( W_xph, K_, w_xph, n_, w_, m_, t );
	// warning: only works for density-density interactions in frequency space; should be reviewed for other cases

	const auto B_dot2 = state.bos_vertex_m_bare_dot( W_xph /*can be whatever*/, K_, W_xph, mp_, 0, n_, t );
	const auto F_dot2 = state.ferm_vertex_m_bare_dot( W_xph /*can be whatever*/, K_, W_xph, mp_, 0, n_, t );

	bare_vertex_mp_p_in += f_n_p_in * state.vertex_m_bare( W_xph, K_, W_xph, mp_, 0, n_, t );
	
	//B^M(q, k, k') + F(q, k, k') = U(k - k') ~  1/(\omega_0^2 + (k-k')^2)

	vertex_dot_p_in_m += f_n_p_in_star * (dstate_over_dt.phi_d_dot(state, W_xph, K_, w_xph, n_, w_, m_, t ) + dstate_over_dt.phi_d_bar_dot(state, W_xph, K_, w_xph, n_, w_, m_, t ) + state.bos_vertex_d_bare_dot( W_xph , K_, w_xph, n_, w_, m_, t ) + state.ferm_vertex_d_bare_dot( W_xph, K_, w_xph, n_, w_, m_, t ) );
	bare_vertex_dot_mp_p_in += f_n_p_in * ( (B_dot2 + F_dot2));
      }
      //std::cout << W_ << "," << w << "," << FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count << "," << w_plus_W_over_2 << "," << 2*(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count-w) <<  std::endl;
	
      dcomplex bubble = 0.0;
      dcomplex bubble_dot = 0.0;
      if (std::abs(W_)-1 < FrequencyDependenceScheme<Model>::BubbleTU_BosonicGrid().positive_freqs_count){
	bubble = (-bubble_ph[W_xph][w_][m_][mp_][0][0][0][0](K_));
	bubble_dot = (-bubble_ph_dot[W_xph][w_][m_][mp_][0][0][0][0](K_));
      }
      
      val += vertex_dot_p_in_m * bubble * bare_vertex_mp_p_in * Gvec[W_][p_in_plus_P_](0,0) +
	vertex_p_in_m * bubble_dot * bare_vertex_mp_p_in * Gvec[W_][p_in_plus_P_](0,0) +
	vertex_p_in_m * bubble * bare_vertex_dot_mp_p_in * Gvec[W_][p_in_plus_P_](0,0) +
	vertex_p_in_m * bubble * bare_vertex_mp_p_in * Gdotvec[W_][p_in_plus_P_](0,0);// * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }

    //for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
    //  val += Gvec[w] * ;
    
    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ) and a matsubara sum with 1/Beta
    val /= Model::MomentumGrid().ptr->get_volume() * Model::MomentumGrid().ptr->get_volume() * BETA * BETA; //todo: thinking about replacing BETA with inverse_frequency_factor 
    
    return -val;  
}
