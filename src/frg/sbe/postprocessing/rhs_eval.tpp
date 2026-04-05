#include <mymath.h>
#include <frequencies/matsubara_space.h>
#include <cmath>
#include <complex>
#include <frg/flows.h>
#include <frg/flows.h>
#include <tu_projections.h>



// ----------------------------------------- Susceptiblities --------------------------------------------
// SUPERCONDUCTIVITY



template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t)
{
    //        __                ___ 
    //       /  \             /|   |\
    //  ~~~~.    .~~~~ + ~~~~. | V | .~~~~ 
    //       \__/             \|___|/

    dcomplex val( 0.0, 0.0 );

    // Introduce help variables
    int W = idx( ISUSC::W );

    int K = idx( ISUSC::K );

    int n_in = idx( ISUSC::n_in);
    int n_out = idx( ISUSC::n_out);

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	if (m_ != n_out) continue;
#endif
	for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	    if (mp_ != n_in) continue;
#endif
	    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
		for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++wp_ )
		{
			val += 1.0/inv_freq_sum_normalisation(t)/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume() *  // Two pairs of equivalent lines -> Factor 1/4
			    bubble_pp[W][w_][n_in][mp_][0][0][0][0](K) * state_vec.vertex_sc(W, K, w_, mp_, wp_, m_, t) *  
			    bubble_pp[W][wp_][m_][n_out][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs2D()[w_][wp_];
	        }
	}
    }

    if (fRGFlowScheme<Model>::G_ChosenFlowName() == G_FlowSchemeName::Interaction)
	val *= t*t;


    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	{
	    val += 1.0/inv_freq_sum_normalisation(t) *//Two equivalent lines -> factor 1/2
		bubble_pp[W][w_][n_in][n_out][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
	}

    return val;  
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_suscbubble_sc( const idx_susc_t<Model>& idx, const gf_bubble_mat_t<Model>& bubble_pp, double t)
{
    //         __ 
    //        /  \
    //   ~~~~.    .~~~~
    //        \__/

    dcomplex val( 0.0, 0.0 );

    // Introduce help variables
    int W = idx( ISUSC::W );

    int K = idx( ISUSC::K );

    int n_in = idx( ISUSC::n_in);
    int n_out = idx( ISUSC::n_out);

    if (W != 0)
	return val;

    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	{
	    val += 1.0/inv_freq_sum_normalisation(t) * bubble_pp[W][w_][n_in][n_out][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
	}
    return val; // ensure normalization matches paper definition
}


// blueprint to evaluate contributions to susceptibiliy by channel 
template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_verttypefunc( const state_t& state_vec, const idx_susc_t<Model>& idx, dcomplex (*verttypefunc)( const state_t&,int,int,int,int,int,int, const double), const gf_bubble_mat_t<Model>& bubble , double t)
{
    dcomplex val( 0.0, 0.0 );

    // Introduce help variables
    int W = idx( ISUSC::W );
    int K = idx( ISUSC::K );
    int n_in = idx( ISUSC::n_in);
    int n_out = idx( ISUSC::n_out);


    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	if (m_ != n_out) continue;
#endif
	for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	    if (mp_ != n_in) continue;
#endif
	    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
		for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++wp_ )
		    {
			val += 1.0/inv_freq_sum_normalisation(t)/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume() * bubble[W][w_][n_in][mp_][0][0][0][0](K) * (*verttypefunc)(state_vec, W, w_, wp_, K, mp_, m_, t ) * 
			    bubble[W][wp_][m_][n_out][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs2D()[w_][wp_];
		    }
	}
    }
    
    if (fRGFlowScheme<Model>::G_ChosenFlowName() == G_FlowSchemeName::Interaction)
	val *= t*t;


    return val;  
}


// blueprint to evaluate contributions to lambda by channel 
template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_verttypefunc( const state_t& state_vec, const idx_lambda_t<Model>& idx, dcomplex (*verttypefunc)( const state_t&,int,int,int,int,int,int, const double), const gf_bubble_mat_t<Model>& bubble, double bubble_sign, double t)
{
    // we emphasise that for D and M channels, one should pass "-bubble_ph", and "+bubble_pp" for the SC channel. 
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	if (mp_ != m) continue;
#endif
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	val += (*verttypefunc)(state_vec, W, w, w_, K, m, mp_, t )  * bubble_sign*bubble[W][w_][mp_][0][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
    }

	// per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ). The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta
	val /= Model::MomentumGrid().ptr->get_volume() * inv_freq_sum_normalisation(t);

    if (fRGFlowScheme<Model>::G_ChosenFlowName() == G_FlowSchemeName::Interaction)
	val *= t*std::sqrt(t); // the lambda is the vertex corresponding to a 3 point function -> interaction flow scales it with sqrt(t)^3
    
    return val;
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_verttypefunc( const state_t& state_vec, const idx_static_lambda_t<Model>& static_idx, dcomplex (*verttypefunc)( const state_t&,int,int,int,int,int,int, const double), const gf_bubble_mat_t<Model>& bubble, double bubble_sign, double t)
{
    idx_lambda_t<Model> idx;
    idx(ILAMBDA::K) = static_idx(ISTATICLAMBDA::K);
    idx(ILAMBDA::m) = static_idx(ISTATICLAMBDA::m);
    idx(ILAMBDA::W) = 0;
    idx(ILAMBDA::w) = 0;
    return eval_lambda_verttypefunc(state_vec, idx, verttypefunc, bubble, bubble_sign, t);
}


// blueprint to evaluate contributions to lambda by channel 
template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_verttypefunc( const state_t& state_vec, const idx_polarisation_t<Model>& idx, const gf_lambda_t<Model> &lambda, const gf_bubble_mat_t<Model>& bubble, double bubble_sign, double t)
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IPOLARISATION::W );
    int K   = idx( IPOLARISATION::K );

    const auto lambda_with_asymp = [&lambda, &state_vec](int W, int K, int w, int mp) -> dcomplex{ 
	if( (unsigned)(W +FrequenciesCount::lambda::POS_BOS) < FrequenciesCount::lambda::BOS && 
	    (unsigned)(w + FrequenciesCount::lambda::POS_FERM ) < (FrequenciesCount::lambda::FERM - (W+100000)%2))
	    return lambda[W][K][w][mp];
	const int w_last = w > 0 ? FrequenciesCount::lambda::POS_FERM : - FrequenciesCount::lambda::POS_FERM;
       
	return state_vec.lambda_sc_bare(W,K,w,mp) + 
	(lambda[W][K][w_last][mp] - state_vec.lambda_sc_bare(W,K,w_last,mp) ) * w_val(w_last) * w_val(w_last) / w_val(w) / w_val(w);  // todo: specialise
    };
    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
    for_freq(int w_, -FrequenciesCount::lambda::POS_FERM, w_ < FrequenciesCount::lambda::POS_FERM, ++w_)
	val += bubble_sign*bubble[W][w_][0][mp_][0][0][0][0](K) * lambda_with_asymp( W, K, w_, mp_ ) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    // per convention, a sum over form factors comes with 1/sqrt(Vol_BZ). The zero form factor need also another factor of 1/sqrt(Vol_BZ). The evaluation of the bubble at m = 0, also yield another factor of 1/sqrt(Vol_BZ). A matsubara sum comes with 1/Beta. 
    val /= Model::MomentumGrid().ptr->get_volume() * inv_freq_sum_normalisation(t); // todo: check if the volume factor is really necessary (or the previous comment is even correct!!)
    
    return val;
}



// DENSITY
template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t)
{
    dcomplex val( 0.0, 0.0 );

    // Introduce help variables
    int W = idx( ISUSC::W );
    int K = idx( ISUSC::K );
    int n_in = idx( ISUSC::n_in);
    int n_out = idx( ISUSC::n_out);


    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	if (m_ != n_out) continue;
#endif
	for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	    if (mp_ != n_in) continue;
#endif
	    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
		for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++wp_ )
		    {
			val += 1.0/inv_freq_sum_normalisation(t)/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume() * //Two internal loops
			    bubble_ph[W][w_][n_in][mp_][0][0][0][0](K) * 
			    state_vec.vertex_d( W, K, w_, mp_, wp_, m_ , t) *  
			    bubble_ph[W][wp_][m_][n_out][0][0][0][0](K) *
			    FrequencyDependenceScheme<Model>::IntegrationWeightVecs2D()[w_][wp_];
		    }
	}
    }

    if (fRGFlowScheme<Model>::G_ChosenFlowName() == G_FlowSchemeName::Interaction)
	val *= t*t;


    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	{
	    val -= 1/inv_freq_sum_normalisation(t) * // One internal loop
		bubble_ph[W][w_][n_in][n_out][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
	}

    return val; // ensure normalization matches paper definition
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_suscvert_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t)
{
    dcomplex val( 0.0, 0.0 );

    // Introduce help variables
    int W = idx( ISUSC::W );
    int K = idx( ISUSC::K );
    int n_in = idx( ISUSC::n_in);
    int n_out = idx( ISUSC::n_out);


    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	if (m_ != n_out) continue;
#endif
	for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	    if (mp_ != n_in) continue;
#endif
	    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
		for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++wp_ )
		    {
			val += 1.0/inv_freq_sum_normalisation(t)/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume()        * //Two internal loops
			    bubble_ph[W][w_][n_in][mp_][0][0][0][0](K) *
			    state_vec.vertex_d(W, K, w_, mp_, wp_, m_, t) *
			    bubble_ph[W][wp_][m_][n_out][0][0][0][0](K) *
			    FrequencyDependenceScheme<Model>::IntegrationWeightVecs2D()[w_][wp_];
		    }
	}
    }

    if (fRGFlowScheme<Model>::G_ChosenFlowName() == G_FlowSchemeName::Interaction)
	val *= t*t;

    return val; // ensure normalization matches paper definition
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_suscbubble_d( const idx_susc_t<Model>& idx, const gf_bubble_mat_t<Model>& bubble_ph, double t)
{
    dcomplex val( 0.0, 0.0 );

    // Introduce help variables
    int W = idx( ISUSC::W );
    int K = idx( ISUSC::K );
    int n_in = idx( ISUSC::n_in);
    int n_out = idx( ISUSC::n_out);

    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	{
	    val -= 1/inv_freq_sum_normalisation(t) * // One internal loop
		bubble_ph[W][w_][n_in][n_out][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
	}

    return val;  // ensure normalization matches paper definition
}


template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_bare( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_bare, bubble_ph, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_nabla_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_nabla_d, bubble_ph, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_M_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_M_d, bubble_ph, t );
}


template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_d_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_d_double_counting, bubble_ph, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_Virr( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_Virr, bubble_ph, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_nabla_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_nabla_m, bubble_ph, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_M_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_M_m, bubble_ph, t );
}


template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_m_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_m_double_counting, bubble_ph, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_nabla_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_nabla_sc, bubble_ph, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_M_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_M_sc, bubble_ph, t );
}


template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_d_contribution_from_sc_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_sc_double_counting, bubble_ph, t );
}


// sc stuff


template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_bare( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_bare, bubble_pp, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_Virr( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_Virr, bubble_pp, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_nabla_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_nabla_d, bubble_pp, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_M_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_M_d, bubble_pp, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_d_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_d_double_counting, bubble_pp, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_nabla_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_nabla_m, bubble_pp, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_M_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_M_m, bubble_pp, t );
}


template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_m_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_m_double_counting, bubble_pp, t );
}


template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_nabla_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_nabla_sc, bubble_pp, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_M_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_M_sc, bubble_pp, t );
}


template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_sc_contribution_from_sc_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_sc_double_counting, bubble_pp, t );
}

template<typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_suscvert_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t)
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertex_sc, bubble_pp, t );
}



// MAGNETIC
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t)
{
    dcomplex val( 0.0, 0.0 );

    // Introduce help variables
    int W = idx( ISUSC::W );
    int K = idx( ISUSC::K );
    int n_in = idx( ISUSC::n_in);
    int n_out = idx( ISUSC::n_out);

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	if (m_ != n_out) continue;
#endif
	for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#ifdef SETMIXED_BUBBLES_TO_ZERO
	    if (mp_ != n_in) continue;
#endif
	    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
		for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++wp_ )
		    {
			val += 1.0/inv_freq_sum_normalisation(t)/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume() * bubble_ph[W][w_][n_in][mp_][0][0][0][0](K) * 
			    state_vec.vertex_m(W, K, w_, mp_, wp_, m_, t) *  
			    bubble_ph[W][wp_][m_][n_out][0][0][0][0](K) *
			    FrequencyDependenceScheme<Model>::IntegrationWeightVecs2D()[w_][wp_];
		    }
    
	}
    }

    if (fRGFlowScheme<Model>::G_ChosenFlowName() == G_FlowSchemeName::Interaction)
	val *= t*t;

    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	{
	    val -= 1/inv_freq_sum_normalisation(t) * // One internal loop
		bubble_ph[W][w_][n_in][n_out][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
	}

    return val;  // ensure normalization matches paper definition 
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_suscbubble_m( const idx_susc_t<Model>& idx, const gf_bubble_mat_t<Model>& bubble_ph, double t)
{
    dcomplex val( 0.0, 0.0 );

    // Introduce help variables
    int W = idx( ISUSC::W );
    int K = idx( ISUSC::K );
    int n_in = idx( ISUSC::n_in);
    int n_out = idx( ISUSC::n_out);

    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	{
	    val -= 1/inv_freq_sum_normalisation(t) * // One internal loop
		bubble_ph[W][w_][n_in][n_out][0][0][0][0](K) *
		FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];
	}

    return val;  // ensure normalization matches paper definition
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_suscvert_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t)
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertex_m, bubble_ph, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_bare( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_bare, bubble_ph, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_Virr( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_Virr, bubble_ph, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_nabla_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_nabla_d, bubble_ph, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_M_d( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_M_d, bubble_ph, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_d_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_d_double_counting, bubble_ph, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_nabla_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_nabla_m, bubble_ph, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_M_m( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_M_m, bubble_ph, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_m_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_m_double_counting, bubble_ph, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_nabla_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_nabla_sc, bubble_ph, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_M_sc( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_M_sc, bubble_ph, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_susc_m_contribution_from_sc_double_counting( const idx_susc_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    return eval_susc_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_sc_double_counting, bubble_ph, t );
}


//hedin vertices
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_contribution_from_1( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    if (idx( ISTATICLAMBDA::m ) == 0)
	return dcomplex(1.0, 0.0);
    return dcomplex(0.0, 0.0);
}


// d
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_bare( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_bare, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_nabla_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_d_contribution_from_nabla_d, bubble_ph, -1.0, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_nabla_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_nabla_sc, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_nabla_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_nabla_m, bubble_ph,-1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_M_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_M_d, bubble_ph, -1.0, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_M_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_M_sc, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_M_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_M_m, bubble_ph,-1.0, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_d_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_d_double_counting, bubble_ph,-1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_sc_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_sc_double_counting, bubble_ph,-1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_m_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_m_double_counting, bubble_ph,-1.0, t );
}


// m
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_bare( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_bare, bubble_ph, -1.0, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_nabla_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_m_contribution_from_nabla_m, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_nabla_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_nabla_sc, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_nabla_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_nabla_d, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_M_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_M_m, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_M_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_M_sc, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_M_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_M_d, bubble_ph, -1.0, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_d_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_d_double_counting, bubble_ph,-1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_sc_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_sc_double_counting, bubble_ph,-1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_m_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_m_double_counting, bubble_ph,-1.0, t );
}



// sc
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_bare( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_bare, bubble_pp, 1.0, t );
}

// it's zero, but we add it anyway
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_nabla_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_sc_contribution_from_nabla_sc, bubble_pp, 1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_nabla_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_nabla_m, bubble_pp, 1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_nabla_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_nabla_d, bubble_pp, 1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_M_sc( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_M_sc, bubble_pp, 1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_M_m( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_M_m, bubble_pp, 1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_M_d( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_M_d, bubble_pp, 1.0, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_sc_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_sc_double_counting, bubble_pp, 1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_d_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_d_double_counting, bubble_pp, 1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_m_double_counting( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_m_double_counting, bubble_pp, 1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_sc_contribution_from_varphi( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_sc, bubble_pp, 1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_d_contribution_from_varphi( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_d, bubble_ph, -1.0, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_lambda_m_contribution_from_varphi( const idx_static_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    return eval_lambda_verttypefunc( state_vec, idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_m, bubble_ph, -1.0, t );
}


// polarisation
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_sc_contribution_lambda( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t, const gf_lambda_t<Model> &lambda)
{
    return eval_polarisation_verttypefunc( state_vec, idx, lambda, bubble_pp, 1.0, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_d_contribution_lambda( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t, const gf_lambda_t<Model> &lambda)
{
    return eval_polarisation_verttypefunc( state_vec, idx, lambda, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_m_contribution_lambda( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t, const gf_lambda_t<Model> &lambda)
{
    return eval_polarisation_verttypefunc( state_vec, idx, lambda, bubble_ph, -1.0, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_contribution_suscvert( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t, const gf_susc_t<Model> &susc_contribution)
{
    dcomplex val(0.0, 0.0);
    
    // Introduce help variables
    int W = idx( IPOLARISATION::W );
    int K = idx( IPOLARISATION::K );

    val += susc_contribution[W][K][0][0][0][0][0][0];
    
    val /= Model::MomentumGrid().ptr->get_volume();

    return val;
}



template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_m_or_d_contribution_from_bubble( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    dcomplex val(0.0, 0.0);
    
    // Introduce help variables
    int W = idx( IPOLARISATION::W );
    int K = idx( IPOLARISATION::K );

    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	val += -bubble_ph[W][w_][0][0][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    val /= Model::MomentumGrid().ptr->get_volume()*inv_freq_sum_normalisation(t);

    return val;
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_sc_contribution_from_bubble( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    dcomplex val(0.0, 0.0);
    
    // Introduce help variables
    int W = idx( IPOLARISATION::W );
    int K = idx( IPOLARISATION::K );

    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	val += bubble_pp[W][w_][0][0][0][0][0][0](K) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_];

    val /= Model::MomentumGrid().ptr->get_volume()*inv_freq_sum_normalisation(t);

    return val;
}

// it's zero, but we add it anyway
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_sc_contribution_from_nabla_sc( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    idx_susc_t<Model> susc_idx;
    susc_idx( ISUSC::W ) = idx( IPOLARISATION::W );
    susc_idx( ISUSC::K ) = idx( IPOLARISATION::K );
    susc_idx( ISUSC::n_in ) = 0;
    susc_idx( ISUSC::n_out ) = 0;    

    return eval_susc_verttypefunc( state_vec, susc_idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_sc_contribution_from_nabla_sc, bubble_pp, t )/ Model::MomentumGrid().ptr->get_volume();
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_d_contribution_from_nabla_d( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    idx_susc_t<Model> susc_idx;
    susc_idx( ISUSC::W ) = idx( IPOLARISATION::W );
    susc_idx( ISUSC::K ) = idx( IPOLARISATION::K );
    susc_idx( ISUSC::n_in ) = 0;
    susc_idx( ISUSC::n_out ) = 0;    

    return eval_susc_verttypefunc( state_vec, susc_idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_d_contribution_from_nabla_d, bubble_ph, t )/ Model::MomentumGrid().ptr->get_volume();
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_m_contribution_from_nabla_m( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    idx_susc_t<Model> susc_idx;
    susc_idx( ISUSC::W ) = idx( IPOLARISATION::W );
    susc_idx( ISUSC::K ) = idx( IPOLARISATION::K );
    susc_idx( ISUSC::n_in ) = 0;
    susc_idx( ISUSC::n_out ) = 0;    

    return eval_susc_verttypefunc( state_vec, susc_idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_m_contribution_from_nabla_m, bubble_ph, t )/ Model::MomentumGrid().ptr->get_volume();
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_sc_contribution_from_varphi( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp , double t)
{
    idx_susc_t<Model> susc_idx;
    susc_idx( ISUSC::W ) = idx( IPOLARISATION::W );
    susc_idx( ISUSC::K ) = idx( IPOLARISATION::K );
    susc_idx( ISUSC::n_in ) = 0;
    susc_idx( ISUSC::n_out ) = 0;    

    return eval_susc_verttypefunc( state_vec, susc_idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_sc, bubble_pp, t )/ Model::MomentumGrid().ptr->get_volume();
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_d_contribution_from_varphi( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    idx_susc_t<Model> susc_idx;
    susc_idx( ISUSC::W ) = idx( IPOLARISATION::W );
    susc_idx( ISUSC::K ) = idx( IPOLARISATION::K );
    susc_idx( ISUSC::n_in ) = 0;
    susc_idx( ISUSC::n_out ) = 0;    

    return eval_susc_verttypefunc( state_vec, susc_idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_d, bubble_ph, t )/ Model::MomentumGrid().ptr->get_volume();
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_polarisation_m_contribution_from_varphi( const idx_polarisation_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph , double t)
{
    idx_susc_t<Model> susc_idx;
    susc_idx( ISUSC::W ) = idx( IPOLARISATION::W );
    susc_idx( ISUSC::K ) = idx( IPOLARISATION::K );
    susc_idx( ISUSC::n_in ) = 0;
    susc_idx( ISUSC::n_out ) = 0;    

    return eval_susc_verttypefunc( state_vec, susc_idx, rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_m, bubble_ph, t )/ Model::MomentumGrid().ptr->get_volume();
}



// self-energy

//------ B G0 phi
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_B_G0_phi_sc( const idx_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( ILAMBDA::W );
    int w_in   = idx( ILAMBDA::w );
    int K      = idx( ILAMBDA::K );
    int n_in   = idx( ILAMBDA::m );

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
        for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_ ){
            for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
		{
		    val += (state_vec.nabla_sc( W, K, w_in, n_in, w_, m_, t ) - TUProjections<Model>::VertexBare()[n_in][m_][0][0][0][0]
#if !defined(SBEa_APPROXIMATION)
			    + state_vec.M_sc( W, K, w_in, n_in, w_, m_ )
#endif
)*
			bubble_pp[W][w_][m_][mp_][0][0][0][0](K) *
			TUProjections<Model>::VertexBare()[mp_][0][0][0][0][0];
		}
	}
        int w = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
        val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_pp( W, t ) * 
	    (state_vec.nabla_sc( W, K, w_in, n_in, w, m_, t ) - TUProjections<Model>::VertexBare()[n_in][m_][0][0][0][0]
#if !defined(SBEa_APPROXIMATION)
 + state_vec.M_sc( W, K, w_in, n_in, w, m_ )
#endif
) *  
	    TUProjections<Model>::VertexBare()[m_][0][0][0][0][0];
	// todo: add orbital dependence
    }
    val *= -1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< -1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 and multiply by minus to compensate minus in proj_vert_bare*/ 
    return val; 
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_B_G0_phi_d( const idx_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( ILAMBDA::W );
    int w_in   = idx( ILAMBDA::w );
    int K      = idx( ILAMBDA::K );
    int n_in   = idx( ILAMBDA::m );
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
        for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_ ){
            for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
		{
		    val +=  bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
                        (state_vec.nabla_d( W, K, w_in, n_in, w_, m_, t ) - TUProjections<Model>::VertexBare()[n_in][m_][0][0][0][0]
#if !defined(SBEa_APPROXIMATION)
 + state_vec.M_d( W, K, w_in, n_in, w_, m_ )
#endif
) *
                        TUProjections<Model>::VertexBare()[mp_][0][0][0][0][0];
		}
	}
        int w = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
        val += inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) * (state_vec.nabla_d( W, K, w_in, n_in, w, m_, t ) - TUProjections<Model>::VertexBare()[n_in][m_][0][0][0][0]
#if !defined(SBEa_APPROXIMATION)
 + state_vec.M_d( W, K, w_in, n_in, w, m_ )
#endif
) *
	    TUProjections<Model>::VertexBare()[m_][0][0][0][0][0];
    }
    
    val *= -1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< -1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 and multiply by minus to compensate minus in proj_vert_bare*/ 
    return val; 
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_B_G0_phi_m( const idx_lambda_t<Model>& idx, const state_t& state_vec, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( ILAMBDA::W );
    int w_in   = idx( ILAMBDA::w );
    int K      = idx( ILAMBDA::K );
    int n_in   = idx( ILAMBDA::m );
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
        for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_ ){
            for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
		{
		    val += (state_vec.nabla_m( W, K, w_in, n_in, w_, m_, t ) + TUProjections<Model>::VertexBare()[n_in][m_][0][0][0][0]
#if !defined(SBEa_APPROXIMATION)
			    + state_vec.M_m( W, K, w_in, n_in, w_, m_ )
#endif
) *  
			bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
			TUProjections<Model>::VertexBare()[mp_][0][0][0][0][0];
		}
	}
	int w = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) * (state_vec.nabla_m( W, K, w_in, n_in, w, m_, t ) + TUProjections<Model>::VertexBare()[n_in][m_][0][0][0][0]
#if !defined(SBEa_APPROXIMATION)
 + state_vec.M_m( W, K, w_in, n_in, w, m_ )
#endif
) *  
	    TUProjections<Model>::VertexBare()[m_][0][0][0][0][0];
    }
    val *= -1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< -1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 and multiply by minus to compensate minus in proj_vert_bare*/ 
    return val; 
}

//------ L(B+B+B,G)
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_rhs_Sig_LBG( const idx_1p_t<Model>& idx, const gf_lambda_t<Model>& B_sc, const gf_lambda_t<Model>& B_d, const gf_lambda_t<Model>& B_m, const gf_1p_mat_t<Model>& Gvec, double t )
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );   /*!< outside frequency. To be distiguished to internal integration frequency */

    
    int p_in = Model::CoarseToFineIdxMap()[idx(I1P::k)];   /*!< outside momentum. To be distiguished to internal integration momentum */
    
    for( int p = 0; p < Model::GetFineMomentaCount(); ++p )
        for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ ){
 
	    int p_plus_p_in = Model::FineMomentumGrid().momentum_plus_momentum_idx_map[p][p_in];

	    int p_minus_p_in = Model::FineMomentumGrid().momentum_plus_momentum_idx_map
		[p]
		[Model::FineMomentumGrid().negative_momentum_idx_map[p_in]];
 
 
            // translation for pp-channel
            int W_pp  = w_ + w + 1; 
            int w_pp  = w - div2_ceil(W_pp);
            int K_pp  = Model::FineToCoarseIdxMap()[p_plus_p_in];
            // translation for ph-channel
            int W_ph = w_ - w;
            int w_ph = w + div2_floor(W_ph);
            int K_ph  = Model::FineToCoarseIdxMap()[p_minus_p_in];
            // translation for xph-channel: same as ph
            int W_xph = W_ph;
            int w_xph = w_ph;
            int K_xph  = K_ph;  //TODO: simplify K_xph=K_ph

            for( int m_=0; m_ < Model::GetMomentumFormFactorsCount(); ++m_) 
		{          
		    auto &f_m_ = Model::GetFormFactorInFineMomentumIdxSpace(m_);
		    auto &f_0 = Model::GetFormFactorInFineMomentumIdxSpace(0);
	   
		    val +=  Gvec[w_][p](0,0) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_]      /**< Minus sign from internal loop */
			* std::conj(f_m_[p_in]) 
                        * f_0[p_in]
			//                * ( B_sc[W_pp][w_pp][K_pp][m_] + 0.5*B_d[W_ph][w_ph][K_ph][m_] - 0.5*B_m[W_ph][w_ph][K_ph][m_] - B_m[W_xph][w_xph][K_xph][m_] );
                        * ( B_sc[W_pp][K_pp][w_pp][m_] + 0.5*B_d[W_ph][K_ph][w_ph][m_] - 0.5*B_m[W_ph][K_ph][w_ph][m_] - B_m[W_xph][K_xph][w_xph][m_] );
		}
        }
    return val *= 1.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount();     /*< 1.0/BETA/PATCH_COUNT -> Normalize freq/momentum summation */ 
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_rhs_Sig_LBG_sc( const idx_1p_t<Model>& idx, const gf_lambda_t<Model>& B_sc, const gf_1p_mat_t<Model>& Gvec, double t )
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );   /*!< outside frequency. To be distiguished to internal integration frequency */
    int p_in = Model::CoarseToFineIdxMap()[idx(I1P::k)];   /*!< outside momentum. To be distiguished to internal integration momentum */
       /*!< outside momentum. To be distiguished to internal integration momentum */
    
    for( int p = 0; p < Model::GetFineMomentaCount(); ++p )
        for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
	    {
		int p_plus_p_in = Model::FineMomentumGrid().momentum_plus_momentum_idx_map[p][p_in];

		// translation for pp-channel
		int W_pp  = w_ + w + 1;
		int w_pp  = w - div2_ceil(W_pp);
	        int K_pp  = Model::FineToCoarseIdxMap()[p_plus_p_in];
		for( int m_=0; m_ < Model::GetMomentumFormFactorsCount(); ++m_)
		    {
			auto &f_m_ = Model::GetFormFactorInFineMomentumIdxSpace(m_);
			auto &f_0 = Model::GetFormFactorInFineMomentumIdxSpace(0);

			val +=  Gvec[w_][p](0,0) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_]      /**< Minus sign from internal loop */
			    * std::conj(f_m_[p_in]) 
			    * f_0[p_in]
			    * B_sc[W_pp][K_pp][w_pp][m_]; 
		    }
	    }
    return val *= 1.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount();     /*< 1.0/BETA/PATCH_COUNT -> Normalize freq/momentum summation */ 
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_rhs_Sig_LBG_d( const idx_1p_t<Model>& idx, const gf_lambda_t<Model>& B_d, const gf_1p_mat_t<Model>& Gvec, double t )
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );   /*!< outside frequency. To be distiguished to internal integration frequency */
    int p_in = Model::CoarseToFineIdxMap()[idx(I1P::k)];   /*!< outside momentum. To be distiguished to internal integration momentum */

    for( int p = 0; p < Model::GetFineMomentaCount(); ++p )
        for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ ){
	    int p_minus_p_in = Model::FineMomentumGrid().momentum_plus_momentum_idx_map
		[p]
		[Model::FineMomentumGrid().negative_momentum_idx_map[p_in]];
	    

            /* translation for ph-channel*/
            int W_ph = w_ - w;
            int w_ph = w + div2_floor(W_ph);
            int K_ph  = Model::FineToCoarseIdxMap()[p_minus_p_in];
            for( int m_=0; m_ < Model::GetMomentumFormFactorsCount(); ++m_)
		{
		    auto &f_m_ = Model::GetFormFactorInFineMomentumIdxSpace(m_);
		    auto &f_0 = Model::GetFormFactorInFineMomentumIdxSpace(0);
		   
		    val +=  Gvec[w_][p](0,0) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_]      /**< Minus sign from internal loop */
                        * std::conj(f_m_[p_in]) 
                        * f_0[p_in]
                        * B_d[W_ph][K_ph][w_ph][m_];
		}
	}
    return val *= 1.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount();     /*< 1.0/BETA/PATCH_COUNT -> Normalize freq/momentum summation */
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_rhs_Sig_LBG_m( const idx_1p_t<Model>& idx, const gf_lambda_t<Model>& B_m, const gf_1p_mat_t<Model>& Gvec, double t )
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w );   /*!< outside frequency. To be distiguished to internal integration frequency */
    int p_in = Model::CoarseToFineIdxMap()[idx(I1P::k)];   /*!< outside momentum. To be distiguished to internal integration momentum */

    for( int p = 0; p < Model::GetFineMomentaCount(); ++p )
        for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ ){
	    
	    int p_minus_p_in = Model::FineMomentumGrid().momentum_plus_momentum_idx_map
		[p]
		[Model::FineMomentumGrid().negative_momentum_idx_map[p_in]];
 
            /* translation for xph-channel*/
            int W_xph = w_ - w;
            int w_xph = w + div2_floor(W_xph);
            int K_xph  = Model::FineToCoarseIdxMap()[p_minus_p_in]; 

            for( int m_=0; m_ < Model::GetMomentumFormFactorsCount(); ++m_)
		{
		    auto &f_m_ = Model::GetFormFactorInFineMomentumIdxSpace(m_);
		    auto &f_0 = Model::GetFormFactorInFineMomentumIdxSpace(0);

		    val +=  Gvec[w_][p](0,0) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_]     /**< Minus sign from internal loop */
        		* std::conj(f_m_[p_in]) 
                        * f_0[p_in]
                        * B_m[W_xph][K_xph][w_xph][m_] ;
		}
	}
    return val *= 1.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount();     /*< 1.0/BETA/PATCH_COUNT -> Normalize freq/momentum summation */
}



template <typename Model, typename state_t > 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_Sig_SDE_integrand_m( const idx_Sig_SDE_integrand_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec, const double t )
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( ISIG_SDE_INTEGRAND::w );
    int k = idx(ISIG_SDE_INTEGRAND::k);
    int W = idx(ISIG_SDE_INTEGRAND::W);
    int K = idx(ISIG_SDE_INTEGRAND::K);


    const unsigned p = Model::GetFineMomentumIdxFromCoarse(k);  
    const unsigned P = Model::GetFineMomentumIdxFromCoarse(K);
    const unsigned p_plus_P = Model::SumFineMomentaIdxes(p, P);
    
    int w_plus_W_over_2 = w + div2_floor( W );
    
    val += Gvec[w + W][p_plus_P](0,0) * state.mom_lambda_m(W, K, w_plus_W_over_2, k) * state.w_m(W, K, t);
    
    return val *= 1.0/Model::GetFineMomentaCount()/std::sqrt(Model::MomentumGrid().ptr->get_volume()); 
}


template <typename Model, typename state_t > 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_Sig_SDE_integrand_d( const idx_Sig_SDE_integrand_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec, const double t )
{
  dcomplex val( 0.0, 0.0 );
  
  int w = idx( ISIG_SDE_INTEGRAND::w );
  int k = idx(ISIG_SDE_INTEGRAND::k);
  int W = idx(ISIG_SDE_INTEGRAND::W);
  int K = idx(ISIG_SDE_INTEGRAND::K);
  
  int p = Model::GetFineMomentumIdxFromCoarse(k);  
  int P = Model::GetFineMomentumIdxFromCoarse(K);  
  
  unsigned p_plus_P = Model::SumFineMomentaIdxes(p, P);
  int w_plus_W_over_2 = w + div2_floor( W );
  
  
  val += Gvec[w + W][p_plus_P](0,0) * state.mom_lambda_d(W, K, w_plus_W_over_2, k) * state.w_d(W, K, t);
  
  val -= 2. * std::sqrt(Model::MomentumGrid().ptr->get_volume()) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*SquareHubbardParams::UINT * Gvec[w + W][p_plus_P](0,0);
  
  
  return val *= 1.0/Model::GetFineMomentaCount()/std::sqrt(Model::MomentumGrid().ptr->get_volume());
  
}


template <typename Model, typename state_t > 
dcomplex rhs_postproc_sbe_t<Model, state_t>::eval_Sig_SDE_integrand_sc( const idx_Sig_SDE_integrand_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec, const double t )
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( ISIG_SDE_INTEGRAND::w );
    int k = idx(ISIG_SDE_INTEGRAND::k);
    int W = idx(ISIG_SDE_INTEGRAND::W);
    int K = idx(ISIG_SDE_INTEGRAND::K);

    int p = Model::GetFineMomentumIdxFromCoarse(k);  
    int P = Model::GetFineMomentumIdxFromCoarse(K);  

    unsigned P_minus_p = Model::SumFineMomentaIdxes(P, Model::GetNegativeFineMomentumIdx(p));
    unsigned K_minus_k = Model::GetCoarseMomentumIdxFromFine(P_minus_p);
    int W_over_2_minus_w = div2_floor( W ) - w - 1;
    val += Gvec[W - w - 1][P_minus_p](0,0) * state.mom_lambda_sc(W, K, W_over_2_minus_w, K_minus_k) * state.w_sc(W, K, t);

    
    return val *= -1.0/Model::GetFineMomentaCount()/std::sqrt(Model::MomentumGrid().ptr->get_volume());

 
}
