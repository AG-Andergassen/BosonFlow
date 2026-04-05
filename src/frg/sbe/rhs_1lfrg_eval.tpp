//#include <translate.h>
#include <mymath.h>
#include <frequencies/matsubara_space.h>
#include <cmath>
#include <complex>
#include <grid.h>
#include <frg/flows.h>
#include <models/square_hubbard.h>


using namespace std; 
using boost::bind;

/**
 * \brief Flow equations for the Hedin vertex lambda in the pp, ph and xph channel
 *
 *
 * for loop over frequency might get fixed boundary
 * orbital dependencies omitted
 *
 */
template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambda_sc( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp, double t )
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
	for_freq(int wp_, negative_frequency_integration_range, wp_ < positive_frequency_integration_range, ++wp_){
	    dcomplex temp( 0.0, 0.0 );
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	    {
		int mpp_ = mp_;
#else
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
		temp += state.lambda_sc( W, K, wp_, mpp_ ) *  
		    bubble_pp[W][wp_][mpp_][mp_][0][0][0][0](K) ;
	    }
	    val += temp * state.B_irreducible_vertex_sc( W, K, wp_, mp_, w, m, t );
	    }
#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
	
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_pp( W, t ) * 
	    state.lambda_sc( W, K, wp_, mp_ ) *  
	    state.B_irreducible_vertex_sc( W, K, wp_, mp_, w, m, t );
#endif 
    }
   
    val *= 1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 
    return val; 
}

   
template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambda_d( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
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
		    bubble_ph[W][wp_][mpp_][mp_][0][0][0][0](K);
	    }
	    val += temp * state.B_irreducible_vertex_d( W, K, wp_, mp_, w, m, t);
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   /**< explicit large frequency contribution -- alternative to weight_vec */

	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_ph( W, t ) *
	    state.lambda_d( W, K, wp_, mp_ ) *
	    state.B_irreducible_vertex_d( W, K, wp_, mp_, w, m, t );
#endif
    }

    val *= -1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 
    return val;
}

template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambda_m( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
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
		    bubble_ph[W][wp_][mpp_][mp_][0][0][0][0](K); 
	    }
	    val += temp * state.B_irreducible_vertex_m( W, K, wp_, mp_, w, m, t);
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   /**< explicit large frequency contribution -- alternative to weight_vec */
 
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_ph( W, t ) *
	    state.lambda_m( W, K, wp_, mp_ ) *  
	    state.B_irreducible_vertex_m( W, K, wp_, mp_, w, m, t );
#endif 
    }

    val *= -1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 
    return val; 
}


/**
 * \brief Flow equations for the screened interaction w in the pp, ph and xph channel
 *
 *
 * for loop over frequency might get fixed boundary
 * orbital dependencies omitted
 *
 */
template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_w_sc( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );
    
    dcomplex w_sc = state.w_sc(W,K,t);
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
	for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_ ){
	    dcomplex temp( 0.0, 0.0 );
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	    {
		int mpp_ = mp_;
#else
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
		temp += bubble_pp[W][wp_][mp_][mpp_][0][0][0][0](K) 
		    * state.lambda_sc( W, K, wp_, mpp_ );
	    }
	    val += state.lambda_sc( W, K, wp_, mp_ ) * temp;
	}
        
#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
	
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_pp( W, t ) * 
	    state.lambda_sc( W, K, wp_, mp_ ) *  
	    state.lambda_sc( W, K, wp_, mp_ );
#endif 
    } 
    
    val *= w_sc * w_sc / inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 

    return val; 
}

template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_w_d( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int K   = idx( ILAMBDA::K );

    dcomplex w_d = state.w_d( W, K, t );
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
	{
	    for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_ )
		{
		    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
		    {
			int mpp_ = mp_;
#else
		    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
			temp += bubble_ph[W][wp_][mp_][mpp_][0][0][0][0](K) * state.lambda_d( W, K, wp_, mpp_ );
		    }
		    val += temp * state.lambda_d( W, K, wp_, mp_ );
		}

#ifndef STATIC_CALCULATION         
	    int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   /**< explicit large frequency contribution -- alternative to weight_vec */

	    val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_ph( W, t ) *
		state.lambda_d( W, K, wp_, mp_ ) *
		state.lambda_d( W, K, wp_, mp_ );
#endif
	}

    val *= -w_d * w_d / inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */

    return val;
}

template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_w_m( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int K   = idx( ILAMBDA::K );

    dcomplex w_m = state.w_m( W, K, t );
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ )
	{
	    for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_ )
		{
		    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
		    {
			int mpp_ = mp_;
#else
		    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
			    temp += bubble_ph[W][wp_][mp_][mpp_][0][0][0][0](K) * state.lambda_m( W, K, wp_, mpp_ ); 
			}
		    val += temp * state.lambda_m( W, K, wp_, mp_ );
		}

#ifndef STATIC_CALCULATION         
	    int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   /**< explicit large frequency contribution -- alternative to weight_vec */

	    val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_ph( W, t ) *
		state.lambda_m( W, K, wp_, mp_ ) *  
		state.lambda_m( W, K, wp_, mp_ ); 
#endif
	}

    val *= -1.0 * w_m * w_m / inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 

    return val; 
}


/**
 * \brief Flow equations for the rest function in the sc, d and m channel. These diagramms are U-irreducible but bubble reducible
 *
 *
 * for loop over frequency might get fixed boundary
 * orbital dependencies omitted
 *
 */

template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_M_sc( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( IM::W );
   int w_in   = idx( IM::w );
   int w_out  = idx( IM::wp);

   int K      = idx( IM::K );
   int n_in   = idx( IM::m );
   int n_out  = idx( IM::mp );

#ifdef SET_OFFDIAGONAL_M_IN_FROMFACTOR_SPACE_TO_ZERO
 if (n_in != n_out)
       return val;
#endif

#ifdef BOSONISE_M_LAZY
   if (n_in != n_out)
       return val;

   // shift w_out to capture a "good" v0 and -v0
   const int sign_of_w_out = (int(w_out >= 0) - int(w_out < 0));
   w_out +=  sign_of_w_out * abs(div2_floor(W)) - int(w_out < 0) *( 1 + ((W+100000)%2));
#endif

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif

   
   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
       for_freq( int w_, negative_frequency_integration_range, w_ < positive_frequency_integration_range, ++w_ )
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
	   for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
	       val += state.B_irreducible_vertex_sc( W, K, w_in, n_in, w_, m_, t ) *
                    bubble_pp[W][w_][m_][mp_][0][0][0][0](K) *
		   state.B_irreducible_vertex_sc( W, K, w_, mp_, w_out, n_out, t );
         }
   val*= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();
   return val;
}

   
template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_M_d( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( IM::W );
   int w_in   = idx( IM::w );
   int w_out  = idx( IM::wp);

   int K      = idx( IM::K );
   int n_in   = idx( IM::m );
   int n_out  = idx( IM::mp );

#ifdef SET_OFFDIAGONAL_M_IN_FROMFACTOR_SPACE_TO_ZERO
 if (n_in != n_out)
       return val;
#endif


#ifdef BOSONISE_M_LAZY
   if (n_in != n_out)
       return val;

   // shift w_out to capture a "good" v0 and -v0
   const int sign_of_w_out = (int(w_out >= 0) - int(w_out < 0));
   w_out +=  sign_of_w_out * abs(div2_floor(W)) - int(w_out < 0) *( 1 + ((W+100000)%2));
#endif

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif

   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
   for_freq( int w_, negative_frequency_integration_range, w_ < positive_frequency_integration_range, ++w_ )

#ifdef SET_MIXED_BUBBLES_TO_ZERO
   {
       int mp_ = m_;
#else
   for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
       val += state.B_irreducible_vertex_d( W, K, w_in, n_in, w_, m_, t ) *
	   bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
	   state.B_irreducible_vertex_d( W, K, w_, mp_, w_out, n_out, t );
   }
   val*= -1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();
   return val;
}

   
template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_M_m( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( IM::W );
   int w_in   = idx( IM::w );
   int w_out  = idx( IM::wp);

   int K      = idx( IM::K );
   int n_in   = idx( IM::m );
   int n_out  = idx( IM::mp );

#ifdef SET_OFFDIAGONAL_M_IN_FROMFACTOR_SPACE_TO_ZERO
 if (n_in != n_out)
       return val;
#endif

#ifdef BOSONISE_M_LAZY
   if (n_in != n_out)
       return val;

   // shift w_out to capture a "good" v0 and -v0
   const int sign_of_w_out = (int(w_out >= 0) - int(w_out < 0));
   w_out +=  sign_of_w_out * abs(div2_floor(W)) - int(w_out < 0) *( 1 + ((W+100000)%2));
#endif

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif


   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
   for_freq( int w_, negative_frequency_integration_range, w_ < positive_frequency_integration_range, ++w_ )
#ifdef SET_MIXED_BUBBLES_TO_ZERO
   {
       int mp_ = m_;
#else
   for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
       val += state.B_irreducible_vertex_m( W, K, w_in, n_in, w_, m_, t ) *
	   bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
	   state.B_irreducible_vertex_m( W, K, w_, mp_, w_out, n_out, t );
   }
   val*= -1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();
   return val;
}
   
// note: does not mention orbital indices (compared to the parquet version)...
template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_Sig_conv( const idx_1p_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Svec, const double t )
{
    dcomplex val( 0.0, 0.0 );

    int w = idx( I1P::w ); 
    int p_in = Model::GetFineMomentumIdxFromCoarse(idx(I1P::k));  
    
    for( int p_ = 0; p_ < Model::GetFineMomentaCount(); ++p_ )
    for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange(), ++w_ )
    {     
	val += - Svec[w_][p_](0,0) * FrequencyDependenceScheme<Model>::IntegrationWeightVecs1D()[w_] *   /**< Minus sign from internal loop */
	    (2. * state.vertex_4pt( w , w_, w, p_in, p_, p_in, t ) -
	     state.vertex_4pt( w_, w, w, p_, p_in, p_in, t ) ); 
    }       
    
    return val *= 1.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount();     /*< 1.0/BETA/PATCH_COUNT -> Normalize freq/momentum summation */                 
}
