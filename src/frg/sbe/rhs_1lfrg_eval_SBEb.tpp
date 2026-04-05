//#include <translate.h>
#include <mymath.h>
#include <frequencies/matsubara_space.h>
#include <cmath>
#include <complex>
#include <frg/flows.h>
#include <models/square_hubbard.h>




template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_SBEb_lambda_sc( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

#ifdef SET_MIXED_BUBBLES_TO_ZERO
    {
	int mp_ = 0; 
#else        
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
	for_freq(int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_){
	    val += bubble_pp[W][wp_][0][mp_][0][0][0][0](K) * state.B_irreducible_vertex_sc( W, K, wp_, mp_, w, m, t );  
	}
#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
	
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_pp( W, t ) * 
	    state.B_irreducible_vertex_sc( W, K, wp_, mp_, w, m, t );
#endif 
    } 
   
    val *= 1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 
    return val; 
}

template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_SBEb_lambda_d( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );


#ifdef SET_MIXED_BUBBLES_TO_ZERO
    {
	int mp_ = 0; 
#else        
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
	for_freq(int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_){
	    val += bubble_ph[W][wp_][0][mp_][0][0][0][0](K) * state.B_irreducible_vertex_d( W, K, wp_, mp_, w, m, t );
	}
#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
	
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_ph( W, t ) * 
	    state.B_irreducible_vertex_d( W, K, wp_, mp_, w, m, t );
#endif 
    } 

    val *= -1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 
    return val;
}

template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_SBEb_lambda_m( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );


#ifdef SET_MIXED_BUBBLES_TO_ZERO
    {
	int mp_ = 0; 
#else        
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
	for_freq(int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_){
	    val += bubble_ph[W][wp_][0][mp_][0][0][0][0](K) * state.B_irreducible_vertex_m( W, K, wp_, mp_, w, m, t );
	}
#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
	
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_ph( W, t ) * 
	    state.B_irreducible_vertex_m( W, K, wp_, mp_, w, m, t );
#endif 
    } 

    val *= -1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 
    return val; 
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_SBEb_w_sc( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );
    
    dcomplex w_sc = state.w_sc(W,K, t);

#ifdef SET_MIXED_BUBBLES_TO_ZERO
    {
	int mp_ = 0; 
#else       
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
	for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_ ){
	    val += bubble_pp[W][wp_][mp_][0][0][0][0][0](K) 
		* (2.0*state.lambda_sc( W, K, wp_, mp_ )- (1.0 ? mp_ == 0 : 0.0) );
	}
        
#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
	
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_pp( W, t ) * 
	    (2.0*state.lambda_sc( W, K, wp_, mp_ )- (1.0 ? mp_ == 0 : 0.0) );
#endif 
    } 
    
    val *= w_sc * w_sc / inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 
    return val; 
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_SBEb_w_d( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int K   = idx( ILAMBDA::K );

    dcomplex w_d = state.w_d( W, K, t );

#ifdef SET_MIXED_BUBBLES_TO_ZERO
    {
	int mp_ = 0; 
#else       
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
	for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_ ){
	    val += bubble_ph[W][wp_][mp_][0][0][0][0][0](K) 
		* (2.0*state.lambda_d( W, K, wp_, mp_ )- (1.0 ? mp_ == 0 : 0.0) );
	}
        
#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
	
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_ph( W, t ) * 
	    (2.0*state.lambda_d( W, K, wp_, mp_ )- (1.0 ? mp_ == 0 : 0.0) );
#endif 
    } 


    val *= -w_d * w_d / inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */
    return val;
}

template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_SBEb_w_m( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int K   = idx( ILAMBDA::K );

    dcomplex w_m = state.w_m( W, K, t );

#ifdef SET_MIXED_BUBBLES_TO_ZERO
    {
	int mp_ = 0; 
#else       
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
	for_freq( int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_ ){
	    val += bubble_ph[W][wp_][mp_][0][0][0][0][0](K) 
		* (2.0*state.lambda_m( W, K, wp_, mp_ )- (1.0 ? mp_ == 0 : 0.0) );
	}
        
#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */
	
	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_SGpGS_ph( W, t ) * 
	    (2.0*state.lambda_m( W, K, wp_, mp_ )- (1.0 ? mp_ == 0 : 0.0) );
#endif 
    } 

    val *= -1.0 * w_m * w_m / inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();     /*< 1.0/BETA/4.0/PI/PI -> Normalize freq summation and momentum integration V_BZ = (2 * PI)^2 */ 
    return val; 
}
