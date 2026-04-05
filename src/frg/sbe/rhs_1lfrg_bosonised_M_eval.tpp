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
 * \brief Flow equations for the rest function in the sc, d and m channel. These diagramms are U-irreducible but bubble reducible
 *
 *
 * for loop over frequency might get fixed boundary
 * orbital dependencies omitted
 *
 */

// nu_0(W_idx)
// minus_nu_0(W_idx)


template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_wM_sc( const idx_w_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state, const gf_bubble_mat_t<Model>& bubble_pp, int sgn, double t)
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( IW_M::W );

   int K      = idx( IW_M::K );
   int n_in   = idx( IW_M::m );
   int n_out  = idx( IW_M::n );
   
   const double sign = double(sgn);

   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
       for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_ )
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
	   for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
	       val += state.B_irreducible_vertex_sc( W, K, nu_0(W), n_in, w_, m_, t ) *
                    bubble_pp[W][w_][m_][mp_][0][0][0][0](K) *
		   state.B_irreducible_vertex_sc( W, K, w_, mp_, nu_0(W), n_out, t );

	       val += sign*state.B_irreducible_vertex_sc( W, K, nu_0(W), n_in, w_, m_, t ) *
                    bubble_pp[W][w_][m_][mp_][0][0][0][0](K) *
		   state.B_irreducible_vertex_sc( W, K, w_, mp_, minus_nu_0(W), n_out, t );
         }
	   val *= 1.0/2.0;
   val*= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();
   return val;
   
}

template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_wM_d( const idx_w_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state, const gf_bubble_mat_t<Model>& bubble_ph, int sgn, double t )
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( IW_M::W );

   int K      = idx( IW_M::K );
   int n_in   = idx( IW_M::m );
   int n_out  = idx( IW_M::n );

   const double sign = double(sgn);


   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
   for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W,2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_ )

#ifdef SET_MIXED_BUBBLES_TO_ZERO
   {
       int mp_ = m_;
#else
   for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
       val += state.B_irreducible_vertex_d( W, K, nu_0(W), n_in, w_, m_, t ) *
	   bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
	   state.B_irreducible_vertex_d( W, K, w_, mp_, nu_0(W), n_out, t );
       val += sign*state.B_irreducible_vertex_d( W, K, nu_0(W), n_in, w_, m_, t ) *
	   bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
	   state.B_irreducible_vertex_d( W, K, w_, mp_, minus_nu_0(W), n_out, t );
   }
   val*= 1/2.0;
   val*= -1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();
   return val;
}


template <typename Model, typename state_t >
    dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_wM_m( const idx_w_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state, const gf_bubble_mat_t<Model>& bubble_ph, int sgn, double t )
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( IW_M::W );
   
   int K      = idx( IW_M::K );
   int n_in   = idx( IW_M::m );
   int n_out  = idx( IW_M::n );

   const double sign = double(sgn);

   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
   for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W,2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_ )
#ifdef SET_MIXED_BUBBLES_TO_ZERO
   {
       int mp_ = m_;
#else
   for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
       val += state.B_irreducible_vertex_m( W, K, nu_0(W), n_in, w_, m_, t ) *
	   bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
	   state.B_irreducible_vertex_m( W, K, w_, mp_, nu_0(W), n_out, t );
       val += sign*state.B_irreducible_vertex_m( W, K, nu_0(W), n_in, w_, m_, t ) *
	   bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
	   state.B_irreducible_vertex_m( W, K, w_, mp_, minus_nu_0(W), n_out, t );
   }
   val*= 1.0/2.0;
   val*= -1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();
   return val;
}

   
template <typename Model, typename state_t >
    dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambdaM_sc( const idx_lambda_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state, const gf_bubble_mat_t<Model>& bubble_pp, int sgn, const gf_w_M_t<Model> &wM_dot_gf, double t)
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( ILAMBDA_M::W );
   int w_in   = idx( ILAMBDA_M::w );
   
   int K      = idx( ILAMBDA_M::K );
   int n_in   = idx( ILAMBDA_M::m );
   int n_out  = idx( ILAMBDA_M::n );

   const double sign = double(sgn);

   dcomplex wM;
   if (sign > 0) 
       wM = state.wM_sc_plus(W, K, n_in, n_out);
   else 
       wM = state.wM_sc_minus(W, K, n_in, n_out);

   //if (std::abs(wM) < 1e-7)
   // return val;

   dcomplex wM_dot = wM_dot_gf[W][K][n_in][n_out];
   

   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
       for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_ )
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
	   for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif

	       val += state.B_irreducible_vertex_sc( W, K, w_in, n_in, w_, m_, t ) *
                    bubble_pp[W][w_][m_][mp_][0][0][0][0](K) *
		   state.B_irreducible_vertex_sc( W, K, w_, mp_, nu_0(W), n_out, t );

	       val += sign*state.B_irreducible_vertex_sc( W, K, w_in, n_in, w_, m_, t ) *
                    bubble_pp[W][w_][m_][mp_][0][0][0][0](K) *
		   state.B_irreducible_vertex_sc( W, K, w_, mp_, minus_nu_0(W), n_out, t );

	   }
	   //val *= 1.0/(2.0 * wM);
     val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();
	   
     //val += - wM_dot*(state.M_sc(W, K, w_in, n_in, nu_0(W), n_out) + sign * state.M_sc(W, K, w_in, n_in, minus_nu_0(W), n_out) )/(2.0 * wM*wM);

     return val;
}

template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambdaM_d( const idx_lambda_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state, const gf_bubble_mat_t<Model>& bubble_ph, int sgn, const gf_w_M_t<Model> &wM_dot_gf, double t)
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( ILAMBDA_M::W );
   int w_in   = idx( ILAMBDA_M::w );
   
   int K      = idx( ILAMBDA_M::K );
   int n_in   = idx( ILAMBDA_M::m );
   int n_out  = idx( ILAMBDA_M::n );

   const double sign = double(sgn);

   dcomplex wM;
   if (sign > 0) 
       wM = state.wM_d_plus(W, K, n_in, n_out);
   else 
       wM = state.wM_d_minus(W, K, n_in, n_out);

   //if (std::abs(wM) < 1e-7)
   //    return val;


   dcomplex wM_dot = wM_dot_gf[W][K][n_in][n_out];
   

   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
       for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_ )
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
	   for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif

	       val += state.B_irreducible_vertex_d( W, K, w_in, n_in, w_, m_, t ) *
                    bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
		   state.B_irreducible_vertex_d( W, K, w_, mp_, nu_0(W), n_out, t );

	       val += sign*state.B_irreducible_vertex_d( W, K, w_in, n_in, w_, m_, t ) *
                    bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
		   state.B_irreducible_vertex_d( W, K, w_, mp_, minus_nu_0(W), n_out, t );

	   }
	   //val *= 1.0/(2.0 * wM);
     val *= -1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();
	   
     //val += - wM_dot*(state.M_d(W, K, w_in, n_in, nu_0(W), n_out) + sign * state.M_d(W, K, w_in, n_in, minus_nu_0(W), n_out) )/(2.0 * wM*wM);

     return val;
}

template <typename Model, typename state_t >
    dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambdaM_m( const idx_lambda_M_t<Model>& idx, const state_frg_sbe_bosonised_M_t<Model>& state, const gf_bubble_mat_t<Model>& bubble_ph, int sgn, const gf_w_M_t<Model> &wM_dot_gf, double t)
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( ILAMBDA_M::W );
   int w_in   = idx( ILAMBDA_M::w );
   
   int K      = idx( ILAMBDA_M::K );
   int n_in   = idx( ILAMBDA_M::m );
   int n_out  = idx( ILAMBDA_M::n );

   const double sign = double(sgn);

   dcomplex wM;
   if (sign > 0) 
       wM = state.wM_m_plus(W, K, n_in, n_out);
   else 
       wM = state.wM_m_minus(W, K, n_in, n_out);
  
   //if (std::abs(wM) < 1e-7)
   //    return val;

   dcomplex wM_dot = wM_dot_gf[W][K][n_in][n_out];

   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ )
       for_freq( int w_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), w_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++w_ )
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
	   for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif

	       val += state.B_irreducible_vertex_m( W, K, w_in, n_in, w_, m_, t ) *
                    bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
		   state.B_irreducible_vertex_m( W, K, w_, mp_, nu_0(W), n_out, t );

	       val += sign*state.B_irreducible_vertex_m( W, K, w_in, n_in, w_, m_, t ) *
                    bubble_ph[W][w_][m_][mp_][0][0][0][0](K) *
		   state.B_irreducible_vertex_m( W, K, w_, mp_, minus_nu_0(W), n_out, t );

	   }
	   //val *= 1.0/(2.0 * wM);
     val *= -1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();
	   
     //val += - wM_dot*(state.M_m(W, K, w_in, n_in, nu_0(W), n_out) + sign * state.M_m(W, K, w_in, n_in, minus_nu_0(W), n_out) )/(2.0 * wM*wM);

     return val;
}
