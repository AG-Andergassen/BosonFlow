#include <frg/flows.h>
#include <frg/symmetries_common.h>
#include <frg/rhs_base.h>
#include <frg/observables_common.h>

template <typename Model, typename state_t>
observer_frg_t<Model, state_t>::observer_frg_t(observables_frg_common_t<Model, state_t> &observables) : m_observables(observables)
{
}

template <typename Model, typename state_t >
bool observer_frg_t<Model, state_t>::operator() ( const state_t& state_vec, const double t, bool is_final )
{
    m_observables.update(state_vec, t, is_final);

    // Calculate maximal coupling
    double max_cpl = get_max_cpl( state_vec ); 
   max_cpl = max_cpl/4.0/M_PI/M_PI; 
    
    return ( max_cpl < MAX_COUPLING ) ? 1 : 0; 
}



/**
 * \brief Returns maximal value of the screened interaction among all channels, prints value and index for each channel
 *
 * used in the condition whether the flow has to be stopped due to divergent couplings
 *
 */
template <typename Model, typename state_t> 
double  observer_frg_t<Model, state_t>::get_max_cpl( const state_t& state )
{
    const dcomplex* max_sig_ptr = max_element( state.gf_Sig().data(), state.gf_Sig().data() + state.gf_Sig().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_sig_distance = std::distance( state.gf_Sig().data(), max_sig_ptr ); 
   std::cout  << " Max Sig " << *max_sig_ptr << " at idx " << state.gf_Sig().get_idx( max_sig_distance ) << std::endl; 


   const dcomplex* max_sc_ptr = max_element( state.gf_w_sc().data(), state.gf_w_sc().data() + state.gf_w_sc().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_sc_distance = std::distance( state.gf_w_sc().data(), max_sc_ptr ); 
   std::cout  << " Max w sc " << *max_sc_ptr << " at idx " << state.gf_w_sc().get_idx( max_sc_distance ) << std::endl; 

   const dcomplex* max_d_ptr = max_element( state.gf_w_d().data(), state.gf_w_d().data() + state.gf_w_d().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_d_distance = std::distance( state.gf_w_d().data(), max_d_ptr ); 
   std::cout  << " Max w d " << *max_d_ptr << " at idx " << state.gf_w_d().get_idx( max_d_distance ) << std::endl; 

   const dcomplex* max_m_ptr = max_element( state.gf_w_m().data(), state.gf_w_m().data() + state.gf_w_m().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_m_distance = std::distance( state.gf_w_m().data(), max_m_ptr ); 
   std::cout << " Max w m " << *max_m_ptr << " at idx " << state.gf_w_m().get_idx( max_m_distance ) << std::endl; 


#ifndef BOSONISE_M 
   const dcomplex* max_M_sc_ptr = max_element( state.gf_M_sc().data(), state.gf_M_sc().data() + state.gf_M_sc().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_M_sc_distance = std::distance( state.gf_M_sc().data(), max_M_sc_ptr ); 
   std::cout  << " Max M sc " << *max_M_sc_ptr << " at idx " << state.gf_M_sc().get_idx( max_M_sc_distance ) << std::endl; 

   const dcomplex* max_M_d_ptr = max_element( state.gf_M_d().data(), state.gf_M_d().data() + state.gf_M_d().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_M_d_distance = std::distance( state.gf_M_d().data(), max_M_d_ptr ); 
   std::cout  << " Max M d " << *max_M_d_ptr << " at idx " << state.gf_M_d().get_idx( max_M_d_distance ) << std::endl; 

   const dcomplex* max_M_m_ptr = max_element( state.gf_M_m().data(), state.gf_M_m().data() + state.gf_M_m().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_M_m_distance = std::distance( state.gf_M_m().data(), max_M_m_ptr ); 
   std::cout << " Max M m " << *max_M_m_ptr << " at idx " << state.gf_M_m().get_idx( max_M_m_distance ) << std::endl; 
#else
   const dcomplex* max_M_sc_ptr = max_element( state.gf_wM_sc_plus().data(), state.gf_wM_sc_plus().data() + state.gf_wM_sc_plus().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_M_sc_distance = std::distance( state.gf_wM_sc_plus().data(), max_M_sc_ptr ); 
   std::cout  << " Max wM sc " << *max_M_sc_ptr << " at idx " << state.gf_wM_sc_plus().get_idx( max_M_sc_distance ) << std::endl; 

   const dcomplex* max_M_d_ptr = max_element( state.gf_wM_d_plus().data(), state.gf_wM_d_plus().data() + state.gf_wM_d_plus().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_M_d_distance = std::distance( state.gf_wM_d_plus().data(), max_M_d_ptr ); 
   std::cout  << " Max wM d " << *max_M_d_ptr << " at idx " << state.gf_wM_d_plus().get_idx( max_M_d_distance ) << std::endl; 

   const dcomplex* max_M_m_ptr = max_element( state.gf_wM_m_plus().data(), state.gf_wM_m_plus().data() + state.gf_wM_m_plus().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_M_m_distance = std::distance( state.gf_wM_m_plus().data(), max_M_m_ptr ); 
   std::cout << " Max wM m " << *max_M_m_ptr << " at idx " << state.gf_wM_m_plus().get_idx( max_M_m_distance ) << std::endl; 
#endif


   const dcomplex* max_lambda_sc_ptr = max_element( state.gf_lambda_sc().data(), state.gf_lambda_sc().data() + state.gf_lambda_sc().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_lambda_sc_distance = std::distance( state.gf_lambda_sc().data(), max_lambda_sc_ptr ); 
   std::cout  << " Max lambda sc " << *max_lambda_sc_ptr << " at idx " << state.gf_lambda_sc().get_idx( max_lambda_sc_distance ) << std::endl; 

   const dcomplex* max_lambda_d_ptr = max_element( state.gf_lambda_d().data(), state.gf_lambda_d().data() + state.gf_lambda_d().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_lambda_d_distance = std::distance( state.gf_lambda_d().data(), max_lambda_d_ptr ); 
   std::cout  << " Max lambda d " << *max_lambda_d_ptr << " at idx " << state.gf_lambda_d().get_idx( max_lambda_d_distance ) << std::endl; 

   const dcomplex* max_lambda_m_ptr = max_element( state.gf_lambda_m().data(), state.gf_lambda_m().data() + state.gf_lambda_m().num_elements(), []( const dcomplex& a, const dcomplex& b )->bool{ return std::abs(a) < std::abs(b); } ); 
   int max_lambda_m_distance = std::distance( state.gf_lambda_m().data(), max_lambda_m_ptr ); 
   std::cout << " Max lambda m " << *max_lambda_m_ptr << " at idx " << state.gf_lambda_m().get_idx( max_lambda_m_distance ) << std::endl; 

   double max_w_lambda = std::max(std::max(std::max(std::max( std::max( std::abs(*max_sc_ptr), std::abs(*max_d_ptr) ), std::abs(*max_m_ptr) ), std::abs(*max_lambda_sc_ptr) ), std::abs(*max_lambda_d_ptr) ), std::abs(*max_lambda_m_ptr) );
   double max_w_lambda_M = std::max(std::max( std::max( std::abs(*max_M_sc_ptr), max_w_lambda), abs(*max_M_d_ptr)), abs(*max_M_m_ptr) );
   return max_w_lambda_M;
}

