
/*******************************************************************************************//** @file
 *  		
 * 	file: 		observer.h
 * 	contents:  	functor to track observables during the flow
 * 
 ****************************************************************************************************/


#pragma once
#include <frg/sbe/state.h>

#include <frg/observables_common.h>
#include <models/concrete_available_models.h>

template <typename Model, typename state_t > 
class observer_frg_t											///< Functor to calculate observables from the state and stores them during the flow/in between the integration steps. It returns false if a quantity in the vertex in the state diverges.
{
 public:
    bool operator() ( const state_t& state_vec, const double t, bool is_final = false ); 					        
    observer_frg_t(observables_frg_common_t<Model, state_t> &observables);

    // SBE
    double get_max_cpl( const state_t& state_vec ); 	///< Find the coupling with the largest absolute value. Returns the absolute value and the position
//sbe end

    
    observables_frg_common_t<Model, state_t> &m_observables;
    
};

#include "../src/frg/observer.tpp"
