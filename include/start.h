#pragma once

#include <boost/numeric/odeint/integrate/integrate_adaptive.hpp>
#include <boost/numeric/odeint/integrate/detail/integrate_adaptive.hpp>
#include <boost/numeric/odeint.hpp>

#include <models/concrete_available_models.h>

#include <frg/sbe/observables.h>
#include <frg/observer.h>

#include <frg/sbe/state.h>
#include <frg/sbe/state_bosonised_M.h>
#include <frg/sbe/rhs_1lfrg.h>
#include <frg/sbe/rhs_1lfrg_with_f.h>
#include <frg/sbe/symmetries.h>
#include <frg/sbe/symmetries_bosonised_M.h>
#include <frg/sbe/interpolators.h>


#include <frg/sbe/rhs_mfrg.h>
#include <frg/sbe/rhs_mfrg_with_f.h>
#include <frg/sbe/symmetries_mfrg.h>

#include <frg/sbe/postprocessing/symmetries.h>
#include <frg/sbe/postprocessing/symmetries.h>

#include <self_consistent/sbe/rhs.h>

#include <file_IO/file_IO_sbe.h>
#include <file_IO/multifile_IO_sbe.h>


using Model = THE_MODEL;

#ifdef BOSONISE_M
    using state_t = state_frg_sbe_bosonised_M_t<Model>;
#else
    using state_t = state_frg_sbe_t<Model>; 
#endif

#ifdef FLOW_EQUATION_METHOD
    #ifdef MULTILOOP
        #ifdef INCLUDE_FREE_ENERGY
            using rhs_t = rhs_sbe_mfrg_with_F_t<Model, state_t>;
        #else
            using rhs_t = rhs_sbe_mfrg_t<Model, state_t>;
        #endif
        using Symmetries = SymmetriesmfRGSBE<Model>;
        using Interpolators = Interpolators1lfRGSBE<Model>;
    #else
        #ifdef INCLUDE_FREE_ENERGY
            using rhs_t = rhs_sbe_1lfrg_with_F_t<Model, state_t>;
        #else
            using rhs_t = rhs_sbe_1lfrg_t<Model, state_t>;
        #endif
        #ifdef BOSONISE_M
            using Symmetries = SymmetriesfRGSBEBosonised_M<Model>;
        #else 
            using Symmetries = Symmetries1lfRGSBE<Model>;
        #endif
        using Interpolators = Interpolators1lfRGSBE<Model>;
    #endif
#elif defined(SELF_CONSISTENT_METHOD)
    using rhs_t = rhs_sbe_selfconsistent_t<Model, state_t>;
    using Symmetries = Symmetries1lfRGSBE<Model>;
    using Interpolators = Interpolators1lfRGSBE<Model>;
#endif

using observables_t = observables_frg_sbe_t<Model, state_t>;
using PostprocessingSymmetries = SymmetriesfRGSBEPostprocessing<Model>;

using FileIO = FileIOSBE<Model, state_t>;
#ifdef MULTIFILE_OUTPUT
    using MultiFileIO = MultiFileIOSBE<Model, state_t>;
#endif


void init_state(state_t &state_vec);

/*
template <typename Model>
void init_observables(observables_t<Model> &observables);
*/

bool integrate_ode(state_t &state_vec, rhs_t &rhs, observer_frg_t<Model, state_t> &observer, double start_time, double final_time, double initial_timestep);


bool solve_selfconsistently(state_t &state_vec, rhs_t &rhs, observer_frg_t<Model, state_t> &observer);


void solve_frg(const resume_details_t &resume_details);

namespace boost {
    namespace numeric {
	namespace odeint {
	    namespace detail {
		/**
		 * Changed integrate adaptive to consider abort when monitored couplings exceed a maximum values
		 * in the vicinity of a specific phase transition
		 */
		template< class Stepper, class Time, class Observer >
		    bool integrate_adaptive_check(
						  Stepper stepper, rhs_t system, state_t &start_state,
						  Time &start_time, Time end_time, Time &dt,
						  Observer observer           /**< controlled_stepper_tag */
						  );
	    } // Namespace detail
	    
	    template< class Stepper, class Time, class Observer >
		bool integrate_adaptive_check(
					      Stepper stepper, rhs_t system, state_t &start_state,
					      Time start_time, Time end_time, Time dt,
					      Observer observer )
	    {
		return detail::integrate_adaptive_check(
							stepper, system, start_state,
							start_time, end_time, dt,
							observer ); //, typename Stepper::stepper_category() );
	    }

	} // Namespace odeint
    } // Namespace numeric
} // Namespace boost
