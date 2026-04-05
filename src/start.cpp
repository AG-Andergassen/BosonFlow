#include <start.h>
#include <params_technical.h>
#include <params_physical.h>
#include <runtime_config.h>

#include <tu_projections.h>

#include <frg/flows.h>

#include <frg/observables_common.h>
#include <frg/observer.h>

#include <frequencies/schemes.h>


using namespace boost::numeric::odeint;


int main ( int argc, char** argv){
    
    // Parse command-line arguments early
    RuntimeConfig::ParseCommandLine(argc, argv);
    
    std::cout << std::scientific << std::setprecision(4);     /**< Set output format */

    Model::Init();

#ifdef LOGARITHMIC_GRIDS
    FrequencyDependenceScheme<Model>::Init(FrequencyDependenceSchemeName::LogarithmicGrids);
#else
    FrequencyDependenceScheme<Model>::Init(FrequencyDependenceSchemeName::UniformGrids);
#endif

    // temporary
    //Model::Test();
    
    Symmetries::Init();
    
    Interpolators::Init();

    TUProjections<Model>::Init();

    G_FlowSchemeName G_flow_name = G_FlowSchemeName::None;
    U_FlowSchemeName U_flow_name = U_FlowSchemeName::None;

    // pick G cutoff
#if defined(INTRP_FLOW) || defined(INVERSE_INTRP_FLOW) || defined(INVERSE_INTRP_FLOW_OMEGA)
    #ifdef INTRP_FLOW
    G_flow_name = G_FlowSchemeName::Interpolating;
    #elif defined(INVERSE_INTRP_FLOW)
    G_flow_name = G_FlowSchemeName::InverseInterpolating;
    #elif defined(INVERSE_INTRP_FLOW_OMEGA)
    G_flow_name = G_FlowSchemeName::InverseInterpolatingOmega;
    #endif
    
    // load DMFT data using MultiFileIO
    fRGFlowScheme<Model>::Set_Ginit_ForInterpolatingFlow( FileIO::read_initial_Ginit_from_DMFT_data(INPUT_REFERENCE_DATA_DIRECTORY) );
#elif defined(INT_FLOW)
    G_flow_name = G_FlowSchemeName::Interaction;
#elif defined(TEMP_FLOW)
    G_flow_name = G_FlowSchemeName::Temperature;
#elif defined(OMEGA_FLOW)
    G_flow_name = G_FlowSchemeName::Omega;
#endif

    // pick U cutoff
#ifdef FREE_U_FLOW
    U_flow_name = U_FlowSchemeName::Free;
#elif defined(ACTUAL_U_FLOW)
    U_flow_name = U_FlowSchemeName::Multiplicative;
#endif

    fRGFlowScheme<Model>::UseFlow(G_flow_name, U_flow_name);

    PostprocessingSymmetries::Init();
    observables_t::SetToTrackAllAvailableObservables();

    resume_details_t resume_details;
    const auto &chosen_flow = fRGFlowScheme<Model>::ChosenFlowParametrizationInfo();
    resume_details.start_time = chosen_flow.t_start;
    resume_details.final_time = chosen_flow.t_final;
    resume_details.init_timestep = chosen_flow.init_t_step;
#ifdef MULTIFILE_OUTPUT
    MultiFileIO::Init(resume_details);
#endif
    
    solve_frg(resume_details);

    Symmetries::CleanUp();

    return EXIT_SUCCESS;
}
// todo: add resume details as argument
void init_state(state_t &state)
{
#if defined(INTRP_FLOW) || defined(INVERSE_INTRP_FLOW) || defined(INVERSE_INTRP_FLOW_OMEGA)
    state.init_bare();
    std::cout << " ... Initialising to DMFT state" << std::endl;
    FileIO::read_initial_state_from_DMFT_data(INPUT_REFERENCE_DATA_DIRECTORY, state);
#else
    std::cout << " ... Initialising to bare values" << std::endl;
    state.init_bare();
#endif

#ifdef AUTORESUME_CALCULATION
    if (MultiFileIO::CurrentScaleIdx() >= 0){ // return to zero
	FileIO::read_state_from_checkpoint_file(MultiFileIO::OutputDirectory() + "/" + std::to_string(MultiFileIO::CurrentScaleIdx()) + ".h5", state); 
	return;
    }
#endif

}

void solve_frg(const resume_details_t &resume_details)
{
#ifdef MULTIFILE_OUTPUT
    MultiFileIO::WritePropertiesFile();
#endif

    state_t state_vec;
    init_state(state_vec);

    rhs_t rhs;

    observables_t observables;

    // the "observer" is a functor which updates the observables every integration step
    observer_frg_t<Model, state_t> observer(observables);
    bool success;

#ifdef MULTIFILE_OUTPUT
    MultiFileIO::WriteCurrentScaleFile(state_vec, observer.m_observables.m_name_data_map);
#endif

#ifdef FLOW_EQUATION_METHOD
    if (!resume_details.is_finished)
	success = integrate_ode(state_vec, rhs, observer, resume_details.start_time, resume_details.final_time, resume_details.init_timestep);
    else{
      std::cout << "Previous run is completed. Will redo post-processing calculations.." << std::endl;
      bool is_final = true; // todo: handle properly
      observer( state_vec, resume_details.start_time, is_final );

    }
#elif defined(SELF_CONSISTENT_METHOD)
    success = solve_selfconsistently(state_vec, rhs, observer);
#endif

    std::string to_be_appended_to_file_name = "";
    if( !success ) 
	to_be_appended_to_file_name = "_DIVERGENT";

#ifdef MULTIFILE_OUTPUT
    MultiFileIO::WriteFinalFile(state_vec, observables.m_name_data_map, to_be_appended_to_file_name);
#else
    FileIO::write_final_file(state_vec, observables.m_name_data_map, to_be_appended_to_file_name);
#endif

}

/*
void init_observables(observables_t<Model> & observables)
{
    // if there's a file input, read observables from it
#ifdef READIN
    std::cout << "Reading observables from file:" << INPUT_FILE_LOC + INPUT_FILE_NAME << std::endl;
    FileIO::read_observables_from_checkpoint_file(INPUT_FILE_LOC + INPUT_FILE_NAME, observables); 
#endif    
}
*/


// Function to adjust the proposed time step to avoid skipping mandatory steps
double adjust_timestep(const std::vector<double>& mandatory_steps, double current_time, double proposed_step, double dt) {
    // If there are no mandatory steps, return the proposed step as-is
    if (mandatory_steps.empty()) {
        return proposed_step;
    }
    std::cout <<"Remaining number of mandatory steps: " << mandatory_steps.size() << std::endl;
    if (std::abs(proposed_step - current_time) < 1e-5)
	return proposed_step;
    
    // Determine the direction of integration
    bool is_forward = (dt > 0);

    // Find the next mandatory step based on the direction of integration
    double next_mandatory_step;
    if (is_forward) {
        // Forward integration: find the first mandatory step > current_time
        auto it = std::upper_bound(mandatory_steps.begin(), mandatory_steps.end(), current_time);
        if (it == mandatory_steps.end()) {
            // No mandatory steps left, return the proposed step
            return proposed_step;
        }
        next_mandatory_step = *it;
    } else {
        // Backward integration: find the first mandatory step < current_time
        auto it = std::lower_bound(mandatory_steps.begin(), mandatory_steps.end(), current_time);
        if (it == mandatory_steps.begin()) {
            // No mandatory steps left, return the proposed step
            return proposed_step;
        }
        next_mandatory_step = *(it - 1);
    }

    // Ensure the solver doesn't skip the next mandatory step
    if (is_forward) {
        return std::min(proposed_step, next_mandatory_step );
    } else {
        return std::max(proposed_step, next_mandatory_step );
    }
}


#ifdef FLOW_EQUATION_METHOD 
bool integrate_ode(state_t &state_vec, rhs_t &rhs, observer_frg_t<Model, state_t> &observer, double start_time, double final_time, double initial_timestep){
    std::cout << " Initializing ODE solver... " << std::endl << std::endl;
 
#ifdef AUTORESUME_CALCULATION
    //get two scales, to get the INIT_T_STEP_SIZE. Get also T_START.
#endif

#ifndef DEBUG_EULER 
    // Integrate ODE with adaptive runge_kutta_cash_karp54:   
    typedef ERR_STEPPER< state_t, double, state_t, double, vector_space_algebra > error_stepper_t;       /**< ERR_STEPPER of the ODE solver. @see params_technical.h */
    std::cout << " Starting ODE solver with " << ERR_STEPPER_STRING  << std::endl;
    bool success = integrate_adaptive_check( make_controlled< error_stepper_t >( ERR_ABS, ERR_REL ), rhs, state_vec, start_time, final_time, initial_timestep, observer);

#else
   
    // Only for debugging, integrate with constant step Euler method. Use FLAG: DEBUG_EULER
    typedef ERR_STEPPER< state_t, double, state_t, double, vector_space_algebra > stepper_t;
    stepper_t stepper; 
    int steps = integrate_const( stepper, rhs, state_vec, start_time, final_time, (final_time-start_time)/NSTEPS, observer); 
    bool success = true; // note: it doesn't stop at divergences, as well as misses several important adjustments to the self-energy 
#endif // DEBUG_EULER

    return success;
}



template< class Stepper, class Time, class Observer >
bool boost::numeric::odeint::detail::integrate_adaptive_check(
    Stepper stepper, rhs_t rhs, state_t &state,
    Time &start_time, Time end_time, Time &dt,
    Observer observer           /**< controlled_stepper_tag */
			      )
{
    typename odeint::unwrap_reference< Observer >::type &obs = observer;            /**< Observer used to track flow observables */

    const size_t max_attempts = 1000;                                               /**< controlling number of attempt for the adaptive integration */
    const char *error_string = "Integrate adaptive : Maximal number of iterations reached. A step size could not be found.";
    size_t count = 0;
    size_t AorB  = 0;
    double initial_filling = 0.0;

    if (RuntimeConfig::HAS_FILLING_OVERRIDE){
        initial_filling = RuntimeConfig::FILLING_OVERRIDE;
        state.adjust_chemical_potential_shift(initial_filling, start_time);
    }
    else {
        initial_filling = state.eval_filling(start_time);
    }

    std::cout << "init filling:"  << initial_filling << std::endl;

    while( less_with_sign( start_time, end_time, dt ) )
    {                            /**< @see boost/numeric/odeint/integrate/detail/less_with_sign.hpp */

	// compute flow observables

    // adjustments to self-energy
#ifdef FIX_FILLING
	double delta_mu_before = state.gf_delta_mu();

	state.adjust_chemical_potential_shift(initial_filling, start_time);

	state.m_d_delta_mu_over_dt = (delta_mu_before - state.gf_delta_mu())/dt;
#endif

	bool is_max_coupling_finite = obs(state, start_time);

        // write current results as output file. How often this is done can be specified in params_technical.h, updatenum. Use flag OUT_UPDATE in makefile. Allows to continue calculations later on.
	if( (count)%UPDATER_FREQ == 0 and count >= UPDATER_START )
        {
#ifdef MULTIFILE_OUTPUT
	    MultiFileIO::IncrementScaleIdx();
	    MultiFileIO::WritePropertiesFile();
	    MultiFileIO::WriteCurrentScaleFile(state, observer.m_observables.m_name_data_map);
#else
#ifdef OUT_UPDATE        
	    FileIO::write_checkpoint_file( state, observer.m_observables.m_name_data_map, count, start_time, end_time, dt, AorB); 
           ++AorB;
#endif // ends OUT_UPDATE
#endif
        }
        if( ! is_max_coupling_finite ){                                      /**< Check that couplings are not exploding before performing the ODE step */
            std::cout << std::endl << " !!!! Vertex exceeded MAX_COUPLING at scale " << start_time << ", ODE solver stopped !!!! " << std::endl << std::endl;
	    bool is_final = true; // todo: handle properly
	    obs( state, start_time, is_final );

            return 0;                                      /**< ODE integration aborted! */
            
        }
        if( less_with_sign( end_time, start_time + dt, dt ) ){                      /**< This condition allows to not exceed the integation upper limit -> \f$ \Lambda_{fin} \f$ = end_time */
            dt = end_time - start_time;
        }


	Time prev_time = start_time;

        size_t trials = 0;
        controlled_step_result res = success;
	
	Time adjusted_start_time_plus_dt = adjust_timestep(MANDATORY_TIMESTEP_PERCENTAGES, double(start_time), double(start_time) + double(dt), double(dt));
	std::cout << "Old dt: "<< dt << std::endl;
	std::cout << "pre_start_time: "<< start_time << " adjusted time: " << adjusted_start_time_plus_dt << std::endl; 
	if (std::abs(adjusted_start_time_plus_dt - start_time) > 1e-5){
	    dt = adjusted_start_time_plus_dt - start_time;
	    std::cout << "New dt: " << dt << " start time: "<< start_time << " adjusted time: " << adjusted_start_time_plus_dt << std::endl;
	}
	do {
            std::cout << "trial is " << trials << std::endl;
#ifdef PRECOMPUTE_STATE_PROJECTIONS
	    state.precalculate_projected_Ms_and_nablas(start_time);
#endif 
	    std::cout <<"OLD start_time and dt" << start_time << ", " << dt << std::endl;
            res = stepper.try_step( rhs, state, start_time, dt );          /**< attempts an integration step with the given dt. If error is too big, it will make dt smaller and return a result which is 'fail'. If error is within tolerance, it updates start_time to start_time + dt, and then updates dt to a new proposed timestep. @return success when the step has been found.  */
	    std::cout <<"NEW start_time and dt" << start_time << ", " << dt << std::endl;
	    /*if (res != fail){
		std::cout << "prev_time" << prev_time << std::endl;
		std::cout << "prev start_time: " << start_time << std::endl;
	    	start_time = adjust_timestep(MANDATORY_TIMESTEP_PERCENTAGES, prev_time, start_time, dt);
		std::cout << "new start_time: " << start_time << std::endl;
		}*/
	    
            ++trials;
        } while( ( res == fail ) && ( trials < max_attempts ) );
	
        if( trials == max_attempts ) 
	    throw std::overflow_error( error_string );
        ++count;

    }

#ifdef FIX_FILLING
        double delta_mu_before = state.gf_delta_mu();
	state.adjust_chemical_potential_shift(initial_filling, start_time);
	state.m_d_delta_mu_over_dt = (delta_mu_before - state.gf_delta_mu())/dt;
#endif

    
    bool is_final = true;

    obs( state, start_time, is_final );

    if (fRGFlowScheme<Model>::G_ChosenFlowName() == G_FlowSchemeName::Temperature){
	state.normalise_Sig_Temperature_Flow(start_time);
	std::cout << "Self-energy was normalized" << std::endl;
    }

    
    std::cout << " ODE solver finished with " << count << " integration steps." << std::endl << std::endl;
    
    return 1;

} // integrate_adaptive_check

#endif // FLOW_EQUATION_METHOD


#ifdef SELF_CONSISTENT_METHOD
bool solve_selfconsistently(state_t &state_vec, rhs_t &rhs, observer_frg_t<Model, state_t>&observer)
{
    state_t old_state = state_vec;
    const double t_final = fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_final;
    
    const double anderson_parameter = 0.9;

    int vertex_iterations_counter = 0;
    int selfenergy_iterations_counter = 0;
    double initial_filling = 0.0;

    if (RuntimeConfig::HAS_FILLING_OVERRIDE){
        initial_filling = RuntimeConfig::FILLING_OVERRIDE;
        state_vec.adjust_chemical_potential_shift(initial_filling, t_final);
    }
    else {
        initial_filling = state_vec.eval_filling(t_final);
    }

    
    bool is_max_coupling_finite = observer(state_vec, t_final);
    double t = t_final;
    
    //for (float shift_filling = 0; shift_filling <= 0.100001; shift_filling += 0.005){
    float shift_filling = 0.0;
    {
      //for (t = 0.1; t <= 1.0; t += 0.01)
	{
    state_vec.adjust_chemical_potential_shift(initial_filling - shift_filling, t);
    double the_norm = 1e50;
    std::cout << "Current t: "<< t << std::endl;
    do {
        rhs.update_bubbles_and_G(state_vec, t);
	//do {
	    vertex_iterations_counter ++;

            #ifdef PRECOMPUTE_STATE_PROJECTIONS
	        state_vec.precalculate_projected_Ms_and_nablas(t);
            #endif 

	    old_state = anderson_parameter*old_state + (1.0 - anderson_parameter)*state_vec;
            
	    rhs.vertex(old_state, state_vec, the_norm);

	    is_max_coupling_finite = observer(state_vec, t);

	    // write current results as output file. How often this is done can be specified in params_technical.h, updatenum. Use flag OUT_UPDATE in makefile. Allows to continue calculations later on.
#ifdef MULTIFILE_OUTPUT
	    MultiFileIO::IncrementScaleIdx();
	    MultiFileIO::WritePropertiesFile();
	    //MultiFileIO::WriteCurrentScaleFile(state_vec, observer.m_observables.m_name_data_map);
#endif
	    
	    std::cout << "Iterations: Vertex - " << vertex_iterations_counter << " Selfenergy - " << selfenergy_iterations_counter << std::endl;
	    std::cout << "Norm: " << norm(old_state - state_vec) << std::endl;

	    //} while ( norm(old_state - state_vec) > 1e-6 );

	selfenergy_iterations_counter ++;
        rhs.selfenergy(old_state, state_vec.gf_Sig());

#ifdef FIX_FILLING
	state_vec.adjust_chemical_potential_shift(initial_filling-shift_filling, t);
#endif

	std::cout << "Iterations: Selfenergy - " << selfenergy_iterations_counter << std::endl;
	std::cout << "Norm: " << norm(old_state - state_vec) << std::endl;

    } while (norm(old_state - state_vec)/norm(state_vec) > 1e-5);
    }
    }
    bool is_final = true;

    is_max_coupling_finite = observer(state_vec, t_final, is_final);
    
    double solution_norm = norm(state_vec);

    if (std::isnan(solution_norm) || std::isinf(solution_norm)){
	return false; 
    }else{
	std::cout << "Converged!" << std::endl;
	return true;
    }
    
    
}

#endif
