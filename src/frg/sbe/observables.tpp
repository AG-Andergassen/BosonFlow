#include <frg/sbe/postprocessing/rhs.h>
#include <frg/sbe/postprocessing/state.h>
#include <frg/sbe/rhs_mfrg.h>


template <typename Model, typename state_t >
void observables_frg_sbe_t<Model, state_t>::SetToTrackAllAvailableObservables()
{
    observables_frg_common_t<Model, state_t>::SetToTrackAllAvailableObservables();

    observables_frg_common_t<Model, state_t>::ObservablesListToTrack().push_back({"Free_Energy", FlowObservableType::COMPLEX, 0});

    observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack().push_back({
	    "w_info",
		{
		    {"Leading_name", FlowObservableType::STRING, 0},
		    {"Max_sc", FlowObservableType::COMPLEX, 0},
		    {"Max_d", FlowObservableType::COMPLEX, 0},
		    {"Max_m", FlowObservableType::COMPLEX, 0},
		    {"Max_sc_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_d_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_m_idx", FlowObservableType::GF_MULTIINDEX, 0}
		}
	});

    observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack().push_back({
	    "S_Wave_Susc_info",
		{
		    {"Susc_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_w_t<Model>::base_t::ndims},
		    {"Susc_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_w_t<Model>::base_t::ndims},
		    {"Susc_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_w_t<Model>::base_t::ndims},

		    {"Leading_name", FlowObservableType::STRING, 0},
		    {"Max_sc", FlowObservableType::COMPLEX, 0},
		    {"Max_d", FlowObservableType::COMPLEX, 0},
		    {"Max_m", FlowObservableType::COMPLEX, 0},
		    {"Max_sc_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_d_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_m_idx", FlowObservableType::GF_MULTIINDEX, 0}
		}
	});

    observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack().push_back({
	    "Postprocessing_Susc_info",
		{
		    {"Susc_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},

#ifdef SPLIT_SUSC_CONTRIBUTIONS
		    {"Susc_sc_vertex_contribution", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_bubble_contribution", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_bubble_contribution", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_bubble_contribution", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},

		    {"Susc_sc_vertex_contribution_from_bare_vertex", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_bare_vertex", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_bare_vertex", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_vertex_contribution_from_nabla_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_nabla_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_nabla_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_vertex_contribution_from_M_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_M_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_M_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_vertex_contribution_from_sc_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_sc_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_sc_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_vertex_contribution_from_nabla_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_nabla_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_nabla_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_vertex_contribution_from_M_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_M_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_M_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_vertex_contribution_from_m_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_m_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_m_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_vertex_contribution_from_nabla_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_nabla_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_nabla_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_vertex_contribution_from_M_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_M_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_M_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_sc_vertex_contribution_from_d_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_d_vertex_contribution_from_d_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
		    {"Susc_m_vertex_contribution_from_d_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_susc_t<Model>::base_t::ndims},
#endif

		    {"Leading_name", FlowObservableType::STRING, 0},
		    {"Max_sc", FlowObservableType::COMPLEX, 0},
		    {"Max_d", FlowObservableType::COMPLEX, 0},
		    {"Max_m", FlowObservableType::COMPLEX, 0},
		    {"Max_sc_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_d_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_m_idx", FlowObservableType::GF_MULTIINDEX, 0}
		}
	});


    observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack().push_back({
	    "Postprocessing_Polarisation_info",
		{
		    {"Polarisation_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_sc_minus_bubble", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_d_minus_bubble", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_m_minus_bubble", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_sc_bubble", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_d_bubble", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_m_bubble", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_sc_contribution_from_nabla_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_m_contribution_from_nabla_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},
		    {"Polarisation_d_contribution_from_nabla_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_polarisation_t<Model>::base_t::ndims},

		    {"Leading_name", FlowObservableType::STRING, 0},
		    {"Max_sc", FlowObservableType::COMPLEX, 0},
		    {"Max_d", FlowObservableType::COMPLEX, 0},
		    {"Max_m", FlowObservableType::COMPLEX, 0},
		    {"Max_sc_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_d_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_m_idx", FlowObservableType::GF_MULTIINDEX, 0}
		}
	});


    observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack().push_back({
	    "Postprocessing_Lambda_info",
		{
		    {"Lambda_sc_minus_1", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_minus_1", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_minus_1", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
#ifdef SPLIT_LAMBDA_CONTRIBUTIONS
		    {"Lambda_sc_contribution_from_bare_vertex", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_bare_vertex", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_bare_vertex", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_nabla_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_nabla_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_nabla_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_sc_contribution_from_sc_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_sc_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_sc_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_sc_contribution_from_nabla_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_sc_contribution_from_nabla_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_nabla_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_sc_contribution_from_m_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_m_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_m_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_nabla_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_sc_contribution_from_nabla_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_nabla_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},  
		    {"Lambda_sc_contribution_from_d_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_d_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_d_double_counting", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_sc_contribution_from_M_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_sc_contribution_from_M_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_sc_contribution_from_M_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_M_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_M_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_d_contribution_from_M_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_M_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_M_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    {"Lambda_m_contribution_from_M_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_static_lambda_t<Model>::base_t::ndims},
		    
#endif

		    {"Leading_name", FlowObservableType::STRING, 0},
		    {"Max_sc", FlowObservableType::COMPLEX, 0},
		    {"Max_d", FlowObservableType::COMPLEX, 0},
		    {"Max_m", FlowObservableType::COMPLEX, 0},
		    {"Max_sc_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_d_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_m_idx", FlowObservableType::GF_MULTIINDEX, 0}
		}
	});


#ifdef CALCULATE_SIG_SDE_INTEGRANDS 
    observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack().push_back({
	    "Postprocessing_SDE_Integrands",
		{
		    {"Sig_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_Sig_SDE_integrand_t<Model>::base_t::ndims},
		    {"Sig_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_Sig_SDE_integrand_t<Model>::base_t::ndims},
		    {"Sig_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_Sig_SDE_integrand_t<Model>::base_t::ndims},
		    {"Leading_name", FlowObservableType::STRING, 0},
		    {"Max_sc", FlowObservableType::COMPLEX, 0},
		    {"Max_d", FlowObservableType::COMPLEX, 0},
		    {"Max_m", FlowObservableType::COMPLEX, 0},
		    {"Max_sc_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_d_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_m_idx", FlowObservableType::GF_MULTIINDEX, 0}
		}
	});
#endif

    
    /*    
        observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack().push_back({
	    "vertex_info",
		{
		    {"vertex_sc", FlowObservableType::COMPLEX_GF_CONTAINER, gf_M_t<Model>::base_t::ndims},
		    {"vertex_d", FlowObservableType::COMPLEX_GF_CONTAINER, gf_M_t<Model>::base_t::ndims},
		    {"vertex_m", FlowObservableType::COMPLEX_GF_CONTAINER, gf_M_t<Model>::base_t::ndims},

		    {"Leading_name", FlowObservableType::STRING, 0},
		    {"Max_sc", FlowObservableType::COMPLEX, 0},
		    {"Max_d", FlowObservableType::COMPLEX, 0},
		    {"Max_m", FlowObservableType::COMPLEX, 0},
		    {"Max_sc_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_d_idx", FlowObservableType::GF_MULTIINDEX, 0},
		    {"Max_m_idx", FlowObservableType::GF_MULTIINDEX, 0}
		}
		});
    */


#ifdef MULTILOOP
    observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack().push_back({
	    "multiloop_info",
		{
		    {"Sig_Iterations_Count", FlowObservableType::DOUBLE, 0},
		    {"Total_Vertex_Corrections_Count", FlowObservableType::DOUBLE, 0},
		    {"Last_Vertex_Corrections_Count", FlowObservableType::DOUBLE, 0},
		    {"Max_Vertex_Corrections_Count", FlowObservableType::DOUBLE, 0}
		}
	});

#ifdef STORE_MULTILOOP_VERTEX_CORRECTIONS

    if (SELFENERGY_ITERATIONS_MAX == 1)
    for (int i = 0; i < 1+EXTRA_LOOP_NUM; i ++){      
      observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack().push_back({
	  "multiloop_corrections_" + std::to_string(i+1),
	  {
	    {"w_sc_dot", FlowObservableType::COMPLEX_GF_CONTAINER, gf_w_t<Model>::base_t::ndims},
	    {"w_d_dot", FlowObservableType::COMPLEX_GF_CONTAINER, gf_w_t<Model>::base_t::ndims},
	    {"w_m_dot", FlowObservableType::COMPLEX_GF_CONTAINER, gf_w_t<Model>::base_t::ndims},
	    {"lambda_sc_dot", FlowObservableType::COMPLEX_GF_CONTAINER, gf_lambda_t<Model>::base_t::ndims},
	    {"lambda_d_dot", FlowObservableType::COMPLEX_GF_CONTAINER, gf_lambda_t<Model>::base_t::ndims},
	    {"lambda_m_dot", FlowObservableType::COMPLEX_GF_CONTAINER, gf_lambda_t<Model>::base_t::ndims}
#ifdef RESTFUNC
	    ,{"M_sc_dot", FlowObservableType::COMPLEX_GF_CONTAINER, gf_M_t<Model>::base_t::ndims},
	    {"M_d_dot", FlowObservableType::COMPLEX_GF_CONTAINER, gf_M_t<Model>::base_t::ndims}
	    {"M_m_dot", FlowObservableType::COMPLEX_GF_CONTAINER, gf_M_t<Model>::base_t::ndims}
#endif
	  }
	});

      
    }
    #endif 
#endif

    
}

template <typename Model, typename state_t >
void observables_frg_sbe_t<Model, state_t>::update(const state_t &state, double t, bool is_final)
{
    observables_frg_common_t<Model, state_t >::update(state, t, is_final);
    
    state_postproc_t<Model> postproc_state;
    rhs_postproc_sbe_t<Model, state_t> postproc_rhs;

    bool include_postprocessing = is_final;

#ifdef POSTPROCESSING_EVERY_STEP
    include_postprocessing = true;
#endif

    if (include_postprocessing)
	postproc_rhs(state, postproc_state, t);

    
    for (auto &group_observables: observables_frg_common_t<Model, state_t >::GroupObservablesListToTrack()){
	std::string group_name = std::get<0>(group_observables);

        if (group_name == "w_info"){
	    auto w_info = this->Eval_max_info(state.gf_w_sc(), state.gf_w_d(), state.gf_w_m());

	    this->append_observable_string(group_name+"/Leading_name", std::get<0>(w_info));
	    this->append_observable_complex(group_name+"/Max_sc", std::get<1>(w_info));
	    this->append_observable_complex(group_name+"/Max_d", std::get<2>(w_info));
	    this->append_observable_complex(group_name+"/Max_m", std::get<3>(w_info));

            this->append_observable_multiindex(group_name+"/Max_sc_idx", std::get<4>(w_info));
	    this->append_observable_multiindex(group_name+"/Max_d_idx", std::get<5>(w_info));
	    this->append_observable_multiindex(group_name+"/Max_m_idx", std::get<6>(w_info));
	}

	if (group_name == "vertex_info"){
	    gf_M_t<Model> vertex_sc(0, 1), vertex_m(0, 1), vertex_d(0, 1);
	    // todo: fix U_FLOW cutoff
	    vertex_sc.init([&state, t](idx_M_t<Model> idx){
		    return state.lambda_sc(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m))*
			state.w_sc(idx(IM::W), idx(IM::K), t)*
			state.lambda_sc(idx(IM::W), idx(IM::K), idx(IM::wp), idx(IM::mp))
#if !defined(SBEa_APPROXIMATION)
			+ state.M_sc(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp))
#endif
		      - fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*state.vertex_sc_bare(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp), t);
		});

	    vertex_d.init([&state, t](idx_M_t<Model> idx){
		    return state.lambda_d(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m))*
			state.w_d(idx(IM::W), idx(IM::K), t)*
			state.lambda_d(idx(IM::W), idx(IM::K), idx(IM::wp), idx(IM::mp))
#if !defined(SBEa_APPROXIMATION)
			+ state.M_d(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp))
#endif
		      - fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*state.vertex_d_bare(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp), t);
		});

	    vertex_m.init([&state, t](idx_M_t<Model> idx){
		    return state.lambda_m(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m))*
			state.w_m(idx(IM::W), idx(IM::K), t)*
			state.lambda_m(idx(IM::W), idx(IM::K), idx(IM::wp), idx(IM::mp))
#if !defined(SBEa_APPROXIMATION)
			+ state.M_m(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp))
#endif
		      - fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*state.vertex_m_bare(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp), t);
		});


	    this->append_observable_complex_gf_container(group_name + "/vertex_sc", vertex_sc);
	    this->append_observable_complex_gf_container(group_name + "/vertex_d", vertex_d);
	    this->append_observable_complex_gf_container(group_name + "/vertex_m", vertex_m);

	    auto vertex_info = this->Eval_max_info(vertex_sc, vertex_d, vertex_m);

	    this->append_observable_string(group_name+"/Leading_name", std::get<0>(vertex_info));
	    this->append_observable_complex(group_name+"/Max_sc", std::get<1>(vertex_info));
	    this->append_observable_complex(group_name+"/Max_d", std::get<2>(vertex_info));
	    this->append_observable_complex(group_name+"/Max_m", std::get<3>(vertex_info));

	    //	    std::cout << "Max Phi_sc " << std::get<1>(vertex_info) << std::endl;
	    //std::cout << "Max Phi_m " << -std::get<3>(vertex_info) << std::endl;


            this->append_observable_multiindex(group_name+"/Max_sc_idx", std::get<4>(vertex_info));
	    this->append_observable_multiindex(group_name+"/Max_d_idx", std::get<5>(vertex_info));
	    this->append_observable_multiindex(group_name+"/Max_m_idx", std::get<6>(vertex_info));


	}
	
	if (group_name == "S_Wave_Susc_info"){
	    // s-wave susceptibilities have the same type as the screened interactions, but we only consider one bosonic frequency (pos_bfreq_count = 0), namely iw = 0
	    gf_w_t<Model> susc_sc(0), susc_d(0), susc_m(0);

	    // note: should be scaled in the plotting script (left to the user)

	    susc_sc.init([&state, t](idx_w_t<Model> idx){
		    double U = fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::B_sc(idx(IW::W), idx(IW::K), 0, 0, 0, 0);
		    return (state.gf_w_sc()(idx) - U)/U/U * Model::MomentumGrid().ptr->get_volume();
		});

	    susc_d.init([&state, t](idx_w_t<Model> idx){
		    double U = fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::B_d(idx(IW::W), idx(IW::K), 0, 0, 0, 0);
		    return (state.gf_w_d()(idx) - U)/U/U * Model::MomentumGrid().ptr->get_volume();
		});

	    susc_m.init([&state, t](idx_w_t<Model> idx){
		    double U = fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::B_m(idx(IW::W), idx(IW::K), 0, 0, 0, 0);     
		    return (state.gf_w_m()(idx) - U)/U/U * Model::MomentumGrid().ptr->get_volume();
		});

	    this->append_observable_complex_gf_container(group_name + "/Susc_sc", susc_sc);
	    this->append_observable_complex_gf_container(group_name + "/Susc_d", susc_d);
	    this->append_observable_complex_gf_container(group_name + "/Susc_m", susc_m);

	    auto susc_info = this->Eval_max_info(susc_sc, susc_d, susc_m);

	    this->append_observable_string(group_name+"/Leading_name", std::get<0>(susc_info));
	    this->append_observable_complex(group_name+"/Max_sc", std::get<1>(susc_info));
	    this->append_observable_complex(group_name+"/Max_d", std::get<2>(susc_info));
	    this->append_observable_complex(group_name+"/Max_m", std::get<3>(susc_info));

            this->append_observable_multiindex(group_name+"/Max_sc_idx", std::get<4>(susc_info));
	    this->append_observable_multiindex(group_name+"/Max_d_idx", std::get<5>(susc_info));
	    this->append_observable_multiindex(group_name+"/Max_m_idx", std::get<6>(susc_info));
	}


	if (group_name == "multiloop_info"){
	  this->append_observable_double(group_name+"/Sig_Iterations_Count", double(rhs_sbe_mfrg_t<Model, state_t>::LatestNumberOfSelfEnergyIterations()));
	  this->append_observable_double(group_name+"/Total_Vertex_Corrections_Count", double(rhs_sbe_mfrg_t<Model, state_t>::LatestTotalNumberOfVertexCorrections()) );
	  this->append_observable_double(group_name+"/Last_Vertex_Corrections_Count", double(rhs_sbe_mfrg_t<Model, state_t>::LatestNumberOfVertexCorrections()) );
	  this->append_observable_double(group_name+"/Max_Vertex_Corrections_Count", double(rhs_sbe_mfrg_t<Model, state_t>::MaxNumberOfVertexCorrections()) );
	}

	for (int i = 0; i < 1+EXTRA_LOOP_NUM; i ++){
	  if (group_name == "multiloop_corrections_" + std::to_string(i+1)){
	    this->append_observable_complex_gf_container(group_name + "/w_sc_dot", rhs_sbe_mfrg_t<Model, state_t>::LatestStateDotMultiloopCorrectionsOrders()[i].gf_w_sc());
	    this->append_observable_complex_gf_container(group_name + "/w_d_dot", rhs_sbe_mfrg_t<Model, state_t>::LatestStateDotMultiloopCorrectionsOrders()[i].gf_w_d());
	    this->append_observable_complex_gf_container(group_name + "/w_m_dot", rhs_sbe_mfrg_t<Model, state_t>::LatestStateDotMultiloopCorrectionsOrders()[i].gf_w_m());
	    this->append_observable_complex_gf_container(group_name + "/lambda_sc_dot", rhs_sbe_mfrg_t<Model, state_t>::LatestStateDotMultiloopCorrectionsOrders()[i].gf_lambda_sc());
	    this->append_observable_complex_gf_container(group_name + "/lambda_d_dot", rhs_sbe_mfrg_t<Model, state_t>::LatestStateDotMultiloopCorrectionsOrders()[i].gf_lambda_d());
	    this->append_observable_complex_gf_container(group_name + "/lambda_m_dot", rhs_sbe_mfrg_t<Model, state_t>::LatestStateDotMultiloopCorrectionsOrders()[i].gf_lambda_m());

#if !defined(SBEa_APPROXIMATION)
	    this->append_observable_complex_gf_container(group_name + "/M_sc_dot", rhs_sbe_mfrg_t<Model, state_t>::LatestStateDotMultiloopCorrectionsOrders()[i].gf_M_sc());
	    this->append_observable_complex_gf_container(group_name + "/M_d_dot", rhs_sbe_mfrg_t<Model, state_t>::LatestStateDotMultiloopCorrectionsOrders()[i].gf_M_d());
	    this->append_observable_complex_gf_container(group_name + "/M_m_dot", rhs_sbe_mfrg_t<Model, state_t>::LatestStateDotMultiloopCorrectionsOrders()[i].gf_M_m());
#endif
	  }
	}

	
	if (include_postprocessing && group_name == "Postprocessing_SDE_Integrands"){
	    this->append_observable_complex_gf_container(group_name + "/Sig_sc", postproc_state.gf_Sig_SDE_integrand_sc());
	    this->append_observable_complex_gf_container(group_name + "/Sig_d", postproc_state.gf_Sig_SDE_integrand_d());
	    this->append_observable_complex_gf_container(group_name + "/Sig_m", postproc_state.gf_Sig_SDE_integrand_m());

	    auto sig_info = this->Eval_max_info(postproc_state.gf_Sig_SDE_integrand_sc(), postproc_state.gf_Sig_SDE_integrand_d(), postproc_state.gf_Sig_SDE_integrand_m());

	    this->append_observable_string(group_name+"/Leading_name", std::get<0>(sig_info));
	    this->append_observable_complex(group_name+"/Max_sc", std::get<1>(sig_info));
	    this->append_observable_complex(group_name+"/Max_d", std::get<2>(sig_info));
	    this->append_observable_complex(group_name+"/Max_m", std::get<3>(sig_info));

            this->append_observable_multiindex(group_name+"/Max_sc_idx", std::get<4>(sig_info));
	    this->append_observable_multiindex(group_name+"/Max_d_idx", std::get<5>(sig_info));
	    this->append_observable_multiindex(group_name+"/Max_m_idx", std::get<6>(sig_info));

	}

	if (include_postprocessing && group_name == "Postprocessing_Susc_info"){

	    this->append_observable_complex_gf_container(group_name + "/Susc_sc", postproc_state.gf_susc_sc());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d", postproc_state.gf_susc_d());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m", postproc_state.gf_susc_m());

	    auto susc_info = this->Eval_max_info(postproc_state.gf_susc_sc(), postproc_state.gf_susc_d(), postproc_state.gf_susc_m());

	    this->append_observable_string(group_name+"/Leading_name", std::get<0>(susc_info));
	    this->append_observable_complex(group_name+"/Max_sc", std::get<1>(susc_info));
	    this->append_observable_complex(group_name+"/Max_d", std::get<2>(susc_info));
	    this->append_observable_complex(group_name+"/Max_m", std::get<3>(susc_info));

            this->append_observable_multiindex(group_name+"/Max_sc_idx", std::get<4>(susc_info));
	    this->append_observable_multiindex(group_name+"/Max_d_idx", std::get<5>(susc_info));
	    this->append_observable_multiindex(group_name+"/Max_m_idx", std::get<6>(susc_info));

#ifdef SPLIT_SUSC_CONTRIBUTIONS
	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution", postproc_state.gf_suscvert_sc());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution", postproc_state.gf_suscvert_d());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution", postproc_state.gf_suscvert_m());
	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_bubble_contribution", postproc_state.gf_suscbubble_sc());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_bubble_contribution", postproc_state.gf_suscbubble_d());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_bubble_contribution", postproc_state.gf_suscbubble_m());

	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_nabla_sc", postproc_state.gf_suscvert_sc_contribution_from_nabla_sc());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_nabla_sc", postproc_state.gf_suscvert_d_contribution_from_nabla_sc());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_nabla_sc", postproc_state.gf_suscvert_m_contribution_from_nabla_sc());	    
	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_M_sc", postproc_state.gf_suscvert_sc_contribution_from_M_sc());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_M_sc", postproc_state.gf_suscvert_d_contribution_from_M_sc());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_M_sc", postproc_state.gf_suscvert_m_contribution_from_M_sc());

	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_sc_double_counting", postproc_state.gf_suscvert_sc_contribution_from_sc_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_sc_double_counting", postproc_state.gf_suscvert_d_contribution_from_sc_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_sc_double_counting", postproc_state.gf_suscvert_m_contribution_from_sc_double_counting());

	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_nabla_d", postproc_state.gf_suscvert_sc_contribution_from_nabla_d());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_nabla_d", postproc_state.gf_suscvert_d_contribution_from_nabla_d());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_nabla_d", postproc_state.gf_suscvert_m_contribution_from_nabla_d());	    
	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_M_d", postproc_state.gf_suscvert_sc_contribution_from_M_d());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_M_d", postproc_state.gf_suscvert_d_contribution_from_M_d());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_M_d", postproc_state.gf_suscvert_m_contribution_from_M_d());

	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_d_double_counting", postproc_state.gf_suscvert_sc_contribution_from_d_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_d_double_counting", postproc_state.gf_suscvert_d_contribution_from_d_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_d_double_counting", postproc_state.gf_suscvert_m_contribution_from_d_double_counting());

	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_nabla_m", postproc_state.gf_suscvert_sc_contribution_from_nabla_m());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_nabla_m", postproc_state.gf_suscvert_d_contribution_from_nabla_m());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_nabla_m", postproc_state.gf_suscvert_m_contribution_from_nabla_m());	    
	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_M_m", postproc_state.gf_suscvert_sc_contribution_from_M_m());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_M_m", postproc_state.gf_suscvert_d_contribution_from_M_m());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_M_m", postproc_state.gf_suscvert_m_contribution_from_M_m());

	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_m_double_counting", postproc_state.gf_suscvert_sc_contribution_from_m_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_m_double_counting", postproc_state.gf_suscvert_d_contribution_from_m_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_m_double_counting", postproc_state.gf_suscvert_m_contribution_from_m_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Susc_sc_vertex_contribution_from_bare_vertex", postproc_state.gf_suscvert_sc_contribution_from_bare());
	    this->append_observable_complex_gf_container(group_name + "/Susc_d_vertex_contribution_from_bare_vertex", postproc_state.gf_suscvert_d_contribution_from_bare());
	    this->append_observable_complex_gf_container(group_name + "/Susc_m_vertex_contribution_from_bare_vertex", postproc_state.gf_suscvert_m_contribution_from_bare());
#endif

	}

	if (include_postprocessing && group_name == "Postprocessing_Polarisation_info"){
	    gf_polarisation_t<Model> polarisation_sc = postproc_state.gf_polarisation_sc_contribution_from_1() + postproc_state.gf_polarisation_sc_contribution_from_varphi();
	    gf_polarisation_t<Model> polarisation_d = postproc_state.gf_polarisation_d_contribution_from_1() + postproc_state.gf_polarisation_d_contribution_from_varphi();
	    gf_polarisation_t<Model> polarisation_m = postproc_state.gf_polarisation_m_contribution_from_1() + postproc_state.gf_polarisation_m_contribution_from_varphi();

	    this->append_observable_complex_gf_container(group_name + "/Polarisation_sc", polarisation_sc);
	    this->append_observable_complex_gf_container(group_name + "/Polarisation_d", polarisation_d);
	    this->append_observable_complex_gf_container(group_name + "/Polarisation_m", polarisation_m);

	    auto polarisation_info = this->Eval_max_info(polarisation_sc, polarisation_d, polarisation_m);

	    this->append_observable_string(group_name+"/Leading_name", std::get<0>(polarisation_info));
	    this->append_observable_complex(group_name+"/Max_sc", std::get<1>(polarisation_info));
	    this->append_observable_complex(group_name+"/Max_d", std::get<2>(polarisation_info));
	    this->append_observable_complex(group_name+"/Max_m", std::get<3>(polarisation_info));

            this->append_observable_multiindex(group_name+"/Max_sc_idx", std::get<4>(polarisation_info));
	    this->append_observable_multiindex(group_name+"/Max_d_idx", std::get<5>(polarisation_info));
	    this->append_observable_multiindex(group_name+"/Max_m_idx", std::get<6>(polarisation_info));



	    this->append_observable_complex_gf_container(group_name + "/Polarisation_sc_minus_bubble", postproc_state.gf_polarisation_sc_contribution_from_varphi());
	    this->append_observable_complex_gf_container(group_name + "/Polarisation_d_minus_bubble", postproc_state.gf_polarisation_d_contribution_from_varphi());
	    this->append_observable_complex_gf_container(group_name + "/Polarisation_m_minus_bubble", postproc_state.gf_polarisation_m_contribution_from_varphi());

	    this->append_observable_complex_gf_container(group_name + "/Polarisation_sc_bubble", postproc_state.gf_polarisation_sc_contribution_from_1());
	    this->append_observable_complex_gf_container(group_name + "/Polarisation_d_bubble", postproc_state.gf_polarisation_d_contribution_from_1());
	    this->append_observable_complex_gf_container(group_name + "/Polarisation_m_bubble", postproc_state.gf_polarisation_m_contribution_from_1());

	    this->append_observable_complex_gf_container(group_name + "/Polarisation_sc_contribution_from_nabla_sc", postproc_state.gf_polarisation_sc_contribution_from_nabla_sc());
	    this->append_observable_complex_gf_container(group_name + "/Polarisation_d_contribution_from_nabla_d", postproc_state.gf_polarisation_d_contribution_from_nabla_d());
	    this->append_observable_complex_gf_container(group_name + "/Polarisation_m_contribution_from_nabla_m", postproc_state.gf_polarisation_m_contribution_from_nabla_m());
	}

	if (include_postprocessing && group_name == "Postprocessing_Lambda_info"){

	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_minus_1", postproc_state.gf_lambda_sc_contribution_from_varphi());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_minus_1", postproc_state.gf_lambda_d_contribution_from_varphi());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_minus_1", postproc_state.gf_lambda_m_contribution_from_varphi());

	    auto lambda_info = this->Eval_max_info(postproc_state.gf_lambda_sc_contribution_from_varphi(), postproc_state.gf_lambda_d_contribution_from_varphi(), postproc_state.gf_lambda_m_contribution_from_varphi());

	    this->append_observable_string(group_name+"/Leading_name", std::get<0>(lambda_info));
	    this->append_observable_complex(group_name+"/Max_sc", std::get<1>(lambda_info));
	    this->append_observable_complex(group_name+"/Max_d", std::get<2>(lambda_info));
	    this->append_observable_complex(group_name+"/Max_m", std::get<3>(lambda_info));

            this->append_observable_multiindex(group_name+"/Max_sc_idx", std::get<4>(lambda_info));
	    this->append_observable_multiindex(group_name+"/Max_d_idx", std::get<5>(lambda_info));
	    this->append_observable_multiindex(group_name+"/Max_m_idx", std::get<6>(lambda_info));

#ifdef SPLIT_LAMBDA_CONTRIBUTIONS
	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_nabla_d", postproc_state.gf_lambda_sc_contribution_from_nabla_d());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_nabla_m", postproc_state.gf_lambda_sc_contribution_from_nabla_m());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_nabla_sc", postproc_state.gf_lambda_d_contribution_from_nabla_sc());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_nabla_m", postproc_state.gf_lambda_d_contribution_from_nabla_m());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_nabla_sc", postproc_state.gf_lambda_m_contribution_from_nabla_sc());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_nabla_d", postproc_state.gf_lambda_m_contribution_from_nabla_d());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_nabla_m", postproc_state.gf_lambda_m_contribution_from_nabla_m());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_nabla_d", postproc_state.gf_lambda_d_contribution_from_nabla_d());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_nabla_sc", postproc_state.gf_lambda_sc_contribution_from_nabla_sc());

	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_M_d", postproc_state.gf_lambda_sc_contribution_from_M_d());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_M_m", postproc_state.gf_lambda_sc_contribution_from_M_m());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_M_sc", postproc_state.gf_lambda_d_contribution_from_M_sc());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_M_m", postproc_state.gf_lambda_d_contribution_from_M_m());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_M_sc", postproc_state.gf_lambda_m_contribution_from_M_sc());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_M_d", postproc_state.gf_lambda_m_contribution_from_M_d());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_M_m", postproc_state.gf_lambda_m_contribution_from_M_m());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_M_d", postproc_state.gf_lambda_d_contribution_from_M_d());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_M_sc", postproc_state.gf_lambda_sc_contribution_from_M_sc());
	    
	    
	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_sc_double_counting", postproc_state.gf_lambda_sc_contribution_from_sc_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_d_double_counting", postproc_state.gf_lambda_sc_contribution_from_d_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_m_double_counting", postproc_state.gf_lambda_sc_contribution_from_m_double_counting());

	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_sc_double_counting", postproc_state.gf_lambda_d_contribution_from_sc_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_d_double_counting", postproc_state.gf_lambda_d_contribution_from_d_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_m_double_counting", postproc_state.gf_lambda_d_contribution_from_m_double_counting());

	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_sc_double_counting", postproc_state.gf_lambda_m_contribution_from_sc_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_d_double_counting", postproc_state.gf_lambda_m_contribution_from_d_double_counting());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_m_double_counting", postproc_state.gf_lambda_m_contribution_from_m_double_counting());

	    
	    this->append_observable_complex_gf_container(group_name + "/Lambda_m_contribution_from_bare_vertex", postproc_state.gf_lambda_m_contribution_from_bare());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_sc_contribution_from_bare_vertex", postproc_state.gf_lambda_sc_contribution_from_bare());
	    this->append_observable_complex_gf_container(group_name + "/Lambda_d_contribution_from_bare_vertex", postproc_state.gf_lambda_d_contribution_from_bare());
#endif

	}
    }

    

    for (auto &observable : observables_frg_common_t<Model, state_t >::ObservablesListToTrack()){
	std::string observable_name = std::get<0>(observable);
	
	if (observable_name == "Free_Energy"){
	  //std::complex<double> value = state.eval_f(t);
          std::complex<double> value = state.gf_f()[0];
	  this->append_observable_complex(observable_name, value);
	}
    }
}
