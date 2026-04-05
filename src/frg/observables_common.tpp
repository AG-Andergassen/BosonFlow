#include <frg/common_frg_gf_types.h>
#include <frg/flows.h>
#include <frequencies/matsubara_space.h>
#include <frg/rhs_base.h>

template <typename Model, typename state_t>
observables_frg_common_t<Model, state_t>::observables_frg_common_t()
{
    auto observables_list = ObservablesListToTrack();
    
    for (auto subgroup : GroupObservablesListToTrack())
    {
	auto subgroup_name = std::get<0>(subgroup);
	auto observables_in_subgroup = std::get<1>(subgroup);
	for (auto o : observables_in_subgroup){
	    observables_list.push_back({subgroup_name + "/" + std::get<0>(o), std::get<1>(o), std::get<2>(o)});
	}
	
    }

    for (auto &observable: observables_list ){
	std::string observable_name = std::get<0>(observable);
	FlowObservableType observable_type = std::get<1>(observable);
	unsigned gf_ndims = std::get<2>(observable);

	switch (observable_type){
	
	case FlowObservableType::INT:
	    m_name_data_map[observable_name] = new std::vector<int>();
	    break;

	case FlowObservableType::DOUBLE:
	    m_name_data_map[observable_name] = new std::vector<double>();
	    break;
	
	case FlowObservableType::COMPLEX:
	    m_name_data_map[observable_name] = new std::vector<std::complex<double> >();
	    break;
	   
	case FlowObservableType::GF_MULTIINDEX:
	    m_name_data_map[observable_name] = new std::vector< std::vector<int> >();
	    break;


	case FlowObservableType::STRING:
	    m_name_data_map[observable_name] = new std::vector<std::string>();
	    break;

	case FlowObservableType::COMPLEX_GF_CONTAINER:

#define MakeGFNameDataMap(X) \
	    if (gf_ndims == X) \
		m_name_data_map[observable_name] = new std::vector<gf<std::complex<double>, X> >()

	    // to be able to include gf observables of ndim up to 21 
	    MakeGFNameDataMap(1); MakeGFNameDataMap(2); MakeGFNameDataMap(3);
	    MakeGFNameDataMap(4); MakeGFNameDataMap(5); MakeGFNameDataMap(6);
	    MakeGFNameDataMap(7); MakeGFNameDataMap(8); MakeGFNameDataMap(9);
	    MakeGFNameDataMap(10); MakeGFNameDataMap(11); MakeGFNameDataMap(12);
	    MakeGFNameDataMap(13); MakeGFNameDataMap(14); MakeGFNameDataMap(15);
	    MakeGFNameDataMap(16); MakeGFNameDataMap(17); MakeGFNameDataMap(18);
	    MakeGFNameDataMap(19); MakeGFNameDataMap(20); MakeGFNameDataMap(21);
	    break;

	default:;
	}
	
    }
}


template <typename Model, typename state_t>
observables_frg_common_t<Model, state_t>::~observables_frg_common_t()
{
    for (auto &observable: ObservablesListToTrack() ){
	std::string observable_name = std::get<0>(observable);
	FlowObservableType observable_type = std::get<1>(observable);
	unsigned gf_ndims = std::get<2>(observable);
       
	switch (observable_type){
	    
	case FlowObservableType::INT:
	    delete static_cast<std::vector<int>*>(m_name_data_map[observable_name]);
	    break;
	case FlowObservableType::DOUBLE:
	    delete static_cast<std::vector<double>*>(m_name_data_map[observable_name]);
	    break;
	
	case FlowObservableType::COMPLEX:
	    delete static_cast<std::vector<std::complex<double> >*>(m_name_data_map[observable_name]);
	    break;

	case FlowObservableType::STRING:
	    delete static_cast<std::vector<std::string >*>(m_name_data_map[observable_name]);
	    break;

	case FlowObservableType::GF_MULTIINDEX:
	    delete static_cast<std::vector<std::vector<int> >*>(m_name_data_map[observable_name]);
	    break;

	    
	case FlowObservableType::COMPLEX_GF_CONTAINER:

#define DeleteGFData(X) \
	    if (gf_ndims == X) \
		delete static_cast<std::vector<gf<std::complex<double>, X> >*>(m_name_data_map[observable_name]);

	    // to be able to include gf observables of ndim up to 21 
	    DeleteGFData(1); DeleteGFData(2); DeleteGFData(3);
	    DeleteGFData(4); DeleteGFData(5); DeleteGFData(6);
	    DeleteGFData(7); DeleteGFData(8); DeleteGFData(9);
	    DeleteGFData(10); DeleteGFData(11); DeleteGFData(12);
	    DeleteGFData(13); DeleteGFData(14); DeleteGFData(15);
	    DeleteGFData(16); DeleteGFData(17); DeleteGFData(18);
	    DeleteGFData(19); DeleteGFData(20); DeleteGFData(21);
	    break;


	default:;
	}
    }
}

template <typename Model, typename state_t>
void observables_frg_common_t<Model, state_t>::append_observable_string(std::string name, std::string value)
{
#ifdef MULTIFILE_OUTPUT
    static_cast<std::vector<std::string>*>(m_name_data_map[name])->clear();
#endif
    static_cast<std::vector<std::string>*>(m_name_data_map[name])->push_back(value);
}


template <typename Model, typename state_t>
void observables_frg_common_t<Model, state_t>::append_observable_double(std::string name, double value)
{
#ifdef MULTIFILE_OUTPUT
    static_cast<std::vector<double>*>(m_name_data_map[name])->clear();
#endif
    static_cast<std::vector<double>*>(m_name_data_map[name])->push_back(value);
}


template <typename Model, typename state_t>
void observables_frg_common_t<Model, state_t>::append_observable_complex(std::string name, std::complex<double> value)
{
#ifdef MULTIFILE_OUTPUT
    static_cast<std::vector<std::complex<double> >*>(m_name_data_map[name])->clear();
#endif
    static_cast<std::vector<std::complex<double> >*>(m_name_data_map[name])->push_back(value); 
}
    

template <typename Model, typename state_t>
void observables_frg_common_t<Model, state_t>::append_observable_multiindex(std::string name, std::vector<int> value)
{
#ifdef MULTIFILE_OUTPUT
    static_cast<std::vector<std::vector<int> >*>(m_name_data_map[name])->clear();
#endif
    static_cast<std::vector<std::vector<int> >*>(m_name_data_map[name])->push_back(value); 
}


template <typename Model, typename state_t> 
template <unsigned N>
void observables_frg_common_t<Model, state_t>::append_observable_complex_gf_container(std::string name, const gf<std::complex<double>, N> &value)
{
#ifdef MULTIFILE_OUTPUT
    static_cast<std::vector< gf<std::complex<double>, N> >*>(m_name_data_map[name])->clear();
#endif
    static_cast<std::vector< gf<std::complex<double>, N> >*>(m_name_data_map[name])->push_back(value);
}



template <typename Model, typename state_t>
void observables_frg_common_t<Model, state_t>::update(const state_t &state, double t, bool is_final)
{
    for (auto &observable : ObservablesListToTrack()){
	std::string observable_name = std::get<0>(observable); 
	
	if (observable_name == "Filling"){
	    append_observable_double(observable_name, Eval_filling(state, t));
	}else if (observable_name == "t"){
	    append_observable_double(observable_name, t);
	}else if (observable_name == "Sig"){
	    append_observable_complex_gf_container<gf_1p_t<Model>::base_t::ndims>(observable_name, state.gf_Sig());
	}else if (observable_name == "delta_mu"){
	    append_observable_double(observable_name, state.gf_delta_mu());
	}else if (observable_name == "delta_mu_dot"){
	    append_observable_double(observable_name, state.m_d_delta_mu_over_dt);
	}else if (observable_name == "Beta"){
	    append_observable_double(observable_name, BETA);
	}else{
	    // ...
	}
    }
}

template <typename Model, typename state_t>
void observables_frg_common_t<Model, state_t>::SetToTrackAllAvailableObservables()
{
    ObservablesListToTrack().push_back({"Filling", FlowObservableType::DOUBLE, 0});

    // if self-energy flow
    ObservablesListToTrack().push_back({"Sig", FlowObservableType::COMPLEX_GF_CONTAINER, gf_1p_t<Model>::base_t::ndims});


    // scale parameters
    ObservablesListToTrack().push_back({"t", FlowObservableType::DOUBLE, 0});

    // chemical potential shift
    ObservablesListToTrack().push_back({"delta_mu", FlowObservableType::DOUBLE, 0});
    ObservablesListToTrack().push_back({"delta_mu_dot", FlowObservableType::DOUBLE, 0});

    ObservablesListToTrack().push_back({"Beta", FlowObservableType::DOUBLE, 0});


    /// observables on special points

    //ObservablesListToTrack().push_back({"Sig_at_w0_", FlowObservableType::COMPLEX, 0});

    // if not static calculation
    //ObservablesListToTrack().push_back({"Sig_at_w1_", FlowObservableType::COMPLEX, 0});

    // observables on special paths
    //...

}




template <typename Model, typename state_t>
double observables_frg_common_t<Model, state_t>::Eval_filling( const state_t& state, double t )
{   
    std::cout << state.gf_delta_mu()<<std::endl;
    double filling = state.eval_filling(t, state.gf_delta_mu());
    std::cout << "Filling: " << filling << std::endl;
    return filling; 
}


template <typename Model, typename state_t>
template <typename gf_type>
std::tuple<std::string, std::complex<double>, std::complex<double>, std::complex<double>, std::vector<int>, std::vector<int>, std::vector<int> > observables_frg_common_t<Model, state_t>::Eval_max_info(const gf_type &gf_sc, const gf_type &gf_d, const gf_type &gf_m)
{
    auto [max_sc, flat_max_idx_sc] = Eval_max_coupling(gf_sc.data(), gf_sc.num_elements());
    auto max_idx_sc = gf_sc.get_idx(flat_max_idx_sc);

    auto [max_d, flat_max_idx_d] = Eval_max_coupling(gf_d.data(), gf_d.num_elements());
    auto max_idx_d = gf_d.get_idx(flat_max_idx_d);

    auto [max_m, flat_max_idx_m] = Eval_max_coupling(gf_m.data(), gf_m.num_elements());
    auto max_idx_m = gf_m.get_idx(flat_max_idx_m);

    std::string leading_str; 
    if (std::abs(std::abs(max_sc) - std::abs(max_m)) < 1e-3 && std::abs(std::abs(max_sc) - std::abs(max_d)) < 1e-3 ){
	leading_str = "";
    }else if (std::abs(std::abs(max_sc) - std::abs(max_m)) < 1e-3 && std::abs(max_sc) > std::abs(max_d)){
	leading_str = "sc,m";
    }else if (std::abs(std::abs(max_sc) - std::abs(max_d)) < 1e-3 && std::abs(max_sc) > std::abs(max_m)){
	leading_str = "sc,d";
    }else if (std::abs(std::abs(max_d) - std::abs(max_m)) < 1e-3  && std::abs(max_d) > std::abs(max_sc)){
	leading_str = "d,m";
    }else if (std::abs(max_sc) - std::abs(max_m) > 1e-3  && std::abs(max_sc) - std::abs(max_d) > 1e-3){
	leading_str = "sc";
    }else if (std::abs(max_m) - std::abs(max_d) > 1e-3){
	leading_str = "m";
    }else{
	leading_str = "d";
    }


    std::vector<int> max_multiindex_sc(max_idx_sc.idx_arr.begin(), max_idx_sc.idx_arr.end());
    std::vector<int> max_multiindex_d(max_idx_d.idx_arr.begin(), max_idx_d.idx_arr.end());
    std::vector<int> max_multiindex_m(max_idx_m.idx_arr.begin(), max_idx_m.idx_arr.end());

    return {leading_str, max_sc, max_d, max_m, max_multiindex_sc, max_multiindex_d, max_multiindex_m};
}

template <typename Model, typename state_t>
std::pair<std::complex<double>, int > observables_frg_common_t<Model, state_t>::Eval_max_coupling( const std::complex<double> *gf_data, unsigned gf_num_elements )
{
    auto comparator = []( const std::complex<double>& a, const std::complex<double>& b )->bool{ return std::abs(a) < std::abs(b);};

    const std::complex<double>* the_max_element_ptr = std::max_element( gf_data, gf_data + gf_num_elements, comparator);
    int flat_index_of_max_element = std::distance(gf_data, the_max_element_ptr);

    return {*the_max_element_ptr, flat_index_of_max_element};
}
