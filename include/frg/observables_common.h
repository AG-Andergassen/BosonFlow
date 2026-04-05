#pragma once

#include <base_gf_types.h>
#include <models/concrete_available_models.h>
#include <string>

enum class FlowObservableType{ GF_MULTIINDEX, INT, DOUBLE, COMPLEX, COMPLEX_GF_CONTAINER, STRING };


// a global container that contains all observables to be tracked and outputted
template <typename Model, typename state_t>
class observables_frg_common_t
{
 public:    
    observables_frg_common_t();

    ~observables_frg_common_t();

    virtual void update(const state_t &state, double t, bool is_final);    

    void append_observable_string(std::string name, std::string value);
    void append_observable_multiindex(std::string name, std::vector<int> value);
    void append_observable_int(std::string name, int value);
    void append_observable_double(std::string name, double value);
    void append_observable_complex(std::string name, std::complex<double> value);    
    template<unsigned N>	
    void append_observable_complex_gf_container(std::string name, const gf<std::complex<double>, N> &value); 



    std::unordered_map<std::string, void* > m_name_data_map;
    
    static void SetObservablesToTrack();

    static void SetToTrackAllAvailableObservables();

    static std::vector< std::tuple<std::string, FlowObservableType, unsigned> > &ObservablesListToTrack()
    {
	static std::vector< std::tuple<std::string, FlowObservableType, unsigned> > observables_list_to_track;
	return observables_list_to_track;
    }

    static std::vector< std::tuple<std::string, std::vector<std::tuple<std::string, FlowObservableType, unsigned> > > > &GroupObservablesListToTrack()
    {
	static std::vector< std::tuple<std::string, std::vector<std::tuple<std::string, FlowObservableType, unsigned> > > > group_observables_list_to_track;
	return group_observables_list_to_track;
    }


    // functions that evaluate the observables
    static double Eval_filling(const state_t& state, double t);

    template< typename gf_type >
    static std::tuple<std::string, std::complex<double>, std::complex<double>, std::complex<double>, std::vector<int>, std::vector<int>, std::vector<int> > Eval_max_info(const gf_type &gf_sc, const gf_type &gf_d, const gf_type &gf_m);

    static std::pair<std::complex<double>, int > Eval_max_coupling( const std::complex<double> *gf_data, unsigned gf_num_elements );


};

// include template implementation
#include "../src/frg/observables_common.tpp"
