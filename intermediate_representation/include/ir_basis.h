#pragma once

#include <julia.h>
#include <iostream>
#include <vector>
#include <complex>
#include <string> 
#include <tuple>

JULIA_DEFINE_FAST_TLS


struct IR_2p_basis_data_t{
    std::vector<int> sampling_freqs_full;

    std::vector<std::vector<std::complex<double> > > E; // interpolation matrix through the IR basis from the sampling frequencies to the entire frequency set
};

struct IR_3p_basis_data_t{
    std::vector<std::tuple<int, int, int> > sampling_freqs_full;

    std::vector<std::vector<std::complex<double> > > E_pp; // interpolation matrix from the IR basis in the pp parametrisation
    std::vector<std::vector<std::complex<double> > > E_ph; // interpolation matrix from the IR basis in the pp parametrisation
    std::vector<std::vector<std::complex<double> > > E_xph; // interpolation matrix from the IR basis in the xph parametrisation
};

struct IR_4p_basis_data_t{
    std::vector<std::tuple<int, int, int, int> > sampling_freqs_full;    

    std::vector<std::vector<std::complex<double> > > E_pp; // interpolation matrix from the IR basis in the pp parametrisation
    std::vector<std::vector<std::complex<double> > > E_ph; // interpolation matrix from the IR basis in the pp parametrisation
    std::vector<std::vector<std::complex<double> > > E_xph; // interpolation matrix from the IR basis in the xph parametrisation

};


class IR_basis_factory_t
{
    IR_basis_factory_t(double beta, double error);

    // e.g. for the self-energy or green's functions
    IR_2p_basis_data_t make_2p_basis(int total_number_of_freqs, double bandwidth, Statistics zeta, bool augment_with_constant = false);
    
    IR_2p_basis_data_t make_2p_basis_with_max_sampling_frequency(int max_sampling_frequency, Statistics zeta, bool augment_with_constant = false);

    // e.g. for hedin vertex or the truncated unity bubble
    IR_3p_basis_data_t make_3p_basis(int total_number_of_freqs_bos, double bandwidth_bos, int total_number_of_freqs_ferm, double bandwidth_ferm, bool augment_with_constant = false);

    IR_3p_basis_data_t make_3p_basis_with_max_sampling_frequency(int max_sampling_frequency_bos, int max_sampling_frequency_ferm, bool augment_with_constant = false);


    // e.g. for the rest function
    IR_4p_basis_data_t make_4p_basis(int total_number_of_freqs_bos, double bandwidth_bos, int total_number_of_freqs_ferm, double bandwidth_ferm, bool augment_with_constant = false);

    IR_4p_basis_data_t make_4p_basis_with_max_sampling_frequency(int max_sampling_frequency_bos, int max_sampling_frequency_ferm, bool augment_with_constant = false);


    double obtain_bandwidth_from_max_sampling_frequency(int max_sampling_frequency, Statistics zeta);

 private:
    void initialise_julia_context_and_headers();

    void destroy_julia_context();

    template<typename T> 
    void get_julia_array(const char* julia_array_name, std::vector<T> &array)
    {
	jl_value_t *julia_var = jl_eval_string(julia_array_name);
	jl_array_t *julia_arr = (jl_array_t*)(julia_var);
	T *data = (T*)jl_array_data(julia_arr);

	size_t size0 = jl_array_dim(julia_arr, 0);

	array.resize(size0);
	for (int i = 0; i < size0; i ++)
	    array[i] = data[i];
    }

    template<typename T>
    void get_julia_2Darray(const char* julia_array_name, std::vector<std::vector<T> > &array)
    {
	jl_value_t *julia_var = jl_eval_string(julia_array_name);
	jl_array_t *julia_arr = (jl_array_t*)(julia_var);
	T *data = (T*)jl_array_data(julia_arr);

	size_t size0 = jl_array_dim(julia_arr, 0);
	size_t size1 = jl_array_dim(julia_arr, 1);

	array.resize(size0, std::vector<T>(size1));

	for (int j = 0; j < size1; j ++)
	    for (int i = 0; i < size0; i ++)
		array[i][j] = data[i + size0*j];
    
    }

    void get_julia_complex2Darray(const char* julia_array_name, std::vector<std::vector<std::complex<double> > > &array)
    {
	std::vector<std::vector<double> > real_part;
	std::vector<std::vector<double> > imag_part;
	get_julia_2Darray(("real.(" + std::string(julia_array_name) + ")").c_str(), real_part);
	get_julia_2Darray(("imag.(" + std::string(julia_array_name) + ")").c_str(), imag_part);

	array.resize(real_part.size(), std::vector<std::complex<double> >(real_part[0].size()));

	for (int i = 0; i < real_part.size(); i ++)
	    for (int j = 0; j < real_part[0].size(); j ++)
		array[i][j] = std::complex<double>(real_part[i][j], imag_part[i][j]);

    }

    template<typename T>
    void get_julia_tuple4array(const char* julia_array_name, std::vector<std::tuple<T, T, T, T>  > &array)
    {
	std::vector<T> firsts, seconds, thirds, fourths;

	get_julia_array(("map(x -> x[1], (" + std::string(julia_array_name) + ")").c_str(), firsts);
	get_julia_array(("map(x -> x[2], (" + std::string(julia_array_name) + ")").c_str(), seconds);
	get_julia_array(("map(x -> x[3], (" + std::string(julia_array_name) + ")").c_str(), thirds);
	get_julia_array(("map(x -> x[4], (" + std::string(julia_array_name) + ")").c_str(), fourths);


	array.resize(firsts.size(), std::tuple<T, T, T, T>);

	for (int i = 0; i < firsts.size(); i ++)
	    array.push_back( std::make_tuple(firsts[i], seconds[i], thirds[i], fourths[i]) );

    }

    template<typename T>
    void get_julia_tuple3array(const char* julia_array_name, std::vector<std::tuple<T, T, T, T>  > &array)
    {
	std::vector<T> firsts, seconds, thirds;

	get_julia_array(("map(x -> x[1], (" + std::string(julia_array_name) + ")").c_str(), firsts);
	get_julia_array(("map(x -> x[2], (" + std::string(julia_array_name) + ")").c_str(), seconds);
	get_julia_array(("map(x -> x[3], (" + std::string(julia_array_name) + ")").c_str(), thirds);

	array.resize(firsts.size(), std::tuple<T, T, T>);

	for (int i = 0; i < firsts.size(); i ++)
	    array.push_back( std::make_tuple(firsts[i], seconds[i], thirds[i]) );

    }


};
