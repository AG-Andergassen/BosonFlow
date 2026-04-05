#pragma once

#include <file_IO/H5Compat.h>

#include <params_physical.h>
#include <params_technical.h> 
#include <ctime>
#include <file_IO/H5Tools.h>
#include <cmath>
#include <frg/flows.h>
#include <tu_projections.h>
#include <frg/observables_common.h>


template <typename Model, typename State>
class FileIOBase
{
 public:
    FileIOBase() = delete; ///< This is not a "type", so no objects of this type may be constructed. Rather, one should use the static methods

 protected:

    static H5std_string construct_name_from_config(); ///< Builds file name from the preset settings of the calculation

    static void write_config( H5::H5File& file );		///< Write configuration to output file
    static void write_params( H5::H5File& file );		///< Write parameters to output file

    static void write_CheckpointInfo( H5::H5File& file, const double count, const double start_time, const double end_time, const double dt);      //<     Write number of current integration step, current t, and dt

    static void write_flow_observables( H5::H5File& file, const std::unordered_map<std::string, void* > &observables_name_data_map );	///<	Write observables tracked during the flow
    static void write_flow_observables_in_subgroup( H5::Group &group, const std::vector<std::tuple< std::string, FlowObservableType, unsigned> > &observables_info, const std::unordered_map<std::string, void* > &observables_name_data_map, std::string group_directory);
    
    
      static void write_bubble( H5::H5File& file, const gf_bubble_mat_t<Model> bubble_pp, const gf_bubble_mat_t<Model> bubble_ph ); //< Write bubble data (pp and ph)
      static void write_projection( H5::H5File& file ); //< Write projection data
      static void write_Gfunc( H5::H5File& file, const gf_1p_mat_real_t<Model> Gvec_real, const gf_1p_mat_real_t<Model> Svec_real ); //< Write Green's functions (G and S) in real space

    static void write_Sig( H5::H5File& file, const State& state_vec );		///<	Write self-energy to output file
};

// include implementations
#include "../src/file_IO/file_IO_base.tpp"
