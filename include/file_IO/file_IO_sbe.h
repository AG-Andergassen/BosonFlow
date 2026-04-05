#pragma once 

#include <file_IO/file_IO_base.h>


template <typename Model, typename State>
class FileIOSBE : public FileIOBase<Model, State>
{
 public:
    FileIOSBE() = delete; ///< This is not a "type", so no objects of this type may be constructed. Rather, one should use the static methods

    static void write_final_file(const State &state_vec, const std::unordered_map<std::string, void* > &observables_name_data_map,  const std::string to_be_appended_to_file_name = ""); ///< Writes state data and observables into the .h5 file at the end of the run
    static void write_checkpoint_file( const State& state_vec, const std::unordered_map<std::string, void* > &observables_name_data_map, size_t count, double start_time, double end_time, double dt, size_t AorB );  ///<    Creates a checkpoint file with state data and observables during the run. UPDATER_FREQ and UPDATER_START control at which integration step = count this happens for the first time and at which intervals after that. AorB is used to alternatingly overwrite files to avoid using too much memory.

    static void read_state_from_checkpoint_file(const std::string input_file_name, State &state);

    static gf_1p_mat_t<Model> read_initial_Ginit_from_DMFT_data(const std::string input_file_name);

    static void read_initial_state_from_DMFT_data(const std::string input_file_name, State &state);

 protected:
    static void write_common_data(H5::H5File &file, const State &state, const std::unordered_map<std::string, void* > &observables_name_data_map);
    
    static void write_w_func( H5::H5File& file, const State& state_vec ); 	
    static void write_lambda_func( H5::H5File& file, const State& state_vec ); 	
    static void write_M_func( H5::H5File& file, const State& state_vec ); 	
};


#include "../src/file_IO/file_IO_sbe.tpp"
//#include "../src/file_IO/file_read_sbe_dmft_data.tpp"
