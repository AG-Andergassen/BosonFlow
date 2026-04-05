
/*******************************************************************************************//** @file
 *  		
 * 	file: 		multifile_IO_sbe.h
 * 	contents:  	Functions to write output to hdf5 file
 * 
 ****************************************************************************************************/


#pragma once

#include <file_IO/H5Compat.h>
#include <frg/sbe/state.h>
#include <frg/observables_common.h>
#include <file_IO/file_IO_sbe.h>

struct resume_details_t
{
    double start_time = 0.0;
    double final_time = 1.0;
    double init_timestep = 0.05;
    double is_finished = false;
};

template <typename Model, typename State>
class MultiFileIOSBE : public FileIOSBE<Model, State>
{
 public:
    MultiFileIOSBE() = delete; ///< This is not a "type", so no objects of this type may be constructed. Rather, one should use the static methods

    static void Init(resume_details_t &resume_details, std::string to_be_appended_to_out_dir_name = ""); // creates output directory

    static void WritePropertiesFile();

    static void IncrementScaleIdx();

    static void WriteCurrentScaleFile(const State &state, const std::unordered_map<std::string, void* > &observables_name_data_map);

    static void WriteFinalFile(const State &state_ptr, const std::unordered_map<std::string, void* > &observables_name_data_map, const std::string to_be_appended_to_file_name = "");


    static int &CurrentScaleIdx()
    {
	static int current_scale_idx = 0;
	return current_scale_idx;
    }

    static double &StartRunTimeSeconds()
    {
	static double start_run_time_seconds = 0;
	return start_run_time_seconds;
    }

    static std::string &OutputDirectory()
    {
	static std::string output_directory = "";
	return output_directory;
    }

};

// include implementations
#include "../src/file_IO/multifile_IO_sbe.tpp"
