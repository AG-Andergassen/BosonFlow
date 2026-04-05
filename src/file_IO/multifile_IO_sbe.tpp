#include <file_IO/H5Tools.h>

#include <params_physical.h>
#include <params_technical.h> 
#include <ctime>
#include <cmath>
#include <models/concrete_available_models.h>
#include <frg/flows.h>
#include <tu_projections.h>
#include <frg/observables_common.h>
#include <runtime_config.h>

#include <sys/stat.h>  // for  stat
#include <sys/types.h> // mkdir

#include <unistd.h> // access

#include <chrono> // chorono::steady_clock::now()

#include <file_IO/file_IO_base.h>

using namespace H5; 
using namespace std; 

template <typename Model, typename State>
void MultiFileIOSBE<Model, State>::Init(resume_details_t &resume_details, std::string to_be_appended_to_out_dir_name)
{
    // create output directory
    mkdir(RuntimeConfig::OUTPUT_DIRECTORY.c_str(), 0777);

    std::string dir_name(RuntimeConfig::OUTPUT_DIRECTORY + "/dat_");

    dir_name = dir_name + FileIOBase<Model, State>::construct_name_from_config();

    dir_name = dir_name + to_be_appended_to_out_dir_name;   

    OutputDirectory() = dir_name;


#ifdef AUTORESUME_CALCULATION 

    if (access((OutputDirectory()+"/Params"+".h5").c_str(), F_OK) != -1){
	H5std_string params_file_name(OutputDirectory()+"/Params"+".h5");
	H5::H5File params_h5file( params_file_name, H5F_ACC_RDONLY );

	read_int_attribute( CurrentScaleIdx(), params_h5file, "ScalesCount" );

	params_h5file.close();

	H5std_string file_name;

	if (access((OutputDirectory()+"/final"+".h5").c_str(), F_OK) != -1){
	    file_name = H5std_string(OutputDirectory() + "/final.h5");
	    resume_details.is_finished = true;
	}else{
	    resume_details.is_finished = false;
	    CurrentScaleIdx() -= 1;
	    if (CurrentScaleIdx() < 0)
		CurrentScaleIdx() = 0;
	
	    std::cout << "Resuming after file number: " << CurrentScaleIdx() << std::endl;
	

	    file_name = H5std_string(OutputDirectory() + "/" + std::to_string(CurrentScaleIdx()) + ".h5");
	
	}

	H5::H5File h5file( file_name, H5F_ACC_RDONLY );
	
	
	H5::Group Flow_obs = h5file.openGroup("/Flow_obs");
	read_double_attribute( resume_details.start_time, Flow_obs, "t" );
	Flow_obs.close();

	h5file.close();

	if (CurrentScaleIdx() > 1){
	    H5std_string file_name2(OutputDirectory() + "/" + std::to_string(CurrentScaleIdx() - 1) + ".h5");
	    H5::H5File h5file2( file_name2, H5F_ACC_RDONLY );
	    double t2;

	    H5::Group Flow_obs2 = h5file2.openGroup("/Flow_obs");
	    read_double_attribute( t2, Flow_obs2, "t" );
	    Flow_obs2.close();

	    resume_details.init_timestep = resume_details.start_time - t2;
	    std::cout << "Timestep: " <<resume_details.init_timestep<< std::endl;
	    h5file2.close();
	}
    }else{
	std::cout << "Couldn't find output directory. Making a new one.." << std::endl;
	mkdir(dir_name.c_str(), 0777);
    }
    /*
    struct stat st;
    if (stat(dir_name, &st) == 0 && S_ISDIR(st.st_mode)){
	// dir exists
	}*/
#else
    mkdir(dir_name.c_str(), 0777);
#endif

    
    std::cout << "The directory for the output is \"" << dir_name << "\"" << std::endl;

    std::cout << "Start recording time for run...";
    StartRunTimeSeconds() = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now().time_since_epoch()).count();

}

template <typename Model, typename State>
void MultiFileIOSBE<Model, State>::WritePropertiesFile()
{
    H5std_string file_name(OutputDirectory()+"/Params"+".h5");

    H5::H5File h5file( file_name, H5F_ACC_TRUNC );
    
    H5::Group general_group( h5file.createGroup("/General") );

    write_double( BETA, general_group, "Beta");
    write( FrequenciesCount::COUNT, general_group, "COUNT");

    // write( T_START, technical_group, "T_START");
    // write( T_END, technical_group, "T_END");
    
    write( CurrentScaleIdx(), h5file, "ScalesCount");

    double time_duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now().time_since_epoch()).count() - StartRunTimeSeconds();
    write(time_duration_seconds, h5file, "RuntimeSeconds");

    H5::Group model_group( h5file.createGroup("/Model") );

    write(Model::GetName(), model_group, "Name");

    auto model_params = Model::GetParamNameValuePairs();

    for( auto par: model_params ){
	std::cout << par.first << par.second << std::endl;
	write_double(par.second, model_group, par.first);
    }
    
    H5::Group special_points_group( model_group.createGroup("/Model/Special_points") );
    for (unsigned i = 0 ; i < Model::SpecialPointsAndPaths().points_names.size(); i ++){
	int idx = Model::SpecialPointsAndPaths().points[i];
	write(idx, special_points_group, Model::SpecialPointsAndPaths().points_names[i]);
    }
    
    H5::Group special_paths_group( model_group.createGroup("/Model/Special_paths") );
    for (unsigned i = 0 ; i < Model::SpecialPointsAndPaths().paths_names.size(); i ++){
	write(Model::SpecialPointsAndPaths().paths[i], special_paths_group, Model::SpecialPointsAndPaths().paths_names[i]);
    }
    

    //    write( count,    group, "_StepCount" );
    //write( start_time,  group, "_start_time" );
    //write( end_time,    group, "_end_time" );
    //write( dt,       group, "_dt" );
    
    
    general_group.close();
    special_points_group.close();
    special_paths_group.close();
    model_group.close();

    h5file.close();
}

template <typename Model, typename State>
void MultiFileIOSBE<Model, State>::IncrementScaleIdx()
{    
    CurrentScaleIdx() += 1;
}


template <typename Model, typename State>
void MultiFileIOSBE<Model, State>::WriteCurrentScaleFile(const State &state, const std::unordered_map<std::string, void* > &observables_name_data_map)
{   
    H5std_string file_name(OutputDirectory() + "/" + std::to_string(CurrentScaleIdx()) + ".h5");
 
    H5::H5File h5file( file_name, H5F_ACC_TRUNC );

    double time_duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now().time_since_epoch()).count() - StartRunTimeSeconds();

    
    std::cout << "Writing time duration..." << std::endl;
    write(time_duration_seconds, h5file, "TimestampSeconds");

    FileIOBase<Model, State>::write_Sig( h5file, state );
    FileIOSBE<Model, State>::write_w_func( h5file, state );
    FileIOSBE<Model, State>::write_lambda_func( h5file, state );
#if !defined(SBEa_APPROXIMATION)
    FileIOSBE<Model, State>::write_M_func( h5file, state );
#endif
    
    FileIOBase<Model, State>::write_flow_observables( h5file, observables_name_data_map );
    
    h5file.close();

    std::cout << "Runtime duration in seconds: " << time_duration_seconds << std::endl;
}

template <typename Model, typename State>
 void MultiFileIOSBE<Model, State>::WriteFinalFile(const State &state, const std::unordered_map<std::string, void* > &observables_name_data_map, const std::string to_be_appended_to_file_name)
{   
    H5std_string file_name(OutputDirectory() + "/final" + to_be_appended_to_file_name + ".h5");
 
    H5::H5File h5file( file_name, H5F_ACC_TRUNC );

    double time_duration_seconds = std::chrono::duration_cast<std::chrono::seconds>(std::chrono::steady_clock::now().time_since_epoch()).count() - StartRunTimeSeconds();
    write(time_duration_seconds, h5file, "TimestampSeconds");

    FileIOBase<Model, State>::write_Sig( h5file, state );
    FileIOSBE<Model, State>::write_w_func( h5file, state );
    FileIOSBE<Model, State>::write_lambda_func( h5file, state );
#if !defined(SBEa_APPROXIMATION)
    FileIOSBE<Model, State>::write_M_func( h5file, state );
#endif
    FileIOBase<Model, State>::write_flow_observables( h5file, observables_name_data_map );

    h5file.close();

    std::cout << "Runtime duration in seconds: " << time_duration_seconds << std::endl;
}
