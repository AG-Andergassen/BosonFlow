

template <typename Model, typename State>
void FileIOSBE<Model, State>::write_final_file(const State &state, const std::unordered_map<std::string, void* > &observables_name_data_map, const std::string to_be_appended_to_file_name)
{

    H5std_string file_name("dat_");
    file_name.append(FileIOBase<Model, State>::construct_name_from_config());
    file_name.append(static_cast<H5std_string>(to_be_appended_to_file_name));
    file_name.append(".h5");

    H5::H5File file( file_name, H5F_ACC_TRUNC );
    
    write_common_data(file, state, observables_name_data_map);
    file.close();
}


template <typename Model, typename State>
void FileIOSBE<Model, State>::write_checkpoint_file( const State& state, const std::unordered_map<std::string, void* > &observables_name_data_map, size_t integration_step_number, double start_time, double end_time, double dt, size_t AorB )
{
   H5std_string file_name("Checkpoint_");
   file_name.append(FileIOBase<Model, State>::construct_name_from_config());
   cout << "AorB = " << AorB << endl; 
   if( AorB%2 == 0 ){ 
       file_name.append("A");
   }
   else{
       file_name.append("B");
   }
   file_name.append(".h5");
   
   cout << "Creating checkpoint of current flow results. " << endl ;
   
   H5::H5File file( file_name, H5F_ACC_TRUNC );
   
   write_common_data(file, state, observables_name_data_map);

   // I'm not yet sure if this line is necessary
   FileIOBase<Model, State>::write_CheckpointInfo( file, static_cast< double >(integration_step_number), start_time, end_time, dt);
   
   file.close();
   cout << "Checkpoint file creation for start_time = " << start_time << " successful." << endl;
}

template <typename Model, typename State>
void FileIOSBE<Model, State>::read_state_from_checkpoint_file(const std::string input_file_name, State &state)
{

    using namespace H5;

    H5::H5File input_file(static_cast<H5std_string>(input_file_name), H5F_ACC_RDONLY);
    std::cout << "Got file" << std::endl;
    std::cout << input_file_name << std::endl;

    H5::Group Sig_group = input_file.openGroup("/Sig") ;
    read( state.gf_Sig(), Sig_group );
    std::cout << "Read self-energy data successfully.." << std::endl;
    Sig_group.close();

    H5::Group w_group =  input_file.openGroup("/w_func");
    read( state.gf_w_sc(), w_group,"_SC");
    read( state.gf_w_d(), w_group,"_D");
    read( state.gf_w_m(), w_group,"_M");
    std::cout << "Read w data successfully.." << std::endl;
    w_group.close();

    H5::Group lambda_group =  input_file.openGroup("/lambda_func");
    read( state.gf_lambda_sc(), lambda_group,"_SC");
    read( state.gf_lambda_d(), lambda_group,"_D");
    read( state.gf_lambda_m(), lambda_group,"_M");
    std::cout << "Read lambda data successfully.." << std::endl;
    lambda_group.close();

#if !defined(SBEa_APPROXIMATION)
    H5::Group M_group =  input_file.openGroup("/M_func");
    read( state.gf_M_sc(), M_group,"_SC");
    read( state.gf_M_d(), M_group,"_D");
    read( state.gf_M_m(), M_group,"_M");
    std::cout << "Read M data successfully.." << std::endl;
    M_group.close();
#endif

    Group Flow_obs = input_file.openGroup("/Flow_obs");
    read_double_attribute( state.gf_delta_mu(), Flow_obs, "delta_mu" );
    read_double_attribute( state.m_d_delta_mu_over_dt, Flow_obs, "delta_mu_dot" );
    Flow_obs.close();
    std::cout << "Read delta_mu successfully.." << std::endl;
    std::cout << state.gf_delta_mu() << std::endl;
    std::cout << "Read d_delta_mu_over_dt successfully.." << std::endl;
    std::cout << state.m_d_delta_mu_over_dt << std::endl;
}


template <typename Model, typename State>
void FileIOSBE<Model, State>::write_common_data(H5::H5File &file, const State &state, const std::unordered_map<std::string, void* > &observables_name_data_map)
{
    FileIOBase<Model, State>::write_config( file );
    FileIOBase<Model, State>::write_params( file );
    FileIOBase<Model, State>::write_projection( file ); // for debugging purposes
    FileIOBase<Model, State>::write_Sig( file, state );
    write_w_func( file, state );
    write_lambda_func( file, state );
#if !defined(SBEa_APPROXIMATION)
    write_M_func( file, state );
#endif
    FileIOBase<Model, State>::write_flow_observables( file, observables_name_data_map );

    H5::Group group( file.createGroup("/Special_points") );
    for (unsigned i = 0 ; i < Model::SpecialPointsAndPaths().points_names.size(); i ++){
	int idx = Model::SpecialPointsAndPaths().points[i];
	write(idx, group, Model::SpecialPointsAndPaths().points_names[i]);
    }
    group.close();

    H5::Group group2( file.createGroup("/Special_paths") );
    for (unsigned i = 0 ; i < Model::SpecialPointsAndPaths().paths_names.size(); i ++){
	write(Model::SpecialPointsAndPaths().paths[i], group2, Model::SpecialPointsAndPaths().paths_names[i]);
    }
    group2.close();
}


template <typename Model, typename State>
void FileIOSBE<Model, State>::write_w_func( H5::H5File& file, const State& state )
{
   H5::Group group( file.createGroup("/w_func") );

   write_complex_valued_gf_object( state.gf_w_sc(), group, "_SC" ); 
   write_complex_valued_gf_object( state.gf_w_m(), group, "_M" ); 
   write_complex_valued_gf_object( state.gf_w_d(), group, "_D" ); 

    write_bosonic_matsubara_space( bosonic_matsubara_space_t( FrequenciesCount::w::POS_BOS, 2.0*M_PI / BETA ), group );
   
   write( Model::MomentumGrid().ptr->get_points_as_std_vectors(), Model::dim, Model::GetRefinedMomentaCount(), group );

   group.close();
}


template <typename Model, typename State>
void FileIOSBE<Model, State>::write_lambda_func( H5::H5File& file, const State& state )
{
   H5::Group group( file.createGroup("/lambda_func") );

   write_complex_valued_gf_object( state.gf_lambda_sc(), group, "_SC" ); 
   write_complex_valued_gf_object( state.gf_lambda_m(), group, "_M" ); 
   write_complex_valued_gf_object( state.gf_lambda_d(), group, "_D" ); 

    write_bosonic_matsubara_space( bosonic_matsubara_space_t( FrequenciesCount::lambda::POS_BOS, 2.0*M_PI / BETA ), group );
    write_fermionic_matsubara_space( fermionic_matsubara_space_t( FrequenciesCount::lambda::POS_FERM, 2.0*M_PI / BETA ), group );
   
   write( Model::MomentumGrid().ptr->get_points_as_std_vectors(), Model::dim, Model::GetRefinedMomentaCount(), group );

   group.close();
}


template <typename Model, typename State>
void FileIOSBE<Model, State>::write_M_func( H5::H5File& file, const State& state )
{
   H5::Group group( file.createGroup("/M_func") );

#ifndef BOSONISE_M
   write_complex_valued_gf_object( state.gf_M_sc(), group, "_SC" );
   write_complex_valued_gf_object( state.gf_M_m(), group, "_M" );
   write_complex_valued_gf_object( state.gf_M_d(), group, "_D" );
#endif

    write_bosonic_matsubara_space( bosonic_matsubara_space_t( FrequenciesCount::M::POS_BOS, 2.0*M_PI / BETA ), group );
    write_fermionic_matsubara_space( fermionic_matsubara_space_t( FrequenciesCount::M::POS_FERM, 2.0*M_PI / BETA ), group );
   
   write( Model::MomentumGrid().ptr->get_points_as_std_vectors(), Model::dim, Model::GetRefinedMomentaCount(), group );

#ifdef BOSONISE_M
   H5::Group subgroup( file.createGroup("/M_func/Bosonised") );

   write_complex_valued_gf_object( state.gf_wM_sc_plus(), subgroup, "_w_SC_plus" );
   write_complex_valued_gf_object( state.gf_wM_m_plus(), subgroup, "_w_M_plus" );
   write_complex_valued_gf_object( state.gf_wM_d_plus(), subgroup, "_w_D_plus" );

   write_complex_valued_gf_object( state.gf_wM_sc_minus(), subgroup, "_w_SC_minus" );
   write_complex_valued_gf_object( state.gf_wM_m_minus(), subgroup, "_w_M_minus" );
   write_complex_valued_gf_object( state.gf_wM_d_minus(), subgroup, "_w_D_minus" );

   write_complex_valued_gf_object( state.gf_lambdaM_sc_plus(), subgroup, "_lambda_SC_plus" );
   write_complex_valued_gf_object( state.gf_lambdaM_m_plus(), subgroup, "_lambda_M_plus" );
   write_complex_valued_gf_object( state.gf_lambdaM_d_plus(), subgroup, "_lambda_D_plus" );

   write_complex_valued_gf_object( state.gf_lambdaM_sc_minus(), subgroup, "_lambda_SC_minus" );
   write_complex_valued_gf_object( state.gf_lambdaM_m_minus(), subgroup, "_lambda_M_minus" );
   write_complex_valued_gf_object( state.gf_lambdaM_d_minus(), subgroup, "_lambda_D_minus" );


   subgroup.close();
#endif

   group.close();
}
