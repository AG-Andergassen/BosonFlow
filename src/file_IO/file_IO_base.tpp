#include <file_IO/H5Tools.h>


template <typename Model, typename State>
H5std_string FileIOBase<Model, State>::construct_name_from_config()
{
    H5std_string file_name("");
    const auto &chosen_flow = fRGFlowScheme<Model>::ChosenFlowParametrizationInfo();

    auto fname_params = Model::GetParamNameValuePairs();
    fname_params.push_back({"BETA", BETA });
    fname_params.push_back({"C", FrequenciesCount::COUNT });
    fname_params.push_back({"T_START", chosen_flow.t_start});

    for( auto parStr: fname_params )
	{
	    file_name.append("_" + parStr.first );
	    string valStr = to_string( ( long double ) parStr.second );
	    valStr.erase( valStr.find_last_not_of('0') + 1, string::npos );           /**< delete trailing zeros */
	    valStr.erase( valStr.find_last_not_of('.') + 1, string::npos );           /**< delete trailing dot */
	    replace( valStr.begin(), valStr.end(), '.', 'p');                         /**< replace dot with p */
	    file_name.append( valStr );
	}

    /**
     *   System informations
     */

    file_name.append("_"+Model::GetName());                       

    /**
     *   Regulator informations
     */
    file_name.append( "_" + fRGFlowScheme<Model>::ChosenFlowSchemeAbbrev() );                /**< Cutoff scheme */


    /**
     *   LOOP-SCHEME informations
     */

    //file_name.append("_CONSTPREFACTCUBECONSTB");

#ifdef FLOW_EQUATION_METHOD
    file_name.append("_FLOWEQN");
#endif

#ifdef SELF_CONSISTENT_METHOD
    file_name.append("_SELFCONSISTENT");
#endif


#ifdef MULTILOOP
    file_name.append("_MULTI" + MULOOP_NUM_STRING + "LOOP");      /**< Loop-number \f$ \ell > 1 \f$ */
#endif

#if defined MULTILOOP
#if defined SELFEN_FLOW || defined SELFEN_SDE || defined SELFEN_SDE_SBE_MAGNETIC || defined SELFEN_SDE_SBE_DENSITY || defined SELFEN_SDE_SBE_SUPERCONDUCTING
    file_name.append("_" + std::to_string(SELFENERGY_ITERATIONS_MAX) + "SELOOP");      /**< Loop-number \f$ \ell > 2 \f$ */
#endif

#ifdef FORCE_CALCULATE_ALL_SELFENERGY_ITERATIONS
    file_name.append("_FORCESEITER");  
#endif

#ifdef FORCE_CALCULATE_ALL_MULTILOOP_CORRECTIONS
    file_name.append("_FORCEMULTILOOPCORR");  
#endif

#endif

#ifdef SET_MIXED_BUBBLES_TO_ZERO
    file_name.append("_NOMIXEDBUBS");
#endif

#ifdef FIX_FILLING
    file_name.append("_FIXFILLING"); 
#endif

#if defined(SBEa_APPROXIMATION)
    file_name.append("_SBEaAPPROX"); 
#endif

#ifdef NO_HEDIN_VERTEX_FLOW
    file_name.append("_NOHEDIN");
#endif

#ifdef SBEb_APPROXIMATION
    file_name.append("_SBEbAPPROX");
#endif    
    
#ifdef UTILIZE_ALGEBRAIC_SYMMETRIES
    file_name.append("_ALGEBSYMM");                              /**< Algebraic symmetries used. */
#endif
#ifdef UTILIZE_LATTICE_SYMMETRIES
    file_name.append("_LATTSYMM");                               /**< Lattice symmetries used. */
#endif

    /**
     *   Formfactor informations. See params_technical.h.
     */

#ifdef BOSONISE_M
    file_name.append("_BOSEREST");
#endif

#ifdef PRECOMPUTE_STATE_PROJECTIONS
    file_name.append("_PRECOMPPROJ");
#endif

#if defined SELFEN_FLOW || defined SELFEN_SDE || defined SELFEN_SDE_SBE_MAGNETIC || defined SELFEN_SDE_SBE_DENSITY || defined SELFEN_SDE_SBE_SUPERCONDUCTING
    file_name.append("_SELFEN");
#endif
#ifdef SELFEN_FLOW
    file_name.append("_FLOW");
#endif
#ifdef SELFEN_SDE
    file_name.append("_SDE");
#endif
#ifdef SELFEN_SDE_SBE_MAGNETIC
    file_name.append("_SDE-SBE-M");
#endif
#ifdef SELFEN_SDE_SBE_DENSITY
    file_name.append("_SDE-SBE-D");
#endif
#ifdef SELFEN_SDE_SBE_SUPERCONDUCTING
    file_name.append("_SDE-SBE-SC");
#endif
#ifdef SET_OFFDIAGONAL_M_IN_FROMFACTOR_SPACE_TO_ZERO
    file_name.append("_RESTFFDIAGZERO");
#endif
#ifdef DEBUG_EULER
    file_name.append("_"+ NSTEPS_STRING+"STEPS");
#endif


    //cout << "The Filename for the output is \"" << file_name << "**.h5\"" << endl;
    return file_name; 
}


template <typename Model, typename State>
void FileIOBase<Model, State>::write_config( H5::H5File& file )
{
   H5::Group group( file.createGroup("/Config"));
    const auto &chosen_flow = fRGFlowScheme<Model>::ChosenFlowParametrizationInfo();

   std::vector<std::pair<std::string, double>> config_scalar_list = 
   { 
      { "POS_FFREQ_COUNT_SIG", FrequenciesCount::Sig::POS_FERM }, 
      { "POS_INT_RANGE", FrequencyDependenceScheme<Model>::PositiveIntegrationRange() }, 
      { "PATCH_COUNT", Model::GetRefinedMomentaCount() }, 
      { "QN_COUNT", Model::GetQuantumNumbersCount() }, 
      { "FFACTOR_COUNT", Model::GetMomentumFormFactorsCount()},
    { "T_START", chosen_flow.t_start }, 
    { "T_FIN", chosen_flow.t_final }, 
    { "INIT_T_STEP", chosen_flow.init_t_step },
      { "ERR_ABS", ERR_ABS},
      { "ERR_REL", ERR_REL},
      { "ERR_LOOP_ABS", ERR_LOOP_ABS},
      { "ERR_LOOP_REL", ERR_LOOP_REL},
   };

   for( auto conf : config_scalar_list )
      write( conf.second, group, conf.first); 
   
   std::vector<std::pair<std::string, std::string>> config_text_list = 
   {
       { "FLOW_SCHEME", fRGFlowScheme<Model>::ChosenFlowSchemeAbbrev() }, 
   }; 
   
   for( auto conf : config_text_list )
      write( conf.second, group, conf.first); 

   group.close();
}

template <typename Model, typename State>

void FileIOBase<Model, State>::write_params( H5::H5File& file )
{

   H5::Group group( file.createGroup("/Params"));

   std::vector<std::pair<std::string, double>> par_lst = Model::GetParamNameValuePairs();
   par_lst.push_back({"BETA", BETA });
    

   for( auto par : par_lst )
      write( par.second, group, par.first); 

   group.close();
}

template <typename Model, typename State>
void FileIOBase<Model, State>::write_flow_observables( H5::H5File& file, const std::unordered_map<std::string, void* > &observables_name_data_map)
{ 
    H5::Group group( file.createGroup("/Flow_obs") );
    write_flow_observables_in_subgroup(group, 
				       observables_frg_common_t<Model, State >::ObservablesListToTrack(),
				       observables_name_data_map,
				       "/Flow_obs/");

    for (auto group_observables_info : observables_frg_common_t<Model, State >::GroupObservablesListToTrack()){
	auto &[name, observables_info] = group_observables_info;

	H5::Group subgroup(group.createGroup("/Flow_obs/" + name));

	write_flow_observables_in_subgroup(subgroup, 
					   observables_info,
					   observables_name_data_map,
					   "/Flow_obs/" + name + "/");
	subgroup.close();

    }

    group.close();
}


template <typename Model, typename State>

void FileIOBase<Model, State>::write_flow_observables_in_subgroup( H5::Group &group, const std::vector<std::tuple< std::string, FlowObservableType, unsigned> > &observables_info, const std::unordered_map<std::string, void* > &observables_name_data_map, std::string group_directory)
{
    for (auto &observable_info: observables_info ){
	std::string observable_name = std::get<0>(observable_info);
	FlowObservableType observable_type = std::get<1>(observable_info);
	// 10 = "/Flow_obs/".size()
	void *data_ptr = observables_name_data_map.at(group_directory.substr(10) + observable_name);
	
	switch (observable_type){
	    /*
	case FlowObservableType::INT:
	    write( *static_cast<std::vector<int>*>(data_ptr), group, observable_name ); 
	    break;
	    */
	case FlowObservableType::DOUBLE:
#ifdef MULTIFILE_OUTPUT
	    if (static_cast<std::vector<double>*>(data_ptr)->size() > 0)
		write( static_cast<std::vector<double>*>(data_ptr)->at(0), group, observable_name ); 	   
#else
	    write( *static_cast<std::vector<double>*>(data_ptr), group, observable_name ); 	   
#endif 
	    break;
	
	    
    
	case FlowObservableType::COMPLEX:
#ifdef MULTIFILE_OUTPUT
	    if (static_cast<std::vector<std::complex<double> >*>(data_ptr)->size() > 0)
		write( static_cast<std::vector<std::complex<double> >*>(data_ptr)->at(0), group, "_" + observable_name ); 	    
#else
	    write( *static_cast<std::vector<std::complex<double> >*>(data_ptr), group, "_" + observable_name );
#endif
	    break;

	case FlowObservableType::STRING:
#ifdef MULTIFILE_OUTPUT
	    if (static_cast<std::vector<std::string>*>(data_ptr)->size() > 0)
		write( static_cast<std::vector<std::string>*>(data_ptr)->at(0), group, observable_name ); 	 
#else
	    write( *static_cast<std::vector<std::string>*>(data_ptr), group, observable_name ); 	
#endif
	    break;

	case FlowObservableType::GF_MULTIINDEX:
#ifdef MULTIFILE_OUTPUT
	    if (static_cast<std::vector<std::vector<int> >*>(data_ptr)->size() > 0)
		write_vec_int( static_cast<std::vector<std::vector<int> >*>(data_ptr)->at(0), group, observable_name );
#else
	    write( *static_cast<std::vector<std::vector<int> >*>(data_ptr), group, observable_name ); 	
#endif   
	    break;

	    

	case FlowObservableType::COMPLEX_GF_CONTAINER:
	    {
	    unsigned gf_ndims = std::get<2>(observable_info);

#ifdef MULTIFILE_OUTPUT

#define WriteGFWithNDIM(X)\
	    if (gf_ndims == X){\
		auto &vector_of_gf_containers = *static_cast<std::vector<gf<std::complex<double>, X> >* >(data_ptr); \
		if (vector_of_gf_containers.size() > 0)\
		    write_complex_valued_gf_object(vector_of_gf_containers.at(0), group, "_" + observable_name); \
	    }

#else

	    H5::Group subGroup(group.createGroup(group_directory + observable_name));

#define WriteGFWithNDIM(X)\
	    if (gf_ndims == X){\
		auto &vector_of_gf_containers = *static_cast<std::vector<gf<std::complex<double>, X> >* >(data_ptr); \
		for (unsigned lam_index = 0; lam_index < vector_of_gf_containers.size(); lam_index ++) \
		    write_complex_valued_gf_object(vector_of_gf_containers[lam_index], subGroup, "_" + to_string(lam_index) ); \
	    }
#endif
	    WriteGFWithNDIM(1); WriteGFWithNDIM(2); WriteGFWithNDIM(3);
	    WriteGFWithNDIM(4); WriteGFWithNDIM(5); WriteGFWithNDIM(6);
	    WriteGFWithNDIM(7); WriteGFWithNDIM(8); WriteGFWithNDIM(9);
	    WriteGFWithNDIM(10); WriteGFWithNDIM(11); WriteGFWithNDIM(12);
	    WriteGFWithNDIM(13); WriteGFWithNDIM(14); WriteGFWithNDIM(15);
	    WriteGFWithNDIM(16); WriteGFWithNDIM(17); WriteGFWithNDIM(18);
	    WriteGFWithNDIM(19); WriteGFWithNDIM(20); WriteGFWithNDIM(21);
	    
#ifndef MULTIFILE_OUTPUT
	    subGroup.close();
#endif
	    }
      	    break;
	default:;
	}		
    }
}


// write functions for FFT(GG), bubble and Gvec_real ---- not used in normal full run
template <typename Model, typename State>
void FileIOBase<Model, State>::write_bubble( H5::H5File& file, const gf_bubble_mat_t<Model> bubble_pp, const gf_bubble_mat_t<Model> bubble_ph )
{
    gf_bubble_t<Model> b_pp;
    gf_bubble_t<Model> b_ph;
	const int pos_w_bfreq_count = FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count;
    // todo: use proper frequency bounds variables
        for(int W=-pos_w_bfreq_count; W<pos_w_bfreq_count; W++)
            for(int w=-(FrequencyDependenceScheme<Model>::PositiveIntegrationRange()+pos_w_bfreq_count/2); w <FrequencyDependenceScheme<Model>::PositiveIntegrationRange()+pos_w_bfreq_count/2; w++)
	 for(int K=0; K<Model::GetRefinedMomentaCount(); K++)
	    for(int m=0; m<Model::GetMomentumFormFactorsCount(); m++)
	       for(int n=0; n<Model::GetMomentumFormFactorsCount(); n++)
		  for(int s1=0; s1<Model::GetQuantumNumbersCount(); s1++)
		  for(int s2=0; s2<Model::GetQuantumNumbersCount(); s2++)
		  for(int s3=0; s3<Model::GetQuantumNumbersCount(); s3++)
		  for(int s4=0; s4<Model::GetQuantumNumbersCount(); s4++)
	       {
	       b_pp[W][w][K][m][n][s1][s2][s3][s4] = bubble_pp[W][w][m][n][s1][s2][s3][s4](K);
	       b_ph[W][w][K][m][n][s1][s2][s3][s4] = bubble_ph[W][w][m][n][s1][s2][s3][s4](K);
	       }


   H5::Group group( file.createGroup("/bubble") );

   write_complex_valued_gf_object( b_pp, group, "_PP" ); 
   write_complex_valued_gf_object( b_ph, group, "_PH" ); 

   group.close();
}

// write projection matrices ---- not used in normal full run
template <typename Model, typename State>
void FileIOBase<Model, State>::write_projection( H5::H5File& file )
{
   H5::Group group( file.createGroup("/Proj") );
   
   write_complex_valued_gf_object( TUProjections<Model>::Matrix_xph_to_pp(), group, "_XPHtoPP" ); 
   write_complex_valued_gf_object( TUProjections<Model>::Matrix_pp_to_ph(), group, "_PPtoPH" ); 
   write_complex_valued_gf_object( TUProjections<Model>::Matrix_ph_to_xph(), group, "_PHtoXPH" ); 

   group.close();
}

template <typename Model, typename State>
void FileIOBase<Model, State>::write_Gfunc( H5::H5File& file, const gf_1p_mat_real_t<Model> Gvec_real, const gf_1p_mat_real_t<Model> Svec_real )
{
    gf_1p_real_t<Model> Gvec_real_plot;
    gf_1p_real_t<Model> Svec_real_plot;
    const int pos_g_ffreq_count = FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count;

    for(int W=-pos_g_ffreq_count; W<pos_g_ffreq_count; W++)
	for(int K=0; K<Model::GetFineMomentaCount(); K++)
	    for(int s1=0; s1<QN_COUNT; s1++)
		for(int s2=0; s2<QN_COUNT; s2++)
		    {
			Gvec_real_plot[W][K][s1][s2] = Gvec_real[W][s1][s2](K);
			Svec_real_plot[W][K][s1][s2] = Svec_real[W][s1][s2](K);
		    }


    H5::Group group( file.createGroup("/Gvec") );

    write_complex_valued_gf_object( Gvec_real_plot, group, "_Gvec" ); 
    write_complex_valued_gf_object( Svec_real_plot, group, "_Svec" ); 
    group.close();
}


template <typename Model, typename State>
void FileIOBase<Model, State>::write_Sig( H5::H5File& file, const State& state_vec )
{
    H5::Group group( file.createGroup("/Sig") );
    write_complex_valued_gf_object( state_vec.gf_Sig(), group ); 
    write_fermionic_matsubara_space( fermionic_matsubara_space_t( FrequenciesCount::Sig::POS_FERM, 2.0*M_PI / BETA ), group );
    
    write( Model::MomentumGrid().ptr->get_points_as_std_vectors(), Model::dim, Model::GetCoarseMomentaCount(), group );

    group.close();
}


