#include <frg/sbe/symmetries_mfrg.h>

#include <algorithm>

template <typename Model>
void SymmetriesmfRGSBE<Model>::Init()
{
    Symmetries1lfRGSBE<Model>::Init();

    if (SymmetriesmfRGSBE<Model>::AreLatticeSymmetriesIncluded() || SymmetriesmfRGSBE<Model>::AreAlgebraicSymmetriesIncluded())
	return;

#ifdef UTILIZE_ALGEBRAIC_SYMMETRIES
    IncludeAlgebraicSymmetries();
#endif
#ifdef UTILIZE_LATTICE_SYMMETRIES
    IncludeLatticeSymmetries();
#endif
    CalculateIndexEquivalenceClasses();
}


template <typename Model>
void SymmetriesmfRGSBE<Model>::CalculateIndexEquivalenceClasses()
{
    // todo: if had been allocated before, should delete first..
    std::cout << " Constructing left corrections of SBE hedin vertex (lambda) left corrections symmetries... " << std::endl;
    /*    
    IdxEquivClasses_lambda_sc_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_lambda_t<Model>(FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1), GroupActions_lambda_sc_LeftCorrections());
    IdxEquivClasses_lambda_d_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_lambda_t<Model>(FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1), GroupActions_lambda_d_LeftCorrections());
    IdxEquivClasses_lambda_m_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_lambda_t<Model>(FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1), GroupActions_lambda_m_LeftCorrections());
    */
    
    auto sampling_freqs_lambda = {
	std::make_pair(int(ILAMBDA::W), FrequencyDependenceScheme<Model>::lambda_BosonicGrid().sampling_freqs), 
	std::make_pair(int(ILAMBDA::w), FrequencyDependenceScheme<Model>::lambda_FermionicGrid().sampling_freqs)
    };    

    IdxEquivClasses_lambda_sc_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_lambda_t<Model>(), GroupActions_lambda_sc_LeftCorrections(), sampling_freqs_lambda);
    IdxEquivClasses_lambda_d_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_lambda_t<Model>(), GroupActions_lambda_d_LeftCorrections(), sampling_freqs_lambda);
    IdxEquivClasses_lambda_m_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_lambda_t<Model>(), GroupActions_lambda_m_LeftCorrections(), sampling_freqs_lambda);


    auto sampling_freqs_M = {
	std::make_pair(int(IM::W), FrequencyDependenceScheme<Model>::M_BosonicGrid().sampling_freqs), 
	std::make_pair(int(IM::w), FrequencyDependenceScheme<Model>::M_FermionicGrid().sampling_freqs)
    };

    std::cout << " Constructing left corrections of SBE rest function (M) left corrections symmetries... " << std::endl;
    /*IdxEquivClasses_M_sc_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_M_t<Model>(FrequencyDependenceScheme<Model>::M_BosonicGrid().positive_freqs_count, FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1), GroupActions_M_sc_LeftCorrections());
    IdxEquivClasses_M_d_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_M_t<Model>(FrequencyDependenceScheme<Model>::M_BosonicGrid().positive_freqs_count, FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1), GroupActions_M_d_LeftCorrections());
    IdxEquivClasses_M_m_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_M_t<Model>(FrequencyDependenceScheme<Model>::M_BosonicGrid().positive_freqs_count, FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, FrequenciesCount::bubble::POS_FERM+1), GroupActions_M_m_LeftCorrections());*/
    IdxEquivClasses_M_sc_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_M_t<Model>(), GroupActions_M_sc_LeftCorrections(), sampling_freqs_M);
    IdxEquivClasses_M_d_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_M_t<Model>(), GroupActions_M_d_LeftCorrections(), sampling_freqs_M);
    IdxEquivClasses_M_m_LeftCorrections_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_M_t<Model>(), GroupActions_M_m_LeftCorrections(), sampling_freqs_M);

}

template <typename Model>
void SymmetriesmfRGSBE<Model>::IncludeAlgebraicSymmetries()
{
    GroupActions_lambda_sc_LeftCorrections().push_back(cconj_timerev_lambda_sc<Model>);
    GroupActions_lambda_d_LeftCorrections().push_back(cconj_timerev_lambda_d<Model>);
    GroupActions_lambda_m_LeftCorrections().push_back(cconj_timerev_lambda_m<Model>);

    // M (left corrections)
    // Left corrections of M have no algebraic symmetries

    SymmetriesmfRGSBE<Model>::AreAlgebraicSymmetriesIncluded() = true;
}

template <typename Model>
void SymmetriesmfRGSBE<Model>::IncludeLatticeSymmetries()
{
    for (unsigned g_idx = 0; g_idx < Model::RealLattice().point_group_ptr->m_group_size; g_idx ++ ){
	
	// for the SBE-rest function M
	GroupActions_lambda_sc_LeftCorrections().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_lambda(g_idx));
	GroupActions_lambda_d_LeftCorrections().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_lambda(g_idx));
	GroupActions_lambda_m_LeftCorrections().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_lambda(g_idx));

	GroupActions_M_sc_LeftCorrections().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_M(g_idx));
	GroupActions_M_d_LeftCorrections().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_M(g_idx));
	GroupActions_M_m_LeftCorrections().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_M(g_idx));
    }

    SymmetriesmfRGSBE<Model>::AreLatticeSymmetriesIncluded() = true;
}


template <typename Model>
void SymmetriesmfRGSBE<Model>::CleanUp()
{
    Symmetries1lfRGSBE<Model>::CleanUp();

    delete IdxEquivClasses_lambda_sc_LeftCorrections_Ptr();
    delete IdxEquivClasses_lambda_d_LeftCorrections_Ptr();
    delete IdxEquivClasses_lambda_m_LeftCorrections_Ptr();

    delete IdxEquivClasses_M_sc_LeftCorrections_Ptr();
    delete IdxEquivClasses_M_d_LeftCorrections_Ptr();
    delete IdxEquivClasses_M_m_LeftCorrections_Ptr();
}


// instantiate class for the particular models
#define Instantiate(MODEL) template class SymmetriesmfRGSBE<MODEL>;
WITH_MODELS(Instantiate)


