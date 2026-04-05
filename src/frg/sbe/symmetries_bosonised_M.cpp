#include <frg/sbe/symmetries_bosonised_M.h>

#include <algorithm>

template <typename Model>
void SymmetriesfRGSBEBosonised_M<Model>::Init()
{
    Symmetries1lfRGSBE<Model>::Init();

    if (SymmetriesfRGSBEBosonised_M<Model>::AreLatticeSymmetriesIncluded() || SymmetriesfRGSBEBosonised_M<Model>::AreAlgebraicSymmetriesIncluded())
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
void SymmetriesfRGSBEBosonised_M<Model>::CalculateIndexEquivalenceClasses()
{
    // todo: if had been allocated before, should delete first..
    std::cout << " Constructing SBE bosonised M objects symmetries... " << std::endl;
    
    auto sampling_freqs_w_M = {
	std::make_pair(int(IW_M::W), FrequencyDependenceScheme<Model>::lambda_BosonicGrid().sampling_freqs)
    };    
    IdxEquivClasses_w_M_sc_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_w_M_t<Model>(), GroupActions_w_M_sc(), sampling_freqs_w_M);
    IdxEquivClasses_w_M_d_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_w_M_t<Model>(), GroupActions_w_M_d(), sampling_freqs_w_M);
    IdxEquivClasses_w_M_m_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_w_M_t<Model>(), GroupActions_w_M_m(), sampling_freqs_w_M);


    auto sampling_freqs_lambda_M = {
	std::make_pair(int(ILAMBDA_M::W), FrequencyDependenceScheme<Model>::M_BosonicGrid().sampling_freqs), 
	std::make_pair(int(ILAMBDA_M::w), FrequencyDependenceScheme<Model>::M_FermionicGrid().sampling_freqs)
    };
    IdxEquivClasses_lambda_M_sc_Ptr() = new symmetry_grp_t<dcomplex,5>(gf_lambda_M_t<Model>(), GroupActions_lambda_M_sc(), sampling_freqs_lambda_M);
    IdxEquivClasses_lambda_M_d_Ptr() = new symmetry_grp_t<dcomplex,5>(gf_lambda_M_t<Model>(), GroupActions_lambda_M_d(), sampling_freqs_lambda_M);
    IdxEquivClasses_lambda_M_m_Ptr() = new symmetry_grp_t<dcomplex,5>(gf_lambda_M_t<Model>(), GroupActions_lambda_M_m(), sampling_freqs_lambda_M);

}

template <typename Model>
void SymmetriesfRGSBEBosonised_M<Model>::IncludeAlgebraicSymmetries()
{
    // none for now... todo: add more

    SymmetriesfRGSBEBosonised_M<Model>::AreAlgebraicSymmetriesIncluded() = true;
}

template <typename Model>
void SymmetriesfRGSBEBosonised_M<Model>::IncludeLatticeSymmetries()
{
    for (unsigned g_idx = 0; g_idx < Model::RealLattice().point_group_ptr->m_group_size; g_idx ++ ){
	
	// for the bosonised M bosonic propagator w_M
	GroupActions_w_M_sc().push_back(SymmetriesfRGSBEBosonised_M<Model>::MakeLatticeSymmetryMap_w_M(g_idx));
	GroupActions_w_M_d().push_back(SymmetriesfRGSBEBosonised_M<Model>::MakeLatticeSymmetryMap_w_M(g_idx));
	GroupActions_w_M_m().push_back(SymmetriesfRGSBEBosonised_M<Model>::MakeLatticeSymmetryMap_w_M(g_idx));

	// for the bosonised M hedin vertex lambda_M
	GroupActions_lambda_M_sc().push_back(SymmetriesfRGSBEBosonised_M<Model>::MakeLatticeSymmetryMap_lambda_M(g_idx));
	GroupActions_lambda_M_d().push_back(SymmetriesfRGSBEBosonised_M<Model>::MakeLatticeSymmetryMap_lambda_M(g_idx));
	GroupActions_lambda_M_m().push_back(SymmetriesfRGSBEBosonised_M<Model>::MakeLatticeSymmetryMap_lambda_M(g_idx));
    }

    SymmetriesfRGSBEBosonised_M<Model>::AreLatticeSymmetriesIncluded() = true;
}


template <typename Model>
void SymmetriesfRGSBEBosonised_M<Model>::CleanUp()
{
    Symmetries1lfRGSBE<Model>::CleanUp();

    delete IdxEquivClasses_w_M_sc_Ptr();
    delete IdxEquivClasses_w_M_d_Ptr();
    delete IdxEquivClasses_w_M_m_Ptr();

    delete IdxEquivClasses_lambda_M_sc_Ptr();
    delete IdxEquivClasses_lambda_M_d_Ptr();
    delete IdxEquivClasses_lambda_M_m_Ptr();

}

template <typename Model>
std::function<operation(idx_w_M_t<Model>& )> SymmetriesfRGSBEBosonised_M<Model>::MakeLatticeSymmetryMap_w_M(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_w_M_t<Model>& w_M_idx){

	w_M_idx(IW_M::K) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][w_M_idx(IW_M::K)];	

	bool with_minus, with_minus2;
	// note that the form-factors are based on the coarse lattice (e.g. bond form factors from the coarse lattice)
	std::tie(w_M_idx(IW_M::m), with_minus) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][w_M_idx(IW_M::m)];

	std::tie(w_M_idx(IW_M::n), with_minus2) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][w_M_idx(IW_M::n)];

	// exclusive OR so that it's only true if one of the two signs , but not both (since minus times a minus is a plus), is a minus
	bool overall_with_minus = with_minus ^ with_minus2;

	
	return operation(overall_with_minus, false );
    };

    return func;
}

template <typename Model>
std::function<operation(idx_lambda_M_t<Model>& )> SymmetriesfRGSBEBosonised_M<Model>::MakeLatticeSymmetryMap_lambda_M(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_lambda_M_t<Model>& lambda_M_idx){

	lambda_M_idx(ILAMBDA_M::K) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][lambda_M_idx(ILAMBDA_M::K)];	

	bool with_minus, with_minus2;
	// note that the form-factors are based on the coarse lattice (e.g. bond form factors from the coarse lattice)
	std::tie(lambda_M_idx(ILAMBDA_M::m), with_minus) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][lambda_M_idx(ILAMBDA_M::m)];

	std::tie(lambda_M_idx(ILAMBDA_M::n), with_minus2) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][lambda_M_idx(ILAMBDA_M::n)];

	// exclusive OR so that it's only true if one of the two signs , but not both (since minus times a minus is a plus), is a minus
	bool overall_with_minus = with_minus ^ with_minus2;

	
	return operation(overall_with_minus, false );
    };

    return func;
}



// instantiate class for the particular models
#define Instantiate(MODEL) template class SymmetriesfRGSBEBosonised_M<MODEL>;
WITH_MODELS(Instantiate)


