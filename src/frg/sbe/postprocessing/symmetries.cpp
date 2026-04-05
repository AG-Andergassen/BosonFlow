#include <frg/sbe/postprocessing/symmetries.h>

#include <algorithm>
#include <utility>

template <typename Model>
void SymmetriesfRGSBEPostprocessing<Model>::Init()
{
    SymmetriesfRGCommon<Model>::Init();

    if (SymmetriesfRGSBEPostprocessing<Model>::AreAlgebraicSymmetriesIncluded() || SymmetriesfRGSBEPostprocessing<Model>::AreLatticeSymmetriesIncluded())
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
void SymmetriesfRGSBEPostprocessing<Model>::CalculateIndexEquivalenceClasses()
{
    // todo: if had been allocated before, should delete first..

    std::cout << " Constructing postprocessing susceptibilities symmetries... " << std::endl;
    IdxEquivClasses_susc_sc_Ptr() = new symmetry_grp_t<dcomplex,8>(gf_susc_t<Model>(), GroupActions_susc_sc());
    IdxEquivClasses_susc_d_Ptr() = new symmetry_grp_t<dcomplex,8>(gf_susc_t<Model>(), GroupActions_susc_d());
    IdxEquivClasses_susc_m_Ptr() = new symmetry_grp_t<dcomplex,8>(gf_susc_t<Model>(), GroupActions_susc_m());

    std::cout << " Constructing postprocessing polarisation symmetries... " << std::endl;
    IdxEquivClasses_polarisation_sc_Ptr() = new symmetry_grp_t<dcomplex,2>(gf_polarisation_t<Model>(), GroupActions_polarisation_sc());
    IdxEquivClasses_polarisation_d_Ptr() = new symmetry_grp_t<dcomplex,2>(gf_polarisation_t<Model>(), GroupActions_polarisation_d());
    IdxEquivClasses_polarisation_m_Ptr() = new symmetry_grp_t<dcomplex,2>(gf_polarisation_t<Model>(), GroupActions_polarisation_m());

    std::cout << " Constructing postprocessing lambda symmetries... " << std::endl;
    IdxEquivClasses_static_lambda_sc_Ptr() = new symmetry_grp_t<dcomplex,2>(gf_static_lambda_t<Model>(), GroupActions_static_lambda_sc());
    IdxEquivClasses_static_lambda_d_Ptr() = new symmetry_grp_t<dcomplex,2>(gf_static_lambda_t<Model>(), GroupActions_static_lambda_d());
    IdxEquivClasses_static_lambda_m_Ptr() = new symmetry_grp_t<dcomplex,2>(gf_static_lambda_t<Model>(), GroupActions_static_lambda_m());

}

template <typename Model>
void SymmetriesfRGSBEPostprocessing<Model>::IncludeAlgebraicSymmetries()
{
    // suscs
    GroupActions_susc_sc().push_back(cconj_timerev_susc_sc<Model>);
    GroupActions_susc_d().push_back(cconj_timerev_susc_d<Model>);
    GroupActions_susc_m().push_back(cconj_timerev_susc_m<Model>);

    GroupActions_susc_sc().push_back(timerev_susc_sc<Model>);
    GroupActions_susc_d().push_back(hmirror_timerev_susc_d<Model>);
    GroupActions_susc_m().push_back(timerev_susc_m<Model>);

    SymmetriesfRGSBEPostprocessing<Model>::AreAlgebraicSymmetriesIncluded() = true;
}

template <typename Model>
void SymmetriesfRGSBEPostprocessing<Model>::IncludeLatticeSymmetries()
{
    for (unsigned g_idx = 0; g_idx < Model::RealLattice().point_group_ptr->m_group_size; g_idx ++ ){
	
	// for the postprocessing susceptibility
	GroupActions_susc_sc().push_back(SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_susc(g_idx));
	GroupActions_susc_d().push_back(SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_susc(g_idx));
	GroupActions_susc_m().push_back(SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_susc(g_idx));

	// for the postprocessing polarisation
	GroupActions_polarisation_sc().push_back(SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_polarisation(g_idx));
	GroupActions_polarisation_d().push_back(SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_polarisation(g_idx));
	GroupActions_polarisation_m().push_back(SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_polarisation(g_idx));

	// for the postprocessing lambda
	GroupActions_static_lambda_sc().push_back(SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_static_lambda(g_idx));
	GroupActions_static_lambda_d().push_back(SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_static_lambda(g_idx));
	GroupActions_static_lambda_m().push_back(SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_static_lambda(g_idx));
    }

    SymmetriesfRGSBEPostprocessing<Model>::AreLatticeSymmetriesIncluded() = true;
}


template <typename Model>
void SymmetriesfRGSBEPostprocessing<Model>::CleanUp()
{
    SymmetriesfRGCommon<Model>::CleanUp();

    delete IdxEquivClasses_susc_sc_Ptr();
    delete IdxEquivClasses_susc_d_Ptr();
    delete IdxEquivClasses_susc_m_Ptr();

    delete IdxEquivClasses_polarisation_sc_Ptr();
    delete IdxEquivClasses_polarisation_d_Ptr();
    delete IdxEquivClasses_polarisation_m_Ptr();

    delete IdxEquivClasses_static_lambda_sc_Ptr();
    delete IdxEquivClasses_static_lambda_d_Ptr();
    delete IdxEquivClasses_static_lambda_m_Ptr();
}


//        LATTICE SYMMETRIES


/**
 * \brief lattice symmetry for the postproessing susceptibilities
 *
 * same for all channels.
 */



template <typename Model>
std::function<operation(idx_susc_t<Model>& )> SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_susc(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_susc_t<Model>& susc_idx){

	susc_idx(ISUSC::K) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][susc_idx(ISUSC::K)];

	bool with_minus, with_minus2;
	// note that the form-factors are based on the coarse lattice (e.g. bond form factors from the coarse lattice)
	std::tie(susc_idx(ISUSC::n_in), with_minus) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][susc_idx(ISUSC::n_in)];

	std::tie(susc_idx(ISUSC::n_out), with_minus2) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][susc_idx(ISUSC::n_out)];

	// exclusive OR so that it's only true if one of the two signs , but not both (since minus times a minus is a plus), is a minus
	bool overall_with_minus = with_minus ^ with_minus2;
	
	return operation(overall_with_minus, false );
    };

    return func;
}


template <typename Model>
std::function<operation(idx_polarisation_t<Model>& )> SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_polarisation(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_polarisation_t<Model>& polarisation_idx){

	polarisation_idx(IPOLARISATION::K) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][polarisation_idx(IPOLARISATION::K)];

	// exclusive OR so that it's only true if one of the two signs , but not both (since minus times a minus is a plus), is a minus
	bool overall_with_minus = false;
	
	return operation(overall_with_minus, false );
    };

    return func;
}


template <typename Model>
std::function<operation(idx_static_lambda_t<Model>& )> SymmetriesfRGSBEPostprocessing<Model>::MakeLatticeSymmetryMap_static_lambda(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_static_lambda_t<Model>& static_lambda_idx){

	static_lambda_idx(ISTATICLAMBDA::K) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][static_lambda_idx(ISTATICLAMBDA::K)];

	// exclusive OR so that it's only true if one of the two signs , but not both (since minus times a minus is a plus), is a minus
	bool overall_with_minus = false;
	
	return operation(overall_with_minus, false );
    };

    return func;
}



/**
 * \brief diagrammatic symmetry relations for the screened interaction w
 */
template<typename Model> operation cconj_timerev_susc_sc( idx_susc_t<Model>& idx )
{
   bool sign= false;
#ifndef STATIC_CALCULATION
   idx(IW::W) *= -1; 
#endif
   return operation( sign, true ); 
}

template<typename Model> operation cconj_timerev_susc_d( idx_susc_t<Model>& idx )
{
   bool sign= false;
#ifndef STATIC_CALCULATION
   idx(IW::W) *= -1; 
#endif
   return operation( sign, true ); 
}

template<typename Model> operation cconj_timerev_susc_m( idx_susc_t<Model>& idx )
{
   bool sign= false;
#ifndef STATIC_CALCULATION
   idx(IW::W) *= -1;
#endif
   return operation( sign, true ); 
}

template<typename Model> operation timerev_susc_sc( idx_susc_t<Model>& idx ){
   bool sign= false;
   std::swap(idx(ISUSC::n_in),idx(ISUSC::n_out));
   return operation( sign, false ); 
}

template<typename Model> operation hmirror_timerev_susc_d( idx_susc_t<Model>& idx )
{
   bool sign= false;
   std::swap(idx(ISUSC::n_in),idx(ISUSC::n_out));
   return operation( sign, false ); 
}
template<typename Model> operation timerev_susc_m( idx_susc_t<Model>& idx ) {
   bool sign= false;
   std::swap(idx(ISUSC::n_in),idx(ISUSC::n_out));
   return operation( sign, false ); 
}

// instantiate class for the particular models
#define Instantiate(MODEL) template class SymmetriesfRGSBEPostprocessing<MODEL>;
WITH_MODELS(Instantiate)


