#include <frg/sbe/symmetries.h>

#include <models/square_hubbard.h>
#include <algorithm>

template <typename Model>
void Symmetries1lfRGSBE<Model>::Init()
{
    SymmetriesfRGCommon<Model>::Init();

    if (Symmetries1lfRGSBE<Model>::AreAlgebraicSymmetriesIncluded() || Symmetries1lfRGSBE<Model>::AreLatticeSymmetriesIncluded())
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
void Symmetries1lfRGSBE<Model>::CalculateIndexEquivalenceClasses()
{
    // todo: if had been allocated before, should delete first..

    auto sampling_freqs_w = {
	std::make_pair(int(IW::W), FrequencyDependenceScheme<Model>::w_BosonicGrid().sampling_freqs)
    };

    std::cout << " Constructing screened-interaction (w) symmetries... " << std::endl;
    IdxEquivClasses_w_sc_Ptr() = new symmetry_grp_t<dcomplex,2>(gf_w_t<Model>(), GroupActions_w_sc(), sampling_freqs_w);
    IdxEquivClasses_w_d_Ptr() = new symmetry_grp_t<dcomplex,2>(gf_w_t<Model>(), GroupActions_w_d(), sampling_freqs_w);
    IdxEquivClasses_w_m_Ptr() = new symmetry_grp_t<dcomplex,2>(gf_w_t<Model>(), GroupActions_w_m(), sampling_freqs_w);

    auto sampling_freqs_lambda = {
	std::make_pair(int(ILAMBDA::W), FrequencyDependenceScheme<Model>::lambda_BosonicGrid().sampling_freqs), 
	std::make_pair(int(ILAMBDA::w), FrequencyDependenceScheme<Model>::lambda_FermionicGrid().sampling_freqs)
    };
    

    std::cout << " Constructing hedin vertex (lambda) symmetries... " << std::endl;
    IdxEquivClasses_lambda_sc_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_lambda_t<Model>(), GroupActions_lambda_sc(), sampling_freqs_lambda);
    IdxEquivClasses_lambda_d_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_lambda_t<Model>(), GroupActions_lambda_d(), sampling_freqs_lambda);
    IdxEquivClasses_lambda_m_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_lambda_t<Model>(), GroupActions_lambda_m(), sampling_freqs_lambda);

    auto sampling_freqs_M = {
	std::make_pair(int(IM::W), FrequencyDependenceScheme<Model>::M_BosonicGrid().sampling_freqs), 
	std::make_pair(int(IM::w), FrequencyDependenceScheme<Model>::M_FermionicGrid().sampling_freqs)
#ifndef BOSONISE_M_LAZY
,
	std::make_pair(int(IM::wp), FrequencyDependenceScheme<Model>::M_FermionicGrid().sampling_freqs)
#endif
    };

    std::cout << " Constructing SBE rest function (M) symmetries... " << std::endl;
    IdxEquivClasses_M_sc_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_M_t<Model>(), GroupActions_M_sc(), sampling_freqs_M);
    IdxEquivClasses_M_d_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_M_t<Model>(), GroupActions_M_d(), sampling_freqs_M);
    IdxEquivClasses_M_m_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_M_t<Model>(), GroupActions_M_m(), sampling_freqs_M);

    
    auto sampling_freqs_projected_nabla = {
	std::make_pair(int(INABLA::W), FrequencyDependenceScheme<Model>::Projected_nabla_BosonicGrid().sampling_freqs), 
	std::make_pair(int(INABLA::w), FrequencyDependenceScheme<Model>::Projected_nabla_FermionicGrid().sampling_freqs),
	std::make_pair(int(INABLA::wp), FrequencyDependenceScheme<Model>::Projected_nabla_FermionicGrid().sampling_freqs)
    };

    
    std::cout << " Constructing projected nabla symmetries... " << std::endl;
    IdxEquivClasses_projected_nabla_pp_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_nabla_t<Model>( FrequencyDependenceScheme<Model>::Projected_nabla_BosonicGrid().positive_freqs_count, FrequencyDependenceScheme<Model>::Projected_nabla_FermionicGrid().positive_freqs_count), GroupActions_M_sc(), sampling_freqs_projected_nabla);
    IdxEquivClasses_projected_nabla_ph_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_nabla_t<Model>(FrequencyDependenceScheme<Model>::Projected_nabla_BosonicGrid().positive_freqs_count, FrequencyDependenceScheme<Model>::Projected_nabla_FermionicGrid().positive_freqs_count), GroupActions_M_d(), sampling_freqs_projected_nabla);
    IdxEquivClasses_projected_nabla_xph_Ptr() = new symmetry_grp_t<dcomplex,6>(gf_nabla_t<Model>(FrequencyDependenceScheme<Model>::Projected_nabla_BosonicGrid().positive_freqs_count, FrequencyDependenceScheme<Model>::Projected_nabla_FermionicGrid().positive_freqs_count), GroupActions_M_m(), sampling_freqs_projected_nabla);
}

template <typename Model>
void Symmetries1lfRGSBE<Model>::IncludeAlgebraicSymmetries()
{
    // w
    GroupActions_w_sc().push_back(hmirror_w_sc<Model>);
    GroupActions_w_d().push_back(hmirror_w_d<Model>);
    GroupActions_w_m().push_back(hmirror_w_m<Model>);

    GroupActions_w_sc().push_back(cconj_w_sc<Model>);
    //GroupActions_w_d().push_back(cconj_w_ph<Model>);
    GroupActions_w_m().push_back(cconj_w_m<Model>);

    GroupActions_w_sc().push_back(timerev_w_sc<Model>);
    //GroupActions_w_d().push_back(timerev_w_ph<Model>);
    GroupActions_w_m().push_back(timerev_w_m<Model>);


    // lambda
    GroupActions_lambda_sc().push_back(cconj_timerev_lambda_sc<Model>);
    GroupActions_lambda_d().push_back(cconj_timerev_lambda_d<Model>);
    GroupActions_lambda_m().push_back(cconj_timerev_lambda_m<Model>);
    
    // M
#ifndef BOSONISE_M_LAZY
    GroupActions_M_sc().push_back(timerev_M_sc<Model>);
    GroupActions_M_d().push_back(hmirror_timerev_M_d<Model>);
    GroupActions_M_m().push_back(timerev_M_m<Model>);
#endif

    GroupActions_M_sc().push_back(cconj_timerev_M_sc<Model>);
    GroupActions_M_d().push_back(cconj_timerev_M_d<Model>);
    GroupActions_M_m().push_back(cconj_timerev_M_m<Model>);
    
    Symmetries1lfRGSBE<Model>::AreAlgebraicSymmetriesIncluded() = true;
}

template <typename Model>
void Symmetries1lfRGSBE<Model>::IncludeLatticeSymmetries()
{
    for (unsigned g_idx = 0; g_idx < Model::RealLattice().point_group_ptr->m_group_size; g_idx ++ ){
	
	// for the screened interaction w
	GroupActions_w_sc().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_w(g_idx));
	GroupActions_w_d().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_w(g_idx));
	GroupActions_w_m().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_w(g_idx));


	// for the hedin vertex lambda
	GroupActions_lambda_sc().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_lambda(g_idx));
	GroupActions_lambda_d().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_lambda(g_idx));
	GroupActions_lambda_m().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_lambda(g_idx));

	// for the SBE-rest function M
	GroupActions_M_sc().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_M(g_idx));
	GroupActions_M_d().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_M(g_idx));
	GroupActions_M_m().push_back(Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_M(g_idx));
    }

    Symmetries1lfRGSBE<Model>::AreLatticeSymmetriesIncluded() = true;
}


template <typename Model>
void Symmetries1lfRGSBE<Model>::CleanUp()
{
    SymmetriesfRGCommon<Model>::CleanUp();

    delete IdxEquivClasses_w_sc_Ptr();
    delete IdxEquivClasses_w_d_Ptr();
    delete IdxEquivClasses_w_m_Ptr();

    delete IdxEquivClasses_lambda_sc_Ptr();
    delete IdxEquivClasses_lambda_d_Ptr();
    delete IdxEquivClasses_lambda_m_Ptr();

    delete IdxEquivClasses_M_sc_Ptr();
    delete IdxEquivClasses_M_d_Ptr();
    delete IdxEquivClasses_M_m_Ptr();

    delete IdxEquivClasses_projected_nabla_pp_Ptr();
    delete IdxEquivClasses_projected_nabla_ph_Ptr();
    delete IdxEquivClasses_projected_nabla_xph_Ptr();
}


//        LATTICE SYMMETRIES


/**
 * \brief lattice symmetry for the screened interaction w
 *
 * same for all channels.
 */
template <typename Model>
std::function<operation(idx_w_t<Model>& )> Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_w(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_w_t<Model>& w_idx){

	w_idx(IW::K) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][w_idx(IW::K)];	
	
	return operation(false, false );
    };

    return func;
}

/**
 * \brief lattice symmetry for the Hedin vertex lambda
 *
 * same for all channels.
 */
template <typename Model>
std::function<operation(idx_lambda_t<Model>& )> Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_lambda(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_lambda_t<Model>& lambda_idx){

	lambda_idx(ILAMBDA::K) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][lambda_idx(ILAMBDA::K)];	

	bool with_minus;
	// note that the form-factors are based on the coarse lattice (e.g. bond form factors from the coarse lattice)
	std::tie(lambda_idx(ILAMBDA::m), with_minus) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][lambda_idx(ILAMBDA::m)];
	
	return operation(with_minus, false );
    };

    return func;
}


/**
 * \brief lattice symmetry for the U-irreducible but two-particle-reducible vertex M
 *
 * same for all channels.
 */

template <typename Model>
std::function<operation(idx_M_t<Model>& )> Symmetries1lfRGSBE<Model>::MakeLatticeSymmetryMap_M(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_M_t<Model>& M_idx){

	M_idx(IM::K) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][M_idx(IM::K)];	

	bool with_minus, with_minus2;
	// note that the form-factors are based on the coarse lattice (e.g. bond form factors from the coarse lattice)
	std::tie(M_idx(IM::m), with_minus) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][M_idx(IM::m)];

	std::tie(M_idx(IM::mp), with_minus2) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][M_idx(IM::mp)];

	// exclusive OR so that it's only true if one of the two signs , but not both (since minus times a minus is a plus), is a minus
	bool overall_with_minus = with_minus ^ with_minus2;

	
	return operation(overall_with_minus, false );
    };

    return func;
}


/**
 * \brief diagrammatic symmetry relations for the Hedin vertex lambda
 */

template<typename Model> operation cconj_timerev_lambda_sc( idx_lambda_t<Model>& idx )
{
   bool sign = false;
#ifndef STATIC_CALCULATION
   idx(ILAMBDA::W) *= -1; 
   idx(ILAMBDA::w) =-idx(ILAMBDA::w)-1-abs(idx(ILAMBDA::W)%2);
#endif
   return operation( sign, true ); // Quantum numbers don't change after cconj and time reversal
}

template<typename Model> operation cconj_timerev_lambda_d( idx_lambda_t<Model>& idx )
{
   bool sign = false;
#ifndef STATIC_CALCULATION
   idx(ILAMBDA::W) *= -1; 
   idx(ILAMBDA::w) =-idx(ILAMBDA::w)-1-abs(idx(ILAMBDA::W)%2);
#endif
   return operation( sign, true ); // Quantum numbers don't change after cconj and time reversal
}


template<typename Model> operation cconj_timerev_lambda_m( idx_lambda_t<Model>& idx )
{
   bool sign = false;
#ifndef STATIC_CALCULATION
   idx(ILAMBDA::W) *= -1; 
   idx(ILAMBDA::w) =-idx(ILAMBDA::w)-1-abs(idx(ILAMBDA::W)%2);
#endif
   return operation( sign, true ); 
}

/**
 * \brief diagrammatic symmetry relations for the screened interaction w
 */
template<typename Model> operation hmirror_w_sc( idx_w_t<Model>& idx ) //< becomes important when additional quantum numbers are involved
{
   bool sign= false;
   return operation( sign, false ); 
}
template<typename Model> operation hmirror_w_d( idx_w_t<Model>& idx )
{
   bool sign= false;
#ifndef STATIC_CALCULATION
   idx(IW::W) *= -1;
#endif 
   idx(IW::K) =Model::GetNegativeCoarseMomentumIdx(idx(IW::K));
   return operation( sign, false ); 
}

template<typename Model> operation hmirror_w_m( idx_w_t<Model>& idx )
{
   bool sign= false;
#ifndef STATIC_CALCULATION
   idx(IW::W) *= -1;
#endif
   idx(IW::K) = Model::GetNegativeCoarseMomentumIdx(idx(IW::K));
   return operation( sign, false ); 
}

template<typename Model> operation cconj_w_sc( idx_w_t<Model>& idx )
{
   bool sign= false;
#ifndef STATIC_CALCULATION
   idx(IW::W) *= -1; 
#endif
   return operation( sign, true ); 
}

template<typename Model> operation cconj_w_ph( idx_w_t<Model>& idx )
{
   bool sign= false;
   idx(IW::K)  =  Model::GetNegativeCoarseMomentumIdx(idx(IW::K));
   return operation( sign, true ); 
}

template<typename Model> operation cconj_w_m( idx_w_t<Model>& idx )
{
   bool sign= false;
#ifndef STATIC_CALCULATION
   idx(IW::W) *= -1;
#endif
   return operation( sign, true ); 
}

template<typename Model> operation timerev_w_sc( idx_w_t<Model>& idx ) //< becomes important when additional quantum numbers are involved
{
   bool sign= false;
   return operation( sign, false ); 
}

template<typename Model> operation timerev_w_ph( idx_w_t<Model>& idx )
{
   bool sign= false;
#ifndef STATIC_CALCULATION
   idx(IW::W) *= -1;
#endif
   idx(IW::K) = Model::GetNegativeCoarseMomentumIdx(idx(IW::K));
   return operation( sign, false ); 
}
template<typename Model> operation timerev_w_m( idx_w_t<Model>& idx ) //< becomes important when additional quantum numbers are involved
{
   bool sign= false;
   return operation( sign, false ); 
}


/**
 * \brief diagrammatic symmetry relations for the U-irreducible but two-particle-reducible vertex M
 */
template<typename Model> operation cconj_timerev_M_sc( idx_M_t<Model>& idx )
{
   bool sign = false;
#ifndef STATIC_CALCULATION
   idx(IM::W) *= -1;
   // w = 0 -> v = first matsubara freq
   idx(IM::w) =-idx(IM::w)-1-(idx(IM::W)+10000)%2;
   idx(IM::wp) =-idx(IM::wp)-1-(idx(IM::W)+10000)%2;
#ifdef BOSONISE_M_LAZY   // when bosonising M, there's only a 0 and a -1 index
   idx(IM::wp) += (idx(IM::W)+10000)%2;
#endif
#endif
   return operation( sign, true );
}

template<typename Model> operation cconj_timerev_M_d( idx_M_t<Model>& idx )
{
   bool sign = false;
#ifndef STATIC_CALCULATION   
   idx(IM::W) *=  -1;
   idx(IM::w) =-idx(IM::w)-1-(idx(IM::W)+10000)%2;
   idx(IM::wp) =-idx(IM::wp)-1-(idx(IM::W)+10000)%2;
#ifdef BOSONISE_M_LAZY   
   idx(IM::wp) += (idx(IM::W)+10000)%2;
#endif
#endif
   return operation( sign, true );
}

template<typename Model> operation cconj_timerev_M_m( idx_M_t<Model>& idx )
{
   bool sign = false;
#ifndef STATIC_CALCULATION
   idx(IM::W) *= -1;
   idx(IM::w) =-idx(IM::w)-1-(idx(IM::W)+10000)%2;
   idx(IM::wp) =-idx(IM::wp)-1-(idx(IM::W)+10000)%2;
#ifdef BOSONISE_M_LAZY   
   idx(IM::wp) += (idx(IM::W)+10000)%2;
#endif
#endif
   return operation( sign, true );
}

template<typename Model> operation timerev_M_sc( idx_M_t<Model>& idx )
{
   bool sign = false;
#ifndef STATIC_CALCULATION
   std::swap(idx(IM::w),idx(IM::wp));
#endif
   std::swap(idx(IM::m),idx(IM::mp));
   return operation( sign, false );
}

template<typename Model> operation timerev_M_m( idx_M_t<Model>& idx )
{
   bool sign = false;
#ifndef STATIC_CALCULATION
   std::swap(idx(IM::w),idx(IM::wp));
#endif
   std::swap(idx(IM::m),idx(IM::mp));
   return operation( sign, false );
}

template<typename Model> operation hmirror_timerev_M_d( idx_M_t<Model>& idx )
{
   bool sign = false;
#ifndef STATIC_CALCULATION
   std::swap(idx(IM::w),idx(IM::wp));
#endif
   std::swap(idx(IM::m),idx(IM::mp));
   return operation( sign, false );
}



// instantiate class for the particular models
#define Instantiate(MODEL) template class Symmetries1lfRGSBE<MODEL>;
WITH_MODELS(Instantiate)


