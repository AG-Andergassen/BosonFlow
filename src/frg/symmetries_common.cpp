#include <frg/symmetries_common.h>

#include <models/concrete_available_models.h>
#include <algorithm>
#include <frequencies/schemes.h>
#include <frg/interpolators_common.h>

template <typename Model>
void SymmetriesfRGCommon<Model>::Init()
{
    if (SymmetriesfRGCommon<Model>::AreAlgebraicSymmetriesIncluded() || SymmetriesfRGCommon<Model>::AreLatticeSymmetriesIncluded())
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
void SymmetriesfRGCommon<Model>::CalculateIndexEquivalenceClasses()
{
    // todo: if had been allocated before, should delete first..
    
    auto sampling_freqs_bubble = {
	std::make_pair(int(IBUBMAT::W), FrequencyDependenceScheme<Model>::BubbleTU_BosonicGrid().sampling_freqs), 
	std::make_pair(int(IBUBMAT::w), FrequencyDependenceScheme<Model>::BubbleTU_FermionicGrid().sampling_freqs)
    };
    

    std::cout << " Constructing pp-bubble symmetries... " << std::endl; 
    IdxEquivClasses_bubblemat_pp_Ptr() = new symmetry_grp_t<MatPatch,8>(gf_bubble_mat_t<Model>(), GroupActions_bubblemat(), sampling_freqs_bubble);

    std::cout << " Constructing ph-bubble symmetries... " << std::endl;
    IdxEquivClasses_bubblemat_ph_Ptr() = new symmetry_grp_t<MatPatch,8>(gf_bubble_mat_t<Model>(), GroupActions_bubblemat(), sampling_freqs_bubble);
    

    auto sampling_freqs_sig = {std::make_pair(int(I1P::w), FrequencyDependenceScheme<Model>::Sig_FermionicGrid().sampling_freqs), };
    std::cout << " Constructing self-energy (sig) symmetries... " << std::endl;
    IdxEquivClasses_sig_Ptr() = new symmetry_grp_t<dcomplex,4>(gf_1p_t<Model>(), GroupActions_sig(), sampling_freqs_sig);

    std::cout << " Constructing self-energy (sig_kmat) symmetries... " << std::endl;
    IdxEquivClasses_sigkmat_Ptr() = new symmetry_grp_t<MatPatch,3>(gf_Sig_kMat_t<Model>(), GroupActions_sigkmat(), sampling_freqs_sig);

    // todo: add Gvec's symmetries

    auto sampling_freqs_G = {std::make_pair(int(I1PMATREAL::w), FrequencyDependenceScheme<Model>::G_FermionicGrid().sampling_freqs), };
    std::cout << " Constructing real Green's function symmetries... " << std::endl;
    IdxEquivClasses_Gvec_real_Ptr() = new symmetry_grp_t<MatReal,3>(gf_1p_mat_real_t<Model>(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true),GroupActions_Gvec_real(), sampling_freqs_G );

    std::cout << " Constructing projection matrix symmetries... " << std::endl;
    IdxEquivClasses_projection_matrix_Ptr() = new symmetry_grp_t<dcomplex, 6>(gf_proj_matrix_t<Model>(),GroupActions_projection_matrix());
}

template <typename Model>
void SymmetriesfRGCommon<Model>::IncludeAlgebraicSymmetries()
{
    //todo: add check so that the following bubble symmetry is added only if form factors in momentum space are real!
    auto swap_bubble = []( idx_bubble_mat_t<Model>& idx ) 
    {
	std::swap(idx(IBUBMAT::m),idx(IBUBMAT::n)); 
	return operation( false, false ); 
    };
    GroupActions_bubblemat().push_back(swap_bubble);

    auto cconj_bubble = []( idx_bubble_mat_t<Model>& idx ){
#ifndef STATIC_CALCULATION
	idx(IBUBMAT::W) *= -1; 
	idx(IBUBMAT::w) = -idx(IBUBMAT::w)-1-abs(idx(IBUBMAT::W)%2);
#endif
	std::swap(idx(IBUBMAT::m),idx(IBUBMAT::n));
	std::swap(idx(IBUBMAT::s1),idx(IBUBMAT::s1p));
	std::swap(idx(IBUBMAT::s2),idx(IBUBMAT::s2p));
	return operation( false, true ); 
    };

    GroupActions_bubblemat().push_back(cconj_bubble);

    ///< Change sign of bosonic frequency + compl conj
    auto cconj_sig = [](idx_1p_t<Model>& idx ){
#ifndef STATIC_CALCULATION
	idx(I1P::w) = -idx(I1P::w) -1;
#endif
	std::swap(idx(I1P::s_in),idx(I1P::s_out));
	return operation( false, true ); 
    };

    // sig
    GroupActions_sig().push_back(cconj_sig);
    
    // --- Sig k_Mat (different indexing)
    ///< Change sign of bosonic frequency + compl conj
    auto cconj_sigkmat = []( idx_Sig_kMat_t<Model>& idx ){
#ifndef STATIC_CALCULATION
	idx(ISIGKMAT::w) = -idx(ISIGKMAT::w) -1;
#endif
	std::swap(idx(ISIGKMAT::s_in),idx(ISIGKMAT::s_out));
	return operation( false, true ); 
    };

    // sigkmat
    GroupActions_sigkmat().push_back(cconj_sigkmat);

    SymmetriesfRGCommon<Model>::AreAlgebraicSymmetriesIncluded() = true;
}

template <typename Model>
void SymmetriesfRGCommon<Model>::IncludeLatticeSymmetries()
{
    for (unsigned g_idx = 0; g_idx < Model::RealLattice().point_group_ptr->m_group_size; g_idx ++ ){
	// sig
	GroupActions_sig().push_back(SymmetriesfRGCommon<Model>::MakeLatticeSymmetryMap_sig(g_idx));

	//projection_matrix
	GroupActions_projection_matrix().push_back(SymmetriesfRGCommon<Model>::MakeLatticeSymmetryMap_projection_matrix(g_idx));
    }

    SymmetriesfRGCommon<Model>::AreLatticeSymmetriesIncluded() = true;
}


template <typename Model>
void SymmetriesfRGCommon<Model>::CleanUp()
{
    delete IdxEquivClasses_bubblemat_pp_Ptr();

    delete IdxEquivClasses_bubblemat_ph_Ptr();

    delete IdxEquivClasses_sig_Ptr();

    delete IdxEquivClasses_sigkmat_Ptr();

    delete IdxEquivClasses_projection_matrix_Ptr();
}


// --- Sig
template <typename Model>
std::function<operation(idx_1p_t<Model>& )> SymmetriesfRGCommon<Model>::MakeLatticeSymmetryMap_sig(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_1p_t<Model>& sig_idx){

	sig_idx(I1P::k) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][sig_idx(I1P::k)];	
	return operation(false, false );
    };

    return func;
}


template <typename Model>
std::function<operation(idx_proj_matrix_t<Model>& )> SymmetriesfRGCommon<Model>::MakeLatticeSymmetryMap_projection_matrix(unsigned point_group_element_idx)
{
    unsigned g_idx = point_group_element_idx;

    auto func = [g_idx](idx_proj_matrix_t<Model>& proj_idx){
	
	proj_idx(PROJ_MATRIX::K_in) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][proj_idx(PROJ_MATRIX::K_in)];	

	proj_idx(PROJ_MATRIX::K_out) =  Model::MomentumGrid().point_group_in_momentum_idx_space[g_idx][proj_idx(PROJ_MATRIX::K_out)];	


	bool with_minus, with_minus2, with_minus3, with_minus4;
	std::tie(proj_idx(PROJ_MATRIX::m), with_minus) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][proj_idx(PROJ_MATRIX::m)];

	std::tie(proj_idx(PROJ_MATRIX::mp), with_minus2) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][proj_idx(PROJ_MATRIX::mp)];

	std::tie(proj_idx(PROJ_MATRIX::n), with_minus3) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][proj_idx(PROJ_MATRIX::n)];

	std::tie(proj_idx(PROJ_MATRIX::np), with_minus4) =  Model::FormFactors().point_group_in_form_factor_idx_space[g_idx][proj_idx(PROJ_MATRIX::np)];

	// exclusive OR so that it's only true if one of the two signs , but not both (since minus times a minus is a plus), is a minus
	bool overall_with_minus = with_minus ^ with_minus2 ^ with_minus3 ^ with_minus4;
	
	return operation(overall_with_minus, false );
    };

    return func;
}


// instantiate class for the particular models
#define Instantiate(MODEL) template class SymmetriesfRGCommon<MODEL>;
WITH_MODELS(Instantiate)


