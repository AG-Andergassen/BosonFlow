#include <mymath.h>
#include <cmath>
#include <complex>
#include <frg/sbe/state.h>
#include <frg/sbe/postprocessing/symmetries.h>

using namespace std; 
using boost::bind;

template <typename Model, typename state_t>
rhs_postproc_sbe_t<Model, state_t>::rhs_postproc_sbe_t()
{ 
}


template <typename Model, typename state_t>
rhs_postproc_sbe_t<Model, state_t>::~rhs_postproc_sbe_t()
{
}


template <typename Model, typename state_t> 
void rhs_postproc_sbe_t<Model, state_t>::operator() ( const state_t& state, state_postproc_t<Model>& postproc_state, const double t )
{

    // Precalculate Green function and single-scale propagator on frequency grid
    gf_1p_mat_t<Model> Gvec( FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true ); 

    rhs_base_t::compute_Gvec(&Gvec, state, t);

    std::cout << " Postprocessing calculation " << std::endl;


    gf_1p_mat_real_t<Model> Gvec_realspace(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true), Svec_realspace(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, true);   
    gf_bubble_mat_t<Model> bubble_pp, bubble_ph;


#pragma omp parallel default(shared) 
    {
	rhs_base_t::compute_Gvec_real(&Gvec_realspace, Gvec);
#pragma omp barrier //barrier to make sure Gvec_realspace is present
	rhs_base_t::compute_GS_bubbles(&bubble_pp, &bubble_ph, Gvec_realspace, Gvec_realspace);
	//rhs_base_t::compute_GG_bubbles(&bubble_pp, &bubble_ph, Gvec_realspace);
    }

    // TODO: after adding factor of 1/2 in the original definition, remove this
    bubble_pp.init([&bubble_pp](const idx_bubble_mat_t<Model> &idx){return 0.5 * bubble_pp(idx);});
    bubble_ph.init([&bubble_ph](const idx_bubble_mat_t<Model> &idx){return 0.5 * bubble_ph(idx);});

    cout << " ... Susceptibilities " << endl;
   
    cout << " ... SC " << endl;
   
#pragma omp parallel default(shared) 
    {
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_susc_sc(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc( idx, state, bubble_pp, t ); } ); 
    
#ifdef SPLIT_SUSC_CONTRIBUTIONS
#pragma omp single
	std::cout << " ... Split contributions to SC" << endl;
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_suscvert_sc( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscbubble_sc(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_suscbubble_sc( idx, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_nabla_d(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_nabla_d( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_M_d(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_M_d( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_d_double_counting(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_d_double_counting( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_nabla_m(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_nabla_m( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_M_m(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_M_m( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_m_double_counting(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_m_double_counting( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_nabla_sc(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_nabla_sc( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_M_sc(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_M_sc( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_sc_double_counting(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_sc_double_counting( idx, state, bubble_pp, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( postproc_state.gf_suscvert_sc_contribution_from_bare(),    [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_bare( idx, state, bubble_pp, t ); } );
#endif

#pragma omp single
	cout << " ... Density " << endl;

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_susc_d(),      [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d( idx, state, bubble_ph, t ); } );

#ifdef SPLIT_SUSC_CONTRIBUTIONS
#pragma omp single
	std::cout << " ... Split contributions to D" << endl;
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscbubble_d(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_suscbubble_d( idx, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_suscvert_d( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_nabla_d(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_nabla_d( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_M_d(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_M_d( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_d_double_counting(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_d_double_counting( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_nabla_m(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_nabla_m( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_M_m(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_M_m( idx, state, bubble_ph, t ); } );

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_m_double_counting(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_m_double_counting( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_nabla_sc(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_nabla_sc( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_M_sc(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_M_sc( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_sc_double_counting(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_sc_double_counting( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_d_Ptr()->init_batched( postproc_state.gf_suscvert_d_contribution_from_bare(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_d_contribution_from_bare( idx, state, bubble_ph, t ); } );
#endif 

#pragma omp single
	cout << " ... Magnetic " << endl;

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_susc_m(),      [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m( idx, state, bubble_ph, t ); } ); 

#ifdef SPLIT_SUSC_CONTRIBUTIONS
#pragma omp single
	std::cout << " ... Split contributions to M" << endl;
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_suscvert_m( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscbubble_m(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_suscbubble_m( idx, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_nabla_d(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_nabla_d( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_M_d(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_M_d( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_d_double_counting(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_d_double_counting( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_nabla_m(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_nabla_m( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_M_m(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_M_m( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_m_double_counting(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_m_double_counting( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_nabla_sc(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_nabla_sc( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_M_sc(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_M_sc( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_sc_double_counting(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_sc_double_counting( idx, state, bubble_ph, t ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_m_Ptr()->init_batched( postproc_state.gf_suscvert_m_contribution_from_bare(),    [&state, &bubble_ph, t]( const idx_susc_t<Model>& idx ){ return eval_susc_m_contribution_from_bare( idx, state, bubble_ph, t ); } );

#pragma omp barrier

#pragma omp single
	cout << " ... Polarisations " << endl;
    
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( postproc_state.gf_polarisation_sc_contribution_from_1(),  [&state, &bubble_pp, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_sc_contribution_from_bubble( idx, state, bubble_pp, t); } );
	
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( postproc_state.gf_polarisation_sc_contribution_from_varphi(),  [&state, &bubble_pp, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_sc_contribution_from_varphi( idx, state, bubble_pp, t); } );

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( postproc_state.gf_polarisation_sc_contribution_from_nabla_sc(),  [&state, &bubble_pp, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_sc_contribution_from_nabla_sc( idx, state, bubble_pp, t); } ); // it's zero, but anyways we include it to retain the "symmetry" of the code
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( postproc_state.gf_polarisation_sc_contribution_from_nabla_m(),  [&state, &bubble_pp, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_pp, t , postproc_state.gf_suscvert_sc_contribution_from_nabla_m() ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( postproc_state.gf_polarisation_sc_contribution_from_nabla_d(),  [&state, &bubble_pp, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_pp, t , postproc_state.gf_suscvert_sc_contribution_from_nabla_d() ); } );


	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_m_Ptr()->init_batched( postproc_state.gf_polarisation_m_contribution_from_1(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_m_or_d_contribution_from_bubble( idx, state, bubble_ph, t); } );
	
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_m_Ptr()->init_batched( postproc_state.gf_polarisation_m_contribution_from_varphi(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_m_contribution_from_varphi( idx, state, bubble_ph, t); } );
	
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_m_Ptr()->init_batched( postproc_state.gf_polarisation_m_contribution_from_nabla_m(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_m_contribution_from_nabla_m( idx, state, bubble_ph, t); } ); 
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_m_Ptr()->init_batched( postproc_state.gf_polarisation_m_contribution_from_nabla_sc(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_ph, t , postproc_state.gf_suscvert_m_contribution_from_nabla_sc() ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_m_Ptr()->init_batched( postproc_state.gf_polarisation_m_contribution_from_nabla_d(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_ph, t , postproc_state.gf_suscvert_m_contribution_from_nabla_d() ); } );

	
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( postproc_state.gf_polarisation_d_contribution_from_1(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_m_or_d_contribution_from_bubble( idx, state, bubble_ph, t); } );
	
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( postproc_state.gf_polarisation_d_contribution_from_varphi(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_d_contribution_from_varphi( idx, state, bubble_ph, t); } );
	
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( postproc_state.gf_polarisation_d_contribution_from_nabla_d(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_d_contribution_from_nabla_d( idx, state, bubble_ph, t); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( postproc_state.gf_polarisation_d_contribution_from_nabla_sc(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_ph, t , postproc_state.gf_suscvert_d_contribution_from_nabla_sc() ); } );
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( postproc_state.gf_polarisation_d_contribution_from_nabla_m(),  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_ph, t , postproc_state.gf_suscvert_d_contribution_from_nabla_m() ); } );

#pragma omp barrier
	auto suscvert_sc_contribution_from_Virr = postproc_state.gf_suscvert_sc_contribution_from_M_m()
	    +postproc_state.gf_suscvert_sc_contribution_from_m_double_counting()
	    +postproc_state.gf_suscvert_sc_contribution_from_M_sc()
	    +postproc_state.gf_suscvert_sc_contribution_from_sc_double_counting()
	    +postproc_state.gf_suscvert_sc_contribution_from_M_d()
	    +postproc_state.gf_suscvert_sc_contribution_from_d_double_counting()
	    +postproc_state.gf_suscvert_sc_contribution_from_bare();

	auto suscvert_m_contribution_from_Virr = postproc_state.gf_suscvert_m_contribution_from_M_m()
	    +postproc_state.gf_suscvert_m_contribution_from_m_double_counting()
	    +postproc_state.gf_suscvert_m_contribution_from_M_sc()
	    +postproc_state.gf_suscvert_m_contribution_from_sc_double_counting()
	    +postproc_state.gf_suscvert_m_contribution_from_M_d()
	    +postproc_state.gf_suscvert_m_contribution_from_d_double_counting()
	    +postproc_state.gf_suscvert_m_contribution_from_bare();

	auto suscvert_d_contribution_from_Virr = postproc_state.gf_suscvert_d_contribution_from_M_m()
	    +postproc_state.gf_suscvert_d_contribution_from_m_double_counting()
	    +postproc_state.gf_suscvert_d_contribution_from_M_sc()
	    +postproc_state.gf_suscvert_d_contribution_from_sc_double_counting()
	    +postproc_state.gf_suscvert_d_contribution_from_M_d()
	    +postproc_state.gf_suscvert_d_contribution_from_d_double_counting()
	    +postproc_state.gf_suscvert_d_contribution_from_bare();

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( postproc_state.gf_polarisation_sc_contribution_from_Virr(),  [&state, &bubble_pp, &suscvert_sc_contribution_from_Virr,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_pp, t , suscvert_sc_contribution_from_Virr); } );

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( postproc_state.gf_polarisation_d_contribution_from_Virr(),  [&state, &bubble_ph, &suscvert_d_contribution_from_Virr,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_ph, t , suscvert_d_contribution_from_Virr); } );

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_m_Ptr()->init_batched( postproc_state.gf_polarisation_m_contribution_from_Virr(),  [&state, &bubble_ph, &suscvert_m_contribution_from_Virr,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_ph, t , suscvert_m_contribution_from_Virr); } );

    }

    std::cout <<"Equality d: "<< norm( 
				      -2.0 * postproc_state.gf_suscvert_d_contribution_from_bare() + 
				      postproc_state.gf_suscvert_d_contribution_from_nabla_d() + 
				      postproc_state.gf_suscvert_d_contribution_from_M_d() +
				      //      postproc_state.gf_suscvert_d_contribution_from_d_double_counting() + 
				      postproc_state.gf_suscvert_d_contribution_from_nabla_sc() +
				      postproc_state.gf_suscvert_d_contribution_from_M_sc() +
				      //   postproc_state.gf_suscvert_d_contribution_from_sc_double_counting() + 
				      postproc_state.gf_suscvert_d_contribution_from_nabla_m() +
				      postproc_state.gf_suscvert_d_contribution_from_M_m() + 
				      //  postproc_state.gf_suscvert_d_contribution_from_m_double_counting()+ 
				      postproc_state.gf_suscbubble_d()  
				      - postproc_state.gf_susc_d()
				       ) << std::endl;

    std::cout <<"Equality sc: "<< norm( 
				       -2.0*postproc_state.gf_suscvert_sc_contribution_from_bare() + 
				       postproc_state.gf_suscvert_sc_contribution_from_nabla_d() +
				       postproc_state.gf_suscvert_sc_contribution_from_M_d() + 
				       //		       postproc_state.gf_suscvert_sc_contribution_from_d_double_counting() + 
				       postproc_state.gf_suscvert_sc_contribution_from_nabla_sc() + 
				       postproc_state.gf_suscvert_sc_contribution_from_M_sc() + 
				       // postproc_state.gf_suscvert_sc_contribution_from_sc_double_counting() + 
				       postproc_state.gf_suscvert_sc_contribution_from_nabla_m() + 
				       postproc_state.gf_suscvert_sc_contribution_from_M_m() + 
				       //postproc_state.gf_suscvert_sc_contribution_from_m_double_counting() + 
				       postproc_state.gf_suscbubble_sc()
				       - postproc_state.gf_susc_sc()) << std::endl;

    std::cout <<"Equality m: "<< norm(
				      -2.0*postproc_state.gf_suscvert_m_contribution_from_bare() + 
				      postproc_state.gf_suscvert_m_contribution_from_nabla_d() + 
				      postproc_state.gf_suscvert_m_contribution_from_M_d() + 
				      //postproc_state.gf_suscvert_m_contribution_from_d_double_counting() + 
				      postproc_state.gf_suscvert_m_contribution_from_nabla_sc() + 
				      postproc_state.gf_suscvert_m_contribution_from_M_sc() + 
				      //postproc_state.gf_suscvert_m_contribution_from_sc_double_counting() + 
				      postproc_state.gf_suscvert_m_contribution_from_nabla_m() + 
				      postproc_state.gf_suscvert_m_contribution_from_M_m() + 
				      //postproc_state.gf_suscvert_m_contribution_from_m_double_counting() + 
				      postproc_state.gf_suscbubble_m()
				      -  postproc_state.gf_susc_m()) << std::endl;

#endif

    cout << " ... Hedin vertices " << endl;
   


#pragma omp parallel default(shared) 
    {
	// lambdas
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_varphi(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_varphi( idx, state, bubble_pp, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_varphi(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_varphi( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_varphi(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_varphi( idx, state, bubble_ph, t ); });

    }
#ifdef SPLIT_LAMBDA_CONTRIBUTIONS    
	// todo: implement symmetries and parallelisation

    #pragma omp parallel default(shared) 
    {

#pragma omp single
	cout << " ... SC " << endl;    
	// it's zero, but we add it anyway
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_bare(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_bare( idx, state, bubble_pp, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_nabla_sc(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_nabla_sc( idx, state, bubble_pp, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_nabla_d(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_nabla_d( idx, state, bubble_pp, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_nabla_m(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_nabla_m( idx, state, bubble_pp, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_M_sc(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_M_sc( idx, state, bubble_pp, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_M_d(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_M_d( idx, state, bubble_pp, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_M_m(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_M_m( idx, state, bubble_pp, t ); });

	
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_sc_double_counting(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_sc_double_counting( idx, state, bubble_pp, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_d_double_counting(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_d_double_counting( idx, state, bubble_pp, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_sc_Ptr()->init_batched(postproc_state.gf_lambda_sc_contribution_from_m_double_counting(), [&state, &bubble_pp, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_m_double_counting( idx, state, bubble_pp, t ); });

#pragma omp barrier


#pragma omp single
	cout << " ... D " << endl;    
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_bare(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_bare( idx, state, bubble_ph, t ); });
    
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_nabla_d(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_nabla_d( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_nabla_m(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_nabla_m( idx, state, bubble_ph, t ); });
    
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_nabla_sc(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_nabla_sc( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_M_d(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_M_d( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_M_m(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_M_m( idx, state, bubble_ph, t ); });
    
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_M_sc(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_M_sc( idx, state, bubble_ph, t ); });

	
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_d_double_counting(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_d_double_counting( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_sc_double_counting(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_sc_double_counting( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_d_Ptr()->init_batched(postproc_state.gf_lambda_d_contribution_from_m_double_counting(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_m_double_counting( idx, state, bubble_ph, t ); });
    
#pragma omp barrier


#pragma omp single
	cout << " ... M " << endl;    
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_bare(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_bare( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_m_double_counting(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_m_double_counting( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_sc_double_counting(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_sc_double_counting( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_d_double_counting(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_d_double_counting( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_nabla_m(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_nabla_m( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_nabla_sc(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_nabla_sc( idx, state, bubble_ph, t ); });
    
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_nabla_d(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_nabla_d( idx, state, bubble_ph, t ); });

	
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_M_m(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_M_m( idx, state, bubble_ph, t ); });

	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_M_sc(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_M_sc( idx, state, bubble_ph, t ); });
    
	SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_static_lambda_m_Ptr()->init_batched(postproc_state.gf_lambda_m_contribution_from_M_d(), [&state, &bubble_ph, t]( const idx_static_lambda_t<Model>& idx ){ return eval_lambda_m_contribution_from_M_d( idx, state, bubble_ph, t ); });
	
    }
    



    std::cout <<"Equality lambda d: "<< norm(  
				      postproc_state.gf_lambda_d_contribution_from_bare() +
				      postproc_state.gf_lambda_d_contribution_from_nabla_d() + 
				      postproc_state.gf_lambda_d_contribution_from_nabla_sc() + 
				      postproc_state.gf_lambda_d_contribution_from_nabla_m() +  
				      postproc_state.gf_lambda_d_contribution_from_sc_double_counting() +
				      postproc_state.gf_lambda_d_contribution_from_d_double_counting() +
				      postproc_state.gf_lambda_d_contribution_from_m_double_counting() +
				      - postproc_state.gf_lambda_d_contribution_from_varphi()
				       ) << std::endl;

    std::cout <<"Equality lambda sc: "<< norm(  
				      postproc_state.gf_lambda_sc_contribution_from_bare() +
				      postproc_state.gf_lambda_sc_contribution_from_nabla_sc() +
				      postproc_state.gf_lambda_sc_contribution_from_nabla_d() + 
				      postproc_state.gf_lambda_sc_contribution_from_nabla_m() + 
				      postproc_state.gf_lambda_sc_contribution_from_sc_double_counting() +
				      postproc_state.gf_lambda_sc_contribution_from_d_double_counting() +
				      postproc_state.gf_lambda_sc_contribution_from_m_double_counting() +
				      - postproc_state.gf_lambda_sc_contribution_from_varphi()
				       ) << std::endl;

    std::cout <<"Equality lambda m: "<< norm(  
					      postproc_state.gf_lambda_m_contribution_from_bare() + 
					      postproc_state.gf_lambda_m_contribution_from_nabla_sc() + 
					      postproc_state.gf_lambda_m_contribution_from_nabla_d() +
					      postproc_state.gf_lambda_m_contribution_from_nabla_m() + 
					      postproc_state.gf_lambda_m_contribution_from_sc_double_counting() +
					      postproc_state.gf_lambda_m_contribution_from_d_double_counting() +
					      postproc_state.gf_lambda_m_contribution_from_m_double_counting()
				      - postproc_state.gf_lambda_m_contribution_from_varphi()
				       ) << std::endl;



#endif // SPLIT_LAMBDA_CONTRIBUTIONS    


    cout << " ... Self-energies integrands (SDE) " << endl;

    
    postproc_state.gf_Sig_SDE_integrand_sc().init([&state, &Gvec,  t]( const idx_Sig_SDE_integrand_t<Model>& idx ){ return eval_Sig_SDE_integrand_sc( idx, state, Gvec, t ); } );
    postproc_state.gf_Sig_SDE_integrand_d().init([&state, &Gvec,  t]( const idx_Sig_SDE_integrand_t<Model>& idx ){ return eval_Sig_SDE_integrand_d( idx, state, Gvec, t ); } );
    postproc_state.gf_Sig_SDE_integrand_m().init([&state, &Gvec,  t]( const idx_Sig_SDE_integrand_t<Model>& idx ){ return eval_Sig_SDE_integrand_m( idx, state, Gvec, t ); } );
    
}



// contributions to d
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertex_d( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
    return state_vec.vertex_d( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_d( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
    return state_vec.B_irreducible_vertex_d( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_bare( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
  return (Model::vertex_local_part_bare(0, 0, n_in, 0, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
  //return (Model::B_d(W, K, w_in, n_in, w_out, n_out) + Model::F_d(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_Virr( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
    return  vertx_d_contribution_from_bare( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_d_contribution_from_M_d( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_d_contribution_from_M_m( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_d_contribution_from_M_sc( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_d_contribution_from_d_double_counting( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_d_contribution_from_m_double_counting( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_d_contribution_from_sc_double_counting( state_vec, W, w_in, w_out, K, n_in, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_nabla_d( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return state_vec.nabla_d( W, K, w_in, n_in, w_out, n_out, t )- 0.5*state_vec.nabla_d_ph_to_xph(W, K, w_in, n_in, w_out, n_out, t);
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_d_contribution_from_nabla_d( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return - 0.5*state_vec.nabla_d_ph_to_xph(W, K, w_in, n_in, w_out, n_out, t);
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_M_d( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
#if !defined(SBEa_APPROXIMATION) 
    return state_vec.M_d( W, K, w_in, n_in, w_out, n_out )  - 0.5 * state_vec.M_d_ph_to_xph(W, K, w_in, n_in, w_out, n_out);
#endif
    return 0.0;
	
}



template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_d_double_counting( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return (- Model::B_d(W, K, w_in, n_in, w_out, n_out) + 1.5*Model::B_m(W, K, w_in, n_in, w_out, n_out) + 0.5*Model::F_d(W, K, w_in, n_in, w_out, n_out) - 1.5*Model::F_m(W, K, w_in, n_in, w_out, n_out) + 2.0*Model::vertex_local_part_bare(W, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_nabla_m( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return -2.0*state_vec.nabla_m_xph_to_ph( W, K, w_in, n_in, w_out, n_out, t ) + 0.5*state_vec.nabla_m_ph_to_xph( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_M_m( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
#if !defined(SBEa_APPROXIMATION)
    return  -2.0*state_vec.M_m_xph_to_ph(W, K, w_in, n_in, w_out, n_out) +0.5*state_vec.M_m_ph_to_xph(W, K, w_in, n_in, w_out, n_out); 
#endif
    return 0.0;
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_m_double_counting( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return (-1.5*(Model::F_d(W, K, w_in, n_in, w_out, n_out) - Model::F_m(W, K, w_in, n_in, w_out, n_out)) - 3.0*Model::vertex_local_part_bare(W, w_in, n_in, w_out, n_out) - 1.5*Model::B_m(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_nabla_sc( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return 2.0*state_vec.nabla_sc_pp_to_ph( W, K, w_in, n_in, w_out, n_out, t ) - state_vec.nabla_sc_pp_to_xph( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_M_sc( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) { 
#if !defined(SBEa_APPROXIMATION)
    return 2.0*state_vec.M_sc_pp_to_ph(W, K, w_in, n_in, w_out, n_out) - state_vec.M_sc_pp_to_xph(W, K, w_in, n_in, w_out, n_out);
#endif
    return 0.0;
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_d_contribution_from_sc_double_counting( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return (- (Model::F_d(W, K, w_in, n_in, w_out, n_out) - Model::F_m(W, K, w_in, n_in, w_out, n_out)) - 2.0*Model::vertex_local_part_bare(W, w_in, n_in, w_out, n_out) - Model::B_m(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}


// Evaluate contributions to superconducting susceptibility by channel
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertex_sc( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
    return state_vec.vertex_sc( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_sc( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
    return state_vec.B_irreducible_vertex_sc( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_bare( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
  return (Model::vertex_local_part_bare(0, 0, n_in, 0, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
  //return (Model::B_sc(W, K, w_in, n_in, w_out, n_out) + Model::F_sc(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_Virr( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
    return  vertx_sc_contribution_from_bare( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_sc_contribution_from_M_d( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_sc_contribution_from_M_m( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_sc_contribution_from_M_sc( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_sc_contribution_from_d_double_counting( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_sc_contribution_from_m_double_counting( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_sc_contribution_from_sc_double_counting( state_vec, W, w_in, w_out, K, n_in, n_out, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_nabla_d( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return 0.5 * state_vec.nabla_d_ph_to_pp( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_M_d( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t) {
#if !defined(SBEa_APPROXIMATION)
    return 0.5*state_vec.M_d_ph_to_pp( W, K, w_in, n_in, w_out, n_out );
#endif
    return 0.0;
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_d_double_counting( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return (- Model::F_sc(W, K, w_in, n_in, w_out, n_out) - Model::vertex_local_part_bare(W, w_in, n_in, w_out, n_out) + 0.5*Model::B_sc(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_nabla_m( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return - state_vec.nabla_m_xph_to_pp( W, K, w_in, n_in, w_out, n_out, t ) -0.5 * state_vec.nabla_m_ph_to_pp( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_M_m( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
#if !defined(SBEa_APPROXIMATION)
    return - state_vec.M_m_xph_to_pp( W, K, w_in, n_in, w_out, n_out ) -0.5 * state_vec.M_m_ph_to_pp( W, K, w_in, n_in, w_out, n_out );
#endif
    return 0.0;
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_m_double_counting( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return (- 1.5*Model::B_sc(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}



template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_nabla_sc( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return state_vec.nabla_sc( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_sc_contribution_from_nabla_sc( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return 0.0;
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_M_sc( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
#if !defined(SBEa_APPROXIMATION)
    return state_vec.M_sc( W, K, w_in, n_in, w_out, n_out );
#endif
    return 0.0;
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_sc_contribution_from_sc_double_counting( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return (- Model::B_sc(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}


// m contributions
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertex_m( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
    return state_vec.vertex_m( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_m( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
    return state_vec.B_irreducible_vertex_m( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_bare( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
  return -(Model::vertex_local_part_bare(0, 0, n_in, 0, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
  //return  (Model::B_m(W, K, w_in, n_in, w_out, n_out) + Model::F_m(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_Virr( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t )
{
    return  vertx_m_contribution_from_bare( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_m_contribution_from_M_d( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_m_contribution_from_M_m( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_m_contribution_from_M_sc( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_m_contribution_from_d_double_counting( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_m_contribution_from_m_double_counting( state_vec, W, w_in, w_out, K, n_in, n_out, t )
	+ vertx_m_contribution_from_sc_double_counting( state_vec, W, w_in, w_out, K, n_in, n_out, t );
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_nabla_d( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return - 0.5 * state_vec.nabla_d_ph_to_xph( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_M_d( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
#if !defined(SBEa_APPROXIMATION) 
    return - 0.5*state_vec.M_d_ph_to_xph( W, K, w_in, n_in, w_out, n_out );
#endif
    return 0.0;
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_d_double_counting( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return (- Model::F_m(W, K, w_in, n_in, w_out, n_out) + Model::vertex_local_part_bare(W, w_in, n_in, w_out, n_out) + 0.5 * Model::B_m(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_nabla_m( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return 0.5 * state_vec.nabla_m_ph_to_xph( W, K, w_in, n_in, w_out, n_out, t ) + state_vec.nabla_m( W, K, w_in, n_in, w_out, n_out, t );
}

// needed for lambda, because lambda_m gets feedback only from nabla_m_ph part
template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::B_irreducible_vertex_m_contribution_from_nabla_m( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return 0.5 * state_vec.nabla_m_ph_to_xph( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_M_m( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
#if !defined(SBEa_APPROXIMATION)
    return 0.5*state_vec.M_m_ph_to_xph( W, K, w_in, n_in, w_out, n_out ) + state_vec.M_m( W, K, w_in, n_in, w_out, n_out );
#endif
    return 0.0;
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_m_double_counting( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return (-1.5*Model::B_m(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}


template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_nabla_sc( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return - state_vec.nabla_sc_pp_to_xph( W, K, w_in, n_in, w_out, n_out, t );
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_M_sc( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) { 
#if !defined(SBEa_APPROXIMATION)
    return - state_vec.M_sc_pp_to_xph( W, K, w_in, n_in, w_out, n_out );
#endif
    return 0.0;
}

template <typename Model, typename state_t> 
dcomplex rhs_postproc_sbe_t<Model, state_t>::vertx_m_contribution_from_sc_double_counting( const state_t& state_vec, int W, int w_in, int w_out, int K, int n_in, int n_out, const double t ) {
    return (-  Model::B_m(W, K, w_in, n_in, w_out, n_out)) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);
}
