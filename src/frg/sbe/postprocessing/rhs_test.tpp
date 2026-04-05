#include <frg/sbe/postprocessing/rhs.h>
#include <mymath.h>
#include <cmath>
#include <complex>
#include <frg/sbe/state.h>
#include <frg/sbe/postprocessing/symmetries.h>

/*

template<typename Model> 
void rhs_postproc_sbe_t<Model>::test( const state_frg_sbe_t<Model>& state, state_postproc_t<Model>& postproc_state, const double t )
{
    // call only after initialising the state

    gf_susc_t<Model> suscvert_sc_contribution_from_Virr;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_susc_sc_Ptr()->init_batched( suscvert_sc_contribution_from_Virr,  [&state, &bubble_pp, t]( const idx_susc_t<Model>& idx ){ return eval_susc_sc_contribution_from_Virr( idx, state, bubble_pp, t ); } );
    
    std::cout << "Equality: " << norm(suscvert_sc_contribution_from_Virr 
		      - postproc_state.gf_suscvert_sc_contribution_from_M_m()
		      - postproc_state.gf_suscvert_sc_contribution_from_m_double_counting()
		      - postproc_state.gf_suscvert_sc_contribution_from_M_sc()
		      - postproc_state.gf_suscvert_sc_contribution_from_sc_double_counting()
		      - postproc_state.gf_suscvert_sc_contribution_from_M_d()
		      - postproc_state.gf_suscvert_sc_contribution_from_d_double_counting()
		      - postproc_state.gf_suscvert_sc_contribution_from_bare()) << std::endl;


    gf_lambda_t<Model> lambda_sc_contribution_from_Virr(0);
    lambda_sc_contribution_from_Virr.init([&state, &bubble_pp, t]( const idx_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_Virr( idx, state, bubble_pp, t ); });

    gf_polarisation_t<Model> polarisation_sc_contribution_from_Virr;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( polarisation_sc_contribution_from_Virr,  [&state, &bubble_pp, &lambda_sc_contribution_from_Virr,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_sc_contribution_lambda( idx, state, bubble_pp, t , lambda_sc_contribution_from_Virr); } );

    gf_polarisation_t<Model> polarisation_sc_contribution_from_Virr2;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( polarisation_sc_contribution_from_Virr2,  [&state, &bubble_pp, &suscvert_sc_contribution_from_Virr,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_pp, t , suscvert_sc_contribution_from_Virr); } );

    std::cout << "Equality of polarisations Virr: " << norm(polarisation_sc_contribution_from_Virr2 - polarisation_sc_contribution_from_Virr) << std::endl;

    // nabla sc
    gf_lambda_t<Model> lambda_sc_contribution_from_nabla_d(0);
    lambda_sc_contribution_from_nabla_d.init([&state, &bubble_pp, t]( const idx_lambda_t<Model>& idx ){ return eval_lambda_sc_contribution_from_nabla_d( idx, state, bubble_pp, t ); });

    gf_polarisation_t<Model> polarisation_sc_contribution_from_nabla_d;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( polarisation_sc_contribution_from_nabla_d,  [&state, &bubble_pp, &lambda_sc_contribution_from_nabla_d,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_sc_contribution_lambda( idx, state, bubble_pp, t , lambda_sc_contribution_from_nabla_d); } );

    gf_polarisation_t<Model> polarisation_sc_contribution_from_nabla_d2;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( polarisation_sc_contribution_from_nabla_d2,  [&state, &bubble_pp, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_pp, t , postproc_state.gf_suscvert_sc_contribution_from_nabla_d() ); } );

    std::cout << "Equality of polarisations sc nabla_d: " << norm(polarisation_sc_contribution_from_nabla_d - polarisation_sc_contribution_from_nabla_d2) << std::endl;


    // nabla d
    gf_lambda_t<Model> lambda_d_contribution_from_nabla_sc(0);
    lambda_d_contribution_from_nabla_sc.init([&state, &bubble_ph, t]( const idx_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_nabla_sc( idx, state, bubble_ph, t ); });

    gf_polarisation_t<Model> polarisation_d_contribution_from_nabla_sc;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( polarisation_d_contribution_from_nabla_sc,  [&state, &bubble_ph, &lambda_d_contribution_from_nabla_sc,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_d_contribution_lambda( idx, state, bubble_ph, t , lambda_d_contribution_from_nabla_sc); } );

    gf_polarisation_t<Model> polarisation_d_contribution_from_nabla_sc2;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( polarisation_d_contribution_from_nabla_sc2,  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_ph, t , postproc_state.gf_suscvert_d_contribution_from_nabla_sc() ); } );

    std::cout << "Equality of polarisations d nabla_sc: " << norm(polarisation_d_contribution_from_nabla_sc - polarisation_d_contribution_from_nabla_sc2) << std::endl;



    // nabla m
    gf_lambda_t<Model> lambda_d_contribution_from_nabla_sc(0);
    lambda_d_contribution_from_nabla_sc.init([&state, &bubble_ph, t]( const idx_lambda_t<Model>& idx ){ return eval_lambda_d_contribution_from_nabla_sc( idx, state, bubble_ph, t ); });

    gf_polarisation_t<Model> polarisation_d_contribution_from_nabla_sc;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( polarisation_d_contribution_from_nabla_sc,  [&state, &bubble_ph, &lambda_d_contribution_from_nabla_sc,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_d_contribution_lambda( idx, state, bubble_ph, t , lambda_d_contribution_from_nabla_sc); } );

    gf_polarisation_t<Model> polarisation_d_contribution_from_nabla_sc2;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_d_Ptr()->init_batched( polarisation_d_contribution_from_nabla_sc2,  [&state, &bubble_ph, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_ph, t , postproc_state.gf_suscvert_d_contribution_from_nabla_sc() ); } );

    std::cout << "Equality of polarisations d nabla_sc: " << norm(polarisation_d_contribution_from_nabla_sc - polarisation_d_contribution_from_nabla_sc2) << std::endl;


    // bubble
    gf_lambda_t<Model> lambda_contribution_from_1(0);
    lambda_contribution_from_1.init([&state, &bubble_pp, t]( const idx_lambda_t<Model>& idx ){ return eval_lambda_contribution_from_1( idx, state, bubble_pp, t ); });

    gf_polarisation_t<Model> polarisation_sc_contribution_from_1;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( polarisation_sc_contribution_from_1,  [&state, &bubble_pp, &lambda_contribution_from_1, t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_sc_contribution_lambda( idx, state, bubble_pp, t , lambda_contribution_from_1); } );

    gf_polarisation_t<Model> polarisation_sc_contribution_from_1_2;
    SymmetriesfRGSBEPostprocessing<Model>::IdxEquivClasses_polarisation_sc_Ptr()->init_batched( polarisation_sc_contribution_from_1_2,  [&state, &bubble_pp, &postproc_state,  t]( const idx_polarisation_t<Model>& idx ){ return eval_polarisation_contribution_suscvert( idx, state, bubble_pp, t , postproc_state.gf_suscbubble_sc() ); } );

    std::cout << "Equality of polarisations 1: " << norm(polarisation_sc_contribution_from_1 - polarisation_sc_contribution_from_1_2) << std::endl;

}

*/




// instantiate class for the particular models
#define Instantiate(MODEL) template class rhs_postproc_sbe_t<MODEL>;
WITH_MODELS(Instantiate)
