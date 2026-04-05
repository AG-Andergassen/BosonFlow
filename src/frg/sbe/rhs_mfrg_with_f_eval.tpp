#include <frg/sbe/rhs_mfrg_with_f.h>
#include <mymath.h> //should be added from rhs_1lfrg.cpp, right?
#include <cmath> //should be added from rhs_1lfrg.cpp, right?
#include <complex> //should be added from rhs_1lfrg.cpp, right?
#include <frg/sbe/symmetries.h> //should be added from rhs_1lfrg.cpp, right?
#include <frg/flows.h> //should be added from rhs_1lfrg.cpp, right?
#include <models/square_hubbard.h>



// free energy
template <typename Model, typename state_t > 
dcomplex rhs_sbe_mfrg_with_F_t<Model, state_t>::eval_F_m( const idx_f_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, double t)
{
    dcomplex val( 0.0, 0.0 );

    //Hatree term of Self energy
    double sig_h = -state.eval_filling(0, 0) * Model::vertex_4pt_bare(00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t);

    for( int p_ = 0; p_ < Model::GetFineMomentaCount(); ++p_ )
    for_freq( int w_, -FrequenciesCount::Integration::POS_RANGE, w_ < FrequenciesCount::Integration::POS_RANGE, ++w_ )
    {
      
      const int k_ = Model::GetCoarseMomentumIdxFromFine(p_);
      dcomplex sig = state.SigMat(w_, k_)(0, 0) + sig_h;
      dcomplex g = Gvec[w_][p_](0,0);
      val += - Svec[w_][p_](0,0)*sig/(1.0 + sig*g);
    }
//If the cut off is not Omega flow there is a 1/nu term that not convergeces which needs to be taken into account like that
#ifndef OMEGA_FLOW
    val += state.eval_filling(0, 0) * fRGFlowScheme<Model>::U_MultiplicativeCutoff(t) * Model::vertex_4pt_bare(00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00, 00)/2;
#endif

    return val *= 2.0/inv_freq_sum_normalisation(t)/Model::GetFineMomentaCount();     /*< Normalize freq/momentum summation; times 2 because of spin */                
}
