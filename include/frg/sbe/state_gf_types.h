
/******************************************************************************************//** @file
 *  		
 * 	file: 		state_gf_types.h
 * 	contents:  	Definitions of containers needed for parquet-decomposition based fRG state 
 *                      ( wrappers around gf container )
 * 
 ****************************************************************************************************/


#pragma once

#include <frg/common_frg_gf_types.h>


/// State types:
/**
 * \brief gf_f_t is the class for the "Free energy"
 * no indices; just a complex number
 */

template<typename Model> class gf_f_t : public gf< dcomplex, 1 > 	
{
   public:
      using base_t = gf< dcomplex, 1 >; 

      gf_f_t():
	 base_t(boost::extents[1])
   {}
      INSERT_COPY_AND_ASSIGN(gf_f_t)
}; 
template<typename Model> using idx_f_t = typename gf_f_t<Model>::idx_t; 

/// State types:
/**
 * \brief gf_w_t is the class for the "screened interaction" according to the BSE
 *
 * 2 bosonic dependencies W and K
 * complex valued
 */
enum class IW{ W, K }; 
template<typename Model> class gf_w_t : public gf< dcomplex, 2 > 	
{
   public:
      using base_t = gf< dcomplex, 2 >; 

      gf_w_t( int pos_bfreq_count_ = FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][Model::GetRefinedMomentaCount()])
   {}
      INSERT_COPY_AND_ASSIGN(gf_w_t)
}; 
template<typename Model> using idx_w_t = typename gf_w_t<Model>::idx_t; 

/**
 * \brief gf_lambda_t is the class for the Hedin vertex according to the SBE
 *
 * 2 bosonic dependencies W and K
 * 2 fermionic dependencies w and m
 * complex valued
 */
enum class ILAMBDA{ W, K, w, m }; 
template<typename Model> class gf_lambda_t : public gf< dcomplex, 4 > 		
{
   public:
      using base_t = gf< dcomplex, 4 >; 

      gf_lambda_t( int pos_bfreq_count_ = FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count, int pos_ffreq_count_ = FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][Model::GetRefinedMomentaCount()][ffreq(pos_ffreq_count_)][Model::GetMomentumFormFactorsCount()])
   {}
      INSERT_COPY_AND_ASSIGN(gf_lambda_t)
}; 
template<typename Model> using idx_lambda_t = typename gf_lambda_t<Model>::idx_t; 

/**
 * \brief gf_M_t is the class for the U-irreducible but 2-particle-reducible vertex according to the SBE
 *
 * 2 bosonic dependencies W and K
 * 4 fermionic dependencies w, m, wp, and mp
 * complex valued
 */
enum class IM{ W, K, w, m, wp, mp };
template<typename Model> class gf_M_t : public gf< dcomplex, 6 >
{
   public:
      using base_t = gf< dcomplex, 6 >;

      gf_M_t( int pos_bfreq_count_ = FrequencyDependenceScheme<Model>::M_BosonicGrid().positive_freqs_count, int pos_ffreq_count_ = FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count ):
#ifndef BOSONISE_M_LAZY
             base_t( boost::extents[bfreq(pos_bfreq_count_)][Model::GetRefinedMomentaCount()][ffreq(pos_ffreq_count_)][Model::GetMomentumFormFactorsCount()][ffreq(pos_ffreq_count_)][Model::GetMomentumFormFactorsCount()])
#else
	     base_t( boost::extents[bfreq(pos_bfreq_count_)][Model::GetRefinedMomentaCount()][ffreq(pos_ffreq_count_)][Model::GetMomentumFormFactorsCount()][ffreq(1)][Model::GetMomentumFormFactorsCount()])
#endif
   {}
      gf_M_t( int pos_bfreq_count_, int pos_ffreq_count_, int pos_ffreq_count2_):
         base_t( boost::extents[bfreq(pos_bfreq_count_)][Model::GetRefinedMomentaCount()][ffreq(pos_ffreq_count_)][Model::GetMomentumFormFactorsCount()][ffreq(pos_ffreq_count2_)][Model::GetMomentumFormFactorsCount()])
   {}
      gf_M_t( int pos_bfreq_count_, int pos_ffreq_count_, int pos_ffreq_count2_, int ff_count1_, int ff_count2_):  
         base_t( boost::extents[bfreq(pos_bfreq_count_)][Model::GetRefinedMomentaCount()][ffreq(pos_ffreq_count_)][ff_count1_][ffreq(pos_ffreq_count2_)][ff_count2_])
   {}

      gf_M_t( int pos_bfreq_count_, int pos_ffreq_count_, int pos_ffreq_count2_, int ff_count1_, int ff_count2_, int bosonic_mom_count):  
         base_t( boost::extents[bfreq(pos_bfreq_count_)][bosonic_mom_count][ffreq(pos_ffreq_count_)][ff_count1_][ffreq(pos_ffreq_count2_)][ff_count2_])
   {}

      INSERT_COPY_AND_ASSIGN(gf_M_t)
};
template<typename Model> using idx_M_t = typename gf_M_t<Model>::idx_t;



/**
 * \brief gf_lambda_postproc_t is the class for the Hedin vertex according to the SBE at W = 0, and w = 1. Stores the postprocessed lambda
 *
 * 2 momentum dependencies K and m
 * complex valued
 */
enum class ISTATICLAMBDA{ K, m }; 
template<typename Model> class gf_static_lambda_t : public gf< dcomplex, 2 > 		
{
   public:
      using base_t = gf< dcomplex, 2 >; 

      gf_static_lambda_t( ):
	 base_t( boost::extents[Model::GetRefinedMomentaCount()][Model::GetMomentumFormFactorsCount()])
   {}
      INSERT_COPY_AND_ASSIGN(gf_static_lambda_t)
}; 
template<typename Model> using idx_static_lambda_t = typename gf_static_lambda_t<Model>::idx_t; 


enum class ISIG_SDE_INTEGRAND{ W, K, w, k };

template<typename Model> class gf_Sig_SDE_integrand_t : public gf< dcomplex, 4 > 		
{
   public:
      using base_t = gf< dcomplex, 4 >; 

      gf_Sig_SDE_integrand_t( ):
	base_t( boost::extents[bfreq(FrequenciesCount::Integration::POS_RANGE)][Model::GetRefinedMomentaCount()][ffreq(3)][Model::GetRefinedMomentaCount()])
   {}
      INSERT_COPY_AND_ASSIGN(gf_Sig_SDE_integrand_t)
}; 
template<typename Model> using idx_Sig_SDE_integrand_t = typename gf_Sig_SDE_integrand_t<Model>::idx_t; 



enum class INABLA{ W, K, w, m, wp, mp };
template<typename Model> using gf_nabla_t = gf_M_t<Model>;
template<typename Model> using idx_nabla_t = typename gf_nabla_t<Model>::idx_t;


// for the self-energy fluctuation diagnostics
enum class ISIG_POSTPROCESSING{ W, Q, w, k, s_in, s_out }; 
template<typename Model> class gf_Sig_postprocessing_t : public  gf< dcomplex, 6 > 		///< Container type for one-particle correlation functions (self-energy), holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 6 >; 

  gf_Sig_postprocessing_t( int pos_ffreq_count_ = FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count, int pos_bfreq_count_ = FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count - FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count, int mom_count = Model::GetCoarseMomentaCount() ): // IRBasis::GetSigPositiveFrequenciesCount()
    base_t( boost::extents[bfreq(pos_bfreq_count_)][mom_count][ffreq(2)][mom_count][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_Sig_postprocessing_t)
}; 
template<typename Model> using idx_Sig_postprocessing_t = typename gf_Sig_postprocessing_t<Model>::idx_t; 

