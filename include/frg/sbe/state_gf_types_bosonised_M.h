
/******************************************************************************************//** @file
 *  		
 * 	file: 		state_gf_types.h
 * 	contents:  	Definitions of containers needed for parquet-decomposition based fRG state 
 *                      ( wrappers around gf container )
 * 
 ****************************************************************************************************/


#pragma once

#include <frg/sbe/state_gf_types.h>


/// State types:


/// State types:
/**
 * \brief gf_w_t is the class for the "screened interaction" according to the BSE
 *
 * 2 bosonic dependencies W and K
 * complex valued
 */
enum class IW_M{ W, K, m, n }; 
template<typename Model> class gf_w_M_t : public gf< dcomplex, 4 > 	
{
   public:
      using base_t = gf< dcomplex, 4 >; 

      gf_w_M_t( int pos_bfreq_count_ = FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][Model::GetRefinedMomentaCount()][Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()])
   {}
      INSERT_COPY_AND_ASSIGN(gf_w_M_t)
}; 
template<typename Model> using idx_w_M_t = typename gf_w_M_t<Model>::idx_t; 

/**
 * \brief gf_lambda_t is the class for the Hedin vertex according to the SBE
 *
 * 2 bosonic dependencies W and K
 * 2 fermionic dependencies w and m
 * complex valued
 */
enum class ILAMBDA_M{ W, K, w, m, n };
template<typename Model> class gf_lambda_M_t : public gf< dcomplex, 5 > 		
{
   public:
      using base_t = gf< dcomplex, 5 >; 

      gf_lambda_M_t( int pos_bfreq_count_ = FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count, int pos_ffreq_count_ = FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][Model::GetRefinedMomentaCount()][ffreq(pos_ffreq_count_)][Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()])
   {}
      INSERT_COPY_AND_ASSIGN(gf_lambda_M_t)
}; 
template<typename Model> using idx_lambda_M_t = typename gf_lambda_M_t<Model>::idx_t; 


