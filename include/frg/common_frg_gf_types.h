#pragma once

#include <base_gf_types.h>
#include <frequencies/schemes.h>

// for the self-energy
enum class I1P{ w, k, s_in, s_out }; 
template<typename Model> class gf_1p_t : public  gf< dcomplex, 4 > 		///< Container type for one-particle correlation functions (self-energy), holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 4 >; 

 gf_1p_t( int pos_freq_count_ = FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count ): // IRBasis::GetSigPositiveFrequenciesCount()
      base_t( boost::extents[ffreq(pos_freq_count_)][Model::GetCoarseMomentaCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_t)
}; 
template<typename Model> using idx_1p_t = typename gf_1p_t<Model>::idx_t; 


/**
 * \brief gf_susc_t is the class for the susceptibility
 *
 * 2 bosonic dependencies W and K and 2 fermionic form factor indices n_in and n_out, and 4 spins
 * complex valued
 */

enum class ISUSC{ W, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out }; 
template<typename Model> class gf_susc_t : public gf< dcomplex, 8 > 		///< Container type for susceptibilities, holds dcomplex
{
 public:
    using base_t = gf< dcomplex, 8 >; 
    
 gf_susc_t( int pos_bfreq_count = 0 /*FrequencyDependenceScheme<Model>::w_BosonicGrid().positive_freqs_count*/ ): 
    base_t( boost::extents[bfreq(pos_bfreq_count)] 
	    [Model::GetRefinedMomentaCount()][Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()]
	    [Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()] )
	{}
    INSERT_COPY_AND_ASSIGN(gf_susc_t)
	}; 
template<typename Model> using idx_susc_t = typename gf_susc_t<Model>::idx_t; 

/**
 * \brief gf_polarisation_t is the class for the postprocessing "polarisation" according to the BSE
 *
 * 2 bosonic dependencies W and K
 * complex valued
 */
enum class IPOLARISATION{ W, K }; 
template<typename Model> class gf_polarisation_t : public gf< dcomplex, 2 > 	
{
   public:
      using base_t = gf< dcomplex, 2 >; 

      gf_polarisation_t( int pos_bfreq_count_ = 0 ):
	 base_t( boost::extents[bfreq(pos_bfreq_count_)][Model::GetRefinedMomentaCount()])
   {}
      INSERT_COPY_AND_ASSIGN(gf_polarisation_t)
}; 
template<typename Model> using idx_polarisation_t = typename gf_polarisation_t<Model>::idx_t; 
