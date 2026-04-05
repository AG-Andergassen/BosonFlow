
/******************************************************************************************//** @file
 *  		
 * 	file: 		base_gf_types.h
 * 	contents:  	Definitions of base containers used around the code ( wrapper around gf container )
 * 
 ****************************************************************************************************/

#pragma once

#include <complex>
#include <Eigen/Core>
#include <vector>
#include <boost/multi_array.hpp>

#include <params_technical.h>
#include <gf.h>

#include <frequencies/schemes.h>

inline extent_range ffreq( int n, bool definite_frequency_dependence = false )
{
    if (definite_frequency_dependence)
	return extent_range( -n, n ); // regardless of whether calculation is static or not

#ifndef STATIC_CALCULATION
    return extent_range( -n, n );   // range [ -n , n )
#else
    return extent_range(0, 1);
#endif
}

inline extent_range bfreq( int n, bool definite_frequency_dependence = false )
{
    if (definite_frequency_dependence)
	return extent_range( -n, n + 1 ); // regardless of whether calculation is static or not

#ifndef STATIC_CALCULATION
    return extent_range( -n, n + 1 );   // range [ -n , n + 1 )
#else
    return extent_range(0, 1);
#endif
}


using dcomplex = std::complex<double>;						///< Complex double type
// todo: pass parameters from Model.. For that, need them to be constexpr
using MatQN = Eigen::Matrix<dcomplex, QN_COUNT, QN_COUNT, Eigen::RowMajor>;	///< Complex matrix representing the discrete quantum number structure
using MatQNQN = Eigen::Matrix<dcomplex, QN_COUNT*QN_COUNT, QN_COUNT*QN_COUNT, Eigen::RowMajor>;	///< Complex matrix representing the discrete quantum number structure of two-particle function

using MatReal = Eigen::Matrix<dcomplex, Eigen::Dynamic, 1>;	///<FFT_DIM * FFT_DIM Complex matrix representing the real grid structure
using MatPatch = Eigen::Matrix<dcomplex, Eigen::Dynamic, 1>;	///< PATCH_COUNT Complex matrix representing the real grid structure

#define INSERT_COPY_AND_ASSIGN(X) 					\
X( const X & gf_obj ):    						\
   base_t( gf_obj )							\
{}       								\
X( X && gf_obj ) noexcept :							\
   base_t( std::move(gf_obj) ) 						\
{}      								\
X & operator=( const X & gf_obj )					\
{									\
   if( this == &gf_obj ) \
        return *this; \
   base_t::operator=( gf_obj ); 					\
   return *this; 							\
} 									\
X & operator=( X && gf_obj ) noexcept						\
{									\
   base_t::operator=( std::move( gf_obj) ); 				\
   return *this; 							\
}\
~ X(){\
}


/**
 * \brief gf_1p_mat_t is the class for the Green's function on the fine momentum grid
 *
 *  1 fermionic dependency w and p
 *  the value is a matrix with QN_COUNT*QN_COUNT entries
 */

// Container and index types
enum class I1PMAT{ w, p }; 
template<typename Model> class gf_1p_mat_t : public gf< MatQN, 2 > 			///< Matrix-valued container type for one-particle correlation functions used to store G and S defined on the fine momentum grid to be fourier transformed for bubble calculation
{
   public:
      using base_t = gf< MatQN, 2 >; 

 gf_1p_mat_t( int pos_freq_count_, int patch_count_, bool definite_frequency_dependence = false ):
      base_t( boost::extents[ffreq(pos_freq_count_, definite_frequency_dependence)][patch_count_] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_mat_t)
}; 
template<typename Model> using idx_1p_mat_t = typename gf_1p_mat_t<Model>::idx_t;  

enum class I1PMATREAL{ w, s_in, s_out }; 
template<typename Model> class gf_1p_mat_real_t : public gf< MatReal, 3 > 			///< Matrix-valued container type for one-particle correlation functions used to store fourier-transformed G and S for bubble calculation using the convolution theorem
{
   public:
      using base_t = gf< MatReal, 3 >; 

 gf_1p_mat_real_t( int pos_freq_count_= FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, bool definite_frequency_dependence = false ):
      base_t( boost::extents[ffreq(pos_freq_count_, definite_frequency_dependence)][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_1p_mat_real_t)
}; 
template<typename Model> using idx_1p_mat_real_t = typename gf_1p_mat_real_t<Model>::idx_t;   

enum class I1PREAL{ w, k, s_in, s_out }; 

template<typename Model> class gf_1p_real_t : public gf< dcomplex, 4 > 			///< Container type for one-particle correlation functions used only for output.cpp to store G and S, holds dcomplex
{
 public:
    using base_t = gf< dcomplex, 4 >; 
      
 gf_1p_real_t( int pos_freq_count_= FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, bool definite_frequency_dependence = false ):
    base_t( boost::extents[ffreq(pos_freq_count_, definite_frequency_dependence)][Model::GetFineMomentaCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()] ) 
   {}


      INSERT_COPY_AND_ASSIGN(gf_1p_real_t)
}; 
template<typename Model> using idx_1p_real_t = typename gf_1p_real_t<Model>::idx_t;   

enum class ISIGKMAT{ w, s_in, s_out }; // only used for parquet code
template<typename Model> class gf_Sig_kMat_t : public  gf< MatPatch, 3 > 		///< Container type for one-particle correlation functions, stores self-energy as a matrix for self-energy SDE flow, holds dcomplex
{
   public:
      using base_t = gf< MatPatch, 3 >; 

 gf_Sig_kMat_t( int pos_freq_count_ = FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count ):
	 base_t( boost::extents[ffreq(pos_freq_count_)][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_Sig_kMat_t)
}; 
template<typename Model> using idx_Sig_kMat_t = typename gf_Sig_kMat_t<Model>::idx_t; 

enum class IBUBMAT{ W, w, m, n, s1, s1p, s2, s2p }; 
template<typename Model> class gf_bubble_mat_t : public gf< MatPatch, 8 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< MatPatch, 8 >; 

 gf_bubble_mat_t( int pos_bfreq_count_= FrequencyDependenceScheme<Model>::BubbleTU_BosonicGrid().positive_freqs_count, int pos_freq_count_= FrequencyDependenceScheme<Model>::BubbleTU_FermionicGrid().positive_freqs_count, bool definite_frequency_dependence = false):
      base_t( boost::extents[bfreq(pos_bfreq_count_, definite_frequency_dependence)][ffreq(pos_freq_count_, definite_frequency_dependence)][Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_bubble_mat_t)
}; 
template<typename Model> using idx_bubble_mat_t = typename gf_bubble_mat_t<Model>::idx_t;

enum class IBUB{ W, w, K, n_in, n_out, s1_in, s2_in, s1_out, s2_out }; 
template<typename Model> class gf_bubble_t : public gf< dcomplex, 9 > 			///< Matrix-valued container type for one-particle correlation functions
{
   public:
      using base_t = gf< dcomplex, 9 >; 

 gf_bubble_t( int pos_bfreq_count_= FrequencyDependenceScheme<Model>::BubbleTU_BosonicGrid().positive_freqs_count, int pos_freq_count_= FrequencyDependenceScheme<Model>::BubbleTU_FermionicGrid().positive_freqs_count, bool definite_frequency_dependence = false):
      base_t( boost::extents[bfreq(pos_bfreq_count_, definite_frequency_dependence)][ffreq(pos_freq_count_, definite_frequency_dependence)][Model::GetRefinedMomentaCount()][Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_bubble_t)
}; 
template<typename Model> using idx_bubble_t = typename gf_bubble_t<Model>::idx_t;

enum class ISUSC_POSTPROC{ K }; 
class gf_susc_postproc_t : public gf< dcomplex, 1 > 		///< Container type for post-processing susceptibilities, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 1 >; 

      gf_susc_postproc_t( ):
	 base_t( boost::extents[3] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_susc_postproc_t)
}; 
using idx_susc_postproc_t = gf_susc_postproc_t::idx_t; 


// todo: move the two classes below to tu_projections.h

// Class for the bare vertex in the form factors + band basis -> only needed for input to the frg flow 
enum class VERT_BARE_FF{n_in, n_out, s1_in, s2_in, s1_out, s2_out }; 
template<typename Model> class gf_vert_bare_ff_t : public gf< dcomplex, 6 > 		///< Container type for the bare 2-particle vertex
{
   public:
      using base_t = gf< dcomplex, 6 >; 

      gf_vert_bare_ff_t():
	 base_t( boost::extents[Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()]
	       [Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()][Model::GetQuantumNumbersCount()] )
   {}
      INSERT_COPY_AND_ASSIGN(gf_vert_bare_ff_t)
}; 
template<typename Model> using idx_vert_bare_ff_t = typename gf_vert_bare_ff_t<Model>::idx_t; 
/*
// Class for the calculating the matrix M_{n, m}(k) = (f_n*f^*_m)(k). Used for writing the bare vertex in a given channel  in ff space; V^0_{m, n}(q) = \int_{k-k'} M_{n, m}(k-k')V^0_{k,k'}(q). Valid if V^0 depends on k-k' instead of k and k' individually, which will be the case for density density interactions.
enum class IFF_CONVOLUTION{ n, m, k}; 

template<typename Model> class gf_ff_convolution_t : public gf< dcomplex, 3 > 		///< Container type for the projection operators, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 3>; 

      gf_proj_matrix_t():
      base_t( boost::extents[Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()][Model::GetRefinedMomentaCount()])
   {}
      INSERT_COPY_AND_ASSIGN(gf_ff_convolution_t)
}; 
template<typename Model> using idx_ff_convolution_t = typename gf_ff_convolution_t<Model>::idx_t; 
*/


// Class for the precalculation of the projection matrices -> only needed for input to the frg flow 
enum class PROJ_MATRIX{ K_in, K_out, m, n, mp, np}; 

template<typename Model> class gf_proj_matrix_t : public gf< dcomplex, 6 > 		///< Container type for the projection operators, holds dcomplex
{
   public:
      using base_t = gf< dcomplex, 6>; 

      gf_proj_matrix_t():
      base_t( boost::extents[Model::GetRefinedMomentaCount()][Model::GetRefinedMomentaCount()]
	       [Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()][Model::GetMomentumFormFactorsCount()])
   {}
      INSERT_COPY_AND_ASSIGN(gf_proj_matrix_t)
}; 
template<typename Model> using idx_proj_matrix_t = typename gf_proj_matrix_t<Model>::idx_t; 


class gf_double_t : public gf< double, 1 > 			///< double-valued gf container; necessary to include scalar state elements
{
   public:
      using base_t = gf< double, 1 >; 

 gf_double_t():
	 base_t( boost::extents[extent_range(0, 1)] ) 
   {}
      INSERT_COPY_AND_ASSIGN(gf_double_t)
}; 

