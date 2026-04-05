
/******************************************************************************************//** @file
 *  		
 * 	file: 		mymath.h
 * 	contents:  	Contains useful mathematical functions / routines
 * 
 ****************************************************************************************************/

#pragma once

#include <base_gf_types.h>
#include <complex>
#include <cmath>

//---------------- MatReal get indices 1d array to 3d array -----------------------------------
// [Deprecated function below]
// it does not really matter what the grid here looks like.
template <typename Model>
inline int Matidx(const int rx, const int ry){
    return (rx+10*Model::GetFFTDim()) % Model::GetFFTDim()+Model::GetFFTDim()*((ry+10*Model::GetFFTDim()) % Model::GetFFTDim());
}


template <typename Model>
inline int GetBubbleFFTRealCoordXFromIdx(const int R_idx) // the shift is to let R = (0, 0) be the middle
{
    return R_idx/Model::GetFFTDim() - Model::GetFFTDim()/2;
}


template <typename Model>
inline int GetBubbleFFTRealCoordYFromIdx(const int R_idx)
{
    return R_idx % Model::GetFFTDim() - Model::GetFFTDim()/2;
}


// ---- Useful functions

inline int div2_ceil( const int W )
{
   return (W - 1000000)/2 + 500000; 
}

inline int div2_floor( const int W )
{
   return (W + 1000000)/2 - 500000; 
}

template <typename T> double sgn( const T val) 
{
   return (T(0) < val) - (val < T(0));
}


// generate_sum_function provides the fitting of a finite matsubara sum in order to get the extension of the sum up to infinity
// The Least-Square fitting uses the following fitting function : sum_{i=0}^{fit_order} a_i (x)^{-i}
// This requires the minimization of the chisquare which has been done in a semi-analytic way (see notes by Georg Roehringer)
// generate_sum_func returns sumfit of a given single-variable function 
std::vector<double> generate_tail_weights( int iMin, int tail_length, int fit_order ); 
void generate_weights( gf<double, 1> &weight_vecs_1d, int iMin, int tail_length, int fit_order ); 
void generate_2d_weights( gf<double, 2> &weight_vecs_2d, int iMin, int tail_length, int fit_order ); 


template< typename value_t >
boost::function< value_t ( boost::function< value_t ( int i ) > ) > generate_sum_func( int iMin, int tail_length, int fit_order )
{
   std::vector<double> tail_weights = generate_tail_weights( iMin, tail_length, fit_order ); 

   return [tail_weights,iMin,tail_length]( boost::function< value_t ( int i ) > func )
   {
      value_t val; 
      for( int i = -iMin; i <= iMin; ++i )
	 val += func(i); 

      for( int i = 1; i < tail_length; ++i )
	 val += ( func(iMin + i) + func(-iMin -i) ) * tail_weights[i]; 

      return val; 
   }; 
}


/// folds integer x on the periodic interval [0, periodicity) 
inline unsigned fold(const int x, const unsigned periodicity)
{
    return (x % periodicity + periodicity) % periodicity;
}

inline double nF(double z, double T)
{
    return 1/(std::exp(z/T) + 1);
}

double d_nF_over_dz(double z, double T); 

double d_nF_over_dT(double z, double T);

double d2_nF_over_dTdz(double z, double T); 



