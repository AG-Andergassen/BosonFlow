
/************************************************************************************************//**
 *  		
 * 	file: 		mymath.cpp
 * 	contents:   	See mymath.h
 * 
 ****************************************************************************************************/


#include <mymath.h>
#include <iostream>
#include <Eigen/Core>
#include <Eigen/Eigenvalues>


std::vector<double> generate_tail_weights( int iMin, int tail_length, int fit_order )
{
   Eigen::MatrixXd MatFit(fit_order,fit_order); 
   MatFit.setZero(); 

   int iMax = tail_length + iMin; 

   for( int i = 0; i < fit_order; ++i )
      for( int j = 0; j < fit_order; ++j )
	 for( int k = iMin; k < iMax; ++k )
	    MatFit(i,j) += 1.0 / pow(k, i + j); 

   Eigen::MatrixXd MatInv = MatFit.inverse(); 

   std::vector<double> w_array(tail_length,0.0); 
   for( int i = 0; i < tail_length; ++i )
      for( int l = 0; l < fit_order; ++l )
	 w_array[i] += MatInv(0,l) / pow(i + iMin, l); 

   std::vector<double> tail_weights(tail_length,0.0);
   for( int i = 0; i < tail_length; ++i )
      for( int l = i; l < tail_length; ++l )
	 tail_weights[i] += w_array[l]; 

   return tail_weights; 
}

void generate_weights( gf<double, 1> &weight_vecs_1d, int iMin, int tail_length, int fit_order ) /**< weight vectors used to correct the matsubara-freq summation on a finite range. Semi-analytical fitting */   
{
   std::vector<double> tail_weights = generate_tail_weights( iMin, tail_length, fit_order ); 

   weight_vecs_1d.init( []( const gf<double, 1>::idx_t& idx ){ return 1.0; } ); 

#ifdef STATIC_CALCULATION         
   return;
#endif

   for( int i = 0; i < tail_length; ++i )
   {
      weight_vecs_1d[iMin + i] = tail_weights[i]; 
      weight_vecs_1d[-iMin - i - 1] = tail_weights[i]; 
   } 
}

void generate_2d_weights( gf<double, 2> &weight_vecs_2d, int iMin, int tail_length, int fit_order ) /**< weight vectors used to correct the double matsubara-freq summation on a finite range. Semi-analytical fitting */ 
{
   std::vector<double> tail_weights = generate_tail_weights( iMin, tail_length, fit_order ); 

   weight_vecs_2d.init( []( const gf<double, 2>::idx_t& idx ){ return 1.0; } ); 
   
#ifdef STATIC_CALCULATION         
   return;
#endif

   for( int i = 0; i < tail_length; ++i )
   {
      for( int j = -iMin - i - 1; j <= iMin + i; ++j )
      {
	 weight_vecs_2d[iMin + i][j] = tail_weights[i]; 
	 weight_vecs_2d[-iMin - i - 1][j] = tail_weights[i]; 
	 weight_vecs_2d[j][iMin + i] = tail_weights[i]; 
	 weight_vecs_2d[j][-iMin - i - 1] = tail_weights[i]; 
      }
   }
}

// derivative of the Fermi-function with respect to temperature (needed for the static SG bubble)
double d_nF_over_dT(double z, double T)
{
    return z / (2 *T * T * (1.0 + std::cosh(z/T)));
}

double d_nF_over_dz(double z, double T)
{
    return - (1.0/T) /(2 + 2 * std::cosh(z/T));
}

double d2_nF_over_dTdz(double z, double T)
{
    const double term = 1.0/(2 *T * T * (1.0 + std::cosh(z/T)));
    return term - z * 2 * T * T * std::sinh(z/T)/term/term;
}
