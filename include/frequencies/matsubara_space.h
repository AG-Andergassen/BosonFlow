#pragma once

#include <vector>
#include <string>
#include <iostream>
#include <params_technical.h>
#include <functional>
#include <mymath.h>
#include <params_physical.h>

/**
 *  \class matsubara_space_t
 *  \brief Base class for the Matsubara frequencies
 */

//enum class Statistics{BOSE = 0, FERMI = 1};

class matsubara_space_t
{
    public:
        matsubara_space_t( unsigned int pos_freq_count, double step_size, Statistics zeta);

        inline const double& operator[]( const int w ){
 	        return grid_points[w + pos_freq_count]; 
        }
        
        inline const double * data() const{
 	        return grid_points.data(); 
        }
        
        unsigned int get_pos_freq_count() const;
        double get_step_size() const;

	void print( std::string fname);	                /**< Prints frequency grid to file */
 
    private:
        std::vector<double> grid_points;  		/**< Values of grid */
        const int 	 pos_freq_count;		    /**< Number of elements in grid */
        const double step_size;			        /**< Step size of equidistant grid */
        
};

/**
 *  \class fermionic_matsubara_space_t
 *  \brief Class for the fermionic Matsubara frequencies
 */

class fermionic_matsubara_space_t : public matsubara_space_t
{	
    public: 	
        fermionic_matsubara_space_t( unsigned int pos_freq_count, double step_size ) : matsubara_space_t(pos_freq_count, step_size, Statistics::FERMI){};	/**< Initialing constructor for equidistant fermionic frequency grid */
 
        void print(){matsubara_space_t::print("f_grid.dat");};	            /**< Prints frequency grid to file "f_grid.dat" */
};

/**
 *  \class bosonic_matsubara_space_t
 *  \brief Class for the bosonic Matsubara frequencies
 */

class bosonic_matsubara_space_t : public matsubara_space_t
{	
    public:
        bosonic_matsubara_space_t( unsigned int pos_freq_count, double step_size ) : matsubara_space_t(pos_freq_count, step_size, Statistics::BOSE){};	/**< Initialing constructor for equidistant bosonic frequency grid */

        void print(){matsubara_space_t::print("bos_grid.dat");};	            /**< Prints frequency grid to file "bos_grid.dat" */
};


// should be moved inside frequency grid class as  methods
/**
 *  \brief Useful inline functions to evaluate Matsubara frequecies given index
 */

inline double w_val( int w_idx )					/**< Return value of fermionic matsubara frequency for a given index */
{
   return M_PI / BETA * ( 2 * w_idx + 1 );
}

inline double W_val( int W_idx )					/**< Return value of bosonic matsubara frequency for a given index */
{
   return 2.0 * M_PI / BETA * W_idx;
}


inline double inv_freq_sum_normalisation( const double t)
{
#ifdef TEMP_FLOW
   return 1.;
#else
   return BETA;
#endif
}

inline double Beta( const double t)
{
#ifdef TEMP_FLOW
   return exp(-1.*LN_10 * t);
#else
   return BETA;
#endif
}

inline double w_val( int w_idx, const double t )					/**< Return value of fermionic matsubara frequency for a given index */
{
   return M_PI / Beta(t) * ( 2 * w_idx + 1 );
}

inline double W_val( int W_idx, const double t )					/**< Return value of bosonic matsubara frequency for a given index */
{
   return 2.0 * M_PI / Beta(t) * W_idx;
}
