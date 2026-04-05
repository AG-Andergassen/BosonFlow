#include <cmath>
#include <complex>
#include <string>
#include <vector>

#ifdef TEMP_FLOW
#include <params_physical.h>
#endif

#pragma once

#ifndef STATIC_CALCULATION
#define for_freq(FREQ_INDEX_NAME, INITIAL_VALUE, COMPARISON, INCREMENTION)  for (FREQ_INDEX_NAME = INITIAL_VALUE; COMPARISON; INCREMENTION)
#else
#define for_freq(FREQ_INDEX_NAME, INITIAL_VALUE, COMPARISON, INCREMENTION)  for (FREQ_INDEX_NAME=0, TEMP=0; TEMP < 1; TEMP ++) 
#endif



/**
 *  General constant definition
 */
const std::complex<double> I( 0.0, 1.0 );                       /**< Imaginary unit */
const double LN_10 = 2.30258509299;                             /**< Natural log of 10 */
const double CHOP_ERR = 1E-16;                                  /**< chop cosine, sine and weight to zero under this value */

namespace FrequenciesIRInfo
{
    const double ERROR = 1e-3;
    const double ENERGYSCALE = 1e6;
}

/**
 *  (Matsubara) frequency setting. All ranges defined as multiples of COUNT
 */
namespace FrequenciesCount
{
#ifdef MICKY_MOUSE
    const int COUNT = 2;
#else
    const int COUNT = COUNT_VAL;
#endif
    const int GRID_MULTIPLIER = 2;

    namespace Sig
    {
	const int POS_FERM = 10 * COUNT;                      /**< Amount of positive frequencies in self-energy grid */
	const int FERM = 2 * POS_FERM;            /**< Amount of frequencies in self-energy grid */
    }

    
    // ------- SBE ----    
    namespace w // screen interaction
    {
	const int POS_BOS = 64 * COUNT;
	const int BOS = 2 * POS_BOS + 1;
    }
    

    namespace lambda // yukawa coupling
    {
	const int POS_FERM = 2 * COUNT;
	const int FERM = 2 * POS_FERM;
	const int POS_BOS = 2 * COUNT;
	const int BOS = 2 * POS_BOS + 1;
    }

    namespace Integration
    {
      const int POS_RANGE =  32*COUNT; //4 * Sig::FERM;      //32*COUNT            /**< Positive range for internal integrations. Also the number of points in the matsubara frequency weight vector */
    }


    namespace bubble // the bubbles (GS and GG bubbles)
    {
	const int POS_BOS = 64 * COUNT;
	const int POS_FERM = Integration::POS_RANGE + w::POS_BOS/2; // 64 count
    }

    namespace M
    {
	// rest function
	const int POS_FERM = 2.0*COUNT;
	const int FERM = 2 * POS_FERM;
	const int POS_BOS = 2 * COUNT;
	const int BOS = 2 * POS_BOS + 1;
    }

    namespace Projected_M
    {
	// rest function
	const int POS_FERM = COUNT;
	const int FERM = 2 * POS_FERM;
	const int POS_BOS = 2 * COUNT;
	const int BOS = 2 * POS_BOS + 1;
    }

    namespace Projected_nabla
    {
	// rest function
	const int POS_FERM = 3 * COUNT;
	const int FERM = 2 * POS_FERM;
	const int POS_BOS = 2 * COUNT;
	const int BOS = 2 * POS_BOS + 1;
    }


  const int POS_1P_RANGE = w::POS_BOS + Integration::POS_RANGE; // 64 * COUNT + 32 * COUNT 

    const int LARGE = 5 * POS_1P_RANGE;  /**< large frequency used for the limit to infinity only for (SBE)*/

    const int TAIL_LENGTH = Integration::POS_RANGE / 5;                      /**< Length of tail used for fitting matsubara sum. @see */
    const int FIT_ORDER = 4;                                        /**< Fit tail function has exponents one lower than this constant */

    static_assert(bubble::POS_BOS == w::POS_BOS, "Bosonic frequencies count of the bubble and the screened interaction (w) must be the same!");
}
/*namespace OutputOptions
{
    const char* OUTPUT_DIRECTORY = "/gpfs/data/fs71925/kfraboul/STORAGE_OUTPUT/frgstationSBE2";
    }*/

// shifts the definition of the bosonic bare interaction by a constant. This is needed in case the the bosonic propagators/screened interactions w would be assigned zero otherwise (which would result in the freezing of the flow of the bosonic degrees of freedom).
const double BOSONIC_INTERACTION_M_SHIFT = 0.0;
const double BOSONIC_INTERACTION_SC_SHIFT =  0.0;
const double BOSONIC_INTERACTION_D_SHIFT = 0.0;
// --- Parquet. Todo : put them in namespace above 

// ----- Output ranges


#ifdef MICKY_MOUSE
const int K_DIM = 4;                                         /**< number of k-points in each k-dimension. e.g. 2D -> (x,y). should be even to contain (pi, pi) */
const int P_IN_K = 1;                                      /**< number of px-points for each ky-point (ratio of fourier-transform grid to coarse grid), symmetries work better if odd */
#else
const int K_DIM = K_DIM_VAL;                                         /**< number of k-points in each k-dimension. e.g. 2D -> (x,y) */
const int P_IN_K = P_IN_K_VAL;//41;                                      /**< number of px-points for each ky-point (ratio of fourier-transform grid to coarse grid), symmetries work better if odd */
#endif

// TODO: make the Hubbard atom (a 0d models) not regard this parameter
const double FORMFACTOR_SHELL_COUNT = FORMFACTOR_SHELL_COUNT_VAL;//1.99; // number of form factor shells to include

const int QN_COUNT = 1;                                         /**< Amount of possible tuples of the discrete quantum numbers */

const int EXTRA_LOOP_NUM = EXTRA_LOOP_NUM_VAL;                                       /**< Number of loops beyond 1-loop. E.g. 0 -> 1 loop + katanin. 1 -> 2 loop. 2 -> 3 loop ..etc*/
const std::string MULOOP_NUM_STRING( std::to_string(EXTRA_LOOP_NUM+1));


/**
 *  \brief precision settings for ODE solver and loop-convergence condition.
 *
 */
const double ERR_ABS = ERR_ABS_VAL;                                    /**< Absolute error for ODE routines */ 
const double ERR_REL = ERR_REL_VAL;                                    /**< Relative error for ODE routines */
const double ERR_LOOP_ABS = ERR_LOOP_ABS_VAL;		///< Absolute error for loop-condition
const double ERR_LOOP_REL = ERR_LOOP_REL_VAL;		///< Relative error for loop-condition

const double MAX_COUPLING = MAX_COUPLING_VAL;                               /**< If this coupling is exceeded in the vertex function the ODE solver will abort */

const int UPDATER_FREQ = 1;                                  /**<For every UPDATER_FREQth integration step an output file is created */
const int UPDATER_START = 0;					/**<Minimum integration step to start creating an output file */

//Self-Energy Correction Iteration
const double ABS_ERROR_SELFENERGY_ITERATIONS = ABS_ERROR_SELFENERGY_ITERATIONS_VAL;				///< Condition for Self-Energy CORRection ITeration loop termination -absolute  
const double ABS_ERROR_VERTEX_MULTILOOP_CORRECTIONS = 1e-4; 

const int SELFENERGY_ITERATIONS_MAX = SELFENERGY_ITERATIONS_MAX_VAL;					///< Maximum Number of Self-Energy CORRection ITeration loops



#ifndef DEBUG_EULER
#define   ERR_STEPPER       runge_kutta_cash_karp54
//#define   ERR_STEPPER       ssp_rk43_error_stepper
const std::string ERR_STEPPER_STRING( "runge_kutta_cash_karp54" );
//const std::string ERR_STEPPER_STRING( "runge_kutta_bogacki_shampine" );

#else
#define 	ERR_STEPPER 			euler
const std::string ERR_STEPPER_STRING( "euler" );
const int NSTEPS=NSTEPS_VAL;
const std::string NSTEPS_STRING( std::to_string(NSTEPS));
#endif

const std::vector<double> MANDATORY_TIMESTEP_PERCENTAGES = {};//{0.602, 0.301, 0.125, 0.0, -0.097, -0.176, -0.243};//{0.2, 0.4, 0.6, 0.8}; // if nonempty, ensures the ODE solver will include the flow steps corresponding to these percentages in the integration of the flow equation


//0.582879 #9 beta = 0.261
//0.1287  #10 beta = 0.743
// -0.056 #11 beta = 1.13
// -0.2422 #12 beta = 1.7466
// -0.30103 # final beta = 2.0

/*********************  INTERPOLATING FLOW ********************/

const int InputUirreducibleVertexPositiveFermFreqCount = 32;
const int InputUirreducibleVertexPositiveBosFreqCount = 32; // "positive" contains the zero frequency

const int Input_lambda_PositiveFermFreqCount = 150;
const int Input_lambda_PositiveBosFreqCount = 150; // "positive" contains the zero frequency

const std::string INPUT_REFERENCE_DATA_DIRECTORY = "./dmft_input/dmft_input_halffilling_U3_Beta5.h5"; 



/*********************      INPUT FILE      ********************/
#ifdef READIN
const std::string  INPUT_FILE_LOC( "./");
const std::string  INPUT_FILE_NAME( "input.h5");
#endif // ends READIN
