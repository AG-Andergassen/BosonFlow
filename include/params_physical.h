
/******************************************************************************************//** @file
 *
 * 	file: 		params_physical.h
 * 	contents:  	Declare global parameters as extern
 *
 ****************************************************************************************************/


#pragma once

#include <vector>
#include <map>
#include <string>
#include <array>
#include <runtime_config.h>

// Choose from: ChainHubbard, SquareHubbard, FCCHubbard, TriangularHubbard, SquareHubbardHolstein, HubbardAtom, AndersonImpurity, AndersonImpurityHolstein, SquareHubbardLongRange, SquareHubbardPeierls (see models/concrete_available_models.h)
// Set via -D THE_MODEL_VAL=<ModelName> in config.mk; falls back to HubbardAtom if unset.
#ifndef THE_MODEL_VAL
    #define THE_MODEL HubbardAtom
#else
    #define THE_MODEL THE_MODEL_VAL
#endif


// Runtime-configurable BETA (set via command-line or defaults to 5.0)
#define BETA RuntimeConfig::BETA


namespace SquareHubbardParams
{
    inline const double& UINT = RuntimeConfig::UINT;     ///< on site interaction

    inline const double& VINT = RuntimeConfig::VINT;     ///< nearest neighbor interaction
  
    inline const double& MU = RuntimeConfig::MU;         ///< chemical potential (0 = half-filling)
  
    inline const double& T_PRIME = RuntimeConfig::T_PRIME; ///< next nearest neighbor hopping
    
    inline const std::map<std::string, double>& REFINE_AT_POINTS = RuntimeConfig::REFINE_AT_POINTS; // e.g. --refine-at idx_00:0.02
}

// (k_s, U) =  (0.1, 19.0423), (1.0, 11.1275), (5.0, 1.4126), (10, 0.38334), (100, 0.00394666), (3, 3.32506), (7, 0.759622), 


namespace SquareHubbardLongRangeParams
{
    inline const double& Q = RuntimeConfig::Q;       /// < long range coupling
    inline const double& K_s = RuntimeConfig::K_s;   /// < "thomas-fermi"-like screening
}


namespace SquareHubbardHolsteinParams
{   
    inline const std::string PHONON_TYPE = "acoustic";//"optical"; //"acoustic" 
    inline const double& g0 = RuntimeConfig::g0;       ///< Holstein e-p coupling
    inline const double& OMEGA0 = RuntimeConfig::OMEGA0; ///< adiabiticity parameter/phonon frequency

    // remaining params are read from SquareHubbardParams
    inline struct V_ph_proxy { operator double() const { return 2.0*RuntimeConfig::g0*RuntimeConfig::g0/RuntimeConfig::OMEGA0; } } V_ph;
}
 
namespace SquareHubbardPeierlsParams
{   
    inline const std::string PHONON_TYPE = "optical";//"optical"; //"acoustic" 
    inline const double& g0 = RuntimeConfig::g0;       ///< Holstein e-p coupling
    inline const double& OMEGA0 = RuntimeConfig::OMEGA0; ///< adiabiticity parameter/phonon frequency

    // remaining params are read from SquareHubbardParams
    inline struct V_ph_proxy { operator double() const { return 2.0*RuntimeConfig::g0*RuntimeConfig::g0/RuntimeConfig::OMEGA0; } } V_ph;
}

namespace TriangularHubbardParams
{
    inline const double& UINT = RuntimeConfig::UINT;
    inline const double& VINT = RuntimeConfig::VINT;

    inline const double& WINT = RuntimeConfig::WINT;
    inline const double& WWINT = RuntimeConfig::WWINT;

    inline const double& MU = RuntimeConfig::MU;

    inline const double& T_PRIME = RuntimeConfig::T_PRIME;
    inline const double& T_PRIME_PRIME = RuntimeConfig::T_PRIME_PRIME;

    inline const std::map<std::string, double>& REFINE_AT_POINTS = RuntimeConfig::REFINE_AT_POINTS;
}

namespace ChainHubbardParams
{
    inline const double& UINT = RuntimeConfig::UINT;
    inline const double& VINT = RuntimeConfig::VINT;
    
    inline const double& MU = RuntimeConfig::MU;
    inline const double& T_PRIME = RuntimeConfig::T_PRIME;
    
    inline const std::map<std::string, double>& REFINE_AT_POINTS = RuntimeConfig::REFINE_AT_POINTS;
}

namespace FCCHubbardParams
{
    inline const double& UINT = RuntimeConfig::UINT;
    
    inline const double& MU = RuntimeConfig::MU;
    inline const double& T_PRIME = RuntimeConfig::T_PRIME;
    
    inline const double& T_PRIME_PRIME = RuntimeConfig::T_PRIME_PRIME;

    inline const std::map<std::string, double>& REFINE_AT_POINTS = RuntimeConfig::REFINE_AT_POINTS;
}


namespace HubbardAtomParams
{
    inline const double& UINT = RuntimeConfig::UINT;
    
    inline const double& MU = RuntimeConfig::MU;
    
    inline const std::map<std::string, double>& REFINE_AT_POINTS = RuntimeConfig::REFINE_AT_POINTS;
}


namespace AndersonImpurityParams
{
    inline const double& UINT = RuntimeConfig::UINT;
    
    inline const double& MU = RuntimeConfig::MU;

    inline const std::string& DOS_TYPE = RuntimeConfig::DOS_TYPE;
    
    inline const int BATH_LEVEL_COUNT = 1;

    // DISCRETE density of states
    inline const std::array<double, BATH_LEVEL_COUNT> E_vals = {1.0, };
    inline const std::array<std::complex<double>, BATH_LEVEL_COUNT> V_vals = {1.0, };

    // BOX density of states (if chosen, see P. Chalupa et al 2022)
    inline const double& D = RuntimeConfig::D;

    inline const double& DELTA_0 = RuntimeConfig::DELTA_0;

    inline const std::map<std::string, double>& REFINE_AT_POINTS = RuntimeConfig::REFINE_AT_POINTS;
}


namespace AndersonImpurityHolsteinParams
{   
    inline const double& g0 = RuntimeConfig::g0;
    inline const double& OMEGA0 = RuntimeConfig::OMEGA0;

    // remaining params are read from SquareHubbardParams
    inline struct V_ph_proxy { operator double() const { return 2.0*RuntimeConfig::g0*RuntimeConfig::g0/RuntimeConfig::OMEGA0; } } V_ph;
}
 
