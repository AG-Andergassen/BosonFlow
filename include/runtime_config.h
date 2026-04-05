#pragma once

#include <map>
#include <string>

namespace RuntimeConfig {
    // Global physical parameters accessible at runtime
    // These are shared across all model namespaces
    extern double BETA;          ///< Inverse temperature (all models)
    extern double UINT;          ///< On-site interaction (all models that use it)
    extern double VINT;          ///< NN interaction (SquareHubbard, TriangularHubbard, ChainHubbard)
    extern double MU;            ///< Chemical potential (all models that use it)
    extern double T_PRIME;       ///< NNN hopping (SquareHubbard, TriangularHubbard, ChainHubbard, FCCHubbard)
    extern double T_PRIME_PRIME; ///< NNNN hopping (TriangularHubbard, FCCHubbard)
    extern double WINT;          ///< NNN interaction (TriangularHubbard)
    extern double WWINT;         ///< NNNN interaction (TriangularHubbard)
    
    // SquareHubbardLongRangeParams
    extern double Q;             ///< Long-range coupling
    extern double K_s;           ///< Thomas-Fermi screening
    
    // Holstein/Peierls/andersonImpurityHolstein phonon params
    extern double g0;            ///< e-p coupling (Holstein variants)
    extern double OMEGA0;        ///< Phonon frequency (Holstein variants)
    
    // AndersonImpurityParams specifics
    extern double D;             ///< Half-bandwidth
    extern double DELTA_0;       ///< Hybridisation strength
    extern std::string DOS_TYPE; ///< Anderson impurity DOS type: BOX or CONST
    extern std::string OUTPUT_DIRECTORY; ///< Base directory for output files (set via positional CLI arg)
    extern bool HAS_FILLING_OVERRIDE; ///< If true, use FILLING_OVERRIDE as target filling
    extern double FILLING_OVERRIDE;   ///< Target filling set via --filling

    // Runtime map of refinement points, populated via --refine-at
    // Example: {"idx_00": 0.02, "idx_pipi": 0.01}
    extern std::map<std::string, double> REFINE_AT_POINTS;
    
    // Parse command line arguments to override defaults
    void ParseCommandLine(int argc, char** argv);
    
    // Print usage information
    void PrintUsage();
}
