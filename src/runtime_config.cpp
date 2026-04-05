#include <runtime_config.h>
#include <iostream>
#include <cstring>
#include <string>
#include <algorithm>
#include <cctype>

#define RUNTIMECONFIG_STRINGIFY_INNER(x) #x
#define RUNTIMECONFIG_STRINGIFY(x) RUNTIMECONFIG_STRINGIFY_INNER(x)

namespace RuntimeConfig {
    // Global physical parameters shared across all models
    // These represent the canonical values; individual model namespaces
    // reference these via const& to avoid copies.
    
    double BETA          = 1.0;
    double UINT          = 2.0;          // On-site interaction U
    double VINT          = 0.0;          // NN interaction V
    double MU            = 0.0;          // Chemical potential
    double T_PRIME       = 0.0;          // NNN hopping
    double T_PRIME_PRIME = 0.0;        // NNNN hopping (FCC default)
    double WINT          = 0.0;          // NNN interaction (triangular)
    double WWINT         = 0.0;          // NNNN interaction (triangular)
    double Q             = 1.0;   // Long-range coupling (2.0 / 0.00394666)
    double K_s           = 0.0;          // Thomas-Fermi screening
    double g0            = 1.0;  // e-p coupling (Holstein, default)
    double OMEGA0        = 1.0;       // Phonon frequency
    double D             = 10.0;         // Half-bandwidth (Anderson)
    double DELTA_0       = 0.63;         // Hybridisation strength (Anderson)
    std::string DOS_TYPE = "CONST";  // Anderson impurity DOS type
    std::string OUTPUT_DIRECTORY = "";
    std::map<std::string, double> REFINE_AT_POINTS = {};
    bool HAS_FILLING_OVERRIDE = false;
    double FILLING_OVERRIDE = 0.0;

    static bool IsFlag(const std::string& arg, const char* long_name) {
        const std::string double_dash = std::string("--") + long_name;
        const std::string single_dash = std::string("-") + long_name;
        return arg == double_dash || arg == single_dash;
    }

    static void FailArgument(const std::string& message) {
        std::cerr << "Error: " << message << "\n\n";
        PrintUsage();
        std::exit(1);
    }

    static double ParseDoubleOrFail(const std::string& flag, const std::string& token) {
        try {
            size_t pos = 0;
            const double value = std::stod(token, &pos);
            if (pos != token.size()) {
                FailArgument("Invalid value '" + token + "' for " + flag + ".");
            }
            return value;
        } catch (...) {
            FailArgument("Invalid value '" + token + "' for " + flag + ".");
        }
        return 0.0;
    }

    static std::string ParseDosTypeOrFail(const std::string& token) {
        std::string normalized = token;
        std::transform(normalized.begin(), normalized.end(), normalized.begin(),
                       [](unsigned char c){ return static_cast<char>(std::toupper(c)); });

        if (normalized != "BOX" && normalized != "CONST") {
            FailArgument("Invalid value '" + token + "' for --dos-type. Allowed values are BOX or CONST.");
        }
        return normalized;
    }
    
    void ParseCommandLine(int argc, char** argv) {
        bool output_dir_set = false;
        for (int i = 1; i < argc; ++i) {
            std::string arg(argv[i]);

            // First positional argument is interpreted as output directory.
            if (!arg.empty() && arg[0] != '-') {
                if (!output_dir_set) {
                    OUTPUT_DIRECTORY = arg;
                    output_dir_set = true;
                } else {
                    FailArgument("Unexpected extra positional argument '" + arg + "'. Only one output directory is allowed.");
                }
                continue;
            }

            if (arg == "--help" || arg == "-h") {
                PrintUsage();
                std::exit(0);
            }

            if (IsFlag(arg, "refine-at-clear")) {
                REFINE_AT_POINTS.clear();
                continue;
            }

            if (IsFlag(arg, "filling")) {
                if (i + 1 >= argc) {
                    FailArgument("--filling requires a numeric value.");
                }
                FILLING_OVERRIDE = ParseDoubleOrFail("--filling", argv[++i]);
                HAS_FILLING_OVERRIDE = true;
                continue;
            }

            if (IsFlag(arg, "dos-type")) {
                if (i + 1 >= argc) {
                    FailArgument("--dos-type requires a value (BOX or CONST).");
                }
                DOS_TYPE = ParseDosTypeOrFail(argv[++i]);
                continue;
            }

            if (IsFlag(arg, "refine-at")) {
                if (i + 1 >= argc) {
                    FailArgument("--refine-at requires an argument in the form '<idx_name>:<percent>'.");
                }

                std::string token(argv[++i]);
                size_t sep_pos = token.find(':');
                if (sep_pos == std::string::npos) {
                    FailArgument("Invalid --refine-at argument '" + token + "'. Use '<idx_name>:<percent>'.");
                }

                const std::string idx_name = token.substr(0, sep_pos);
                const std::string value_token = token.substr(sep_pos + 1);

                if (idx_name.empty()) {
                    FailArgument("--refine-at index name cannot be empty.");
                }

                const double percent = ParseDoubleOrFail("--refine-at", value_token);
                if (percent < 0.0 || percent > 1.0) {
                    FailArgument("--refine-at percent must be between 0 and 1 (inclusive). Got '" + value_token + "'.");
                }

                REFINE_AT_POINTS[idx_name] = percent;
                continue;
            }

            // Generic parameter parsing
            #define PARSE_PARAM(flag_name, var) \
                if (IsFlag(arg, flag_name)) { \
                    if (i + 1 >= argc) { \
                        FailArgument(std::string("Missing value for ") + arg + "."); \
                    } \
                    var = ParseDoubleOrFail(arg, argv[++i]); \
                    continue; \
                }

            PARSE_PARAM("beta", BETA)
            PARSE_PARAM("u", UINT)
            PARSE_PARAM("u-prime", VINT)
            PARSE_PARAM("mu", MU)
            PARSE_PARAM("t-prime", T_PRIME)
            PARSE_PARAM("t-prime-prime", T_PRIME_PRIME)
            PARSE_PARAM("wint", WINT)
            PARSE_PARAM("wwint", WWINT)
            PARSE_PARAM("q", Q)
            PARSE_PARAM("k-s", K_s)
            PARSE_PARAM("g0", g0)
            PARSE_PARAM("omega0", OMEGA0)
            PARSE_PARAM("d", D)
            PARSE_PARAM("delta-0", DELTA_0)

            #undef PARSE_PARAM

            if (!arg.empty() && arg[0] == '-') {
                FailArgument("Unknown parameter '" + arg + "'.");
            }
        }

        if (!output_dir_set) {
            FailArgument("Missing required positional argument <output-dir>.");
        }
    }
    
    void PrintUsage() {
    #ifdef THE_MODEL_VAL
        const char* chosen_model = RUNTIMECONFIG_STRINGIFY(THE_MODEL_VAL);
    #else
        const char* chosen_model = "HubbardAtom";
    #endif

        std::cout << "\n=== BosonFlow Runtime Configuration ===\n\n"
                  << "Active compile-time model: " << chosen_model << "\n\n"
                  << "Global:\n"
                  << "  <output-dir>                     Required first positional output directory\n"
                  << "  --beta <value>                   Inverse temperature (default: 5.0)\n"
                  << "  --filling <value>                Target filling for initial chemical-potential adjustment\n"
                  << "  --dos-type <BOX|CONST>           Anderson impurity DOS/hybridisation type (default: CONST; BOX is not compatible with T-flow)\n"
                  << "  --refine-at <idx>:<percent>      Add/override refinement point (percent in [0,1], repeatable)\n"
                  << "  --refine-at-clear                Clear all refinement points\n"
                  << "  --help, -h                       Show this help message\n\n";

        const std::string model = chosen_model;
        std::cout << model << " parameters:\n";

        if (model == "SquareHubbard") {
            std::cout << "  --u <value>                      On-site interaction U (default: 2.0)\n"
                      << "  --u-prime <value>                Nearest-neighbor interaction V (default: 0.0)\n"
                      << "  --mu <value>                     Chemical potential (default: 0.0)\n"
                      << "  --t-prime <value>                Next-nearest-neighbor hopping t' (default: 0.0)\n\n";
        } else if (model == "TriangularHubbard") {
            std::cout << "  --u <value>                      On-site interaction U\n"
                      << "  --u-prime <value>                Nearest-neighbor interaction V\n"
                      << "  --wint <value>                   Next-nearest-neighbor interaction W\n"
                      << "  --wwint <value>                  Next-next-nearest-neighbor interaction W'\n"
                      << "  --mu <value>                     Chemical potential\n"
                      << "  --t-prime <value>                Next-nearest-neighbor hopping t'\n"
                      << "  --t-prime-prime <value>          Next-next-nearest-neighbor hopping t''\n\n";
        } else if (model == "ChainHubbard") {
            std::cout << "  --u <value>                      On-site interaction U\n"
                      << "  --u-prime <value>                Nearest-neighbor interaction V\n"
                      << "  --mu <value>                     Chemical potential\n"
                      << "  --t-prime <value>                Next-nearest-neighbor hopping t'\n\n";
        } else if (model == "FCCHubbard") {
            std::cout << "  --u <value>                      On-site interaction U\n"
                      << "  --mu <value>                     Chemical potential\n"
                      << "  --t-prime <value>                Next-nearest-neighbor hopping t'\n"
                      << "  --t-prime-prime <value>          Next-next-nearest-neighbor hopping t'' (default: -0.13)\n\n";
        } else if (model == "HubbardAtom") {
            std::cout << "  --u <value>                      On-site interaction U\n"
                      << "  --mu <value>                     Chemical potential\n\n";
        } else if (model == "AndersonImpurity") {
            std::cout << "  --u <value>                      On-site interaction U\n"
                      << "  --mu <value>                     Chemical potential\n"
                      << "  --d <value>                      Half-bandwidth D (default: 10.0)\n"
                      << "  --delta-0 <value>                Hybridisation strength Delta_0 (default: 0.63)\n\n";
        } else if (model == "SquareHubbardLongRange") {
            std::cout << "  --u <value>                      On-site interaction U\n"
                      << "  --u-prime <value>                Nearest-neighbor interaction V\n"
                      << "  --mu <value>                     Chemical potential\n"
                      << "  --t-prime <value>                Next-nearest-neighbor hopping t'\n"
                      << "  --q <value>                      Long-range coupling Q (default: 507.11)\n"
                      << "  --k-s <value>                    Thomas-Fermi screening k_s (default: 0.0)\n\n";
        } else if (model == "SquareHubbardHolstein" || model == "SquareHubbardPeierls") {
            std::cout << "  --u <value>                      On-site interaction U\n"
                      << "  --u-prime <value>                Nearest-neighbor interaction V\n"
                      << "  --mu <value>                     Chemical potential\n"
                      << "  --t-prime <value>                Next-nearest-neighbor hopping t'\n"
                      << "  --g0 <value>                     Electron-phonon coupling g0 (default: 15.81)\n"
                      << "  --omega0 <value>                 Phonon frequency omega0 (default: 1000.0)\n\n";
        } else if (model == "AndersonImpurityHolstein") {
            std::cout << "  --u <value>                      On-site interaction U\n"
                      << "  --mu <value>                     Chemical potential\n"
                      << "  --d <value>                      Half-bandwidth D (if BOX DOS)\n"
                      << "  --delta-0 <value>                Hybridisation strength Delta_0\n"
                      << "  --g0 <value>                     Electron-phonon coupling g0\n"
                      << "  --omega0 <value>                 Phonon frequency omega0\n\n";
        } else {
            std::cout << "  --u <value>                      On-site interaction U\n"
                      << "  --mu <value>                     Chemical potential\n\n";
        }

        std::cout << "Note: physical parameters are shared runtime variables; names are repeated per model for readability.\n\n"
                  << "Examples:\n"
                  << "  ./bin/run ./my/output/directory --beta 10.0 --u 3.0 --mu 0.2\n"
                  << "  ./bin/run ./results --beta 5.0 --u 2.0 --u-prime 1.5 --mu 0.5\n"
                  << "  ./bin/run ./results --beta 8.0 --u 2.0 --refine-at idx_00:0.02 --refine-at idx_pipi:0.015\n"
                  << "  ./bin/run --help\n\n";
    }
}
