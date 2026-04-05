# config.mk -- fRG feature flags
# Toggle features by commenting / uncommenting lines.
# Included by the main makefile.

# -- Method (pick one) --------------------------------------------------------
#CFLAGS += -D FLOW_EQUATION_METHOD
CFLAGS += -D SELF_CONSISTENT_METHOD


# -- Model selection ----------------------------------------------------------
# Choose from: HubbardAtom, AndersonImpurity, AndersonImpurityHolstein,
#              SquareHubbard, SquareHubbardHolstein, SquareHubbardPeierls,
#              SquareHubbardLongRange, TriangularHubbard, ChainHubbard, FCCHubbard
CFLAGS += -D THE_MODEL_VAL=HubbardAtom


# -- Output path settings ------------------------------------------------------
CFLAGS += -D OUTPUT_DIRECTORY_VAL=\"./\"  # the default value

# -- Grid / frequency / form-facItor resolution --------------------------------
CFLAGS += -D COUNT_VAL=5  
CFLAGS += -D K_DIM_VAL=16
CFLAGS += -D P_IN_K_VAL=5
CFLAGS += -D FORMFACTOR_SHELL_COUNT_VAL=1.0  # For models without spatial structure (e.g. Anderson impurity, Hubbard atom), this must be 1.0

# -- G cutoff scheme (pick one) -------------------------------------------------
#CFLAGS += -D INT_FLOW               # Interaction flow
#CFLAGS += -D EBERL_FLOW             # Eberlein flow (WARNING: implement G_latt correctly)
CFLAGS += -D OMEGA_FLOW              # Omega flow -- Salmhofer/Gierig (arXiv:1208.6131)
#CFLAGS += -D TEMP_FLOW              # Temperature flow (DOI:10.1103/PhysRevB.64.184516); rescale V after flow
#CFLAGS += -D INTRP_FLOW        
#CFLAGS += -D INVERSE_INTRP_FLOW
#CFLAGS += -D INVERSE_INTRP_FLOW_OMEGA

# -- Enable for model specific U-flow -------------------------------------------------
#CFLAGS += -D ACTUAL_U_FLOW          # bare interaction is interpolatied simply U(t) = tU, t goes from 0 to 1 
#CFLAGS += -D FREE_U_FLOW            # bare interaction is flowing in a model dependent fashion


# -- Output & I/O -------------------------------------------------------------
CFLAGS += -D MULTIFILE_OUTPUT     # if disabled, only final of the flow/iteration is outputted
#CFLAGS += -D AUTORESUME_CALCULATION  # Resume from existing files if present
#CFLAGS += -D OUT_UPDATE              # Write data updatenum times during flow
#CFLAGS += -D POSTPROCESSING_EVERY_STEP  # Compute post-processing quantities and fluctuation diagnostics at every step (not recommended)

# -- Self-energy (pick one) --------------------------------------------------------------
# If MULTILOOP is used, SELFEN_FLOW (standard flow) can be switched with one of the SDE options;
# the right-hand side of the flow equation will be initialized to the standard one.
CFLAGS += -D SELFEN_FLOW                        # Loop-ordered fRG self-energy flow
#CFLAGS += -D SELFEN_SDE_CONVENTIONAL_MAGNETIC   # SDE (conventional, magnetic channel)
#CFLAGS += -D SELFEN_SDE_CONVENTIONAL_DENSITY    # SDE (conventional, density channel)
#CFLAGS += -D SELFEN_SDE_SBE_MAGNETIC             # SDE via SBE formalism, magnetic channel
#CFLAGS += -D SELFEN_SDE_SBE_DENSITY             # SDE via SBE formalism, density channel
#CFLAGS += -D SELFEN_SDE_SBE_SUPERCONDUCTING     # SDE via SBE formalism, SC channel

#CFLAGS += -D FIX_FILLING             # Fix filling to initial value (set by chemical potential)
#CFLAGS += -D PIN_SELFENERGY_TAIL_TO_ZERO # Pin the high-frequency tail of the self-energy to zero (often not recommended)


# -- Multiloop settings ----------------------------------------
# vertex
#CFLAGS += -D MULTILOOP               # Multiloop corrections l > 2; set loop_number in const.h
CFLAGS += -D EXTRA_LOOP_NUM_VAL=0    # Number of additional loops beyond 1-loop; 0 only includes the Katanin correction
CFLAGS += -D ERR_LOOP_ABS_VAL=1e-5   # Absolute error threshold for multiloop convergence (per loop)
CFLAGS += -D ERR_LOOP_REL_VAL=1e-4   # Relative error threshold for multiloop convergence (per loop)
#CFLAGS += -D FORCE_CALCULATE_ALL_MULTILOOP_CORRECTIONS # Calculate all multiloop corrections up to EXTRA_LOOP_NUM_VAL, even if convergence is reached earlier (not recommended)
#CFLAGS += -D STORE_MULTILOOP_VERTEX_CORRECTIONS # Store multiloop corrections separately in the .hdf5 output

CFLAGS += -D SELFENERGY_ITERATIONS_MAX_VAL=100 # If multiloop is on and an SDE self-energy is used: Maximum number of iterations for self-consistent SDE solution (if SELFEN_SDE is used)
CFLAGS += -D ABS_ERROR_SELFENERGY_ITERATIONS_VAL=1e-4 # Absolute convergence threshold for self-consistent SDE iterations
#CFLAGS += -D FORCE_CALCULATE_ALL_SELFENERGY_ITERATIONS # Perform the maximum number of iterations for the self-consistent SDE solution, even if convergence is reached earlier (not recommended)
#CFLAGS += -D NO_KATANIN             # Perform multiloop without Katanin correction


# -- Fluctuation diagnostics ----------------------------------------
CFLAGS += -D SPLIT_SUSC_CONTRIBUTIONS # Split the contributions to susceptibilities into different channels in postprocessing and output them separately in the .hdf5 file
CFLAGS += -D SPLIT_LAMBDA_CONTRIBUTIONS # Split the contributions to the Hedin vertex into different channels in postprocessing and output them separately in the .hdf5 file
#CFLAGS += -D CALCULATE_SIG_SDE_INTEGRANDS # Calculate the integrands of the SDE for the self-energy during the flow and output them in the .hdf5 file for fluctuation diagnostics



# -- SBE related ----------------------------------------
#CFLAGS += -D F_NONZERO              # Use extended parameters
#CFLAGS += -D SBEa_APPROXIMATION  # Enable SBEa approximation (disables rest-function flow)
#CFLAGS += -D NO_HEDIN_VERTEX_FLOW   # Neglect the flow of the Hedin vertex (set it to the bare value, i.e. 1)
#CFLAGS += -D SBEb_APPROXIMATION	 # Neglect the rest function flow and remove its feedback completely
#CFLAGS += -D BOSONISE_M             # Full bosonic propagator / Yukawa coupling implementation
#CFLAGS += -D BOSONISE_M_LAZY        # Lazy: compute necessary entries of M for the approximate reconstruction
#CFLAGS += -D BARE_INTERACTION_NOT_DENSITY_DENSITY # If bare interaction is not purely density-density, the double-counting terms in the B+F splitting are more complicated (e.g. in the presence of SSH phonons)

# -- Performance & numerics ---------------------------------------------------
#CFLAGS += -D PRECOMPUTE_STATE_PROJECTIONS # Precompute the state projections (form-factor integrals) for all combinations of external momenta; can lead to high memory usage 
#CFLAGS += -D PRECOMPUTE_MOMENTUM_ARITHMETIC  # Fast index arithmetic; high memory cost for large grids
#CFLAGS += -D ITERATE_DYSON_EQUATION_FOR_BOSONIC_SC_PROPAGATOR # Instead of w = B/(1-P*B), iterate over w = B + BPw (recommended for strong coupling and/or near instabilities)
#CFLAGS += -D ITERATE_DYSON_EQUATION_FOR_BOSONIC_D_PROPAGATOR
#CFLAGS += -D ITERATE_DYSON_EQUATION_FOR_BOSONIC_M_PROPAGATOR


# -- General approximations --------------------------------------------------------
#CFLAGS += -D SET_MIXED_BUBBLES_TO_ZERO           # TODO: make work for multiloop / self-consistent SBE
#CFLAGS += -D SET_OFFDIAGONAL_M_IN_FROMFACTOR_SPACE_TO_ZERO
CFLAGS += -D RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX  # Recommended with PRECOMPUTE_STATE_PROJECTIONS


# -- Form factors & symmetries ------------------------------------------------
CFLAGS += -D UTILIZE_ALGEBRAIC_SYMMETRIES  # Use algebraic (time-reversal, crossing  ..etc symmetry) symmetries (highly recommended)
CFLAGS += -D UTILIZE_LATTICE_SYMMETRIES    # Use lattice point-group symmetries (highly recommended)


# -- ODE solver settings -------------------------------------------------------
CFLAGS += -D ERR_ABS_VAL=1e-5  # Absolute error threshold for ODE solver (per step)
CFLAGS += -D ERR_REL_VAL=1e-4 # Relative error threshold for ODE solver (per step)
CFLAGS += -D MAX_COUPLING_VAL=1e+4 # Divergence criterion:If the maximum coupling exceeds this value during the flow, the flow will be stopped (to prevent divergence and numerical instability)


# -- Debug & diagnostics ------------------------------------------------------
#CFLAGS += -D DEBUG_EULER            # Euler integration over lambda with fixed step size
CFLAGS += -D NSTEPS_VAL=100          # Used only when DEBUG_EULER is enabled
#CFLAGS += -D STATIC_CALCULATION     # fRG with static vertices (freq-independent Sigma and V; Sigma not fully correct)
CFLAGS += -D MICKY_MOUSE            # Sets very low technical parameters, for test runs


# -- Experimental / WIP -------------------------------------------------------
#CFLAGS += -D SU2_UNRESOLVED         # WIP: full form of fRG equations
#CFLAGS += -D INCLUDE_FREE_ENERGY
#CFLAGS += -D LOGARITHMIC_GRIDS     # Logarithmic frequency grids
                                    # Default (no flag): uniform grids

