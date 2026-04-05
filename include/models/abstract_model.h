#pragma once

#include <real_space_lattice.h>
#include <momentum_space.h>
#include <delta_function_vector.h>
#include <lattice_point_group.h>
#include <form_factor_container.h>

#include <fftw3.h> /**< Needed for bubbles calculations */
#include <string>

#include <params_physical.h>

#include <params_technical.h>


struct special_points_and_paths_t
{
    // a collection of the indices of special points (e.g. symmetry points like gamma, M) and paths (e.g. Fermi-surface, some paths touching Gamma and M in the BZ ..etc) in the momentum space of the model. Should be useful for computing observables, or plotting the data of a run.
    std::vector<unsigned> points;
    std::vector<std::vector<unsigned> > paths;

    // names to identify the points/paths by; useful for outputting. Must not contain special characters other than "_"
    std::vector<std::string> points_names;
    std::vector<std::string> paths_names;
};

template <unsigned dim>
struct momentum_space_data_t
{
    momentum_space_t<dim> *ptr{nullptr};

    // volume weights of momentum space points
    std::vector<double> volume_weights;

    // Use: k_prime_idx = point_group_in_momentum_idx_space[g_idx][k_idx]. The result k_prime is k after being acted upon by g
    std::vector<std::vector<unsigned> > point_group_in_momentum_idx_space;

    // form factors in momentum idx space. To access form factor use e.g.: auto mom_space_ff = form_factors_in_momentum_idx_space[shell_idx][ff_idx]
    std::vector<std::vector<ff_in_mom_idx_space_t> > form_factors_in_momentum_idx_space;

    // Maps containing results of operations in momentum idx space

    // Use: minus_k_idx = negative_momentum_idx_map[k_idx]
    std::vector<unsigned> negative_momentum_idx_map;

    // Use: k1_plus_k2_idx = momentum_plus_momentum_idx_map[k1][k2]
    std::vector<std::vector<unsigned> > momentum_plus_momentum_idx_map;

    // results of computing sin(p) + i cos(p) with p in the FINE LATTICE momentum idx space
    std::vector<complex_vector_t<dim> > exp_ik_idx_map;

    // results of computing sin(pr) + i cos(pr) with p in the FINE LATTICE momentum idx space and r is in the real lattice
    std::vector<std::vector<std::complex<double> > > exp_ikr_idx_map;
};

template <unsigned dim>
struct real_lattice_data_t
{
    real_space_lattice_t<dim> *ptr{nullptr};
    lattice_point_group_t<dim> *point_group_ptr{nullptr};
    std::vector<coord_units_t<dim> > unit_points;
};

template <unsigned dim, typename ff_tuple_idx_t>
struct form_factors_data_t
{
    form_factor_container_t<dim> *container_ptr{nullptr}; 

    // Use: auto [n_prime_idx, has_minus] = point_group_in_form_factor_idx_space[g_idx][n_idx]
    // The has_minus bool variable is true if the resulting form factor after mapping through g picks a minus.
    std::vector<std::vector<std::tuple<unsigned, bool> > > point_group_in_form_factor_idx_space;  

    std::map<ff_tuple_idx_t, unsigned> tuple_index_to_flat_index_map; 
    std::vector<ff_tuple_idx_t> flat_index_to_tuple_index_map;

    // conj(f)*f weights that appear in the fourier transform expression. Access by [R][m][n] with R real lattice index, m and n form factor indices
    std::vector<std::vector<std::vector<dcomplex> > > fourier_transform_weights;
    std::vector<unsigned > fourier_transform_weights_r_support_indices;
    std::vector<unsigned> real_space_to_real_fourier_transform_idx_map;

    // calculates  M_{nm}(k) = convolution(conj(f_n), f_m)(k). Used for writing the bare vertex in a given channel  in ff space; V^0_{m, n}(q) = \int_{k-k'} M_{n, m}(k-k')V^0_{k,k'}(q). Valid if V^0 depends on k-k' instead of k and k' individually, which will be the case for density density interactions.
    std::vector<std::vector<std::vector<dcomplex> > > ff_convolution; 
};

struct fftw_data_t
{
    fftw_complex *input, *output;        /**< input output objects for the FFTW */
    fftw_plan backward_plan, forward_plan;                /**< plan for FFTW. p: plan Gvec and Svec-> FFT (BACKWARD). p_b: plan bubble -> FFT (FORWARD) */   

};

// since all methods are static, this class definition only serves as a template/tutorial of what a complete model must defined
template<unsigned dim_>
class AbstractModel
{
 public: 
    static const unsigned dim = dim_;

    static const std::string GetName(); // name of the model
    
    // the parameter name-value pairs returned here will be included in the output file
    static const std::vector<std::pair<std::string, double>> GetParamNameValuePairs();
    
    // a model should have all these methods
    static const unsigned GetQuantumNumbersCount();
    static const unsigned GetCoarseMomentaCount();
    static const unsigned GetRefinedMomentaCount();
    static const unsigned GetFineMomentaCount();
    static const unsigned GetFFTDim();
    static const unsigned GetMomentumFormFactorsCount();

    // a model should specify a dispersion (chemical potential included)
    static double E(const int idx_p, const int idx_w); 

    // dispersion as a function of p coordinate. Used only when adaptive integration over BZ is needed (e.g. calculation of static bubble)
    static double E_of_p_coord(const coord_t<dim> p_coord); 


    // a model should specify its bare vertex
    static double vertex_4pt_bare( 
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out,
	const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out,
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       ); // bare vertex/interaction of the model

    // a model should define its real space, momentum space, and an (optionally) fine grained momentum space (for bubble computation)
    static real_lattice_data_t<dim> &RealLattice();
    static momentum_space_data_t<dim> &MomentumGrid();
    static momentum_space_data_t<dim> &FineMomentumGrid();

    // a model should define an instance of this structure, which may be empty, if there are no special points nor paths
    static special_points_and_paths_t &SpecialPointsAndPaths();

    // a model should define the following accessor functions
    static coord_t<2> GetMomentumCoordFromCoarseIdx(const int idx_k);
    static coord_t<2> GetMomentumCoordFromFineIdx(const int idx_p);
    static unsigned &MomentaCountUnrefined(); // will be different from the size of the original momentum grid in case there additional refinement on top
    
    // a model should def its form factors and the following accessor functions
    static form_factors_data_t<dim, shell_ff_idx_t> &FormFactors();
    static delta_function_vector_t<dim> GetFormFactor(const unsigned n);
    static ff_in_mom_idx_space_t &GetFormFactorInMomentumIdxSpace(const unsigned n);
    static ff_in_mom_idx_space_t &GetFormFactorInFineMomentumIdxSpace(const unsigned n);

  
    // fourier transform library (FFTW) related objects
    static fftw_data_t &FFTW(unsigned thread_idx);
    static std::vector<fftw_data_t> &FFTWs();

    	
    // a model should define maps between the coarse and the fine momentum spaces (if given)
    static std::vector<unsigned> &CoarseToFineIdxMap();
    static std::vector<unsigned> &FineToCoarseIdxMap();
};
