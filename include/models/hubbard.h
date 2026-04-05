#pragma once

#include <models/abstract_model.h>

template<unsigned dim_>
class Hubbard : public AbstractModel<dim_>
{
 public:
    Hubbard() = delete;

    static inline const unsigned GetQuantumNumbersCount(){return 1/*QN_COUNT*/;};
    static inline const unsigned GetCoarseMomentaCount(){return MomentaCountUnrefined();};
    static inline const unsigned GetRefinedMomentaCount(){return MomentumGrid().ptr->m_points.size();};
    static inline const unsigned GetFineMomentaCount(){return FineMomentumGrid().ptr->m_points.size();};
    static inline const unsigned GetFFTDim(){return static_cast<unsigned>(std::sqrt(FineMomentumGrid().ptr->m_points.size()));};
    static inline const unsigned GetMomentumFormFactorsCount(){return FormFactors().flat_index_to_tuple_index_map.size();};
    

    static void Init(
		     basis_t<dim_> basis, 
		     std::vector<std::vector<matrix_t<dim_> > > point_group, 
		     std::map<std::string, coord_t<dim_> > special_points_coords, 
		     std::map<std::string, std::vector<coord_t<dim_> > > special_paths_coords, 
		     std::map<std::string, double> refine_at_points = {}, 
		     std::vector< delta_function_vector_t<dim_> > extra_form_factors = {} 
		     );

    static real_lattice_data_t<dim_> &RealLattice()
    {
	static real_lattice_data_t<dim_> real_lattice;
	return real_lattice;
    }

    static momentum_space_data_t<dim_> &MomentumGrid()
    {
	static momentum_space_data_t<dim_> momentum_grid;
	return momentum_grid;
    }

    static momentum_space_data_t<dim_> &FineMomentumGrid()
    {
	static momentum_space_data_t<dim_> fine_momentum_grid;
	return fine_momentum_grid;
    }


    static form_factors_data_t<dim_, shell_ff_idx_t> &FormFactors()
    {
	static form_factors_data_t<dim_, shell_ff_idx_t> form_factors;
	return form_factors;
    }

    static unsigned &MomentaCountUnrefined()
    {
	static unsigned momenta_count_unrefined;
	return momenta_count_unrefined;
    }

    static std::vector<unsigned> &CoarseToFineIdxMap()
    {
	static std::vector<unsigned> coarse_to_fine_idx_map;
	return coarse_to_fine_idx_map;
    }

    static std::vector<unsigned> &FineToCoarseIdxMap()
    {
	static std::vector<unsigned> fine_to_coarse_idx_map;
	return fine_to_coarse_idx_map;
    }

    // fourier transform library (FFTW) related objects
    static fftw_data_t &FFTW(unsigned thread_idx)
    {
        return FFTWs()[thread_idx];
    }

    static std::vector<fftw_data_t> &FFTWs()
    {
	static std::vector<fftw_data_t> fftw_data;
	return fftw_data;
    }

    static special_points_and_paths_t &SpecialPointsAndPaths()
    {
	static special_points_and_paths_t special_points_and_paths;
	return special_points_and_paths;
    }

    static coord_t<dim_> GetMomentumCoordFromCoarseIdx(const int idx_k);

    static coord_t<dim_> GetMomentumCoordFromFineIdx(const int idx_p);


    static delta_function_vector_t<dim_> GetFormFactor(const unsigned n);

    static ff_in_mom_idx_space_t &GetFormFactorInMomentumIdxSpace(const unsigned n);
    
    static ff_in_mom_idx_space_t &GetFormFactorInFineMomentumIdxSpace(const unsigned n);

  static unsigned SumFineMomentaIdxes(const unsigned p1_idx, const unsigned p2_idx);

  static unsigned GetNegativeFineMomentumIdx(const unsigned p_idx);
  
  static unsigned GetNegativeCoarseMomentumIdx(const unsigned k_idx);

  static unsigned GetFineMomentumIdxFromCoarse(const unsigned p_idx);

  static unsigned GetCoarseMomentumIdxFromFine(const unsigned p_idx);
  
 private:
    static void InitLattices(basis_t<dim_> basis);
    static void InitLatticesPointGroup(std::vector<std::vector<matrix_t<dim_> > > point_group);
    
    static void InitFormFactors(std::vector< delta_function_vector_t<dim_> > extra_form_factors);
    
    static double InitMomentumSymmetryMaps();

    static void InitFormFactorsSymmetryMaps();

    static void InitFormFactorsFourierTransformWeights();

    // creates the coarse-to-fine momentum index map, and returns the maximum error/deviation that was found
    static double InitCoarseToFineMomentumIdxMap();

    static double InitFineToCoarseMomentumIdxMap();
    
    static double InitPrecalculatedMomentumIdxArithmeticMaps();

    static void InitPrecalculatedTrigonometricMaps();

    static void InitSpecialPointsAndPaths(std::map<std::string, coord_t<dim_> > special_points_coords, std::map<std::string, std::vector<coord_t<dim_> > > special_paths_coords);

    static void RefineLattices(std::map<std::string, coord_t<dim_> > special_points_coords, std::map<std::string, std::vector<coord_t<dim_> > > special_paths_coords, std::map<std::string, double> refine_at_points);


};

