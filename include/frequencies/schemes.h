

/*******************************************************************************************//** @file
 *  		
 * 	file: 		ir_basis.h
 * 	contents:  	Definition of projection-related quantities 
 * 
****************************************************************************************************/
#pragma once

#include <vector>
#include <gf.h>
#include <memory>

enum class FrequencyDependenceSchemeName{UniformGrids, LogarithmicGrids, IRBasis};

enum class Statistics{BOSE = 0, FERMI = 1};

struct frequency_grid_data_t
{
    frequency_grid_data_t(Statistics statistics_):statistics(statistics_){}

    Statistics statistics;
    int positive_freqs_count;
    int freqs_count;
    
    std::vector<int> sampling_freqs;
    std::vector<int> freqs_left_sampling_neighbors;
    std::vector<std::vector<int> > freqs_4_sampling_neighbors; 
};

struct frequency_ranges_sbe_t
{
    unsigned G_fermionic_positive_freqs_count;
    unsigned Sig_fermionic_positive_freqs_count;
    unsigned w_bosonic_positive_freqs_count;
    unsigned lambda_bosonic_positive_freqs_count;
    unsigned lambda_fermionic_positive_freqs_count;
    unsigned M_bosonic_positive_freqs_count;
    unsigned M_fermionic_positive_freqs_count;
    unsigned bubble_tu_fermionic_positive_freqs_count;
    unsigned bubble_tu_bosonic_positive_freqs_count;
    unsigned positive_integration_range;
    unsigned projected_nabla_fermionic_positive_freqs_count;
    unsigned projected_nabla_bosonic_positive_freqs_count;
    unsigned projected_M_fermionic_positive_freqs_count;
    unsigned projected_M_bosonic_positive_freqs_count;

    unsigned G_fermionic_freqs_count;
    unsigned Sig_fermionic_freqs_count;
    unsigned w_bosonic_freqs_count;
    unsigned lambda_bosonic_freqs_count;
    unsigned lambda_fermionic_freqs_count;
    unsigned M_bosonic_freqs_count;
    unsigned M_fermionic_freqs_count;
    unsigned bubble_tu_fermionic_freqs_count;
    unsigned bubble_tu_bosonic_freqs_count;

    unsigned projected_nabla_fermionic_freqs_count;
    unsigned projected_nabla_bosonic_freqs_count;
    unsigned projected_M_fermionic_freqs_count;
    unsigned projected_M_bosonic_freqs_count;


};

template <typename Model> 
class FrequencyDependenceScheme
{
 public:
    static void Init(FrequencyDependenceSchemeName freq_scheme);

    static void InitSamplingFrequencies();

    static void InitFrequencyGrids();

    static const frequency_ranges_sbe_t GetFrequencyRanges();

    static FrequencyDependenceSchemeName &ChosenSchemeName(){
	static FrequencyDependenceSchemeName chosen_scheme_name;
	return chosen_scheme_name;
    }
    
    static frequency_grid_data_t &G_FermionicGridRef(){
	static frequency_grid_data_t G_fermionic_freqs(Statistics::FERMI);
	return G_fermionic_freqs;
    }

    static const frequency_grid_data_t &G_FermionicGrid(){
	return G_FermionicGridRef();
    }

    static frequency_grid_data_t &Sig_FermionicGridRef(){
	static frequency_grid_data_t Sig_fermionic_freqs(Statistics::FERMI);
	return Sig_fermionic_freqs;
    }
    static const frequency_grid_data_t &Sig_FermionicGrid(){
	return Sig_FermionicGridRef();
    }

    static frequency_grid_data_t &w_BosonicGridRef(){
	static frequency_grid_data_t w_bosonic_freqs(Statistics::BOSE);
	return w_bosonic_freqs;
    }
    static const frequency_grid_data_t &w_BosonicGrid(){
	return w_BosonicGridRef();
    }

    static frequency_grid_data_t &lambda_FermionicGridRef(){
	static frequency_grid_data_t lambda_fermionic_freqs(Statistics::FERMI);
	return lambda_fermionic_freqs;
    }
    static const frequency_grid_data_t &lambda_FermionicGrid(){
	return lambda_FermionicGridRef();
    }

    static frequency_grid_data_t &lambda_BosonicGridRef(){
	static frequency_grid_data_t lambda_bosonic_freqs(Statistics::BOSE);
	return lambda_bosonic_freqs;
    }
    static const frequency_grid_data_t &lambda_BosonicGrid(){
	return lambda_BosonicGridRef();
    }

    static frequency_grid_data_t &M_FermionicGridRef(){
	static frequency_grid_data_t M_fermionic_freqs(Statistics::FERMI);
	return M_fermionic_freqs;
    }
    static const frequency_grid_data_t &M_FermionicGrid(){
	return M_FermionicGridRef();
    } 

    static frequency_grid_data_t &M_BosonicGridRef(){
	static frequency_grid_data_t M_bosonic_freqs(Statistics::BOSE);
	return M_bosonic_freqs;
    }
    static const frequency_grid_data_t &M_BosonicGrid(){
	return M_BosonicGridRef();
    }

    static frequency_grid_data_t &Projected_nabla_FermionicGridRef(){
	static frequency_grid_data_t projected_nabla_fermionic_freqs(Statistics::FERMI);
	return projected_nabla_fermionic_freqs;
    }
    static const frequency_grid_data_t &Projected_nabla_FermionicGrid(){
	return Projected_nabla_FermionicGridRef();
    }

    static frequency_grid_data_t &Projected_nabla_BosonicGridRef(){
	static frequency_grid_data_t projected_nabla_bosonic_freqs(Statistics::BOSE);
	return projected_nabla_bosonic_freqs;
    }
    static const frequency_grid_data_t &Projected_nabla_BosonicGrid(){
	return Projected_nabla_BosonicGridRef();
    }

    static frequency_grid_data_t &Projected_M_FermionicGridRef(){
	static frequency_grid_data_t projected_M_fermionic_freqs(Statistics::FERMI);
	return projected_M_fermionic_freqs;
    }
    static const frequency_grid_data_t &Projected_M_FermionicGrid(){
	return Projected_M_FermionicGridRef();
    }

    static frequency_grid_data_t &Projected_M_BosonicGridRef(){
	static frequency_grid_data_t projected_M_bosonic_freqs(Statistics::BOSE);
	return projected_M_bosonic_freqs;
    }
    static const frequency_grid_data_t &Projected_M_BosonicGrid(){
	return Projected_M_BosonicGridRef();
    }

    static frequency_grid_data_t &BubbleTU_FermionicGridRef(){
	static frequency_grid_data_t bubble_fermionic_freqs(Statistics::FERMI);
	return bubble_fermionic_freqs;
    }
    static const frequency_grid_data_t &BubbleTU_FermionicGrid(){
	return BubbleTU_FermionicGridRef();
    }
 

    static frequency_grid_data_t &BubbleTU_BosonicGridRef(){
	static frequency_grid_data_t bubble_bosonic_freqs(Statistics::BOSE);
	return bubble_bosonic_freqs;
    }    
    static const frequency_grid_data_t &BubbleTU_BosonicGrid(){
	return BubbleTU_BosonicGridRef();
    }    

    static int &PositiveIntegrationRangeRef(){
	static int positive_integration_range;
	return positive_integration_range;
    }
    static const int PositiveIntegrationRange(){
        return PositiveIntegrationRangeRef();
    }
    
    

    static std::unique_ptr<gf< double, 1 > > &IntegrationWeightVecs1DPtr(){
	static std::unique_ptr<gf< double, 1 > > integration_weight_vecs_1d_ptr;
	return integration_weight_vecs_1d_ptr;
    }

    static gf< double, 1 > &IntegrationWeightVecs1D(){
	return *IntegrationWeightVecs1DPtr();
    }

    static std::unique_ptr<gf< double, 2 > > &IntegrationWeightVecs2DPtr(){
	static std::unique_ptr<gf< double, 2 > > integration_weight_vecs_2d_ptr;
	return integration_weight_vecs_2d_ptr;
    }
    
    static gf< double, 2 > &IntegrationWeightVecs2D(){
	return *IntegrationWeightVecs2DPtr();
    }


    static bool is_in_frequency_box_2_particle(const int W, const int w, const int wp, const unsigned bosonic_positive_frequencies_count, const unsigned fermionic_positive_frequencies_count);

 private:
    static void ConstructIRBasis();

    static void InitFrequencyGridData(frequency_grid_data_t &frequency_grid, int positive_freqs_count);

    static std::vector<int> GenerateLogSpacedFrequencies(int upper_bound, int number_of_positive_freqs, Statistics zeta);
    
    static void ConstructSamplingNeighbors(frequency_grid_data_t &frequency_grid);
    
};
