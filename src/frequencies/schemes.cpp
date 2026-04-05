#include <frequencies/schemes.h>
#include <frequencies/matsubara_space.h> // for w_val and W_val
#include <params_technical.h>
#include <models/concrete_available_models.h>
#include <boost/math/interpolators/pchip.hpp>
#include <mymath.h>
#include <gf_interpolator.h>
//#include <frequencies/ir_basis.h>

template <typename Model>
void FrequencyDependenceScheme<Model>::Init(const FrequencyDependenceSchemeName freq_scheme)
{
    ChosenSchemeName() = freq_scheme;

    InitFrequencyGrids();
    
    if (ChosenSchemeName() == FrequencyDependenceSchemeName::UniformGrids){
	//...
    }
    else if (ChosenSchemeName() == FrequencyDependenceSchemeName::IRBasis){
	//...
    }
    else {
	//...
    } 

    IntegrationWeightVecs1DPtr().reset(new gf<double, 1>(boost::extents[ffreq(PositiveIntegrationRange())]));
    generate_weights( IntegrationWeightVecs1D(), PositiveIntegrationRange()-FrequenciesCount::TAIL_LENGTH, FrequenciesCount::TAIL_LENGTH, FrequenciesCount::FIT_ORDER );     /**< to correct finite Matsubara freq summation.*/ 

    IntegrationWeightVecs2DPtr().reset(new gf<double, 2>(boost::extents[ffreq(PositiveIntegrationRange())][ffreq(PositiveIntegrationRange())]) );    
    generate_2d_weights( IntegrationWeightVecs2D(), PositiveIntegrationRange()-FrequenciesCount::TAIL_LENGTH, FrequenciesCount::TAIL_LENGTH, FrequenciesCount::FIT_ORDER ); /**< to correct double finite Matsubara freq summation.*/

}


template <typename Model>
void FrequencyDependenceScheme<Model>::InitFrequencyGrids()
{
    InitFrequencyGridData(G_FermionicGridRef(), FrequenciesCount::POS_1P_RANGE);
   
    InitFrequencyGridData(Sig_FermionicGridRef(), FrequenciesCount::Sig::POS_FERM);

    InitFrequencyGridData(w_BosonicGridRef(), FrequenciesCount::w::POS_BOS);

    InitFrequencyGridData(lambda_FermionicGridRef(), FrequenciesCount::lambda::POS_FERM);
    InitFrequencyGridData(lambda_BosonicGridRef(), FrequenciesCount::lambda::POS_BOS);

    InitFrequencyGridData(M_FermionicGridRef(), FrequenciesCount::M::POS_FERM);
    InitFrequencyGridData(M_BosonicGridRef(), FrequenciesCount::M::POS_BOS);

    InitFrequencyGridData(BubbleTU_FermionicGridRef(), FrequenciesCount::bubble::POS_FERM);
    InitFrequencyGridData(BubbleTU_BosonicGridRef(), FrequenciesCount::bubble::POS_BOS);

    InitFrequencyGridData(Projected_nabla_FermionicGridRef(), FrequenciesCount::Projected_nabla::POS_FERM);
    InitFrequencyGridData(Projected_nabla_BosonicGridRef(), FrequenciesCount::Projected_nabla::POS_BOS);


    InitFrequencyGridData(Projected_M_FermionicGridRef(), FrequenciesCount::Projected_M::POS_FERM);
    InitFrequencyGridData(Projected_M_BosonicGridRef(), FrequenciesCount::Projected_M::POS_BOS);
    
    if (ChosenSchemeName() == FrequencyDependenceSchemeName::LogarithmicGrids)
	PositiveIntegrationRangeRef() = FrequenciesCount::GRID_MULTIPLIER*FrequenciesCount::Integration::POS_RANGE;
    else
	PositiveIntegrationRangeRef() = FrequenciesCount::Integration::POS_RANGE;

} 

template <typename Model>
void FrequencyDependenceScheme<Model>::ConstructSamplingNeighbors(frequency_grid_data_t &frequency_grid)
{
    if (frequency_grid.sampling_freqs.size() < 4){
	throw("Number of frequencies is too low to perform interpolation!");
    }
    for (int freq = frequency_grid.sampling_freqs[0]; freq < frequency_grid.sampling_freqs.back(); freq ++)
	for (int j = 0; j < frequency_grid.sampling_freqs.size(); j++){
	    int sampling_freq_index = frequency_grid.sampling_freqs.size() - 1 - j;
	    int sampling_freq = frequency_grid.sampling_freqs[sampling_freq_index];
	    if ( sampling_freq <= freq ){
		frequency_grid.freqs_left_sampling_neighbors.push_back(sampling_freq);
		    
		if (sampling_freq < 0){
		    int offset = (sampling_freq_index - 1) < 0 ? 1 : 0;
		    frequency_grid.freqs_4_sampling_neighbors.push_back(
                               {frequency_grid.sampling_freqs[sampling_freq_index-1+offset], 
				frequency_grid.sampling_freqs[sampling_freq_index+0+offset], 
				frequency_grid.sampling_freqs[sampling_freq_index+1+offset],
				frequency_grid.sampling_freqs[sampling_freq_index+2+offset] }
					     );
		}else{
		    int offset = (sampling_freq_index + 1) >= frequency_grid.sampling_freqs.size() ? 1 : 0;
		    frequency_grid.freqs_4_sampling_neighbors.push_back(
                               {frequency_grid.sampling_freqs[sampling_freq_index-2+offset], 
				frequency_grid.sampling_freqs[sampling_freq_index-1+offset], 
				frequency_grid.sampling_freqs[sampling_freq_index+0+offset],
				frequency_grid.sampling_freqs[sampling_freq_index+1+offset] }
					     );
		    
		}
	     
		break;
	    }
	}
}


template <typename Model>
void FrequencyDependenceScheme<Model>::InitFrequencyGridData(frequency_grid_data_t &frequency_grid, int positive_freqs_count)
{
    if (ChosenSchemeName() == FrequencyDependenceSchemeName::LogarithmicGrids){
	frequency_grid.positive_freqs_count = FrequenciesCount::GRID_MULTIPLIER*positive_freqs_count;
    }else if (ChosenSchemeName() == FrequencyDependenceSchemeName::UniformGrids){
	frequency_grid.positive_freqs_count = positive_freqs_count;
    }else{
	//....
	frequency_grid.positive_freqs_count = positive_freqs_count;
    }

    frequency_grid.freqs_count = 2 * frequency_grid.positive_freqs_count;

    // if bosonic, add one for the 0th matsubara frequency
    if (frequency_grid.statistics == Statistics::BOSE) 
	frequency_grid.freqs_count += 1;

    if (ChosenSchemeName() == FrequencyDependenceSchemeName::LogarithmicGrids){
	if (frequency_grid.statistics == Statistics::FERMI) // to avoid going out of bounds
	    frequency_grid.sampling_freqs = GenerateLogSpacedFrequencies(frequency_grid.positive_freqs_count - 1, 
									 frequency_grid.positive_freqs_count/FrequenciesCount::GRID_MULTIPLIER, 
									 frequency_grid.statistics);
	else
	    frequency_grid.sampling_freqs = GenerateLogSpacedFrequencies(frequency_grid.positive_freqs_count, 
									 frequency_grid.positive_freqs_count/FrequenciesCount::GRID_MULTIPLIER, 
									 frequency_grid.statistics);
	ConstructSamplingNeighbors(frequency_grid);
    }else{ // if other schemes
	for (int i = -frequency_grid.positive_freqs_count; i < frequency_grid.freqs_count; i ++)
	    frequency_grid.sampling_freqs.push_back(i);
    }
}

template <typename Model>
std::vector<int> FrequencyDependenceScheme<Model>::GenerateLogSpacedFrequencies(int upper_bound, int number_of_positive_freqs, Statistics zeta)
{
    // Define the range [a, b]
    int a = 1; // Lower bound
    int b = upper_bound; // Upper bound

    // Number of logarithmically spaced integers to generate
    int N;
    if (zeta == Statistics::BOSE){
	N = number_of_positive_freqs;
    }else{
	N = number_of_positive_freqs - 1; 
    }

    // Calculate the logarithmic spacing factor
    double factor = std::pow(static_cast<double>(b) / a, 1.0 / (N - 1));

    // Initialize the vector to store the generated integers
    std::vector<int> logSpacedIntegers;

    // Generate N logarithmically spaced integers
    for (int i = 0; i < N; i++) {
	int value = static_cast<int>(a * std::pow(factor, i));
	if (logSpacedIntegers.size() == 0 || logSpacedIntegers.back() < value)
	    logSpacedIntegers.push_back(value);
	else
	    logSpacedIntegers.push_back(logSpacedIntegers.back()+1);
    }

    std::vector<int> logspaced_frequencies;
    if (zeta == Statistics::BOSE){
	for (int i = 0; i < logSpacedIntegers.size(); i ++) {
	    logspaced_frequencies.push_back(-logSpacedIntegers[logSpacedIntegers.size() - 1 - i]);
	}
	logspaced_frequencies.push_back(0);
	for (int i = 0; i < logSpacedIntegers.size(); i ++) {
	    logspaced_frequencies.push_back(logSpacedIntegers[i]);
	}
    } else { // if zeta == Statistics::FERMI
	for (int i = 0; i < logSpacedIntegers.size(); i ++) {
	    logspaced_frequencies.push_back(-logSpacedIntegers[logSpacedIntegers.size() - 1 - i]-1);
	}
	logspaced_frequencies.push_back(-1);
	logspaced_frequencies.push_back(0);
	for (int i = 0; i < logSpacedIntegers.size(); i ++) {
	    logspaced_frequencies.push_back(logSpacedIntegers[i]);
	}
    }
    return logspaced_frequencies;
}

template <typename Model>
const frequency_ranges_sbe_t FrequencyDependenceScheme<Model>::GetFrequencyRanges()
{
    frequency_ranges_sbe_t ranges;
    ranges.G_fermionic_positive_freqs_count = G_FermionicGrid().positive_freqs_count;
    ranges.Sig_fermionic_positive_freqs_count = Sig_FermionicGrid().positive_freqs_count;
    ranges.w_bosonic_positive_freqs_count = w_BosonicGrid().positive_freqs_count;
    ranges.lambda_bosonic_positive_freqs_count = lambda_BosonicGrid().positive_freqs_count;
    ranges.lambda_fermionic_positive_freqs_count = lambda_FermionicGrid().positive_freqs_count;
    ranges.M_bosonic_positive_freqs_count = M_BosonicGrid().positive_freqs_count;
    ranges.M_fermionic_positive_freqs_count = M_FermionicGrid().positive_freqs_count;
    
    ranges.projected_M_bosonic_positive_freqs_count = Projected_M_BosonicGrid().positive_freqs_count;
    ranges.projected_M_fermionic_positive_freqs_count = Projected_M_FermionicGrid().positive_freqs_count;

    ranges.projected_nabla_bosonic_positive_freqs_count = Projected_nabla_BosonicGrid().positive_freqs_count;
    ranges.projected_nabla_fermionic_positive_freqs_count = Projected_nabla_FermionicGrid().positive_freqs_count;


    ranges.bubble_tu_fermionic_positive_freqs_count = BubbleTU_FermionicGrid().positive_freqs_count;
    ranges.bubble_tu_bosonic_positive_freqs_count = BubbleTU_BosonicGrid().positive_freqs_count;
    ranges.positive_integration_range = PositiveIntegrationRange();

    ranges.G_fermionic_freqs_count = G_FermionicGrid().freqs_count;
    ranges.Sig_fermionic_freqs_count = Sig_FermionicGrid().freqs_count;
    ranges.w_bosonic_freqs_count = w_BosonicGrid().freqs_count;
    ranges.lambda_bosonic_freqs_count = lambda_BosonicGrid().freqs_count;
    ranges.lambda_fermionic_freqs_count = lambda_FermionicGrid().freqs_count;
    ranges.M_bosonic_freqs_count = M_BosonicGrid().freqs_count;
    ranges.M_fermionic_freqs_count = M_FermionicGrid().freqs_count;

    ranges.projected_M_bosonic_freqs_count = Projected_M_BosonicGrid().freqs_count;
    ranges.projected_M_fermionic_freqs_count = Projected_M_FermionicGrid().freqs_count;

    ranges.projected_nabla_bosonic_freqs_count = Projected_nabla_BosonicGrid().freqs_count;
    ranges.projected_nabla_fermionic_freqs_count = Projected_nabla_FermionicGrid().freqs_count;

    ranges.bubble_tu_fermionic_freqs_count = BubbleTU_FermionicGrid().freqs_count;
    ranges.bubble_tu_bosonic_freqs_count = BubbleTU_BosonicGrid().freqs_count;


    return ranges;
}


template <typename Model>
bool FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(const int W, const int w, const int wp, const unsigned bosonic_positive_frequencies_count, const unsigned fermionic_positive_frequencies_count)
{
	const unsigned pos_bos_count = bosonic_positive_frequencies_count;
	const unsigned bos_count = 2 * pos_bos_count + 1;

	const unsigned pos_ferm_count = fermionic_positive_frequencies_count;
	const unsigned ferm_count = 2 * pos_ferm_count;


	return (unsigned)( W + pos_bos_count ) < ( bos_count ) &&
	    (unsigned)( w + pos_ferm_count ) < ( ferm_count - fold(W,2)) &&
	    (unsigned)( wp + pos_ferm_count ) < ( ferm_count -fold(W,2) );
};


// instantiate class for the particular models
#define Instantiate(MODEL) template class FrequencyDependenceScheme<MODEL>;
WITH_MODELS(Instantiate)
