#include <lattice2D.h>
#include <delta_function_vector.h>

class irreps_form_factors_t
{
 public:
    irreps_form_factors_t();

    delta_function_vector_t get(unsigned irrep_idx, unsigned shell_idx, unsigned bond_idx);
    
    delta_function_vector_t project_delta_function(delta_function_vector_t delta_func, unsigned mu);

    std::vector<std::vector<delta_function_vector_t> > make_form_factors_from_shells(std::map<double, std::vector<coord_t> > shells, unsigned mu);

    double compare_idx_space_form_factors(std::vector<std::complex<double> > ff1, std::vector<std::complex<double> > ff2);

    std::vector<std::tuple<unsigned, unsigned, unsigned> > remove_vanishing_and_repeated_form_factors();

    std::vector<std::tuple<unsigned, unsigned, unsigned> > get_point_group_in_form_factor_idx_space();
};
