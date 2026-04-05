#pragma once

#include <primitive_zone.h>


using ff_in_mom_idx_space_t = std::vector< std::complex<double> > ;

using ff_idx_t = unsigned;
using shell_ff_idx_t = std::tuple<unsigned, ff_idx_t>;


template<unsigned dim>
class form_factor_container_t
{
 public: 
    form_factor_container_t(unsigned shells_count);

    ff_in_mom_idx_space_t & 
	get_in_momentum_idx_space(shell_ff_idx_t shell_ff_idx );

    ff_in_mom_idx_space_t &
	get_in_momentum_idx_space(unsigned shell_idx, unsigned ff_idx);
    
    // return delta_function_vector_t form factor corresponding to a shell and a bond index.
    delta_function_vector_t<dim> 
	get(unsigned shell_idx, unsigned ff_idx);

    delta_function_vector_t<dim> 
	get(shell_ff_idx_t shell_ff_idx );

    // get a list of all viable tuple indices (e.g. to iterate over)
    std::vector<shell_ff_idx_t> 
	get_indices();

        
    // taking a point group, and a set of form-factor shells, this method returns a map which instead of mapping momentum index values, will map from form-factor index values to a form-factor index values (upto a sign or complex conjugation specified by a string). 
    std::vector<std::map<shell_ff_idx_t, std::tuple<ff_idx_t, std::string> > > 
	make_point_group_in_shells_form_factor_idx_space(const std::vector<std::vector<unsigned> > &point_group_in_mom_idx_space, const std::vector<std::vector< ff_in_mom_idx_space_t > > &form_factors_shells_in_momentum_idx_space);

    // taking a point group and form-factors in one shell, this method returns a map which instead of mapping momentum index values, will map between the form-factor index values.
    std::vector<std::map<unsigned, std::tuple<unsigned, std::string> > > 
	with_shell_make_point_group_in_form_factor_idx_space(const std::vector<ff_in_mom_idx_space_t > &form_factors_in_momentum_idx_space, const std::vector<std::vector<unsigned> > &point_group_in_mom_idx_space);

    std::vector<std::vector<ff_in_mom_idx_space_t> > 
	make_form_factors_in_momentum_idx_space(const momentum_space_t<dim> &momentum_space);


    unsigned m_shells_count;
    std::vector<std::vector<delta_function_vector_t<dim> > > m_form_factors_shells;

    const double epsilon = 1e-10;

    
protected:
    ff_in_mom_idx_space_t 
        make_form_factor_in_momentum_idx_space(const momentum_space_t<dim> &momentum_space, delta_function_vector_t<dim> form_factor);
    
    double 
        compare_idx_space_form_factors(ff_in_mom_idx_space_t ff1, ff_in_mom_idx_space_t ff2);
    

};

#include "../src/form_factor_container.tpp"
