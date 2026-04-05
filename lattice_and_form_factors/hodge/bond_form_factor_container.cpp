#include <bond_form_factor_container.h>

template <unsigned dim>
bond_form_factor_container_t<dim>::bond_form_factor_container(lattice_t<dim> lattice, unsigned shells_count): m_lattice(lattice), m_shells_count(shells_count)
{
    m_form_factors_shells = make_form_factors_from_shells(m_lattice.m_shells);
    m_form_factors_shells_in_momentum_idx_space = make_form_factors_in_momentum_idx_space(m_form_factors_shells);
}


template <unsigned dim>
std::vector<std::vector<delta_function_vector_t<dim> > > bond_form_factor_container_t<dim>::make_form_factors_from_shells(std::map<double, std::vector<coord_t<dim> > > shells, unsigned shells_count)
{
    std::vector<std::vector<delta_function_vector_t<dim> > > form_factors_shells;
    unsigned shell_counter = 0;
    for (auto &shell_iter : shells){
	if (shell_counter == shells_count)
	    break;
	std::vector<delta_function_vector_t<dim> > form_factors_in_shell;
	for (coord_t<dim> &bond_vec : shell_iter->second){
	    // for every bond in the shell, we define a localised form factor in that shell (a delta function)
	    delta_function_vector_t<dim> bond({bond_vec}, {1.0});
	    form_factors_in_shell.push_back(bond);
	}
	form_factors_shells.push_back(form_factors_in_shell);

	shell_counter ++;
    }
}

template <unsigned dim>
std::vector<std::complex<double> > bond_form_factor_container_t<dim>::make_form_factor_in_mom_idx_space(delta_function_vector_t<dim> form_factor)
{
    std::vector<std::complex<double> > idx_form_factor;
    for (double mom : m_lattice.m_BZ_mesh)
	idx_form_factor.push_back(form_factor.evaluate_fourier_transform(mom));
    return idx_form_factor;
}

std::vector<std::vector<std::vector<std::complex<double> > > > bond_form_factor_container_t<dim>::make_form_factors_in_momentum_idx_space(std::vector<std::vector<delta_function_vector_t<dim> > > form_factors_shells)
{
    std::vector<std::vector<std::vector<std::complex<double> > > > FFs;

    for (unsigned shell_idx = 0; shell_idx < m_shells_count; shell_idx++){
	std::vector<std::vector< std::complex<double> > > shells_FFs;
	for (unsigned bond_idx = 0; bond_idx < m_lattice.m_BZ_mesh; bond_idx++){
	    shells_FFs.push_back(make_form_factor_in_mom_idx_space(form_factors_shells[shell_idx][bond_idx]));
	}
	FFs.push_back(shells_FFs);
    }
    return FFs;
}

template <unsigned dim>
delta_function_vector_t<dim> bond_form_factor_container_t<dim>::get(unsigned shell_idx, unsigned bond_idx)
{
    return m_form_factors_shells[shell_idx][bond_idx];
}
