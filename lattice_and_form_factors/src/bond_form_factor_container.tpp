#include <bond_form_factor_container.h>

template <unsigned dim>
bond_form_factor_container_t<dim>::bond_form_factor_container_t(real_space_lattice_t<dim> lattice, unsigned shells_count ): form_factor_container_t<dim>(shells_count), m_lattice(lattice)
{
    form_factor_container_t<dim>::m_form_factors_shells = make_form_factors_from_lattice_shells(m_lattice.m_shells, form_factor_container_t<dim>::m_shells_count);
}


template <unsigned dim>
std::vector<std::vector<delta_function_vector_t<dim> > > bond_form_factor_container_t<dim>::make_form_factors_from_lattice_shells(std::map<double, std::vector<coord_t<dim> > > shells, unsigned shells_count)
{
    std::vector<std::vector<delta_function_vector_t<dim> > > form_factors_shells;
    unsigned shell_counter = 0;
    for (auto shell_iter = shells.begin(); shell_iter != shells.end(); shell_iter ++){
	if (shell_counter == shells_count)
	    break;
	std::vector<delta_function_vector_t<dim> > form_factors_in_shell;
	for (const coord_t<dim> &bond_vec : shell_iter->second){
	    // for every bond in the shell, we define a localised form factor in that shell (a delta function)
	    delta_function_vector_t<dim> bond({bond_vec}, {1.0});
	    form_factors_in_shell.push_back(bond);
	}
	form_factors_shells.push_back(form_factors_in_shell);

	shell_counter ++;
    }
    return form_factors_shells;
}
