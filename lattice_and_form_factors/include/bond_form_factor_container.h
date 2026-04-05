#pragma once

#include <form_factor_container.h>
#include <real_space_lattice.h>
#include <delta_function_vector.h>


template <unsigned dim>
class bond_form_factor_container_t : public form_factor_container_t<dim>
{
 public:
    bond_form_factor_container_t(real_space_lattice_t<dim> lattice, unsigned shells_count);

    std::vector<std::vector<delta_function_vector_t<dim> > > 
	make_form_factors_from_lattice_shells(std::map<double, std::vector<coord_t<dim> > > shells, unsigned shells_count);

    real_space_lattice_t<dim> m_lattice;
};

#include "../src/bond_form_factor_container.tpp"
