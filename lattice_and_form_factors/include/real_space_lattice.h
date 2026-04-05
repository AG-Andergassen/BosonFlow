#pragma once 

#include <grid.h>

template <unsigned dim>
class real_space_lattice_t : public grid_t<dim>
{
 public:
    real_space_lattice_t(basis_t<dim> real_basis, unsigned real_lattice_range);

    // convert coord to the lattice basis
    coord_units_t<dim> make_lattice_coord_in_lattice_basis(const coord_t<dim> coord) const;

    // get a subset of the lattice up to only shell shells_count
    std::vector<coord_t<dim> > get_subset(unsigned shells_count) const;

    // get a subset of the lattice indices that correspond to vectors up to only shell shells_count
    std::vector<unsigned > get_subset_indices(unsigned shells_count) const;

    // calculate reciprocal lattice basis vectors from real lattice basis vectors
    basis_t<dim> make_reciprocal_lattice_basis() const;

    basis_t<dim> m_basis;
    matrix_t<dim> m_basis_matrix;
    matrix_t<dim> m_basis_matrix_inv;

    std::vector<double> m_norms_squared;
    std::map<double, std::vector<coord_t<dim> > > m_shells;
    std::map<double, std::vector<unsigned > > m_shells_indices;

};

#include "../src/real_space_lattice.tpp"
