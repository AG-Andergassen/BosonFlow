#pragma once

#include <grid.h>

template <unsigned dim>
class momentum_space_t : public grid_t<dim>
{
 public:
    momentum_space_t(const std::vector<coord_t<dim> > &points);

    // span a fundamental primitive cell in the lattice. Returns the mesh and the minimum interpoint spacing (needed for computing errors)
    std::tuple<std::vector<coord_t<dim> >, double> span_patched_fundamental_cell(basis_t<dim> reciprocal_basis, unsigned points_per_momentum_cell_edge) const;

    // for debug purposes, spews out a string with a python-style list of the points on the momentum_mesh
    std::string make_debug_mesh_python_list(std::string var_name) const;

    // calculates exp(ik) = cos k + i sin k on the (default) momentum lattice
    std::vector< complex_vector_t<dim> > precalculate_exp_ik() const;

    // calculates exp(ikr) = cos k.r + i sin k.r on the (default) momentum lattice and given set of real lattice vectors 
    std::vector<std::vector< std::complex<double> > > precalculate_exp_ikr(std::vector<coord_t<dim> > rs) const;

    // calculates the corresponding indices to taking the negative of momentum coordinates. Returns a map in index space (from index of momentum, to index of the negative of that momentum) as well as the error.
    std::tuple<std::vector<unsigned>, double> precalculate_negative_indices(const unsigned grid_points_count) const;
    
    // calculates the corresponding indices to summing up two momentum coordinates. Returns a 2D sum map in index space (from the two indices of the momenta, to the index of the sum of the aforementioned momenta) as well as the error
    std::tuple<std::vector<std::vector<unsigned> >, double> precalculate_sum_of_indices(const unsigned grid_points_count) const;

    double get_volume() const;

    double m_volume;
};

#include "../src/momentum_space.tpp"
