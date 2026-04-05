#pragma once

#include <Eigen/Core>
#include <Eigen/LU> // for Eigen::Matrix::inverse()

#include <momentum_space.h>

using namespace Eigen;


template <unsigned dim>
class lattice_point_group_t
{
 public: 
    lattice_point_group_t(const std::vector<std::vector<matrix_t<dim> > > &group_conjugacy_classes_matrices, std::vector<std::string> element_names = {});
    
    lattice_point_group_t(const std::vector<matrix_t<dim> > &group_matrices, std::vector<std::string> element_names = {});
    

    // Generates a map for the group action between the BZ indices instead of 2D coordinates. The error here refers to the maximum misalignment of points after the symmetry transformation, roughly in units of lattice spacing. Returns the maximum error.
    std::tuple<std::vector<std::vector<unsigned> >, double> make_group_in_momentum_idx_space(momentum_space_t<dim> &momentum_space, double error_threshold = 0.1);

    // private:
    std::vector<matrix_t<dim> > m_group_matrices;
    unsigned m_group_size;

    std::vector<std::string> m_element_names;
};

#include "../src/lattice_point_group.tpp"
