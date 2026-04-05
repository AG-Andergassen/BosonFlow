#include <lattice_point_group.h>

template <unsigned dim>
lattice_point_group_t<dim>::lattice_point_group_t(const std::vector<std::vector<matrix_t<dim> > > &group_conjugacy_classes_matrices, std::vector<std::string> element_names)
{
    m_element_names = element_names;
    for (auto conj_class : group_conjugacy_classes_matrices)
	m_group_matrices.insert(m_group_matrices.end(), conj_class.begin(), conj_class.end());
    
    m_group_size = m_group_matrices.size();
}


template <unsigned dim>
lattice_point_group_t<dim>::lattice_point_group_t(const std::vector<matrix_t<dim> > &group_matrices, std::vector<std::string> element_names): m_group_matrices(group_matrices)
{
    m_element_names = element_names;
    m_group_size = m_group_matrices.size();
}


template <unsigned dim>
std::tuple<std::vector<std::vector<unsigned> >, double> lattice_point_group_t<dim>::make_group_in_momentum_idx_space(momentum_space_t<dim> &momentum_space, double error_threshold)
{
    double max_error;
    int error_idx;
    matrix_t<dim> error_g;
    bool is_first = true;
    
    std::vector<std::vector<unsigned> > group_in_momentum_idx_space;

    for (unsigned g_idx = 0; g_idx < m_group_size; g_idx ++){
	matrix_t<dim> g = m_group_matrices[g_idx];
	std::vector<unsigned> g_in_idx_space;
	for (unsigned mom_idx = 0; mom_idx < momentum_space.m_points.size(); mom_idx ++){
	    vector_t<dim> new_mom_vec = g.inverse() * momentum_space.m_points[mom_idx];
	    auto [new_idx, error] = momentum_space.get_idx_from_coord(new_mom_vec);
	    g_in_idx_space.push_back(new_idx);
	    
	    if (is_first || error > max_error){
		is_first = false;
		max_error = error;
		error_g = g;
	    }

	    if (error > error_threshold)
		throw std::invalid_argument("momentum grid does not respect symmetry with given error threshold! The error is " + std::to_string(error) + " with group element of index: " + std::to_string(g_idx));
	    
	}
	group_in_momentum_idx_space.push_back(g_in_idx_space);
    }
    
    return std::make_tuple(group_in_momentum_idx_space, max_error);
}

