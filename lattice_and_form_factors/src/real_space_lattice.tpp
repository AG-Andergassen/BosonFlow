#include <real_space_lattice.h>


template <unsigned dim>
real_space_lattice_t<dim>::real_space_lattice_t(basis_t<dim> real_basis, unsigned real_lattice_range) : grid_t<dim>(real_basis, real_lattice_range) 
{
    m_basis = real_basis;
    
    for (unsigned i = 0; i < dim; i ++)
	m_basis_matrix.col(i) = real_basis[i];

    // todo: if fails, should have throw an error that basis vectors aren't linearly independent.
    m_basis_matrix_inv = m_basis_matrix.inverse();

    m_norms_squared = grid_t<dim>::compute_norms_squared(grid_t<dim>::m_points);

    // initialises also shells_indices
    m_shells = grid_t<dim>::make_shells(grid_t<dim>::m_points, m_norms_squared, m_shells_indices);
}

template <unsigned dim>
coord_units_t<dim> real_space_lattice_t<dim>::make_lattice_coord_in_lattice_basis(const coord_t<dim> coord) const
{
    coord_t<dim> transformed_coord = m_basis_matrix_inv * coord;
    coord_units_t<dim> new_coord;
    for (unsigned i = 0; i < dim; i ++)
	new_coord(i) = static_cast<int>(std::round(transformed_coord(i)));
    return new_coord; 
}
    
template <unsigned dim>
std::vector<coord_t<dim> > real_space_lattice_t<dim>::get_subset(unsigned shells_count) const
{
    std::vector<coord_t<dim> > pts;
    unsigned i = 0;
    for (auto itr = m_shells.begin(); itr != m_shells.end(); itr ++){
	if (i >= shells_count)
	    break;
	i++;

	const auto &shell = itr->second;
	for (auto pt : shell)
	    pts.push_back(pt);	
    }    
    return pts;
}

template <unsigned dim>
std::vector< unsigned > real_space_lattice_t<dim>::get_subset_indices(unsigned shells_count) const
{
    std::vector<unsigned > indices;
    unsigned i = 0;
    for (auto itr = m_shells_indices.begin(); itr != m_shells_indices.end(); itr ++){
	if (i >= shells_count)
	    break;
	i++;

	const auto &shell_indices = itr->second;
	for (auto idx : shell_indices)
	    indices.push_back(idx);	
    }    
    return indices;
}


template <unsigned dim>
basis_t<dim> real_space_lattice_t<dim>::make_reciprocal_lattice_basis() const
{
    basis_t<dim> reciprocal_basis;

    if (dim == 0){
	reciprocal_basis[0](0) = 1.0;
	return reciprocal_basis;
    }
    double prefactor = 2 * M_PI;// / compute_enclosed_volume(real_basis);
    
    // evaluates $k_j = ^\star (\bigwedge_{j \neq i} r_j)$
    std::vector<unsigned> indices;
    for (unsigned i = 0; i < dim; i++)
	indices.push_back(i);

    std::vector<std::tuple<std::vector<unsigned>, int> > S = grid_t<dim>::make_permutations(indices);
    

    for (unsigned j = 0; j < dim; j ++){ 
	coord_t<dim> k_j = coord_t<dim>::Zero();

	for (const auto &sigma : S){ 
	    int sign_of_sigma = std::get<1>(sigma);
	    std::vector<unsigned> permuted_indices = std::get<0>(sigma);
	    coord_t<dim> e_l = coord_t<dim>::Zero(); 
	    e_l(permuted_indices[0]) = 1;

	    double val = 1;
	    for (unsigned i = 0; i < dim; i ++){
		if (i == j)
		    continue;
		if (i == 0)
		    val *= m_basis[i](permuted_indices[j]);
		else
		    val *= m_basis[i](permuted_indices[i]);
	    }
	    k_j += e_l * sign_of_sigma * val;
	    
	}

	reciprocal_basis[j] =  k_j * prefactor / k_j.dot(m_basis[j]) ;
    }

    return reciprocal_basis;
}
