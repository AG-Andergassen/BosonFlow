#include <momentum_space.h>
#include <limits>

template <unsigned dim>
momentum_space_t<dim>::momentum_space_t(const std::vector<coord_t<dim> > &points): grid_t<dim>(points)
{
    // for a momentum space that is just a discrete set of points, the volume is just the number of elements.
    m_volume = grid_t<dim>::m_size;
}

template <unsigned dim>
std::tuple<std::vector<coord_t<dim> >, double > momentum_space_t<dim>::span_patched_fundamental_cell(basis_t<dim> reciprocal_basis, unsigned points_per_momentum_cell_edge) const
{
    double min_spacing = std::numeric_limits<double>::max();

    for (unsigned i = 0; i < reciprocal_basis.size(); i ++) 
	min_spacing = std::min(std::sqrt(reciprocal_basis[i].dot(reciprocal_basis[i]))/(points_per_momentum_cell_edge+1), min_spacing);

    if (points_per_momentum_cell_edge == 0){
	std::vector<coord_t<dim> > mom_grid = {coord_t<dim>::Zero()};
	return std::make_tuple(mom_grid, min_spacing);
    }

    std::vector<coord_t<dim> > cell;

    auto coeffs_container = grid_t<dim>::make_integer_coeffs(points_per_momentum_cell_edge/2,  (points_per_momentum_cell_edge%2) );

    for (auto coeffs : coeffs_container){
	coord_t<dim> vec = coord_t<dim>::Zero();
	for (unsigned i = 0; i < vec.size(); i ++){
	    vec += (reciprocal_basis[i]*coeffs[i]) / points_per_momentum_cell_edge;
	    vec += (reciprocal_basis[i] * (points_per_momentum_cell_edge/2))/points_per_momentum_cell_edge; // shift the cell so that the zero point is at the corner. Note that this factor multplying the recirpocal basis is not strictly 0.5 for an odd points_per_momentum_cell_edge
	}
	cell.push_back( vec);
    }


    return std::make_tuple(cell, min_spacing);
}

template <unsigned dim>
std::vector< complex_vector_t<dim> > momentum_space_t<dim>::precalculate_exp_ik() const
{
    std::vector< complex_vector_t<dim> > exp_iks;

    for (auto pt : grid_t<dim>::m_points){
	complex_vector_t<dim> value_at_pt;
	for (unsigned i = 0; i < value_at_pt.size(); i ++){
	    value_at_pt(i) = std::complex( std::cos(pt(i)), std::sin(pt(i)) );
	}
	exp_iks.push_back(value_at_pt);
    }
    return exp_iks;
}

template <unsigned dim>
std::vector<std::vector< std::complex<double> > > momentum_space_t<dim>::precalculate_exp_ikr(std::vector<coord_t<dim> > rs) const
{
    std::vector<std::vector< std::complex<double> > > exp_ikrs;

    for (auto k_pt : grid_t<dim>::m_points){
	std::vector<std::complex<double> > with_rs;
	for (auto r_pt : rs){
	    double phase = 0;
	    std::complex<double> val;
	    for (unsigned i = 0; i < k_pt.size(); i ++){
		phase += k_pt(i) * r_pt(i);
	    }
	    with_rs.push_back(std::complex( std::cos(phase), std::sin(phase)) );
	}
	exp_ikrs.push_back(with_rs);
    }
    return exp_ikrs;
}

template <unsigned dim>
double momentum_space_t<dim>::get_volume() const
{
    return m_volume;
}


template <unsigned dim>
std::tuple<std::vector<unsigned>, double> momentum_space_t<dim>::precalculate_negative_indices(const unsigned grid_points_count) const
{
    double max_error = -1;

    std::vector<unsigned> neg_mom_idxes;
    neg_mom_idxes.resize(grid_t<dim>::m_points.size());

    for (unsigned mom_idx = 0; mom_idx < grid_t<dim>::m_points.size(); mom_idx++){
	auto [new_idx, error] = this->get_idx_from_coord(-grid_t<dim>::m_points[mom_idx], grid_points_count);
	max_error = std::max(max_error, error);
	neg_mom_idxes[mom_idx] = new_idx;
    }

    return std::make_tuple(neg_mom_idxes, max_error);
}

template <unsigned dim>
std::tuple<std::vector<std::vector<unsigned> >, double> momentum_space_t<dim>::precalculate_sum_of_indices(const unsigned grid_points_count) const
{
    double max_error = -1;

    std::vector<std::vector<unsigned> > sum_moms_idxes;
    // todo: think about whether to take m_size or m_points.size()...
    
    sum_moms_idxes.resize(grid_t<dim>::m_points.size());
    for (unsigned i = 0; i < grid_t<dim>::m_points.size(); i++){
	sum_moms_idxes[i].resize(grid_t<dim>::m_points.size());
    }
    
#pragma omp parallel default(shared)
    {
    #pragma omp for schedule(nonmonotonic:dynamic) nowait
    for (unsigned mom1_idx = 0; mom1_idx < grid_t<dim>::m_points.size(); mom1_idx++)
    for (unsigned mom2_idx = 0; mom2_idx < grid_t<dim>::m_points.size(); mom2_idx++){
	auto [new_idx, error] = this->get_idx_from_coord(grid_t<dim>::m_points[mom1_idx] + grid_t<dim>::m_points[mom2_idx], grid_points_count);
#pragma omp critical
	max_error = std::max(max_error, error);
	
	sum_moms_idxes[mom1_idx][mom2_idx] = new_idx;
    }
    }
	
    return std::make_tuple(sum_moms_idxes, max_error);
}
