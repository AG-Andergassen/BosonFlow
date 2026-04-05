#include <grid.h>

#define USE_KNN_SEARCH

template <unsigned dim>
grid_t<dim>::grid_t(const std::vector<coord_t<dim> > &points): m_points(points), m_grid_kdtree_adaptor(m_points), m_points_kdtree((dim == 0) ? 1 : dim, m_grid_kdtree_adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10))
{
    m_size = m_points.size();

    populate_search_bins();
}

template <unsigned dim>
grid_t<dim>::grid_t(basis_t<dim> basis, unsigned range):m_grid_kdtree_adaptor(m_points), m_points_kdtree((dim == 0) ? 1 : dim, m_grid_kdtree_adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10))
{
    m_points = span_lattice(basis, range);
    m_size = m_points.size();

    populate_search_bins();
}


template <unsigned dim>
std::vector<coord_t<dim> > grid_t<dim>::span_lattice(basis_t<dim> basis, unsigned range ) const
{
    std::vector<coord_t<dim> > span;
    
    auto coeffs_container = make_integer_coeffs(range);

    for (auto coeffs : coeffs_container){
	coord_t<dim> r = coord_t<dim>::Zero();
	for (int i = 0; i < dim; i ++)
	    r += basis[i] * coeffs[i];

	span.push_back(r);
    }
       
    return span;
}

template <unsigned dim>
int_basis_coeffs_t<dim> grid_t<dim>::make_integer_coeffs(unsigned range, bool is_closed_interval) const
{ 
    // choose whether to make integers in [-range, range) or integers in [-range, range]
    unsigned interval_size = is_closed_interval ? 2 * range+1 : 2 * range;

    int coeffs_container_size = std::pow(interval_size, dim);

    int_basis_coeffs_t<dim> coeffs_container;

    // encode the coefficients as the digits of a "number in base $interval_size$". The ith coefficient is then read as the ith digit of this number.
    for (int coeffs_number = 0; coeffs_number < coeffs_container_size; coeffs_number ++){
	int coeffs_number_mutable = coeffs_number;
        typename int_basis_coeffs_t<dim>::value_type coeffs;
	for (int i = coeffs.size() - 1; i >= 0; i--){
	    // extract the first digit and set it to be the ith component
	    coeffs[i] = coeffs_number_mutable % interval_size;
	    
	    // throw away the first digit
	    coeffs_number_mutable -= coeffs[i];
	    coeffs_number_mutable /= interval_size;

	    // 'recenter' the extracted coeffecient so that it goes in [-range, range] (or [-range, range) )
	    coeffs[i] -= range;
	}
	// push the coefficients into the set of coefficients
	coeffs_container.push_back(coeffs);
    }
    return coeffs_container;
}
/*
// 0-dimensional case
template <>
std::vector< std::array<int, 1> > grid_t<0>::make_integer_coeffs(unsigned range, bool is_closed_interval) const
{ 
    std::vector<std::array<int, 1> > coeffs_container;
    std::array<int, 1> coeffs = {0.0,};
    coeffs_container.push_back(coeffs);
    return coeffs_container;
}
*/

template <unsigned dim>
std::vector<double> grid_t<dim>::compute_norms_squared(const std::vector<coord_t<dim> > &vecs) const
{
    std::vector<double> norms_squared;
    for (const coord_t<dim> &vec : vecs){
	double norm_squared = vec.dot(vec);
	norms_squared.push_back(norm_squared);
    } 
    return norms_squared;
}


template <unsigned dim>
std::map<double, std::vector<coord_t<dim> > > grid_t<dim>::make_shells(const std::vector<coord_t<dim> > &span, const std::vector<double> &norms_squared, std::map<double, std::vector<unsigned> > &shells_indices, double max_norm_squared) const
{
    std::map<double, std::vector<coord_t<dim> > > shells;

    if (span.size() == 0)
    	throw std::invalid_argument("The set of points \"span\" given to make shells of is empty.");
    
    for (unsigned i =0; i < norms_squared.size(); i ++){
	bool key_exists = false;
	double coord_norm_squared = norms_squared[i];
	if (coord_norm_squared > max_norm_squared)
	    continue;
	for (auto itr = shells.begin(); itr != shells.end(); itr ++){
	    double norm_squared = itr->first;
	    // if a shell with this distance is already there, add lattice vector to shell. 
	    if (std::abs(coord_norm_squared - norm_squared) < epsilon){
		key_exists = true;
		shells[norm_squared].push_back(span[i]);
		shells_indices[norm_squared].push_back(i);
		break;
	    }
	}
	// if no shell with that distance exists, create a new shell with this lattice vector in it
	if (! key_exists){
	    shells[norms_squared[i]] = {span[i]};
	    shells_indices[norms_squared[i]] = {i};
	}
    }
    return shells;
}


template <unsigned dim>
std::map<double, std::vector<coord_t<dim> > > grid_t<dim>::make_shells(const std::vector<coord_t<dim> > &span, double max_norm_squared) const
{
    std::map<double, std::vector<unsigned> > shells_indices;
    return make_shells(span, shells_indices, max_norm_squared);
}

template <unsigned dim>
std::map<double, std::vector<coord_t<dim> > > grid_t<dim>::make_shells(const std::vector<coord_t<dim> > &span, std::map<double, std::vector<unsigned> > &shells_indices, double max_norm_squared) const
{
    std::vector<double> norms_squared;

    norms_squared = compute_norms_squared(span);
    return make_shells(span, norms_squared, shells_indices, max_norm_squared);
}

template <unsigned dim>
std::map<double, std::vector<coord_t<dim> > > grid_t<dim>::make_shells(const std::vector<coord_t<dim> > &span, std::map<double, std::vector<unsigned> > &shells_indices, const coord_t<dim> center, double max_norm_squared) const
{
    std::vector<double> norms_squared;
    
    std::vector<coord_t<dim> > centered_span = span;
    for (unsigned i = 0; i < centered_span.size(); i ++)
	centered_span[i] = centered_span[i] - center;

    norms_squared = compute_norms_squared(centered_span);
    return make_shells(centered_span, norms_squared, shells_indices, max_norm_squared);
}

template <unsigned dim>
std::vector<std::tuple<std::vector<unsigned>, int> > grid_t<dim>::make_permutations(const std::vector<unsigned> &ls) const
{
    if (ls.size() <= 1)
	return { {ls, +1} };
    
    std::vector<std::tuple<std::vector<unsigned>, int> > permutations;
    for (unsigned i = 0; i < ls.size(); i ++){
	std::vector<unsigned> ls_copy = ls;
	unsigned ele = ls_copy[i];
	ls_copy.erase(ls_copy.begin()+i);
	int perm_sign = std::pow(-1, i);
	std::vector<std::tuple<std::vector<unsigned>, int> > sub_perms = make_permutations(ls_copy);
	for (auto &sub_perm : sub_perms){
	    std::vector<unsigned> new_permutation = {ele};
	    new_permutation.insert(new_permutation.end(), std::get<0>(sub_perm).begin(), std::get<0>(sub_perm).end() );
	    int new_perm_sign = perm_sign * std::get<1>(sub_perm);
	    permutations.push_back({new_permutation, new_perm_sign});
	}
    }
    return permutations;
}


template <unsigned dim>
double grid_t<dim>::compute_enclosed_volume(basis_t<dim> vecs) const
{
    matrix_t<dim> cols;
    for (unsigned i = 0; i < vecs.size(); i ++)
    for (unsigned j = 0; j < vecs[0].size(); j ++){
	cols(i, j) = vecs[i](j);
    }
    return std::abs(cols.determinant());
}



template <unsigned dim>
std::unordered_map<coord_t<dim> , unsigned, CoordHash<coord_t<dim> > > grid_t<dim>::get_coord_to_idx_map(const std::vector<coord_t<dim> > &mesh)
{
    std::unordered_map<coord_t<dim>, unsigned, CoordHash<coord_t<dim> > > coord_to_idx;
    for (unsigned k_idx = 0; k_idx < mesh.size(); k_idx++){
	coord_to_idx[mesh[k_idx]] = k_idx;
    }
    return coord_to_idx;
}



template <unsigned dim>
std::tuple<unsigned, double> grid_t<dim>::get_idx_from_coord(coord_t<dim> coord) const
{
  return grid_t<dim>::get_idx_from_coord(coord, m_points.size());
}

template <unsigned dim>
std::tuple<unsigned, double> grid_t<dim>::get_idx_from_coord(coord_t<dim> coord, unsigned points_count) const
{
  
#ifdef USE_KNN_SEARCH

    size_t ret_index;
    double dist_squared;
    nanoflann::KNNResultSet<double> resultSet(1);
    resultSet.init(&ret_index, &dist_squared);
    
    nanoflann::SearchParameters params;
    params.eps = 0.0;  // Exact search
    
    m_points_kdtree.findNeighbors(resultSet, coord.data(), params);
    return {static_cast<unsigned>(ret_index), std::sqrt(dist_squared)};

#else

    unsigned idx;
    double deviation_squared;
    bool is_first = true;

    // old implementation
    for (unsigned i = 0; i < points_count; i ++){
	coord_t<dim> potential_coord = m_points[i];
	
	double norm_squared = (potential_coord - coord).dot(potential_coord - coord);
	if (is_first || norm_squared < deviation_squared){
	    is_first = false;
	    deviation_squared = norm_squared;
	    idx = i;
	    if (deviation_squared <= epsilon)
	      break;
	}
    }
    
    return {idx, std::sqrt(deviation_squared)};
#endif
}


template<unsigned dim>
void grid_t<dim>::populate_search_bins() {
    m_points_kdtree.buildIndex();
    return;
}

template<unsigned dim>
unsigned grid_t<dim>::get_coord_search_bin_idx(const coord_t<dim>& coord) const {
    if (m_search_bins_of_coord_indices.empty()) {
        throw std::runtime_error("Bins are not populated; call populate_bins first.");
    }


    // Array to hold the min coordinate values for each dimension
    std::array<double, dim> bin_sizes, min_coords;
    for (int d = 0; d < dim; ++d) {
        bin_sizes[d] = (m_max_coord - m_min_coord) / m_num_bins_per_dim; // Example for 0th coord
    }

    unsigned bin_idx = 0;
    unsigned stride = 1;

    for (int d = 0; d < dim; ++d) {
        unsigned coord_bin_idx = static_cast<unsigned>(
						       std::floor((coord[d] - min_coords[d]) / bin_sizes[d]));

        // Clamp to bounds
        coord_bin_idx = std::min(coord_bin_idx, m_num_bins_per_dim - 1);

        bin_idx += coord_bin_idx * stride;
        stride *= m_num_bins_per_dim;
    }

    return bin_idx;
}


template <unsigned dim>
std::string grid_t<dim>::make_debug_grid_python_list(std::string var_name)
{
    std::string output; 
    output += var_name + "= [";

    for (unsigned mom_idx = 0; mom_idx < m_points.size(); mom_idx ++){
	coord_t<dim> pt = m_points.at(mom_idx);
	output += "[";
	for (unsigned i = 0; i < dim; i++){
	    output += std::to_string(pt(i));
	    if (i != pt.size() - 1)
		output += ",";
	}
	output += "]";
	if (mom_idx != m_points.size() - 1)
	    output += ",";
    }

    output += "]\n";
        
    return output;
}


template <unsigned dim>
std::vector<double> grid_t<dim>::get_points_as_std_vectors() const
{
    std::vector<double> pts_data;
    
    for (auto p: m_points){
	for (unsigned i = 0; i < p.size(); i++){
	    pts_data.push_back(p(i));
	}
    }
    return pts_data;
}
