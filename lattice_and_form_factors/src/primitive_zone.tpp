#include <primitive_zone.h>

//#include <ranges>  // requires C++20

template <unsigned dim>
primitive_zone_t<dim>::primitive_zone_t(): momentum_space_t<dim>({})
{
}


template <unsigned dim>
primitive_zone_t<dim>::primitive_zone_t(basis_t<dim> reciprocal_basis, unsigned points_per_momentum_cell_edge): momentum_space_t<dim>({})
{
    m_reciprocal_basis = reciprocal_basis;
    m_points_per_momentum_cell_edge = points_per_momentum_cell_edge;
    std::tie(grid_t<dim>::m_points, m_momentum_mesh_min_interpoint_spacing) = momentum_space_t<dim>::span_patched_fundamental_cell(m_reciprocal_basis, points_per_momentum_cell_edge);
    grid_t<dim>::m_size = grid_t<dim>::m_points.size();
    momentum_space_t<dim>::m_volume = grid_t<dim>::compute_enclosed_volume(reciprocal_basis);

    // since grid is homogeneous, assign equal weights to all points. They must sum up to 1
    for (auto p: grid_t<dim>::m_points)
	m_volume_weights.push_back(1.0/grid_t<dim>::m_points.size());
    m_volume_weights_unrefined = m_volume_weights;
    
    m_points_count_unrefined = grid_t<dim>::m_points.size();


    for (unsigned i = 0; i < m_reciprocal_basis.size(); i ++)
	m_reciprocal_basis_matrix.col(i) = m_reciprocal_basis[i];
    
    m_reciprocal_basis_matrix_inv = m_reciprocal_basis_matrix.inverse();

    // for binary search in the grid
    grid_t<dim>::populate_search_bins();
}



template <unsigned dim>
std::tuple<unsigned, double> primitive_zone_t<dim>::get_idx_from_coord(coord_t<dim> coord) const
{
  return get_idx_from_coord(coord, grid_t<dim>::m_points.size());
}

template <unsigned dim>
std::tuple<unsigned, double> primitive_zone_t<dim>::get_idx_from_coord_unrefined(coord_t<dim> coord) const
{
  return get_idx_from_coord(coord, m_points_count_unrefined);
}

template <unsigned dim>
  std::tuple<unsigned, double> primitive_zone_t<dim>::get_idx_from_coord(coord_t<dim> coord, unsigned grid_points_count) const
{ 
    std::vector<std::tuple<unsigned, double> > idxes_and_deviations;
    bool is_first = true;
    
    unsigned smallest_deviation_idx;
    double smallest_deviation;

    for (unsigned i = 0; i < coord.size(); i++){
	double reciprocal_basis_length_squared = m_reciprocal_basis[i].dot(m_reciprocal_basis[i]);
	double normalised_coord_in_reciprocal_basis = coord.dot(m_reciprocal_basis[i])/reciprocal_basis_length_squared;
	// coarsely reduce large winding number in this reciprocal lattice direction if too large > 1
	while ( std::abs(normalised_coord_in_reciprocal_basis) > 0.5 ){
	    int sign = (normalised_coord_in_reciprocal_basis > 0) ? 1 : ((normalised_coord_in_reciprocal_basis < 0) ? -1 : 0);
	    coord = coord - m_reciprocal_basis[i] * sign;
	    normalised_coord_in_reciprocal_basis = coord.dot(m_reciprocal_basis[i])/reciprocal_basis_length_squared;
	}
    }

    // the "finer" check for small winding number periodicity (in all dimension, limited to winding numbers 1, 0, -1)
    auto winding_coeffs_container = grid_t<dim>::make_integer_coeffs(1);
    
    // we try "rewinding" the coordinate with windings (1, 0, -1), then take the coordinates that conindices with a point in the brillouin zone
    for (auto winding_coeffs : winding_coeffs_container){
	coord_t<dim> winded_coord = coord;
	for (int i = 0; i < coord.size(); i ++)
	    winded_coord += m_reciprocal_basis[i] * winding_coeffs[i];
	auto [idx, deviation] = grid_t<dim>::get_idx_from_coord(winded_coord, grid_points_count);
	if (is_first){
	    is_first = false;
	    smallest_deviation_idx = idx;
	    smallest_deviation = deviation;
	    if (smallest_deviation < grid_t<dim>::epsilon)
		break;
	    
	    continue;
	}
	if (deviation < smallest_deviation){
	    smallest_deviation_idx = idx;
	    smallest_deviation = deviation;
	}
    }

    return std::make_tuple(smallest_deviation_idx, smallest_deviation/m_momentum_mesh_min_interpoint_spacing);
}

template <unsigned dim>
std::vector<unsigned> primitive_zone_t<dim>::get_path_indices(std::vector<coord_t<dim> > points, bool is_closed)
{
    return get_path_indices(points, is_closed, grid_t<dim>::m_points.size());
}

template <unsigned dim>
std::vector<unsigned> primitive_zone_t<dim>::get_path_indices_unrefined(std::vector<coord_t<dim> > points, bool is_closed)
{
    return get_path_indices(points, is_closed, m_points_count_unrefined);
}

template <unsigned dim>
std::vector<unsigned> primitive_zone_t<dim>::get_path_indices(std::vector<coord_t<dim> > points, bool is_closed, unsigned grid_points_count)
{
  if (points.size() < 2){
    return {std::get<0>(get_idx_from_coord(points[0], grid_points_count))};
  }

  std::vector<unsigned> path;
  
  // todo: rename m_momentum_mesh_min_interpoint_spacing to m_momentum_mesh_typical_interpoint_spacing ... Since now refining should make the min distance smaller. And it doesn't work well if that changes for the error calculation in get_idx_from_coord
  
  // get a safe path discretisation so as to not miss any points
  double min_dist;
  
  {

    auto smallestNonZeroNormSquared = [&](const std::vector<coord_t<dim> >& coords, double eps) {
      std::vector<double> norms;
      for (const auto& coord : coords)
	  norms.push_back(coord.dot(coord));
        
      std::sort(norms.begin(), norms.end()); // Sort norms

      // Remove near-duplicates
      norms.erase(std::unique(norms.begin(), norms.end(), [&](double a, double b) { return std::fabs(a - b) < eps; }), norms.end());
      
      return norms;
    };
    auto norms_squared = smallestNonZeroNormSquared(this->m_points, grid_t<dim>::epsilon);


    // auto shells = this->make_shells(this->m_points, smallest_norm_squared + epsilon);
    //std::vector<double> norms_squared_hierarchy;
    //for (auto itr = shells.begin(); itr != shells.end(); itr ++){
    //  norms_squared_hierarchy.push_back(itr->first);
    //}
    min_dist = std::sqrt(norms_squared[1]) /2.0; //(std::sqrt(norms_squared_hierarchy[1]) - std::sqrt(norms_squared_hierarchy[0])) / 2.0;
    for (unsigned i = 2; i < norms_squared.size(); i ++){
      double min_dist_candidate = (std::sqrt(norms_squared[i]) - std::sqrt(norms_squared[i-1])) / 2.0;
      if (min_dist > min_dist_candidate)
	  min_dist = min_dist_candidate;
    }
  }
  
  unsigned until = is_closed ? points.size() : points.size() - 1;

  for (unsigned i = 0; i < until; i ++){
    coord_t<dim> start_pt = points[i];
    coord_t<dim> end_pt = points[(i+1)%points.size()];

    coord_t<dim> d = (start_pt - end_pt);

    double max_N = std::round(std::sqrt(d.dot(d))/min_dist);

    path.push_back(std::get<0>(get_idx_from_coord(start_pt)));
    for (int n = 1; n <= max_N; n++){
      unsigned pt_idx_candidate = std::get<0>(get_idx_from_coord(start_pt * (1.0 - n/max_N) + end_pt * (n / max_N), grid_points_count));
      if (pt_idx_candidate != path.back())
	path.push_back(pt_idx_candidate);
    }
  }
  return path;
}

// TODO: support symmetries
template <unsigned dim>
void primitive_zone_t<dim>::refine_from_finer_grid(momentum_space_t<dim> *fine_grid_ptr_, coord_t<dim> center_of_refinement_, float refine_percent)
{

  // project center to *unrefined grid*. This is important for the volume weights to be calculated better (perfectly). 
  //coord_t<dim> center_of_refinement = this->m_points[std::get<0>(get_idx_from_coord_unrefined(center_of_refinement_))];
  coord_t<dim> center_of_refinement = fine_grid_ptr_->m_points[std::get<0>(fine_grid_ptr_->get_idx_from_coord(center_of_refinement_))];
 

  // recenter grid so that pt 0 is at the center. 
  auto fine_grid = *fine_grid_ptr_;
  coord_t<dim> fine_grid_center_coord;
  for (int i = 0; i < fine_grid_center_coord.size(); i++){
      fine_grid_center_coord(i) = 0.0;
      for (int j = 0; j < fine_grid_center_coord.size(); j ++)
	  fine_grid_center_coord(i) += m_reciprocal_basis[j][i]/2.0;
  }
  
  fine_grid_ptr_->populate_search_bins();

  // project to actual center on fine grid
  fine_grid_center_coord = fine_grid_ptr_->m_points[std::get<0>(fine_grid_ptr_->get_idx_from_coord(fine_grid_center_coord))];

  // grid adjusted so that e.g. the coordinates go from -pi to pi instead of from 0 to 2pi
  for (int i = 0; i <  fine_grid.m_points.size(); i ++)
	  fine_grid.m_points[i] -= fine_grid_center_coord;
  // for get_idx_from_coord to work fast assuming no periodicity handling is needed
  std::cout << "constructing kd tree" << std::endl;

  fine_grid.populate_search_bins();

  std::cout << "shifting complete" << std::endl;
  
  

  // note: this fine_grid_ptr points to a grid that has manually been shifted.. this makes using get_idx_from_coord on it sometimes return wrong results (periodicity handling gets broken)
  auto fine_grid_ptr = &fine_grid;


  auto smallestNormSquareds = [&](const std::vector<coord_t<dim> >& coords, double p, double eps) {
      std::vector<double> norm_squareds;
      for (const auto& coord : coords)
	  norm_squareds.push_back(coord.dot(coord));
        
      std::sort(norm_squareds.begin(), norm_squareds.end()); // Sort norms

      // Remove near-duplicates
      norm_squareds.erase(std::unique(norm_squareds.begin(), norm_squareds.end(), [&](double a, double b) {          return std::fabs(a - b) < eps; }), norm_squareds.end());
      
      size_t count = static_cast<size_t>(p * norm_squareds.size() ); // p% of the size
      return std::vector<double>(norm_squareds.begin(), norm_squareds.begin() + std::min(count, norm_squareds.size()));
  };
  auto smallest_norm_squareds = smallestNormSquareds(fine_grid_ptr->m_points, refine_percent, 1e-12);

  
  std::map<double, std::vector<unsigned> > fine_shells_indices;

  
  auto fine_shells = this->make_shells(fine_grid_ptr->m_points, fine_shells_indices, smallest_norm_squareds.back());
  std::vector<double> norms_squared_hierarchy;
  for (auto itr = fine_shells.begin(); itr != fine_shells.end(); itr ++){
    norms_squared_hierarchy.push_back(itr->first);
  }

  
  // pick a percentage of the fine_shells, then take size 
  unsigned hierarchy_size = norms_squared_hierarchy.size() - 1; //std::max(0.0f, std::round(norms_squared_hierarchy.size()*refine_percent - 1));
  std::cout << norms_squared_hierarchy.back() << ",," << smallest_norm_squareds.back() << std::endl;
  
  double max_refined_norm_squared = smallest_norm_squareds.back();// norms_squared_hierarchy[hierarchy_size];

  std::vector<unsigned> refined_pts_indices;
  for (unsigned shell_idx = 0; shell_idx <= hierarchy_size ; shell_idx ++){
    for (unsigned pt_idx = 0; pt_idx < fine_shells_indices[norms_squared_hierarchy[shell_idx]].size(); pt_idx ++){
      refined_pts_indices.push_back(fine_shells_indices[norms_squared_hierarchy[shell_idx]][pt_idx]);
    }
  }
  
  // get the norms of the coarse shells around refinement
  std::map<double, std::vector<unsigned> > coarse_shells_indices;
  std::vector<coord_t<dim> > coarse_grid_translated_extended;
  {
    auto winding_coeffs_container = grid_t<dim>::make_integer_coeffs(1);
    
    for (auto winding_coeffs : winding_coeffs_container){
      coord_t<dim> winded_coord = coord_t<dim>::Zero();
      for (auto pt:this->m_points){
      for (int i = 0; i < m_reciprocal_basis.size(); i ++)
	  winded_coord += m_reciprocal_basis[i] * winding_coeffs[i];
	
	coarse_grid_translated_extended.push_back(pt - center_of_refinement + winded_coord);
      }
    }
  }
  
  auto coarse_shells = this->make_shells(coarse_grid_translated_extended, coarse_shells_indices, max_refined_norm_squared);
  std::vector<double> coarse_norms_squared_hierarchy;
  for (auto itr = coarse_shells.begin(); itr != coarse_shells.end(); itr ++){
    coarse_norms_squared_hierarchy.push_back(itr->first);
  }
  
      
  // get max coarse shell that will be refined
  unsigned max_coarse_shell = coarse_norms_squared_hierarchy.size() - 1;
  /*  for (max_coarse_shell = 0; max_coarse_shell < coarse_norms_squared_hierarchy.size() && coarse_norms_squared_hierarchy[max_coarse_shell] <= max_refined_norm_squared; max_coarse_shell ++) {
  }
  std::cout << max_coarse_shell << "," << coarse_norms_squared_hierarchy.size() <<  std::endl;*/
  std::cout << coarse_norms_squared_hierarchy[max_coarse_shell] << ",," << max_refined_norm_squared << std::endl;
  
  // get which of the points coming from the refined lattice are really new
  std::vector<unsigned> existing_pts_fine_indices;
  std::vector<unsigned> existing_pts_coarse_indices;
  
  std::vector<unsigned > coarse_pts_in_refined_region_indices;
  for (unsigned shell_idx = 0; shell_idx < max_coarse_shell+1; shell_idx ++){
    for (unsigned pt_idx = 0; pt_idx < coarse_shells[coarse_norms_squared_hierarchy[shell_idx]].size(); pt_idx ++){
      double error;
      unsigned idx;
      std::tie(idx, error) = fine_grid_ptr->get_idx_from_coord(coarse_shells[coarse_norms_squared_hierarchy[shell_idx]][pt_idx]);
      if (error < grid_t<dim>::epsilon){
	existing_pts_fine_indices.push_back(idx);
	//here, we're also wrapping to reduced primitive/brillouin zone (update: this naive wrapping is problematic!!)
	existing_pts_coarse_indices.push_back(coarse_shells_indices[coarse_norms_squared_hierarchy[shell_idx]][pt_idx] % this->m_points.size());
      }
      coarse_pts_in_refined_region_indices.push_back(coarse_shells_indices[coarse_norms_squared_hierarchy[shell_idx]][pt_idx] % this->m_points.size() );
    }
  }
  

  std::vector<coord_t<dim> > new_pts;
  for (unsigned idx : refined_pts_indices){
    // if it doesn't already exist, add it
    if (std::find(existing_pts_fine_indices.begin(), existing_pts_fine_indices.end(), idx) == existing_pts_fine_indices.end()){
      new_pts.push_back(fine_grid_ptr->m_points[idx]);
    }
  }

  double total_coarse_refined_region_weight = 0;
  for (auto idx : coarse_pts_in_refined_region_indices){
    total_coarse_refined_region_weight += m_volume_weights[idx];
  }
  
  double new_refined_region_pts_weight = total_coarse_refined_region_weight/(refined_pts_indices.size());
  
  for (auto idx : existing_pts_coarse_indices){
    this->m_volume_weights[idx] = new_refined_region_pts_weight;
  }
  std::cout << new_pts.size() << std::endl;
  for (auto pt: new_pts){
    this->m_points.push_back(fine_grid_ptr_->m_points[std::get<0>(fine_grid_ptr_->get_idx_from_coord(pt + center_of_refinement))]);
    this->m_volume_weights.push_back(new_refined_region_pts_weight);
  }
  
  this->populate_search_bins();
  
}


template <unsigned dim>
void primitive_zone_t<dim>::refine_from_finer_grid_along_path(momentum_space_t<dim> *fine_grid_ptr, std::vector<coord_t<dim>> points, double refine_percent, bool is_closed_path)
{
  if (points.size() == 0)
      return;
  auto path = get_path_indices_unrefined(points, is_closed_path);
  for (auto pt_idx : path)
    refine_from_finer_grid(fine_grid_ptr, this->m_points[pt_idx], refine_percent);
}
