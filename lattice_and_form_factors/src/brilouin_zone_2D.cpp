#include <brillouin_zone_2D.h>


brillouin_zone_2D_t::brillouin_zone_2D_t(std::array<coord_t<2>, 2> reciprocal_basis, std::array<coord_t<2>, 2> bz_grid_basis, unsigned points_per_momentum_cell_edge)
{
    primitive_zone_t<2>::m_reciprocal_basis = reciprocal_basis;
    primitive_zone_t<2>::m_points_per_momentum_cell_edge = points_per_momentum_cell_edge;
    std::tie(m_points, m_momentum_mesh_min_interpoint_spacing) = span_patched_BZ_wigner_seitz(m_reciprocal_basis, bz_grid_basis, points_per_momentum_cell_edge);
    grid_t<2>::m_size = grid_t<2>::m_points.size();
    momentum_space_t<2>::m_volume = compute_enclosed_volume(reciprocal_basis);

    // since grid is homogeneous, assign equal weights to all points. They must sum up to 1
    for (auto p: grid_t<2>::m_points)
      primitive_zone_t<2>::m_volume_weights.push_back(1.0/grid_t<2>::m_points.size());

    primitive_zone_t<2>::m_volume_weights_unrefined = m_volume_weights;
    
    primitive_zone_t<2>::m_points_count_unrefined = grid_t<2>::m_points.size();

    
    for (unsigned i = 0; i < 2; i ++)
	primitive_zone_t<2>::m_reciprocal_basis_matrix.col(i) = primitive_zone_t<2>::m_reciprocal_basis[i];
    
    primitive_zone_t<2>::m_reciprocal_basis_matrix_inv = primitive_zone_t<2>::m_reciprocal_basis_matrix.inverse();
    
    populate_search_bins();
}


coord_t<2> brillouin_zone_2D_t::get_orthogonal_vec_2D(coord_t<2> vec)
{
    return {vec(1), - vec(0)};
}

coord_t<2> brillouin_zone_2D_t::find_intercept(coord_t<2> a, coord_t<2> c, coord_t<2> b, coord_t<2> d)
{
    //find intercept of the (braag) lines given by
    // (x(s), y(s)) = c + s*a
    // (x'(s), y'(s)) = d + s*b 

    double determinant = b(1) * a(0) - b(0) * a(1);

    if (determinant == 0){
	throw std::invalid_argument("No intercept of parallel lines!"); 
    }
    double s = (c(1) * a(0) - c(0) * a(1) + d(0) * a(1) - d(1) * a(0))/determinant;

    coord_t<2> intercept = d + b * s;
    return intercept;
}

std::tuple<coord_t<2>, coord_t<2> > brillouin_zone_2D_t::get_braag_line(coord_t<2> vec)
{
    // the line bisecting a (reciprocal) lattice vector
    // (x(s), y(s)) = s * a + c
    coord_t<2> c = vec*0.5;
    coord_t<2> a = get_orthogonal_vec_2D(vec);
    return std::make_tuple(a, c);
}


bool brillouin_zone_2D_t::check_if_convex(const std::vector<coord_t<2> > &vertices)
{
    int orientation_counter = 0;
    for (unsigned i = 0; i < vertices.size(); i ++){
	coord_t<2> a = vertices[i];
	coord_t<2> b;
	if (i + 1 == vertices.size())
	    b = vertices[0];
	else
	    b = vertices[i + 1];
	double cross_product_with_z = a(0) * b(1) - a(1) * b(0);
	if (cross_product_with_z > 0)
	    orientation_counter += 1;
	else
	    orientation_counter -= 1;
    }
    return std::abs(orientation_counter) == vertices.size();
}

std::vector<coord_t<2> > brillouin_zone_2D_t::get_wigner_seitz_vertices(std::array<coord_t<2>, 2> basis)
{
    std::vector<coord_t<2> > lattice = span_lattice(basis);
    std::map<double, std::vector<coord_t<2> > > shells = make_shells(lattice);

    auto shells_itr = shells.begin();
    shells_itr++; // nearest neighbohr
    std::vector<coord_t<2> > closest_neighbohrs = shells_itr->second;

    if (closest_neighbohrs.size() == 2){ // in case of a non-homogeneous lattice, one really needs also take the next nearest neighbors to form the braag lines which intersection will form  a wigner-seitz cell
	shells_itr++;
	std::vector<coord_t<2> > next_nearest_neighbohrs = shells_itr->second;
        closest_neighbohrs.insert(closest_neighbohrs.end(), next_nearest_neighbohrs.begin(), next_nearest_neighbohrs.end());
    }

    std::vector<unsigned> vertices_indices;
    for (unsigned i = 0; i < closest_neighbohrs.size(); i++)
	vertices_indices.push_back(i);
    
    std::vector<std::tuple<std::vector<unsigned>, int> > vertex_idx_permutations = make_permutations(vertices_indices);
    
    // try to find the permutation which forms an oriented convex polygon
    for (std::tuple<std::vector<unsigned>, int> &idxes_with_signs: vertex_idx_permutations){
	auto &[idxes, permutation_sign] = idxes_with_signs; // we don't need the permutation sign
	std::vector<coord_t<2> > permuted_closest_neighbohrs;
	for (unsigned idx: idxes)
	    permuted_closest_neighbohrs.push_back(closest_neighbohrs[idx]);

	if (check_if_convex(permuted_closest_neighbohrs)){
	    // pick this permutation
	    closest_neighbohrs = permuted_closest_neighbohrs;
	    break;
	}
    }
    
    // form the braag lines and find their intersections (these will be the vertices of the wigner-seitz cell)
    std::vector<coord_t<2> > vertices;
    for (unsigned i = 0; i < closest_neighbohrs.size(); i ++){
	unsigned i_minus_1 = (i == 0) ? closest_neighbohrs.size()-1 : i - 1;
        
	auto [a, c] = get_braag_line(closest_neighbohrs[i]);
	auto [b, d] = get_braag_line(closest_neighbohrs[i_minus_1]);
	
	coord_t<2> intercept = find_intercept(a, c, b, d);
	vertices.push_back(intercept);
    }

    return vertices;
}

std::tuple<std::vector<coord_t<2> >, double> brillouin_zone_2D_t::span_patched_BZ_wigner_seitz(std::array<coord_t<2>, 2> reciprocal_basis, std::array<coord_t<2>, 2> bz_grid_basis, unsigned points_per_momentum_cell_edge){
    // This method has been tested thoroughly only for a regular lattice.
    std::vector<coord_t<2> > vertices = get_wigner_seitz_vertices(reciprocal_basis);

    auto edges = form_mesh_edges_from_vertices(vertices, points_per_momentum_cell_edge);
    
    // make a flattened edges array
    std::vector<coord_t<2> > edges_flat;
    for (auto &edge : edges){ 
	edges_flat.insert(edges_flat.end(), edge.begin(), edge.end());
    }
    if (edges_flat.size() % 2 != 0)
	throw std::invalid_argument("Boundary should have an even number of points!");

    // since the BZ is a torus, half the boundary points are repeated points
    edges_flat.resize(int(edges_flat.size()/2));
    edges_flat.erase(edges_flat.begin());
    

    std::vector<coord_t<2> > BZ;
    
    // append the boundary points of the mesh to the BZ
    BZ.insert(BZ.end(), edges_flat.begin(), edges_flat.end());

    double interpoint_spacing = std::sqrt((edges[0][0] - edges[0][1]).dot(edges[0][0] - edges[0][1]));

    std::vector<coord_t<2> > BZ_insides = span_mesh_strictly_in_polygon(vertices, bz_grid_basis, interpoint_spacing);
    BZ.insert(BZ.end(), BZ_insides.begin(), BZ_insides.end());

    return std::make_tuple(BZ, interpoint_spacing);
}

bool brillouin_zone_2D_t::is_inside_polygon(const std::vector<coord_t<2> > &vertices, coord_t<2> pt)
{
    // adopted from: https://wrf.ecse.rpi.edu/Research/Short_Notes/pnpoly.html   
    int j = 0;
    bool c = false;
    for (int i = 0; i < vertices.size()-1; i ++){
	j = i + 1;
	if ( ( ( vertices[i](1) > pt(1) ) != ( vertices[j](1) > pt(1) ) ) && ( pt(0) < ( vertices[j](0) - vertices[i](0) ) * ( pt(1) - vertices[i](1) )/ ( vertices[j](1) - vertices[i](1) ) + vertices[i](0)) )
	    c = !c;
    }
    return c;
}

bool brillouin_zone_2D_t::is_strictly_inside_polygon(const std::vector<coord_t<2> > &vertices, coord_t<2> pt)
{
    coord_t<2> pt1 = {pt(0) + epsilon, pt(1)+epsilon};
    coord_t<2> pt2 = {pt(0) - epsilon, pt(1)+epsilon};
    coord_t<2> pt3 = {pt(0) + epsilon, pt(1)-epsilon};
    coord_t<2> pt4 = {pt(0) - epsilon, pt(1)-epsilon};
    return is_inside_polygon(vertices, pt1) && is_inside_polygon(vertices, pt2) && is_inside_polygon(vertices, pt3) && is_inside_polygon(vertices, pt4);
}

std::vector<std::vector<coord_t<2> > > brillouin_zone_2D_t::form_mesh_edges_from_vertices(const std::vector<coord_t<2> > &vertices, unsigned points_per_momentum_cell_edge)
{
    std::vector<std::vector<coord_t<2> > > edges;
    
    // make perimeter
    for (unsigned i = 0; i < vertices.size(); i ++){
	unsigned i_minus_1 = (i == 0) ? vertices.size()-1 : i - 1;
	coord_t<2> spacing_vec = vertices[i] - vertices[i_minus_1];
	double polygon_side_distance = std::sqrt(spacing_vec.dot(spacing_vec) );
	
	int N = points_per_momentum_cell_edge;
	std::vector<coord_t<2> > pts;
	for (unsigned j = 0; j < N; j++){
	    double s = static_cast<double>(j)/static_cast<double>(N);
	    coord_t<2> new_pt = vertices[i]*s + vertices[i_minus_1]*(1.0-s);
	    pts.push_back(new_pt);
	}
	edges.push_back(pts);
    }

    return edges;
}


std::vector<coord_t<2> > brillouin_zone_2D_t::span_mesh_strictly_in_polygon(const std::vector<coord_t<2> > &vertices, const std::array<coord_t<2>, 2> &normalised_basis, double interpoint_spacing)
{
    std::vector<coord_t<2> > span;
    
    std::vector<coord_t<2> > polygon_vertices = vertices;
    // augment the vertices array by repeating its first point towards the end. Necessary for the call to is_inside_polygon that will follow to work properly.
    polygon_vertices.push_back(vertices[0]);

    coord_t<2> basis_1 = normalised_basis[0] * interpoint_spacing;
    coord_t<2> basis_2 = normalised_basis[1] * interpoint_spacing;

    int range = 100;
    for (int b = -range; b <= range; b++)
    for (int a = -range; a <= range; a++){
	coord_t<2> r = basis_1 * a + basis_2 * b;
	if (is_strictly_inside_polygon(polygon_vertices, r))
	    span.push_back(r);
    }

    return span;
}

