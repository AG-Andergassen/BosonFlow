#pragma once
#include <primitive_zone.h>


class brillouin_zone_2D_t : public primitive_zone_t<2>
{
 public:
    // for the bz_grid_basis, pass the real lattice basis for the BZ to have the same symmetries as the original lattice
    brillouin_zone_2D_t(std::array<coord_t<2>, 2> reciprocal_basis, std::array<coord_t<2>, 2> bz_grid_basis, unsigned points_per_momentum_cell_edge);

    // get a 2D vector orthogonal to a 2D vector 
    coord_t<2> get_orthogonal_vec_2D(coord_t<2> vec);

    // find the the point of intersection of the two lines given by (x(s), y(s)) = c + s*a  and (x'(s), y'(s)) = d + s*b
    coord_t<2> find_intercept(coord_t<2> a, coord_t<2> c, coord_t<2> b, coord_t<2> d);

    // get the "braag" line corresponding to a lattice vector. Returns a tuple {a, c} for a line given by (x(s), y(s)) = c + s*a 
    std::tuple<coord_t<2>, coord_t<2> > get_braag_line(coord_t<2> vec);

    // check if an ordered set of vertices form an oriented convex polygon
    bool check_if_convex(const std::vector<coord_t<2> > &vertices);

    // span the unique unit cell in the reciprocal space formed by braag planes; the Brillouin zone
    std::tuple<std::vector<coord_t<2> >, double> span_patched_BZ_wigner_seitz(std::array<coord_t<2>, 2> reciprocal_basis, std::array<coord_t<2>, 2> bz_grid_basis, unsigned points_per_momentum_cell_edge);
    
    // check if point lies within a polygon formed by vertices
    bool is_inside_polygon(const std::vector<coord_t<2> > &vertices, coord_t<2> pt);

    // check if point lies within a polygon formed by vertices (excluding boundary)
    bool is_strictly_inside_polygon(const std::vector<coord_t<2> > &vertices, coord_t<2> pt);

    std::vector<coord_t<2> > get_wigner_seitz_vertices(std::array<coord_t<2>, 2> basis);

    std::vector<std::vector<coord_t<2> > > form_mesh_edges_from_vertices(const std::vector<coord_t<2> > &vertices, unsigned points_per_momentum_cell_edge);

    // span a mesh for a polygon formed by given vertices
    std::vector<coord_t<2> > span_mesh_strictly_in_polygon(const std::vector<coord_t<2> > &vertices, const std::array<coord_t<2>, 2> &normalised_basis, double interpoint_spacing);
};
