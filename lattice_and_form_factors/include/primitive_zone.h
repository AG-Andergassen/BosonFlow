#pragma once
#include <momentum_space.h>


template <unsigned dim>
class primitive_zone_t : public momentum_space_t<dim>
{
 public:
    primitive_zone_t(basis_t<dim> reciprocal_basis, unsigned points_per_momentum_cell_edge);
  
    std::tuple<unsigned, double> get_idx_from_coord(coord_t<dim> coord) const override;

    std::tuple<unsigned, double> get_idx_from_coord_unrefined(coord_t<dim> coord) const;

    std::vector<unsigned> get_path_indices(std::vector<coord_t<dim> > points, bool is_closed = false);

    std::vector<unsigned> get_path_indices_unrefined(std::vector<coord_t<dim> > points, bool is_closed = false);


    // refine this grid based on a regularly spaced fine grid
    void refine_from_finer_grid(momentum_space_t<dim> *fine_grid_ptr, coord_t<dim> center, float refine_percent);

    void refine_from_finer_grid_along_path(momentum_space_t<dim> *fine_grid_ptr, std::vector<coord_t<dim>> points, double refine_percent, bool is_closed_path);
  
    unsigned m_points_per_momentum_cell_edge;

    double m_momentum_mesh_min_interpoint_spacing;
    basis_t<dim> m_reciprocal_basis;

    std::vector<double> m_volume_weights; // what fraction of "volume" does every point/patch represent

    std::vector<double> m_volume_weights_unrefined; // what fraction of "volume" does every point/patch represent, before any refinement of grid

    // number of grid points that represent the unrefined grid
    unsigned m_points_count_unrefined;

    matrix_t<dim> m_reciprocal_basis_matrix;
    matrix_t<dim> m_reciprocal_basis_matrix_inv;


 protected:
    primitive_zone_t();

    std::tuple<unsigned, double> get_idx_from_coord(coord_t<dim> coord, unsigned points_count) const override;

    std::vector<unsigned> get_path_indices(std::vector<coord_t<dim> > points, bool is_closed, unsigned grid_points_count);

};


#include "../src/primitive_zone.tpp"
