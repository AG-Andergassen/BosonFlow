#pragma once 

#include <tuple>
#include <vector>
#include <string>
#include <map>
#include <Eigen/Core>
#include <iostream>
#include <Eigen/LU>
#include <unordered_map>
#include <type_traits>
#include <boost/geometry.hpp>
#include <boost/geometry/index/rtree.hpp>
#include <nanoflann.hpp>

template <unsigned dim1, unsigned dim2>
    //template <typename std::enable_if<dim1 != 0>::type>
struct matrix_type
{
    using T = Eigen::Matrix<double, dim1, dim2>;

    constexpr bool func(){
	static_assert(dim1 != 0, "can't call this with dim1 = 0...");
    
	return true;
    }
};

template <>
struct matrix_type<0, 1>
{
    using T = Eigen::Matrix<double, 1, 1>;
};

template <>
struct matrix_type<0, 0>
{
    using T = Eigen::Matrix<double, 1, 1>;
};


template <unsigned dim>
using vector_t = typename matrix_type<dim, 1>::T;


template <unsigned dim>
struct complex_vector_type
{
    using T = Eigen::Matrix<std::complex<double>, dim, 1>;
};

template <>
struct complex_vector_type<0>
{
    using T = Eigen::Matrix<std::complex<double>, 1, 1>;
};

template <unsigned dim>
using complex_vector_t = typename complex_vector_type<dim>::T;

template <unsigned dim>
class coord_t : public vector_t<dim>
{
   using vector_t<dim>::vector_t;

 public:
   // needed for boost's vornoi diagrams
    double get(int i) const {return this->operator()(i);}

    void set(int i, double value) const {this->operator()(i) = value;}

    static std::size_t dimension(){ return dim; }

//double &operator[](unsigned i ) {return this->operator()(i);}
};


//using coord_t = vector_t<dim>;

template <unsigned dim>
using rtree_point_t= boost::geometry::model::point<double, dim, boost::geometry::cs::cartesian>; 

// wrapper for nanoflann to construct a kdtree for fast nearest neighbor search
template <unsigned dim>
struct GridNanoflannAdaptor {
    const std::vector<coord_t<dim> >& points;
    
    // Constructor
    GridNanoflannAdaptor(const std::vector<coord_t<dim> >& pts) : points(pts) {}
    
    // Returns the number of data points
    inline size_t kdtree_get_point_count() const {
	return points.size();
    }
    
    // Returns a specific dimension of a point
    inline double kdtree_get_pt(const size_t idx, const size_t i) const {
	return points[idx](i);
    }
    
    // Optional bounding-box computation
    template <class BBOX>
    bool kdtree_get_bbox(BBOX& /*bb*/) const {
        return false; // Return true if bounding box computation is enabled
    }
};


template <unsigned dim>
struct int_vector_type
{
    using T = Eigen::Matrix<int, dim, 1>;
};

template <>
struct int_vector_type<0>
{
    using T = Eigen::Matrix<int, 1, 1>;
};

template <unsigned dim>
using coord_units_t = typename int_vector_type<dim>::T;


template <unsigned dim>
using matrix_t = typename matrix_type<dim, dim>::T;


template <unsigned dim>
struct basis_type
{
    using T = std::array<coord_t<dim>, dim>;
    using T_coeffs = std::vector< std::array<int, dim> >;
};

template <>
struct basis_type<0>
{
    using T = std::array<coord_t<0>, 1>;
    using T_coeffs = std::vector< std::array<int, 1> >;
};

template <unsigned dim>
using basis_t = typename basis_type<dim>::T;

template <unsigned dim>
using int_basis_coeffs_t = typename basis_type<dim>::T_coeffs; 

template <class T>
class CoordHash;

template<unsigned dim>
class CoordHash< Eigen::Matrix<double, dim, 1> >
{
 public:
    std::size_t operator()(Eigen::Matrix<double, dim, 1> const& s) const 
	{
	    std::size_t h = (unsigned&)(s(0));
	    for (unsigned i = 1; i <= s.rows(); i ++)
		h = h ^ ((unsigned&)(s(i)) << 1);
	    return h;
	}
};

/*
namespace boost {
    namespace polygon {

template <>
    struct geometry_concept<coord_t<2>> {
    typedef point_concept type;  // Tells Boost this is a point
 };

template <>
    struct point_traits<coord_t<2>> {
    typedef double coordinate_type;

    static inline coordinate_type get(const coord_t<2>& point, int index) {
        return point.get(index);  // Access x or y via index
    }
 };

    }  // namespace polygon
}  // namespace boost
*/

//typedef boost::polygon::voronoi_diagram<double> VoronoiDiagram;

template <unsigned dim>
class grid_t
{
 public: 


    grid_t(const std::vector<coord_t<dim> > &points);
    
    grid_t(basis_t<dim> basis, unsigned range);
    
    std::vector<double> get_points_as_std_vectors() const;

    std::vector<coord_t<dim> > span_lattice(basis_t<dim> basis, unsigned range = 5 ) const;
    
    std::vector<double> compute_norms_squared(const std::vector<coord_t<dim> > &vecs) const;

    std::map<double, std::vector<coord_t<dim> > > make_shells(const std::vector<coord_t<dim> > &span, const std::vector<double> &norms_squared, std::map<double, std::vector<unsigned> > &shells_indices, double max_norm_squared = std::numeric_limits<double>::max()) const;

    std::map<double, std::vector<coord_t<dim> > > make_shells(const std::vector<coord_t<dim> > &span, double max_norm_squared = std::numeric_limits<double>::max()) const;

    std::map<double, std::vector<coord_t<dim> > > make_shells(const std::vector<coord_t<dim> > &span, std::map<double, std::vector<unsigned> > &shells_indices, double max_norm_squared = std::numeric_limits<double>::max()) const;

    std::map<double, std::vector<coord_t<dim> > > make_shells(const std::vector<coord_t<dim> > &span, std::map<double, std::vector<unsigned> > &shells_indices, const coord_t<dim> center, double max_norm_squared = std::numeric_limits<double>::max()) const;

    // get volume made by a (dim) number of vectors
    double compute_enclosed_volume(basis_t<dim> vecs) const;

    // returns all permutations of elements in a given list of (unsigned integers) together with the sign of the permutation
    std::vector<std::tuple<std::vector<unsigned>, int> > make_permutations(const std::vector<unsigned> &ls) const;

    // returns all dim-dimensional integer vectors in which entries are inside [-range, range]. It's done in row major order (to be compatable with FFTW format). i.e. if the i-th index is the origin (0, 0, 0 ... 0, 0), then the (i+1)-th index is (0, 0, 0 ... 0, 1)
    int_basis_coeffs_t<dim> make_integer_coeffs(unsigned range, bool is_closed_interval = true) const;

    // get a mapping between coordinates lying in a grid and their indices
    std::unordered_map<coord_t<dim>, unsigned, CoordHash<coord_t<dim> > > get_coord_to_idx_map(const std::vector<coord_t<dim> > &mesh);
    
    
    // get the index corresponding to a coordinate in a grid. Override e.g. to add periodicity
    virtual std::tuple<unsigned, double>  get_idx_from_coord(coord_t<dim> coord) const;

    // for debugging purposes, out put m_points as a python list.
    std::string make_debug_grid_python_list(std::string var_name);


    std::vector<coord_t<dim> > m_points; // the grid points in the canonical basis
    //std::vector<coord_units_t<dim> > m_points_units; // the grid points in the lattice basis (integer coords)
    unsigned m_size;

    const double epsilon = 1e-10;

    void populate_search_bins();
    unsigned get_coord_search_bin_idx(const coord_t<dim>& coord) const;


    grid_t(const grid_t& other_grid) : m_points(other_grid.m_points), m_size(other_grid.m_size), m_grid_kdtree_adaptor(m_points), m_points_kdtree((dim == 0) ? 1 : dim, m_grid_kdtree_adaptor, nanoflann::KDTreeSingleIndexAdaptorParams(10))
    {
    };

protected:
    // get the index corresponding to a coordinate in a grid, for only a given number of points. Override e.g. to add periodicity

    virtual std::tuple<unsigned, double>  get_idx_from_coord(coord_t<dim> coord, unsigned points_count) const;

    std::vector<std::vector<unsigned>> m_search_bins_of_coord_indices;
    double m_min_coord;  // Minimum value of the 0th coordinate in the grid
    double m_max_coord;  // Maximum value of the 0th coordinate in the grid

    unsigned m_num_bins_per_dim;

    std::vector<double> m_min_coords;
    std::vector<double> m_max_coords;

    // grid of indices
    std::vector<std::uint16_t> m_dense_search_container;

    //boost::geometry::index::rtree<std::pair<coord_t<dim>, unsigned>, boost::geometry::index::kd_tree<dim>> rtree({});
    
    //static constexpr std::size_t adjusted_dim = (dim == 0) ? 1 : dim;
    //static constexpr std::size_t adjusted_dim = std::conditional_t<(dim == 3), std::integral_constant<std::size_t, dim + 1>, std::integral_constant<std::size_t, dim>>::value;

    boost::geometry::index::rtree<std::pair<rtree_point_t< (dim == 0) ? 1 : dim >, unsigned>, boost::geometry::index::quadratic<16>> m_points_rtree;


    using KDTree = nanoflann::KDTreeSingleIndexAdaptor<
	nanoflann::L2_Simple_Adaptor<double, GridNanoflannAdaptor<dim> >,
        GridNanoflannAdaptor<dim>,
	((dim == 0) ? 1 : dim) >;


    GridNanoflannAdaptor<dim> m_grid_kdtree_adaptor;
    KDTree m_points_kdtree;

};

#include <grid.cpp>
