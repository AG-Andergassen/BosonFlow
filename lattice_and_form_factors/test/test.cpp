#include <real_space_lattice.h>
#include <brillouin_zone_2D.h>
#include <delta_function_vector.h>
#include <lattice_point_group.h>
#include <bond_form_factor_container.h>
#include <cmath>
#include <iostream>

void test_delta_function_vector()
{
    delta_function_vector_t<2> delta({{0, 0} }, {1});
    std::cout << (delta + delta).get_pretty_string() << std::endl;
    std::cout << delta.evaluate({0, 0}) << "  " << delta.evaluate({0, 0.1}) << std::endl;
}

void test_bond_form_factors()
{
    std::array<coord_t<2>, 2> basis;
    basis[0] = {1, 0};
    basis[1] = {0, 1};

    real_space_lattice_t<2> square_lattice(basis, 3);

    std::vector<std::vector<matrix_t<2>> > C4v;
    std::vector<std::string> names;

    {
	matrix_t<2> e, c4, c4inv, c2, sig_x, sig_y, sig_d, sig_dp;
    
	e << 1, 0,
	    0, 1;
    
	c4 << 0, 1,
	    -1, 0;
    
	//std::cout << c4;

	c4inv = c4.inverse();
	c2 = c4 * c4;
    
	sig_x << -1, 0,
	    0, 1;
	sig_y << 1, 0,
	    0, -1;
	sig_d << 0, 1,
	    1, 0;
	sig_dp << 0, -1,
	    -1, 0;

	C4v = {{e}, {c4, c4inv}, {c2}, {sig_x, sig_y}, {sig_d, sig_dp}};
	names = {"e", "c4", "c4inv", "c2", "sig_x", "sig_y", "sig_d", "sig_dp"};
    }

    std::vector<std::vector<int> > C4v_character_table = 
	{ // [e][c4][c2][sig_x][sig_d]
	    { 1,  1,  1,   1  ,   1  }, // A1
	    { 1,  1,  1,  -1  ,  -1  }, // B2
	    { 1, -1,  1,   1  ,  -1  }, // B1
	    { 1, -1,  1,  -1  ,   1  }, // B2
	    { 2,  0, -2,   0  ,   0  }  // E
	};

    
    lattice_point_group_t<2> sq_point_group(C4v, names); 

    primitive_zone_t<2> square_lattice_PZ(square_lattice.make_reciprocal_lattice_basis(), 5);

    //brillouin_zone_2D_t square_lattice_PZ(square_lattice.make_reciprocal_lattice_basis(), square_lattice.m_basis, 5);

    auto [group_in_mom_idx_space ,err] = sq_point_group.make_group_in_momentum_idx_space(square_lattice_PZ, 1e-10);
    std::cout << "Symmetry errors: " << err << std::endl; 

    bond_form_factor_container_t<2> square_bond_ffs(square_lattice, 3);

    auto ff_shells_in_mom_idx_space = square_bond_ffs.make_form_factors_in_momentum_idx_space(square_lattice_PZ);

    std::cout << "Acting on the form factor: " << square_bond_ffs.get(1, 1).get_pretty_string() << std::endl;
    
    std::vector<std::map<shell_ff_idx_t, std::tuple<ff_idx_t, std::string> > > pt_grp = square_bond_ffs.make_point_group_in_shells_form_factor_idx_space(group_in_mom_idx_space, ff_shells_in_mom_idx_space);

    for (unsigned g_idx = 0; g_idx < sq_point_group.m_group_size; g_idx++){
	unsigned result = std::get<0>(pt_grp[g_idx].at(std::make_tuple(1, 1)));
	std::cout << "Action of " << sq_point_group.m_element_names[g_idx] << " gives "; 
	std::cout << square_bond_ffs.get(1, result).get_pretty_string() << " upto operation " << std::get<1>(pt_grp[g_idx].at(std::make_tuple(1, 1))) << std::endl;
    }
}

void test_lattice_symmetries()
{

    std::vector<std::vector<matrix_t<2>> > C4v;

    {
	matrix_t<2> e, c4, c4inv, c2, sig_x, sig_y, sig_d, sig_dp;
    
	e << 1, 0,
	    0, 1;
    
	c4 << 0, 1,
	    -1, 0;
    
	//std::cout << c4;

	c4inv = c4.inverse();
	c2 = c4 * c4;
    
	sig_x << -1, 0,
	    0, 1;
	sig_y << 1, 0,
	    0, -1;
	sig_d << 0, 1,
	    1, 0;
	sig_dp << 0, -1,
	    -1, 0;

	C4v = {{e}, {c4, c4inv}, {c2}, {sig_x, sig_y}, {sig_d, sig_dp}};
    }

    std::vector<std::vector<int> > C4v_character_table = 
	{ // [e][c4][c2][sig_x][sig_d]
	    { 1,  1,  1,   1  ,   1  }, // A1
	    { 1,  1,  1,  -1  ,  -1  }, // B2
	    { 1, -1,  1,   1  ,  -1  }, // B1
	    { 1, -1,  1,  -1  ,   1  }, // B2
	    { 2,  0, -2,   0  ,   0  }  // E
	};

    
    std::array<coord_t<2>, 2> basis;
    basis[0] = {1, 0};
    basis[1] = {0, 1};

    real_space_lattice_t<2> square_lattice(basis, 3);

    primitive_zone_t<2> square_lattice_PZ(square_lattice.make_reciprocal_lattice_basis(), 5);


    lattice_point_group_t<2> sq_point_group(C4v); 

    auto [sq_group_in_mom_idx_space ,err] = sq_point_group.make_group_in_momentum_idx_space(square_lattice_PZ, 1e-10);
    std::cout << "Symmetry errors: " << err << std::endl; 

    
    
    std::vector<std::vector<matrix_t<2>> > C6v;

    {
	matrix_t<2> e, c6, c6inv, c2, c3, c3inv, sig_y, sig_y2, sig_y3, sig_d, sig_d2, sig_d3;

	e << 1, 0,
	    0, 1;
    
	c6 << std::cos(M_PI/3.0) , -std::sin(M_PI/3.0),
	    std::sin(M_PI/3.0), std::cos(M_PI/3.0);
	c6inv = c6.inverse();

	c3 = c6*c6;
	c3inv = c3.inverse();
    
	c2 << std::cos(M_PI), -std::sin(M_PI), 
	    std::sin(M_PI), std::cos(M_PI);

	sig_y << 1, 0, 
	    0, -1;
	sig_y2 = sig_y * c3;
	sig_y3 = sig_y2 * c3;
    
	sig_d = sig_y * c6;
	sig_d2 = sig_d * c3;
	sig_d3 = sig_d2 * c3;

	C6v = {{e}, {c6, c6inv}, {c3, c3inv}, {c2}, {sig_y, sig_y2, sig_y3}, {sig_d, sig_d2, sig_d3}};
    }

    std::array<coord_t<2>, 2> tri_basis;
    tri_basis[0] = {1, 0};
    tri_basis[1] = {-0.5, std::sqrt(3)/2.0};

    real_space_lattice_t<2> triangular_lattice(tri_basis, 3);
    
	//lattice_point_group_t<2> tri_point_group(triangular_lattice, C4v); 

	//err = tri_point_group.make_group_in_momentum_idx_space(1e-10);
	//std::cout << "Symmetry errors: " << err << std::endl; 


    lattice_point_group_t<2> tri_point_group(C6v); 

    brillouin_zone_2D_t tri_lattice_BZ(triangular_lattice.make_reciprocal_lattice_basis(), triangular_lattice.m_basis, 5);

    auto [tri_group_in_mom_idx_space ,err2] = tri_point_group.make_group_in_momentum_idx_space(tri_lattice_BZ, 1e-10);
    std::cout << "Symmetry errors: " << err2 << std::endl; 
    
}

void test_lattice()
{
    std::array<coord_t<1>, 1> basis1;
    basis1[0](0) = 1;

    real_space_lattice_t<1> chain_lattice(basis1, 3);

    primitive_zone_t<1> chain_BZ(chain_lattice.make_reciprocal_lattice_basis(), 9);

    std::cout << "Number of BZ pts (chain): " << chain_BZ.m_points.size() << std::endl;
    for (coord_t<1> p : chain_BZ.m_points){
    	std::cout << "p: " << p(0) << std::endl;
    }


    std::array<coord_t<2>, 2> basis;
    basis[0] = {1, 0};
    basis[1] = {0, 1};

    real_space_lattice_t<2> square_lattice(basis, 3);

    primitive_zone_t<2> square_BZ(square_lattice.make_reciprocal_lattice_basis(), 5);

    auto [neg_k_idxes, neg_k_idxes_error] = square_BZ.precalculate_negative_indices();
    auto [sum_k_idxes, sum_k_idxes_error] = square_BZ.precalculate_sum_of_indices();

    std::cout << "Number of BZ pts (square): " << square_BZ.m_points.size() << std::endl;
    for (unsigned p_idx = 0; p_idx < square_BZ.m_points.size(); p_idx ++){
	coord_t<2> p = square_BZ.m_points[p_idx];
	coord_t<2> minus_p = square_BZ.m_points[neg_k_idxes[p_idx]];
	
    	std::cout << "p: " << p(0) << " " << p(1) << std::endl;
	std::cout << "-p: " << minus_p(0) << " " << minus_p(1) << std::endl;
    }

    {
	coord_t<2> p1 = square_BZ.m_points[12];
	coord_t<2> p2 = square_BZ.m_points[9];
	coord_t<2> sum1 = p1 + p2;
	coord_t<2> sum2 = square_BZ.m_points[sum_k_idxes[12][9]];
	std::cout << sum1(0) << " " << sum1(1)<< " should be equal to " << sum2(0) << " " << sum2(1) << std::endl;
    }
    

    std::array<coord_t<3>, 3> basis3;
    basis3[0] = {1, 0, 0};
    basis3[1] = {0, 1, 0};
    basis3[2] = {0, 0, 1};

    real_space_lattice_t<3> cube_lattice(basis3, 3);

    primitive_zone_t<3> cube_BZ(cube_lattice.make_reciprocal_lattice_basis(), 5);

    std::cout << "Number of BZ pts (cube): " << cube_BZ.m_points.size() << std::endl;
    for (coord_t<3> p : cube_BZ.m_points){
    	std::cout << "p: " << p(0) << " " << p(1) << " " << p(2) << std::endl;
    }


    //lattice_2D_t square_lattice({{1, 0}, {0, 1}}, 1.0);
    std::array<coord_t<2>, 2> tri_basis;
    tri_basis[0] = {1, 0};
    tri_basis[1] = {-0.5, std::sqrt(3)/2.0};

    real_space_lattice_t<2> triangular_lattice(tri_basis, 3);

    brillouin_zone_2D_t triangular_BZ(triangular_lattice.make_reciprocal_lattice_basis(), triangular_lattice.m_basis, 5); 

    std::cout << "Number of BZ pts: " << triangular_BZ.m_points.size() << std::endl;
    for (coord_t<2> p : triangular_BZ.m_points){
    	std::cout << "p: " << p(0) << " " << p(1) << std::endl;
    }

}



int main(int argc, char *argv[])
{
    test_lattice();
    test_delta_function_vector();
    test_lattice_symmetries();
    test_bond_form_factors();
    return 0;
}
