#include <models/triangular_hubbard.h>
#include <params_physical.h>

void TriangularHubbard::Init()
{
    // basis for triangular lattice
    // Remark: The basis vectors in reciprocal space are: b1=(2*M_PI,2*M_PI/sqrt(3)), b2=(0,4*M_PI/sqrt(3))
    std::array<coord_t<2>, 2> basis;
    basis[0] = {1, 0};
    basis[1] = {-0.5, std::sqrt(3)/2.0};

    // define point group of triangular lattice
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

    /*    // define character table of triangular lattice
    std::vector<std::vector<int> > C4v_character_table = 
	{ // [e][c4][c2][sig_x][sig_d]
	    { 1,  1,  1,   1  ,   1  }, // A1
	    { 1,  1,  1,  -1  ,  -1  }, // B2
	    { 1, -1,  1,   1  ,  -1  }, // B1
	    { 1, -1,  1,  -1  ,   1  }, // B2
	    { 2,  0, -2,   0  ,   0  }  // E
	};
    */
    

    // special points
    std::map<std::string, coord_t<2> > special_points_coords;
    special_points_coords["idx_00"] = {0., 0.};
    special_points_coords["idx_Gamma"] = {0., 0.};
    special_points_coords["idx_M"] = {M_PI, M_PI/std::sqrt(3.0)};
    special_points_coords["idx_Mprime"] = {0., 2.*M_PI/std::sqrt(3.0)};
    special_points_coords["idx_Mprimeprime"] = {-M_PI, M_PI/std::sqrt(3.0)};
    special_points_coords["idx_VanHove"] = {M_PI, M_PI/std::sqrt(3.0)};
    special_points_coords["idx_K"] = {2.*M_PI/3., 2.*M_PI/std::sqrt(3.0)};
    special_points_coords["idx_Kprime"] = {4.*M_PI/3., 0.};


    // special paths
    std::map<std::string, std::vector<coord_t<2> > > special_paths_coords;
    special_paths_coords["path_K_Gamma_M_Kprime"] = {{2.*M_PI/3., 2.*M_PI/std::sqrt(3.)}, {0., 0.}, {M_PI, M_PI/std::sqrt(3.)}, {4.*M_PI/3., 0.}};

    Hubbard<2>::Init(basis, C6v, special_points_coords, special_paths_coords, TriangularHubbardParams::REFINE_AT_POINTS);

    std::cout << "Precomputing E(p) for the model..." << std::endl;
    PrecomputeE();

    std::cout << "Precomputing the potential for the model..." << std::endl;
    PrecomputePotential();

}

double TriangularHubbard::E_of_p_coord(const coord_t<2> p_coord)
{
    const double px = p_coord(1), py = p_coord(0);
    
    const double cos_px = std::cos(px), cos_py = std::cos(py);
    
    
    // for t1:
    const double cos_px_o_2 = std::cos(px/2.0);
    const double cos_sqrt_3_py_o_2 = std::cos(std::sqrt(3.0) * py/2.0 );

    // for t2
    const double cos_3_px_o_2 = std::cos(3.0 * px/2.0 );
    const double cos_sqrt_3_py = std::cos(std::sqrt(3.0) * py );


    // for t3
    const double cos_2_px = std::cos(2.0 * px);

    //const double t1 = -2.5 / 0.5, t2 = 1, t3 = 0.25/0.5;
    const double t1 = 1.0, t2 = TriangularHubbardParams::T_PRIME, t3 = TriangularHubbardParams::T_PRIME_PRIME; 

    const double mu = TriangularHubbardParams::MU;

    return - 2.0 * t1 * ( cos_px + 2*cos_px_o_2 * cos_sqrt_3_py_o_2 )
	  - 2.0 * t2 * (2*cos_3_px_o_2 * cos_sqrt_3_py_o_2 + cos_sqrt_3_py)
	+ 2.0 * t3 * (cos_2_px - 2*cos_px * cos_sqrt_3_py) - mu;

}

void TriangularHubbard::PrecomputeE()
{
    PrecomputedE().resize(GetFineMomentaCount());

    for (unsigned p_idx = 0; p_idx < GetFineMomentaCount(); p_idx ++){
	coord_t<2> p = FineMomentumGrid().ptr->m_points[p_idx];
	PrecomputedE()[p_idx] = E_of_p_coord(p);
    }
}

double TriangularHubbard::Potential_of_K_coord(const coord_t<2> K_coord)
{
    const double Kx = K_coord(1), Ky = K_coord(0);
    
    const double cos_Kx = std::cos(Kx), cos_Ky = std::cos(Ky);    
    
    // for V:
    const double cos_Kx_o_2 = std::cos(Kx/2.0);
    const double cos_sqrt_3_Ky_o_2 = std::cos(std::sqrt(3.0) * Ky/2.0 );

    // for W
    const double cos_3_Kx_o_2 = std::cos(3.0 * Kx/2.0 );
    const double cos_sqrt_3_Ky = std::cos(std::sqrt(3.0) * Ky );


    // for WW
    const double cos_2_Kx = std::cos(2.0 * Kx);

    const double U = TriangularHubbardParams::UINT, V = TriangularHubbardParams::VINT, W = TriangularHubbardParams::WINT, WW = TriangularHubbardParams::WWINT; 

    return + 2.0 * V * ( cos_Kx + 2*cos_Kx_o_2 * cos_sqrt_3_Ky_o_2 )
	  + 2.0 * W * (2*cos_3_Kx_o_2 * cos_sqrt_3_Ky_o_2 + cos_sqrt_3_Ky)
	- 2.0 * WW * (cos_2_Kx - 2*cos_Kx * cos_sqrt_3_Ky);
}

void TriangularHubbard::PrecomputePotential()
{
    PrecomputedPotential().resize(GetFineMomentaCount());

    for (unsigned K_idx = 0; K_idx < GetFineMomentaCount(); K_idx ++){
	coord_t<2> K = FineMomentumGrid().ptr->m_points[K_idx];
	PrecomputedPotential()[K_idx] = Potential_of_K_coord(K);
    }
}

double TriangularHubbard::E(const int idx_p, const int idx_w)
{
    return PrecomputedE()[idx_p];
}

double TriangularHubbard::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       )
{
# ifdef F_NONZERO
    int idx_K = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));

    return -TriangularHubbardParams::UINT-PrecomputedPotential()[idx_K];
#endif
    return -TriangularHubbardParams::UINT;
}


double TriangularHubbard::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * TriangularHubbardParams::UINT;
    else
	return 0.;
}

double TriangularHubbard::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * TriangularHubbardParams::UINT;
    else
	return 0.;
}


double TriangularHubbard::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0){
	return - MomentumGrid().ptr->get_volume() * ( 2 * PrecomputedPotential()[idx_K] + TriangularHubbardParams::UINT);
    }else
	return 0.;
}


double TriangularHubbard::F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 6){ // if first shell form factors
	return - MomentumGrid().ptr->get_volume() * TriangularHubbardParams::VINT;
    }
    return 0.0;
}

double TriangularHubbard::F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 6){
	return MomentumGrid().ptr->get_volume() * TriangularHubbardParams::VINT;
    }
    return 0.0;
}


double TriangularHubbard::F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 6){
	return MomentumGrid().ptr->get_volume() * TriangularHubbardParams::VINT;
    }
    return 0.0;
}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double TriangularHubbard::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * TriangularHubbardParams::UINT;
    else
	return 0.;
}

double TriangularHubbard::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return -TriangularHubbardParams::UINT;
}
