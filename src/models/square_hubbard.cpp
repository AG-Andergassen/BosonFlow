#include <models/square_hubbard.h>
#include <primitive_zone.h>
#include <params_physical.h>

void SquareHubbard::Init()
{
    // basis for square lattice
    std::array<coord_t<2>, 2> basis;
    basis[0] = {1, 0};
    basis[1] = {0, 1};

    // define point group of square lattice
    std::vector<std::vector<matrix_t<2>> > C4v;

    {
	matrix_t<2> e, c4, c4inv, c2, sig_x, sig_y, sig_d, sig_dp;
    
	e << 1, 0,
	    0, 1;
    
	c4 << 0, 1,
	    -1, 0;

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
	std::vector<std::string> names = {"e", "c4", "c4inv", "c2", "sig_x", "sig_y", "sig_d", "sig_dp"};
    }

    // define character table of square lattice
    std::vector<std::vector<int> > C4v_character_table = 
	{ // [e][c4][c2][sig_x][sig_d]
	    { 1,  1,  1,   1  ,   1  }, // A1
	    { 1,  1,  1,  -1  ,  -1  }, // B2
	    { 1, -1,  1,   1  ,  -1  }, // B1
	    { 1, -1,  1,  -1  ,   1  }, // B2
	    { 2,  0, -2,   0  ,   0  }  // E
	};

    // special points
    std::map<std::string, coord_t<2> > special_points_coords;
    special_points_coords["idx_00"] = {0, 0};
    special_points_coords["idx_pipi"] = {M_PI, M_PI};
    special_points_coords["idx_Gamma"] = {0, 0};
    special_points_coords["idx_X"] =  {M_PI, 0};
    special_points_coords["idx_M"] =  {M_PI, M_PI};
    special_points_coords["idx_VanHove"] = {M_PI, M_PI};

    // special paths
    std::map<std::string, std::vector<coord_t<2> > > special_paths_coords;
    special_paths_coords["path_fermi_surface"] = {{0, M_PI}, {M_PI, 0}, {0, -M_PI}, {-M_PI, 0}};
    special_paths_coords["path_Gamma_X_M"] = {{0, 0}, {M_PI, 0}, {M_PI, M_PI}};
    special_paths_coords["path_Gamma_M_X"] = {{0, 0}, {M_PI, M_PI}, {M_PI, 0}};
	
    Hubbard<2>::Init(basis, C4v, special_points_coords, special_paths_coords, SquareHubbardParams::REFINE_AT_POINTS);

    std::cout << "Precomputing E(p) for the model..." << std::endl;
    PrecomputeE();    
}

double SquareHubbard::E_of_p_coord(const coord_t<2> p_coord)
{
    const double px = p_coord(0), py = p_coord(1);
    
    const double cos_px = std::cos(px), cos_py = std::cos(py);
    
    return - 2.0 * ( cos_px + cos_py ) - 4.0 * (cos_px * cos_py ) * SquareHubbardParams::T_PRIME - SquareHubbardParams::MU;
}

void SquareHubbard::PrecomputeE()
{
    PrecomputedE().resize(GetFineMomentaCount());

    for (unsigned p_idx = 0; p_idx < GetFineMomentaCount(); p_idx ++){
	coord_t<2> p = FineMomentumGrid().ptr->m_points[p_idx];
	PrecomputedE()[p_idx] = E_of_p_coord(p);
    }
}


double SquareHubbard::E(const int idx_p, const int idx_w)
{
    return PrecomputedE()[idx_p];
}

double SquareHubbard::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       )
{
# ifdef F_NONZERO
    double val = 0;
    val += -SquareHubbardParams::UINT;
    int idx_K = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));
    
    const complex_vector_t<2> e_iK = FineMomentumGrid().exp_ik_idx_map[idx_K];
    
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    const double cos_Ky  = e_iK(1).real();
    
    val += - 2 * SquareHubbardParams::VINT * ( cos_Kx + cos_Ky);
    return val;
#endif
    return -SquareHubbardParams::UINT;
}



double SquareHubbard::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * SquareHubbardParams::UINT;
    else
	return 0.;
}

double SquareHubbard::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * SquareHubbardParams::UINT;
    else
	return 0.;
}


double SquareHubbard::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){

	const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
	// real part of e^{ip} is cos p
	const double cos_Kx  = e_iK(0).real();
	const double cos_Ky  = e_iK(1).real();
        
	val += - MomentumGrid().ptr->get_volume() * SquareHubbardParams::UINT;
        val += - MomentumGrid().ptr->get_volume() * 4 * SquareHubbardParams::VINT * (cos_Kx + cos_Ky);
	return val;
    }else
	return 0.;
}


double SquareHubbard::F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (FORMFACTOR_SHELL_COUNT == 1.5) // projection into antisymmetric d-wave ff would just give zero.
	return 0;
    if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 4){ // if first shell form factors. TODO: Warning, this function makes no sense with a d-wave form factor!!
    	return - MomentumGrid().ptr->get_volume() * SquareHubbardParams::VINT;
    }
    return 0.0;
}

double SquareHubbard::F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (FORMFACTOR_SHELL_COUNT == 1.5)
	return 0;
    if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 4){
	return MomentumGrid().ptr->get_volume() * SquareHubbardParams::VINT;
    }
    return 0.0;
}


double SquareHubbard::F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (FORMFACTOR_SHELL_COUNT == 1.5)
	return 0;
    if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 4){
	return MomentumGrid().ptr->get_volume() * SquareHubbardParams::VINT;
    }
    return 0.0;
}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double SquareHubbard::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * SquareHubbardParams::UINT;
    else
	return 0.;
}

double SquareHubbard::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return -SquareHubbardParams::UINT;
}



double SquareHubbard::vertex_DC(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out
			       )
{
    double val = 0;
    val -= -2.0*SquareHubbardParams::UINT;
    return val;
}


double SquareHubbard::vertex_DC_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (idx_m == 0 && idx_mp == 0){
        return -2.0 * SquareHubbard::vertex_local_part_bare(idx_W, idx_w, idx_m, idx_wp, idx_mp);
    
    }else{
	return 0;
    }
}

double SquareHubbard::vertex_DC_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (idx_m == 0 && idx_mp == 0){
        return -2.0 * SquareHubbard::vertex_local_part_bare(idx_W, idx_w, idx_m, idx_wp, idx_mp);
    }else{
	return 0;
    }
}

double SquareHubbard::vertex_DC_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (idx_m == 0 && idx_mp == 0){
        return 2.0 * SquareHubbard::vertex_local_part_bare(idx_W, idx_w, idx_m, idx_wp, idx_mp);
    }else{
	return 0;
    }
}


/*
double SquareHubbard::vertex_real_space_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_r1, const int idx_r2, const int idx_r3, const int idx_r4, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    if (idx_r1 == idx_r2 && idx_r1 == idx_r3 && idx_r1 == idx_r4){
	return -SquareHubbardParams::UINT;
    }
    return 0;
}
*/
