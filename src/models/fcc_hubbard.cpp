#include <models/fcc_hubbard.h>
#include <primitive_zone.h>
#include <params_physical.h>

void FCCHubbard::Init()
{
#ifdef F_NONZERO
    throw std::invalid_argument("FCC Hubbard model with extended interactions not implemented (turn off the F_NONZERO flag).");
#endif

    // basis for square lattice
    std::array<coord_t<3>, 3> basis;
    basis[0] = {1, 0, 0};
    basis[1] = {0, 1, 0};
    basis[2] = {0, 0, 1};

    // define point group of square lattice
    std::vector<std::vector<matrix_t<3>> > C4v;

    {
	matrix_t<3> e, c4x, c4y, c4z, c4xinv, c4yinv, c4zinv, c2x, c2y, c2z, sig_x, sig_y, sig_z, sig_dx, sig_dy, sig_dz, sig_dpx, sig_dpy, sig_dpz;
    
	e << 1, 0, 0,
	    0, 1, 0,
	    0, 0, 1;
    
	c4z << 0, 1, 0,
	    -1, 0, 0,
	    0, 0, 1;
	c4x << 1, 0, 0,
	    0, 0, 1,
	    0 ,-1, 0;
	c4y << 0, 0, 1,
	    0, 1, 0,
	    -1 ,0, 0;

	c4xinv = c4x.inverse();
	c4yinv = c4y.inverse();
	c4zinv = c4z.inverse();
	
	c2x = c4x * c4x;
	c2y = c4y * c4y;
	c2z = c4z * c4z;

	sig_x << -1, 0, 0,
	    0, 1, 0,
	    0, 0, 1;
	sig_y << 1, 0, 0,
	    0, -1, 0,
	    0, 0, 1;
	sig_z << 1, 0, 0,
	    0, 1, 0,
	    0, 0, -1;

	sig_dz << 0, 1, 0,
	    1, 0, 0,
	    0, 0, 1;
	sig_dy << 0, 0, 1,
	    0, 1, 0,
	    1, 0, 0;
	sig_dx << 1, 0, 0,
	    0, 0, 1,
	    0, 1, 0;

	sig_dpz << 0, -1, 0,
	    -1, 0, 0,
	    0, 0, 1;
	sig_dpy << 0, 0, -1,
	    0, 1, 0,
	    -1, 0, 0;
	sig_dpx << 1, 0, 0,
	    0, 0, -1,
	    0, -1, 0;
	// Todo: think of the division to equivalence classes (this will do for now) 
	C4v = {{e}, {c4x, c4xinv, c4y, c4yinv, c4z, c4zinv}, {c2x, c2y, c2z}, {sig_x, sig_y, sig_z}, {sig_dx, sig_dy, sig_dz, sig_dpx, sig_dpy, sig_dpz}};
    }
   
    // special points
    std::map<std::string, coord_t<3> > special_points_coords;
    special_points_coords["idx_000"] = {0, 0, 0};
	special_points_coords["idx_L"] = {M_PI, M_PI, M_PI};
	special_points_coords["idx_W"] = {2*M_PI, M_PI, 0};
	special_points_coords["idx_X"] = {2*M_PI, 0, 0};
	special_points_coords["idx_X1"] =  {M_PI, 0, 0};

    // special paths
    std::map<std::string, std::vector<coord_t<3> > > special_paths_coords;
    // none added yet...

    Hubbard<3>::Init(basis, C4v, special_points_coords, special_paths_coords, FCCHubbardParams::REFINE_AT_POINTS);
    
    std::cout << "Precomputing E(p) for the model..." << std::endl;
    PrecomputeE();

}


double FCCHubbard::E_of_p_coord(const coord_t<3> p_coord)
{
    	double px = p_coord(0), py = p_coord(1), pz = p_coord(2);

	double E_of_p = -4 * (std::cos(px/2.0)*std::cos(py/2.0) + std::cos(px/2.0)*std::cos(pz/2.0) + std::cos(py/2.0)*std::cos(pz/2.0))+
	    2*FCCHubbardParams::T_PRIME * (std::cos(px) + std::cos(py) + std::cos(pz)) + 4*FCCHubbardParams::T_PRIME_PRIME*(std::cos(px)*std::cos(py) + std::cos(px)*std::cos(pz) + std::cos(py)*std::cos(pz)) - FCCHubbardParams::MU;

	return E_of_p;
}


void FCCHubbard::PrecomputeE()
{
    PrecomputedE().resize(GetFineMomentaCount());

    for (unsigned p_idx = 0; p_idx < GetFineMomentaCount(); p_idx ++){
	coord_t<3> p = FineMomentumGrid().ptr->m_points[p_idx];
	PrecomputedE()[p_idx] = E_of_p_coord(p);
    }
}


double FCCHubbard::E(const int idx_p, const int idx_w)
{
    /*
-4 (cos kx/2 cos ky/2 + cos kx/2 cos kz/2 + cos ky/2 cos kz/2)+
2T'(cos kx + cos ky + cos kz) + 4T'' (cos kx cosky + coskx cos kz + cos ky cos kz)

T''_c = 1/8 + T' /4
     */
    return PrecomputedE()[idx_p];
}



double FCCHubbard::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       )
{
    return -FCCHubbardParams::UINT;
}

double FCCHubbard::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * FCCHubbardParams::UINT;
    else
	return 0.;
}

double FCCHubbard::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * FCCHubbardParams::UINT;
    else
	return 0.;
}


double FCCHubbard::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0)
        return - MomentumGrid().ptr->get_volume() * FCCHubbardParams::UINT;
    else
	return 0.;
}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double FCCHubbard::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * FCCHubbardParams::UINT;
    else
	return 0.;
}

double FCCHubbard::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return -FCCHubbardParams::UINT;
}
