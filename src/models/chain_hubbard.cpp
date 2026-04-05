#include <models/chain_hubbard.h>
#include <primitive_zone.h>
#include <params_physical.h>

void ChainHubbard::Init()
{
    //#ifdef F_NONZERO
    //throw std::invalid_argument("Chain Hubbard model with extended interactions not implemented (turn off the F_NONZERO flag).");
    //#endif


    std::array<coord_t<1>, 1> basis;
    basis[0](0) = 1.0; // lattice spacing

    // define 'point group' of chain
    std::vector<std::vector<matrix_t<1>> > Z2;

    {
	matrix_t<1> e, sig_x; // just identity and reflection
    
	e << 1;
    
        sig_x << -1;
    
	Z2 = {{e}, {sig_x}};
	std::vector<std::string> names = {"e", "sig_x"};
    }

    // special points
    std::map<std::string, coord_t<1> > special_points_coords;
    coord_t<1> temp;
    temp(0) = 0;
    special_points_coords["idx_0"] = temp;
    temp(0) = M_PI;
    special_points_coords["idx_pi"] = temp;
    temp(0) = -M_PI;
    special_points_coords["idx_minus_pi"] = temp;
    
    // special paths
    std::map<std::string, std::vector<coord_t<1> > > special_paths_coords;
    // there are none..

    Hubbard<1>::Init(basis, Z2, special_points_coords, special_paths_coords, ChainHubbardParams::REFINE_AT_POINTS);

    std::cout << "Precomputing E(p) for the model..." << std::endl;
    PrecomputeE();
}


double ChainHubbard::E_of_p_coord(const coord_t<1> p_coord)
{
    const double px = p_coord(0);
    
    const double cos_px = std::cos(px);
    
    return - 2.0 * cos_px - 4.0 *cos_px * ChainHubbardParams::T_PRIME - ChainHubbardParams::MU;
}

void ChainHubbard::PrecomputeE()
{
    PrecomputedE().resize(GetFineMomentaCount());

    for (unsigned p_idx = 0; p_idx < GetFineMomentaCount(); p_idx ++){
	coord_t<1> p = FineMomentumGrid().ptr->m_points[p_idx];
	PrecomputedE()[p_idx] = E_of_p_coord(p);
    }
}


double ChainHubbard::E(const int idx_p, const int idx_w)
{
    return PrecomputedE()[idx_p];
}

double ChainHubbard::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       )
{

    double val = 0;
    val += -ChainHubbardParams::UINT;
    int idx_K = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));
    
    const complex_vector_t<1> e_iK = FineMomentumGrid().exp_ik_idx_map[idx_K];
    
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    
    val += - 2 * ChainHubbardParams::VINT * cos_Kx;
    return val;

}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double ChainHubbard::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * ChainHubbardParams::UINT;
    else
	return 0.;
}

double ChainHubbard::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * ChainHubbardParams::UINT;
    else
	return 0.;
}

double ChainHubbard::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * ChainHubbardParams::UINT;
    else
	return 0.;
}


double ChainHubbard::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){

	const complex_vector_t<1> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
	// real part of e^{ip} is cos p
	const double cos_Kx  = e_iK(0).real();
        
	val += - MomentumGrid().ptr->get_volume() * ChainHubbardParams::UINT;
        val += - MomentumGrid().ptr->get_volume() * 4 * ChainHubbardParams::VINT * cos_Kx;
	return val;

    }else
	return 0.;
}


double ChainHubbard::F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (FORMFACTOR_SHELL_COUNT == 1.5) // projection into antisymmetric d-wave ff would just give zero.
	return 0;
    if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 2){ 
    	return - MomentumGrid().ptr->get_volume() * ChainHubbardParams::VINT;
    }
    return 0.0;
}

double ChainHubbard::F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (FORMFACTOR_SHELL_COUNT == 1.5)
	return 0;
    if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 2){
	return MomentumGrid().ptr->get_volume() * ChainHubbardParams::VINT;
    }
    return 0.0;
}


double ChainHubbard::F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (FORMFACTOR_SHELL_COUNT == 1.5)
	return 0;
    if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 2){
	return MomentumGrid().ptr->get_volume() * ChainHubbardParams::VINT;
    }
    return 0.0;
}


double ChainHubbard::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return -ChainHubbardParams::UINT;
}
