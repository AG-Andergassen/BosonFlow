#include <models/anderson_impurity.h>
#include <primitive_zone.h>
#include <params_physical.h>
#include <frequencies/matsubara_space.h>
#include <stdexcept>

void AndersonImpurity::Init()
{
    // zero dimensional basis is just "1"
    basis_t<0> basis;
    basis[0](0) = 1;

    Hubbard<0>::Init(basis, {}, {}, {}, AndersonImpurityParams::REFINE_AT_POINTS);
}


std::complex<double> AndersonImpurity::E(const int idx_p, const int idx_w)
{
    double w = w_val(idx_w);

    #ifdef TEMP_FLOW
    // todo: let these depend on scale 
    #endif
    
    std::complex<double> Delta(0.0, 0.0);

    if (AndersonImpurityParams::DOS_TYPE == "DISCRETE"){
	for (int n = 0; n < AndersonImpurityParams::BATH_LEVEL_COUNT; n ++){
	    Delta += (std::conj(AndersonImpurityParams::V_vals[n]) * AndersonImpurityParams::V_vals[n])/(I*w  + AndersonImpurityParams::E_vals[n]);
	}
    }else if (AndersonImpurityParams::DOS_TYPE == "BOX"){
	Delta = (I *2.0*AndersonImpurityParams::DELTA_0/M_PI) * std::atan(AndersonImpurityParams::D/w);
    }else if (AndersonImpurityParams::DOS_TYPE == "CONST"){ // broadband limit (D -> \infty)
      double sign = w > 0 ? +1 : -1;
        Delta = I * AndersonImpurityParams::DELTA_0 * sign;
    } else {
        throw std::invalid_argument("Unsupported DOS_TYPE='" + AndersonImpurityParams::DOS_TYPE + "'. Allowed values are BOX or CONST.");
    }
    return - AndersonImpurityParams::MU - Delta;
}


double AndersonImpurity::E_of_p_coord(const coord_t<0> p_coord)
{
    // this function is used only in the computation of the static bubble, so...
    throw std::invalid_argument("The Anderson impurity model is not supported with a static calculation. ");
    return 0;
}


double AndersonImpurity::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       )
{
    return -AndersonImpurityParams::UINT;
}



double AndersonImpurity::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * AndersonImpurityParams::UINT;
    else
	return 0.;
}

double AndersonImpurity::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * AndersonImpurityParams::UINT;
    else
	return 0.;
}


double AndersonImpurity::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){        
	val += - MomentumGrid().ptr->get_volume() * AndersonImpurityParams::UINT;
	return val;
    }else
	return 0.;
}


double AndersonImpurity::F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    return 0.0;
}

double AndersonImpurity::F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    return 0.0;
}


double AndersonImpurity::F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    return 0.0;
}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double AndersonImpurity::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    return - MomentumGrid().ptr->get_volume() * AndersonImpurityParams::UINT;
}

double AndersonImpurity::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return -AndersonImpurityParams::UINT;
}
