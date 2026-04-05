#include <models/chain_hubbard.h>
#include <primitive_zone.h>
#include <params_physical.h>


double ChainHubbard::U_of_t(const double t)
{
    return t * ChainHubbardParams::UINT * std::cos(2*M_PI*t * 5);
}

double ChainHubbard::V_of_t(const double t)
{
    return t * ChainHubbardParams::VINT * std::sin(2*M_PI*t * 5);
}


double ChainHubbard::U_of_t_dot(const double t)
{
    return ChainHubbardParams::UINT * (std::cos(2*M_PI*t * 5) - (2*M_PI* 5) * t * std::sin(2*M_PI*t * 5));
}

double ChainHubbard::V_of_t_dot(const double t)
{
    return ChainHubbardParams::VINT * (std::sin(2*M_PI*t * 5) + (2*M_PI* 5) * t * std::cos(2*M_PI*t * 5)); 
}


double ChainHubbard::vertex_4pt_bare_of_t(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t 
			       )
{

    double val = 0;
    val += -U_of_t(t);
    int idx_K = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));
    
    const complex_vector_t<1> e_iK = FineMomentumGrid().exp_ik_idx_map[idx_K];
    
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    
    val += - 2 * V_of_t(t) * cos_Kx;
    return val;

}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double ChainHubbard::vertex_local_part_bare_of_t(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * U_of_t(t);
    else
	return 0.;
}

double ChainHubbard::B_sc_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t)
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * U_of_t(t);
    else
	return 0.;
}

double ChainHubbard::B_m_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * U_of_t(t);
    else
	return 0.;
}


double ChainHubbard::B_d_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){

	const complex_vector_t<1> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
	// real part of e^{ip} is cos p
	const double cos_Kx  = e_iK(0).real();
        
	val += - MomentumGrid().ptr->get_volume() * U_of_t(t);
        val += - MomentumGrid().ptr->get_volume() * 4 * V_of_t(t) * cos_Kx;
	return val;

    }else
	return 0.;
}


double ChainHubbard::vertex_4pt_local_part_bare_of_t(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t)
{
    return -U_of_t(t);
}



// the derivatives

double ChainHubbard::vertex_4pt_bare_of_t_dot(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t 
			       )
{

    double val = 0;
    val += -U_of_t_dot(t);
    int idx_K = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));
    
    const complex_vector_t<1> e_iK = FineMomentumGrid().exp_ik_idx_map[idx_K];
    
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    
    val += - 2 * V_of_t_dot(t) * cos_Kx;
    return val;
}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double ChainHubbard::vertex_local_part_bare_of_t_dot(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * U_of_t_dot(t);
    else
	return 0.;
}

double ChainHubbard::B_sc_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t)
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * U_of_t_dot(t);
    else
	return 0.;
}

double ChainHubbard::B_m_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * U_of_t_dot(t);
    else
	return 0.;
}


double ChainHubbard::B_d_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){

	const complex_vector_t<1> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
	// real part of e^{ip} is cos p
	const double cos_Kx  = e_iK(0).real();
        
	val += - MomentumGrid().ptr->get_volume() * U_of_t_dot(t);
        val += - MomentumGrid().ptr->get_volume() * 4 * V_of_t_dot(t) * cos_Kx;
	return val;

    }else
	return 0.;
}


double ChainHubbard::vertex_4pt_local_part_bare_of_t_dot(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t)
{
    return -U_of_t_dot(t);
}

