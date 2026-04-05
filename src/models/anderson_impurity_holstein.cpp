#include <models/anderson_impurity_holstein.h>
#include <primitive_zone.h>
#include <params_physical.h>
#include <frequencies/matsubara_space.h>


void AndersonImpurityHolstein::Init()
{
    AndersonImpurity::Init();
}


double AndersonImpurityHolstein::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out
			       )
{
# ifdef F_NONZERO
    double val = 0;
    val += -(AndersonImpurityParams::UINT - AndersonImpurityHolsteinParams::V_ph);
    
    const int idx_W = idx_w1_out - idx_w1_in;
    
    //if (std::abs(idx_W) < FrequenciesCount::LARGE){
    
    // TODO: with the temperature flow, the second argument must be the flow scale
    //double T = exp( t * LN_10 ) ;  
    //double Sqrt_T = exp( t * LN_10 * 0.5 ) ;  

    const double W = W_val(idx_W);// * BETA * T;

    double Omega0 = AndersonImpurityHolsteinParams::OMEGA0;
        
    // electron phonon coupling
    val -= AndersonImpurityHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);
    //}

    return val;
#endif
    return -AndersonImpurityParams::UINT;
}



double AndersonImpurityHolstein::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (AndersonImpurityParams::UINT - AndersonImpurityHolsteinParams::V_ph);
    else
	return 0.;
}

double AndersonImpurityHolstein::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * (AndersonImpurityParams::UINT - AndersonImpurityHolsteinParams::V_ph);
    else
	return 0.;
}


double AndersonImpurityHolstein::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){
        
	val += - MomentumGrid().ptr->get_volume() * (AndersonImpurityParams::UINT - AndersonImpurityHolsteinParams::V_ph);

	//if (std::abs(idx_W) < FrequenciesCount::LARGE){
	    // TODO: with the temperature flow, the second argument must be the flow scale
	const double W = W_val(idx_W, 0);
	double Omega0 = AndersonImpurityHolsteinParams::OMEGA0;
	
	// electron phonon coupling
	// TODO: think seriously about whether a small imaginary part is needed (so far can't think of anything)
	val += - MomentumGrid().ptr->get_volume() * 2.0 * AndersonImpurityHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);
	//#pragma omp single
	//std::cout << "val and W: " << 2.0 * g0 * g0 * Omega0/(-W*W - Omega0*Omega0) << ", " << W << Omega0 << g0 << std::endl;
	//}
	
	return val;
    }else
	return 0.;
}

// TODO: Postprocessing susceptibilities require these fermionic interactions, and so far they are not correct from this model
double AndersonImpurityHolstein::F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;
	    
        const double W = W_val(idx_W_ph, 0);
	double Omega0 = AndersonImpurityHolsteinParams::OMEGA0;
	
	val += -MomentumGrid().ptr->get_volume() * AndersonImpurityHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);
    }
    return val;
}

double AndersonImpurityHolstein::F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;
	
	//if (std::abs(idx_W_ph) < FrequenciesCount::LARGE){
	// TODO: with the temperature flow, the second argument must be the flow scale
	const double W = W_val(idx_W_ph, 0);
	double Omega0 = AndersonImpurityHolsteinParams::OMEGA0;
	
	val += MomentumGrid().ptr->get_volume() * AndersonImpurityHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);
	    //}
    }

    return val;	
}


double AndersonImpurityHolstein::F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;

	//if (std::abs(idx_W_ph) < FrequenciesCount::LARGE){
	// TODO: with the temperature flow, the second argument must be the flow scale
	const double W = W_val(idx_W_ph, 0);
	double Omega0 = AndersonImpurityHolsteinParams::OMEGA0;
	
	val += MomentumGrid().ptr->get_volume() * AndersonImpurityHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);
    }
    return val;	
}


double AndersonImpurityHolstein::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (AndersonImpurityParams::UINT - AndersonImpurityHolsteinParams::V_ph);
    else
	return 0.;
}


double AndersonImpurityHolstein::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return -(AndersonImpurityParams::UINT - AndersonImpurityHolsteinParams::V_ph);
}
