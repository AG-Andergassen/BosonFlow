#include <models/anderson_impurity_holstein.h>
#include <primitive_zone.h>
#include <params_physical.h>
#include <frequencies/matsubara_space.h>


double AndersonImpurityHolstein::U_of_t(const double t)
{
#ifdef TEMP_FLOW
    return AndersonImpurityParams::UINT;
#else
    return t*AndersonImpurityParams::UINT;
#endif
}

double AndersonImpurityHolstein::V_ph_of_t(const double t)
{
#ifdef TEMP_FLOW
    return AndersonImpurityHolsteinParams::V_ph;
#else
    return t*t*AndersonImpurityHolsteinParams::V_ph;
#endif
}

double AndersonImpurityHolstein::U_of_t_dot(const double t)
{
#ifdef TEMP_FLOW
    return 0.0;
#else
    return AndersonImpurityParams::UINT;
#endif
}

double AndersonImpurityHolstein::V_ph_of_t_dot(const double t)
{
#ifdef TEMP_FLOW
    return 0.0;
#else
    return 2.0*t*AndersonImpurityHolsteinParams::V_ph;
#endif 
}

double AndersonImpurityHolstein::L_of_t(const int idx_W, const double t)
{
    double Omega0 = AndersonImpurityHolsteinParams::OMEGA0;

#ifdef TEMP_FLOW
    const double T = std::exp( t * LN_10 );
    const double W = 2*M_PI*idx_W*T;// * BETA * T;
#else
    const double W = W_val(idx_W);
#endif

    return W*W/(W*W + Omega0*Omega0);
}


double AndersonImpurityHolstein::L_of_t_dot(const int idx_W, const double t)
{
    double Omega0 = AndersonImpurityHolsteinParams::OMEGA0;

#ifdef TEMP_FLOW
    const double T = std::exp( t * LN_10 );
    const double dTdt = LN_10 * T;

    const double W = 2*M_PI*idx_W*T;

    return 2.0* W*Omega0*Omega0 * 2*M_PI* dTdt* idx_W /((W*W + Omega0*Omega0)*(W*W + Omega0*Omega0));

#else
    return 0.0;


#endif
    
}


double AndersonImpurityHolstein::vertex_4pt_bare_of_t(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t 
			       )
{
# ifdef F_NONZERO
    double val = 0;
    val += -(U_of_t(t) - V_ph_of_t(t));
    
    const int idx_W = idx_w1_out - idx_w1_in;
    
    val -= V_ph_of_t(t) * L_of_t(idx_W, t);
    
    return val;
#endif
    return -U_of_t(t);
}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double AndersonImpurityHolstein::vertex_local_part_bare_of_t(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (U_of_t(t) - V_ph_of_t(t));
    else
	return 0.;
}

double AndersonImpurityHolstein::B_sc_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t)
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (U_of_t(t) - V_ph_of_t(t));
    else
	return 0.;
}

double AndersonImpurityHolstein::B_m_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * (U_of_t(t) - V_ph_of_t(t));
    else
	return 0.;
}


double AndersonImpurityHolstein::B_d_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){
        
	val += - MomentumGrid().ptr->get_volume() * (U_of_t(t) - V_ph_of_t(t));

        val += - MomentumGrid().ptr->get_volume() * 2.0 * V_ph_of_t(t) * L_of_t(idx_W, t);
        
	return val;
    }else
	return 0.;
}


double AndersonImpurityHolstein::F_sc_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t  )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;
	val += -MomentumGrid().ptr->get_volume() * V_ph_of_t(t) * L_of_t(idx_W_ph, t);
    }
    return val;
}

double AndersonImpurityHolstein::F_m_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t  )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;
        	
	val += MomentumGrid().ptr->get_volume() * V_ph_of_t(t) * L_of_t(idx_W_ph, t);
    }

    return val;	
}


double AndersonImpurityHolstein::F_d_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t  )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;

	const double W = W_val(idx_W_ph, 0);
	double Omega0 = AndersonImpurityHolsteinParams::OMEGA0;
	
	val += MomentumGrid().ptr->get_volume() * V_ph_of_t(t) * L_of_t(idx_W_ph, t);
    }
    return val;	
}


double AndersonImpurityHolstein::vertex_4pt_local_part_bare_of_t(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t)
{
    return -(U_of_t(t) - V_ph_of_t(t));
}



// the derivatives

double AndersonImpurityHolstein::vertex_4pt_bare_of_t_dot(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t 
			       )
{
# ifdef F_NONZERO
    double val = 0;
    val += -(U_of_t_dot(t) - V_ph_of_t_dot(t));
    
    const int idx_W = idx_w1_out - idx_w1_in;
    
        
    val -= V_ph_of_t_dot(t) * L_of_t(idx_W, t) + V_ph_of_t(t) * L_of_t_dot(idx_W, t);

    return val;
#endif
    return -U_of_t_dot(t);

}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double AndersonImpurityHolstein::vertex_local_part_bare_of_t_dot(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (U_of_t_dot(t) - V_ph_of_t_dot(t));
    else
	return 0.;

}

double AndersonImpurityHolstein::B_sc_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t)
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (U_of_t_dot(t) - V_ph_of_t_dot(t));
    else
	return 0.;
}

double AndersonImpurityHolstein::B_m_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * (U_of_t_dot(t) - V_ph_of_t_dot(t));
    else
	return 0.;

}


double AndersonImpurityHolstein::B_d_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){
        
	val += - MomentumGrid().ptr->get_volume() * (U_of_t_dot(t) - V_ph_of_t_dot(t));

	val += - MomentumGrid().ptr->get_volume() * 2.0 * (V_ph_of_t_dot(t) * L_of_t(idx_W, t) + V_ph_of_t(t) * L_of_t_dot(idx_W, t));
	//#pragma omp single
	//std::cout << "val and W: " << 2.0 * g0 * g0 * Omega0/(-W*W - Omega0*Omega0) << ", " << W << Omega0 << g0 << std::endl;
	//}
	
	return val;
    }else
	return 0.;
}

double AndersonImpurityHolstein::F_sc_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t  )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;
	val += -MomentumGrid().ptr->get_volume() * (V_ph_of_t_dot(t) * L_of_t(idx_W_ph, t) + V_ph_of_t(t) * L_of_t_dot(idx_W_ph, t));
    }
    return val;
}

double AndersonImpurityHolstein::F_m_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t  )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;
        
	val += MomentumGrid().ptr->get_volume() * (V_ph_of_t(t) * L_of_t_dot(idx_W_ph, t) + V_ph_of_t_dot(t) * L_of_t(idx_W_ph, t));
    }

    return val;	
}


double AndersonImpurityHolstein::F_d_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t  )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;
	
	val += MomentumGrid().ptr->get_volume() * (V_ph_of_t(t) * L_of_t_dot(idx_W_ph, t) + V_ph_of_t_dot(t) * L_of_t(idx_W_ph, t));
    }
    return val;	
}




double AndersonImpurityHolstein::vertex_4pt_local_part_bare_of_t_dot(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t)
{
    return -(U_of_t_dot(t) - V_ph_of_t_dot(t));
}

