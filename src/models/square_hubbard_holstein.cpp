#include <models/square_hubbard_holstein.h>
#include <primitive_zone.h>
#include <params_physical.h>
#include <frequencies/matsubara_space.h>


void SquareHubbardHolstein::Init()
{
    SquareHubbard::Init();

    if (SquareHubbardHolsteinParams::PHONON_TYPE != "optical")
        PrecomputeAcousticPhononDispersion();
}

void SquareHubbardHolstein::PrecomputeAcousticPhononDispersion()
{
    PrecomputedAcousticPhononDispersion().resize(GetFineMomentaCount());

    for (unsigned p_idx = 0; p_idx < GetFineMomentaCount(); p_idx ++){
	coord_t<2> p_over_two = FineMomentumGrid().ptr->m_points[p_idx]/2.0;
	
        PrecomputedAcousticPhononDispersion()[p_idx] = 0.5*(std::abs(std::sin(p_over_two(0))) + std::abs(std::sin(p_over_two(1))));
	if (PrecomputedAcousticPhononDispersion()[p_idx] == 0.0)
	    PrecomputedAcousticPhononDispersion()[p_idx] = 1e-10; // very small number
    }
}


double SquareHubbardHolstein::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out
			       )
{
# ifdef F_NONZERO
    double val = 0;
    val += -(SquareHubbardParams::UINT - SquareHubbardHolsteinParams::V_ph);
    
    int idx_K = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));
    
    const complex_vector_t<2> e_iK = FineMomentumGrid().exp_ik_idx_map[idx_K];
    
    
    const int idx_W = idx_w1_out - idx_w1_in;
    

    const double W = W_val(idx_W);// * BETA * T;

    double Omega0;
    if (SquareHubbardHolsteinParams::PHONON_TYPE == "optical")
      Omega0 = SquareHubbardHolsteinParams::OMEGA0;
    else{
      // |sin K_x| + |sin K_y|
      Omega0 = SquareHubbardHolsteinParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[idx_K];
    }
    // electron phonon coupling
    val -= SquareHubbardHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);


    // extended part of the interaction:
    
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    const double cos_Ky  = e_iK(1).real();
    
    val += - 2 * SquareHubbardParams::VINT * ( cos_Kx + cos_Ky);
    
    return val;
#endif
    return -SquareHubbardParams::UINT;
}



double SquareHubbardHolstein::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (SquareHubbardParams::UINT - SquareHubbardHolsteinParams::V_ph);
    else
	return 0.;
}

double SquareHubbardHolstein::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * (SquareHubbardParams::UINT - SquareHubbardHolsteinParams::V_ph);
    else
	return 0.;
}


double SquareHubbardHolstein::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){
	const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
	
	val += - MomentumGrid().ptr->get_volume() * (SquareHubbardParams::UINT - SquareHubbardHolsteinParams::V_ph);

	// real part of e^{ip} is cos p
	const double cos_Kx  = e_iK(0).real();
	const double cos_Ky  = e_iK(1).real();
	val += - MomentumGrid().ptr->get_volume() * 4 * SquareHubbardParams::VINT * (cos_Kx + cos_Ky);

	const double W = W_val(idx_W, 0);
	double Omega0;
	if (SquareHubbardHolsteinParams::PHONON_TYPE == "optical")
	  Omega0 = SquareHubbardHolsteinParams::OMEGA0;
	else{
	  // |sin K_x| + |sin K_y|
	  Omega0 = SquareHubbardHolsteinParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[GetFineMomentumIdxFromCoarse(idx_K)];
	}
	    
	// electron phonon coupling
	val += - MomentumGrid().ptr->get_volume() * 2.0 * SquareHubbardHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);
        
    }
	
    return val;
}

// TODO: Postprocessing susceptibilities require these fermionic interactions, and so far they are not correct from this model
double SquareHubbardHolstein::F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;
	    
	//if (std::abs(idx_W) < FrequenciesCount::LARGE){
	    // TODO: with the temperature flow, the second argument must be the flow scale
	    const double W = W_val(idx_W_ph, 0);
	    double Omega0;
	    if (SquareHubbardHolsteinParams::PHONON_TYPE == "optical")
	      Omega0 = SquareHubbardHolsteinParams::OMEGA0;
	    else{
	      // |sin K_x| + |sin K_y|
	      Omega0 = SquareHubbardHolsteinParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[GetFineMomentumIdxFromCoarse(idx_K)];
	    }
	    val += -MomentumGrid().ptr->get_volume() * SquareHubbardHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);
	    //}
    }else if (FORMFACTOR_SHELL_COUNT == 1.5)
	val += 0;
    else if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 4){
	val += -MomentumGrid().ptr->get_volume() * SquareHubbardParams::VINT;
    }
    return val;
}

double SquareHubbardHolstein::F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;
	
	//if (std::abs(idx_W_ph) < FrequenciesCount::LARGE){
	    // TODO: with the temperature flow, the second argument must be the flow scale
	    const double W = W_val(idx_W_ph, 0);
	    double Omega0;
	    if (SquareHubbardHolsteinParams::PHONON_TYPE == "optical")
	      Omega0 = SquareHubbardHolsteinParams::OMEGA0;
	    else{
	      // |sin K_x| + |sin K_y|
	      Omega0 = SquareHubbardHolsteinParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[GetFineMomentumIdxFromCoarse(idx_K)];
	    }
	    val += MomentumGrid().ptr->get_volume() * SquareHubbardHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);
	    //}
    }else if (FORMFACTOR_SHELL_COUNT == 1.5)
	val += 0;
    else if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 4){
	val += MomentumGrid().ptr->get_volume() * SquareHubbardParams::VINT;
    }
    return val;	
}


double SquareHubbardHolstein::F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;

    if (idx_m == 0 && idx_m == idx_mp){
	const int idx_W_ph = idx_wp - idx_w;

	//if (std::abs(idx_W_ph) < FrequenciesCount::LARGE){
	    // TODO: with the temperature flow, the second argument must be the flow scale
	    const double W = W_val(idx_W_ph, 0);
	    double Omega0;
	    if (SquareHubbardHolsteinParams::PHONON_TYPE == "optical")
	      Omega0 = SquareHubbardHolsteinParams::OMEGA0;
	    else{
	      // |sin K_x| + |sin K_y|
	      //int idx_K = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));
	      Omega0 = SquareHubbardHolsteinParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[GetFineMomentumIdxFromCoarse(idx_K)];
	    }
	    val += MomentumGrid().ptr->get_volume() * SquareHubbardHolsteinParams::V_ph * W*W/(W*W + Omega0*Omega0);
	    //}
    }else if (FORMFACTOR_SHELL_COUNT == 1.5)
	val += 0;
    else if(idx_m == idx_mp && 1 <= idx_m && idx_m <= 4){
	val += MomentumGrid().ptr->get_volume() * SquareHubbardParams::VINT;
    }
    return val;	
}


double SquareHubbardHolstein::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (SquareHubbardParams::UINT - SquareHubbardHolsteinParams::V_ph);
    else
	return 0.;
}


double SquareHubbardHolstein::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return -(SquareHubbardParams::UINT - SquareHubbardHolsteinParams::V_ph);
}


double SquareHubbardHolstein::vertex_DC(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out
			       )
{
    double val = 0;
    val -= -2.0*(SquareHubbardParams::UINT - SquareHubbardHolsteinParams::V_ph);
    return val;
}


double SquareHubbardHolstein::vertex_DC_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (idx_m == 0 && idx_mp == 0){
        return -2.0 * SquareHubbardHolstein::vertex_local_part_bare(idx_W, idx_w, idx_m, idx_wp, idx_mp);
    
    }else{
	return 0;
    }
}

double SquareHubbardHolstein::vertex_DC_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (idx_m == 0 && idx_mp == 0){
        return -2.0 * SquareHubbardHolstein::vertex_local_part_bare(idx_W, idx_w, idx_m, idx_wp, idx_mp);
    }else{
	return 0;
    }
}

double SquareHubbardHolstein::vertex_DC_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if (idx_m == 0 && idx_mp == 0){
        return 2.0 * SquareHubbardHolstein::vertex_local_part_bare(idx_W, idx_w, idx_m, idx_wp, idx_mp);
    }else{
	return 0;
    }
}

