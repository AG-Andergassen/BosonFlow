#include <models/square_hubbard_peierls.h>
#include <primitive_zone.h>
#include <params_physical.h>
#include <frequencies/matsubara_space.h>


void SquareHubbardPeierls::Init()
{
    SquareHubbard::Init();

    if (SquareHubbardPeierlsParams::PHONON_TYPE != "optical")
        PrecomputeAcousticPhononDispersion();
}

void SquareHubbardPeierls::PrecomputeAcousticPhononDispersion()
{
    PrecomputedAcousticPhononDispersion().resize(GetFineMomentaCount());

    for (unsigned p_idx = 0; p_idx < GetFineMomentaCount(); p_idx ++){
	coord_t<2> p_over_two = FineMomentumGrid().ptr->m_points[p_idx]/2.0;
	
        PrecomputedAcousticPhononDispersion()[p_idx] = 0.5*(std::abs(std::sin(p_over_two(0))) + std::abs(std::sin(p_over_two(1))));
	if (PrecomputedAcousticPhononDispersion()[p_idx] == 0.0)
	    PrecomputedAcousticPhononDispersion()[p_idx] = 1e-10; // very small number
    }
}


double SquareHubbardPeierls::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out
			       )
{
# ifdef F_NONZERO
    double val = 0;
    val += -SquareHubbardParams::UINT;
    
    int idx_K_ph = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));
    int idx_W_ph = idx_w1_out - idx_w1_in;

    int idx_K_xph = SumFineMomentaIdxes(idx_p1_in, GetNegativeFineMomentumIdx(idx_p2_out));
    //int idx_W_xph = idx_w1_in - idx_w2_out;

    int idx_K_pp = SumFineMomentaIdxes(idx_p1_out, idx_p2_out);
    //int idx_W_pp = idx_w1_out + idx_w2_out + 1;
    
    const complex_vector_t<2> e_iK_ph = FineMomentumGrid().exp_ik_idx_map[idx_K_ph];
    const complex_vector_t<2> e_iK_xph = FineMomentumGrid().exp_ik_idx_map[idx_K_xph];
    const complex_vector_t<2> e_iK_pp = FineMomentumGrid().exp_ik_idx_map[idx_K_pp];
    
    // real part of e^{ip} is cos p
    const double cos_Kx_ph  = e_iK_ph(0).real();
    const double cos_Ky_ph  = e_iK_ph(1).real();

    const double cos_Kx_xph  = e_iK_xph(0).real();
    const double cos_Ky_xph  = e_iK_xph(1).real();
    
    const double cos_Kx_pp  = e_iK_pp(0).real();
    const double cos_Ky_pp  = e_iK_pp(1).real();
        
    
    const double W_ph = W_val(idx_W_ph);

    double Omega0;
    if (SquareHubbardPeierlsParams::PHONON_TYPE == "optical")
	Omega0 = SquareHubbardPeierlsParams::OMEGA0;
    else{
	// |sin K_x| + |sin K_y|
	Omega0 = SquareHubbardPeierlsParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[idx_K_ph];
    }
    double V_R = -SquareHubbardPeierlsParams::V_ph  + SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(W_ph*W_ph + Omega0*Omega0);
    
    val -= V_R * 2.0 *(cos_Kx_pp + cos_Ky_pp + cos_Kx_xph + cos_Ky_xph );
    
    return val;
#endif
    return -SquareHubbardParams::UINT;
}



double SquareHubbardPeierls::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0){ 
	double val = 0;
	const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
	// real part of e^{ip} is cos p
	const double cos_Kx  = e_iK(0).real();
	const double cos_Ky  = e_iK(1).real();
        
	val += - MomentumGrid().ptr->get_volume() * SquareHubbardParams::UINT;
        
	// electron phonon coupling
	val += + MomentumGrid().ptr->get_volume() * 2.0 * SquareHubbardPeierlsParams::V_ph * (cos_Kx + cos_Ky);
	
	return val;    
    }
    else
	return 0.;
}

double SquareHubbardPeierls::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0){ 
	double val = 0;

	const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
	// real part of e^{ip} is cos p
	const double cos_Kx  = e_iK(0).real();
	const double cos_Ky  = e_iK(1).real();
        
	val += MomentumGrid().ptr->get_volume() * SquareHubbardParams::UINT;
            
	// electron phonon coupling
	val += - MomentumGrid().ptr->get_volume() * 2.0 * SquareHubbardPeierlsParams::V_ph * (cos_Kx + cos_Ky);
	
	return val;	
    }else
	return 0.;
}


double SquareHubbardPeierls::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0){
	double val = 0.0;
	const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
	// real part of e^{ip} is cos p
	const double cos_Kx  = e_iK(0).real();
	const double cos_Ky  = e_iK(1).real();
        
	val += - MomentumGrid().ptr->get_volume() * SquareHubbardParams::UINT;
        val += - MomentumGrid().ptr->get_volume() * 4 * SquareHubbardParams::VINT * (cos_Kx + cos_Ky);
	    
	// electron phonon coupling
	val += - MomentumGrid().ptr->get_volume() * 2.0 * SquareHubbardPeierlsParams::V_ph * (cos_Kx + cos_Ky);
	
	return val;
    }else
	return 0.;
}

// TODO: Postprocessing susceptibilities require these fermionic interactions, and so far they are not correct from this model
double SquareHubbardPeierls::F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;
    
    const int idx_W_ph = idx_wp - idx_w;

    const double W = W_val(idx_W_ph, 0);
    double Omega0;
    if (SquareHubbardPeierlsParams::PHONON_TYPE == "optical")
	Omega0 = SquareHubbardPeierlsParams::OMEGA0;
    else{
	// |sin K_x| + |sin K_y|
	//Omega0 = SquareHubbardPeierlsParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[GetFineMomentumIdxFromCoarse(idx_K)];
    }

    const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];

    const double cos_Kx  = e_iK(0).real();
    const double cos_Ky  = e_iK(1).real();

    const double sin_Kx  = e_iK(0).imag();
    const double sin_Ky  = e_iK(1).imag();


    if (idx_m == 0 && idx_m == idx_mp){
	val += -MomentumGrid().ptr->get_volume() * 2.0*SquareHubbardPeierlsParams::V_ph * W*W/(W*W + Omega0*Omega0) * (cos_Kx + cos_Ky);

    }
    else if ((idx_mp == 3 || idx_mp == 4)  && idx_mp == idx_m ){ // check if d-wave or extended s-wave
      val += MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * (cos_Kx + cos_Ky);
    }
    else if (idx_mp == idx_m){ // diagonal px or py component
	double cos_;
	if (idx_m == 1)
	    cos_ = cos_Kx;
	else
	    cos_ = cos_Ky;
	
	// factor of 2.0 comes from the sqrt(2) that multiplies the p-wave ff
	val += -MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * cos_;
	val *= 2.0;
	
    }else if ((idx_m == 1 && idx_mp == 3) || (idx_m == 3 && idx_mp == 1)){ // px-d component
        val += MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * sin_Kx;
	val *= std::sqrt(2.0);
    }else if ((idx_m == 2 && idx_mp == 3) || (idx_m == 3 && idx_mp == 2)){ // py-d component
        val += -MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * sin_Ky;
	val *= std::sqrt(2.0);
    }else if ((idx_m == 1 && idx_mp == 4) || (idx_m == 4 && idx_mp == 1)){ // px-d component
        val += MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * sin_Kx;
	val *= std::sqrt(2.0);
    }else if ((idx_m == 2 && idx_mp == 4) || (idx_m == 4 && idx_mp == 2)){ // py-d component
        val += MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * sin_Ky;
	val *= std::sqrt(2.0);
    }else if ((idx_m == 3 && idx_mp == 4) || (idx_m == 4 && idx_mp == 3)){ // esw-d component
      val += MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * (cos_Kx - cos_Ky);
    }
    return val;
}

double SquareHubbardPeierls::F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
  // note: \mathcal{F}^{M}(q) = - \mathcal{F}^{SC}(-q)
  
    double val = 0.0;
    
    const int idx_W_ph = idx_wp - idx_w;

    const double W = W_val(idx_W_ph, 0);
    double Omega0;
    if (SquareHubbardPeierlsParams::PHONON_TYPE == "optical")
	Omega0 = SquareHubbardPeierlsParams::OMEGA0;
    else{
	// |sin K_x| + |sin K_y|
	//Omega0 = SquareHubbardPeierlsParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[GetFineMomentumIdxFromCoarse(idx_K)];
    }

    const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];

    const double cos_Kx  = e_iK(0).real();
    const double cos_Ky  = e_iK(1).real();

    const double sin_Kx  = e_iK(0).imag();
    const double sin_Ky  = e_iK(1).imag();


    if (idx_m == 0 && idx_m == idx_mp){
	val += MomentumGrid().ptr->get_volume() * 2.0*SquareHubbardPeierlsParams::V_ph * W*W/(W*W + Omega0*Omega0) * (cos_Kx + cos_Ky);

    }
    else if ((idx_mp == 3 || idx_mp == 4)  && idx_mp == idx_m ){ // check if d-wave or extended s-wave
      val += -MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * (cos_Kx + cos_Ky);
    }
    else if (idx_mp == idx_m){ // diagonal px or py component
	double cos_;
	if (idx_m == 1)
	    cos_ = cos_Kx;
	else
	    cos_ = cos_Ky;
	
	// factor of 2.0 comes from the sqrt(2) that multiplies the p-wave ff
	val += MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * cos_;
	val *= 2.0;
	
    }else if ((idx_m == 1 && idx_mp == 3) || (idx_m == 3 && idx_mp == 1)){ // px-d component
        val += MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * sin_Kx;
	val *= std::sqrt(2.0);
    }else if ((idx_m == 2 && idx_mp == 3) || (idx_m == 3 && idx_mp == 2)){ // py-d component
        val += -MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * sin_Ky;
	val *= std::sqrt(2.0);
    }else if ((idx_m == 1 && idx_mp == 4) || (idx_m == 4 && idx_mp == 1)){ // px-esw component
        val += MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * sin_Kx;
	val *= std::sqrt(2.0);
    }else if ((idx_m == 2 && idx_mp == 4) || (idx_m == 4 && idx_mp == 2)){ // py-esw component
        val += MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * sin_Ky;
	val *= std::sqrt(2.0);
    }
    if ((idx_m == 3 && idx_mp == 4) || (idx_m == 4 && idx_mp == 3)){ // esw-d component
      val += -MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W*W + Omega0*Omega0) * (cos_Kx - cos_Ky);
    }
    return val;
}


double SquareHubbardPeierls::F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
  double val = 0.0;
    
  const int idx_W_xph = idx_wp - idx_w;

  const double W_xph = W_val(idx_W_xph, 0);
  const double W_ph = W_val(idx_W, 0);

  double Omega0;
  if (SquareHubbardPeierlsParams::PHONON_TYPE == "optical")
    Omega0 = SquareHubbardPeierlsParams::OMEGA0;
  else{
    // |sin K_x| + |sin K_y|
    //Omega0 = SquareHubbardPeierlsParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[GetFineMomentumIdxFromCoarse(idx_K)];
  }

  const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];

  const double cos_Kx  = e_iK(0).real();
  const double cos_Ky  = e_iK(1).real();

  const double sin_Kx  = e_iK(0).imag();
  const double sin_Ky  = e_iK(1).imag();


    if (idx_m == 0 && idx_m == idx_mp){
	val += -F_sc( idx_W, idx_K, idx_w, idx_m, idx_wp, idx_mp );
    }else if ((idx_mp == 3 || idx_mp == 4)  && idx_mp == idx_m ){ // check if d-wave or extended s-wave
      val += -MomentumGrid().ptr->get_volume()*(-SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W_ph*W_ph + Omega0*Omega0) * (cos_Kx + cos_Ky + 2) + 0.5 * SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W_xph*W_xph + Omega0*Omega0) * (cos_Kx + cos_Ky));
    }else if (idx_m == idx_mp){ // px + py component
	double cos_;
	if (idx_m == 1)
	    cos_ = cos_Kx;
	else
	    cos_ = cos_Ky;
	
	val += -MomentumGrid().ptr->get_volume()*(-0.5*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W_xph*W_xph + Omega0*Omega0) * cos_ - SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(W_ph*W_ph + Omega0*Omega0) * (1.0 - cos_));
 	val *= 2.0;
    }else if ((idx_m == 1 && idx_mp == 3) || (idx_m == 3 && idx_mp == 1)){ // px-d component
      val += -MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * (Omega0*Omega0/(W_ph*W_ph + Omega0*Omega0) * sin_Kx - 0.5 * Omega0*Omega0/(W_xph*W_xph + Omega0*Omega0) * sin_Kx);
	val *= std::sqrt(2.0);
    }else if ((idx_m == 2 && idx_mp == 3) || (idx_m == 3 && idx_mp == 2)){ // py-d component
        val += -MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * (-Omega0*Omega0/(W_ph*W_ph + Omega0*Omega0) * sin_Ky  + 0.5 * Omega0*Omega0/(W_xph*W_xph + Omega0*Omega0) * sin_Ky);
	val *= std::sqrt(2.0);
    }else if ((idx_m == 1 && idx_mp == 4) || (idx_m == 4 && idx_mp == 1)){ // px-esw component
      val += -MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * (Omega0*Omega0/(W_ph*W_ph + Omega0*Omega0) * sin_Kx - 0.5 * Omega0*Omega0/(W_xph*W_xph + Omega0*Omega0) * sin_Kx);
	val *= std::sqrt(2.0);
    }else if ((idx_m == 2 && idx_mp == 4) || (idx_m == 4 && idx_mp == 2)){ // py-esw component
        val += -MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * (Omega0*Omega0/(W_ph*W_ph + Omega0*Omega0) * sin_Ky  - 0.5 * Omega0*Omega0/(W_xph*W_xph + Omega0*Omega0) * sin_Ky);
	val *= std::sqrt(2.0);
    }else if ((idx_m == 3 && idx_mp == 4) || (idx_m == 4 && idx_mp == 3)){ // d-esw component
      val += -MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * (-Omega0*Omega0/(W_ph*W_ph + Omega0*Omega0)   + 0.5 * Omega0*Omega0/(W_xph*W_xph + Omega0*Omega0) )* (cos_Kx - cos_Ky);
    }
    return val;	
}


double SquareHubbardPeierls::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * SquareHubbardParams::UINT;
    else
	return 0.;
}


double SquareHubbardPeierls::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return -SquareHubbardParams::UINT;
}


double SquareHubbardPeierls::vertex_DC(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out
			       )
{
    double val = 0;
    val -= -2.0*SquareHubbardParams::UINT;
    
    int idx_K_ph = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));
    int idx_W_ph = idx_w1_out - idx_w1_in;

    int idx_K_xph = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p2_in));
    //int idx_W_xph = idx_w1_in - idx_w2_out;

    int idx_K_pp = SumFineMomentaIdxes(idx_p1_in, idx_p2_in);
    //int idx_W_pp = idx_w1_out + idx_w2_out + 1;
    
    const complex_vector_t<2> e_iK_ph = FineMomentumGrid().exp_ik_idx_map[idx_K_ph];
    const complex_vector_t<2> e_iK_xph = FineMomentumGrid().exp_ik_idx_map[idx_K_xph];
    const complex_vector_t<2> e_iK_pp = FineMomentumGrid().exp_ik_idx_map[idx_K_pp];
    
    // real part of e^{ip} is cos p
    const double cos_Kx_ph  = e_iK_ph(0).real();
    const double cos_Ky_ph  = e_iK_ph(1).real();

    const double cos_Kx_xph  = e_iK_xph(0).real();
    const double cos_Ky_xph  = e_iK_xph(1).real();
    
    const double cos_Kx_pp  = e_iK_pp(0).real();
    const double cos_Ky_pp  = e_iK_pp(1).real();
        
    
    const double W_ph = W_val(idx_W_ph);

    double Omega0;
    if (SquareHubbardPeierlsParams::PHONON_TYPE == "optical")
	Omega0 = SquareHubbardPeierlsParams::OMEGA0;
    else{
	// |sin K_x| + |sin K_y|
	Omega0 = SquareHubbardPeierlsParams::OMEGA0 * PrecomputedAcousticPhononDispersion()[idx_K_ph];
    }
    
    val -= 2.0*SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(W_ph*W_ph + Omega0*Omega0) *(cos_Kx_pp + cos_Ky_pp + cos_Kx_xph + cos_Ky_xph );
    
    return val;
}



double SquareHubbardPeierls::vertex_DC_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    const double cos_Ky  = e_iK(1).real();

    const double sin_Kx  = e_iK(0).imag();
    const double sin_Ky  = e_iK(1).imag();
        
    const double W_ph = W_val(idx_wp - idx_w);

    double Omega0;
    if (SquareHubbardPeierlsParams::PHONON_TYPE == "optical")
	Omega0 = SquareHubbardPeierlsParams::OMEGA0;
    else{
	// |sin K_x| + |sin K_y|
	//throw std::exception(std::string("Acoustic phonon for Peierls model not implemented yet"));// needs to handle projection into s-wave form factor
    }


    if (idx_m == 0 && idx_mp == 0){
      return -2.0 * SquareHubbardPeierls::vertex_local_part_bare(idx_W, idx_w, idx_m, idx_wp, idx_mp) - MomentumGrid().ptr->get_volume() * (cos_Kx + cos_Ky)*2* SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(W_ph*W_ph + Omega0*Omega0);
    }else if((idx_m == 3 || idx_m == 4) && idx_mp == idx_m){
      return -0.5*MomentumGrid().ptr->get_volume() * SquareHubbardPeierlsParams::V_ph* W_ph*W_ph/(W_ph*W_ph + Omega0 * Omega0) * (cos_Kx + cos_Ky);
    }else if (idx_m == idx_mp){ // diagonal px or py component
      double cos_;
      if (idx_m == 1)
	cos_ = cos_Kx;
      else
	cos_ = cos_Ky;
      // factor of two to account for p-wave ff normalisation
      return + 2.0*MomentumGrid().ptr->get_volume() * 0.5 * SquareHubbardPeierlsParams::V_ph* W_ph*W_ph/(W_ph*W_ph + Omega0 * Omega0) * cos_;
    }else if ((idx_m == 1 && idx_mp == 3) || (idx_m == 3 && idx_mp == 1)){ // px-d component
      return std::sqrt(2.0) * -MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(Omega0*Omega0 + W_ph*W_ph) * sin_Kx;
    }else if ((idx_m == 2 && idx_mp == 3) || (idx_m == 3 && idx_mp == 2)){ // py-d component
      return std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(Omega0*Omega0 + W_ph*W_ph) * sin_Ky;
    }else if ((idx_m == 1 && idx_mp == 4) || (idx_m == 4 && idx_mp == 1)){ // px-esw component
      return std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * (1 + Omega0*Omega0/(Omega0*Omega0 + W_ph*W_ph)) * sin_Kx;
    }else if ((idx_m == 2 && idx_mp == 4) || (idx_m == 4 && idx_mp == 2)){ // py-esw component
      return std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * (1 + Omega0*Omega0/(Omega0*Omega0 + W_ph*W_ph)) * sin_Ky;
    }else if ((idx_m == 3 && idx_mp == 4) || (idx_m == 4 && idx_mp == 3)){ // d-esw component
      return MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * (1 + Omega0*Omega0/(Omega0*Omega0 + W_ph*W_ph)) * (cos_Kx - cos_Ky);
    }else{
	return 0;
    }
}

double SquareHubbardPeierls::vertex_DC_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    const double cos_Ky  = e_iK(1).real();

    const double sin_Kx  = e_iK(0).imag();
    const double sin_Ky  = e_iK(1).imag();
        
    const double W_xph = W_val(idx_wp - idx_w);
    const double W_ph = W_val(idx_W);
	
	
    double Omega0;
    if (SquareHubbardPeierlsParams::PHONON_TYPE == "optical")
	Omega0 = SquareHubbardPeierlsParams::OMEGA0;
    else{
	// |sin K_x| + |sin K_y|
	//throw std::exception(std::string("Acoustic phonon for Peierls model not implemented yet"));// needs to handle projection into s-wave form factor
    }


    if (idx_m == 0 && idx_mp == 0){
        return -vertex_DC_sc( idx_W, idx_K, idx_w, idx_m, idx_wp, idx_mp );
    }
    else if ((idx_mp == 3 || idx_mp == 4) && idx_mp == idx_m){
      return -MomentumGrid().ptr->get_volume() * SquareHubbardPeierlsParams::V_ph * (-0.5 * W_xph*W_xph/(W_xph*W_xph + Omega0 * Omega0) * (cos_Kx + cos_Ky) + W_ph*W_ph/(W_ph*W_ph + Omega0*Omega0)*(2.0 + cos_Kx + cos_Ky));
    }
    else if (idx_m == idx_mp){
	double cos_;
	if (idx_m == 1)
	    cos_ = cos_Kx;
	else
	    cos_ = cos_Ky;
	// factor of 2 p-wave ff normalisation
	return - 2.0*MomentumGrid().ptr->get_volume() * SquareHubbardPeierlsParams::V_ph * (0.5 * W_xph*W_xph/(W_xph*W_xph + Omega0 * Omega0) * cos_ + W_ph*W_ph/(W_ph*W_ph + Omega0*Omega0)*(1.0 - cos_));
    }else if ((idx_m == 1 && idx_mp == 3) || (idx_m == 3 && idx_mp == 1)){ // px-d component
      return std::sqrt(2.0) * -MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * (-W_ph*W_ph/(Omega0*Omega0 + W_ph*W_ph) * sin_Kx + 0.5 * W_xph*W_xph/(Omega0*Omega0 + W_xph*W_xph) * sin_Kx);
    }else if ((idx_m == 2 && idx_mp == 3) || (idx_m == 3 && idx_mp == 2)){ // py-d component
      return std::sqrt(2.0) * -MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * (W_ph*W_ph/(Omega0*Omega0 + W_ph*W_ph) * sin_Ky - 0.5 * W_xph*W_xph/(Omega0*Omega0 + W_xph*W_xph) * sin_Ky);
    }else if ((idx_m == 1 && idx_mp == 4) || (idx_m == 4 && idx_mp == 1)){ // px-esw component
      return -std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * ( Omega0*Omega0/(Omega0*Omega0 + W_ph*W_ph) - 0.5 * (1.0 + Omega0*Omega0/(Omega0*Omega0 + W_xph*W_xph)))*sin_Kx;
    }else if ((idx_m == 2 && idx_mp == 4) || (idx_m == 4 && idx_mp == 2)){ // py-esw component
      return -std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * ( Omega0*Omega0/(Omega0*Omega0 + W_ph*W_ph) - 0.5 * (1.0 + Omega0*Omega0/(Omega0*Omega0 + W_xph*W_xph)))*sin_Ky;
    }else if ((idx_mp == 4 && idx_m == 3) || (idx_mp == 3 && idx_m == 4)){
      return -MomentumGrid().ptr->get_volume() * SquareHubbardPeierlsParams::V_ph * (  0.5* (1 + Omega0*Omega0/(W_xph*W_xph + Omega0*Omega0)) - Omega0*Omega0/(W_ph*W_ph + Omega0*Omega0) )*(cos_Kx - cos_Ky)  ; 
    }else{
	return 0;
    }
}

/*
else if ((idx_m == 1 && idx_mp == 3) || (idx_m == 3 && idx_mp == 1)){ // px-d component
      return std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * 0.5 * (1 + Omega0*Omega0/(Omega0*Omega0 + W_xph*W_xph) )*sin_Kx;
    }else if ((idx_m == 2 && idx_mp == 3) || (idx_m == 3 && idx_mp == 2)){ // py-d component
      return -std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*SquareHubbardPeierlsParams::V_ph * 0.5 * (1 + Omega0*Omega0/(Omega0*Omega0 + W_xph*W_xph) )*sin_Ky;
    }
*/
double SquareHubbardPeierls::vertex_DC_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
	
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    const double cos_Ky  = e_iK(1).real();

    const double sin_Kx  = e_iK(0).imag();
    const double sin_Ky  = e_iK(1).imag();
        
    const double W_ph = W_val(idx_wp - idx_w);
	
    double Omega0;
    if (SquareHubbardPeierlsParams::PHONON_TYPE == "optical")
	Omega0 = SquareHubbardPeierlsParams::OMEGA0;
    else{
	// |sin K_x| + |sin K_y|
	//throw std::exception (std::string("Acoustic phonon for Peierls model not implemented yet") );// needs to handle projection into s-wave form factor
    }

    if (idx_m == 0 && idx_mp == 0){
      return -vertex_DC_sc( idx_W, idx_K, idx_w, idx_m, idx_wp, idx_mp );
    }else if((idx_m == 3 || idx_m == 4) && idx_mp == idx_m){
      return MomentumGrid().ptr->get_volume() * 0.5 * SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(W_ph*W_ph + Omega0 * Omega0) * (cos_Kx + cos_Ky);
    }else if (idx_m == idx_mp){ // diagonal px or py component
	double cos_;
	if (idx_m == 1)
	    cos_ = cos_Kx;
	else
	    cos_ = cos_Ky;
	// factor of 2 p-wave ff normalisation
	return - 2.0*MomentumGrid().ptr->get_volume() * 0.5 * SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(W_ph*W_ph + Omega0 * Omega0) * cos_;
    }else if ((idx_m == 1 && idx_mp == 3) || (idx_m == 3 && idx_mp == 1)){ // px-d component
        return std::sqrt(2.0) * -MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(Omega0*Omega0 + W_ph*W_ph) * sin_Kx;
    }else if ((idx_m == 2 && idx_mp == 3) || (idx_m == 3 && idx_mp == 2)){ // py-d component
        return std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(Omega0*Omega0 + W_ph*W_ph) * sin_Ky;
    }else if ((idx_m == 1 && idx_mp == 4) || (idx_m == 4 && idx_mp == 1)){ // px-esw component
      return -std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * (W_ph*W_ph/(Omega0*Omega0 + W_ph*W_ph)) * sin_Kx;
    }else if ((idx_m == 2 && idx_mp == 4) || (idx_m == 4 && idx_mp == 2)){ // py-esw component
      return -std::sqrt(2.0) * MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * (W_ph*W_ph/(Omega0*Omega0 + W_ph*W_ph)) * sin_Ky;
    }else if ((idx_m == 3 && idx_mp == 4) || (idx_m == 4 && idx_mp == 3)){ // d-esw component
      return MomentumGrid().ptr->get_volume()*0.5*SquareHubbardPeierlsParams::V_ph * (W_ph*W_ph/(Omega0*Omega0 + W_ph*W_ph)) * (cos_Kx - cos_Ky);
    }else{
	return 0;
    }
}


void SquareHubbardPeierls::Test()
{
  std::vector<int> W_idxes = {0, 0, -6, 100000, 2};
  std::vector<int> w_idxes = {0, 5, 1, 10, 6};
  std::vector<int> wp_idxes = {0, -3, 0, 10, 2};
  std::vector<int> q_idxes = {0, 13, 12, 5, 15};
  
  auto V_sc = [](double qx, double qy, double kx, double ky, double kpx, double kpy, int idx_W, int idx_w, int idx_wp){
      const int idx_W_xph = idx_wp - idx_w;

      const double W_xph = W_val(idx_W_xph, 0);
      const double W_ph = W_val(idx_W, 0);

      const double Omega0 = SquareHubbardPeierlsParams::OMEGA0;

      return MomentumGrid().ptr->get_volume()* 2*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(Omega0*Omega0 + W_xph*W_xph)*
	(std::cos(qx) + std::cos(qy) + std::cos(kx + kpx - qx) + std::cos(ky + kpy - qy));
	;

    };

  auto V_d = [](double qx, double qy, double kx, double ky, double kpx, double kpy, int idx_W, int idx_w, int idx_wp){
      const int idx_W_xph = idx_wp - idx_w;

      const double W_xph = W_val(idx_W_xph, 0);
      const double W_ph = W_val(idx_W, 0);

      const double Omega0 = SquareHubbardPeierlsParams::OMEGA0;

      return MomentumGrid().ptr->get_volume()* (
						4*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(Omega0*Omega0 + W_ph*W_ph)*
						(std::cos(kx - kpx) + std::cos(ky - kpy) + std::cos(kx + kpx + qx) + std::cos(ky + kpy + qy)) -2*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(Omega0*Omega0 + W_xph*W_xph) * (std::cos(qx) + std::cos(qy) + std::cos(kx + kpx + qx) +  + std::cos(ky + kpy + qy))
						) ;
	;

    };

  auto V_m = [](double qx, double qy, double kx, double ky, double kpx, double kpy, int idx_W, int idx_w, int idx_wp){
      const int idx_W_xph = idx_wp - idx_w;

      const double W_xph = W_val(idx_W_xph, 0);
      const double W_ph = W_val(idx_W, 0);

      const double Omega0 = SquareHubbardPeierlsParams::OMEGA0;

      return -MomentumGrid().ptr->get_volume()* 2*SquareHubbardPeierlsParams::V_ph * Omega0*Omega0/(Omega0*Omega0 + W_xph*W_xph)*
	(std::cos(qx) + std::cos(qy) + std::cos(kx + kpx + qx) + std::cos(ky + kpy + qy));
	;

    };


    auto V_m_dc = [](double qx, double qy, double kx, double ky, double kpx, double kpy, int idx_W, int idx_w, int idx_wp){
      const int idx_W_xph = idx_wp - idx_w;

      const double W_xph = W_val(idx_W_xph, 0);
      const double W_ph = W_val(idx_W, 0);

      const double Omega0 = SquareHubbardPeierlsParams::OMEGA0;

      return MomentumGrid().ptr->get_volume()* 2*SquareHubbardPeierlsParams::V_ph * W_xph*W_xph/(Omega0*Omega0 + W_xph*W_xph)*
	(std::cos(qx) + std::cos(qy) + std::cos(kx + kpx + qx) + std::cos(ky + kpy + qy));
	;

    };


    auto V_d_dc = [](double qx, double qy, double kx, double ky, double kpx, double kpy, int idx_W, int idx_w, int idx_wp){
      const int idx_W_xph = idx_wp - idx_w;

      const double W_xph = W_val(idx_W_xph, 0);
      const double W_ph = W_val(idx_W, 0);

      const double Omega0 = SquareHubbardPeierlsParams::OMEGA0;

      return MomentumGrid().ptr->get_volume()* (2*SquareHubbardPeierlsParams::V_ph * W_xph*W_xph/(Omega0*Omega0 + W_xph*W_xph)*
						(std::cos(qx) + std::cos(qy) + std::cos(kx + kpx + qx) + std::cos(ky + kpy + qy))
						- 4 * SquareHubbardPeierlsParams::V_ph * W_ph*W_ph/(W_ph*W_ph + Omega0*Omega0) * (std::cos(kx + kpx + qx) + std::cos(ky + kpy + qy) + std::cos(kx - kpx) + std::cos(ky - kpy) ) 
						);
	;

    };

    auto V_sc_dc = [](double qx, double qy, double kx, double ky, double kpx, double kpy, int idx_W, int idx_w, int idx_wp){
      const int idx_W_xph = idx_wp - idx_w;

      const double W_xph = W_val(idx_W_xph, 0);
      const double W_ph = W_val(idx_W, 0);

      const double Omega0 = SquareHubbardPeierlsParams::OMEGA0;

      return -MomentumGrid().ptr->get_volume()* 2*SquareHubbardPeierlsParams::V_ph * W_xph*W_xph/(Omega0*Omega0 + W_xph*W_xph)*
	(std::cos(qx) + std::cos(qy) + std::cos(kx + kpx - qx) + std::cos(ky + kpy - qy));
	;

    };
 
  
    auto f_swave = [](double kx, double ky){
      return 1.0;
    };

    auto f_dwave = [](double kx, double ky){
      return std::cos(kx) - std::cos(ky);
    };

    auto f_pxwave = [](double kx, double ky){
      return std::sin(kx) * std::sqrt(2);
    };

    auto f_pywave = [](double kx, double ky){
      return std::sin(ky) * std::sqrt(2);
    };

    auto f_eswave = [](double kx, double ky){
      return std::cos(kx) + std::cos(ky);
    };

    std::vector<std::function<double(double, double)>> formfactors = {
								      f_swave, f_pxwave, f_pywave, f_dwave, f_eswave
    };
    
    int f1 = 0;
    int f2 = 0;
    
    
    {
      std::vector<double> vals = {0.0, 0.0, 0.0, 0.0, 0.0};
      for (auto k : MomentumGrid().ptr->m_points)
	for (auto kp : MomentumGrid().ptr->m_points){
	  for (int i = 0; i < vals.size(); i++ )
	    vals[i] += V_m_dc(MomentumGrid().ptr->m_points[q_idxes[i]](0),
			    MomentumGrid().ptr->m_points[q_idxes[i]](1),
			    k(0), k(1), kp(0), kp(1), W_idxes[i], w_idxes[i], wp_idxes[i]) *
	      formfactors[f1](k(0), k(1)) *
	      formfactors[f2](kp(0), kp(1)) /
	      MomentumGrid().ptr->m_points.size() / MomentumGrid().ptr->m_points.size();
	}

    for (int i = 0; i < vals.size(); i++)
      std::cout << "...Vertex calculated: " <<//F_d
        vertex_DC_m( W_idxes[i], q_idxes[i], w_idxes[i], f1, wp_idxes[i], f2 ) << ", Analytical: " << vals[i] << std::endl;
    }

    /*auto ff_px = GetFormFactorInMomentumIdxSpace(1);
    auto ff_py = GetFormFactorInMomentumIdxSpace(2);
    auto ff_dx2_y2 = GetFormFactorInMomentumIdxSpace(3);
    

    
    std::cout << "px" << std::endl;
    for (auto ff_val : ff_px){
      std::cout << "," << ff_val;
    }
    std::cout << std::endl;

    

    std::cout << "py" << std::endl;
    for (auto ff_val : ff_py){
      std::cout << "," << ff_val;
    }
    std::cout << std::endl;

    std::cout << "d-wave" << std::endl;
    for (auto ff_val : ff_dx2_y2){
      std::cout << "," << ff_val;
    }
    */
    std::cout << std::endl;
    
}
