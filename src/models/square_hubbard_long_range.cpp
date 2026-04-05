#include <models/square_hubbard_long_range.h>
#include <primitive_zone.h>
#include <params_physical.h>
#include <frequencies/matsubara_space.h>




double SquareHubbardLongRange::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out
			       )
{
# ifdef F_NONZERO
    double val = 0;
    
    int idx_K = SumFineMomentaIdxes(idx_p1_out, GetNegativeFineMomentumIdx(idx_p1_in));
    
    const complex_vector_t<2> e_iK = FineMomentumGrid().exp_ik_idx_map[idx_K];
    
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    const double cos_Ky  = e_iK(1).real();

    const double K_s_squared = SquareHubbardLongRangeParams::K_s* SquareHubbardLongRangeParams::K_s;

    const double A = cos_Kx + cos_Ky - 3 - K_s_squared;

    val =  - SquareHubbardLongRangeParams::Q/(2* std::sqrt(A*A - 1));
#endif
    return 0;
}

//void SquareHubbardLongRange::PrecomputeV0_local()
//{
  
//}


double SquareHubbardLongRange::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (SquareHubbardParams::UINT);
    else
	return 0.;
}

double SquareHubbardLongRange::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{

  if(idx_m==0 && idx_mp==0) 
	return -MomentumGrid().ptr->get_volume() * (-SquareHubbardParams::UINT);
    else
	return 0.;
}


double SquareHubbardLongRange::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
  if(!(idx_m==0 && idx_mp==0)){ // todo: check why this even depends on ff indices ...
    return 0.0;
  }

  double V0_bar = 0;

# ifdef F_NONZERO // todo: do something about these flags

    const complex_vector_t<2> e_iK = MomentumGrid().exp_ik_idx_map[idx_K];
    
    // real part of e^{ip} is cos p
    const double cos_Kx  = e_iK(0).real();
    const double cos_Ky  = e_iK(1).real();

    const double K_s_squared = SquareHubbardLongRangeParams::K_s* SquareHubbardLongRangeParams::K_s;

    const double A = cos_Kx + cos_Ky - 3 - K_s_squared;

    V0_bar =  (SquareHubbardLongRangeParams::Q/(2* std::sqrt(A*A - 1)) - SquareHubbardParams::UINT);
#endif

    
    return - MomentumGrid().ptr->get_volume() *(SquareHubbardParams::UINT + 2*V0_bar);
}

// TODO: Postprocessing susceptibilities require these fermionic interactions, and so far they are not correct from this model
double SquareHubbardLongRange::F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    return 0.0;
}

double SquareHubbardLongRange::F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;
    return val;	
}


double SquareHubbardLongRange::F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    //todo: adjust
    double val = 0.0;
    return val;	
}


double SquareHubbardLongRange::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * (SquareHubbardParams::UINT);
    else
	return 0.;
}


double SquareHubbardLongRange::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return - SquareHubbardParams::UINT;
}
