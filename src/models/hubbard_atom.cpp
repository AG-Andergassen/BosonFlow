#include <models/hubbard_atom.h>
#include <primitive_zone.h>
#include <params_physical.h>

//#define B_PLUS_F

void HubbardAtom::Init()
{
    // in 0-dimensions, 1 is the only basis element
    basis_t<0> basis;
    basis[0](0) = 1;

    Hubbard<0>::Init(basis, {}, {}, {}, HubbardAtomParams::REFINE_AT_POINTS);

    std::cout << "Volume of BZ: " << MomentumGrid().ptr->get_volume() << std::endl;
    std::cout << "Number of points: " << MomentumGrid().ptr->m_points.size() << std::endl;
    std::cout << "Number of points: " << FineMomentumGrid().ptr->m_points.size() << std::endl;
}


double HubbardAtom::E(const int idx_p, const int idx_w)
{
    return - HubbardAtomParams::MU;
}

double HubbardAtom::E_of_p_coord(const coord_t<0> p_coord)
{
    return - HubbardAtomParams::MU;
}


double HubbardAtom::vertex_4pt_bare(
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, 			  
        const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, 
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       )
{
    return -HubbardAtomParams::UINT;
}



double HubbardAtom::B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return - MomentumGrid().ptr->get_volume() * HubbardAtomParams::UINT
	    #ifdef B_PLUS_F
	    *0.5
            #endif

;
    else
	return 0.;
}

double HubbardAtom::B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    if(idx_m==0 && idx_mp==0) 
	return MomentumGrid().ptr->get_volume() * HubbardAtomParams::UINT
#ifdef B_PLUS_F
	    *0.5
            #endif
;
    else
	return 0.;
}


double HubbardAtom::B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    double val = 0.0;
    if(idx_m==0 && idx_mp==0){        
	val += - MomentumGrid().ptr->get_volume() * HubbardAtomParams::UINT
#ifdef B_PLUS_F
	    *0.5
            #endif
;
	return val;
    }else
	return 0.;
}


double HubbardAtom::F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
#ifdef B_PLUS_F
    return B_sc(idx_W, idx_K, idx_w, idx_m, idx_wp, idx_mp);
#endif
    return 0.0;
}

double HubbardAtom::F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
#ifdef B_PLUS_F
    return B_m(idx_W, idx_K, idx_w, idx_m, idx_wp, idx_mp);
#endif

    return 0.0;
}


double HubbardAtom::F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
#ifdef B_PLUS_F
    return B_d(idx_W, idx_K, idx_w, idx_m, idx_wp, idx_mp);
#endif

    return 0.0;
}


// needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
double HubbardAtom::vertex_local_part_bare(const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp )
{
    return - MomentumGrid().ptr->get_volume() * HubbardAtomParams::UINT
#ifdef B_PLUS_F
	    *0.25
            #endif
;
}

double HubbardAtom::vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out)
{
    return -HubbardAtomParams::UINT
#ifdef B_PLUS_F
	    *0.25
            #endif
;
}
