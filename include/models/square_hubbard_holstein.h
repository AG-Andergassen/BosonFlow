#pragma once

#include <models/square_hubbard.h>

class SquareHubbardHolstein : public SquareHubbard
{
 public:
    SquareHubbardHolstein() = delete;
    
    static const std::string GetName(){return "Square Hubbard-Holstein";};

    // the parameter name-value pairs returned here will be included in the output file
    static const std::vector<std::pair<std::string, double>> GetParamNameValuePairs()
    {
        return  {{ "g0", SquareHubbardHolsteinParams::g0 }, { "OMEGA0", SquareHubbardHolsteinParams::OMEGA0 }, { "U", SquareHubbardParams::UINT }, { "V", SquareHubbardParams::VINT }, { "TP", SquareHubbardParams::T_PRIME}, {"Mu", SquareHubbardParams::MU}, {"KDIM", K_DIM}, {"FineMoms", GetFineMomentaCount()}, {"RefdMoms", GetRefinedMomentaCount()}, {"FFShellCount", FORMFACTOR_SHELL_COUNT}};
    };

    static void Init();
  
    static std::vector<double> &PrecomputedAcousticPhononDispersion()
    {
	static std::vector<double> precomputed_acoustic_phonon_dispersion;
	return precomputed_acoustic_phonon_dispersion;
    }


    static void PrecomputeAcousticPhononDispersion();

    static double vertex_4pt_bare( 
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out,
	const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out,
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       ); // bare vertex/interaction of the model

    // used for extended parameters
    static double B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    // used for extended parameters
    static double F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double vertex_local_part_bare( const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out);

    static double vertex_DC( 
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out,
	const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out,
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       ); 
  
    static double vertex_DC_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double vertex_DC_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double vertex_DC_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

};
