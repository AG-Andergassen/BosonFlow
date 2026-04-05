#pragma once

#include <models/abstract_model.h>
#include <models/hubbard.h>

class TriangularHubbard : public Hubbard<2>
{
 public:
    TriangularHubbard() = delete;
    
    static const std::string GetName(){return "Triangular Hubbard";};

    // the parameter name-value pairs returned here will be included in the output file
    static const std::vector<std::pair<std::string, double>> GetParamNameValuePairs()
    {
        return {{ "U", TriangularHubbardParams::UINT }, { "V", TriangularHubbardParams::VINT }, { "W", TriangularHubbardParams::WINT }, { "WW", TriangularHubbardParams::WWINT }, { "TP", TriangularHubbardParams::T_PRIME}, { "TPP", TriangularHubbardParams::T_PRIME_PRIME}, {"Mu", TriangularHubbardParams::MU}, {"KDIM", K_DIM}, {"FineMoms", GetFineMomentaCount()}, {"RefdMoms", GetRefinedMomentaCount()}, {"FFShellCount", FORMFACTOR_SHELL_COUNT}};
    };

    static double E(const int idx_k, const int idx_w);

    static double E_of_p_coord(const coord_t<2> p_coord);

    static double Potential_of_K_coord(const coord_t<2> K_coord);

    static double vertex_4pt_bare( 
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out,
	const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out,
	const int s1_in, const int s2_in, const int s1_out, const int s2_out 
			       ); // bare vertex/interaction of the model

    // used for extended parameters
    static double B_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double B_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double B_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double F_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double F_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double F_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    // needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
    static double vertex_local_part_bare( const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    // needed to prevent double counting in the calculation of the full vertex (for the full s.e. flow
    static double vertex_4pt_local_part_bare(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out);


    static void Init();

    static double vertex_DC_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double vertex_DC_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double vertex_DC_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

 protected:    
    static void PrecomputeE();

    static std::vector<double> &PrecomputedE()
    {
	static std::vector<double> precomputed_E;
	return precomputed_E;
    }

    static void PrecomputePotential();

    static std::vector<double> &PrecomputedPotential()
    {
	static std::vector<double> precomputed_potential;
	return precomputed_potential;
    }

};
