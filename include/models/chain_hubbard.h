#pragma once

#include <models/abstract_model.h>
#include <models/hubbard.h>

class ChainHubbard : public Hubbard<1>
{
 public:
    ChainHubbard() = delete;
    
    static const std::string GetName(){return "Hubbard Chain";};

    // the parameter name-value pairs returned here will be included in the output file
    static const std::vector<std::pair<std::string, double>> GetParamNameValuePairs()
    {
        return {{ "U", ChainHubbardParams::UINT }, { "TP", ChainHubbardParams::T_PRIME}, {"Mu", ChainHubbardParams::MU}, {"FineMoms", GetFineMomentaCount()}, {"RefdMoms", GetRefinedMomentaCount()}, {"FFShellCount", FORMFACTOR_SHELL_COUNT}};
    };

    static double E(const int idx_k, const int idx_w);

    static double E_of_p_coord(coord_t<1> p_coord);


    //    constexpr unsigned get_innerbox_bosonic_frequencies_count();
    //    constexpr unsigned get_innerbox_fermionic_frequencies_count();

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


    // used for Free_U flow
    static double U_of_t(const double t);
    static double V_of_t(const double t);
    static double U_of_t_dot(const double t);
    static double V_of_t_dot(const double t);

    static double B_sc_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double B_d_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double B_m_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double vertex_4pt_bare_of_t( 
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out,
	const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out,
	const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t 
			       ); // bare vertex/interaction of the model

    // needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
    static double vertex_local_part_bare_of_t( const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    // needed to prevent double counting in the calculation of the full vertex (for the full s.e. flow
    static double vertex_4pt_local_part_bare_of_t(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t);

    static double B_sc_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double B_d_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double B_m_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double vertex_4pt_bare_of_t_dot( 
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out,
	const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out,
	const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t 
			       ); // bare vertex/interaction of the model

    // needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
    static double vertex_local_part_bare_of_t_dot( const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    // needed to prevent double counting in the calculation of the full vertex (for the full s.e. flow
    static double vertex_4pt_local_part_bare_of_t_dot(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t);

	static double vertex_DC_sc( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double vertex_DC_d( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static double vertex_DC_m( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp );

    static void Init();

 protected:    
    static void PrecomputeE();

    static std::vector<double> &PrecomputedE()
    {
	static std::vector<double> precomputed_E;
	return precomputed_E;
    }

};
