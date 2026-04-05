#pragma once

#include <models/anderson_impurity.h>

class AndersonImpurityHolstein : public AndersonImpurity
{
 public:
    AndersonImpurityHolstein() = delete;
    
    static const std::string GetName(){return "Anderson Impurity-Holstein";};

    // the parameter name-value pairs returned here will be included in the output file
    static const std::vector<std::pair<std::string, double>> GetParamNameValuePairs()
    {
        return  {{ "g0", AndersonImpurityHolsteinParams::g0 }, { "OMEGA0", AndersonImpurityHolsteinParams::OMEGA0 }, { "U", AndersonImpurityParams::UINT }, {"Mu", AndersonImpurityParams::MU}, {"KDIM", K_DIM}, {"FineMoms", GetFineMomentaCount()}, {"RefdMoms", GetRefinedMomentaCount()}, {"FFShellCount", FORMFACTOR_SHELL_COUNT}, {"BATH_LEVEL_COUNT", AndersonImpurityParams::BATH_LEVEL_COUNT}};
    };

    static void Init();
  
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

    // used for Free_U flow
    static double U_of_t(const double t);
    static double V_ph_of_t(const double t);
    static double U_of_t_dot(const double t);
    static double V_ph_of_t_dot(const double t);

    // used when Free_U + Temperature flow flag is turned on
    static double L_of_t(const int idx_W, const double t);
    static double L_of_t_dot(const int idx_W, const double t);


    static double B_sc_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double B_d_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double B_m_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double F_sc_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double F_d_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double F_m_of_t( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );


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

    static double F_sc_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double F_d_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    static double F_m_of_t_dot( const int idx_W, const int idx_K, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );


    static double vertex_4pt_bare_of_t_dot( 
	const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out,
	const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out,
	const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t 
			       ); // bare vertex/interaction of the model

    // needed to prevent double counting in the calculation of sbe B_irreducible_vertex_d method
    static double vertex_local_part_bare_of_t_dot( const int idx_W, const int idx_w, const int idx_m, const int idx_wp, const int idx_mp, const double t );

    // needed to prevent double counting in the calculation of the full vertex (for the full s.e. flow
    static double vertex_4pt_local_part_bare_of_t_dot(const int idx_w1_in, const int idx_w2_in, const int idx_w1_out, const int idx_w2_out, const int idx_p1_in, const int idx_p2_in, const int idx_p1_out, const int idx_p2_out, const int s1_in, const int s2_in, const int s1_out, const int s2_out, const double t);



};
