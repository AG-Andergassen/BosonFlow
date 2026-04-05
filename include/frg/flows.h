#pragma once

#include <base_gf_types.h>
#include <models/concrete_available_models.h>
#include <string>

enum class G_FlowSchemeName{None, Interaction, Temperature, Omega, Eberlein, Interpolating, InverseInterpolating, InverseInterpolatingOmega};

enum class U_FlowSchemeName{None, Multiplicative, Free};



// contains all objects to do with cutoff schemes (Lambda dependent G, Single scale propagator S ..etc)
template <typename Model>
class fRGFlowScheme
{
 public:    
    struct FlowParamInfo {
        double t_start;
        double t_final;
        double init_t_step;
    };

    static void UseFlow(G_FlowSchemeName G_flow_name, U_FlowSchemeName U_flow_name);

    static std::string &ChosenFlowSchemeAbbrev()
    {
	static std::string chosen_flow_scheme_abbrev("UNKNOWN_FLOW");
	return chosen_flow_scheme_abbrev;
    }

    static void Set_Ginit_ForInterpolatingFlow(gf_1p_mat_t<Model> Ginit);

    static G_FlowSchemeName &G_ChosenFlowName()
    {
	static G_FlowSchemeName G_chosen_flow_name;
	return G_chosen_flow_name;
    }

    static U_FlowSchemeName &U_ChosenFlowName()
    {
	static U_FlowSchemeName U_chosen_flow_name;
	return U_chosen_flow_name;
    }

    static FlowParamInfo &ChosenFlowParametrizationInfo()
    {
	static FlowParamInfo chosen_flow{0.0, 1.0, 0.05};
	return chosen_flow;
    }


    /**
     *  \brief Definition fRG functions
     */
    static MatQN G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );	            /**< Scale-dependent Greens function, introduce regulator here! */
    static MatQN G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );	            /**< Scale-dependent Greens function, introduce regulator here! */
    static MatQN S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt );             /**< Single-scale propagator */

    static dcomplex static_SG_pp_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_n, const double t);

    static dcomplex static_SG_ph_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_n, const double t);

    static double U_MultiplicativeCutoff( const double t ); // scale dependent (local) bare interaction

    static double U_MultiplicativeCutoff_dot( const double t ); // scale derivative of (local) bare interaction

/**
 * \brief Asymptotics (for large fermionic frequencies) of functions eval_diag_bubble_pp, eval_diag_bubble_ph, eval_diag_bubble_GG_pp and eval_diag_bubble_GG_ph
 */
    static dcomplex asymptotic_SG_pp( const int idx_W, const double t );	                                            /**< Give estimate for the integral 1/2/PI * S(W/2-w-1,Lam) * G(W/2+w,Lam) for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE ) */
    static dcomplex asymptotic_SGpGS_pp( const int idx_W, const double t );	                                        /**< Give estimate for the integral 1/2/PI * [ S(w-W/2,Lam) * G(w+W/2,Lam) + S(w-W/2,Lam) * G(w+W/2,Lam) ] for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE ) */

    static dcomplex asymptotic_SG_ph( const int idx_W, const double t );	                                            /**< Give estimate for the integral 1/2/PI * S(w-W/2,Lam) * G(W/2+w,Lam) for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE ) */
    static dcomplex asymptotic_SGpGS_ph( const int idx_W, const double t );	                                        /**< Give estimate for the integral 1/2/PI * [ S(w-W/2,Lam) * G(w+W/2,Lam) + S(w-W/2,Lam) * G(w+W/2,Lam) ] for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE ) */

    static dcomplex asymptotic_GG_pp( const int idx_W, const double t );	                                            /**< Give estimate for the integral 1/2/PI * G(W/2-w-1,Lam) * G(W/2+w,Lam) for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE ) */
    static dcomplex asymptotic_GG_ph( const int idx_W, const double t );	                                            /**< Give estimate for the integral 1/2/PI * G(w-W/2,Lam) * G(w+W/2,Lam) for the outside region N \ [ -POS_INT_RANGE, POS_INT_RANGE ) */

    static double inv_freq_sum_normalisation( const double t);

    static double LocalBareInteraction(const double t);
    
    static double Beta( const double t);

    static gf_1p_mat_t<Model> &InterpolatingFlow_Ginit()
    {
	static gf_1p_mat_t<Model> interpolating_flow_Ginit(FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count, Model::GetFineMomentaCount(), true);
	return interpolating_flow_Ginit;
    };
    
 private:
    // shorthand for functors with args: w, k, t, selfEn
    using G_func_t = std::function<MatQN (const int, const int, const double, const MatQN&, const double)>; 
    using S_func_t = std::function<MatQN (const int, const int, const double, const MatQN&, const double, const double)>; 
    using Ucutoff_func_t = std::function<double (const double)>; 
    

    using asymptotic_G_func_t = std::function<dcomplex(const int, const double)>;
    using static_bubble_func_t = std::function<dcomplex(const coord_t<Model::dim>, const coord_t<Model::dim>, const int, const int, const double)>;


    static G_func_t &Chosen_G()
    {
	static G_func_t chosen_G;
	return chosen_G;
    }
    static G_func_t &Chosen_G_latt()
    {
	static G_func_t chosen_G_latt;
	return chosen_G_latt;
    }
    static S_func_t &Chosen_S()
    {
	static S_func_t chosen_S;
	return chosen_S;
    }

    static static_bubble_func_t &Chosen_static_SG_pp_integrand()
    {
	static static_bubble_func_t chosen_static_SG_pp_integrand;
	return chosen_static_SG_pp_integrand;
    }

    static static_bubble_func_t &Chosen_static_SG_ph_integrand()
    {
	static static_bubble_func_t chosen_static_SG_ph_integrand;
	return chosen_static_SG_ph_integrand;
    }


    static asymptotic_G_func_t &Chosen_asymptotic_SG_pp()
    {
	static asymptotic_G_func_t chosen_asymptotic_SG_pp;
	return chosen_asymptotic_SG_pp;
    }
    static asymptotic_G_func_t &Chosen_asymptotic_GG_pp()
    {
	static asymptotic_G_func_t chosen_asymptotic_GG_pp;
	return chosen_asymptotic_GG_pp;
    }
    static Ucutoff_func_t &Chosen_U_MultiplicativeCutoff()
    {
	static Ucutoff_func_t chosen_U_multiplicative_cutoff;
	return chosen_U_multiplicative_cutoff;
    }
    static Ucutoff_func_t &Chosen_U_MultiplicativeCutoff_dot()
    {
	static Ucutoff_func_t chosen_U_multiplicative_cutoff_dot;
	return chosen_U_multiplicative_cutoff_dot;
    }


    static MatQN InterpolatingFlow_G(const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN InverseInterpolatingFlow_G(const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN InverseInterpolatingOmegaFlow_G(const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN TemperatureFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN InteractionFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN OmegaFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN EberleinFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );

    static MatQN InterpolatingFlow_G_latt(const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN InverseInterpolatingFlow_G_latt(const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN InverseInterpolatingOmegaFlow_G_latt(const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN TemperatureFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN InteractionFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN OmegaFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN EberleinFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );

    static MatQN InterpolatingFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt );
    static MatQN InverseInterpolatingFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt );
    static MatQN InverseInterpolatingOmegaFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt );
    static MatQN TemperatureFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt );
    static MatQN InteractionFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt );
    static MatQN OmegaFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt );
    static MatQN EberleinFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt );

    static dcomplex InterpolatingFlow_asymptotic_SG_pp(const int idx_W, const double t);
    static dcomplex InverseInterpolatingFlow_asymptotic_SG_pp(const int idx_W, const double t);
    static dcomplex InverseInterpolatingOmegaFlow_asymptotic_SG_pp(const int idx_W, const double t);
    static dcomplex TemperatureFlow_asymptotic_SG_pp(const int idx_W, const double t);
    static dcomplex InteractionFlow_asymptotic_SG_pp(const int idx_W, const double t);
    static dcomplex OmegaFlow_asymptotic_SG_pp(const int idx_W, const double t);
    static dcomplex EberleinFlow_asymptotic_SG_pp(const int idx_W, const double t);
    
    
    static dcomplex InterpolatingFlow_asymptotic_GG_pp(const int idx_W, const double t);
    static dcomplex InverseInterpolatingFlow_asymptotic_GG_pp(const int idx_W, const double t);
    static dcomplex InverseInterpolatingOmegaFlow_asymptotic_GG_pp(const int idx_W, const double t);
    static dcomplex TemperatureFlow_asymptotic_GG_pp(const int idx_W, const double t);
    static dcomplex InteractionFlow_asymptotic_GG_pp(const int idx_W, const double t);
    static dcomplex OmegaFlow_asymptotic_GG_pp(const int idx_W, const double t);
    static dcomplex EberleinFlow_asymptotic_GG_pp(const int idx_W, const double t);
    
    static dcomplex TemperatureFlow_static_SG_pp_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_mp, const double t);
    static dcomplex TemperatureFlow_static_SG_ph_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_mp, const double t);
    static dcomplex InteractionFlow_static_SG_pp_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_mp, const double t);
    static dcomplex InteractionFlow_static_SG_ph_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_mp, const double t);

    static double Actual_UFlow_MultiplicativeCutoff(const double t);
    static double Actual_UFlow_MultiplicativeCutoff_dot(const double t);

    static double U_TrivialCutoff(const double t);
    static double U_TrivialCutoff_dot(const double t);
    
    static MatQN Nonflowing_G(const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN Nonflowing_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu );
    static MatQN Nonflowing_S(const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt );
    static dcomplex Nonflowing_asymptotic_SG_pp(const int idx_W, const double t);
    static dcomplex Nonflowing_asymptotic_GG_pp(const int idx_W, const double t);
};

