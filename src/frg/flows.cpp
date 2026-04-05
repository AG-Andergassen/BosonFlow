#include <frg/flows.h>
#include <frequencies/matsubara_space.h> // for w_val and W_val
#include <Eigen/LU> // for matrix inverse
#include <mymath.h>


template <typename Model>
void fRGFlowScheme<Model>::UseFlow(const G_FlowSchemeName G_flow, const U_FlowSchemeName U_flow)
{
    G_ChosenFlowName() = G_flow;
    U_ChosenFlowName() = U_flow;

    FlowParamInfo flow_parameter;
    flow_parameter.t_start = 0.0;
    flow_parameter.t_final = 1.0;
    flow_parameter.init_t_step = 0.05;

    if (G_ChosenFlowName() == G_FlowSchemeName::None && U_ChosenFlowName() == U_FlowSchemeName::None){
	throw std::invalid_argument("An fRG flow requires a cutoff in either G or U!");
    }

    switch (G_flow)
    {
    case G_FlowSchemeName::None:
std::cout << "G none" << std::endl;
    flow_parameter.t_start = 0.0;
    flow_parameter.t_final = 1.0;
    flow_parameter.init_t_step = 0.05;
    ChosenFlowSchemeAbbrev() = "NO_GFL";
	Chosen_G() = Nonflowing_G;
	Chosen_G_latt() = Nonflowing_G_latt;
	Chosen_S() = Nonflowing_S;
	Chosen_asymptotic_SG_pp() = Nonflowing_asymptotic_SG_pp;
	Chosen_asymptotic_GG_pp() = Nonflowing_asymptotic_GG_pp;

	// not implemented
	Chosen_static_SG_pp_integrand() = InteractionFlow_static_SG_pp_integrand;
	Chosen_static_SG_ph_integrand() = InteractionFlow_static_SG_ph_integrand;
	break;
    case G_FlowSchemeName::Interpolating:
    flow_parameter.t_start = 0.0;
    flow_parameter.t_final = 1.0;
    flow_parameter.init_t_step = 0.01;
	ChosenFlowSchemeAbbrev() = "INTRP";
	Chosen_G() = InterpolatingFlow_G;
	Chosen_G_latt() = InterpolatingFlow_G_latt;
	Chosen_S() = InterpolatingFlow_S;
	Chosen_asymptotic_SG_pp() = InterpolatingFlow_asymptotic_SG_pp;
	Chosen_asymptotic_GG_pp() = InterpolatingFlow_asymptotic_GG_pp;
	// not implemented
	Chosen_static_SG_pp_integrand() = InteractionFlow_static_SG_pp_integrand;
	Chosen_static_SG_ph_integrand() = InteractionFlow_static_SG_ph_integrand;
	break;
    case G_FlowSchemeName::InverseInterpolating:
    flow_parameter.t_start = 0.0;
    flow_parameter.t_final = 1.0;
    flow_parameter.init_t_step = 0.01;
	ChosenFlowSchemeAbbrev() = "INVERSE_INTRP";
	Chosen_G() = InverseInterpolatingFlow_G;
	Chosen_G_latt() = InverseInterpolatingFlow_G_latt;
	Chosen_S() = InverseInterpolatingFlow_S;
	Chosen_asymptotic_SG_pp() = InverseInterpolatingFlow_asymptotic_SG_pp;
	Chosen_asymptotic_GG_pp() = InverseInterpolatingFlow_asymptotic_GG_pp;
	// not implemented
	Chosen_static_SG_pp_integrand() = InteractionFlow_static_SG_pp_integrand;
	Chosen_static_SG_ph_integrand() = InteractionFlow_static_SG_ph_integrand;
	break;
    case G_FlowSchemeName::InverseInterpolatingOmega:
    flow_parameter.t_start = 5.0;
    flow_parameter.t_final = -1.0;
    flow_parameter.init_t_step = -1.0;
	ChosenFlowSchemeAbbrev() = "INVERSE_INTRP_OMEGA";
	Chosen_G() = InverseInterpolatingOmegaFlow_G;
	Chosen_G_latt() = InverseInterpolatingOmegaFlow_G_latt;
	Chosen_S() = InverseInterpolatingOmegaFlow_S;
	Chosen_asymptotic_SG_pp() = OmegaFlow_asymptotic_SG_pp;
	Chosen_asymptotic_GG_pp() = OmegaFlow_asymptotic_GG_pp;
	// not implemented
	Chosen_static_SG_pp_integrand() = InteractionFlow_static_SG_pp_integrand;
	Chosen_static_SG_ph_integrand() = InteractionFlow_static_SG_ph_integrand;
	break;
    case G_FlowSchemeName::Interaction:
    flow_parameter.t_start = 0.0;
    flow_parameter.t_final = 1.0;
    flow_parameter.init_t_step = 0.25;
	ChosenFlowSchemeAbbrev() = "INTFL";
	Chosen_G() = InteractionFlow_G;
	Chosen_G_latt() = InteractionFlow_G_latt;
	Chosen_S() = InteractionFlow_S;
	Chosen_asymptotic_SG_pp() = InteractionFlow_asymptotic_SG_pp;
	Chosen_asymptotic_GG_pp() = InteractionFlow_asymptotic_GG_pp;
	Chosen_static_SG_pp_integrand() = InteractionFlow_static_SG_pp_integrand;
	Chosen_static_SG_ph_integrand() = InteractionFlow_static_SG_ph_integrand;
	break;
    case G_FlowSchemeName::Temperature:
    flow_parameter.t_start = 10.0;
    flow_parameter.t_final = log(1.0/BETA)/LN_10;
    flow_parameter.init_t_step = -2.0;
	ChosenFlowSchemeAbbrev() = "TEMPFL";
	Chosen_G() = TemperatureFlow_G;
	Chosen_G_latt() = TemperatureFlow_G_latt;
	Chosen_S() = TemperatureFlow_S;
	Chosen_asymptotic_SG_pp() = TemperatureFlow_asymptotic_SG_pp;
	Chosen_asymptotic_GG_pp() = TemperatureFlow_asymptotic_GG_pp;
	Chosen_static_SG_pp_integrand() = TemperatureFlow_static_SG_pp_integrand;
	Chosen_static_SG_ph_integrand() = TemperatureFlow_static_SG_ph_integrand;
	break;
    case G_FlowSchemeName::Omega:
    flow_parameter.t_start = 5.0;
    flow_parameter.t_final = -1.0;
    flow_parameter.init_t_step = -1.0;
	ChosenFlowSchemeAbbrev() = "OMFL";
	Chosen_G() = OmegaFlow_G;
	Chosen_G_latt() = OmegaFlow_G_latt;
	Chosen_S() = OmegaFlow_S;
	Chosen_asymptotic_SG_pp() = OmegaFlow_asymptotic_SG_pp;
	Chosen_asymptotic_GG_pp() = OmegaFlow_asymptotic_GG_pp;
	break;
    case G_FlowSchemeName::Eberlein:
    flow_parameter.t_start = 10.0;
    flow_parameter.t_final = -10.0;
    flow_parameter.init_t_step = -2.0;
	ChosenFlowSchemeAbbrev() = "EBFL";
	Chosen_G() = EberleinFlow_G;
	Chosen_G_latt() = EberleinFlow_G_latt;
	Chosen_S() = EberleinFlow_S;
	Chosen_asymptotic_SG_pp() = EberleinFlow_asymptotic_SG_pp;
	Chosen_asymptotic_GG_pp() = EberleinFlow_asymptotic_GG_pp;
#ifdef STATIC_CALCULATION
	throw std::invalid_argument("Static fRG with Eberlein flow not implemented.");
#endif
	break;
    default:
	;
    }
    
    // TODO: kind of redundant..
    switch (U_flow){
    case U_FlowSchemeName::None:
	Chosen_U_MultiplicativeCutoff() = U_TrivialCutoff;
	Chosen_U_MultiplicativeCutoff_dot() = U_TrivialCutoff_dot;
	break;
    case U_FlowSchemeName::Multiplicative:
	std::cout << "U multiplicative" << std::endl;
    if (G_flow == G_FlowSchemeName::None) {
        ChosenFlowSchemeAbbrev() = "ACTUFL";
    } else {
        ChosenFlowSchemeAbbrev() += "_ACTUFL";
    }
	Chosen_U_MultiplicativeCutoff() = Actual_UFlow_MultiplicativeCutoff;
	Chosen_U_MultiplicativeCutoff_dot() = Actual_UFlow_MultiplicativeCutoff_dot;
	break;
    case U_FlowSchemeName::Free:
    if (G_flow == G_FlowSchemeName::None) {
        ChosenFlowSchemeAbbrev() = "FREEUFL";
    } else {
        ChosenFlowSchemeAbbrev() += "_FREEUFL";
    }
	Chosen_U_MultiplicativeCutoff() = Actual_UFlow_MultiplicativeCutoff;
	Chosen_U_MultiplicativeCutoff_dot() = Actual_UFlow_MultiplicativeCutoff_dot;
	break;
    default:
	;
    }

    if (G_flow == G_FlowSchemeName::None && U_flow != U_FlowSchemeName::None) {
    flow_parameter.t_start = 1e-8;
    flow_parameter.t_final = 1.0;
    flow_parameter.init_t_step = 0.05;
    }

    ChosenFlowParametrizationInfo() = flow_parameter;

}


template <typename Model>
void fRGFlowScheme<Model>::Set_Ginit_ForInterpolatingFlow(gf_1p_mat_t<Model> Ginit)
{
    InterpolatingFlow_Ginit() = Ginit;
}

template <typename Model>
MatQN fRGFlowScheme<Model>::G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    return Chosen_G()(idx_w, idx_p, t, selfEn, delta_mu);
}


template <typename Model>
MatQN fRGFlowScheme<Model>::G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    return Chosen_G_latt()(idx_w, idx_p, t, selfEn, delta_mu);
}

template <typename Model>
MatQN fRGFlowScheme<Model>::S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt )
{
    return Chosen_S()(idx_w, idx_p, t, selfEn, delta_mu, d_delta_mu_over_dt);
}

template <typename Model>
double fRGFlowScheme<Model>::U_MultiplicativeCutoff( const double t )
{
    return Chosen_U_MultiplicativeCutoff()(t);
}

template <typename Model>
double fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot( const double t )
{
    return Chosen_U_MultiplicativeCutoff_dot()(t);
}


template <typename Model>
dcomplex fRGFlowScheme<Model>::static_SG_pp_integrand( const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_n, const double t )
{
    return Chosen_static_SG_pp_integrand()(p_coord, q_coord, idx_m, idx_n, t);
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::static_SG_ph_integrand( const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_n, const double t )
{
    return Chosen_static_SG_ph_integrand()(p_coord, q_coord, idx_m, idx_n, t);
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::asymptotic_SG_pp(const int idx_W, const double t)
{
    return 0.0;
    return Chosen_asymptotic_SG_pp()(idx_W, t);
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::asymptotic_GG_pp(const int idx_W, const double t)
{
    return 0.0;
    return Chosen_asymptotic_GG_pp()(idx_W, t);
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::asymptotic_SGpGS_pp( const int idx_W, const double t )
{
    return 0.0;
    return 2.0 * asymptotic_SG_pp( idx_W, t );
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::asymptotic_SG_ph( const int idx_W, const double t )
{
    return 0.0;
    return -asymptotic_SG_pp( idx_W, t );
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::asymptotic_SGpGS_ph( const int idx_W, const double t )
{
    return 0.0;
    return 2.0 * asymptotic_SG_ph( idx_W, t );
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::asymptotic_GG_ph( const int idx_W, const double t )
{
    return 0.0;
    return -asymptotic_GG_pp( idx_W, t );
}


// Flow dependent

template <typename Model>
MatQN fRGFlowScheme<Model>::InterpolatingFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0_lattice_inv, G0;
    double w = w_val(idx_w);
    G0_lattice_inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    G0 << t * G0_lattice_inv.inverse() + (1-t) * InterpolatingFlow_Ginit()[idx_w][idx_p];
    return ( G0.inverse() - selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::InterpolatingFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    return InterpolatingFlow_G(idx_w, idx_p, t, selfEn, delta_mu); 
}

template <typename Model>
MatQN fRGFlowScheme<Model>::InterpolatingFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt )
{
    MatQN G0_lattice_inv, G0, G, G0inv, G0_lattice;
    double w = w_val(idx_w);
    G0_lattice_inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    G0_lattice << G0_lattice_inv.inverse();
    G0 << t * G0_lattice_inv.inverse() + (1-t) * InterpolatingFlow_Ginit()[idx_w][idx_p];
    G0inv << G0.inverse();
    G << ( G0inv - selfEn ).inverse();

    MatQN d_delta_mu_over_dt_mat;
    d_delta_mu_over_dt_mat << d_delta_mu_over_dt;

    return G * G0inv * (G0_lattice - InterpolatingFlow_Ginit()[idx_w][idx_p] - t * G0_lattice * d_delta_mu_over_dt_mat * G0_lattice ) * G0inv * G ;
}


template <typename Model>
dcomplex fRGFlowScheme<Model>::InterpolatingFlow_asymptotic_SG_pp(const int idx_W, const double t)
{
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    if( idx_W == 0 ) 
	return (BETA*t)/(2.0*pow(M_PI,2.0)*PIR);
   
    return (BETA*t*atanh(W/(2.0*PIR)))/(pow(M_PI,2.0)*W);
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::InterpolatingFlow_asymptotic_GG_pp(const int idx_W, const double t)
{
    return t * InterpolatingFlow_asymptotic_SG_pp( idx_W, t );
}

// inverse interpolating flow

template <typename Model>
MatQN fRGFlowScheme<Model>::InverseInterpolatingFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0_lattice_inv, G0_inv;
    double w = w_val(idx_w);
    G0_lattice_inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    G0_inv << t * G0_lattice_inv + (1-t) * InterpolatingFlow_Ginit()[idx_w][idx_p].inverse();
    return ( G0_inv - selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::InverseInterpolatingFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    return InverseInterpolatingFlow_G(idx_w, idx_p, t, selfEn, delta_mu); 
}

template <typename Model>
MatQN fRGFlowScheme<Model>::InverseInterpolatingFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt )
{
  MatQN G0_lattice_inv, G0_inv, G;
  double w = w_val(idx_w);
  G0_lattice_inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
  G0_inv << t * G0_lattice_inv + (1-t) * InterpolatingFlow_Ginit()[idx_w][idx_p].inverse();
  G <<  ( G0_inv - selfEn ).inverse();
  MatQN d_delta_mu_over_dt_mat;
  d_delta_mu_over_dt_mat << d_delta_mu_over_dt;
  
  return - G * ( (G0_lattice_inv + t*d_delta_mu_over_dt_mat) - InterpolatingFlow_Ginit()[idx_w][idx_p].inverse() ) * G;
}


template <typename Model>
dcomplex fRGFlowScheme<Model>::InverseInterpolatingFlow_asymptotic_SG_pp(const int idx_W, const double t)
{
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    if( idx_W == 0 ) 
	return (BETA*t)/(2.0*pow(M_PI,2.0)*PIR);
   
    return (BETA*t*atanh(W/(2.0*PIR)))/(pow(M_PI,2.0)*W);
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::InverseInterpolatingFlow_asymptotic_GG_pp(const int idx_W, const double t)
{
    return t * InverseInterpolatingFlow_asymptotic_SG_pp( idx_W, t );
}


// inverse interpolating OMEGA flow

template <typename Model>
MatQN fRGFlowScheme<Model>::InverseInterpolatingOmegaFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0_lattice_inv, G0_inv;
    double w = w_val(idx_w);
    G0_lattice_inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    const auto &flow = ChosenFlowParametrizationInfo();
    constexpr double AA = 1.0;
    constexpr double AAA = 0.0;
    const double A = -1.0 * (AA * exp(LN_10 * flow.t_final) + AAA) / flow.t_final / flow.t_final / flow.t_final;
    const double Lam = AA * exp( t * LN_10 ) + A * t * t * t + AAA;
    double regulator = w*w / ( Lam*Lam + w*w );

    G0_inv << regulator * G0_lattice_inv + (1-regulator) * InterpolatingFlow_Ginit()[idx_w][idx_p].inverse();
    return ( G0_inv - selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::InverseInterpolatingOmegaFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    return InverseInterpolatingOmegaFlow_G(idx_w, idx_p, t, selfEn, delta_mu); 
}

template <typename Model>
MatQN fRGFlowScheme<Model>::InverseInterpolatingOmegaFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt )
{
  MatQN G0_lattice_inv, G0_inv, G;
  double w = w_val(idx_w);
  G0_lattice_inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    const auto &flow = ChosenFlowParametrizationInfo();
        constexpr double AA = 1.0;
        constexpr double AAA = 0.0;
        const double A = -1.0 * (AA * exp(LN_10 * flow.t_final) + AAA) / flow.t_final / flow.t_final / flow.t_final;
  
        double Lam = AA * exp( t * LN_10 ) + A * t * t * t + AAA;
        double dLam_over_dt = AA * LN_10 * exp( t * LN_10 ) + 3. * A * t * t;
  double regulator = w*w / ( Lam*Lam + w*w );
  double d_regulator_over_dLam = -2 *Lam*w*w/((Lam*Lam + w*w)*(Lam*Lam + w*w));
  double d_regulator_over_dt = d_regulator_over_dLam * dLam_over_dt;
  G0_inv << regulator * G0_lattice_inv + (1-regulator) * InterpolatingFlow_Ginit()[idx_w][idx_p].inverse();
  
  G <<  ( G0_inv - selfEn ).inverse();
  MatQN d_delta_mu_over_dt_mat;
  d_delta_mu_over_dt_mat << d_delta_mu_over_dt;

  // -G d_t(G_0^{-1}) G
  return - G * ( (d_regulator_over_dt*G0_lattice_inv + regulator*d_delta_mu_over_dt_mat) - d_regulator_over_dt*InterpolatingFlow_Ginit()[idx_w][idx_p].inverse() ) * G;
}


template <typename Model>
dcomplex fRGFlowScheme<Model>::InverseInterpolatingOmegaFlow_asymptotic_SG_pp(const int idx_W, const double t)
{
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    if( idx_W == 0 ) 
	return (BETA*t)/(2.0*pow(M_PI,2.0)*PIR);
   
    return (BETA*t*atanh(W/(2.0*PIR)))/(pow(M_PI,2.0)*W);
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::InverseInterpolatingOmegaFlow_asymptotic_GG_pp(const int idx_W, const double t)
{
    return t * InverseInterpolatingOmegaFlow_asymptotic_SG_pp( idx_W, t );
}



// interaction flow

template <typename Model>
MatQN fRGFlowScheme<Model>::InteractionFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    return t * ( G0inv - t * selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::InteractionFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    return ( G0inv - t * selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::InteractionFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;

    MatQN G = ( G0inv - t * selfEn ).inverse();
    MatQN d_delta_mu_over_dt_mat;
    d_delta_mu_over_dt_mat << d_delta_mu_over_dt;

    return G * (G0inv -  t*d_delta_mu_over_dt_mat) * G;
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::InteractionFlow_asymptotic_SG_pp(const int idx_W, const double t)
{
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    if( idx_W == 0 ) 
	return (BETA*t)/(2.0*pow(M_PI,2.0)*PIR);
   
    return (BETA*t*atanh(W/(2.0*PIR)))/(pow(M_PI,2.0)*W);
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::InteractionFlow_asymptotic_GG_pp(const int idx_W, const double t)
{
    return t * InteractionFlow_asymptotic_SG_pp( idx_W, t );
}


template <typename Model>
dcomplex fRGFlowScheme<Model>::InteractionFlow_static_SG_ph_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_n, const double t)
{
    // todo: add static self-energy
    double T = 1.0/BETA;

    double E_of_q_plus_p = Model::E_of_p_coord(q_coord + p_coord);
    double E_of_q = Model::E_of_p_coord(q_coord);

    double denominator = E_of_q_plus_p - E_of_q;

    double kernel;
    if (std::abs(denominator) > 1e-3){
	kernel = (nF(E_of_q_plus_p, T) - nF(E_of_q, T))/ denominator;
    }else{
	kernel = std::signbit(denominator) * d_nF_over_dz(E_of_q, T);
    }   
    
    return 2 * t * kernel * 
	Model::GetFormFactor(idx_n).evaluate_fourier_transform(p_coord)*
	std::conj(Model::GetFormFactor(idx_m).evaluate_fourier_transform(p_coord));
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::InteractionFlow_static_SG_pp_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_n, const double t)
{
    // todo: add static self-energy
    double T = 1.0/BETA;
    double E_of_q_plus_p = Model::E_of_p_coord(q_coord + p_coord);
    double E_of_minus_q = Model::E_of_p_coord(-q_coord);
    
    double denominator = E_of_q_plus_p + E_of_minus_q;

    double kernel;
    if (std::abs(denominator) > 1e-3){
	kernel = (nF(E_of_q_plus_p, T) + nF(E_of_minus_q, T))/ denominator;
    }else{
	kernel = (nF(E_of_q_plus_p, T) + nF(E_of_minus_q, T))/ (1e-3);
    }   

    return - 2 * t * kernel * 
	Model::GetFormFactor(idx_n).evaluate_fourier_transform(p_coord)*
	std::conj(Model::GetFormFactor(idx_m).evaluate_fourier_transform(p_coord));
}




template <typename Model>
MatQN fRGFlowScheme<Model>::TemperatureFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    double T = exp( t * LN_10 );  
    double Sqrt_T = exp( t * LN_10 * 0.5 );  

    MatQN G0inv;
    double w = w_val(idx_w) * BETA * T;
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;

    return Sqrt_T  * ( G0inv - Sqrt_T * selfEn ).inverse();    
}

template <typename Model>
MatQN fRGFlowScheme<Model>::TemperatureFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu)
{
    double T = exp( t * LN_10 ) ;  
    double Sqrt_T = exp( t * LN_10 * 0.5 ) ;  

    MatQN G0inv;
    double w = w_val(idx_w) * BETA * T;
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    
    return ( G0inv - Sqrt_T * selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::TemperatureFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt )
{
    double T = exp( t * LN_10 ) ;  
    double Tpow3o2 = exp( t * LN_10 * 1.5 ) ;  
    double Sqrt_T = exp( t * LN_10 * 0.5 ) ;  
    double dTdt = LN_10 * T;
    MatQN GMat = TemperatureFlow_G( idx_w, idx_p, t, selfEn, delta_mu);
    MatQN G0inv_plus;

    
    G0inv_plus << I * w_val(idx_w) * BETA * T + Model::E(idx_p, idx_w) - delta_mu;

    MatQN d_delta_mu_over_dt_mat;
    d_delta_mu_over_dt_mat << d_delta_mu_over_dt;

    return GMat * (- 0.5 *  G0inv_plus * dTdt/Tpow3o2 - Sqrt_T*d_delta_mu_over_dt_mat ) * GMat ; 
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::TemperatureFlow_asymptotic_SG_pp(const int idx_W, const double t)
{
    double T = exp( t * LN_10 ) ;  
    double Sqrt_T = exp( t * LN_10 * 0.5 ) ;  
    double dTdt = LN_10 * T;
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    //if( idx_W == 0 ) return dTdt/(-4.0*pow(M_PI,2.0)*PIR*T);
    if( idx_W == 0 ) return dTdt/(-4.0*pow(M_PI,2.0)*PIR*T*T);

    //return (LN_10*atanh(W/(2.0*PIR)))/(-2.0*pow(M_PI,2.0)*W);
    return (LN_10*atanh(W/(2.0*PIR)))/(-2.0*pow(M_PI,2.0)*W*T);
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::TemperatureFlow_asymptotic_GG_pp(const int idx_W, const double t)
{
    double T = exp( t * LN_10 ) ;  
    double Sqrt_T = exp( t * LN_10 * 0.5 ) ;  
    double dTdt = LN_10 * T;

    return -2. * T * TemperatureFlow_asymptotic_SG_pp( idx_W, t )/dTdt;
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::TemperatureFlow_static_SG_ph_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_n, const double t)
{
    // todo: add static self-energy
    double T = exp( t * LN_10 ) ;  
    double Sqrt_T = exp( t * LN_10 * 0.5 ) ;  
    double dTdt = LN_10 * T;

    double E_of_q_plus_p = Model::E_of_p_coord(q_coord + p_coord);
    double E_of_q = Model::E_of_p_coord(q_coord);

    double denominator = E_of_q_plus_p - E_of_q;

    double kernel;
    if (std::abs(denominator) > 1e-3){
	kernel = (d_nF_over_dT(E_of_q_plus_p, T) - d_nF_over_dT(E_of_q, T))/ denominator;
    }else{
	kernel = std::signbit(denominator) * d2_nF_over_dTdz(E_of_q, T);
    }   
    
    return kernel * 
	Model::GetFormFactor(idx_n).evaluate_fourier_transform(p_coord)*
	std::conj(Model::GetFormFactor(idx_m).evaluate_fourier_transform(p_coord));
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::TemperatureFlow_static_SG_pp_integrand(const coord_t<Model::dim> p_coord, const coord_t<Model::dim> q_coord, const int idx_m, const int idx_n, const double t)
{
    // todo: add static self-energy
    double T = exp( t * LN_10 ) ;  
    double Sqrt_T = exp( t * LN_10 * 0.5 ) ;  
    double dTdt = LN_10 * T;

    double E_of_q_plus_p = Model::E_of_p_coord(q_coord + p_coord);
    double E_of_minus_q = Model::E_of_p_coord(-q_coord);
    
    double denominator = E_of_q_plus_p + E_of_minus_q;

    double kernel;
    if (std::abs(denominator) > 1e-3){
	kernel = (d_nF_over_dT(E_of_q_plus_p, T) + d_nF_over_dT(E_of_minus_q, T))/ denominator;
    }else{
	kernel = std::signbit(denominator) * d2_nF_over_dTdz(E_of_minus_q, T);
    }   

    return - kernel * 
	Model::GetFormFactor(idx_n).evaluate_fourier_transform(p_coord)*
	std::conj(Model::GetFormFactor(idx_m).evaluate_fourier_transform(p_coord));
}


template <typename Model>
MatQN fRGFlowScheme<Model>::OmegaFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;

    const auto &flow = ChosenFlowParametrizationInfo();
    constexpr double AA = 1.0;
    constexpr double AAA = 0.0;
    const double A = -1.0 * (AA * exp(LN_10 * flow.t_final) + AAA) / flow.t_final / flow.t_final / flow.t_final;
    double Lam = AA * exp( t * LN_10 ) + A * t * t * t + AAA;
    double regulator = w*w / ( Lam*Lam + w*w );

    return ( G0inv / regulator - selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::OmegaFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;

    const auto &flow = ChosenFlowParametrizationInfo();
    constexpr double AA = 1.0;
    constexpr double AAA = 0.0;
    const double A = -1.0 * (AA * exp(LN_10 * flow.t_final) + AAA) / flow.t_final / flow.t_final / flow.t_final;
    double Lam = AA * exp( t * LN_10 ) + A * t * t * t + AAA;
    double regulator = w*w / ( Lam*Lam + w*w );

    return ( G0inv - regulator * selfEn ).inverse();

}

template <typename Model>
MatQN fRGFlowScheme<Model>::OmegaFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt  )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;

    const auto &flow = ChosenFlowParametrizationInfo();
    constexpr double AA = 1.0;
    constexpr double AAA = 0.0;
    const double A = -1.0 * (AA * exp(LN_10 * flow.t_final) + AAA) / flow.t_final / flow.t_final / flow.t_final;
    double Lam = AA * exp( t * LN_10 ) + A * t * t * t + AAA;
    double dLamdt = AA * LN_10 * exp( t * LN_10 ) + 3. * A * t * t;
    MatQN GMat = OmegaFlow_G( idx_w, idx_p, t, selfEn, delta_mu );

    double regulator = w*w / ( Lam*Lam + w*w );
    
    MatQN d_delta_mu_over_dt_mat;
    d_delta_mu_over_dt_mat << d_delta_mu_over_dt;

    return GMat * (- 2.0 * Lam / w / w * G0inv * dLamdt -  regulator*d_delta_mu_over_dt_mat) * GMat; 
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::OmegaFlow_asymptotic_SG_pp(const int idx_W, const double t)
{
    const auto &flow = ChosenFlowParametrizationInfo();
    constexpr double AA = 1.0;
    constexpr double AAA = 0.0;
    const double A = -1.0 * (AA * exp(LN_10 * flow.t_final) + AAA) / flow.t_final / flow.t_final / flow.t_final;
    double Lam = AA * exp( t * LN_10 ) + A * t * t * t + AAA;
    double dLamdt = AA * LN_10 * exp( t * LN_10 ) + 3. * A * t * t;
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    if ( Lam < 0.1 ){
	if ( idx_W == 0 )
        return (-(pow(BETA,3.0)*Lam)/(12.0*pow(M_PI,4.0)*pow(PIR,3.0))) * dLamdt;
	return (-(pow(BETA,3.0)*Lam*(-16.0*pow(PIR,3.0)*W + 12.0*PIR*pow(W,3.0) + pow(-4.0*pow(PIR,2.0) + pow(W,2.0),2.0)*atanh((4.0*PIR*W)/(4.0*pow(PIR,2.0) + pow(W,2.0)))))/
        (4.0*pow(M_PI,4.0)*pow(W,3.0)*pow(-4.0*pow(PIR,2.0) + pow(W,2.0),2.0))) * dLamdt;
    }

    return (-(-4.0*M_PI*PIR*BETA*Lam*(pow(M_PI,2.0)*(4.0*pow(PIR,2.0) - 5.0*pow(W,2.0)) - pow(BETA,2.0)*pow(Lam,2.0))*(pow(M_PI,2.0)*pow(W,2.0) + pow(BETA,2.0)*pow(Lam,2.0)) +
          (-(M_PI*W) + BETA*Lam)*(M_PI*W + BETA*Lam)*(pow(M_PI,2.0)*pow(-2.0*PIR + W,2.0) + pow(BETA,2.0)*pow(Lam,2.0))*(pow(M_PI,2.0)*pow(2.0*PIR + W,2.0) + pow(BETA,2.0)*pow(Lam,2.0))*
          (atan((BETA*Lam)/(2.0*M_PI*PIR - M_PI*W)) + atan((BETA*Lam)/(2.0*M_PI*PIR + M_PI*W))) + 2.0*M_PI*W*BETA*Lam*
          (pow(M_PI,4.0)*pow(-4.0*pow(PIR,2.0) + pow(W,2.0),2.0) + 2.0*pow(M_PI,2.0)*(4.0*pow(PIR,2.0) + pow(W,2.0))*pow(BETA,2.0)*pow(Lam,2.0) + pow(BETA,4.0)*pow(Lam,4.0))*
          atanh((4.0*pow(M_PI,2.0)*PIR*W)/(pow(M_PI,2.0)*(4.0*pow(PIR,2.0) + pow(W,2.0)) + pow(BETA,2.0)*pow(Lam,2.0))))/
        (8.0*M_PI*pow(BETA,2.0)*pow((pow(M_PI,2.0)*pow(W,2.0))/pow(BETA,2.0) + pow(Lam,2.0),2.0)*(pow(M_PI,2.0)*pow(-2.0*PIR + W,2.0) + pow(BETA,2.0)*pow(Lam,2.0))*(pow(M_PI,2.0)*pow(2.0*PIR + W,2.0) + pow(BETA,2.0)*pow(Lam,2.0)))
	    ) * dLamdt;
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::OmegaFlow_asymptotic_GG_pp(const int idx_W, const double t)
{
    const auto &flow = ChosenFlowParametrizationInfo();
    constexpr double AA = 1.0;
    constexpr double AAA = 0.0;
    const double A = -1.0 * (AA * exp(LN_10 * flow.t_final) + AAA) / flow.t_final / flow.t_final / flow.t_final;
    double Lam = AA * exp( t * LN_10 ) + A * t * t * t + AAA;
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    if ( idx_W == 0 ){
        if ( Lam < 0.1 ) return BETA/(2.0*pow(M_PI,2.0)*PIR) - (pow(BETA,3.0)*pow(Lam,2.0))/(12.0*pow(M_PI,4.0)*pow(PIR,3.0));

        return (2.0*M_PI*PIR*BETA*Lam + (4.0*pow(M_PI,2.0)*pow(PIR,2.0) + pow(BETA,2.0)*pow(Lam,2.0))*atan((BETA*Lam)/(2.0*M_PI*PIR)))/(8.0*pow(M_PI,3.0)*pow(PIR,2.0)*Lam + 2.0*M_PI*pow(BETA,2.0)*pow(Lam,3.0));
    }

    return (BETA*(M_PI*W*BETA*Lam*(atan((BETA*Lam)/(2.0*M_PI*PIR - M_PI*W)) + atan((BETA*Lam)/(2.0*M_PI*PIR + M_PI*W))) +
	       (2.0*pow(M_PI,2.0)*pow(W,2.0) + pow(BETA,2.0)*pow(Lam,2.0))*atanh((4.0*pow(M_PI,2.0)*PIR*W)/(pow(M_PI,2.0)*(4.0*pow(PIR,2.0) + pow(W,2.0)) + pow(BETA,2.0)*pow(Lam,2.0)))))/
	(4.0*pow(M_PI,2.0)*(pow(M_PI,2.0)*pow(W,3.0) + W*pow(BETA,2.0)*pow(Lam,2.0)));
}

template <typename Model>
MatQN fRGFlowScheme<Model>::EberleinFlow_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    double w = w_val(idx_w);

    double Lam = exp( t * LN_10 );
    double wt = sgn( w ) * sqrt( w*w + Lam*Lam );

    MatQN G0inv;    
    G0inv << I*wt - Model::E(idx_p, idx_w) + delta_mu;


    return ( G0inv - selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::EberleinFlow_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;

    return ( G0inv - selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::EberleinFlow_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt  )
{
    // TODO: d_delta_mu_over_dt corrections not yet included
    double w = w_val(idx_w);
    double Lam = exp( t * LN_10 );
    MatQN GMat = EberleinFlow_G( idx_w, idx_p, t, selfEn, delta_mu );
    
    return GMat * (- I * sgn( w ) * Lam / sqrt( w*w + Lam*Lam ) *  LN_10 * Lam )*GMat;
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::EberleinFlow_asymptotic_SG_pp(const int idx_W, const double t)
{
    double Lam = exp( t * LN_10 );
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    return 0.0;  /**< TODO: find an analytic experession for the S-G bubble in high-frequency regime */

}

template <typename Model>
dcomplex fRGFlowScheme<Model>::EberleinFlow_asymptotic_GG_pp(const int idx_W, const double t)
{
    double Lam = exp( t * LN_10 );
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    return 0.0;  /**< TODO: find an analytic experession for the G-G bubble in high-frequency regime */
}


// Actual_UFlow
template <typename Model>
MatQN fRGFlowScheme<Model>::Nonflowing_G( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    return ( G0inv - selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::Nonflowing_G_latt( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;
    return ( G0inv - selfEn ).inverse();
}

template <typename Model>
MatQN fRGFlowScheme<Model>::Nonflowing_S( const int idx_w, const int idx_p, const double t, const MatQN& selfEn, const double delta_mu, const double d_delta_mu_over_dt )
{
    MatQN G0inv;
    double w = w_val(idx_w);
    G0inv << I*w - Model::E(idx_p, idx_w) + delta_mu;

    MatQN G = ( G0inv - selfEn ).inverse();
    MatQN d_delta_mu_over_dt_mat;
    d_delta_mu_over_dt_mat << d_delta_mu_over_dt;

    return G * (t*d_delta_mu_over_dt_mat) * G;
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::Nonflowing_asymptotic_SG_pp(const int idx_W, const double t)
{
    return 0;
}

template <typename Model>
dcomplex fRGFlowScheme<Model>::Nonflowing_asymptotic_GG_pp(const int idx_W, const double t)
{
    double PIR = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(idx_W/2);
    double W = idx_W;

    if( idx_W == 0 ) 
	return (BETA)/(2.0*pow(M_PI,2.0)*PIR);

    return (BETA*atanh(W/(2.0*PIR)))/(pow(M_PI,2.0)*W);
}

template <typename Model>
double fRGFlowScheme<Model>::U_TrivialCutoff(const double t)
{
    return 1.0;
}

template <typename Model>
double fRGFlowScheme<Model>::U_TrivialCutoff_dot(const double t)
{
    return 0.0;
}

template <typename Model>
double fRGFlowScheme<Model>::Actual_UFlow_MultiplicativeCutoff(const double t)
{
    const auto &flow = ChosenFlowParametrizationInfo();
    const double t_range = flow.t_final - flow.t_start;
    const double t_normalized = (t - flow.t_start) / t_range;
    return t_normalized;
}

template <typename Model>
double fRGFlowScheme<Model>::Actual_UFlow_MultiplicativeCutoff_dot(const double t)
{
    const auto &flow = ChosenFlowParametrizationInfo();
    const double t_range = flow.t_final - flow.t_start;
    const double dt_normalized_dt = 1.0 / t_range;
    return dt_normalized_dt;
}

// Actual_UFlow end
template <typename Model>
double fRGFlowScheme<Model>::LocalBareInteraction(const double t)
{
    if (G_ChosenFlowName() == G_FlowSchemeName::Interaction)
	return t*t*std::real(Model::vertex_4pt_local_part_bare(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
    else
	return std::real(Model::vertex_4pt_local_part_bare(0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0, 0));
}
    

// instantiate class for the particular models
#define Instantiate(MODEL) template class fRGFlowScheme<MODEL>;
WITH_MODELS(Instantiate)
