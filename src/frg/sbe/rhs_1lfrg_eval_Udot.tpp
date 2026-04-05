template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_w_sc_Udot_part( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    const auto bare_B_dot = state.bos_vertex_sc_bare_dot(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, t);

    dcomplex left_part( 0.0, 0.0 ), right_part( 0.0, 0.0);
    dcomplex F_part(0.0, 0.0);

    const int lambda_bos_freq_count = FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count;

#ifdef SET_MIXED_BUBBLES_TO_ZERO
    {
	int mpp_ = 0;
#else
    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif    
        for_freq(int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_){
 
	    left_part += state.lambda_sc( W, K, wp_, mpp_ ) *  
		bubble_GG_pp[W][wp_][mpp_][  0  ][0][0][0][0](K) ;
	    right_part += bubble_GG_pp[W][wp_][  0  ][mpp_][0][0][0][0](K) *
	        state.lambda_sc( W, K, wp_, mpp_ );
	}
	if ( W  <= lambda_bos_freq_count && -lambda_bos_freq_count <= W) 
	for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    F_part += m_lambda_sc_left_Fdot_part[W][K][wp_][mpp_] * bubble_GG_pp[W][wp_][  0  ][mpp_][0][0][0][0](K) * state.lambda_sc( W, K, wp_, mpp_ );
	}
#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */

	left_part +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_pp( W, t ) * 
	    state.lambda_sc( W, K, wp_, mpp_ );

	right_part +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_pp( W, t ) * 
	    state.lambda_sc( W, K, wp_, mpp_ );
#endif 

    }

    left_part *=  1.0 / inv_freq_sum_normalisation(t);
    right_part *=  1.0 / inv_freq_sum_normalisation(t);
    F_part *= 1.0 / inv_freq_sum_normalisation(t);

    const auto VolBZ = Model::MomentumGrid().ptr->get_volume();
    
    dcomplex w_sc = state.w_sc(W,K,t);
    val += bare_B_dot +
         w_sc*left_part*bare_B_dot/VolBZ + bare_B_dot*right_part*w_sc/VolBZ +
	 w_sc*left_part*bare_B_dot*right_part*w_sc/VolBZ/VolBZ+
	 w_sc*F_part*w_sc/VolBZ;
	; 



#ifdef F_NONZERO
    //    
#endif
    
    return val; 
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_w_d_Udot_part( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );

    const auto bare_B_dot = state.bos_vertex_d_bare_dot(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, t) ;

    dcomplex left_part( 0.0, 0.0 ), right_part( 0.0, 0.0);
    dcomplex F_part(0.0, 0.0);

    const int lambda_bos_freq_count = FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count;

#ifdef SET_MIXED_BUBBLES_TO_ZERO
    {
	int mpp_ = 0;
#else
    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif    
        for_freq(int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_){
 
	    left_part += state.lambda_d( W, K, wp_, mpp_ ) *  
		bubble_GG_ph[W][wp_][mpp_][  0  ][0][0][0][0](K) ;
	    right_part += bubble_GG_ph[W][wp_][  0  ][mpp_][0][0][0][0](K) * state.lambda_d( W, K, wp_, mpp_ );
	}
	if ( W  <= lambda_bos_freq_count && -lambda_bos_freq_count <= W)
	for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    F_part += m_lambda_d_left_Fdot_part[W][K][wp_][mpp_] * bubble_GG_ph[W][wp_][  0  ][mpp_][0][0][0][0](K) * state.lambda_d( W, K, wp_, mpp_ );
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */

	left_part +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) * 
	    state.lambda_d( W, K, wp_, mpp_ );

	right_part +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) * 
	    state.lambda_d( W, K, wp_, mpp_ );
#endif 

    }

    left_part *=  1.0 / inv_freq_sum_normalisation(t);
    right_part *=  1.0 / inv_freq_sum_normalisation(t);
    F_part *= 1.0 / inv_freq_sum_normalisation(t);

    const auto VolBZ = Model::MomentumGrid().ptr->get_volume();
    
    
    dcomplex w_d = state.w_d(W,K,t);
    
    val += -bare_B_dot +
      (w_d*left_part*bare_B_dot + bare_B_dot*right_part*w_d)/VolBZ -
           w_d*left_part*bare_B_dot*right_part*w_d/VolBZ/VolBZ -
	   w_d* F_part * w_d / VolBZ; 
    val *= -1.0;
   
    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_w_m_Udot_part( const idx_w_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( IW::W );
    int K   = idx( IW::K );
    dcomplex F_part(0.0, 0.0);

    const int lambda_bos_freq_count = FrequencyDependenceScheme<Model>::lambda_BosonicGrid().positive_freqs_count;

    const auto bare_B_dot = state.bos_vertex_m_bare_dot(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, t);

    
    dcomplex left_part( 0.0, 0.0 ), right_part( 0.0, 0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
    {
	int mpp_ = 0;
#else
    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif    
        for_freq(int wp_, -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2), wp_ < FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2), ++wp_){
 
	    left_part += state.lambda_m( W, K, wp_, mpp_ ) *  
		bubble_GG_ph[W][wp_][mpp_][  0  ][0][0][0][0](K) ;
	    right_part += bubble_GG_ph[W][wp_][  0  ][mpp_][0][0][0][0](K) * state.lambda_m( W, K, wp_, mpp_ );
	}
	if ( W  <= lambda_bos_freq_count && -lambda_bos_freq_count <= W)
	for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    F_part += m_lambda_m_left_Fdot_part[W][K][wp_][mpp_] * bubble_GG_ph[W][wp_][  0  ][mpp_][0][0][0][0](K) * state.lambda_m( W, K, wp_, mpp_ );
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);    /**< explicit large frequency contribution -- alternative to weight_vec */

	left_part +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) * 
	    state.lambda_m( W, K, wp_, mpp_ );

	right_part +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) * 
	    state.lambda_m( W, K, wp_, mpp_ );
#endif 

    }

    left_part *=  1.0 / inv_freq_sum_normalisation(t);
    right_part *=  1.0 / inv_freq_sum_normalisation(t);
    F_part  *=  1.0 / inv_freq_sum_normalisation(t);

    const auto VolBZ = Model::MomentumGrid().ptr->get_volume();
    
    
    dcomplex w_m = state.w_m(W,K,t);
    
    val += bare_B_dot  -
          (w_m*left_part*bare_B_dot + bare_B_dot*right_part*w_m)/VolBZ +
           w_m*left_part*bare_B_dot*right_part*w_m/VolBZ/VolBZ + 
	   w_m*F_part*w_m/VolBZ;

    
    return val; 
}

    

template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambda_sc_Udot_part( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    dcomplex F_part_2l( 0.0, 0.0 );
    dcomplex F_part_1l = m_lambda_sc_left_Fdot_part[W][K][w][m];

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
    
#ifdef SET_MIXED_BUBBLES_TO_ZERO
    int mp_ = m_;
    {
#else
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif    
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    F_part_2l += m_lambda_sc_left_Fdot_part[W][K][wp_][mp_] * bubble_GG_pp[W][wp_][  mp_  ][m_][0][0][0][0](K) * state.B_irreducible_vertex_sc( W, K, wp_, m_, w, m, t );
	}
    }
    }

    F_part_2l *= 1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    val += F_part_1l + F_part_2l;

    return val; 
}




template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambda_d_Udot_part( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    dcomplex F_part_2l( 0.0, 0.0 );
    dcomplex F_part_1l = m_lambda_d_left_Fdot_part[W][K][w][m];

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
    
#ifdef SET_MIXED_BUBBLES_TO_ZERO
    int mp_ = m_;
    {
#else
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif    
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    F_part_2l += m_lambda_d_left_Fdot_part[W][K][wp_][mp_] * bubble_GG_ph[W][wp_][  mp_  ][m_][0][0][0][0](K) * state.B_irreducible_vertex_d( W, K, wp_, m_, w, m, t );
	}

    }
    }

    F_part_2l *= 1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    val += - F_part_1l + F_part_2l;

    return val; 
}



template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambda_m_Udot_part( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

    dcomplex F_part_2l( 0.0, 0.0 );
    dcomplex F_part_1l = m_lambda_m_left_Fdot_part[W][K][w][m];

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
    
#ifdef SET_MIXED_BUBBLES_TO_ZERO
    int mp_ = m_;
    {
#else
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif    
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::lambda_FermionicGrid().positive_freqs_count, ++wp_ ){
	    F_part_2l += m_lambda_m_left_Fdot_part[W][K][wp_][mp_] * bubble_GG_ph[W][wp_][  mp_  ][m_][0][0][0][0](K) * state.B_irreducible_vertex_m( W, K, wp_, m_, w, m, t );
	}
    }
    }

    F_part_2l *= 1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    val += -F_part_1l + F_part_2l;

    return val; 
}




template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_M_sc_Udot_part( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( IM::W );
   int w_in   = idx( IM::w );
   int w_out  = idx( IM::wp);

   int K      = idx( IM::K );
   int n_in   = idx( IM::m );
   int n_out  = idx( IM::mp );

   dcomplex F_part_1l = m_M_sc_left_Fdot_part[W][K][w_in][n_in][w_out][n_out] + m_M_sc_left_Fdot_part[W][K][w_out][n_out][w_in][n_in];

   dcomplex F_part_2l( 0.0, 0.0 );

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
    
#ifdef SET_MIXED_BUBBLES_TO_ZERO
    int mp_ = m_;
    {
#else
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif    
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, ++wp_ ){
	    F_part_2l += m_M_sc_left_Fdot_part[W][K][w_in][n_in][wp_][mp_] * bubble_GG_pp[W][wp_][  mp_  ][m_][0][0][0][0](K) * state.B_irreducible_vertex_sc( W, K, wp_, m_, w_out, n_out, t );
	}

    }
    }

    F_part_2l *= 1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    val += F_part_1l + F_part_2l;

    return val; 
}



template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_M_d_Udot_part( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
   dcomplex val( 0.0, 0.0 );

   int W      = idx( IM::W );
   int w_in   = idx( IM::w );
   int w_out  = idx( IM::wp);

   int K      = idx( IM::K );
   int n_in   = idx( IM::m );
   int n_out  = idx( IM::mp );

   dcomplex F_part_1l = m_M_d_left_Fdot_part[W][K][w_in][n_in][w_out][n_out] + m_M_d_left_Fdot_part[W][K][w_out][n_out][w_in][n_in];

   dcomplex F_part_2l( 0.0, 0.0 );

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
    
#ifdef SET_MIXED_BUBBLES_TO_ZERO
    int mp_ = m_;
    {
#else
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif    
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, ++wp_ ){
	    F_part_2l += m_M_d_left_Fdot_part[W][K][w_in][n_in][wp_][mp_] * bubble_GG_ph[W][wp_][  mp_  ][m_][0][0][0][0](K) * state.B_irreducible_vertex_d( W, K, wp_, m_, w_out, n_out, t );
	}

    }
    }

    F_part_2l *= 1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    val += -F_part_1l + F_part_2l;

    return val; 
}
   

template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_M_m_Udot_part( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{ // todo: bosonize_M_lazy
   dcomplex val( 0.0, 0.0 );

   int W      = idx( IM::W );
   int w_in   = idx( IM::w );
   int w_out  = idx( IM::wp);

   int K      = idx( IM::K );
   int n_in   = idx( IM::m );
   int n_out  = idx( IM::mp );

   dcomplex F_part_1l = m_M_m_left_Fdot_part[W][K][w_in][n_in][w_out][n_out] + m_M_m_left_Fdot_part[W][K][w_out][n_out][w_in][n_in];

   dcomplex F_part_2l( 0.0, 0.0 );

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
    
#ifdef SET_MIXED_BUBBLES_TO_ZERO
    int mp_ = m_;
    {
#else
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif    
        for_freq( int wp_, -FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, wp_ < FrequencyDependenceScheme<Model>::M_FermionicGrid().positive_freqs_count, ++wp_ ){
	    F_part_2l += m_M_m_left_Fdot_part[W][K][w_in][n_in][wp_][mp_] * bubble_GG_ph[W][wp_][  mp_  ][m_][0][0][0][0](K) * state.B_irreducible_vertex_m( W, K, wp_, m_, w_out, n_out, t );
	}
    }
    }

    F_part_2l *= 1.0/ inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    val += -F_part_1l + F_part_2l;

    return val; 
}


// F_dot part only there when EXTENDED flag is on (for an extended model)
template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambda_sc_left_Fdot_part( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif

    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
	for_freq( int wp_, negative_frequency_integration_range, wp_ < positive_frequency_integration_range, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	    {
		int mpp_ = mp_;
#else
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
		temp += state.lambda_sc( W, K, wp_, mpp_ ) *
		    bubble_GG_pp[W][wp_][mpp_][mp_][0][0][0][0](K);
	    }
	    val += temp * state.ferm_vertex_sc_bare_dot(W, K, wp_, mp_, w, m, t);
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_pp( W, t ) *
	    state.lambda_sc( W, K, wp_, mp_ ) *
	    state.ferm_vertex_sc_bare_dot(W, K, wp_, mp_, w, m, t);
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambda_d_left_Fdot_part( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif

    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
	for_freq( int wp_, negative_frequency_integration_range, wp_ < positive_frequency_integration_range, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	    {
		int mpp_ = mp_;
#else
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
		temp += state.lambda_d( W, K, wp_, mpp_ ) *
		    bubble_GG_ph[W][wp_][mpp_][mp_][0][0][0][0](K);
	    }
	    val += temp * state.ferm_vertex_d_bare_dot(W, K, wp_, mp_, w, m, t);
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) *
	    state.lambda_d( W, K, wp_, mp_ ) *
	    state.ferm_vertex_d_bare_dot(W, K, wp_, mp_, w, m, t);
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}


template <typename Model, typename state_t > 
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_lambda_m_left_Fdot_part( const idx_lambda_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W   = idx( ILAMBDA::W );
    int w   = idx( ILAMBDA::w );

    int K   = idx( ILAMBDA::K );
    int m   = idx( ILAMBDA::m );

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif
    
    for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
	for_freq( int wp_, negative_frequency_integration_range, wp_ < positive_frequency_integration_range, ++wp_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	    {
		int mpp_ = mp_;
#else
	    for( int mpp_ = 0; mpp_ < Model::GetMomentumFormFactorsCount(); ++mpp_ ){
#endif
		temp += state.lambda_m( W, K, wp_, mpp_ ) *
		    bubble_GG_ph[W][wp_][mpp_][mp_][0][0][0][0](K);
	    }
	    val += temp * state.ferm_vertex_m_bare_dot(W, K, wp_, mp_, w, m, t);
	}

#ifndef STATIC_CALCULATION         
	int wp_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

	val +=  inv_freq_sum_normalisation(t) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) *
	    state.lambda_m( W, K, wp_, mp_ ) *
	     state.ferm_vertex_m_bare_dot(W, K, wp_, mp_, w, m, t);
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}




template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_M_sc_left_Fdot_part( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_pp, double t )
{
   dcomplex val( 0.0, 0.0 );
    
   int W      = idx( IM::W );
   int w_in   = idx( IM::w );
   int w_out  = idx( IM::wp);

   int K      = idx( IM::K );
   int n_in   = idx( IM::m );
   int n_out  = idx( IM::mp );

#ifdef SET_OFFDIAGONAL_M_IN_FROMFACTOR_SPACE_TO_ZERO
 if (n_in != n_out)
       return val;
#endif

#ifdef BOSONISE_M_LAZY
   if (n_in != n_out)
       return val;

   // shift w_out to capture a "good" v0 and -v0
   const int sign_of_w_out = (int(w_out >= 0) - int(w_out < 0));
   w_out +=  sign_of_w_out * abs(div2_floor(W)) - int(w_out < 0) *( 1 + ((W+100000)%2));
#endif

#ifdef SET_OFFDIAGONAL_M_IN_FROMFACTOR_SPACE_TO_ZERO
 if (n_in != n_out)
       return val;
#endif

#ifdef BOSONISE_M_LAZY
   if (n_in != n_out)
       return val;

   // shift w_out to capture a "good" v0 and -v0
   const int sign_of_w_out = (int(w_out >= 0) - int(w_out < 0));
   w_out +=  sign_of_w_out * abs(div2_floor(W)) - int(w_out < 0) *( 1 + ((W+100000)%2));
#endif


#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
   const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
   const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
   const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
   const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif
    
   for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
       for_freq( int w_, negative_frequency_integration_range, w_ < positive_frequency_integration_range, ++w_ ){
	   dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	{
	    int mp_ = m_;
#else
        for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
                temp += state.B_irreducible_vertex_sc( W, K, w_in, n_in, w_, mp_, t ) * bubble_GG_pp[W][w_][mp_][m_][0][0][0][0](K) ;   
            }
	    val += temp * state.ferm_vertex_sc_bare_dot(W, K, w_, m_, w_out, n_out, t);
        }

#ifndef STATIC_CALCULATION         
	int w_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

        val += inv_freq_sum_normalisation(t) * state.B_irreducible_vertex_sc( W, K, w_in, n_in, w_, m_, t ) * fRGFlowScheme<Model>::asymptotic_GG_pp( W, t ) * state.ferm_vertex_sc_bare_dot(W, K, w_, m_, w_out, n_out, t);
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

   return val;
}



template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_M_d_left_Fdot_part( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
    dcomplex val( 0.0, 0.0 );

    int W      = idx( IM::W );
    int w_in   = idx( IM::w );
    int w_out  = idx( IM::wp);

    int K      = idx( IM::K );
    int n_in   = idx( IM::m );
    int n_out  = idx( IM::mp );


#ifdef SET_OFFDIAGONAL_M_IN_FROMFACTOR_SPACE_TO_ZERO
 if (n_in != n_out)
       return val;
#endif


#ifdef BOSONISE_M_LAZY
   if (n_in != n_out)
       return val;

   // shift w_out to capture a "good" v0 and -v0
   const int sign_of_w_out = (int(w_out >= 0) - int(w_out < 0));
   w_out +=  sign_of_w_out * abs(div2_floor(W)) - int(w_out < 0) *( 1 + ((W+100000)%2));
#endif

#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif
    
    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
        for_freq( int w_, negative_frequency_integration_range, w_ < positive_frequency_integration_range, ++w_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
            for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
                temp += state.B_irreducible_vertex_d( W, K, w_in, n_in, w_, mp_, t ) * bubble_GG_ph[W][w_][mp_][m_][0][0][0][0](K) ;   
            }
	    val += temp * state.ferm_vertex_d_bare_dot(W, K, w_, m_, w_out, n_out, t);
        }

#ifndef STATIC_CALCULATION         
	int w_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

        val += inv_freq_sum_normalisation(t) * state.B_irreducible_vertex_d( W, K, w_in, n_in, w_, m_, t ) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t ) * state.ferm_vertex_d_bare_dot(W, K, w_, m_, w_out, n_out, t);
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}
   

template <typename Model, typename state_t >
dcomplex rhs_sbe_1lfrg_t<Model, state_t>::eval_M_m_left_Fdot_part( const idx_M_t<Model>& idx, const state_t& state, const gf_bubble_mat_t<Model>& bubble_GG_ph, double t )
{
   dcomplex val( 0.0, 0.0 );
    
   int W      = idx( IM::W );
   int w_in   = idx( IM::w );
   int w_out  = idx( IM::wp);

   int K      = idx( IM::K );
   int n_in   = idx( IM::m );
   int n_out  = idx( IM::mp );

   
#ifdef SET_OFFDIAGONAL_M_IN_FROMFACTOR_SPACE_TO_ZERO
 if (n_in != n_out)
       return val;
#endif

#ifdef BOSONISE_M_LAZY
   if (n_in != n_out)
       return val;

   // shift w_out to capture a "good" v0 and -v0
   const int sign_of_w_out = (int(w_out >= 0) - int(w_out < 0));
   w_out +=  sign_of_w_out * abs(div2_floor(W)) - int(w_out < 0) *( 1 + ((W+100000)%2));
#endif

    
#ifndef RESTRICT_FREQUENCY_SUMS_TO_STATE_PROJECTIONS_BOX
    const int negative_frequency_integration_range = -FrequencyDependenceScheme<Model>::PositiveIntegrationRange() - abs(W/2) - fold(W, 2);
    const int positive_frequency_integration_range = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);
#else //precomputed state projections are zero outside of a frequency box
    const int negative_frequency_integration_range = -state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count - 01;
    const int positive_frequency_integration_range = state.m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count + 01;
#endif

    for( int m_ = 0; m_ < Model::GetMomentumFormFactorsCount(); ++m_ ){
      for_freq( int w_, negative_frequency_integration_range, w_ < positive_frequency_integration_range, ++w_ ){
	    dcomplex temp = (0.0,0.0);
#ifdef SET_MIXED_BUBBLES_TO_ZERO
	   {
		int mp_ = m_;
#else
            for( int mp_ = 0; mp_ < Model::GetMomentumFormFactorsCount(); ++mp_ ){
#endif
                temp += state.B_irreducible_vertex_m( W, K, w_in, n_in, w_, mp_, t ) * bubble_GG_ph[W][w_][mp_][m_][0][0][0][0](K) ;   
            }
	    val += temp * state.ferm_vertex_m_bare_dot(W, K, w_, m_, w_out, n_out, t);
	}

#ifndef STATIC_CALCULATION         
	int w_ = FrequencyDependenceScheme<Model>::PositiveIntegrationRange() + abs(W/2);   //< explicit large frequency contribution -- alternative to weight_vec

        val += inv_freq_sum_normalisation(t) * state.B_irreducible_vertex_m( W, K, w_in, n_in, w_, m_, t ) * fRGFlowScheme<Model>::asymptotic_GG_ph( W, t )  * state.ferm_vertex_m_bare_dot(W, K, w_, m_, w_out, n_out, t);
#endif
    }

    val *= 1.0/inv_freq_sum_normalisation(t)/Model::MomentumGrid().ptr->get_volume();

    return val;
}

