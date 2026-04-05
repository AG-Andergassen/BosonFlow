// Implementation for template
#include <frg/flows.h>
#include <frg/symmetries_common.h>
#include <boost/math/tools/roots.hpp>
#include <frequencies/matsubara_space.h>

// --- SELF ENERGY ---
template <typename Model, typename... OtherTypes >
dcomplex state_vector_frg_base_t<Model, OtherTypes... >::Sig( int w, int k, int s_in, int s_out ) const{
#ifdef STATIC_CALCULATION         
	return 0.0; 
#endif

    if ( w < -FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count || w > FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count - 1 ) 
	return 0.0; 
    
    return gf_Sig()[w][k][s_in][s_out]; 
}


template <typename Model, typename... OtherTypes >
MatQN state_vector_frg_base_t<Model, OtherTypes... >::SigMat( int w, int k ) const{
    // todo: for static_calculation with self-energy, consider something else..
#ifdef STATIC_CALCULATION         
	return MatQN::Zero(); 
#endif

    if ( w < -FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count || w > FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count - 1 ) 
	return MatQN::Zero(); 
    //todo: implement proper spin dependence
    return Eigen::Map<const MatQN>( &(gf_Sig()[w][k][0][0]) );  
}


template <typename Model, typename... OtherTypes >
MatQN state_vector_frg_base_t<Model, OtherTypes... >::SigMat_big( const idx_1p_mat_t<Model>& idx ) const{
    // todo: for static_calculation with self-energy, consider something else..
#ifdef STATIC_CALCULATION         
	return MatQN::Zero(); 
#endif

    int w = idx(I1PMAT::w); 
    int idx_k = Model::FineToCoarseIdxMap()[ idx(I1PMAT::p) ];
    if ( w < -FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count || w > FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count - 1 ) 
	return MatQN::Zero(); 
    
    return Eigen::Map<const MatQN>( &(gf_Sig()[w][idx_k][0][0]) );  
}


template <typename Model, typename... OtherTypes >
void state_vector_frg_base_t<Model, OtherTypes... >::init_bare()
{
    m_initialised = true;
    gf_delta_mu() = 0.0;

    gf_Sig().init(  [](const idx_1p_t<Model> &idx){return 0.0;}  );

}

template <typename Model, typename... OtherTypes >
void state_vector_frg_base_t<Model, OtherTypes... >::init_zero()
{
    m_initialised = true;  

    gf_Sig().init(  [](const idx_1p_t<Model> &idx){return 0.0;}  );
    gf_delta_mu() = 0.0;
}


template <typename Model, typename... OtherTypes >
void state_vector_frg_base_t<Model, OtherTypes... >::normalise_Sig_Temperature_Flow(const double t)
{
    gf_1p_t<Model> Sig_normalized;
    SymmetriesfRGCommon<Model>::IdxEquivClasses_sig_Ptr()->init( Sig_normalized, [this,t]( const idx_1p_t<Model>& idx ){ 
	    double Sqrt_T = exp( t * LN_10 * 0.5 );
	    return Sqrt_T * this->Sig( idx(I1P::w), idx(I1P::k), idx(I1P::s_in), idx(I1P::s_out));
	} );
    gf_Sig() = Sig_normalized;
}


template <typename Model, typename... OtherTypes >
dcomplex state_vector_frg_base_t<Model, OtherTypes... >::calculate_RE_Sig_VH()
{
    int W = 0;
    int K = -1;
    for (unsigned i = 0; i < Model::SpecialPointsAndPaths().points_names.size(); i ++){
	if (Model::SpecialPointsAndPaths().points_names[i] == "idx_VanHove")
	    K = Model::SpecialPointsAndPaths().points[i];
    }

    if (K == -1){
	K = 0;
	std::cout << "Warning: flag VAN_HOVE is on, but model does not define a VanHove idx in its list of special points.. Will pin the fermi-surface at idx-0 instead." << std::endl;
    }

    int S_in = 0;
    int S_out = 0;
    dcomplex RE_Sig_VH = this->Sig(W,K,S_in,S_out).real();
    return RE_Sig_VH;
}


template <typename Model, typename... OtherTypes >
bool state_vector_frg_base_t<Model, OtherTypes... >::adjust_Sig_VH()
{
    dcomplex RE_Sig_VH = this->calculate_RE_Sig_VH();
    gf_1p_t<Model> Sig_corr;
    SymmetriesfRGCommon<Model>::IdxEquivClasses_sig_Ptr()->init( Sig_corr, [this,RE_Sig_VH]( const idx_1p_t<Model>& idx ){ 
	    return this->Sig(idx( I1P::w ),idx(I1P::k),idx( I1P::s_in ),idx( I1P::s_out ))-RE_Sig_VH;; 
	} );	
    this->gf_Sig() = Sig_corr;
    return 1;
}


template <typename Model, typename... OtherTypes >
double state_vector_frg_base_t<Model, OtherTypes... >::eval_filling(const double t, const double delta_mu) const
{
   dcomplex val( 0.0, 0.0 );
   
   for( int w = -FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count; w < FrequencyDependenceScheme<Model>::G_FermionicGrid().positive_freqs_count; ++w ){
   for( int p = 0; p < Model::GetFineMomentaCount(); ++p )
   for(int s_in =0; s_in < Model::GetQuantumNumbersCount(); ++s_in)
   for(int s_out =0; s_out < Model::GetQuantumNumbersCount(); ++s_out){
       int k = Model::GetCoarseMomentumIdxFromFine(p);
       val += fRGFlowScheme<Model>::G_latt(w, p, t, SigMat(w, k), delta_mu)(s_out, s_in)/Beta(t)/double(Model::GetFineMomentaCount());
   }
   }
   
   if( std::abs( std::imag( val )) > 0.00001 )
       std::cout << " CAUTION, nonvanishing imaginary part of filling: " << std::imag( val ) << std::endl; 
   double filling = std::real( val ) + 0.5;
   return filling; 
}


template <typename Model, typename... OtherTypes >
double state_vector_frg_base_t<Model, OtherTypes... >::eval_filling(const double t) const
{
    return this->eval_filling(t, this->gf_delta_mu());
}


template <typename T>
int newton_find_root(std::function<T(T)> func, T &x_guess, T tolerance, int maxIterations) {
    T prevX = x_guess;
    T f_of_prevX = func(x_guess);
    if (std::abs(f_of_prevX) < tolerance)
	return 0;

    T prev_prevX = x_guess - 1e-2;
    T f_of_prev_prevX = func(prev_prevX);
    
    for (int i = 0; i < maxIterations; ++i) {
	x_guess = prevX - f_of_prevX * (prevX - prev_prevX)/(f_of_prevX - f_of_prev_prevX);

	prev_prevX = prevX;
	f_of_prev_prevX = f_of_prevX;
	
	prevX = x_guess;
	f_of_prevX = func(prevX);

	if (std::abs(f_of_prevX) < tolerance)
	    return i;
    }

    return maxIterations + 1;
}

template <typename Model, typename... OtherTypes >
void state_vector_frg_base_t<Model, OtherTypes... >::adjust_chemical_potential_shift(const double desired_filling, const double t)
{
#ifndef STATIC_CALCULATION 
    double sig_high_frequency_asymptote = std::real(gf_Sig()[-FrequencyDependenceScheme<Model>::Sig_FermionicGrid().positive_freqs_count][0][0][0]);
    
    // we shift away the high frequency asymptote
    SymmetriesfRGCommon<Model>::IdxEquivClasses_sig_Ptr()->init( this->gf_Sig(), [this, sig_high_frequency_asymptote]( const idx_1p_t<Model>& idx ){ 
	return this->Sig(idx( I1P::w ),idx(I1P::k),idx( I1P::s_in ),idx( I1P::s_out ))-sig_high_frequency_asymptote; 
      } );	
#endif

#ifdef PIN_SELFENERGY_TAIL_TO_ZERO
    // we absorb it into the chemical potential shift
    gf_delta_mu() -= sig_high_frequency_asymptote;
    std::cout << "... Shifted high frequency asymptote of Sig by: " << sig_high_frequency_asymptote << std::endl;
#endif
    auto n_of_delta_mu = [t, desired_filling, this](double delta_mu_) { 
      return this->eval_filling(t, delta_mu_) - desired_filling; 
    };

    // Set up the fixed-point iteration solver
    int maxIterations = 50;
    double tolerance = 1e-5; // require agreement upto 5 decimal places
    double delta_mu = gf_delta_mu();
    int iterations;
    
    for (int i = 0; i < 10; i++){
      iterations = newton_find_root<double>(n_of_delta_mu, delta_mu, tolerance, maxIterations);
      if (iterations >= maxIterations){
	delta_mu = gf_delta_mu() + std::rand() / RAND_MAX * 0.02 - 0.01;
      }else{
	gf_delta_mu() = delta_mu;
	break;
      }
    }

    if (iterations >= maxIterations){
	std::cout << "Warning!: fixed point iteration for the chemical potential did not converge!";
	std::cerr << "Warning!: fixed point iteration for the chemical potential did not converge!";
    }

    std::cout << "... delta_mu: " << gf_delta_mu() << " (" << iterations <<  " iterations)"<< std::endl;
}





