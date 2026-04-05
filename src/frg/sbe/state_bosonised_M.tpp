#include <mymath.h>
#include <tu_projections.h>

#include <iostream>

#define SQU(a) ((a)*(a))

template<typename Model, typename... OtherTypes > 
void state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::init_bare(){
    //initialise state with vertices being the bare couplings

    state_base_2_t::init_bare();

    gf_wM_sc_plus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_d_plus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_m_plus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_sc_minus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_d_minus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_m_minus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});

    gf_lambdaM_sc_plus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_d_plus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_m_plus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_sc_minus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_d_minus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_m_minus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
}

template<typename Model, typename... OtherTypes > 
void state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::init_zero(){
    //initialise state vector with zeroes

    state_base_2_t::init_zero();

    gf_wM_sc_plus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_d_plus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_m_plus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_sc_minus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_d_minus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});
    gf_wM_m_minus().init( [this](const idx_w_M_t<Model> &idx){ return 0.0;});

    gf_lambdaM_sc_plus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_d_plus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_m_plus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_sc_minus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_d_minus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
    gf_lambdaM_m_minus().init( [this](const idx_lambda_M_t<Model> &idx){return 0.0;}  );
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::M_sc( const int W, const int K, const int w, const int m, const int wp, const int mp) const
{
   dcomplex val(0.0,0.0);

   if( (unsigned)( W + this->m_frequency_ranges.M_bosonic_positive_freqs_count ) < (this->m_frequency_ranges.M_bosonic_freqs_count)   &&
       (unsigned)(w + this->m_frequency_ranges.M_fermionic_positive_freqs_count) < (this->m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2) &&
       (unsigned)(wp + this->m_frequency_ranges.M_fermionic_positive_freqs_count) < (this->m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2))
   {   
       auto wM_plus = wM_sc_plus(W, K, m, mp);
       if (std::abs(wM_plus) > 10e-7) 
	   val += (lambdaM_sc_plus(W, K, w, m, mp)/wM_plus/2.0) * wM_sc_plus(W, K, m, mp) * (lambdaM_sc_plus(W, K, wp, m, mp)/wM_plus/2.0);
       
       auto wM_minus = wM_sc_minus(W, K, m, mp);
       if (std::abs(wM_minus) > 10e-7)
	   val += (lambdaM_sc_minus(W, K, w, m, mp)/wM_minus/2.0) * wM_sc_minus(W, K, m, mp) * (lambdaM_sc_minus(W, K, wp, m, mp)/wM_minus/2.0);
   }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::M_d( const int W, const int K, const int w, const int m, const int wp, const int mp) const
{
   dcomplex val(0.0,0.0);
   if( (unsigned)( W + this->m_frequency_ranges.M_bosonic_positive_freqs_count ) < (this->m_frequency_ranges.M_bosonic_freqs_count)   &&
       (unsigned)(w + this->m_frequency_ranges.M_fermionic_positive_freqs_count) < (this->m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2) &&
       (unsigned)(wp + this->m_frequency_ranges.M_fermionic_positive_freqs_count) < (this->m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2))
   {   
       auto wM_plus = wM_d_plus(W, K, m, mp);
       if (std::abs(wM_plus) > 10e-7) 
	   val += (lambdaM_d_plus(W, K, w, m, mp)/wM_plus/2.0) * wM_d_plus(W, K, m, mp) * (lambdaM_d_plus(W, K, wp, m, mp)/wM_plus/2.0);
       
       auto wM_minus = wM_d_minus(W, K, m, mp);
       if (std::abs(wM_minus) > 10e-7)
	   val += (lambdaM_d_minus(W, K, w, m, mp)/wM_minus/2.0) * wM_d_minus(W, K, m, mp) * (lambdaM_d_minus(W, K, wp, m, mp)/wM_minus/2.0);
   }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::M_m( const int W, const int K, const int w, const int m, const int wp, const int mp) const
{
   dcomplex val(0.0,0.0);
   if( (unsigned)( W + this->m_frequency_ranges.M_bosonic_positive_freqs_count ) < (this->m_frequency_ranges.M_bosonic_freqs_count)   &&
       (unsigned)(w + this->m_frequency_ranges.M_fermionic_positive_freqs_count) < (this->m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2) &&
       (unsigned)(wp + this->m_frequency_ranges.M_fermionic_positive_freqs_count) < (this->m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2))
   {   
       auto wM_plus = wM_m_plus(W, K, m, mp);
       if (std::abs(wM_plus) > 10e-7) 
	   val += (lambdaM_m_plus(W, K, w, m, mp)/wM_plus/2.0) * wM_m_plus(W, K, m, mp) * (lambdaM_m_plus(W, K, wp, m, mp)/wM_plus/2.0);
       
       auto wM_minus = wM_m_minus(W, K, m, mp);
       if (std::abs(wM_minus) > 10e-7)
	   val += (lambdaM_m_minus(W, K, w, m, mp)/wM_minus/2.0) * wM_m_minus(W, K, m, mp) * (lambdaM_m_minus(W, K, wp, m, mp)/wM_minus/2.0);
   }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::wM_sc_plus( const int W, const int K, const int m, const int mp) const{
#ifdef STATIC_CALCULATION
    return gf_wM_sc_plus()[0][K][m][mp];
#endif

 if( (unsigned)(W + this->m_frequency_ranges.w_bosonic_positive_freqs_count) < this->m_frequency_ranges.w_bosonic_freqs_count )  // TODO: if( likely())
     return gf_wM_sc_plus()[W][K][m][mp];
 return 0.0;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::wM_d_plus( const int W, const int K, const int m, const int mp) const{
#ifdef STATIC_CALCULATION
    return gf_wM_d_plus()[0][K][m][mp];
#endif

 if( (unsigned)(W + this->m_frequency_ranges.w_bosonic_positive_freqs_count) < this->m_frequency_ranges.w_bosonic_freqs_count )
     return gf_wM_d_plus()[W][K][m][mp];
 return 0.0;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::wM_m_plus( const int W, const int K, const int m, const int mp) const{
#ifdef STATIC_CALCULATION
    return gf_wM_m_plus()[0][K][m][mp];
#endif

 if( (unsigned)(W + this->m_frequency_ranges.w_bosonic_positive_freqs_count) < this->m_frequency_ranges.w_bosonic_freqs_count )
     return gf_wM_m_plus()[W][K][m][mp];
 return 0.0;
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::wM_sc_minus( const int W, const int K, const int m, const int mp) const{
#ifdef STATIC_CALCULATION
    return gf_wM_sc_minus()[0][K][m][mp];
#endif

 if( (unsigned)(W + this->m_frequency_ranges.w_bosonic_positive_freqs_count) < this->m_frequency_ranges.w_bosonic_freqs_count )  // TODO: if( likely())
     return gf_wM_sc_minus()[W][K][m][mp];
 return 0.0;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::wM_d_minus( const int W, const int K, const int m, const int mp) const{
#ifdef STATIC_CALCULATION
    return gf_wM_d_minus()[0][K][m][mp];
#endif

 if( (unsigned)(W + this->m_frequency_ranges.w_bosonic_positive_freqs_count) < this->m_frequency_ranges.w_bosonic_freqs_count )
     return gf_wM_d_minus()[W][K][m][mp];
 return 0.0;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::wM_m_minus( const int W, const int K, const int m, const int mp) const{
#ifdef STATIC_CALCULATION
    return gf_wM_m_minus()[0][K][m][mp];
#endif

 if( (unsigned)(W + this->m_frequency_ranges.w_bosonic_positive_freqs_count) < this->m_frequency_ranges.w_bosonic_freqs_count )
     return gf_wM_m_minus()[W][K][m][mp];
 return 0.0;
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::lambdaM_sc_plus( const int W, const int K, const int w, const int m, const int n) const
{
#ifdef STATIC_CALCULATION
    return gf_lambdaM_sc_plus()[0][K][0][m][n];
#endif

   if( (unsigned)(W +this->m_frequency_ranges.lambda_bosonic_positive_freqs_count) < this->m_frequency_ranges.lambda_bosonic_freqs_count && 
     (unsigned)(w + this->m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (this->m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
       return gf_lambdaM_sc_plus()[W][K][w][m][n];
   return 0.0; 
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::lambdaM_d_plus( const int W, const int K, const int w, const int m, const int n) const
{
#ifdef STATIC_CALCULATION
    return gf_lambdaM_d_plus()[0][K][0][m][n];
#endif

   if( (unsigned)(W +this->m_frequency_ranges.lambda_bosonic_positive_freqs_count) < this->m_frequency_ranges.lambda_bosonic_freqs_count && 
     (unsigned)(w + this->m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (this->m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
       return gf_lambdaM_d_plus()[W][K][w][m][n];
   return 0.0; 
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::lambdaM_m_plus( const int W, const int K, const int w, const int m, const int n) const
{
#ifdef STATIC_CALCULATION
    return gf_lambdaM_m_plus()[0][K][0][m][n];
#endif

   if( (unsigned)(W +this->m_frequency_ranges.lambda_bosonic_positive_freqs_count) < this->m_frequency_ranges.lambda_bosonic_freqs_count && 
     (unsigned)(w + this->m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (this->m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
   {
       return gf_lambdaM_m_plus()[W][K][w][m][n];
   }
   return 0.0;
}   


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::lambdaM_sc_minus( const int W, const int K, const int w, const int m, const int n) const
{
#ifdef STATIC_CALCULATION
    return gf_lambdaM_sc_minus()[0][K][0][m][n];
#endif

   if( (unsigned)(W +this->m_frequency_ranges.lambda_bosonic_positive_freqs_count) < this->m_frequency_ranges.lambda_bosonic_freqs_count && 
     (unsigned)(w + this->m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (this->m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
       return gf_lambdaM_sc_minus()[W][K][w][m][n];
   return 0.0; 
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::lambdaM_d_minus( const int W, const int K, const int w, const int m, const int n) const
{
#ifdef STATIC_CALCULATION
    return gf_lambdaM_d_minus()[0][K][0][m][n];
#endif

   if( (unsigned)(W +this->m_frequency_ranges.lambda_bosonic_positive_freqs_count) < this->m_frequency_ranges.lambda_bosonic_freqs_count && 
     (unsigned)(w + this->m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (this->m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
       return gf_lambdaM_d_minus()[W][K][w][m][n];
   return 0.0; 
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_bosonised_M_t<Model, OtherTypes...>::lambdaM_m_minus( const int W, const int K, const int w, const int m, const int n) const
{
#ifdef STATIC_CALCULATION
    return gf_lambdaM_m_minus()[0][K][0][m][n];
#endif

   if( (unsigned)(W +this->m_frequency_ranges.lambda_bosonic_positive_freqs_count) < this->m_frequency_ranges.lambda_bosonic_freqs_count && 
     (unsigned)(w + this->m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (this->m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
   {
       return gf_lambdaM_m_minus()[W][K][w][m][n];
   }
   return 0.0;
}   
