#pragma once

#include <frg/sbe/state_gf_types.h>
#include <arithmetic_tuple.h>
#include <models/concrete_available_models.h>


template <typename Model, typename... OtherTypes >
class state_vector_frg_base_t : public ReaK::arithmetic_tuple< \
  gf_1p_t<Model>,  // self-energy
    gf_double_t, // delta_mu
  OtherTypes... >
{
 public:
    using base_t = ReaK::arithmetic_tuple< \
	gf_1p_t<Model>,  // self-energy
	gf_double_t,   // delta_mu
	OtherTypes... >;

    using Sig_t = gf_1p_t<Model>; 
    
 state_vector_frg_base_t():
    base_t()
        {}; 
    
 state_vector_frg_base_t( const state_vector_frg_base_t<Model, OtherTypes... >& state_vec ):
    base_t( state_vec )
    {}; 
    
 state_vector_frg_base_t( state_vector_frg_base_t<Model, OtherTypes... >&& state_vec ):
    base_t( std::move(state_vec) )
	{}; 
    
    state_vector_frg_base_t<Model, OtherTypes... >& operator=( const state_vector_frg_base_t<Model, OtherTypes... >& state_vec ){
	if( this == &state_vec)
	    return *this;
	base_t::operator=( state_vec );
	return *this;
    }
    
    state_vector_frg_base_t<Model, OtherTypes... >& operator=( state_vector_frg_base_t<Model, OtherTypes... >&& state_vec ){
	base_t::operator=( std::move( state_vec ) );
	return *this;
    }
    

    
    /**
     * \brief inline functions which return the \f$ n \f$ tuple element. Types specified above.
     *
     * The first inline function has not const such that it can be updated during the ODE integration
     * The second const inline function prevents any modification of the object inside the function
     */

    inline Sig_t& gf_Sig()			        { return std::get<0>( *this ); }    
    inline const Sig_t& gf_Sig() const	        { return std::get<0>( *this ); }  

    inline double& gf_delta_mu()			        { return std::get<1>( *this )[0]; }    
    inline const double& gf_delta_mu() const	        { return std::get<1>( *this )[0]; }  

    static constexpr int offset_idx(unsigned idx)
    {
	return idx + 2;
    }

    dcomplex Sig( int w, int k, int s_in, int s_out ) const; 		/**< Return self-energy */
 
    MatQN SigMat( int w, int k ) const; 				                /**< Return self-energy quantum number matrix for specific momentum and frequency given the current state vector */
    
    MatQN SigMat_big( const idx_1p_mat_t<Model>& idx ) const; 	    /**< Return self-energy quantum number matrix for specific momentum and frequency given the current state vector. Finer momentum resolution for the FT bubble calculation. */
    
    void adjust_chemical_potential_shift(const double desired_filling, const double t);

    void normalise_Sig_Temperature_Flow(const double t);

    bool adjust_Sig_VH();

    dcomplex calculate_RE_Sig_VH();

    void init_bare();          /*< initialization to bare quantities */

    void init_zero();          /*< initialization to zero */
    

    double eval_filling( const double t, const double delta_mu_) const;

    double eval_filling( const double t) const;

    double m_d_delta_mu_over_dt = 0.0; // used if filling is to be kept fixed. The derivative of the chemical potential shift with respect to the RG scale.

    bool m_initialised = false;
};


// implementation follows for template
#include "../src/frg/state_base.tpp"
