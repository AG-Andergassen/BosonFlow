#pragma once

#include <frg/sbe/state.h>
#include <frg/sbe/state_gf_types_bosonised_M.h>

#include <boost/numeric/odeint.hpp>
#include <complex>
#include <iostream>


template <typename Model, typename... OtherTypes >
class state_frg_sbe_bosonised_M_t: public state_frg_sbe_t<Model,
  gf_w_M_t<Model>, gf_w_M_t<Model>, gf_w_M_t<Model>,  // screened interaction +
  gf_w_M_t<Model>, gf_w_M_t<Model>, gf_w_M_t<Model>,  // screened interaction -
  gf_lambda_M_t<Model>, gf_lambda_M_t<Model>, gf_lambda_M_t<Model>,  // hedin vertices + 
  gf_lambda_M_t<Model>, gf_lambda_M_t<Model>, gf_lambda_M_t<Model>,  // hedin vertices -
    OtherTypes...> 
  /**< Type of the state vector of the ODE solver. @see arithmetic_tuple.h for ReaK wrapper */
{
 public:
    using state_base_2_t = state_frg_sbe_t<Model,
  gf_w_M_t<Model>, gf_w_M_t<Model>, gf_w_M_t<Model>,  // screened interaction +
  gf_w_M_t<Model>, gf_w_M_t<Model>, gf_w_M_t<Model>,  // screened interaction -
  gf_lambda_M_t<Model>, gf_lambda_M_t<Model>, gf_lambda_M_t<Model>,  // hedin vertices + 
	gf_lambda_M_t<Model>, gf_lambda_M_t<Model>, gf_lambda_M_t<Model>, OtherTypes... >;  //hedin vertices -

    using w_M_t = gf_w_M_t<Model>; 
    using lambda_M_t = gf_lambda_M_t<Model>; 
    /**
     * \brief inline functions which return the \f$ n \f$ tuple element. Types specified above.
     *
     * The first inline function has not const such that it can be updated during the ODE integration
     * The second const inline function prevents any modification of the object inside the function
     */
         
    inline w_M_t& gf_wM_d_plus()			        { return std::get<state_base_2_t::offset_idx(0)>( *this ); }    
    inline const w_M_t& gf_wM_d_plus() const	        { return std::get<state_base_2_t::offset_idx(0)>( *this ); }    
 
    inline w_M_t& gf_wM_m_plus()			        { return std::get<state_base_2_t::offset_idx(1)>( *this ); }    
    inline const w_M_t& gf_wM_m_plus() const	        { return std::get<state_base_2_t::offset_idx(1)>( *this ); }    
 
    inline w_M_t& gf_wM_sc_plus()			        { return std::get<state_base_2_t::offset_idx(2)>( *this ); }    
    inline const w_M_t& gf_wM_sc_plus() const	        { return std::get<state_base_2_t::offset_idx(2)>( *this ); }    

    inline w_M_t& gf_wM_d_minus()			        { return std::get<state_base_2_t::offset_idx(3)>( *this ); }    
    inline const w_M_t& gf_wM_d_minus() const	        { return std::get<state_base_2_t::offset_idx(3)>( *this ); }    
 
    inline w_M_t& gf_wM_m_minus()			        { return std::get<state_base_2_t::offset_idx(4)>( *this ); }    
    inline const w_M_t& gf_wM_m_minus() const	        { return std::get<state_base_2_t::offset_idx(4)>( *this ); }    
 
    inline w_M_t& gf_wM_sc_minus()			        { return std::get<state_base_2_t::offset_idx(5)>( *this ); }    
    inline const w_M_t& gf_wM_sc_minus() const	        { return std::get<state_base_2_t::offset_idx(5)>( *this ); }    
  
    inline lambda_M_t& gf_lambdaM_d_plus()			{ return std::get<state_base_2_t::offset_idx(6)>( *this ); }    
    inline const lambda_M_t& gf_lambdaM_d_plus() const	{ return std::get<state_base_2_t::offset_idx(6)>( *this ); }    
 
    inline lambda_M_t& gf_lambdaM_m_plus()		        { return std::get<state_base_2_t::offset_idx(7)>( *this ); }    
    inline const lambda_M_t& gf_lambdaM_m_plus() const	{ return std::get<state_base_2_t::offset_idx(7)>( *this ); }    
 
    inline lambda_M_t& gf_lambdaM_sc_plus()			{ return std::get<state_base_2_t::offset_idx(8)>( *this ); }    
    inline const lambda_M_t& gf_lambdaM_sc_plus() const	{ return std::get<state_base_2_t::offset_idx(8)>( *this ); }    
    
    inline lambda_M_t& gf_lambdaM_d_minus()			{ return std::get<state_base_2_t::offset_idx(9)>( *this ); }    
    inline const lambda_M_t& gf_lambdaM_d_minus() const	{ return std::get<state_base_2_t::offset_idx(9)>( *this ); }    
 
    inline lambda_M_t& gf_lambdaM_m_minus()		        { return std::get<state_base_2_t::offset_idx(10)>( *this ); }    
    inline const lambda_M_t& gf_lambdaM_m_minus() const	{ return std::get<state_base_2_t::offset_idx(10)>( *this ); }    
 
    inline lambda_M_t& gf_lambdaM_sc_minus()			{ return std::get<state_base_2_t::offset_idx(11)>( *this ); }    
    inline const lambda_M_t& gf_lambdaM_sc_minus() const	{ return std::get<state_base_2_t::offset_idx(11)>( *this ); }    
    
               
 state_frg_sbe_bosonised_M_t():
    state_base_2_t()
        {}; 
 
 state_frg_sbe_bosonised_M_t( const state_frg_sbe_bosonised_M_t<Model, OtherTypes...>& state_vec ):
    state_base_2_t( state_vec )
    {}; 
 
 state_frg_sbe_bosonised_M_t( state_frg_sbe_bosonised_M_t<Model, OtherTypes...>&& state_vec ):
    state_base_2_t( std::move(state_vec) )
	{}; 
    
    state_frg_sbe_bosonised_M_t<Model, OtherTypes...>& operator=( const state_frg_sbe_bosonised_M_t<Model, OtherTypes...>& state_vec ){
	if( this == &state_vec)
	    return *this;
	state_base_2_t::operator=( state_vec );
	return *this;
    }
 
    state_frg_sbe_bosonised_M_t<Model, OtherTypes...>& operator=( state_frg_sbe_bosonised_M_t<Model, OtherTypes...>&& state_vec ){
	state_base_2_t::operator=( std::move( state_vec ) );
	return *this;
    }
    
    /**
     *  \brief Initialization of state_frg_sbe_t<Model> state-vector
     */
    
    static constexpr int offset_idx(const unsigned idx)
    {
	// 12: number of new components introduced by this state
	return state_base_2_t::offset_idx(idx + 12);
    }
 
    virtual void init_bare();         /*< initialization to the state bare quantities */

    virtual void init_zero(); /*< initialise this state to zero *>/

    /**
     *   \brief state_frg_sbe_t<Model> member functions. 
     *   
     *   \param set of indices which specify the specific GF
     *   \return the correspondent GF value for that index set
     */     
    
    dcomplex wM_sc_plus( const int W, const int K, const int m, const int mp) const;	/**< Return screened interaction w in the U-sc-reducible channel */
    dcomplex wM_m_plus( const int W, const int K, const int m, const int mp) const;	/**< Return screened interaction w in the U-m-reducible channel */
    dcomplex wM_d_plus( const int W, const int K, const int m, const int mp) const; 	/**< Return screened interaction w in the U-d-reducible channel */
    dcomplex wM_sc_minus( const int W, const int K, const int m, const int mp) const;	/**< Return screened interaction w in the U-sc-reducible channel */
    dcomplex wM_m_minus( const int W, const int K, const int m, const int mp) const;	/**< Return screened interaction w in the U-m-reducible channel */
    dcomplex wM_d_minus( const int W, const int K, const int m, const int mp) const; 	/**< Return screened interaction w in the U-d-reducible channel */

    // todo: fix comments
    dcomplex lambdaM_sc_plus( const int W, const int K, const int w, const int m, const int n) const; 	/**< Return Hedin vertex lambda in the U-sc_reducible channel */
    dcomplex lambdaM_m_plus( const int W, const int K, const int w, const int m, const int n) const; 	/**< Return Hedin vertex lambda in the U-m_reducible channel */
    dcomplex lambdaM_d_plus( const int W, const int K, const int w, const int m, const int n) const; 	/**< Return Hedin vertex lambda in the U-d_reducible channel */
    dcomplex lambdaM_sc_minus( const int W, const int K, const int w, const int m, const int n) const; 	/**< Return Hedin vertex lambda in the U-sc_reducible channel */
    dcomplex lambdaM_m_minus( const int W, const int K, const int w, const int m, const int n) const; 	/**< Return Hedin vertex lambda in the U-m_reducible channel */
    dcomplex lambdaM_d_minus( const int W, const int K, const int w, const int m, const int n) const; 	/**< Return Hedin vertex lambda in the U-d_reducible channel */


    virtual dcomplex M_sc( const int W, const int K, const int w, const int m, const int wp, const int mp) const;       /**< Return U-sc-irreducible but pp-reducilbe vertex */
    virtual dcomplex M_m( const int W, const int K, const int w, const int m, const int wp, const int mp) const;        /**< Return U-m-irreducible  but xph-reducible vertex */
    virtual dcomplex M_d( const int W, const int K, const int w, const int m, const int wp, const int mp) const;        /**< Return U-d-irreducible but ph-reducible vertex */
}; 


namespace boost {
    namespace numeric {
	namespace odeint {
	    template<typename Model>
		struct vector_space_norm_inf< state_frg_sbe_bosonised_M_t<Model> >
		{
		    typedef double result_type;
		    double operator()( const state_frg_sbe_bosonised_M_t<Model> &state_vec ) const
		    {
			using namespace std; 
			return norm( state_vec ); 
		    }
		};
	}
    }
}

// include implementation
#include <../src/frg/sbe/state_bosonised_M.tpp>
