#pragma once

#include <frg/state_base.h>
// todo: check why/whether to include the relevant gf_types header
#include <boost/numeric/odeint.hpp>
#include <complex>
#include <iostream>

template <typename gf_type>
struct projected_gf_t
{
    bool are_initialised = false;
    bool are_calculated = false;

    gf_type *sc_pp_to_ph_ptr;
    gf_type *sc_pp_to_xph_ptr;
    gf_type *d_ph_to_xph_ptr;
    gf_type *d_ph_to_pp_ptr;
    gf_type *m_ph_to_pp_ptr;
    gf_type *m_ph_to_xph_ptr;
    gf_type *m_xph_to_pp_ptr;
    gf_type *m_xph_to_ph_ptr;

    projected_gf_t<gf_type>& copy_if_initialised( const projected_gf_t<gf_type>& projected_gf_container ){
	if( this == &projected_gf_container)
	    return *this;
	if (projected_gf_container.are_initialised & this->are_initialised){
	    *sc_pp_to_ph_ptr = *projected_gf_container.sc_pp_to_ph_ptr;
	    *sc_pp_to_xph_ptr = *projected_gf_container.sc_pp_to_xph_ptr;
	    *d_ph_to_xph_ptr = *projected_gf_container.d_ph_to_xph_ptr;
	    *d_ph_to_pp_ptr = *projected_gf_container.d_ph_to_pp_ptr;
	    *m_ph_to_pp_ptr = *projected_gf_container.m_ph_to_pp_ptr;
	    *m_ph_to_xph_ptr = *projected_gf_container.m_ph_to_xph_ptr;
	    *m_xph_to_pp_ptr = *projected_gf_container.m_xph_to_pp_ptr;
	    *m_xph_to_ph_ptr = *projected_gf_container.m_xph_to_ph_ptr;
	}
	return *this;
    }
 
    
};

template <typename Model, typename... OtherTypes >
class state_frg_sbe_t: public state_vector_frg_base_t<Model,
  gf_w_t<Model>, gf_w_t<Model>, gf_w_t<Model>,  // screened interaction
  gf_lambda_t<Model>, gf_lambda_t<Model>, gf_lambda_t<Model>,  // hedin vertices
  gf_M_t<Model>, gf_M_t<Model>, gf_M_t<Model>, // related to rest functions
    gf_f_t<Model>,
    OtherTypes...> 
  /**< Type of the state vector of the ODE solver. @see arithmetic_tuple.h for ReaK wrapper */
{
 public:
    using state_base_t = state_vector_frg_base_t<
	Model, gf_w_t<Model>, gf_w_t<Model>, gf_w_t<Model>,  // screened interactioon
	gf_lambda_t<Model>, gf_lambda_t<Model>, gf_lambda_t<Model>,  // hedin vertices
        gf_M_t<Model>, gf_M_t<Model>, gf_M_t<Model>, // rest funcs related
	gf_f_t<Model>, OtherTypes... >;  // the free-energy

    using base_t = typename state_base_t::base_t;

    using w_t = gf_w_t<Model>; 
    using lambda_t = gf_lambda_t<Model>; 
    using M_t = gf_M_t<Model>;
    using f_t = gf_f_t<Model>;
    /**
     * \brief inline functions which return the \f$ n \f$ tuple element. Types specified above.
     *
     * The first inline function has not const such that it can be updated during the ODE integration
     * The second const inline function prevents any modification of the object inside the function
     */
         
    inline w_t& gf_w_d()			        { return std::get<state_base_t::offset_idx(0)>( *this ); }    
    inline const w_t& gf_w_d() const	        { return std::get<state_base_t::offset_idx(0)>( *this ); }    
 
    inline w_t& gf_w_m()			        { return std::get<state_base_t::offset_idx(1)>( *this ); }    
    inline const w_t& gf_w_m() const	        { return std::get<state_base_t::offset_idx(1)>( *this ); }    
 
    inline w_t& gf_w_sc()			        { return std::get<state_base_t::offset_idx(2)>( *this ); }    
    inline const w_t& gf_w_sc() const	        { return std::get<state_base_t::offset_idx(2)>( *this ); }    
 
    inline lambda_t& gf_lambda_d()			{ return std::get<state_base_t::offset_idx(3)>( *this ); }    
    inline const lambda_t& gf_lambda_d() const	{ return std::get<state_base_t::offset_idx(3)>( *this ); }    
 
    inline lambda_t& gf_lambda_m()		        { return std::get<state_base_t::offset_idx(4)>( *this ); }    
    inline const lambda_t& gf_lambda_m() const	{ return std::get<state_base_t::offset_idx(4)>( *this ); }    
 
    inline lambda_t& gf_lambda_sc()			{ return std::get<state_base_t::offset_idx(5)>( *this ); }    
    inline const lambda_t& gf_lambda_sc() const	{ return std::get<state_base_t::offset_idx(5)>( *this ); }    
 
    inline M_t& gf_M_d()                            { return std::get<state_base_t::offset_idx(6)>( *this ); }
    inline const M_t& gf_M_d() const                { return std::get<state_base_t::offset_idx(6)>( *this ); }

    inline M_t& gf_M_m()                            { return std::get<state_base_t::offset_idx(7)>( *this ); }
    inline const M_t& gf_M_m() const                { return std::get<state_base_t::offset_idx(7)>( *this ); }

    inline M_t& gf_M_sc()                           { return std::get<state_base_t::offset_idx(8)>( *this ); }
    inline const M_t& gf_M_sc() const               { return std::get<state_base_t::offset_idx(8)>( *this ); }

    inline f_t& gf_f()                           { return std::get<state_base_t::offset_idx(9)>( *this ); }
    inline const f_t& gf_f() const               { return std::get<state_base_t::offset_idx(9)>( *this ); }


    projected_gf_t<gf_nabla_t<Model> > m_projected_nablas;
    projected_gf_t<gf_M_t<Model> > m_projected_Ms;

    bool m_force_w_dot_asymptotics_to_zero = false;
  
    const frequency_ranges_sbe_t m_frequency_ranges;
  
 state_frg_sbe_t():
    state_base_t(), m_frequency_ranges(FrequencyDependenceScheme<Model>::GetFrequencyRanges())
        {}; 
 
 state_frg_sbe_t( const state_frg_sbe_t<Model, OtherTypes...>& state_vec ):
    state_base_t( state_vec ), m_frequency_ranges(FrequencyDependenceScheme<Model>::GetFrequencyRanges())
    {}; 
 
 state_frg_sbe_t( state_frg_sbe_t<Model, OtherTypes...>&& state_vec ):
    state_base_t( std::move(state_vec) ), m_frequency_ranges(FrequencyDependenceScheme<Model>::GetFrequencyRanges())
	{}; 
    
    state_frg_sbe_t<Model, OtherTypes...>& operator=( const state_frg_sbe_t<Model, OtherTypes...>& state_vec ){
	if( this == &state_vec)
	    return *this;

	//m_projected_nablas = state_vec.m_projected_nablas;
	//m_projected_Ms = state_vec.m_projected_nablas;
	m_force_w_dot_asymptotics_to_zero = state_vec.m_force_w_dot_asymptotics_to_zero;
	
	state_base_t::operator=( state_vec );

        #ifdef PRECOMPUTE_STATE_PROJECTIONS	
	    initialise_projected_Ms_and_nablas();
	    m_projected_nablas.copy_if_initialised(state_vec.m_projected_nablas);
	    m_projected_Ms.copy_if_initialised(state_vec.m_projected_Ms);
        #endif
	return *this;
    }
 
    state_frg_sbe_t<Model, OtherTypes...>& operator=( state_frg_sbe_t<Model, OtherTypes...>&& state_vec ){

      state_base_t::operator=( std::move( state_vec ) );
	m_force_w_dot_asymptotics_to_zero = state_vec.m_force_w_dot_asymptotics_to_zero;

	return *this;
    }

    ~state_frg_sbe_t();
  
  
    /**
     *  \brief Initialization of state_frg_sbe_t<Model> state-vector
     */
    
    static constexpr int offset_idx(const unsigned idx)
    {
	// 10: number of new components introduced by this state
	return state_base_t::offset_idx(idx + 10);
    }

    using Uirreducible_vertex_t = gf_M_t<Model>;     

    static Uirreducible_vertex_t &InterpolatingFlow_Uirreducible_d()
    {
	static Uirreducible_vertex_t interpolating_flow_Uirreducible_d(InputUirreducibleVertexPositiveBosFreqCount-1,
								       InputUirreducibleVertexPositiveFermFreqCount-1,
								       InputUirreducibleVertexPositiveFermFreqCount-1, 1, 1, 1);
	return interpolating_flow_Uirreducible_d;
    };

    static Uirreducible_vertex_t &InterpolatingFlow_Uirreducible_sc()
    {
	static Uirreducible_vertex_t interpolating_flow_Uirreducible_sc(InputUirreducibleVertexPositiveBosFreqCount-1,
									InputUirreducibleVertexPositiveFermFreqCount-1,
									InputUirreducibleVertexPositiveFermFreqCount-1, 1, 1, 1); // todo: make more flexible for possibility of momentum dependent starting point
	return interpolating_flow_Uirreducible_sc;
    };

    static Uirreducible_vertex_t &InterpolatingFlow_Uirreducible_m()
    {
	static Uirreducible_vertex_t interpolating_flow_Uirreducible_m(InputUirreducibleVertexPositiveBosFreqCount-1,
								       InputUirreducibleVertexPositiveFermFreqCount-1,
								       InputUirreducibleVertexPositiveFermFreqCount-1, 1, 1, 1);
	return interpolating_flow_Uirreducible_m;
    };

    static lambda_t &InterpolatingFlow_lambda_d()
    {
	static lambda_t interpolating_flow_lambda_d(Input_lambda_PositiveBosFreqCount-1,
								       Input_lambda_PositiveFermFreqCount-1);
	return interpolating_flow_lambda_d;
    };

    static lambda_t &InterpolatingFlow_lambda_sc()
    {
	static lambda_t interpolating_flow_lambda_sc(Input_lambda_PositiveBosFreqCount-1,
								       Input_lambda_PositiveFermFreqCount-1);
	return interpolating_flow_lambda_sc;
    };

    static lambda_t &InterpolatingFlow_lambda_m()
    {
	static lambda_t interpolating_flow_lambda_m(Input_lambda_PositiveBosFreqCount-1,
								       Input_lambda_PositiveFermFreqCount-1);
	return interpolating_flow_lambda_m;
    };
 
    virtual void init_bare();         /*< initialization to the state bare quantities */

    virtual void init_zero(); /*< initialise this state to zero >*/

    void initialise_projected_Ms_and_nablas();

    void precalculate_projected_Ms_and_nablas(const double t);

    void precalculate_projected_Ms_and_nablas_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const double t);


    /**
     *   \brief state_frg_sbe_t<Model> member functions. 
     *   
     *   \param set of indices which specify the specific GF
     *   \return the correspondent GF value for that index set
     */ 

    dcomplex vert_bare( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out) const;

    dcomplex vertex_DC( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out, const double t) const;

  
    dcomplex vertex_sc_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex vertex_d_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex vertex_m_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex vertex_DC_sc( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex vertex_DC_d( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex vertex_DC_m( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

  dcomplex vertex_DC_sc_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex vertex_DC_d_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex vertex_DC_m_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

  
    dcomplex local_vertex_4pt_bare( const int W, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex local_vertex_4pt_bare_dot( const int W, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex bos_vertex_sc_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex bos_vertex_d_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex bos_vertex_m_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex bos_vertex_sc_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex bos_vertex_d_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex bos_vertex_m_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;


    dcomplex ferm_vertex_sc_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex ferm_vertex_d_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex ferm_vertex_m_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex ferm_vertex_sc_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex ferm_vertex_d_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex ferm_vertex_m_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    
    dcomplex lambda_sc_bare( const int W, const int K, const int w, const int m ) const;
    
    dcomplex lambda_d_bare( const int W, const int K, const int w, const int m ) const;

    dcomplex lambda_m_bare( const int W, const int K, const int w, const int m ) const;
          

    dcomplex nabla_sc( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const; 	/**< Return the sc single boson exchanges */
    dcomplex nabla_m( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const; 	/**< Return the m single boson exchanges  */
    dcomplex nabla_d( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const; 	/**< Return the d single boson exchanges  */

    dcomplex vertex_4pt( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out, const double t) const; /**< Return the one-particle reducible vertex in the fermionic notation. Used in the conventional flow of the self-energy @\ref rhs::eval_Sig_conv() */
    

    dcomplex vertex_sc( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const; 	/**< Return vertex in sc channel parametrisation */
    dcomplex vertex_m( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const; 	/**< Return vertex in m channel parametrisation  */
    dcomplex vertex_d( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const; 	/**< Return vertex in d channel parametrisation  */
    
    
    dcomplex mom_lambda( dcomplex (state_frg_sbe_t<Model, OtherTypes...>::*lambda_func)( const int, const int, const int, const int) const, const int W, const int K, const int w, const int p) const;  /**< Return the Hedin vertex lambda depending where the form factor dependence is expanded on momenta. The lambda_func defines the channel */
    
    dcomplex mom_lambda_sc(const int W, const int K, const int w, const int p) const;  /**< Return the sc Hedin vertex lambda depending where the form factor dependence is expanded on momenta. */
    
    dcomplex mom_lambda_d(const int W, const int K, const int w, const int p) const;  /**< Return the d Hedin vertex lambda depending where the form factor dependence is expanded on momenta. */
    
    dcomplex mom_lambda_m(const int W, const int K, const int w, const int p) const;  /**< Return the sc Hedin vertex lambda depending where the form factor dependence is expanded on momenta. */

    dcomplex mom_lambda_sc_dot(const int W, const int K, const int w, const int p) const;  /**< Return the scale derivative of the sc Hedin vertex lambda depending where the form factor dependence is expanded on momenta. */
    
    dcomplex mom_lambda_d_dot(const int W, const int K, const int w, const int p) const;  /**< Return the scale derivative of the d Hedin vertex lambda depending where the form factor dependence is expanded on momenta. */
    
    dcomplex mom_lambda_m_dot(const int W, const int K, const int w, const int p) const;  /**< Return the scale derivative of the sc Hedin vertex lambda depending where the form factor dependence is expanded on momenta. */
   

       
    dcomplex mom_nabla_sc( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out, const double t ) const; /**< Return Hedin vertex in the superconducting channel in the purely fermionic notation. Used in @\ref state::vertex() */
    dcomplex mom_nabla_d( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out, const double t ) const;  /**< Return Hedin vertex in the density channel in the purely fermionic notation. Used in @\ref state::vertex() */
    dcomplex mom_nabla_m_ph( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out, const double t ) const; /**< Return Hedin vertex in the magnetic channel (in ph notation) in the purely fermionic notation. Used in @\ref state::vertex() */
    dcomplex mom_nabla_m( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out, const double t ) const; /**< Return Hedin vertex in the magnetic channel in the purely fermionic notation. Used in @\ref state::vertex() */

#if !defined(SBEa_APPROXIMATION)
    dcomplex mom_M_sc( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out ) const; /**< Return rest function in the superconducting channel in the purely fermionic notation. Used in @\ref state::vertex() */
    dcomplex mom_M_d( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out ) const;  /**< Return rest function in the denisty channel in the purely fermionic notation. Used in @\ref state::vertex() */
    dcomplex mom_M_m_ph( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out ) const; /**< Return rest function in the magnetic channel (in ph notation) in the purely fermionic notation. Used in @\ref state::vertex() */
    dcomplex mom_M_m( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out ) const;  /**< Return rest function in the magnetic channel in the purely fermionic notation. Used in @\ref state::vertex() */
# endif // !defined(SBEa_APPROXIMATION)
 
    dcomplex w_sc( const int W, const int K, const double t) const;	/**< Return screened interaction w in the U-sc-reducible channel */
    dcomplex w_m( const int W, const int K, const double t) const;	/**< Return screened interaction w in the U-m-reducible channel */
    dcomplex w_d( const int W, const int K, const double t) const; 	/**< Return screened interaction w in the U-d-reducible channel */
 
    dcomplex lambda_sc( const int W, const int K, const int w, const int m) const; 	/**< Return Hedin vertex lambda in the U-sc_reducible channel */
    dcomplex lambda_m( const int W, const int K, const int w, const int m) const; 	/**< Return Hedin vertex lambda in the U-m_reducible channel */
    dcomplex lambda_d( const int W, const int K, const int w, const int m) const; 	/**< Return Hedin vertex lambda in the U-d_reducible channel */

    dcomplex w_sc_dot( const int W, const int K, const double t) const;     /**< Return screened interaction w in the U-sc-reducible channel */
    dcomplex w_m_dot( const int W, const int K, const double t) const;      /**< Return screened interaction w in the U-m-reducible channel */
    dcomplex w_d_dot( const int W, const int K, const double t) const;      /**< Return screened interaction w in the U-d-reducible channel */

    dcomplex lambda_sc_dot( const int W, const int K, const int w, const int m) const;      /**< Return Hedin vertex lambda in the U-sc_reducible channel */
    dcomplex lambda_m_dot( const int W, const int K, const int w, const int m) const;       /**< Return Hedin vertex lambda in the U-m_reducible channel */
    dcomplex lambda_d_dot( const int W, const int K, const int w, const int m) const;       /**< Return Hedin vertex lambda in the U-d_reducible channel */

    virtual dcomplex M_sc( const int W, const int K, const int w, const int m, const int wp, const int mp) const;       /**< Return U-sc-irreducible but pp-reducilbe vertex */
    virtual dcomplex M_m( const int W, const int K, const int w, const int m, const int wp, const int mp) const;        /**< Return U-m-irreducible  but xph-reducible vertex */
    virtual dcomplex M_d( const int W, const int K, const int w, const int m, const int wp, const int mp) const;        /**< Return U-d-irreducible but ph-reducible vertex */
       
    dcomplex B_irreducible_vertex_sc( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const; 	/**< Return U-sc-irreducible vertex */
    dcomplex B_irreducible_vertex_m( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const; 	/**< Return U-m-irreducible vertex */
    dcomplex B_irreducible_vertex_d( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const; 	/**< Return U-d-irreducible vertex */
      
    dcomplex nabla_d_ph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_m_ph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_m_xph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_sc_pp_to_ph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_m_xph_to_ph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_sc_pp_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_d_ph_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_m_ph_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

#if !defined(SBEa_APPROXIMATION)
    virtual dcomplex M_d_ph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp ) const;
    virtual dcomplex M_m_ph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp ) const;
    virtual dcomplex M_m_xph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp ) const;
    virtual dcomplex M_sc_pp_to_ph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const;
    virtual dcomplex M_m_xph_to_ph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const;
    virtual dcomplex M_sc_pp_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const;
    virtual dcomplex M_d_ph_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const;
    virtual dcomplex M_m_ph_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const;       
# endif //!defined(SBEa_APPROXIMATION)   

    dcomplex phi_sc_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const;
    dcomplex phi_d_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const;
    dcomplex phi_m_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const;
  
    dcomplex phi_sc_bar_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const;
    dcomplex phi_d_bar_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const;
    dcomplex phi_m_bar_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const;

    dcomplex nabla_sc_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_d_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_m_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

  
    dcomplex nabla_d_ph_to_pp_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_m_ph_to_pp_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_m_xph_to_pp_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_sc_pp_to_ph_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_m_xph_to_ph_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_sc_pp_to_xph_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_d_ph_to_xph_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;
    dcomplex nabla_m_ph_to_xph_dot( const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const;

    dcomplex eval_f(const double t) const;

}; 


namespace boost {
    namespace numeric {
	namespace odeint {
	    template<typename Model>
		struct vector_space_norm_inf< state_frg_sbe_t<Model> >
		{
		    typedef double result_type;
		    double operator()( const state_frg_sbe_t<Model> &state_vec ) const
		    {
			using namespace std; 
			return norm( state_vec ); 
		    }
		};
	}
    }
}


// include implementation
#include <../src/frg/sbe/state.tpp>
