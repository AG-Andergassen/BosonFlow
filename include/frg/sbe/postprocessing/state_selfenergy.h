#pragma once

#include <frg/state_base.h>

// polarisation sc contr 1
// polarisation d or m contr 1
// polarisation sc contr d nabla
// polarisation sc contr m nabla
// polarisation d contr sc nabla
// polarisation d contr m nabla
// polarisation m contr d nabla
// polarisation m contr sc nabla

// polarisation sc contr d Virr
// polarisation sc contr m Virr
// polarisation sc contr sc Virr
// polarisation d contr sc Virr
// polarisation d contr m Virr
// polarisation d contr d Virr
// polarisation m contr d Virr
// polarisation m contr sc Virr
// polarisation m contr m Virr


template <typename Model>
class state_postproc_t: public state_vector_frg_base_t<Model,
    gf_1p_t<Model>, gf_1p_t<Model>, gf_1p_t<Model>, gf_1p_t<Model>,  // self-energy contributions
    gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>,   // sc susceptibility contributions
    gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, // d susceptibility contributions
    gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>,  gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, // m susceptibility contributions

    gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>, // polarisation contribution from varphi
    gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>, // polarisation contribution from 1
    gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>, // nabla contributions
    gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>,
    gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>,
    gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>, // Virr contributions
    

    gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, // lambda from contribution from nablas
    gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>,
    gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>,
    gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, // lambda from contribution from double counting terms (-Vbos)
    gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>,
    gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>,
    gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, // Virr contributions
    gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model> // varphi contributions
    
> 
  /**< Type of the state vector of the ODE solver. @see arithmetic_tuple.h for ReaK wrapper */
{
 public:
    using state_base_t = state_vector_frg_base_t<Model,
	gf_1p_t<Model>, gf_1p_t<Model>, gf_1p_t<Model>, gf_1p_t<Model>,  // self-energy contributions
	gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>,   // sc susc
	gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>,  // d susc
	gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, gf_susc_t<Model>, // m susc

	gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>, // polarisation contribution from varphi
	gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>, // polarisation contribution from 1
	gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>, // nabla contributions
	gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>,
	gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>,
	gf_polarisation_t<Model>, gf_polarisation_t<Model>, gf_polarisation_t<Model>, // Vbare contributions

	gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, // nabla contributions
	gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>,
	gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>,
	gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, // nabla contributions from double counting terms
	gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>,
	gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>,
	gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, // Vbare contributions
	gf_static_lambda_t<Model>, gf_static_lambda_t<Model>, gf_static_lambda_t<Model> // varphi contributions

     >; 

    using base_t = typename state_base_t::base_t;

    using Sig_t = gf_1p_t<Model>;
    using susc_t = gf_susc_t<Model>;
    using polarisation_t = gf_polarisation_t<Model>;
    using lambda_t = gf_lambda_t<Model>;
    using static_lambda_t = gf_static_lambda_t<Model>;
    
    /**
     *  \brief Initialization of state_postproc_t<Model> state-vector
     */
    
    inline Sig_t& gf_Sig_Lam2Pi()			        { return std::get<state_base_t::offset_idx(0)>( *this ); }    
    inline const Sig_t& gf_Sig_Lam2Pi() const	        { return std::get<state_base_t::offset_idx(0)>( *this ); }    
    
    inline Sig_t& gf_Sig_sc()			        { return std::get<state_base_t::offset_idx(1)>( *this ); }    
    inline const Sig_t& gf_Sig_sc() const	        { return std::get<state_base_t::offset_idx(1)>( *this ); }    

    inline Sig_t& gf_Sig_d()			        { return std::get<state_base_t::offset_idx(2)>( *this ); }    
    inline const Sig_t& gf_Sig_d() const	        { return std::get<state_base_t::offset_idx(2)>( *this ); }    
   
    inline Sig_t& gf_Sig_m()			        { return std::get<state_base_t::offset_idx(3)>( *this ); }    
    inline const Sig_t& gf_Sig_m() const	        { return std::get<state_base_t::offset_idx(3)>( *this ); }    


    // sc   

    inline susc_t& gf_susc_sc()			        { return std::get<state_base_t::offset_idx(4)>( *this ); }    
    inline const susc_t& gf_susc_sc() const	        { return std::get<state_base_t::offset_idx(4)>( *this ); }    

    inline susc_t& gf_suscvert_sc()			        { return std::get<state_base_t::offset_idx(5)>( *this ); }    
    inline const susc_t& gf_suscvert_sc() const	        { return std::get<state_base_t::offset_idx(5)>( *this ); }    

    inline susc_t& gf_suscbubble_sc()			        { return std::get<state_base_t::offset_idx(6)>( *this ); }    
    inline const susc_t& gf_suscbubble_sc() const	        { return std::get<state_base_t::offset_idx(6)>( *this ); }    

    inline susc_t& gf_suscvert_sc_contribution_from_nabla_d()			        { return std::get<state_base_t::offset_idx(7)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_nabla_d() const	        { return std::get<state_base_t::offset_idx(7)>( *this ); }    

    inline susc_t& gf_suscvert_sc_contribution_from_M_d()			        { return std::get<state_base_t::offset_idx(8)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_M_d() const	        { return std::get<state_base_t::offset_idx(8)>( *this ); }    


    inline susc_t& gf_suscvert_sc_contribution_from_d_double_counting()			        { return std::get<state_base_t::offset_idx(9)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_d_double_counting() const	        { return std::get<state_base_t::offset_idx(9)>( *this ); }    


    inline susc_t& gf_suscvert_sc_contribution_from_nabla_m()			        { return std::get<state_base_t::offset_idx(10)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_nabla_m() const	        { return std::get<state_base_t::offset_idx(10)>( *this ); }    

    inline susc_t& gf_suscvert_sc_contribution_from_M_m()			        { return std::get<state_base_t::offset_idx(11)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_M_m() const	        { return std::get<state_base_t::offset_idx(11)>( *this ); }    


    inline susc_t& gf_suscvert_sc_contribution_from_m_double_counting()			        { return std::get<state_base_t::offset_idx(12)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_m_double_counting() const	        { return std::get<state_base_t::offset_idx(12)>( *this ); }    

    inline susc_t& gf_suscvert_sc_contribution_from_nabla_sc()			        { return std::get<state_base_t::offset_idx(13)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_nabla_sc() const	        { return std::get<state_base_t::offset_idx(13)>( *this ); }    

    inline susc_t& gf_suscvert_sc_contribution_from_M_sc()			        { return std::get<state_base_t::offset_idx(14)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_M_sc() const	        { return std::get<state_base_t::offset_idx(14)>( *this ); }    


    inline susc_t& gf_suscvert_sc_contribution_from_sc_double_counting()			        { return std::get<state_base_t::offset_idx(15)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_sc_double_counting() const	        { return std::get<state_base_t::offset_idx(15)>( *this ); }    


    inline susc_t& gf_suscvert_sc_contribution_from_bare()			        { return std::get<state_base_t::offset_idx(16)>( *this ); }    
    inline const susc_t& gf_suscvert_sc_contribution_from_bare() const	        { return std::get<state_base_t::offset_idx(16)>( *this ); }    

    // d
    inline susc_t& gf_susc_d()			        { return std::get<state_base_t::offset_idx(17)>( *this ); }    
    inline const susc_t& gf_susc_d() const	        { return std::get<state_base_t::offset_idx(17)>( *this ); }    

    inline susc_t& gf_suscvert_d()			        { return std::get<state_base_t::offset_idx(18)>( *this ); }    
    inline const susc_t& gf_suscvert_d() const	        { return std::get<state_base_t::offset_idx(18)>( *this ); }    

    inline susc_t& gf_suscbubble_d()			        { return std::get<state_base_t::offset_idx(19)>( *this ); }    
    inline const susc_t& gf_suscbubble_d() const	        { return std::get<state_base_t::offset_idx(19)>( *this ); }    

    ///

    inline susc_t& gf_suscvert_d_contribution_from_nabla_d()			        { return std::get<state_base_t::offset_idx(20)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_nabla_d() const	        { return std::get<state_base_t::offset_idx(20)>( *this ); }    

    inline susc_t& gf_suscvert_d_contribution_from_M_d()			        { return std::get<state_base_t::offset_idx(21)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_M_d() const	        { return std::get<state_base_t::offset_idx(21)>( *this ); }    


    inline susc_t& gf_suscvert_d_contribution_from_d_double_counting()			        { return std::get<state_base_t::offset_idx(22)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_d_double_counting() const	        { return std::get<state_base_t::offset_idx(22)>( *this ); }    


    inline susc_t& gf_suscvert_d_contribution_from_nabla_m()			        { return std::get<state_base_t::offset_idx(23)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_nabla_m() const	        { return std::get<state_base_t::offset_idx(23)>( *this ); }    

    inline susc_t& gf_suscvert_d_contribution_from_M_m()			        { return std::get<state_base_t::offset_idx(24)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_M_m() const	        { return std::get<state_base_t::offset_idx(24)>( *this ); }    


    inline susc_t& gf_suscvert_d_contribution_from_m_double_counting()			        { return std::get<state_base_t::offset_idx(25)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_m_double_counting() const	        { return std::get<state_base_t::offset_idx(25)>( *this ); }    


    inline susc_t& gf_suscvert_d_contribution_from_nabla_sc()			        { return std::get<state_base_t::offset_idx(26)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_nabla_sc() const	        { return std::get<state_base_t::offset_idx(26)>( *this ); }    

    inline susc_t& gf_suscvert_d_contribution_from_M_sc()			        { return std::get<state_base_t::offset_idx(27)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_M_sc() const	        { return std::get<state_base_t::offset_idx(27)>( *this ); }    


    inline susc_t& gf_suscvert_d_contribution_from_sc_double_counting()			        { return std::get<state_base_t::offset_idx(28)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_sc_double_counting() const	        { return std::get<state_base_t::offset_idx(28)>( *this ); }    


    inline susc_t& gf_suscvert_d_contribution_from_bare()			        { return std::get<state_base_t::offset_idx(29)>( *this ); }    
    inline const susc_t& gf_suscvert_d_contribution_from_bare() const	        { return std::get<state_base_t::offset_idx(29)>( *this ); }    

    ///
    

    // m   
    inline susc_t& gf_susc_m()			        { return std::get<state_base_t::offset_idx(30)>( *this ); }    
    inline const susc_t& gf_susc_m() const	        { return std::get<state_base_t::offset_idx(30)>( *this ); }    

    inline susc_t& gf_suscvert_m()			        { return std::get<state_base_t::offset_idx(31)>( *this ); }    
    inline const susc_t& gf_suscvert_m() const	        { return std::get<state_base_t::offset_idx(31)>( *this ); }    

    inline susc_t& gf_suscbubble_m()			        { return std::get<state_base_t::offset_idx(32)>( *this ); }    
    inline const susc_t& gf_suscbubble_m() const	        { return std::get<state_base_t::offset_idx(32)>( *this ); }    

    inline susc_t& gf_suscvert_m_contribution_from_nabla_d()			        { return std::get<state_base_t::offset_idx(33)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_nabla_d() const	        { return std::get<state_base_t::offset_idx(33)>( *this ); }    

    inline susc_t& gf_suscvert_m_contribution_from_M_d()			        { return std::get<state_base_t::offset_idx(34)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_M_d() const	        { return std::get<state_base_t::offset_idx(34)>( *this ); }    


    inline susc_t& gf_suscvert_m_contribution_from_d_double_counting()			        { return std::get<state_base_t::offset_idx(35)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_d_double_counting() const	        { return std::get<state_base_t::offset_idx(35)>( *this ); }    


    inline susc_t& gf_suscvert_m_contribution_from_nabla_m()			        { return std::get<state_base_t::offset_idx(36)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_nabla_m() const	        { return std::get<state_base_t::offset_idx(36)>( *this ); }    

    inline susc_t& gf_suscvert_m_contribution_from_M_m()			        { return std::get<state_base_t::offset_idx(37)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_M_m() const	        { return std::get<state_base_t::offset_idx(37)>( *this ); }    


    inline susc_t& gf_suscvert_m_contribution_from_m_double_counting()			        { return std::get<state_base_t::offset_idx(38)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_m_double_counting() const	        { return std::get<state_base_t::offset_idx(38)>( *this ); }    

    inline susc_t& gf_suscvert_m_contribution_from_nabla_sc()			        { return std::get<state_base_t::offset_idx(39)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_nabla_sc() const	        { return std::get<state_base_t::offset_idx(39)>( *this ); }    

    inline susc_t& gf_suscvert_m_contribution_from_M_sc()			        { return std::get<state_base_t::offset_idx(40)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_M_sc() const	        { return std::get<state_base_t::offset_idx(40)>( *this ); }    


    inline susc_t& gf_suscvert_m_contribution_from_sc_double_counting()			        { return std::get<state_base_t::offset_idx(41)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_sc_double_counting() const	        { return std::get<state_base_t::offset_idx(41)>( *this ); }    

    inline susc_t& gf_suscvert_m_contribution_from_bare()			        { return std::get<state_base_t::offset_idx(42)>( *this ); }    
    inline const susc_t& gf_suscvert_m_contribution_from_bare() const	        { return std::get<state_base_t::offset_idx(42)>( *this ); }    


    inline polarisation_t& gf_polarisation_sc_contribution_from_varphi()			        { return std::get<state_base_t::offset_idx(43)>( *this ); }    
    inline const polarisation_t& gf_polarisation_sc_contribution_from_varphi() const	        { return std::get<state_base_t::offset_idx(43)>( *this ); }    

    inline polarisation_t& gf_polarisation_d_contribution_from_varphi()			        { return std::get<state_base_t::offset_idx(44)>( *this ); }    
    inline const polarisation_t& gf_polarisation_d_contribution_from_varphi() const	        { return std::get<state_base_t::offset_idx(44)>( *this ); }    

    inline polarisation_t& gf_polarisation_m_contribution_from_varphi()			        { return std::get<state_base_t::offset_idx(45)>( *this ); }    
    inline const polarisation_t& gf_polarisation_m_contribution_from_varphi() const	        { return std::get<state_base_t::offset_idx(45)>( *this ); }    


    inline polarisation_t& gf_polarisation_sc_contribution_from_1()			        { return std::get<state_base_t::offset_idx(46)>( *this ); }    
    inline const polarisation_t& gf_polarisation_sc_contribution_from_1() const	        { return std::get<state_base_t::offset_idx(46)>( *this ); }    

    inline polarisation_t& gf_polarisation_d_contribution_from_1()			        { return std::get<state_base_t::offset_idx(47)>( *this ); }    
    inline const polarisation_t& gf_polarisation_d_contribution_from_1() const	        { return std::get<state_base_t::offset_idx(47)>( *this ); }    

    inline polarisation_t& gf_polarisation_m_contribution_from_1()			        { return std::get<state_base_t::offset_idx(48)>( *this ); }    
    inline const polarisation_t& gf_polarisation_m_contribution_from_1() const	        { return std::get<state_base_t::offset_idx(48)>( *this ); }    


    inline polarisation_t& gf_polarisation_sc_contribution_from_nabla_sc()			        { return std::get<state_base_t::offset_idx(49)>( *this ); }    
    inline const polarisation_t& gf_polarisation_sc_contribution_from_nabla_sc() const	        { return std::get<state_base_t::offset_idx(49)>( *this ); }    

    inline polarisation_t& gf_polarisation_sc_contribution_from_nabla_d()			        { return std::get<state_base_t::offset_idx(50)>( *this ); }    
    inline const polarisation_t& gf_polarisation_sc_contribution_from_nabla_d() const	        { return std::get<state_base_t::offset_idx(50)>( *this ); }    

    inline polarisation_t& gf_polarisation_sc_contribution_from_nabla_m()			        { return std::get<state_base_t::offset_idx(51)>( *this ); }    
    inline const polarisation_t& gf_polarisation_sc_contribution_from_nabla_m() const	        { return std::get<state_base_t::offset_idx(51)>( *this ); }    


    inline polarisation_t& gf_polarisation_d_contribution_from_nabla_d()			        { return std::get<state_base_t::offset_idx(52)>( *this ); }    
    inline const polarisation_t& gf_polarisation_d_contribution_from_nabla_d() const	        { return std::get<state_base_t::offset_idx(52)>( *this ); }    

    inline polarisation_t& gf_polarisation_d_contribution_from_nabla_sc()			        { return std::get<state_base_t::offset_idx(53)>( *this ); }    
    inline const polarisation_t& gf_polarisation_d_contribution_from_nabla_sc() const	        { return std::get<state_base_t::offset_idx(53)>( *this ); }    

    inline polarisation_t& gf_polarisation_d_contribution_from_nabla_m()			        { return std::get<state_base_t::offset_idx(54)>( *this ); }    
    inline const polarisation_t& gf_polarisation_d_contribution_from_nabla_m() const	        { return std::get<state_base_t::offset_idx(54)>( *this ); }    


    inline polarisation_t& gf_polarisation_m_contribution_from_nabla_m()			        { return std::get<state_base_t::offset_idx(55)>( *this ); }    
    inline const polarisation_t& gf_polarisation_m_contribution_from_nabla_m() const	        { return std::get<state_base_t::offset_idx(55)>( *this ); }    

    inline polarisation_t& gf_polarisation_m_contribution_from_nabla_sc()			        { return std::get<state_base_t::offset_idx(56)>( *this ); }    
    inline const polarisation_t& gf_polarisation_m_contribution_from_nabla_sc() const	        { return std::get<state_base_t::offset_idx(56)>( *this ); }    

    inline polarisation_t& gf_polarisation_m_contribution_from_nabla_d()			        { return std::get<state_base_t::offset_idx(57)>( *this ); }    
    inline const polarisation_t& gf_polarisation_m_contribution_from_nabla_d() const	        { return std::get<state_base_t::offset_idx(57)>( *this ); }    

    inline polarisation_t& gf_polarisation_sc_contribution_from_Virr()			        { return std::get<state_base_t::offset_idx(58)>( *this ); }    
    inline const polarisation_t& gf_polarisation_sc_contribution_from_Virr() const	        { return std::get<state_base_t::offset_idx(58)>( *this ); }    

    inline polarisation_t& gf_polarisation_d_contribution_from_Virr()			        { return std::get<state_base_t::offset_idx(59)>( *this ); }    
    inline const polarisation_t& gf_polarisation_d_contribution_from_Virr() const	        { return std::get<state_base_t::offset_idx(59)>( *this ); }    

    inline polarisation_t& gf_polarisation_m_contribution_from_Virr()			        { return std::get<state_base_t::offset_idx(60)>( *this ); }    
    inline const polarisation_t& gf_polarisation_m_contribution_from_Virr() const	        { return std::get<state_base_t::offset_idx(60)>( *this ); }    


    // lambda
    inline static_lambda_t& gf_lambda_sc_contribution_from_nabla_d()			        { return std::get<state_base_t::offset_idx(61)>( *this ); }    
    inline const static_lambda_t& gf_lambda_sc_contribution_from_nabla_d() const	        { return std::get<state_base_t::offset_idx(61)>( *this ); }    

    inline static_lambda_t& gf_lambda_sc_contribution_from_nabla_m()			        { return std::get<state_base_t::offset_idx(62)>( *this ); }    
    inline const static_lambda_t& gf_lambda_sc_contribution_from_nabla_m() const	        { return std::get<state_base_t::offset_idx(62)>( *this ); }    

    inline static_lambda_t& gf_lambda_sc_contribution_from_nabla_sc()			        { return std::get<state_base_t::offset_idx(63)>( *this ); }    
    inline const static_lambda_t& gf_lambda_sc_contribution_from_nabla_sc() const	        { return std::get<state_base_t::offset_idx(63)>( *this ); }    



    inline static_lambda_t& gf_lambda_d_contribution_from_nabla_d()			        { return std::get<state_base_t::offset_idx(64)>( *this ); }    
    inline const static_lambda_t& gf_lambda_d_contribution_from_nabla_d() const	        { return std::get<state_base_t::offset_idx(64)>( *this ); }    

    inline static_lambda_t& gf_lambda_d_contribution_from_nabla_sc()			        { return std::get<state_base_t::offset_idx(65)>( *this ); }    
    inline const static_lambda_t& gf_lambda_d_contribution_from_nabla_sc() const	        { return std::get<state_base_t::offset_idx(65)>( *this ); }    

    inline static_lambda_t& gf_lambda_d_contribution_from_nabla_m()			        { return std::get<state_base_t::offset_idx(66)>( *this ); }    
    inline const static_lambda_t& gf_lambda_d_contribution_from_nabla_m() const	        { return std::get<state_base_t::offset_idx(66)>( *this ); }    


    inline static_lambda_t& gf_lambda_m_contribution_from_nabla_m()			        { return std::get<state_base_t::offset_idx(67)>( *this ); }    
    inline const static_lambda_t& gf_lambda_m_contribution_from_nabla_m() const	        { return std::get<state_base_t::offset_idx(67)>( *this ); }    

    inline static_lambda_t& gf_lambda_m_contribution_from_nabla_sc()			        { return std::get<state_base_t::offset_idx(68)>( *this ); }    
    inline const static_lambda_t& gf_lambda_m_contribution_from_nabla_sc() const	        { return std::get<state_base_t::offset_idx(68)>( *this ); }    

    inline static_lambda_t& gf_lambda_m_contribution_from_nabla_d()			        { return std::get<state_base_t::offset_idx(69)>( *this ); }    
    inline const static_lambda_t& gf_lambda_m_contribution_from_nabla_d() const	        { return std::get<state_base_t::offset_idx(69)>( *this ); }    


    // double counting contributions to lambda
    inline static_lambda_t& gf_lambda_sc_contribution_from_sc_double_counting()	        { return std::get<state_base_t::offset_idx(70)>( *this ); }    
    inline const static_lambda_t& gf_lambda_sc_contribution_from_sc_double_counting() const	        { return std::get<state_base_t::offset_idx(70)>( *this ); }    

    inline static_lambda_t& gf_lambda_sc_contribution_from_d_double_counting()	        { return std::get<state_base_t::offset_idx(71)>( *this ); }    
    inline const static_lambda_t& gf_lambda_sc_contribution_from_d_double_counting() const	        { return std::get<state_base_t::offset_idx(71)>( *this ); }    

    inline static_lambda_t& gf_lambda_sc_contribution_from_m_double_counting()	        { return std::get<state_base_t::offset_idx(72)>( *this ); }    
    inline const static_lambda_t& gf_lambda_sc_contribution_from_m_double_counting() const	        { return std::get<state_base_t::offset_idx(72)>( *this ); }    

    inline static_lambda_t& gf_lambda_d_contribution_from_d_double_counting()	        { return std::get<state_base_t::offset_idx(73)>( *this ); }    
    inline const static_lambda_t& gf_lambda_d_contribution_from_d_double_counting() const	        { return std::get<state_base_t::offset_idx(73)>( *this ); }    

    inline static_lambda_t& gf_lambda_d_contribution_from_sc_double_counting()	        { return std::get<state_base_t::offset_idx(74)>( *this ); }    
    inline const static_lambda_t& gf_lambda_d_contribution_from_sc_double_counting() const	        { return std::get<state_base_t::offset_idx(74)>( *this ); } 

    inline static_lambda_t& gf_lambda_d_contribution_from_m_double_counting()	        { return std::get<state_base_t::offset_idx(75)>( *this ); }    
    inline const static_lambda_t& gf_lambda_d_contribution_from_m_double_counting() const	        { return std::get<state_base_t::offset_idx(75)>( *this ); }    

    inline static_lambda_t& gf_lambda_m_contribution_from_m_double_counting()	        { return std::get<state_base_t::offset_idx(76)>( *this ); }    
    inline const static_lambda_t& gf_lambda_m_contribution_from_m_double_counting() const	        { return std::get<state_base_t::offset_idx(76)>( *this ); }    

    inline static_lambda_t& gf_lambda_m_contribution_from_sc_double_counting()	        { return std::get<state_base_t::offset_idx(77)>( *this ); }    
    inline const static_lambda_t& gf_lambda_m_contribution_from_sc_double_counting() const	        { return std::get<state_base_t::offset_idx(77)>( *this ); } 

    inline static_lambda_t& gf_lambda_m_contribution_from_d_double_counting()	        { return std::get<state_base_t::offset_idx(78)>( *this ); }    
    inline const static_lambda_t& gf_lambda_m_contribution_from_d_double_counting() const	        { return std::get<state_base_t::offset_idx(78)>( *this ); }    




    inline static_lambda_t& gf_lambda_sc_contribution_from_bare()			        { return std::get<state_base_t::offset_idx(79)>( *this ); }    
    inline const static_lambda_t& gf_lambda_sc_contribution_from_bare() const	        { return std::get<state_base_t::offset_idx(79)>( *this ); }    

    inline static_lambda_t& gf_lambda_d_contribution_from_bare()			        { return std::get<state_base_t::offset_idx(80)>( *this ); }    
    inline const static_lambda_t& gf_lambda_d_contribution_from_bare() const	        { return std::get<state_base_t::offset_idx(80)>( *this ); }    

    inline static_lambda_t& gf_lambda_m_contribution_from_bare()			        { return std::get<state_base_t::offset_idx(81)>( *this ); }    
    inline const static_lambda_t& gf_lambda_m_contribution_from_bare() const	        { return std::get<state_base_t::offset_idx(81)>( *this ); }    

    inline static_lambda_t& gf_lambda_sc_contribution_from_varphi()			        { return std::get<state_base_t::offset_idx(82)>( *this ); }    
    inline const static_lambda_t& gf_lambda_sc_contribution_from_varphi() const	        { return std::get<state_base_t::offset_idx(82)>( *this ); }    

    inline static_lambda_t& gf_lambda_d_contribution_from_varphi()			        { return std::get<state_base_t::offset_idx(83)>( *this ); }    
    inline const static_lambda_t& gf_lambda_d_contribution_from_varphi() const	        { return std::get<state_base_t::offset_idx(83)>( *this ); }    

    inline static_lambda_t& gf_lambda_m_contribution_from_varphi()			        { return std::get<state_base_t::offset_idx(84)>( *this ); }    
    inline const static_lambda_t& gf_lambda_m_contribution_from_varphi() const	        { return std::get<state_base_t::offset_idx(84)>( *this ); }    


    
 state_postproc_t():
    state_base_t()
        {}; 
    
 state_postproc_t( const state_postproc_t<Model>& state_vec ):
    state_base_t( state_vec )
    {}; 
    
 state_postproc_t( state_postproc_t<Model>&& state_vec ):
    state_base_t( std::move(state_vec) )
	{}; 
    
    state_postproc_t<Model>& operator=( const state_postproc_t<Model>& state_vec ){
	if( this == &state_vec)
	    return *this;
	state_base_t::operator=( state_vec );
    }
 
    state_postproc_t<Model>& operator=( state_postproc_t<Model>&& state_vec ){
	state_base_t::operator=( std::move( state_vec ) );
    }
    
    static constexpr int offset_idx(const unsigned idx)
    {
	// number of new components introduced by this state vector
	return state_base_t::offset_idx(idx + 85);
    }

    void init_0();         /*< initialization to zero or bare quantities */

    dcomplex susc_sc( const int W, const int K, const int n_in, const int n_out ) const; 

    dcomplex susc_d(  const int W, const int K, const int n_in, const int n_out ) const;

    dcomplex susc_m( const int W, const int K, const int n_in, const int n_out ) const;

};
