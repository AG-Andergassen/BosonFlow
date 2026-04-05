#include <mymath.h>
#include <models/square_hubbard.h>
#include <tu_projections.h>
#include <frg/flows.h>
#include <frg/sbe/symmetries.h>

#include <iostream>

#define SQU(a) ((a)*(a))

template<typename Model, typename... OtherTypes > 
void state_frg_sbe_t<Model, OtherTypes...>::init_bare(){
    //initialise state with vertices being the bare couplings

    state_base_t::init_bare();

    gf_w_sc().init( [this](const idx_w_t<Model> &idx) -> dcomplex {
        return bos_vertex_sc_bare(idx(IW::W), idx(IW::K), FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_start);
	});
    gf_w_d().init( [this](const idx_w_t<Model> &idx) -> dcomplex {
        return bos_vertex_d_bare(idx(IW::W), idx(IW::K), FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_start);
	});
    gf_w_m().init( [this](const idx_w_t<Model> &idx) -> dcomplex {
        return bos_vertex_m_bare(idx(IW::W), idx(IW::K), FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, fRGFlowScheme<Model>::ChosenFlowParametrizationInfo().t_start);
	});

    gf_lambda_sc().init( [this](const idx_lambda_t<Model> &idx){
	    return this->lambda_sc_bare(idx(ILAMBDA::W), idx(ILAMBDA::K), idx(ILAMBDA::w), idx(ILAMBDA::m));
	}  );
    gf_lambda_d().init(  [this](const idx_lambda_t<Model> &idx){
	    return this->lambda_d_bare(idx(ILAMBDA::W), idx(ILAMBDA::K), idx(ILAMBDA::w), idx(ILAMBDA::m));
	}   );
    gf_lambda_m().init(  [this](const idx_lambda_t<Model> &idx){
	    return this->lambda_m_bare(idx(ILAMBDA::W), idx(ILAMBDA::K), idx(ILAMBDA::w), idx(ILAMBDA::m));
	}   );

    gf_M_sc().init( [](const idx_M_t<Model> &idx){return 0;} );
    gf_M_d().init( [](const idx_M_t<Model> &idx){return 0;} );
    gf_M_m().init( [](const idx_M_t<Model> &idx){return 0;} );

    gf_f().init( [](const idx_f_t<Model> &idx){return 0;} );
}

template<typename Model, typename... OtherTypes > 
void state_frg_sbe_t<Model, OtherTypes...>::init_zero(){
    //initialise state vector with zeroes

    state_base_t::init_zero();

    gf_w_sc().init( [](const idx_w_t<Model> &idx){ return 0; });
    gf_w_d().init( [](const idx_w_t<Model> &idx){ return 0; });
    gf_w_m().init( [](const idx_w_t<Model> &idx){ return 0; });
    
    gf_lambda_sc().init( [](const idx_lambda_t<Model> &idx){ return 0; });
    gf_lambda_d().init( [](const idx_lambda_t<Model> &idx){ return 0; });
    gf_lambda_m().init( [](const idx_lambda_t<Model> &idx){ return 0; });
    
    gf_M_sc().init( [](const idx_M_t<Model> &idx){return 0;} );
    gf_M_d().init( [](const idx_M_t<Model> &idx){return 0;} );
    gf_M_m().init( [](const idx_M_t<Model> &idx){return 0;} );

    gf_f().init( [](const idx_f_t<Model> &idx){return 0;} );
}

template<typename Model, typename... OtherTypes > 
void state_frg_sbe_t<Model, OtherTypes...>::initialise_projected_Ms_and_nablas()
{
    if (!m_projected_nablas.are_initialised){
	m_projected_nablas.sc_pp_to_ph_ptr = new gf_nabla_t<Model>(m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count);
	m_projected_nablas.sc_pp_to_xph_ptr = new gf_nabla_t<Model>(m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count);

	m_projected_nablas.d_ph_to_xph_ptr = new gf_nabla_t<Model>(m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count);
	m_projected_nablas.d_ph_to_pp_ptr = new gf_nabla_t<Model>(m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count);

	m_projected_nablas.m_ph_to_xph_ptr = new gf_nabla_t<Model>(m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count);
	m_projected_nablas.m_ph_to_pp_ptr = new gf_nabla_t<Model>(m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count);
	m_projected_nablas.m_xph_to_ph_ptr = new gf_nabla_t<Model>(m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count);
	m_projected_nablas.m_xph_to_pp_ptr = new gf_nabla_t<Model>(m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count);
	m_projected_nablas.are_initialised = true;
    }

#if !defined(SBEa_APPROXIMATION)
    if (!m_projected_Ms.are_initialised){
	m_projected_Ms.sc_pp_to_ph_ptr = new gf_M_t<Model>(m_frequency_ranges.projected_M_bosonic_positive_freqs_count, m_frequency_ranges.projected_M_fermionic_positive_freqs_count);
	m_projected_Ms.sc_pp_to_xph_ptr = new gf_M_t<Model>(m_frequency_ranges.projected_M_bosonic_positive_freqs_count, m_frequency_ranges.projected_M_fermionic_positive_freqs_count);

	m_projected_Ms.d_ph_to_xph_ptr = new gf_M_t<Model>(m_frequency_ranges.projected_M_bosonic_positive_freqs_count, m_frequency_ranges.projected_M_fermionic_positive_freqs_count);
	m_projected_Ms.d_ph_to_pp_ptr = new gf_M_t<Model>(m_frequency_ranges.projected_M_bosonic_positive_freqs_count, m_frequency_ranges.projected_M_fermionic_positive_freqs_count);

	m_projected_Ms.m_ph_to_xph_ptr = new gf_M_t<Model>(m_frequency_ranges.projected_M_bosonic_positive_freqs_count, m_frequency_ranges.projected_M_fermionic_positive_freqs_count);
	m_projected_Ms.m_ph_to_pp_ptr = new gf_M_t<Model>(m_frequency_ranges.projected_M_bosonic_positive_freqs_count, m_frequency_ranges.projected_M_fermionic_positive_freqs_count);
	m_projected_Ms.m_xph_to_ph_ptr = new gf_M_t<Model>(m_frequency_ranges.projected_M_bosonic_positive_freqs_count, m_frequency_ranges.projected_M_fermionic_positive_freqs_count);
	m_projected_Ms.m_xph_to_pp_ptr = new gf_M_t<Model>(m_frequency_ranges.projected_M_bosonic_positive_freqs_count, m_frequency_ranges.projected_M_fermionic_positive_freqs_count);

	m_projected_Ms.are_initialised = true;
    }
#endif 

}


template<typename Model, typename... OtherTypes > 
state_frg_sbe_t<Model, OtherTypes...>::~state_frg_sbe_t()
{
  if (m_projected_Ms.are_initialised){
    	delete m_projected_Ms.sc_pp_to_ph_ptr;
	delete m_projected_Ms.sc_pp_to_xph_ptr;

	delete m_projected_Ms.d_ph_to_xph_ptr;
	delete m_projected_Ms.d_ph_to_pp_ptr;

	delete m_projected_Ms.m_ph_to_xph_ptr;
	delete m_projected_Ms.m_ph_to_pp_ptr;
	delete m_projected_Ms.m_xph_to_ph_ptr;
	delete m_projected_Ms.m_xph_to_pp_ptr;
  }
  
  if (m_projected_nablas.are_initialised){
    	delete m_projected_nablas.sc_pp_to_ph_ptr;
	delete m_projected_nablas.sc_pp_to_xph_ptr;

	delete m_projected_nablas.d_ph_to_xph_ptr;
	delete m_projected_nablas.d_ph_to_pp_ptr;

	delete m_projected_nablas.m_ph_to_xph_ptr;
	delete m_projected_nablas.m_ph_to_pp_ptr;
	delete m_projected_nablas.m_xph_to_ph_ptr;
	delete m_projected_nablas.m_xph_to_pp_ptr;
  }
}

template<typename Model, typename... OtherTypes > 
void state_frg_sbe_t<Model, OtherTypes...>::precalculate_projected_Ms_and_nablas(const double t )
{
    initialise_projected_Ms_and_nablas();

    std::cout << "Precalculating projected nabla functions... " << std::endl;
    m_projected_nablas.are_calculated = false;


#pragma omp parallel default(shared) 
    {
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_ph_Ptr()->init_batched(*m_projected_nablas.sc_pp_to_ph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_sc_pp_to_ph(idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_xph_Ptr()->init_batched(*m_projected_nablas.sc_pp_to_xph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_sc_pp_to_xph(idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });

    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_xph_Ptr()->init_batched(*m_projected_nablas.m_ph_to_xph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_m_ph_to_xph(idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_ph_Ptr()->init_batched(*m_projected_nablas.m_xph_to_ph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_m_xph_to_ph(idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_pp_Ptr()->init_batched(*m_projected_nablas.m_ph_to_pp_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_m_ph_to_pp(idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_pp_Ptr()->init_batched(*m_projected_nablas.m_xph_to_pp_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_m_xph_to_pp(idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_xph_Ptr()->init_batched(*m_projected_nablas.d_ph_to_xph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_d_ph_to_xph(idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_pp_Ptr()->init_batched(*m_projected_nablas.d_ph_to_pp_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_d_ph_to_pp(idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    }
    
    m_projected_nablas.are_calculated = true;


    m_projected_Ms.are_calculated = false;
#if !defined(SBEa_APPROXIMATION)
    std::cout << "Precalculating projected M functions... " << std::endl;
    
#pragma omp parallel default(shared) 
    {
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched(*m_projected_Ms.sc_pp_to_ph_ptr, [&](idx_M_t<Model> idx){return this->M_sc_pp_to_ph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched(*m_projected_Ms.sc_pp_to_xph_ptr, [&](idx_M_t<Model> idx){return this->M_sc_pp_to_xph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });

    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched(*m_projected_Ms.m_ph_to_xph_ptr, [&](idx_M_t<Model> idx){return this->M_m_ph_to_xph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched(*m_projected_Ms.m_xph_to_ph_ptr, [&](idx_M_t<Model> idx){return this->M_m_xph_to_ph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched(*m_projected_Ms.m_ph_to_pp_ptr, [&](idx_M_t<Model> idx){return this->M_m_ph_to_pp(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched(*m_projected_Ms.m_xph_to_pp_ptr, [&](idx_M_t<Model> idx){return this->M_m_xph_to_pp(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched(*m_projected_Ms.d_ph_to_xph_ptr, [&](idx_M_t<Model> idx){return this->M_d_ph_to_xph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched(*m_projected_Ms.d_ph_to_pp_ptr, [&](idx_M_t<Model> idx){return this->M_d_ph_to_pp(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    }
    
    m_projected_Ms.are_calculated = true;
# endif

}

template<typename Model, typename... OtherTypes > 
void state_frg_sbe_t<Model, OtherTypes...>::precalculate_projected_Ms_and_nablas_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const double t )
{
    initialise_projected_Ms_and_nablas();

    std::cout << "Precalculating projected nabla functions... " << std::endl;
    m_projected_nablas.are_calculated = false;

    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_ph_Ptr()->init_batched(*m_projected_nablas.sc_pp_to_ph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_sc_pp_to_ph_dot(state, idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_xph_Ptr()->init_batched(*m_projected_nablas.sc_pp_to_xph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_sc_pp_to_xph_dot(state, idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });

    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_xph_Ptr()->init_batched(*m_projected_nablas.m_ph_to_xph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_m_ph_to_xph_dot(state, idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_ph_Ptr()->init_batched(*m_projected_nablas.m_xph_to_ph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_m_xph_to_ph_dot(state, idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_pp_Ptr()->init_batched(*m_projected_nablas.m_ph_to_pp_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_m_ph_to_pp_dot(state, idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_pp_Ptr()->init_batched(*m_projected_nablas.m_xph_to_pp_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_m_xph_to_pp_dot(state, idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_xph_Ptr()->init_batched(*m_projected_nablas.d_ph_to_xph_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_d_ph_to_xph_dot(state, idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_projected_nabla_pp_Ptr()->init_batched(*m_projected_nablas.d_ph_to_pp_ptr, [&](idx_nabla_t<Model> idx){return this->nabla_d_ph_to_pp_dot(state, idx(INABLA::W), idx(INABLA::K), idx(INABLA::w), idx(INABLA::m), idx(INABLA::wp), idx(INABLA::mp), t); });


    
    m_projected_nablas.are_calculated = true;

    
    m_projected_Ms.are_calculated = false;
#if !defined(SBEa_APPROXIMATION)
    std::cout << "Precalculating projected M functions... " << std::endl;

#pragma omp parallel default(shared) 
    {
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched(*m_projected_Ms.sc_pp_to_ph_ptr, [&](idx_M_t<Model> idx){return this->M_sc_pp_to_ph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched(*m_projected_Ms.sc_pp_to_xph_ptr, [&](idx_M_t<Model> idx){return this->M_sc_pp_to_xph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });

    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched(*m_projected_Ms.m_ph_to_xph_ptr, [&](idx_M_t<Model> idx){return this->M_m_ph_to_xph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr()->init_batched(*m_projected_Ms.m_xph_to_ph_ptr, [&](idx_M_t<Model> idx){return this->M_m_xph_to_ph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched(*m_projected_Ms.m_ph_to_pp_ptr, [&](idx_M_t<Model> idx){return this->M_m_ph_to_pp(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched(*m_projected_Ms.m_xph_to_pp_ptr, [&](idx_M_t<Model> idx){return this->M_m_xph_to_pp(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr()->init_batched(*m_projected_Ms.d_ph_to_xph_ptr, [&](idx_M_t<Model> idx){return this->M_d_ph_to_xph(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr()->init_batched(*m_projected_Ms.d_ph_to_pp_ptr, [&](idx_M_t<Model> idx){return this->M_d_ph_to_pp(idx(IM::W), idx(IM::K), idx(IM::w), idx(IM::m), idx(IM::wp), idx(IM::mp)); });
    }

    m_projected_Ms.are_calculated = true;
# endif

}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_sc_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{

    if(m==0 && mp==0)
	return Model::vertex_4pt_bare( w + div2_ceil(W), div2_floor(W) - w - 1, wp + div2_ceil(W), -00,
				   0, 0, 0, -00,
				   -00, -00, -00, -00) * Model::MomentumGrid().ptr->get_volume(); // this is a volume factor
    else
       return 0.;
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_d_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
    if(m==0 && mp==0) 
	return Model::vertex_4pt_bare( w - div2_floor(W), wp + div2_ceil(W), w + div2_ceil(W), -00,
				   0, 0, 0, -00,
				   -00, -00, -00, -00) * Model::MomentumGrid().ptr->get_volume();
    else
       return 0.;
}

// used in conventional SDE
template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_m_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{ 
    return bos_vertex_m_bare( W, K, w, m, wp, mp, t ) + ferm_vertex_m_bare( W, K, w, m, wp, mp, t );

    if(m==0 && mp==0) // w2, w4, w3, w1
	return -Model::vertex_4pt_bare( w - div2_floor(W), wp + div2_ceil(W), wp - div2_floor(W), -00,
				    0, 0, 0, -00,
				    -00,-00,-00,-00) * Model::MomentumGrid().ptr->get_volume();
    else
       return 0.;
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::local_vertex_4pt_bare( const int W, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::vertex_local_part_bare(W, w, m, wp, mp);
#else
    return Model::vertex_local_part_bare_of_t(W, w, m, wp, mp, t);
#endif
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::local_vertex_4pt_bare_dot( const int W, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*Model::vertex_local_part_bare(W, w, m, wp, mp);
#else
    return Model::vertex_local_part_bare_of_t_dot(W, w, m, wp, mp, t);
#endif
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::bos_vertex_sc_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(Model::B_sc(W, K, w, m, wp, mp) - BOSONIC_INTERACTION_SC_SHIFT * Model::MomentumGrid().ptr->get_volume());
#else
    return Model::B_sc_of_t(W, K, w, m, wp, mp, t) - fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*BOSONIC_INTERACTION_SC_SHIFT * Model::MomentumGrid().ptr->get_volume();
#endif
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::bos_vertex_d_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(Model::B_d(W, K, w, m, wp, mp) - BOSONIC_INTERACTION_D_SHIFT * Model::MomentumGrid().ptr->get_volume());
#else
    return Model::B_d_of_t(W, K, w, m, wp, mp, t) - fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*BOSONIC_INTERACTION_D_SHIFT * Model::MomentumGrid().ptr->get_volume();
#endif
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::bos_vertex_m_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{ 
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(Model::B_m(W, K, w, m, wp, mp) - BOSONIC_INTERACTION_M_SHIFT * Model::MomentumGrid().ptr->get_volume());
#else
    return Model::B_m_of_t(W, K, w, m, wp, mp, t) - fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*BOSONIC_INTERACTION_M_SHIFT * Model::MomentumGrid().ptr->get_volume();
#endif
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::bos_vertex_sc_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*(Model::B_sc(W, K, w, m, wp, mp) - BOSONIC_INTERACTION_SC_SHIFT * Model::MomentumGrid().ptr->get_volume());
#else
    return Model::B_sc_of_t_dot(W, K, w, m, wp, mp, t) - fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*BOSONIC_INTERACTION_SC_SHIFT * Model::MomentumGrid().ptr->get_volume();
#endif
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::bos_vertex_d_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*(Model::B_d(W, K, w, m, wp, mp) - BOSONIC_INTERACTION_D_SHIFT * Model::MomentumGrid().ptr->get_volume());
#else
    return Model::B_d_of_t_dot(W, K, w, m, wp, mp, t) - fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*BOSONIC_INTERACTION_D_SHIFT * Model::MomentumGrid().ptr->get_volume();
#endif
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::bos_vertex_m_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{ 
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*(Model::B_m(W, K, w, m, wp, mp) - BOSONIC_INTERACTION_M_SHIFT * Model::MomentumGrid().ptr->get_volume());
#else
    return Model::B_m_of_t_dot(W, K, w, m, wp, mp, t) - fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*BOSONIC_INTERACTION_M_SHIFT * Model::MomentumGrid().ptr->get_volume();
#endif
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::ferm_vertex_sc_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(Model::F_sc(W, K, w, m, wp, mp));
#else
    return Model::F_sc_of_t(W, K, w, m, wp, mp, t);
#endif
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::ferm_vertex_d_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(Model::F_d(W, K, w, m, wp, mp));
#else
    return Model::F_d_of_t(W, K, w, m, wp, mp, t);
#endif
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::ferm_vertex_m_bare( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{ 
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(Model::F_m(W, K, w, m, wp, mp));
#else
    return Model::F_m_of_t(W, K, w, m, wp, mp, t);
#endif
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::ferm_vertex_sc_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*Model::F_sc(W, K, w, m, wp, mp);
#else
    return Model::F_sc_of_t_dot(W, K, w, m, wp, mp, t);
#endif
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::ferm_vertex_d_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*Model::F_d(W, K, w, m, wp, mp);
#else
    return Model::F_d_of_t_dot(W, K, w, m, wp, mp, t);
#endif
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::ferm_vertex_m_bare_dot( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{ 
#ifndef FREE_U_FLOW
    return fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*Model::F_m(W, K, w, m, wp, mp);
#else
    return Model::F_m_of_t_dot(W, K, w, m, wp, mp, t);
#endif
}



template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::lambda_sc_bare( const int W, const int K, const int w, const int m ) const
{
    // useful for frequency asymptotics

#if defined(INTRP_FLOW) || defined(INVERSE_INTRP_FLOW) || defined(INVERSE_INTRP_FLOW_OMEGA)
    // if e.g. DMF2RG is used, we can assume the (better) DMFT lambda is the bare value
    if( (unsigned)( W + Input_lambda_PositiveBosFreqCount - 1 ) < ( 2 * (InputUirreducibleVertexPositiveBosFreqCount - 1) + 1 )   &&
	(unsigned)( w + Input_lambda_PositiveFermFreqCount - 1 ) < ( 2 * (Input_lambda_PositiveFermFreqCount - 1) - (W+100000)%2) ){
        return InterpolatingFlow_lambda_sc()[W][K][w][m];
    }
#endif

    if( m==0 ) 
       return dcomplex(1.,0.);
    else
       return 0.;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::lambda_d_bare( const int W, const int K, const int w, const int m ) const
{ 

#if defined(INTRP_FLOW) || defined(INVERSE_INTRP_FLOW) || defined(INVERSE_INTRP_FLOW_OMEGA)
    // if e.g. DMF2RG is used, we can assume the (better) DMFT lambda is the bare value
    if( (unsigned)( W + Input_lambda_PositiveBosFreqCount - 1 ) < ( 2 * (InputUirreducibleVertexPositiveBosFreqCount - 1) + 1 )   &&
	(unsigned)( w + Input_lambda_PositiveFermFreqCount - 1 ) < ( 2 * (Input_lambda_PositiveFermFreqCount - 1) - (W+100000)%2) ){
        return InterpolatingFlow_lambda_d()[W][K][w][m];
    }
#endif

    if( m==0 ) 
       return dcomplex(1.,0.);
    else
       return 0.;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::lambda_m_bare( const int W, const int K, const int w, const int m ) const
{ 
#if defined(INTRP_FLOW) || defined(INVERSE_INTRP_FLOW) || defined(INVERSE_INTRP_FLOW_OMEGA)
    // if e.g. DMF2RG is used, we can assume the (better) DMFT lambda is the bare value
    if( (unsigned)( W + Input_lambda_PositiveBosFreqCount - 1 ) < ( 2 * (InputUirreducibleVertexPositiveBosFreqCount - 1) + 1 )   &&
	(unsigned)( w + Input_lambda_PositiveFermFreqCount - 1 ) < ( 2 * (Input_lambda_PositiveFermFreqCount - 1) - (W+100000)%2) ){
        return InterpolatingFlow_lambda_m()[W][K][w][m];
    }
#endif

    if( m==0 ) 
       return dcomplex(1.,0.);
    else
       return 0.;
}


// --- VERTEX FUNCTION IN THE PURELY FERMIONIC NOTATION ON THE FINE MOMENTUM GRID---

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_DC( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out, const double t ) const
{
#ifdef BARE_INTERACTION_NOT_DENSITY_DENSITY // models with bare interaction not pure density-density have double counting terms that are not simply -2U_{loc}
  #ifdef FREE_U_FLOW
      return Model::vertex_DC_of_t(w1_in, w2_in, w1_out, -00, p1_in, p2_in, p1_out, -00, -00, -00, -00, -00, t);
  #else
      return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::vertex_DC(w1_in, w2_in, w1_out, -00, p1_in, p2_in, p1_out, -00, -00, -00, -00, -00);
  #endif
#else
  dcomplex bare_U = local_vertex_4pt_bare(0, 0, 0, 0, 0, t )/Model::MomentumGrid().ptr->get_volume();
  return - 2. * bare_U;
#endif
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_DC_sc( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifdef BARE_INTERACTION_NOT_DENSITY_DENSITY
  #ifdef FREE_U_FLOW
         return Model::vertex_DC_sc_of_t(W, K, w, m, wp, mp, t);
  #else
	 return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::vertex_DC_sc(W, K, w, m, wp, mp);
  #endif
#else
     return - 2. * local_vertex_4pt_bare(0, 0, 0, 0, 0, t );
#endif
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_DC_d( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifdef BARE_INTERACTION_NOT_DENSITY_DENSITY
  #ifdef FREE_U_FLOW
         return Model::vertex_DC_d_of_t(W, K, w, m, wp, mp, t);
  #else
	 return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::vertex_DC_d(W, K, w, m, wp, mp);
  #endif
#else
     return - 2. * local_vertex_4pt_bare(0, 0, 0, 0, 0, t );
#endif
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_DC_m( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const
{
#ifdef BARE_INTERACTION_NOT_DENSITY_DENSITY
  #ifdef FREE_U_FLOW
         return Model::vertex_DC_m_of_t(W, K, w, m, wp, mp, t);
  #else
	 return fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::vertex_DC_m(W, K, w, m, wp, mp);
  #endif
#else
     return + 2. * local_vertex_4pt_bare(0, 0, 0, 0, 0, t );
#endif
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_4pt( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out, const double t ) const
{

    dcomplex val = mom_nabla_sc( w1_in, w2_in, w1_out, p1_in, p2_in, p1_out, t)
       + 0.5*mom_nabla_d( w1_in, w2_in, w1_out, p1_in, p2_in, p1_out, t)
       - 0.5*mom_nabla_m_ph( w1_in, w2_in, w1_out, p1_in, p2_in, p1_out, t)
       - mom_nabla_m( w1_in, w2_in, w1_out, p1_in, p2_in, p1_out, t)
      +vertex_DC( w1_in, w2_in, w1_out, p1_in, p2_in, p1_out, t );
      
    val -= fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(- BOSONIC_INTERACTION_SC_SHIFT - 1.5*BOSONIC_INTERACTION_M_SHIFT - 0.5*BOSONIC_INTERACTION_D_SHIFT);

#if defined(INTRP_FLOW) || defined(INVERSE_INTRP_FLOW) || defined(INVERSE_INTRP_FLOW_OMEGA)
    // todo: .. make the number of freqs a variable 

   // superconducting
   //   int W = w1_in + w2_in + 1;
   //int w    = w1_in - div2_ceil(W);
   //int wp   = w1_out - div2_ceil(W);

   // magnetic 
   int W = w2_in - w1_out;
   int w = w1_in + div2_floor(W);
   int wp =  w2_in - div2_ceil(W);
   
   if( (unsigned)( W + InputUirreducibleVertexPositiveBosFreqCount - 1 ) < ( 2*(InputUirreducibleVertexPositiveBosFreqCount-1) + 1 )   &&
       (unsigned)( w + InputUirreducibleVertexPositiveFermFreqCount - 1 ) < ( 2 *(InputUirreducibleVertexPositiveFermFreqCount-1 ) - (W+100000)%2) &&
       (unsigned)( wp + InputUirreducibleVertexPositiveFermFreqCount - 1 ) < ( 2 * (InputUirreducibleVertexPositiveFermFreqCount-1 ) - (W+100000)%2) ){
       
     //int p1_in_plus_p2_in = Model::SumFineMomentaIdxes(p1_in, p2_in);
       
     //int K = Model::GetCoarseMomentumIdxFromFine(p1_in_plus_p2_in);
       
     // magnetic
     
     int p2_in_minus_p1_out = Model::SumFineMomentaIdxes(p2_in, Model::GetNegativeFineMomentumIdx(p1_out));
     
     int K = Model::GetCoarseMomentumIdxFromFine(p2_in_minus_p1_out);
   
       
     // add the U-irreducible vertex of the correlated starting point
     //val += InterpolatingFlow_Uirreducible_sc()[W][K][w][0][wp][0]/Model::MomentumGrid().ptr->get_volume();
     
     // remark: we assume the input only an s-wave component (like DMFT), so we evaluate at m = 0
     // magnetic
     
     //val += InterpolatingFlow_Uirreducible_m()[W][K][w][0][wp][0]/Model::MomentumGrid().ptr->get_volume();
     val += InterpolatingFlow_Uirreducible_m()[W][0][w][0][wp][0]/Model::MomentumGrid().ptr->get_volume();
     
     
   }
#endif
   
   
#if !defined(SBEa_APPROXIMATION)
   val +=  mom_M_sc( w1_in, w2_in, w1_out, p1_in, p2_in, p1_out) + 0.5* mom_M_d( w1_in, w2_in, w1_out, p1_in, p2_in, p1_out) - 0.5* mom_M_m_ph( w1_in, w2_in, w1_out, p1_in, p2_in, p1_out) - mom_M_m( w1_in, w2_in, w1_out, p1_in, p2_in, p1_out);
# endif //!defined(SBEa_APPROXIMATION)
   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_sc( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const{
    return lambda_sc(W, K, w, m) * w_sc(W, K, t) * lambda_sc(W, K, wp, mp);
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_d( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const{
    return lambda_d(W, K, w, m) * w_d(W, K, t) * lambda_d(W, K, wp, mp);
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_m( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const{
    return lambda_m(W, K, w, m) * w_m(W, K, t) * lambda_m(W, K, wp, mp);
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_sc( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const{
    return nabla_sc(W, K, w, m, wp, mp, t) + B_irreducible_vertex_sc(W, K, w, m, wp, mp, t);
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_d( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const{
    return nabla_d(W, K, w, m, wp, mp, t) + B_irreducible_vertex_d(W, K, w, m, wp, mp, t);
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::vertex_m( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const{
    return nabla_m(W, K, w, m, wp, mp, t) + B_irreducible_vertex_m(W, K, w, m, wp, mp, t);
}




/**
 * \brief two-particle reducible vertices in the purely fermionic notation on the fine momentum grid
 */
template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_lambda( dcomplex (state_frg_sbe_t<Model, OtherTypes...>::*lambda_func)( const int, const int, const int, const int) const , const int W, const int K, const int w, const int p) const
{
    dcomplex temp( 0.0, 0.0 );
    for (int m=0; m< Model::GetMomentumFormFactorsCount(); ++m)
	temp += std::conj(Model::GetFormFactorInFineMomentumIdxSpace(m)[p])  * ((this->*lambda_func)( W, K, w, m));

    return temp;
} 

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_lambda_sc( const int W, const int K, const int w, const int p) const
{
    return mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_sc, W, K, w, p);
} 

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_lambda_d( const int W, const int K, const int w, const int p) const
{
    return mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_d, W, K, w, p);
} 

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_lambda_m( const int W, const int K, const int w, const int p) const
{
    return mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_m, W, K, w, p);
} 

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_lambda_sc_dot( const int W, const int K, const int w, const int p) const
{
    return mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_sc_dot, W, K, w, p);
} 

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_lambda_d_dot( const int W, const int K, const int w, const int p) const
{
    return mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_d_dot, W, K, w, p);
} 

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_lambda_m_dot( const int W, const int K, const int w, const int p) const
{
    return mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_m_dot, W, K, w, p);
} 


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_nabla_sc( const int w1_in, const int w2_in, const int w1_out, const int p1_in, const int p2_in, const int p1_out, const double t ) const
{    
   int W_pp = w1_in + w2_in + 1;
   int w    = w1_in - div2_ceil(W_pp);
   int wp   = w1_out - div2_ceil(W_pp);
 
   int p1_in_plus_p2_in = Model::SumFineMomentaIdxes(p1_in, p2_in);

   int K_pp = Model::GetCoarseMomentumIdxFromFine(p1_in_plus_p2_in);

   dcomplex temp_l=mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_sc, W_pp, K_pp, w, p1_in);
   dcomplex temp_r=mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_sc, W_pp, K_pp, wp, p1_out);

   return temp_l * w_sc( W_pp, K_pp, t) * temp_r; 
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_nabla_d( int w1_in, int w2_in, int w1_out, int p1_in, int p2_in, int p1_out, const double t ) const
{
   int W_ph = w1_out - w1_in;
   int w = w1_in + div2_floor(W_ph);
   int wp =  w2_in - div2_ceil(W_ph);

   int p1_out_minus_p1_in = Model::SumFineMomentaIdxes(p1_out,Model::GetNegativeFineMomentumIdx(p1_in));
 
   int K_ph = Model::GetCoarseMomentumIdxFromFine(p1_out_minus_p1_in);
 
   int p1_in_plus_p2_in_minus_p1_out = Model::SumFineMomentaIdxes
     (Model::SumFineMomentaIdxes(p1_in, p2_in),
      Model::GetNegativeFineMomentumIdx(p1_out));
 
   dcomplex temp_l=mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_d, W_ph, K_ph, w, p1_in);
   dcomplex temp_r=mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_d, W_ph, K_ph, wp, p1_in_plus_p2_in_minus_p1_out);

   return temp_l * w_d( W_ph, K_ph, t) * temp_r; 
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_nabla_m_ph( int w1_in, int w2_in, int w1_out, int p1_in, int p2_in, int p1_out, const double t ) const
{
   int W_ph = w1_out - w1_in;
   int w = w1_in + div2_floor(W_ph);
   int wp =  w2_in - div2_ceil(W_ph);

   int p1_out_minus_p1_in = Model::SumFineMomentaIdxes(
	 p1_out,
	 Model::GetNegativeFineMomentumIdx(p1_in));
 
   int K_ph = Model::GetCoarseMomentumIdxFromFine(p1_out_minus_p1_in);
   
   int p1_in_plus_p2_in_minus_p1_out = Model::SumFineMomentaIdxes(
       Model::SumFineMomentaIdxes(p1_in, p2_in),
       Model::GetNegativeFineMomentumIdx(p1_out));

   dcomplex temp_l=mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_m, W_ph, K_ph, w, p1_in);
   dcomplex temp_r=mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_m, W_ph, K_ph, wp, p1_in_plus_p2_in_minus_p1_out);


   return temp_l * w_m( W_ph, K_ph, t) * temp_r;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_nabla_m( int w1_in, int w2_in, int w1_out, int p1_in, int p2_in, int p1_out, const double t ) const{
    
    dcomplex val( 0.0, 0.0 );
    
    int W_xph = w2_in - w1_out;
    int w = w1_in + div2_floor(W_xph);
    int wp =  w2_in - div2_ceil(W_xph);
 
    int p2_in_minus_p1_out = Model::SumFineMomentaIdxes(
       p2_in,
       Model::GetNegativeFineMomentumIdx(p1_out));

    int K_xph = Model::GetCoarseMomentumIdxFromFine(p2_in_minus_p1_out);
    
    dcomplex temp_l=mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_m, W_xph, K_xph, w, p1_in);
    dcomplex temp_r=mom_lambda( &state_frg_sbe_t<Model, OtherTypes...>::lambda_m, W_xph, K_xph, wp, p1_out);

    return temp_l * w_m( W_xph, K_xph, t) * temp_r; 
}

#if !defined(SBEa_APPROXIMATION)
template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_M_sc( int w1_in, int w2_in, int w1_out, int p1_in, int p2_in, int p1_out ) const
{
   int W_pp = w1_in + w2_in + 1;
   int w    = w1_in - div2_ceil(W_pp);
   int wp   = w1_out - div2_ceil(W_pp);

   dcomplex val( 0.0, 0.0 );
   if((unsigned)(W_pp +m_frequency_ranges.M_bosonic_positive_freqs_count) > m_frequency_ranges.M_bosonic_freqs_count ||
      (unsigned)(w + m_frequency_ranges.M_fermionic_positive_freqs_count ) > (m_frequency_ranges.M_fermionic_freqs_count - (W_pp+100000)%2) ||
      (unsigned)(wp + m_frequency_ranges.M_fermionic_positive_freqs_count ) > (m_frequency_ranges.M_fermionic_freqs_count - (W_pp+100000)%2))
      return val;

   
   int p1_in_plus_p2_in = Model::SumFineMomentaIdxes(p1_in, p2_in);

   int K_pp = Model::GetCoarseMomentumIdxFromFine(p1_in_plus_p2_in);


   for (int m=0; m< Model::GetMomentumFormFactorsCount(); ++m)
   for (int mp=0; mp< Model::GetMomentumFormFactorsCount(); ++mp){
       auto &f_m = Model::GetFormFactorInFineMomentumIdxSpace(m);
       auto &f_mp = Model::GetFormFactorInFineMomentumIdxSpace(mp);
	   
       val += std::conj(f_m[p1_in]) * f_mp[p1_out] * M_sc( W_pp, K_pp, w, m, wp, mp );
   }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_M_d( int w1_in, int w2_in, int w1_out, int p1_in, int p2_in, int p1_out ) const
{
   int W_ph = w1_out - w1_in;
   int w = w1_in + div2_floor(W_ph);
   int wp =  w2_in - div2_ceil(W_ph);

   dcomplex val( 0.0, 0.0 );
   if((unsigned)(W_ph +m_frequency_ranges.M_bosonic_positive_freqs_count) > m_frequency_ranges.M_bosonic_freqs_count ||
      (unsigned)(w + m_frequency_ranges.M_fermionic_positive_freqs_count ) > (m_frequency_ranges.M_fermionic_freqs_count - (W_ph+100000)%2) ||
      (unsigned)(w + m_frequency_ranges.M_fermionic_positive_freqs_count ) > (m_frequency_ranges.M_fermionic_freqs_count - (W_ph+100000)%2))
      return val;

   int p1_out_minus_p1_in = Model::SumFineMomentaIdxes(
     p1_out,
     Model::GetNegativeFineMomentumIdx(p1_in));
   
   int K_ph = Model::GetCoarseMomentumIdxFromFine(p1_out_minus_p1_in);

   int p1_in_plus_p2_in_minus_p1_out = Model::SumFineMomentaIdxes
     (Model::SumFineMomentaIdxes(p1_in, p2_in),
      Model::GetNegativeFineMomentumIdx(p1_out));

   for (int m=0; m< Model::GetMomentumFormFactorsCount(); ++m)
   for (int mp=0; mp< Model::GetMomentumFormFactorsCount(); ++mp){
       auto &f_m = Model::GetFormFactorInFineMomentumIdxSpace(m);
       auto &f_mp = Model::GetFormFactorInFineMomentumIdxSpace(mp);
	  
       val += std::conj(f_m[p1_in]) * f_mp[p1_in_plus_p2_in_minus_p1_out]  * M_d( W_ph, K_ph, w, m, wp, mp );
   }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_M_m_ph( int w1_in, int w2_in, int w1_out, int p1_in, int p2_in, int p1_out ) const
{
   int W_ph = w1_out - w1_in;
   int w = w1_in + div2_floor(W_ph);
   int wp =  w2_in - div2_ceil(W_ph);

   dcomplex val( 0.0, 0.0 );
   if((unsigned)(W_ph +m_frequency_ranges.M_bosonic_positive_freqs_count) > m_frequency_ranges.M_bosonic_freqs_count ||
      (unsigned)(w + m_frequency_ranges.M_fermionic_positive_freqs_count ) > (m_frequency_ranges.M_fermionic_freqs_count - (W_ph+100000)%2) ||
      (unsigned)(w + m_frequency_ranges.M_fermionic_positive_freqs_count ) > (m_frequency_ranges.M_fermionic_freqs_count - (W_ph+100000)%2))
      return val;

   
   int p1_out_minus_p1_in = Model::SumFineMomentaIdxes
     (p1_out,
      Model::GetNegativeFineMomentumIdx(p1_in));
   
   int K_ph = Model::GetCoarseMomentumIdxFromFine(p1_out_minus_p1_in);

   int p1_in_plus_p2_in_minus_p1_out = Model::SumFineMomentaIdxes
     (Model::SumFineMomentaIdxes(p1_in,p2_in),
      Model::GetNegativeFineMomentumIdx(p1_out));

   for (int m=0; m< Model::GetMomentumFormFactorsCount(); ++m)
   for (int mp=0; mp< Model::GetMomentumFormFactorsCount(); ++mp){
       auto &f_m = Model::GetFormFactorInFineMomentumIdxSpace(m);
       auto &f_mp = Model::GetFormFactorInFineMomentumIdxSpace(mp);
       
       val += std::conj(f_m[p1_in]) * f_mp[p1_in_plus_p2_in_minus_p1_out] *  M_m( W_ph, K_ph, w, m, wp, mp );
   }

   return val;
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::mom_M_m( int w1_in, int w2_in, int w1_out, int p1_in, int p2_in, int p1_out ) const
{
   int W_xph = w2_in - w1_out;
   int w = w1_in + div2_floor(W_xph);
   int wp =  w2_in - div2_ceil(W_xph);

   dcomplex val( 0.0, 0.0 );
   if((unsigned)(W_xph +m_frequency_ranges.M_bosonic_positive_freqs_count) > m_frequency_ranges.M_bosonic_freqs_count ||
      (unsigned)(w + m_frequency_ranges.M_fermionic_positive_freqs_count ) > (m_frequency_ranges.M_fermionic_freqs_count - (W_xph+100000)%2) ||
      (unsigned)(wp + m_frequency_ranges.M_fermionic_positive_freqs_count ) > (m_frequency_ranges.M_fermionic_freqs_count - (W_xph+100000)%2))
      return val;

   int p2_in_minus_p1_out = Model::SumFineMomentaIdxes(
         p2_in,
	 Model::GetNegativeFineMomentumIdx(p1_out));


   int K_xph = Model::GetCoarseMomentumIdxFromFine(p2_in_minus_p1_out);

   for (int m=0; m< Model::GetMomentumFormFactorsCount(); ++m)
   for (int mp=0; mp< Model::GetMomentumFormFactorsCount(); ++mp){
       auto &f_m = Model::GetFormFactorInFineMomentumIdxSpace(m);
       auto &f_mp = Model::GetFormFactorInFineMomentumIdxSpace(mp);
       
       val += std::conj(f_m[p1_in]) * f_mp[p1_out] * M_m( W_xph, K_xph, w, m, wp, mp );
   }
   
   return val;
}
# endif


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::w_sc( const int W, const int K, const double t) const{
#ifdef STATIC_CALCULATION
    return gf_w_sc()[0][K];
#endif

 if( (unsigned)(W + m_frequency_ranges.w_bosonic_positive_freqs_count) < m_frequency_ranges.w_bosonic_freqs_count )  // TODO: if( likely())
     return gf_w_sc()[W][K];
 return bos_vertex_sc_bare(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, t);
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::w_d( const int W, const int K, const double t) const{
#ifdef STATIC_CALCULATION
    return gf_w_d()[0][K];
#endif

 if( (unsigned)(W + m_frequency_ranges.w_bosonic_positive_freqs_count) < m_frequency_ranges.w_bosonic_freqs_count )
     return gf_w_d()[W][K];
 return bos_vertex_d_bare(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, t);
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::w_m( const int W, const int K, const double t) const{
#ifdef STATIC_CALCULATION
    return gf_w_m()[0][K];
#endif

 if( (unsigned)(W + m_frequency_ranges.w_bosonic_positive_freqs_count) < m_frequency_ranges.w_bosonic_freqs_count )
     return gf_w_m()[W][K];
 return bos_vertex_m_bare(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, t);
}

// Remark: the "dot" quantities are used in the occasion the state object is a _derivative_ of the state. The distinction is needed because e.g. w and w_dot have different asymptotics (one decays to the bare interaction, the other decays to zero)
template<typename Model, typename... OtherTypes >
dcomplex state_frg_sbe_t<Model, OtherTypes...>::w_sc_dot( const int W, const int K, const double t) const{
#ifdef STATIC_CALCULATION
    return gf_w_sc()[0][K];
#endif

 if( (unsigned)(W + m_frequency_ranges.w_bosonic_positive_freqs_count) < m_frequency_ranges.w_bosonic_freqs_count )  // TODO: if( likely())
     return gf_w_sc()[W][K];

 if (m_force_w_dot_asymptotics_to_zero)
   return 0.0;
 
 return bos_vertex_sc_bare_dot(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, t);
}

template<typename Model, typename... OtherTypes >
dcomplex state_frg_sbe_t<Model, OtherTypes...>::w_d_dot( const int W, const int K, const double t) const{
#ifdef STATIC_CALCULATION
    return gf_w_d()[0][K];
#endif

 if( (unsigned)(W + m_frequency_ranges.w_bosonic_positive_freqs_count) < m_frequency_ranges.w_bosonic_freqs_count )
     return gf_w_d()[W][K];

 if (m_force_w_dot_asymptotics_to_zero) // for multiloop correction objects, this would be set to true, to not overcount asymptotics
   return 0.0;
 
 
 return bos_vertex_d_bare_dot(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, t);
}

template<typename Model, typename... OtherTypes >
dcomplex state_frg_sbe_t<Model, OtherTypes...>::w_m_dot( const int W, const int K, const double t) const{
#ifdef STATIC_CALCULATION
    return gf_w_m()[0][K];
#endif

 if( (unsigned)(W + m_frequency_ranges.w_bosonic_positive_freqs_count) < m_frequency_ranges.w_bosonic_freqs_count )
     return gf_w_m()[W][K];
 
 //std::cout << fRGFlowScheme<Model>::U_MultiplicativeCutoff_dot(t)*(Model::B_m(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0) + BOSONIC_INTERACTION_M_SHIFT * Model::MomentumGrid().ptr->get_volume()) << std::endl;

 if (m_force_w_dot_asymptotics_to_zero)
   return 0.0;
 
 
 return bos_vertex_m_bare_dot(W, K, FrequenciesCount::LARGE, 0, FrequenciesCount::LARGE, 0, t);
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::lambda_sc( const int W, const int K, const int w, const int m) const
{
#ifdef STATIC_CALCULATION
    return gf_lambda_sc()[0][K][0][m];
#endif

   if( (unsigned)(W +m_frequency_ranges.lambda_bosonic_positive_freqs_count) < m_frequency_ranges.lambda_bosonic_freqs_count && 
     (unsigned)(w + m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
       return gf_lambda_sc()[W][K][w][m];
   return this->lambda_sc_bare(W,K,w,m); 
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::lambda_d( const int W, const int K, const int w, const int m) const
{
#ifdef STATIC_CALCULATION
    return gf_lambda_d()[0][K][0][m];
#endif

   if( (unsigned)(W +m_frequency_ranges.lambda_bosonic_positive_freqs_count) < m_frequency_ranges.lambda_bosonic_freqs_count && 
     (unsigned)(w + m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
       return gf_lambda_d()[W][K][w][m];
   return this->lambda_d_bare(W,K,w,m); 
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::lambda_m( const int W, const int K, const int w, const int m) const
{
#ifdef STATIC_CALCULATION
    return gf_lambda_m()[0][K][0][m];
#endif
   if( (unsigned)(W +m_frequency_ranges.lambda_bosonic_positive_freqs_count) < m_frequency_ranges.lambda_bosonic_freqs_count && 
     (unsigned)(w + m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
   {
       return gf_lambda_m()[W][K][w][m];
   }
   return this->lambda_m_bare(W,K,w,m);
}   

template<typename Model, typename... OtherTypes >
dcomplex state_frg_sbe_t<Model, OtherTypes...>::lambda_sc_dot( const int W, const int K, const int w, const int m) const
{
#ifdef STATIC_CALCULATION
    return gf_lambda_sc()[0][K][0][m];
#endif

   if( (unsigned)(W +m_frequency_ranges.lambda_bosonic_positive_freqs_count) < m_frequency_ranges.lambda_bosonic_freqs_count &&
     (unsigned)(w + m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
       return gf_lambda_sc()[W][K][w][m];
   return 0.0;
}

template<typename Model, typename... OtherTypes >
dcomplex state_frg_sbe_t<Model, OtherTypes...>::lambda_d_dot( const int W, const int K, const int w, const int m) const
{
#ifdef STATIC_CALCULATION
    return gf_lambda_d()[0][K][0][m];
#endif

   if( (unsigned)(W +m_frequency_ranges.lambda_bosonic_positive_freqs_count) < m_frequency_ranges.lambda_bosonic_freqs_count &&
     (unsigned)(w + m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
       return gf_lambda_d()[W][K][w][m];
   return 0.0;
}

template<typename Model, typename... OtherTypes >
dcomplex state_frg_sbe_t<Model, OtherTypes...>::lambda_m_dot( const int W, const int K, const int w, const int m) const
{
#ifdef STATIC_CALCULATION
    return gf_lambda_m()[0][K][0][m];
#endif

   if( (unsigned)(W +m_frequency_ranges.lambda_bosonic_positive_freqs_count) < m_frequency_ranges.lambda_bosonic_freqs_count &&
     (unsigned)(w + m_frequency_ranges.lambda_fermionic_positive_freqs_count ) < (m_frequency_ranges.lambda_fermionic_freqs_count - (W+100000)%2))
   {
       return gf_lambda_m()[W][K][w][m];
   }
   return 0.0;
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_sc( const int W, const int K, const int w, const int m, const int wp, const int mp) const
{
#ifdef STATIC_CALCULATION
    return gf_M_sc()[0][K][0][m][0][mp];
#endif


   dcomplex val(0.0,0.0);
   if( (unsigned)( W + m_frequency_ranges.M_bosonic_positive_freqs_count ) < (m_frequency_ranges.M_bosonic_freqs_count)   &&
       (unsigned)(w + m_frequency_ranges.M_fermionic_positive_freqs_count ) < (m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2) &&
       (unsigned)(wp + m_frequency_ranges.M_fermionic_positive_freqs_count ) < (m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2))
   {
#ifdef BOSONISE_M_LAZY
       const int v0 = 0 + std::abs(div2_floor(W));
       //const int minus_v0 = -1 - div2_floor(W) - (W % 2);
       if (v0 >= m_frequency_ranges.M_fermionic_positive_freqs_count)
	   return val;
       const dcomplex w_of_M_sc_plus = 0.5*(gf_M_sc()[W][K][v0][m][0][mp] + gf_M_sc()[W][K][v0][m][-1][mp]);
       const dcomplex w_of_M_sc_minus = 0.5*(gf_M_sc()[W][K][v0][m][0][mp] - gf_M_sc()[W][K][v0][m][-1][mp]);

       const dcomplex lambda_of_M_sc_plus_left_times_w = 0.5*(gf_M_sc()[W][K][w][m][0][mp] + gf_M_sc()[W][K][w][m][-1][mp]);
       const dcomplex lambda_of_M_sc_minus_left_times_w = 0.5*(gf_M_sc()[W][K][w][m][0][mp] - gf_M_sc()[W][K][w][m][-1][mp]);

       const dcomplex lambda_of_M_sc_plus_right = 0.5*(gf_M_sc()[W][K][wp][m][0][mp] + gf_M_sc()[W][K][wp][m][-1][mp])/w_of_M_sc_plus;
       const dcomplex lambda_of_M_sc_minus_right = 0.5*(gf_M_sc()[W][K][wp][m][0][mp] - gf_M_sc()[W][K][wp][m][-1][mp])/w_of_M_sc_minus;
       
       val += (abs(w_of_M_sc_plus) < 1e-7) ? 0 : lambda_of_M_sc_plus_left_times_w * lambda_of_M_sc_plus_right;
       val += (abs(w_of_M_sc_minus) < 1e-7) ? 0 : lambda_of_M_sc_minus_left_times_w * lambda_of_M_sc_minus_right;
#else
       val = gf_M_sc()[W][K][w][m][wp][mp];
#endif
   }
   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_d( const int W, const int K, const int w, const int m, const int wp, const int mp) const
{
#ifdef STATIC_CALCULATION
    return gf_M_d()[0][K][0][m][0][mp];;
#endif

   dcomplex val(0.0,0.0);
   if( (unsigned)( W + m_frequency_ranges.M_bosonic_positive_freqs_count ) < (m_frequency_ranges.M_bosonic_freqs_count)   &&
       (unsigned)(w + m_frequency_ranges.M_fermionic_positive_freqs_count) < (m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2) &&
       (unsigned)(wp + m_frequency_ranges.M_fermionic_positive_freqs_count) < (m_frequency_ranges.M_fermionic_freqs_count- (W+100000)%2))
   {       
#ifdef BOSONISE_M_LAZY
       const int v0 = 0 + std::abs(div2_floor(W));
       //const int minus_v0 = -1 - div2_floor(W) - (W % 2);
       if (v0 >= m_frequency_ranges.M_fermionic_positive_freqs_count)
	   return val;
       const dcomplex w_of_M_d_plus = 0.5*(gf_M_d()[W][K][v0][m][0][mp] + gf_M_d()[W][K][v0][m][-1][mp]);
       const dcomplex w_of_M_d_minus = 0.5*(gf_M_d()[W][K][v0][m][0][mp] - gf_M_d()[W][K][v0][m][-1][mp]);

       const dcomplex lambda_of_M_d_plus_left_times_w = 0.5*(gf_M_d()[W][K][w][m][0][mp] + gf_M_d()[W][K][w][m][-1][mp]);
       const dcomplex lambda_of_M_d_minus_left_times_w = 0.5*(gf_M_d()[W][K][w][m][0][mp] - gf_M_d()[W][K][w][m][-1][mp]);

       const dcomplex lambda_of_M_d_plus_right = 0.5*(gf_M_d()[W][K][wp][m][0][mp] + gf_M_d()[W][K][wp][m][-1][mp])/w_of_M_d_plus;
       const dcomplex lambda_of_M_d_minus_right = 0.5*(gf_M_d()[W][K][wp][m][0][mp] - gf_M_d()[W][K][wp][m][-1][mp])/w_of_M_d_minus;
       
       val += (w_of_M_d_plus.real() == 0) ? 0 : lambda_of_M_d_plus_left_times_w * lambda_of_M_d_plus_right;
       val += (w_of_M_d_minus.real() == 0) ? 0 : lambda_of_M_d_minus_left_times_w * lambda_of_M_d_minus_right;
#else
       val = gf_M_d()[W][K][w][m][wp][mp];
#endif

   }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_m( const int W, const int K, const int w, const int m, const int wp, const int mp) const
{
#ifdef STATIC_CALCULATION
    return gf_M_m()[0][K][0][m][0][mp];;
#endif

   dcomplex val(0.0,0.0);
   if( (unsigned)( W + m_frequency_ranges.M_bosonic_positive_freqs_count ) < (m_frequency_ranges.M_bosonic_freqs_count)   &&
       (unsigned)(w + m_frequency_ranges.M_fermionic_positive_freqs_count) < (m_frequency_ranges.M_fermionic_freqs_count - (W+100000)%2) &&
       (unsigned)(wp + m_frequency_ranges.M_fermionic_positive_freqs_count) < (m_frequency_ranges.M_fermionic_freqs_count -(W+100000)%2) )
   {
#ifdef BOSONISE_M_LAZY
       const int v0 = 0 + std::abs(div2_floor(W));
       //const int minus_v0 = -1 - div2_floor(W) - (W % 2);
       if (v0 >= m_frequency_ranges.M_fermionic_positive_freqs_count)
	   return val;
       const dcomplex w_of_M_m_plus = 0.5*(gf_M_m()[W][K][v0][m][0][mp] + gf_M_m()[W][K][v0][m][-1][mp]);
       const dcomplex w_of_M_m_minus = 0.5*(gf_M_m()[W][K][v0][m][0][mp] - gf_M_m()[W][K][v0][m][-1][mp]);

       const dcomplex lambda_of_M_m_plus_left_times_w = 0.5*(gf_M_m()[W][K][w][m][0][mp] + gf_M_m()[W][K][w][m][-1][mp]);
       const dcomplex lambda_of_M_m_minus_left_times_w = 0.5*(gf_M_m()[W][K][w][m][0][mp] - gf_M_m()[W][K][w][m][-1][mp]);

       const dcomplex lambda_of_M_m_plus_right = 0.5*(gf_M_m()[W][K][wp][m][0][mp] + gf_M_m()[W][K][wp][m][-1][mp])/w_of_M_m_plus;
       const dcomplex lambda_of_M_m_minus_right = 0.5*(gf_M_m()[W][K][wp][m][0][mp] - gf_M_m()[W][K][wp][m][-1][mp])/w_of_M_m_minus;
       
       val += (w_of_M_m_plus.real() == 0) ? 0 : lambda_of_M_m_plus_left_times_w * lambda_of_M_m_plus_right;
       val += (w_of_M_m_minus.real() == 0) ? 0 : lambda_of_M_m_minus_left_times_w * lambda_of_M_m_minus_right;
#else
       val = gf_M_m()[W][K][w][m][wp][mp];
#endif
   }

   return val;
}

/**
 * \brief U-irreducible vertex in channels pp, xph and ph
 *
 */

// todo: check if it should be Model::B_sc local bare part here really is needed
template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::B_irreducible_vertex_sc( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
    dcomplex val(0, 0);

    // todo: treat properly the shifted interactions when the flag is used
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count)){
#endif
	val = 0.5*nabla_d_ph_to_pp(W, K, w, m, wp, mp, t) - 0.5*nabla_m_ph_to_pp(W, K, w, m, wp, mp, t) - nabla_m_xph_to_pp(W, K, w, m, wp, mp, t) +vertex_DC_sc(W, K, w, m, wp, mp, t);
	    
    // shifted value from bosonic propagator gf_w necessary if gf_w would be initialised to zero (since gf_w = 0 => gf_w does not flow) 
    val += (m == 0 && mp == 0) ? fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(BOSONIC_INTERACTION_SC_SHIFT + 0.5*BOSONIC_INTERACTION_D_SHIFT + 1.5*BOSONIC_INTERACTION_M_SHIFT) * Model::MomentumGrid().ptr->get_volume() : 0;

#ifdef PRECOMPUTE_STATE_PROJECTIONS // warning: not compatible with flexible U flow
    }else{
	val = fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::F_sc( W, K, w, m, wp, mp );
    }
#endif

#if !defined(SBEa_APPROXIMATION)
   if( FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.M_bosonic_positive_freqs_count, m_frequency_ranges.M_fermionic_positive_freqs_count) ){
       val += M_sc(W, K, w, m, wp, mp) + 0.5*M_d_ph_to_pp(W, K, w, m, wp, mp) -0.5*M_m_ph_to_pp(W, K, w, m, wp, mp) - M_m_xph_to_pp(W, K, w, m, wp, mp);
   }
# endif

#if defined(INTRP_FLOW) || defined(INVERSE_INTRP_FLOW) || defined(INVERSE_INTRP_FLOW_OMEGA)
    if( FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, InputUirreducibleVertexPositiveBosFreqCount - 1, InputUirreducibleVertexPositiveFermFreqCount - 1) ){
	// add the U-irreducible vertex of the correlated starting point
      if (m == 0 && mp == 0)
	val += InterpolatingFlow_Uirreducible_sc()[W][0][w][0][wp][0];
    }
#endif

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::B_irreducible_vertex_d( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const
{
    dcomplex val(0,0);
    
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count))
#endif    
	val = 2.0*nabla_sc_pp_to_ph(W, K, w, m, wp, mp, t) - nabla_sc_pp_to_xph(W, K, w, m, wp, mp, t) - 2.0*nabla_m_xph_to_ph(W, K, w, m, wp, mp, t) + 0.5*nabla_m_ph_to_xph(W, K, w, m, wp, mp, t) - 0.5*nabla_d_ph_to_xph(W, K, w, m, wp, mp, t) + vertex_DC_d(W, K, w, m, wp, mp, t);
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    else
	val = fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::F_d( W, K, w, m, wp, mp );
#endif


    // shifted value from bosonic propagator gf_w necessary if gf_w would be initialised to zero (since gf_w = 0 => gf_w does not flow) 
    val += (m == 0 && mp == 0) ?  fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(BOSONIC_INTERACTION_SC_SHIFT + 0.5*BOSONIC_INTERACTION_D_SHIFT + 1.5*BOSONIC_INTERACTION_M_SHIFT) * Model::MomentumGrid().ptr->get_volume() : 0;

#if !defined(SBEa_APPROXIMATION)
    if( FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.M_bosonic_positive_freqs_count, m_frequency_ranges.M_fermionic_positive_freqs_count) )
   {
      val += M_d(W, K, w, m, wp, mp)-0.5*M_d_ph_to_xph(W, K, w, m, wp, mp) -2.0*M_m_xph_to_ph(W, K, w, m, wp, mp) +0.5*M_m_ph_to_xph(W, K, w, m, wp, mp) + 2.0*M_sc_pp_to_ph(W, K, w, m, wp, mp) - M_sc_pp_to_xph(W, K, w, m, wp, mp);
   }
# endif

#if defined(INTRP_FLOW) || defined(INVERSE_INTRP_FLOW) || defined(INVERSE_INTRP_FLOW_OMEGA)
    if (FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, InputUirreducibleVertexPositiveBosFreqCount - 1, InputUirreducibleVertexPositiveFermFreqCount - 1)){
    // add the U-irreducible vertex of the correlated starting point
      if (m == 0 && mp == 0)
	val += InterpolatingFlow_Uirreducible_d()[W][0][w][0][wp][0];
}
#endif

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::B_irreducible_vertex_m( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const
{
    dcomplex val(0.0, 0.0);


#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count))
#endif
	val = -nabla_sc_pp_to_xph(W, K, w, m, wp, mp, t) - 0.5*nabla_d_ph_to_xph(W, K, w, m, wp, mp, t) + 0.5*nabla_m_ph_to_xph(W, K, w, m, wp, mp, t) + vertex_DC_m(W, K, w, m, wp, mp, t);
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    else
	val = fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*Model::F_m( W, K, w, m, wp, mp );
#endif

    // shifted value from bosonic propagator gf_w necessary if gf_w would be initialised to zero (since gf_w = 0 => gf_w does not flow) 
    val -= (m == 0 && mp == 0) ? fRGFlowScheme<Model>::U_MultiplicativeCutoff(t)*(BOSONIC_INTERACTION_SC_SHIFT + 0.5*BOSONIC_INTERACTION_D_SHIFT + 1.5*BOSONIC_INTERACTION_M_SHIFT) * Model::MomentumGrid().ptr->get_volume() : 0;

#if defined(INTRP_FLOW) || defined(INVERSE_INTRP_FLOW) || defined(INVERSE_INTRP_FLOW_OMEGA)
    if( FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, InputUirreducibleVertexPositiveBosFreqCount - 1, InputUirreducibleVertexPositiveFermFreqCount - 1) ){
    // add the U-irreducible vertex of the correlated starting point
      if (m == 0 && mp == 0)
	val += InterpolatingFlow_Uirreducible_m()[W][0][w][0][wp][0];
    }
#endif


#if !defined(SBEa_APPROXIMATION)
   if( FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.M_bosonic_positive_freqs_count, m_frequency_ranges.M_fermionic_positive_freqs_count) )
   {
      val += M_m(W, K, w, m, wp, mp) -M_sc_pp_to_xph(W, K, w, m, wp, mp) - 0.5*M_d_ph_to_xph(W, K, w, m, wp, mp) + 0.5*M_m_ph_to_xph(W, K, w, m, wp, mp);
   }
# endif
   return val;
}
/** 
 * \brief Projection of the U-reducible vertices nabla to different channels
 *
 * The multi-boson exchange M should be included in the next step
 */

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_d_ph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{

#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.d_ph_to_pp_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_mod(w + div2_ceil( W ) + div2_floor( W_ph ));
   int wp_mod(div2_floor( W ) - wp + div2_floor( W_ph ) - 1);

   auto &Matrix_ph_to_pp = TUProjections<Model>::Matrix_ph_to_pp();
   const int refined_momenta_count = Model::GetRefinedMomentaCount();
   const int form_factors_count = Model::GetMomentumFormFactorsCount();
   auto &volume_weights = Model::MomentumGrid().volume_weights;

   for(int Q = 0; Q<refined_momenta_count; ++Q)
      {
      dcomplex temp(0.0,0.0);
      for(int n = 0; n<form_factors_count; ++n)
         {
         dcomplex temp2(0.0,0.0);
         for(int np = 0; np<form_factors_count; ++np)
	     {
		temp2 += Matrix_ph_to_pp[K][Q][m][mp][n][np] * lambda_d(W_ph, Q, wp_mod, np); 
            }
         temp += lambda_d(W_ph, Q, w_mod, n) * temp2;
         }
      val += temp * w_d(W_ph,Q, t) * volume_weights[Q];
      }
   
   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_m_ph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{

#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.m_ph_to_pp_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_mod(w + div2_ceil( W ) + div2_floor( W_ph ));
   int wp_mod(div2_floor( W ) - wp + div2_floor( W_ph ) - 1);

   auto &Matrix_ph_to_pp = TUProjections<Model>::Matrix_ph_to_pp();
   const int refined_momenta_count = Model::GetRefinedMomentaCount();
   const int form_factors_count = Model::GetMomentumFormFactorsCount();
   auto &volume_weights = Model::MomentumGrid().volume_weights;    

   for(int Q = 0; Q<refined_momenta_count; ++Q)
      {
      dcomplex temp(0.0,0.0);
      for(int n = 0; n<form_factors_count; ++n)
         {
         dcomplex temp2(0.0,0.0);
         for(int np = 0; np<form_factors_count; ++np)
            {
            temp2 += Matrix_ph_to_pp[K][Q][m][mp][n][np] * lambda_m(W_ph, Q, wp_mod, np);
            }
         temp += lambda_m(W_ph, Q, w_mod, n) * temp2;
         }
      val += temp * w_m(W_ph,Q, t) * volume_weights[Q];
      }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_m_xph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.m_xph_to_pp_ptr)[W][K][w][m][wp][mp];
#endif 
    
   dcomplex val( 0.0, 0.0 );

   int W_xph = - ( W + 100000 ) % 2 - wp - w - 1;
   int w_mod(w + div2_ceil( W ) + div2_floor( W_xph ));
   int wp_mod(wp + div2_ceil( W ) + div2_floor( W_xph ));

   auto &Matrix_xph_to_pp = TUProjections<Model>::Matrix_xph_to_pp();
   const int refined_momenta_count = Model::GetRefinedMomentaCount();
   const int form_factors_count = Model::GetMomentumFormFactorsCount();
   auto &volume_weights = Model::MomentumGrid().volume_weights;
    
   for(int Q = 0; Q<refined_momenta_count; ++Q)
      {
      dcomplex temp(0.0,0.0);
      for(int n = 0; n<form_factors_count; ++n)
         {
         dcomplex temp2(0.0,0.0);
         for(int np = 0; np<form_factors_count; ++np)
            {
            temp2 += Matrix_xph_to_pp[K][Q][m][mp][n][np] * lambda_m(W_xph, Q, wp_mod, np); 
            }
         temp += lambda_m(W_xph, Q, w_mod, n) * temp2;
         }
      val += temp * w_m(W_xph,Q, t) * volume_weights[Q];
      }
   
    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_sc_pp_to_ph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.sc_pp_to_ph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_pp = w + wp + ( W + 100000 ) % 2 + 1;
   int w_mod( w - div2_floor( W ) - div2_ceil( W_pp ));
   int wp_mod(w + div2_ceil( W ) - div2_ceil( W_pp ));
   
   auto &Matrix_pp_to_ph = TUProjections<Model>::Matrix_pp_to_ph();
   const int refined_momenta_count = Model::GetRefinedMomentaCount();
   const int form_factors_count = Model::GetMomentumFormFactorsCount();
   auto &volume_weights = Model::MomentumGrid().volume_weights;
     
   for(int Q = 0; Q<refined_momenta_count; ++Q)
      {
      dcomplex temp(0.0,0.0);
      for(int n = 0; n<form_factors_count; ++n)
         {
         dcomplex temp2(0.0,0.0);
         for(int np = 0; np<form_factors_count; ++np)
            {
            temp2 += Matrix_pp_to_ph[K][Q][m][mp][n][np] * lambda_sc(W_pp, Q, wp_mod, np); 
            }
         temp += lambda_sc(W_pp, Q, w_mod, n) * temp2;
         }
      val += temp * w_sc(W_pp,Q, t) * volume_weights[Q];
      }
   
    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_m_xph_to_ph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.m_xph_to_ph_ptr)[W][K][w][m][wp][mp];
#endif 
    
   dcomplex val( 0.0, 0.0 );

   int W_xph = wp - w;
   int w_mod(w - div2_floor( W ) + div2_floor( W_xph ));
   int wp_mod(wp + div2_ceil( W ) - div2_ceil( W_xph ));

   auto &Matrix_xph_to_ph = TUProjections<Model>::Matrix_xph_to_ph();
   const int refined_momenta_count = Model::GetRefinedMomentaCount();
   const int form_factors_count = Model::GetMomentumFormFactorsCount();
   auto &volume_weights = Model::MomentumGrid().volume_weights;
    
   for(int Q = 0; Q<refined_momenta_count; ++Q)
      {
      dcomplex temp(0.0,0.0);
      for(int n = 0; n<form_factors_count; ++n)
         {
         dcomplex temp2(0.0,0.0);
         for(int np = 0; np<form_factors_count; ++np)
            {
            temp2 += Matrix_xph_to_ph[K][Q][m][mp][n][np] * lambda_m(W_xph, Q, wp_mod, np); 
            }
         temp += lambda_m(W_xph, Q, w_mod, n) * temp2;
         }
      val += temp * w_m(W_xph,Q, t) * volume_weights[Q];
      }
   
    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_sc_pp_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.sc_pp_to_xph_ptr)[W][K][w][m][wp][mp];
#endif 
    
   dcomplex val( 0.0, 0.0 );

   int W_pp = w + wp + ( W + 100000 ) % 2 + 1;
   int w_mod(w - div2_floor( W ) - div2_ceil( W_pp ));
   int wp_mod(wp - div2_floor( W ) - div2_ceil( W_pp ));

   auto &Matrix_pp_to_xph = TUProjections<Model>::Matrix_pp_to_xph();
   const int refined_momenta_count = Model::GetRefinedMomentaCount();
   const int form_factors_count = Model::GetMomentumFormFactorsCount();
   auto &volume_weights = Model::MomentumGrid().volume_weights;
    
   for(int Q = 0; Q<refined_momenta_count; ++Q)
      {
      dcomplex temp(0.0,0.0);
      for(int n = 0; n<form_factors_count; ++n)
         {
         dcomplex temp2(0.0,0.0);
         for(int np = 0; np<form_factors_count; ++np)
            {
            temp2 += Matrix_pp_to_xph[K][Q][m][mp][n][np] * lambda_sc(W_pp, Q, wp_mod, np); 
            }
         temp += lambda_sc(W_pp, Q, w_mod, n) * temp2;
         }
      val += temp * w_sc(W_pp,Q, t) * volume_weights[Q];
      }
   
    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_d_ph_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.d_ph_to_xph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_mod(w - div2_floor( W ) + div2_floor( W_ph ));
   int wp_mod(wp + div2_ceil( W ) - div2_ceil( W_ph ));

   auto &Matrix_ph_to_xph = TUProjections<Model>::Matrix_ph_to_xph();
   const int refined_momenta_count = Model::GetRefinedMomentaCount();
   const int form_factors_count = Model::GetMomentumFormFactorsCount();
   auto &volume_weights = Model::MomentumGrid().volume_weights;
    
   for(int Q = 0; Q<refined_momenta_count; ++Q)
      {
      dcomplex temp(0.0,0.0);
      for(int n = 0; n<form_factors_count; ++n)
         {
         dcomplex temp2(0.0,0.0);
         for(int np = 0; np<form_factors_count; ++np)
            {
            temp2 += Matrix_ph_to_xph[K][Q][m][mp][n][np] * lambda_d(W_ph, Q, wp_mod, np); 
            }
         temp += lambda_d(W_ph, Q, w_mod, n) * temp2;
         }
      val += temp * w_d(W_ph,Q, t) * volume_weights[Q];
      }
   
    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_m_ph_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.m_ph_to_xph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_mod(w - div2_floor( W ) + div2_floor( W_ph ));
   int wp_mod(wp + div2_ceil( W ) - div2_ceil( W_ph ));

   auto &Matrix_ph_to_xph = TUProjections<Model>::Matrix_ph_to_xph();
   const int refined_momenta_count = Model::GetRefinedMomentaCount();
   const int form_factors_count = Model::GetMomentumFormFactorsCount();
   auto &volume_weights = Model::MomentumGrid().volume_weights;

   for(int Q = 0; Q< refined_momenta_count; ++Q)
      {
      dcomplex temp(0.0,0.0);
      for(int n = 0; n < form_factors_count ; ++n)
         {
         dcomplex temp2(0.0,0.0);
         for(int np = 0; np < form_factors_count; ++np)
            {
            temp2 += Matrix_ph_to_xph[K][Q][m][mp][n][np] * lambda_m(W_ph, Q, wp_mod, np);
            }
         temp += lambda_m(W_ph, Q, w_mod, n) * temp2;
         }
      val += temp * w_m(W_ph,Q, t) * volume_weights[Q];
      }

    return val;
}

#if !defined(SBEa_APPROXIMATION)
template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_d_ph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_Ms.are_calculated)
	return (*m_projected_Ms.d_ph_to_pp_ptr)[W][K][w][m][wp][mp];
#endif 

    dcomplex val( 0.0, 0.0 );

    int W_ph = wp - w;
    int w_in_mod = w + div2_ceil( W ) + div2_floor( W_ph );
    int w_out_mod = div2_floor( W ) - wp + div2_floor( W_ph ) - 1;

    for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
        for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q){
                val+= TUProjections<Model>::Matrix_ph_to_pp()[K][Q][m][mp][n][np] *
                      M_d( W_ph, Q, w_in_mod, n, w_out_mod, np ) * Model::MomentumGrid().volume_weights[Q];
            }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_m_ph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_Ms.are_calculated)
	return (*m_projected_Ms.m_ph_to_pp_ptr)[W][K][w][m][wp][mp];
#endif 

    dcomplex val( 0.0, 0.0 );

    int W_ph = wp - w;
    int w_in_mod = w + div2_ceil( W ) + div2_floor( W_ph );
    int w_out_mod = div2_floor( W ) - wp + div2_floor( W_ph ) - 1;

    for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
        for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q){
                val+= TUProjections<Model>::Matrix_ph_to_pp()[K][Q][m][mp][n][np] *
                      M_m( W_ph, Q, w_in_mod, n, w_out_mod, np ) * Model::MomentumGrid().volume_weights[Q];
            }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_m_xph_to_pp( const int W, const int K, const int w, const int m, const int wp, const int mp ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_Ms.are_calculated)
	return (*m_projected_Ms.m_xph_to_pp_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val(0.0,0.0);

   int W_xph = - ( W + 100000 ) % 2 - wp - w - 1;
   int w_in_mod(w + div2_ceil( W ) + div2_floor( W_xph ));
   int w_out_mod(wp + div2_ceil( W ) + div2_floor( W_xph ));

    for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
        for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q){
                val+= TUProjections<Model>::Matrix_xph_to_pp()[K][Q][m][mp][n][np] *
                          M_m( W_xph, Q, w_in_mod, n, w_out_mod, np ) * Model::MomentumGrid().volume_weights[Q];
                }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_sc_pp_to_ph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_Ms.are_calculated)
	return (*m_projected_Ms.sc_pp_to_ph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_pp = w + wp + ( W + 100000 ) % 2 + 1;
   int w_in_mod( w - div2_floor( W ) - div2_ceil( W_pp )), w_out_mod(w + div2_ceil( W ) - div2_ceil( W_pp ));

   for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
       for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
           for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
                   val+= TUProjections<Model>::Matrix_pp_to_ph()[K][Q][m][mp][n][np] *
                         M_sc( W_pp, Q, w_in_mod, n, w_out_mod, np ) * Model::MomentumGrid().volume_weights[Q];

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_m_xph_to_ph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_Ms.are_calculated)
	return (*m_projected_Ms.m_xph_to_ph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_xph = wp - w;
   int w_in_mod(w - div2_floor( W ) + div2_floor( W_xph )), w_out_mod(wp + div2_ceil( W ) - div2_ceil( W_xph ));

   for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
       for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
           for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
                   val+= TUProjections<Model>::Matrix_xph_to_ph()[K][Q][m][mp][n][np] *
                         M_m( W_xph, Q, w_in_mod, n, w_out_mod, np ) * Model::MomentumGrid().volume_weights[Q];

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_sc_pp_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_Ms.are_calculated)
	return (*m_projected_Ms.sc_pp_to_xph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_pp = w + wp + ( W + 100000 ) % 2 + 1;
   int w_in_mod(w - div2_floor( W ) - div2_ceil( W_pp )), w_out_mod(wp - div2_floor( W ) - div2_ceil( W_pp ));

   for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
       for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
           for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q){
                   val+=  TUProjections<Model>::Matrix_pp_to_xph()[K][Q][m][mp][n][np] *
                          M_sc( W_pp, Q, w_in_mod, n, w_out_mod, np ) * Model::MomentumGrid().volume_weights[Q];
            }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_d_ph_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_Ms.are_calculated)
	return (*m_projected_Ms.d_ph_to_xph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_in_mod(w - div2_floor( W ) + div2_floor( W_ph ));
   int w_out_mod(wp + div2_ceil( W ) - div2_ceil( W_ph ));

   for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
       for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
           for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
                   val+= TUProjections<Model>::Matrix_ph_to_xph()[K][Q][m][mp][n][np] *
                         M_d(W_ph, Q, w_in_mod, n, w_out_mod, np ) * Model::MomentumGrid().volume_weights[Q];

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::M_m_ph_to_xph( const int W, const int K, const int w, const int m, const int wp, const int mp ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_Ms.are_calculated)
	return (*m_projected_Ms.m_ph_to_xph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_in_mod(w - div2_floor( W ) + div2_floor( W_ph ));
   int w_out_mod(wp + div2_ceil( W ) - div2_ceil( W_ph ));

   for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
       for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
           for(int Q = 0; Q< Model::GetRefinedMomentaCount(); ++Q)
                   val+= TUProjections<Model>::Matrix_ph_to_xph()[K][Q][m][mp][n][np] *
                         M_m(W_ph, Q, w_in_mod, n, w_out_mod, np ) * Model::MomentumGrid().volume_weights[Q];

   return val;
}
# endif //!defined(SBEa_APPROXIMATION)


// needed for multiloop. See equations (A.9) in C. Hille's thesis


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::phi_sc_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const
{
  
  dcomplex lambda_sc_dot_l = lambda_sc_dot(W, K, w, m);
  dcomplex lambda_sc_dot_r = lambda_sc_dot(W, K, wp, mp);
  dcomplex w_sc_dot_ = w_sc_dot(W, K, t);
  dcomplex lambda_sc_l = state.lambda_sc(W, K, w, m);
  dcomplex lambda_sc_r = state.lambda_sc(W, K, wp, mp);
  dcomplex w_sc_ = state.w_sc(W, K, t);
  
  dcomplex nabla_sc_dot = lambda_sc_l * w_sc_ * lambda_sc_dot_r +
                         lambda_sc_dot_l * w_sc_ * lambda_sc_r +
                         lambda_sc_l * w_sc_dot_ * lambda_sc_dot_r;

  dcomplex M_sc_dot = M_sc(W, K, w, m, wp, mp); // its asymptotics are also zero, so we didn't need a special method "M_sc_dot" as opposed to lambda and w
  

  return nabla_sc_dot + M_sc_dot - bos_vertex_sc_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t);
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::phi_sc_bar_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const 
{
    dcomplex val (0.0, 0.0);

#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count)){
#endif
      
	val = 0.5*nabla_d_ph_to_pp_dot(state, W, K, w, m, wp, mp, t)-
	    0.5*nabla_m_ph_to_pp_dot(state, W, K, w, m, wp, mp, t)-
	    nabla_m_xph_to_pp_dot(state, W, K, w, m, wp, mp, t );
#ifdef PRECOMPUTE_STATE_PROJECTIONS    
    }
	// todo: add appropriate asymptotics for ACTUAL_U flow
#endif
	
#if !defined(SBEa_APPROXIMATION)
    if( FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.M_bosonic_positive_freqs_count, m_frequency_ranges.M_fermionic_positive_freqs_count)){
	val += 0.5*M_d_ph_to_pp( W, K, w, m, wp, mp )-
	    0.5*M_m_ph_to_pp( W, K, w, m, wp, mp )-
	    M_m_xph_to_pp( W, K, w, m, wp, mp );
    }
# endif

    // the following is added in the case there's a flowing interaction

    // phi_sc_bar_dot ~ 0.5*dot(nabla_d + M_d - U_d) - 1.5*dot(nabla_m + M_m - U_m)

    // note: for 1l state_dot objects, w_dot decays not to 0 but to \dot{B}. For such objects, m_force_w_dot_asymptotics_to_zero variable will be set to true 
    if (!m_force_w_dot_asymptotics_to_zero) // we must subtract it back
	val -= 2.0 * local_vertex_4pt_bare_dot(W, w, m, wp, mp, t) + ferm_vertex_sc_bare_dot(W, K, w, m, wp, mp, t);//-local_vertex_4pt_bare_dot(W, w, m, wp, mp, t) + bos_vertex_m_bare_dot(W, K, w, m, wp, mp, t) - ferm_vertex_sc_bare_dot(W, K, w, m, wp, mp, t);
	//      val += -0.5*bos_vertex_d_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t) + 1.5 * bos_vertex_m_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t);

    
    return val;
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::phi_d_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const
{
  
  dcomplex lambda_d_dot_l = lambda_d_dot(W, K, w, m);
  dcomplex lambda_d_dot_r = lambda_d_dot(W, K, wp, mp);
  dcomplex w_d_dot_ = w_d_dot(W, K, t);
  dcomplex lambda_d_l = state.lambda_d(W, K, w, m);
  dcomplex lambda_d_r = state.lambda_d(W, K, wp, mp);
  dcomplex w_d_ = state.w_d(W, K, t);
  
  dcomplex nabla_d_dot = lambda_d_l * w_d_ * lambda_d_dot_r +
                         lambda_d_dot_l * w_d_ * lambda_d_r +
                         lambda_d_l * w_d_dot_ * lambda_d_dot_r;

  dcomplex M_d_dot = M_d(W, K, w, m, wp, mp); // its asymptotics are also zero, so we didn't need a special method "M_d_dot" as opposed to lambda and w

  return nabla_d_dot + M_d_dot - bos_vertex_d_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t);
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::phi_d_bar_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const
{
    dcomplex val (0.0, 0.0);

#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count)){
#endif    
	val = -2.*nabla_m_xph_to_ph_dot(state, W, K, w, m, wp, mp, t )+
	    2.*nabla_sc_pp_to_ph_dot(state, W, K, w, m, wp, mp, t )-
	    0.5*nabla_d_ph_to_xph_dot(state, W, K, w, m, wp, mp, t )+
	    0.5*nabla_m_ph_to_xph_dot(state, W, K, w, m, wp, mp, t )-
	    nabla_sc_pp_to_xph_dot(state, W, K, w, m, wp, mp, t );
#ifdef PRECOMPUTE_STATE_PROJECTIONS    
    }
    // todo: add appropriate asymptotics for ACTUAL_U flow
#endif
    

#if !defined(SBEa_APPROXIMATION)
    if( FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.M_bosonic_positive_freqs_count, m_frequency_ranges.M_fermionic_positive_freqs_count)){ 
	val += -2.*M_m_xph_to_ph( W, K, w, m, wp, mp )+
	    2.*M_sc_pp_to_ph( W, K, w, m, wp, mp )-
	    0.5*M_d_ph_to_xph( W, K, w, m, wp, mp )+
	    0.5*M_m_ph_to_xph( W, K, w, m, wp, mp )-
	    M_sc_pp_to_xph( W, K, w, m, wp, mp );
    }
# endif
    
    
    if (!m_force_w_dot_asymptotics_to_zero) // we must subtract it back
	val -= 2.0 * local_vertex_4pt_bare_dot(W, w, m, wp, mp, t) + ferm_vertex_d_bare_dot(W, K, w, m, wp, mp, t);//- 3.0 * local_vertex_4pt_bare_dot(W, w, m, wp, mp, t) + bos_vertex_sc_bare_dot(W, K, w, m, wp, mp, t) + ferm_vertex_m_bare_dot(W, K, w, m, wp, mp, t) - 2.0 * ferm_vertex_d_bare_dot(W, K, w, m, wp, mp, t);
	//val += -bos_vertex_sc_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t) + 1.5 * bos_vertex_m_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t) + 0.5 * bos_vertex_d_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t);
                          
    
    return val;
}    

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::phi_m_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const
{
  
  dcomplex lambda_m_dot_l = lambda_m_dot(W, K, w, m);
  dcomplex lambda_m_dot_r = lambda_m_dot(W, K, wp, mp);
  dcomplex w_m_dot_ = w_m_dot(W, K, t);
  
  dcomplex lambda_m_l = state.lambda_m(W, K, w, m);
  dcomplex lambda_m_r = state.lambda_m(W, K, wp, mp);
  dcomplex w_m_ = state.w_m(W, K, t);
  
  dcomplex nabla_m_dot = lambda_m_dot_l * w_m_ * lambda_m_r +
                         lambda_m_l * w_m_dot_ * lambda_m_r +
                         lambda_m_l * w_m_ * lambda_m_dot_r;
                         
  dcomplex M_m_dot = M_m(W, K, w, m, wp, mp); // its asymptotics are also zero, so we didn't need a special method "M_m_dot" as opposed to lambda and w
  
  return nabla_m_dot + M_m_dot - bos_vertex_m_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t);
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::phi_m_bar_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t) const
{
    dcomplex val (0.0, 0.0);

#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.projected_nabla_bosonic_positive_freqs_count, m_frequency_ranges.projected_nabla_fermionic_positive_freqs_count)){
#endif
	val = -0.5*nabla_d_ph_to_xph_dot(state, W, K, w, m, wp, mp, t )+
	    0.5*nabla_m_ph_to_xph_dot(state, W, K, w, m, wp, mp, t )-
	    nabla_sc_pp_to_xph_dot(state, W, K, w, m, wp, mp, t );
#ifdef PRECOMPUTE_STATE_PROJECTIONS    
    }
    // todo: add appropriate asymptotics for ACTUAL_U flow
#endif
#if !defined(SBEa_APPROXIMATION)
    if( FrequencyDependenceScheme<Model>::is_in_frequency_box_2_particle(W, w, wp, m_frequency_ranges.M_bosonic_positive_freqs_count, m_frequency_ranges.M_fermionic_positive_freqs_count)){ 
	    val += -0.5*M_d_ph_to_xph( W, K, w, m, wp, mp )+
	    0.5*M_m_ph_to_xph( W, K, w, m, wp, mp )-
	    M_sc_pp_to_xph( W, K, w, m, wp, mp );
    }
# endif

    // the following is added in the case there's a flowing interaction                                                                                                                                                          

    if (!m_force_w_dot_asymptotics_to_zero) // we must subtract it back
	val -= -2.0*local_vertex_4pt_bare_dot(W, w, m, wp, mp, t) + ferm_vertex_m_bare_dot(W, K, w, m, wp, mp, t);//local_vertex_4pt_bare_dot(W, w, m, wp, mp, t) + bos_vertex_sc_bare_dot(W, K, w, m, wp, mp, t) + 0.5 * (ferm_vertex_d_bare_dot(W, K, w, m, wp, mp, t) - ferm_vertex_m_bare_dot(W, K, w, m, wp, mp, t));
	//val += bos_vertex_sc_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t) + 0.5 * bos_vertex_d_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t) - 0.5 * bos_vertex_m_bare_dot(W, K, FrequenciesCount::LARGE, m, FrequenciesCount::LARGE, mp, t);
    
    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_d_ph_to_pp_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.d_ph_to_pp_ptr)[W][K][w][m][wp][mp];
#endif 


   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_mod(w + div2_ceil( W ) + div2_floor( W_ph ));
   int wp_mod(div2_floor( W ) - wp + div2_floor( W_ph ) - 1);

   for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
      {
      dcomplex temp1(0.0,0.0);
      dcomplex temp2(0.0,0.0);
      dcomplex temp3(0.0,0.0);
      for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
         {
         dcomplex temp4(0.0,0.0);
         dcomplex temp5(0.0,0.0);
         for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
	     {
		temp4 += TUProjections<Model>::Matrix_ph_to_pp()[K][Q][m][mp][n][np] * state.lambda_d(W_ph, Q, wp_mod, np); 
		temp5 += TUProjections<Model>::Matrix_ph_to_pp()[K][Q][m][mp][n][np] * lambda_d_dot(W_ph, Q, wp_mod, np); 
            }
         temp1 += state.lambda_d(W_ph, Q, w_mod, n) * temp5;
         temp2 += lambda_d_dot(W_ph, Q, w_mod, n) * temp4;
         temp3 += state.lambda_d(W_ph, Q, w_mod, n) * temp4;
         }
      val += ((temp1+temp2) * state.w_d(W_ph,Q, t) + temp3 * w_d_dot(W_ph, Q, t)) * Model::MomentumGrid().volume_weights[Q];
      }
   
   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_m_ph_to_pp_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.m_ph_to_pp_ptr)[W][K][w][m][wp][mp];
#endif 


   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_mod(w + div2_ceil( W ) + div2_floor( W_ph ));
   int wp_mod(div2_floor( W ) - wp + div2_floor( W_ph ) - 1);

   for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
      {
      dcomplex temp1(0.0,0.0);
      dcomplex temp2(0.0,0.0);
      dcomplex temp3(0.0,0.0);
      for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
         {
         dcomplex temp4(0.0,0.0);
         dcomplex temp5(0.0,0.0);
         for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            {
            temp4 += TUProjections<Model>::Matrix_ph_to_pp()[K][Q][m][mp][n][np] * state.lambda_m(W_ph, Q, wp_mod, np);
            temp5 += TUProjections<Model>::Matrix_ph_to_pp()[K][Q][m][mp][n][np] * lambda_m_dot(W_ph, Q, wp_mod, np);
            }
         temp1 += state.lambda_m(W_ph, Q, w_mod, n) * temp5;
         temp2 += lambda_m_dot(W_ph, Q, w_mod, n) * temp4;
         temp3 += state.lambda_m(W_ph, Q, w_mod, n) * temp4;
         }
      val += ((temp1+temp2) * state.w_m(W_ph,Q, t) + temp3 * w_m_dot(W_ph,Q, t)) * Model::MomentumGrid().volume_weights[Q];
      }

   return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_m_xph_to_pp_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{

#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.m_xph_to_pp_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_xph = - ( W + 100000 ) % 2 - wp - w - 1;
   int w_mod(w + div2_ceil( W ) + div2_floor( W_xph ));
   int wp_mod(wp + div2_ceil( W ) + div2_floor( W_xph ));
    
   for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
      {
      dcomplex temp1(0.0,0.0);
      dcomplex temp2(0.0,0.0);
      dcomplex temp3(0.0,0.0);
      for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
         {
         dcomplex temp4(0.0,0.0);
         dcomplex temp5(0.0,0.0);
         for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            {
            temp4 += TUProjections<Model>::Matrix_xph_to_pp()[K][Q][m][mp][n][np] * state.lambda_m(W_xph, Q, wp_mod, np); 
            temp5 += TUProjections<Model>::Matrix_xph_to_pp()[K][Q][m][mp][n][np] * lambda_m_dot(W_xph, Q, wp_mod, np); 
            }
         temp1 += state.lambda_m(W_xph, Q, w_mod, n) * temp5;
         temp2 += lambda_m_dot(W_xph, Q, w_mod, n) * temp4;
         temp3 += state.lambda_m(W_xph, Q, w_mod, n) * temp4;
         }
      val += ((temp1+temp2) * state.w_m(W_xph,Q, t) + temp3 * w_m_dot(W_xph,Q, t)) * Model::MomentumGrid().volume_weights[Q];
      }

    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_sc_pp_to_ph_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.sc_pp_to_ph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_pp = w + wp + ( W + 100000 ) % 2 + 1;
   int w_mod( w - div2_floor( W ) - div2_ceil( W_pp ));
   int wp_mod(w + div2_ceil( W ) - div2_ceil( W_pp ));
    
   for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
      {
      dcomplex temp1(0.0,0.0);
      dcomplex temp2(0.0,0.0);
      dcomplex temp3(0.0,0.0);
      for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
         {
         dcomplex temp4(0.0,0.0);
         dcomplex temp5(0.0,0.0);
         for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            {
            temp4 += TUProjections<Model>::Matrix_pp_to_ph()[K][Q][m][mp][n][np] * state.lambda_sc(W_pp, Q, wp_mod, np); 
            temp5 += TUProjections<Model>::Matrix_pp_to_ph()[K][Q][m][mp][n][np] * lambda_sc_dot(W_pp, Q, wp_mod, np); 
            }
         temp1 += state.lambda_sc(W_pp, Q, w_mod, n) * temp5;
         temp2 += lambda_sc_dot(W_pp, Q, w_mod, n) * temp4;
         temp3 += state.lambda_sc(W_pp, Q, w_mod, n) * temp4;
         }
      val += ((temp1+temp2) * state.w_sc(W_pp,Q, t) + temp3 * w_sc_dot(W_pp,Q, t)) * Model::MomentumGrid().volume_weights[Q];
      }
   
    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_m_xph_to_ph_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.m_xph_to_ph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_xph = wp - w;
   int w_mod(w - div2_floor( W ) + div2_floor( W_xph ));
   int wp_mod(wp + div2_ceil( W ) - div2_ceil( W_xph ));
    
   for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
      {
      dcomplex temp1(0.0,0.0);
      dcomplex temp2(0.0,0.0);
      dcomplex temp3(0.0,0.0);
      for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
         {
         dcomplex temp4(0.0,0.0);
         dcomplex temp5(0.0,0.0);
         for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            {
            temp4 += TUProjections<Model>::Matrix_xph_to_ph()[K][Q][m][mp][n][np] * state.lambda_m(W_xph, Q, wp_mod, np); 
            temp5 += TUProjections<Model>::Matrix_xph_to_ph()[K][Q][m][mp][n][np] * lambda_m_dot(W_xph, Q, wp_mod, np); 
            }
         temp1 += state.lambda_m(W_xph, Q, w_mod, n) * temp5;
         temp2 += lambda_m_dot(W_xph, Q, w_mod, n) * temp4;
         temp3 += state.lambda_m(W_xph, Q, w_mod, n) * temp4;
         }
      val += ((temp1+temp2) * state.w_m(W_xph,Q, t) + temp3 * w_m_dot(W_xph,Q, t)) * Model::MomentumGrid().volume_weights[Q];
      }
   
    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_sc_pp_to_xph_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
    
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.sc_pp_to_xph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_pp = w + wp + ( W + 100000 ) % 2 + 1;
   int w_mod(w - div2_floor( W ) - div2_ceil( W_pp ));
   int wp_mod(wp - div2_floor( W ) - div2_ceil( W_pp ));
    
   for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
      {
      dcomplex temp1(0.0,0.0);
      dcomplex temp2(0.0,0.0);
      dcomplex temp3(0.0,0.0);
      for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
         {
         dcomplex temp4(0.0,0.0);
         dcomplex temp5(0.0,0.0);
         for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            {
            temp4 += TUProjections<Model>::Matrix_pp_to_xph()[K][Q][m][mp][n][np] * state.lambda_sc(W_pp, Q, wp_mod, np); 
            temp5 += TUProjections<Model>::Matrix_pp_to_xph()[K][Q][m][mp][n][np] * lambda_sc_dot(W_pp, Q, wp_mod, np); 
            }
         temp1 += state.lambda_sc(W_pp, Q, w_mod, n) * temp5;
         temp2 += lambda_sc_dot(W_pp, Q, w_mod, n) * temp4;
         temp3 += state.lambda_sc(W_pp, Q, w_mod, n) * temp4;
         }
      val += ((temp1+temp2) * state.w_sc(W_pp,Q, t) + temp3 * w_sc_dot(W_pp,Q, t)) * Model::MomentumGrid().volume_weights[Q];
      }
   
    return val;
}

template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_d_ph_to_xph_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.d_ph_to_xph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_mod(w - div2_floor( W ) + div2_floor( W_ph ));
   int wp_mod(wp + div2_ceil( W ) - div2_ceil( W_ph ));
    
   for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
      {
      dcomplex temp1(0.0,0.0);
      dcomplex temp2(0.0,0.0);
      dcomplex temp3(0.0,0.0);
      for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
         {
         dcomplex temp4(0.0,0.0);
         dcomplex temp5(0.0,0.0);
         for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            {
            temp4 += TUProjections<Model>::Matrix_ph_to_xph()[K][Q][m][mp][n][np] * state.lambda_d(W_ph, Q, wp_mod, np); 
            temp5 += TUProjections<Model>::Matrix_ph_to_xph()[K][Q][m][mp][n][np] * lambda_d_dot(W_ph, Q, wp_mod, np); 
            }
         temp1 += state.lambda_d(W_ph, Q, w_mod, n) * temp5;
         temp2 += lambda_d_dot(W_ph, Q, w_mod, n) * temp4;
         temp3 += state.lambda_d(W_ph, Q, w_mod, n) * temp4;
         }
      val += ((temp1+temp2) * state.w_d(W_ph,Q, t) + temp3 * w_d_dot(W_ph,Q, t)) * Model::MomentumGrid().volume_weights[Q];
      }
   
    return val;
}


template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::nabla_m_ph_to_xph_dot(const state_frg_sbe_t<Model, OtherTypes...> & state, const int W, const int K, const int w, const int m, const int wp, const int mp, const double t ) const{
#ifdef PRECOMPUTE_STATE_PROJECTIONS
    if (m_projected_nablas.are_calculated)
	return (*m_projected_nablas.m_ph_to_xph_ptr)[W][K][w][m][wp][mp];
#endif 

   dcomplex val( 0.0, 0.0 );

   int W_ph = wp - w;
   int w_mod(w - div2_floor( W ) + div2_floor( W_ph ));
   int wp_mod(wp + div2_ceil( W ) - div2_ceil( W_ph ));

   for(int Q = 0; Q<Model::GetRefinedMomentaCount(); ++Q)
      {
      dcomplex temp1(0.0,0.0);
      dcomplex temp2(0.0,0.0);
      dcomplex temp3(0.0,0.0);
      for(int n = 0; n<Model::GetMomentumFormFactorsCount(); ++n)
         {
         dcomplex temp4(0.0,0.0);
         dcomplex temp5(0.0,0.0);
         for(int np = 0; np<Model::GetMomentumFormFactorsCount(); ++np)
            {
            temp4 += TUProjections<Model>::Matrix_ph_to_xph()[K][Q][m][mp][n][np] * state.lambda_m(W_ph, Q, wp_mod, np);
            temp5 += TUProjections<Model>::Matrix_ph_to_xph()[K][Q][m][mp][n][np] * lambda_m_dot(W_ph, Q, wp_mod, np);
            }
         temp1 += state.lambda_m(W_ph, Q, w_mod, n) * temp5;
         temp2 += lambda_m_dot(W_ph, Q, w_mod, n) * temp4;
         temp3 += state.lambda_m(W_ph, Q, w_mod, n) * temp4;
         }
      val += ((temp1+temp2) * state.w_m(W_ph,Q, t) + temp3 * w_m_dot(W_ph,Q, t)) * Model::MomentumGrid().volume_weights[Q];
      }

    return val;
}



//Curently not used
template<typename Model, typename... OtherTypes > 
dcomplex state_frg_sbe_t<Model, OtherTypes...>::eval_f(const double t) const
{
  
   dcomplex val( 0.0, 0.0 );
   // todo: see what to do in case of frequency dependence for free energy
   for( int p_ = 0; p_ < Model::GetFineMomentaCount(); ++p_ ){
       val += std::log(1.0+std::exp(-Beta(t)*(Model::E(p_, 0) - this->m_delta_mu)));
   }
   val *= -2.0/Beta(t)/Model::GetFineMomentaCount();     /*< Normalize freq/momentum summation and times 2 because of spin*/
   std::cout << "Non interacting part of free-energy: "<< val <<std::endl;

   // add correlated part
   val += this->gf_f()[0];
   std::cout << "Full free-energy: "<< val <<std::endl;


   if( std::abs( std::imag( val )) > 0.00001 )
       std::cout << " CAUTION, nonvanishing imaginary part of the free-energy: " << std::imag( val ) << std::endl; 
   
   return val; 
}
