
#pragma once

#include <frg/sbe/state_gf_types.h>
#include <frg/interpolators_common.h>

template <typename Model>
class Interpolators1lfRGSBE : public InterpolatorsfRGCommon<Model>
{
 public:
    static void Init();
    
    static void CleanUp();


    // interpolators for w sc/d/m
    static gf_interpolator_t<dcomplex,2>* &w_sc_Ptr(){
	static gf_interpolator_t<dcomplex,2> *interpolator_w_sc_ptr(nullptr);
	return interpolator_w_sc_ptr;
    }
    static gf_interpolator_t<dcomplex,2>* &w_d_Ptr(){
	static gf_interpolator_t<dcomplex,2> *interpolator_w_d_ptr(nullptr);
	return interpolator_w_d_ptr;
    }
    static gf_interpolator_t<dcomplex,2>* &w_m_Ptr(){
	static gf_interpolator_t<dcomplex,2> *interpolator_w_m_ptr(nullptr);
	return interpolator_w_m_ptr;
    }

    // interpolators for lambda sc/d/m
    static gf_interpolator_t<dcomplex,4>* &lambda_sc_Ptr(){
	static gf_interpolator_t<dcomplex,4> *interpolator_lambda_sc_ptr(nullptr);
	return interpolator_lambda_sc_ptr;
    }
    static gf_interpolator_t<dcomplex,4>* &lambda_d_Ptr(){
	static gf_interpolator_t<dcomplex,4> *interpolator_lambda_d_ptr(nullptr);
	return interpolator_lambda_d_ptr;
    }
    static gf_interpolator_t<dcomplex,4>* &lambda_m_Ptr(){
	static gf_interpolator_t<dcomplex,4> *interpolator_lambda_m_ptr(nullptr);
	return interpolator_lambda_m_ptr;
    }

    // interpolators for M sc/d/m
    static gf_interpolator_t<dcomplex,6>* &M_sc_Ptr(){
	static gf_interpolator_t<dcomplex,6> *interpolator_M_sc_ptr(nullptr);
	return interpolator_M_sc_ptr;
    }
    static gf_interpolator_t<dcomplex,6>* &M_d_Ptr(){
	static gf_interpolator_t<dcomplex,6> *interpolator_M_d_ptr(nullptr);
	return interpolator_M_d_ptr;
    }
    static gf_interpolator_t<dcomplex,6>* &M_m_Ptr(){
	static gf_interpolator_t<dcomplex,6> *interpolator_M_m_ptr(nullptr);
	return interpolator_M_m_ptr;
    }
};
