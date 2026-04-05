#pragma once

#include <frg/sbe/state_gf_types.h>
#include <gf_interpolator.h>
#include <models/concrete_available_models.h>

template <typename Model>
class InterpolatorsfRGCommon
{
 public:
    static void Init();
    
    static void CleanUp();

    // interpolator for bubbles
    static gf_interpolator_t<MatPatch, 8>* &bubblemat_pp_Ptr(){
	static gf_interpolator_t<MatPatch, 8> *interpolator_bubblemat_pp_ptr(nullptr);
	return interpolator_bubblemat_pp_ptr;
    }

    static gf_interpolator_t<MatPatch, 8>* &bubblemat_ph_Ptr(){
	static gf_interpolator_t<MatPatch, 8> *interpolator_bubblemat_ph_ptr(nullptr);
	return interpolator_bubblemat_ph_ptr;
    }

    // interpolator for sig
    static gf_interpolator_t<dcomplex, 4>* &sig_Ptr(){
	static gf_interpolator_t<dcomplex, 4> *interpolator_sig_ptr(nullptr);
	return interpolator_sig_ptr;
    }

    // interpolator for the MatPatch-valued sig
    static gf_interpolator_t<MatPatch,3>* &sigkmat_Ptr(){
	static gf_interpolator_t<MatPatch,3> *interpolator_sigkmat_ptr(nullptr);
	return interpolator_sigkmat_ptr;
    }

    // interpolator for the Green function in real space
    static gf_interpolator_t<MatReal,3>* &Gvec_real_Ptr(){
	static gf_interpolator_t<MatReal,3> *interpolator_Gvec_real_ptr(nullptr);
	return interpolator_Gvec_real_ptr;
    }

    // interpolator for the projection matrices
    static gf_interpolator_t<dcomplex, 6>* &projection_matrix_Ptr(){
	static gf_interpolator_t<dcomplex, 6> *interpolator_projection_matrix_ptr(nullptr);
	return interpolator_projection_matrix_ptr;
    }

};

