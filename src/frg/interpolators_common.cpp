#include <frg/interpolators_common.h>
#include <frg/symmetries_common.h>
#include <models/concrete_available_models.h>
#include <algorithm>

template <typename Model>
void InterpolatorsfRGCommon<Model>::Init()
{
    // todo: if had been allocated before, should delete first..
   
    std::cout << " Constructing self-energy (sig) interpolators... " << std::endl;
    sig_Ptr() = new gf_interpolator_t<dcomplex,4>(*SymmetriesfRGCommon<Model>::IdxEquivClasses_sig_Ptr());
 
    std::cout << " Constructing pp-bubble interpolators... " << std::endl; 
    bubblemat_pp_Ptr() = new gf_interpolator_t<MatPatch,8>(*SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_pp_Ptr());

    std::cout << " Constructing ph-bubble interpolators... " << std::endl;
    bubblemat_ph_Ptr() = new gf_interpolator_t<MatPatch,8>(*SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_ph_Ptr());
   


    std::cout << " Constructing self-energy (sig_kmat) interpolators... " << std::endl;
    sigkmat_Ptr() = new gf_interpolator_t<MatPatch,3>(*SymmetriesfRGCommon<Model>::IdxEquivClasses_sigkmat_Ptr());

    // todo: add Gvec's interpolators

    // todo: add Gvec real interpolators
    std::cout << " Constructing real Green's function interpolators... " << std::endl;
    Gvec_real_Ptr() = new gf_interpolator_t<MatReal,3>(*SymmetriesfRGCommon<Model>::IdxEquivClasses_Gvec_real_Ptr());

    std::cout << " Constructing projection matrix interpolators... " << std::endl;
    projection_matrix_Ptr() = new gf_interpolator_t<dcomplex, 6>(*SymmetriesfRGCommon<Model>::IdxEquivClasses_projection_matrix_Ptr());
}

template <typename Model>
void InterpolatorsfRGCommon<Model>::CleanUp()
{
    delete bubblemat_pp_Ptr();

    delete bubblemat_ph_Ptr();

    delete sig_Ptr();

    delete sigkmat_Ptr();

    delete projection_matrix_Ptr();
}


// instantiate class for the particular models
#define Instantiate(MODEL) template class InterpolatorsfRGCommon<MODEL>;
WITH_MODELS(Instantiate)


