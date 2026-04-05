#include <frg/sbe/interpolators.h>
#include <frg/sbe/symmetries.h>
#include <models/concrete_available_models.h>
#include <algorithm>

template <typename Model>
void Interpolators1lfRGSBE<Model>::Init()
{
    InterpolatorsfRGCommon<Model>::Init();

    std::cout << " Constructing screened-interaction (w) interpolators... " << std::endl;
    w_sc_Ptr() = new gf_interpolator_t<dcomplex,2>(*Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_sc_Ptr());
    w_d_Ptr() = new gf_interpolator_t<dcomplex,2>(*Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_d_Ptr());
    w_m_Ptr() = new gf_interpolator_t<dcomplex,2>(*Symmetries1lfRGSBE<Model>::IdxEquivClasses_w_m_Ptr());


    std::cout << " Constructing hedin vertex (lambda) interpolators... " << std::endl;
    lambda_sc_Ptr() = new gf_interpolator_t<dcomplex,4>(*Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_sc_Ptr());
    lambda_d_Ptr() = new gf_interpolator_t<dcomplex,4>(*Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_d_Ptr());
    lambda_m_Ptr() = new gf_interpolator_t<dcomplex,4>(*Symmetries1lfRGSBE<Model>::IdxEquivClasses_lambda_m_Ptr());


    std::cout << " Constructing SBE rest function (M) interpolators... " << std::endl;
    M_sc_Ptr() = new gf_interpolator_t<dcomplex,6>(*Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_sc_Ptr());
    M_d_Ptr() = new gf_interpolator_t<dcomplex,6>(*Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_d_Ptr());
    M_m_Ptr() = new gf_interpolator_t<dcomplex,6>(*Symmetries1lfRGSBE<Model>::IdxEquivClasses_M_m_Ptr());
}




template <typename Model>
void Interpolators1lfRGSBE<Model>::CleanUp()
{
    InterpolatorsfRGCommon<Model>::CleanUp();

    delete w_sc_Ptr();
    delete w_d_Ptr();
    delete w_m_Ptr();

    delete lambda_sc_Ptr();
    delete lambda_d_Ptr();
    delete lambda_m_Ptr();

    delete M_sc_Ptr();
    delete M_d_Ptr();
    delete M_m_Ptr();
}



// instantiate class for the particular models
#define Instantiate(MODEL) template class Interpolators1lfRGSBE<MODEL>;
WITH_MODELS(Instantiate)


