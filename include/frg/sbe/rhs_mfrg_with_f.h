#pragma once

#include <frg/sbe/rhs_mfrg.h>

template <typename Model, typename state_t >
class rhs_sbe_mfrg_with_F_t : public rhs_sbe_mfrg_t<Model, state_t>
{
public:

    rhs_sbe_mfrg_with_F_t(){}

    using rhs_base_t = ::rhs_base_t<Model, state_t >;

    void operator() ( const state_t& old_state, state_t &dstate_over_dt, const double t );
  
protected:

  void rhs_F_m(gf_f_t<Model> *new_F_ptr, const state_t &old_state, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, const double t);

private:

    // ------- Left corrections for vertex corrections -----------

    // free-energy
    dcomplex eval_F_m( const idx_f_t<Model>& idx, const state_t& state, const gf_1p_mat_t<Model>& Gvec, const gf_1p_mat_t<Model>& Svec, double t);

};

// include implementation
#include <../src/frg/sbe/rhs_mfrg_with_f.tpp> // logic for calculation of right hand side (the free energy)
#include <../src/frg/sbe/rhs_mfrg_with_f_eval.tpp> // (todo: not needed)

