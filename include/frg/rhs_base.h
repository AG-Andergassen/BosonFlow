#pragma once

#include <frg/state_base.h>



template <typename Model, typename State>
class rhs_base_t
{
    virtual void operator() (const State &old_state, State &dstate_over_dt, const double t) = 0;

 protected:
    void compute_Gvec(gf_1p_mat_t<Model> *Gvec_ptr, const State &old_state); // Gvec without a scale
    void compute_Gvec(gf_1p_mat_t<Model> *Gvec_ptr, const State &old_state, const double t);
    void compute_Svec(gf_1p_mat_t<Model> *Svec_ptr, const State &old_state, const double t);
    void compute_Gvec_real(gf_1p_mat_real_t<Model> *Gvec_real_ptr, const gf_1p_mat_t<Model>& Gvec);
    void compute_Svec_real(gf_1p_mat_real_t<Model> *Svec_real_ptr, const gf_1p_mat_t<Model>& Svec);
    void compute_GS_bubbles(gf_bubble_mat_t<Model> *bubble_GS_pp_ptr, gf_bubble_mat_t<Model> *bubble_GS_ph_ptr, const gf_1p_mat_real_t<Model> &Gvec_real, const gf_1p_mat_real_t<Model> &Svec_real);
    void compute_GS_bubbles_trad(gf_bubble_mat_t<Model> *bubble_GS_pp_ptr, gf_bubble_mat_t<Model> *bubble_GS_ph_ptr, const gf_1p_mat_t<Model> &Gvec, const gf_1p_mat_t<Model> &Svec_real);
    void compute_GG_bubbles(gf_bubble_mat_t<Model> *bubble_GG_pp_ptr, gf_bubble_mat_t<Model> *bubble_GG_ph_ptr, const gf_1p_mat_real_t<Model> &Gvec_real);

    void compute_static_GS_bubbles(gf_bubble_mat_t<Model> *bubble_GS_pp_ptr, gf_bubble_mat_t<Model> *bubble_GS_ph_ptr, const double t);

    struct static_bubble_argument_data_t
    {
	double t; // scale
	unsigned idx_m, idx_n; // form factor indices. The m form factor is the one conjugated.
	coord_t<Model::dim> q_coord; // bosonic momentum  
	// TODO: self energy
    };
    

 public:
    
    static MatReal eval_Gvec_real( const idx_1p_mat_real_t<Model>& idx, const gf_1p_mat_t<Model>& Gvec);

    static MatPatch eval_diag_bubble_pp( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_real_t<Model>& Svec_real, const gf_1p_mat_real_t<Model>& Gvec_real);

    static MatPatch eval_diag_bubble_ph( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_real_t<Model>& Svec_real, const gf_1p_mat_real_t<Model>& Gvec_real);

    static MatPatch eval_diag_bubble_GG_pp( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_real_t<Model>& Gvec_real);

    static MatPatch eval_diag_bubble_GG_ph( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_real_t<Model>& Gvec_real);
    
    static MatPatch eval_diag_bubble_pp_trad( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_t<Model>& G, const gf_1p_mat_t<Model>& S );

    static MatPatch eval_diag_bubble_ph_trad( const idx_bubble_mat_t<Model>& idx, const gf_1p_mat_t<Model>& G, const gf_1p_mat_t<Model>& S );

    static MatPatch eval_diag_bubble_pp_static( const idx_bubble_mat_t<Model>& idx, const double t);

    static MatPatch eval_diag_bubble_ph_static( const idx_bubble_mat_t<Model>& idx, const double t);

    // for post-processing
    static MatPatch eval_Sig_Lam2PI_kMat( const idx_Sig_kMat_t<Model>& idx, const gf_1p_mat_real_t<Model>&  Gvec_real, const gf_1p_mat_t<Model>& Gvec, const double t );


    static int	
	eval_static_bubble_pp_integrand(
					unsigned p_coord_dim, // (Model::dim) 
					const double *p_coord_data, // the argument
					void *args_data, // to be cast to type static_bubble_argument_data_t
					unsigned fdim, // always 2, since static bubble is complex valued 
					double *fval); // output (real and imaginary part)
    
    static int eval_static_bubble_ph_integrand(unsigned p_coord_dim, const double *p_coord_data, void *args_data, unsigned fdim, double *fval); 

    static void print_current_scale(const double t);
};

// implementation of template
#include "../src/frg/rhs_base.tpp"
#include "../src/frg/rhs_base_static.tpp"


