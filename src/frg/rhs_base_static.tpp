#include "cubature.h"
#if defined(PCUBATURE)
#  define cubature pcubature
#else
#  define cubature hcubature
#endif





template <typename Model, typename State>
void rhs_base_t<Model, State>::compute_static_GS_bubbles(gf_bubble_mat_t<Model> *bubble_GS_pp_ptr, gf_bubble_mat_t<Model> *bubble_GS_ph_ptr, const double t){
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_pp_Ptr()->init_batched( (*bubble_GS_pp_ptr), [t]( const idx_bubble_mat_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_diag_bubble_pp_static( idx, t ); } );
    SymmetriesfRGCommon<Model>::IdxEquivClasses_bubblemat_ph_Ptr()->init_batched( (*bubble_GS_ph_ptr), [t]( const idx_bubble_mat_t<Model>& idx ){ return rhs_base_t<Model, State>::eval_diag_bubble_ph_static( idx, t ); } );
}



template <typename Model, typename State>
int rhs_base_t<Model, State>::eval_static_bubble_pp_integrand(unsigned dim, const double *p_coord_data, void *args_data, unsigned fdim, double *fval)
{
    coord_t<Model::dim> p_coord;
    
    for (unsigned d = 0; d < Model::dim; d++)
	p_coord[d] = p_coord_data[d];

    // scale from cube
    p_coord = static_cast<primitive_zone_t<Model::dim>*>(Model::MomentumGrid().ptr)->m_reciprocal_basis_matrix * p_coord;

    auto args = *static_cast<static_bubble_argument_data_t*>(args_data);

    std::complex<double> result = fRGFlowScheme<Model>::static_SG_pp_integrand(p_coord, args.q_coord, args.idx_m, args.idx_n, args.t);

    fval[0] = result.real();
    fval[1] = result.imag();

    return 0;
}

template <typename Model, typename State>
int rhs_base_t<Model, State>::eval_static_bubble_ph_integrand(unsigned dim, const double *p_coord_data, void *args_data, unsigned fdim, double *fval)
{
    coord_t<Model::dim> p_coord;
    
    for (unsigned d = 0; d < Model::dim; d++)
	p_coord[d] = p_coord_data[d];

    // scale from cube
    p_coord = static_cast<primitive_zone_t<Model::dim>*>(Model::MomentumGrid().ptr)->m_reciprocal_basis_matrix * p_coord;

    auto args = *static_cast<static_bubble_argument_data_t*>(args_data);

    std::complex<double> result = fRGFlowScheme<Model>::static_SG_ph_integrand(p_coord, args.q_coord, args.idx_m, args.idx_n, args.t);

    fval[0] = result.real();
    fval[1] = result.imag();

    return 0;
}



template <typename Model, typename State>
MatPatch rhs_base_t<Model, State>::eval_diag_bubble_pp_static( const idx_bubble_mat_t<Model>& idx, const double t)
{
    int W   = idx( IBUBMAT::W );    int w   = idx( IBUBMAT::w ); // ignored
    int m   = idx( IBUBMAT::m );
    int n   = idx( IBUBMAT::n );
    int s4  = idx( IBUBMAT::s1 );
    int s4p = idx( IBUBMAT::s1p); 
    int s3  = idx( IBUBMAT::s2 );
    int s3p = idx( IBUBMAT::s2p );

    double xmin[Model::dim];
    double xmax[Model::dim];


    unsigned fdim = 2; // complex valued function
    unsigned dim = Model::dim; // dimension of domain

    size_t max_evaluations = 100000; // if 0,  unlimited number of function evaluations are allowed

    double required_absolute_error = 1e-5; // ignore absolute error
    double required_relative_error = 1e-3;
    error_norm norm = ERROR_PAIRED; // measure L2 norm of vector of two consequtive numbers (i.e. complex numbers)


    MatPatch val_patch;
    val_patch.setZero(Model::GetRefinedMomentaCount());

    for (unsigned q_idx = 0; q_idx < Model::GetRefinedMomentaCount(); q_idx ++){

	double result[2];
	double error[2];
	
	static_bubble_argument_data_t args;
	args.t = t;
	args.q_coord = Model::MomentumGrid().ptr->m_points[q_idx];
	args.idx_m = m; // one that should be conjugated
	args.idx_n = n;

	for (unsigned d = 0; d < Model::dim; d++){
	    xmin[d] = -0.5;
	    xmax[d] = +0.5;
	}

	// perform adaptive integration
	do{
        cubature(fdim, rhs_base_t<Model, State>::eval_static_bubble_pp_integrand, static_cast<void*>(&args), dim, xmin, xmax, max_evaluations, required_absolute_error, required_relative_error, norm, result, error);
	    for (unsigned d = 0; d < Model::dim; d++){
		unsigned int seed = (d + omp_get_thread_num()) ^ q_idx;
		double rand_shift =  (double)rand_r(&seed)/(RAND_MAX);
		xmin[d] = -0.5 + rand_shift;
		xmax[d] = +0.5 + rand_shift;
	    }
	} while (false);//(isinf(result[0])); // check if result is infinity (there's pole at the boundary, shift the integration range 
	val_patch(q_idx) =  std::complex<double>(result[0], result[1]);
	
    }
    return val_patch * Model::MomentumGrid().ptr->m_volume;
}



// todo: eliminate this piece of code repetition, where pp and ph could be instead given as a (e.g. template) argument
template <typename Model, typename State>
MatPatch rhs_base_t<Model, State>::eval_diag_bubble_ph_static( const idx_bubble_mat_t<Model>& idx, const double t)
{
    int W   = idx( IBUBMAT::W );    int w   = idx( IBUBMAT::w ); // ignored
    int m   = idx( IBUBMAT::m );
    int n   = idx( IBUBMAT::n );
    int s4  = idx( IBUBMAT::s1 );
    int s4p = idx( IBUBMAT::s1p); 
    int s3  = idx( IBUBMAT::s2 );
    int s3p = idx( IBUBMAT::s2p );

    double xmin[Model::dim];
    double xmax[Model::dim];


    unsigned fdim = 2; // complex valued function
    unsigned dim = Model::dim; // dimension of domain

    size_t max_evaluations = 100000; // unlimited number of function evaluations are allowed

    double required_absolute_error = 1e-5; // ignore absolute error
    double required_relative_error = 1e-3;
    error_norm norm = ERROR_PAIRED; // measure L2 norm of vector of two consequtive numbers (i.e. complex numbers)


    MatPatch val_patch;
    val_patch.setZero(Model::GetRefinedMomentaCount());

    for (unsigned q_idx = 0; q_idx < Model::GetRefinedMomentaCount(); q_idx ++){
	//std::cout << q_idx << std::endl;

	double result[2];
	double error[2];
	
	static_bubble_argument_data_t args;
	args.t = t;
	args.q_coord = Model::MomentumGrid().ptr->m_points[q_idx];
	args.idx_m = m; // one that should be conjugated
	args.idx_n = n;

	for (unsigned d = 0; d < Model::dim; d++){
	    xmin[d] = -0.5;
	    xmax[d] = +0.5;
	}

	// perform adaptive integration
	do{
        cubature(fdim, rhs_base_t<Model, State>::eval_static_bubble_ph_integrand, static_cast<void*>(&args), dim, xmin, xmax, max_evaluations, required_absolute_error, required_relative_error, norm, result, error);
	    for (unsigned d = 0; d < Model::dim; d++){
		unsigned int seed = (d + omp_get_thread_num()) ^ q_idx;
		double rand_shift =  (double)rand_r(&seed)/(RAND_MAX);
		xmin[d] = -0.5 + rand_shift;
		xmax[d] = +0.5 + rand_shift;
	    }
	} while (false);//isinf(result[0])); // check if result is infinity (there's pole at the boundary, shift the integration range 
	val_patch(q_idx) =  std::complex<double>(result[0], result[1]);
	
    }

    // include jacobian factor
    return val_patch * Model::MomentumGrid().ptr->m_volume;
}


