#include <models/hubbard.h>
#include <bond_form_factor_container.h>
#include <omp.h>

template<unsigned dim_>
void Hubbard<dim_>::Init(basis_t<dim_> basis, 
			 std::vector<std::vector<matrix_t<dim_> > > point_group_matrices, 
			 std::map<std::string, coord_t<dim_> > special_points_coords, 
			 std::map<std::string, std::vector<coord_t<dim_> > > special_paths_coords, 
			 std::map<std::string, double>  refine_at_points, 
			 std::vector< delta_function_vector_t<dim_> > extra_form_factors)
{
    std::cout << "Initialising lattice..." << std::endl;
    InitLattices(basis);

    std::cout << "Applying momentum space refinement(s)..." << std::endl;
    RefineLattices(special_points_coords, special_paths_coords, refine_at_points);

    std::cout << "Defining special points and paths..." << std::endl;
    InitSpecialPointsAndPaths(special_points_coords, special_paths_coords);

    std::cout << "Initialising lattice point group..." << std::endl;
    InitLatticesPointGroup(point_group_matrices);

    InitPrecalculatedTrigonometricMaps();

    std::cout << "Initialising momentum idx space symmetry maps..." << std::endl;
    double mom_symmetry_maps_error = InitMomentumSymmetryMaps();

    std::cout << "Building Coarse to Fine grid momentum idx map..." << std::endl;
    double mom_coarse_to_fine_map_error = InitCoarseToFineMomentumIdxMap();

    std::cout << "Building Fine to Coarse grid momentum idx map..." << std::endl;
    double mom_fine_to_coarse_map_error = InitFineToCoarseMomentumIdxMap();

#ifdef PRECOMPUTE_MOMENTUM_ARITHMETIC
    std::cout << "Initialising momentum idx arithmetic maps..." << std::endl;
    double mom_arithmetic_maps_error = InitPrecalculatedMomentumIdxArithmeticMaps();
#endif
    
    std::cout << "Initialising form factors..." << std::endl;
    InitFormFactors(extra_form_factors);

    std::cout << "Initialising form factors fourier transform weights..." << std::endl;

    InitFormFactorsFourierTransformWeights();

    std::cout << "Initialising form factors idx symmetry maps (from momentum space symmetry maps)..." << std::endl;
    InitFormFactorsSymmetryMaps();

    static int thread_count;
    
    #pragma omp parallel
    {
        #pragma omp single
        thread_count = omp_get_num_threads();
    }

    FFTWs().resize(thread_count);
    
    std::cout << "Number of threads available: " << thread_count << std::endl;
    
    
    int mom_count[dim_];
    for (unsigned i = 0; i < dim_; i ++)
	mom_count[i] = K_DIM * P_IN_K; 

    // Todo: check fftw_plan_many

    // fourier transform on the fine momentum grid
    // since that code runs in parallel, every thread will needs its own copy
    for (unsigned thread_num = 0; thread_num < thread_count; thread_num++){
	FFTW(thread_num).input = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * GetFineMomentaCount());
	FFTW(thread_num).output = (fftw_complex*) fftw_malloc(sizeof(fftw_complex) * GetFineMomentaCount() );
	FFTW(thread_num).forward_plan = fftw_plan_dft(dim_, mom_count, FFTW(thread_num).input, FFTW(thread_num).output, FFTW_FORWARD, FFTW_MEASURE); 
	FFTW(thread_num).backward_plan = fftw_plan_dft(dim_, mom_count, FFTW(thread_num).input, FFTW(thread_num).output, FFTW_BACKWARD, FFTW_MEASURE);
    }

    // todo: add recommendation that mom_count is a power of two.    
}

template<unsigned dim_>
void Hubbard<dim_>::InitLattices(basis_t<dim_> basis)
{
    // make square lattice with 1 nearest neighbohrs in real space. For the projection matrices to be constructed correctly, this should go up be at least 2 times the maximal range of the form factors. We just choose 10 for safety
    RealLattice().ptr = new real_space_lattice_t<dim_>(basis,10);

    // make brillouin zone with K_DIM points along each BZ edge
    MomentumGrid().ptr = new primitive_zone_t<dim_>(RealLattice().ptr->make_reciprocal_lattice_basis(), K_DIM);

    MomentaCountUnrefined() = MomentumGrid().ptr->m_points.size();

    MomentumGrid().volume_weights = dynamic_cast<primitive_zone_t<dim_>*>(MomentumGrid().ptr)->m_volume_weights;
    

    // this is another BZ that is more fine (with P_IN_K * K_DIM points along each BZ edge). Desired for a more careful evaluation of the fourier transform in the computation of the bare bubbles (using the convolution theorem)
    FineMomentumGrid().ptr = new primitive_zone_t<dim_>(RealLattice().ptr->make_reciprocal_lattice_basis(), P_IN_K * K_DIM);

    FineMomentumGrid().volume_weights = dynamic_cast<primitive_zone_t<dim_> *>(FineMomentumGrid().ptr)->m_volume_weights;
}

template<unsigned dim_>
void Hubbard<dim_>::InitLatticesPointGroup(std::vector<std::vector<matrix_t<dim_> > > point_group_matrices)
{
    RealLattice().point_group_ptr = new lattice_point_group_t<dim_>(point_group_matrices); 
}

template<unsigned dim_>
double Hubbard<dim_>::InitMomentumSymmetryMaps()
{
    double err;
    
    std::tie(MomentumGrid().point_group_in_momentum_idx_space, err) = RealLattice().point_group_ptr->make_group_in_momentum_idx_space(*(MomentumGrid().ptr), 1e-10);
    std::cout << "The maximum error in applying the symmetry group to the default momentum grid is: " << err << std::endl;

    std::tie(FineMomentumGrid().point_group_in_momentum_idx_space, err) = RealLattice().point_group_ptr->make_group_in_momentum_idx_space(*(FineMomentumGrid().ptr), 1e-10);
    std::cout << "The maximum error in applying the symmetry group to the fine (fourier transform) momentum grid is: " << err << std::endl;

    //        std::vector<std::map<shell_ff_idx_t, std::tuple<ff_idx_t, std::string> > > pt_grp = square_bond_ffs.make_point_group_in_shells_form_factor_idx_space(sq_point_group);

    return err;
}

template<unsigned dim_>
void Hubbard<dim_>::InitFormFactorsSymmetryMaps()
{
    if (FormFactors().container_ptr == nullptr)
	throw "Form factors were not initialised! They should be initialised first before the form factor symmetry maps can be.";

    // get point group in form-factor index space
    std::vector<std::map<shell_ff_idx_t, std::tuple<ff_idx_t, std::string> > > pt_grp = FormFactors().container_ptr->make_point_group_in_shells_form_factor_idx_space(MomentumGrid().point_group_in_momentum_idx_space, MomentumGrid().form_factors_in_momentum_idx_space);
    
    for (std::map<shell_ff_idx_t, std::tuple<ff_idx_t, std::string> > &g : pt_grp){
	std::vector<std::tuple<unsigned, bool> > g_in_ff_idx_space;	
	for (unsigned ff_flat_idx = 0; ff_flat_idx < FormFactors().flat_index_to_tuple_index_map.size(); ff_flat_idx ++){
	    auto tuple_idx = FormFactors().flat_index_to_tuple_index_map[ff_flat_idx];
	    auto [bond_idx, symm_operation] = g[tuple_idx];
	    
	    auto ff_prime_flat_idx = FormFactors().tuple_index_to_flat_index_map[std::make_tuple(std::get<0>(tuple_idx), bond_idx)];

	    // g_in_ff_idx_space[3][1] == (1, true)
	    g_in_ff_idx_space.push_back(std::make_tuple(ff_prime_flat_idx, symm_operation == "minus"));
	}
	FormFactors().point_group_in_form_factor_idx_space.push_back(g_in_ff_idx_space);
    }
}

template<unsigned dim_>
void Hubbard<dim_>::InitFormFactorsFourierTransformWeights(){
    if (FormFactors().container_ptr == nullptr)
	throw "Form factors were not initialised! They should be initialised first before the form factor weights can be.";

    for (coord_t<dim_> R: RealLattice().ptr->m_points){
	std::vector<std::vector<dcomplex> > fxf;
	for (unsigned m = 0; m < GetMomentumFormFactorsCount(); m++){
	    std::vector<dcomplex> fxf_of_m;
	    for (unsigned n  = 0; n < GetMomentumFormFactorsCount(); n++){
		dcomplex weight(0.0, 0.0);
		for (coord_t<dim_> Rp : RealLattice().ptr->m_points){
		    weight += std::conj(GetFormFactor(m).evaluate(Rp - R)) * 
			GetFormFactor(n).evaluate(Rp);
		}
		fxf_of_m.push_back(weight);
	    }
	    fxf.push_back(fxf_of_m);
	}
	FormFactors().fourier_transform_weights.push_back(fxf);
    }
    // remark: the +3 was previously +1.. todo: understand this.Ans reason is that one additional shell only contains the diagonal neighbor rather than the one further away which is next next nearest neighbor
    FormFactors().fourier_transform_weights_r_support_indices = RealLattice().ptr->get_subset_indices( (FormFactors().container_ptr->m_shells_count - 1) * 2 + 3); // rest would be zero

    FormFactors().real_space_to_real_fourier_transform_idx_map.resize(RealLattice().ptr->m_points.size());
    for (auto r_idx : FormFactors().fourier_transform_weights_r_support_indices){
    auto vec = (2.0 * M_PI/(K_DIM*P_IN_K)) * RealLattice().ptr->m_basis_matrix_inv.transpose() * RealLattice().ptr->m_basis_matrix_inv * RealLattice().ptr->m_points[r_idx];
	FormFactors().real_space_to_real_fourier_transform_idx_map[r_idx] =  std::get<0>(dynamic_cast<primitive_zone_t<dim_>*>(FineMomentumGrid().ptr)->get_idx_from_coord(vec));
    }

    // TODO: move to appropriate separate function

    FormFactors().ff_convolution.resize(GetMomentumFormFactorsCount());
    for (unsigned m = 0; m < GetMomentumFormFactorsCount(); m++){
      FormFactors().ff_convolution[m].resize(GetMomentumFormFactorsCount());
      for (unsigned n = 0; n < GetMomentumFormFactorsCount(); n ++){
	FormFactors().ff_convolution[m][n].resize(GetCoarseMomentaCount());
	for (unsigned k_idx; k_idx < GetCoarseMomentaCount(); k_idx++){
	  FormFactors().ff_convolution[m][n][k_idx] = 0.0;
	  for (coord_t<dim_> R : RealLattice().ptr->m_points){
	    const auto exp_ikR = std::exp( std::complex<double>(0.0, MomentumGrid().ptr->m_points[k_idx].dot(R) ) );
	    FormFactors().ff_convolution[m][n][k_idx] += exp_ikR * std::conj(GetFormFactor(m).evaluate(R)) * GetFormFactor(n).evaluate(R);
	  }
	}
      }
    }
}

template<unsigned dim_>
void Hubbard<dim_>::InitPrecalculatedTrigonometricMaps()
{
    MomentumGrid().exp_ik_idx_map = MomentumGrid().ptr->precalculate_exp_ik();
    FineMomentumGrid().exp_ik_idx_map = FineMomentumGrid().ptr->precalculate_exp_ik();
    
    // needed for the fourier transform for computation of bubble
    MomentumGrid().exp_ikr_idx_map = MomentumGrid().ptr->precalculate_exp_ikr(RealLattice().ptr->m_points);
    FineMomentumGrid().exp_ikr_idx_map = FineMomentumGrid().ptr->precalculate_exp_ikr(RealLattice().ptr->m_points);
}

template<unsigned dim_>
double Hubbard<dim_>::InitPrecalculatedMomentumIdxArithmeticMaps()
{
    
    double err1; 
    std::tie(MomentumGrid().negative_momentum_idx_map, err1) = MomentumGrid().ptr->precalculate_negative_indices(GetCoarseMomentaCount());
    std::cout << "The maximum error in precalculating the negative-momentum indices on the (coarse) momentum grid is: " << err1 << std::endl;

    double err2;

    std::tie(FineMomentumGrid().negative_momentum_idx_map, err2) = FineMomentumGrid().ptr->precalculate_negative_indices(GetFineMomentaCount());
    
    std::cout << "The maximum error in precalculating the negative-momentum indices on the fine/fourier-transform grid is: " << err2 << std::endl;

    double err3;
    std::tie(MomentumGrid().momentum_plus_momentum_idx_map, err3) = MomentumGrid().ptr->precalculate_sum_of_indices(GetCoarseMomentaCount());
    std::cout << "The maximum error in precalculating the momentum sums indices on the (coarse) momentum grid is: " << err3 << std::endl;

    double err4;
    std::tie(FineMomentumGrid().momentum_plus_momentum_idx_map, err4) = FineMomentumGrid().ptr->precalculate_sum_of_indices(GetFineMomentaCount());
    std::cout << "The maximum error in precalculating the momentum sums indices on the fine/fourier-transform grid is: " << err4 << std::endl;

    // return the biggest error
    return std::max(std::max(std::max(err1, err2), err3), err4);
}

template<unsigned dim_>
void Hubbard<dim_>::InitSpecialPointsAndPaths(std::map<std::string, coord_t<dim_> > special_points_coords, std::map<std::string, std::vector<coord_t<dim_> > > special_paths_coords)
{
    // special points
    for (auto itr = special_points_coords.begin(); itr != special_points_coords.end(); itr ++){
	std::string pt_name = itr->first;
	coord_t<dim_> pt_coord = itr->second;

	SpecialPointsAndPaths().points_names.push_back(pt_name);
	SpecialPointsAndPaths().points.push_back(std::get<0>(MomentumGrid().ptr->get_idx_from_coord(pt_coord)));
    }

    // special paths
    for (auto itr = special_paths_coords.begin(); itr != special_paths_coords.end(); itr ++){
	std::string path_name = itr->first;
	std::vector<coord_t<dim_> > path_pts_coords = itr->second;
	std::vector<unsigned> path_idxes = dynamic_cast<primitive_zone_t<dim_>*>(MomentumGrid().ptr)->get_path_indices(path_pts_coords, true); // "true" for we consider the is closed

	SpecialPointsAndPaths().paths_names.push_back(path_name);
	SpecialPointsAndPaths().paths.push_back(path_idxes);
    }
}

template<unsigned dim_>	
void Hubbard<dim_>::RefineLattices(std::map<std::string, coord_t<dim_> > special_points_coords, std::map<std::string, std::vector<coord_t<dim_> > > special_paths_coords, std::map<std::string, double> refine_at_points)
{
    // refine lattice
    for (auto itr = refine_at_points.begin(); itr != refine_at_points.end(); itr++){
	std::cout << "Refining at " << itr->first << "..." << std::endl;
	
	coord_t<dim_> refine_coord = special_points_coords[itr->first];
	double refine_percent = itr->second;
	
	dynamic_cast<primitive_zone_t<dim_>*>(MomentumGrid().ptr)->refine_from_finer_grid(FineMomentumGrid().ptr, refine_coord, refine_percent);
    }
    
    /*
    dynamic_cast<primitive_zone_t<dim_>*>(MomentumGrid().ptr)->refine_from_finer_grid_along_path(FineMomentumGrid().ptr, { {0, M_PI}, {M_PI, 0}, {0, -M_PI}, {-M_PI, 0} }, 0.007, true);
    */

    /*for (unsigned p_idx = 0; p_idx < K_DIM*K_DIM; p_idx ++)
      MomentumGrid().ptr->m_points[p_idx] = RealLattice().ptr->m_basis_matrix.transpose() * MomentumGrid().ptr->m_points[p_idx];*/

    MomentumGrid().volume_weights = dynamic_cast<primitive_zone_t<dim_>*>(MomentumGrid().ptr)->m_volume_weights;

    std::cout << MomentumGrid().ptr->make_debug_grid_python_list("array") << std::endl;
    double sum = 0;
    for (auto w: MomentumGrid().volume_weights)
      sum += w;
    std::cout << "Weights sum up to: " << sum << ", and you're happy if this is 1.0" << std::endl;
}


template<unsigned dim_>
void Hubbard<dim_>::InitFormFactors(std::vector< delta_function_vector_t<dim_> > extra_form_factors) // todo: implement extra form factor dependence
{
    
    const unsigned form_factor_shells = static_cast<unsigned int>(FORMFACTOR_SHELL_COUNT);

	// make bond form factors.
	FormFactors().container_ptr = new bond_form_factor_container_t<dim_>(*RealLattice().ptr, form_factor_shells);
	
	std::cout << "Form factor count: " << FORMFACTOR_SHELL_COUNT << std::endl;
	if (FORMFACTOR_SHELL_COUNT == 1.3){ // p-wave form factors
	    FormFactors().container_ptr->m_form_factors_shells.push_back({});
	    //FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }, std::vector<std::complex<double> >{-std::complex<double>(0, 1)/2.0, std::complex<double>(0, 1)/2.0, -std::complex<double>(0, 1)/2.0, std::complex<double>(0, 1)/2.0} ));
	    //FormFactors().container_ptr->m_shells_count++;
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {-1, 0}, {+1, 0}, {0, -1}, {0, +1} }, std::vector<std::complex<double> >{-std::sqrt(2.0)*std::complex<double>(0, 1)/2.0, std::sqrt(2.0)*std::complex<double>(0, 1)/2.0, 0.0, 0.0} )); // sin k_x form factor
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {-1, 0}, {+1, 0}, {0, -1}, {0, +1} }, std::vector<std::complex<double> >{0.0, 0.0, std::sqrt(2.0)*-std::complex<double>(0, 1)/2.0, std::sqrt(2.0)*std::complex<double>(0, 1)/2.0} )); // sin k_y form factor
	    FormFactors().container_ptr->m_shells_count++;
	}

	if (FORMFACTOR_SHELL_COUNT == 1.5){ // d-wave form factors
	    FormFactors().container_ptr->m_form_factors_shells.push_back({});
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }, {1.0/2, 1.0/2, -1.0/2, -1.0/2} ));
	    FormFactors().container_ptr->m_shells_count++;
	}
	if (FORMFACTOR_SHELL_COUNT == 1.7){
	    FormFactors().container_ptr->m_form_factors_shells.push_back({});
	    //FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }, std::vector<std::complex<double> >{-std::complex<double>(0, 1)/2.0, std::complex<double>(0, 1)/2.0, -std::complex<double>(0, 1)/2.0, std::complex<double>(0, 1)/2.0} ));
	    //FormFactors().container_ptr->m_shells_count++;
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {-1, 0}, {+1, 0}, {0, -1}, {0, +1} }, std::vector<std::complex<double> >{-std::sqrt(2.0)*std::complex<double>(0, 1)/2.0, std::sqrt(2.0)*std::complex<double>(0, 1)/2.0, 0.0, 0.0} )); // (1/sqrt2) * sin k_x form factor
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {-1, 0}, {+1, 0}, {0, -1}, {0, +1} }, std::vector<std::complex<double> >{0.0, 0.0, std::sqrt(2.0)*-std::complex<double>(0, 1)/2.0, std::sqrt(2.0)*std::complex<double>(0, 1)/2.0} )); // (1/sqrt2) * sin k_y form factor
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }, {1.0/2, 1.0/2, -1.0/2, -1.0/2} )); // cos k_x - cos k_y
	    FormFactors().container_ptr->m_shells_count++;
	}

	if (FORMFACTOR_SHELL_COUNT == 1.99){
	  std::cout << "yolo" << std::endl;
	    FormFactors().container_ptr->m_form_factors_shells.push_back({});
	    //FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }, std::vector<std::complex<double> >{-std::complex<double>(0, 1)/2.0, std::complex<double>(0, 1)/2.0, -std::complex<double>(0, 1)/2.0, std::complex<double>(0, 1)/2.0} ));
	    //FormFactors().container_ptr->m_shells_count++;
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {-1, 0}, {+1, 0}, {0, -1}, {0, +1} }, std::vector<std::complex<double> >{-std::sqrt(2.0)*std::complex<double>(0, 1)/2.0, std::sqrt(2.0)*std::complex<double>(0, 1)/2.0, 0.0, 0.0} )); // (1/sqrt2) * sin k_x form factor
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {-1, 0}, {+1, 0}, {0, -1}, {0, +1} }, std::vector<std::complex<double> >{0.0, 0.0, std::sqrt(2.0)*-std::complex<double>(0, 1)/2.0, std::sqrt(2.0)*std::complex<double>(0, 1)/2.0} )); // (1/sqrt2) * sin k_y form factor
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }, {1.0/2, 1.0/2, -1.0/2, -1.0/2} )); // cos k_x - cos k_y
	    FormFactors().container_ptr->m_form_factors_shells[1].push_back(delta_function_vector_t<dim_>({ {1, 0}, {-1, 0}, {0, 1}, {0, -1} }, {1.0/2, 1.0/2, 1.0/2, 1.0/2} )); // cos k_x + cos k_y
	    FormFactors().container_ptr->m_shells_count++;
	}
    
     
    auto ff_tuple_indices = FormFactors().container_ptr->get_indices();
    
    FormFactors().flat_index_to_tuple_index_map.reserve(ff_tuple_indices.size() );

    for (unsigned ff_flat_idx = 0; ff_flat_idx < ff_tuple_indices.size(); ff_flat_idx ++){
	    std::cout << ff_flat_idx << std::endl;
	    FormFactors().tuple_index_to_flat_index_map[ff_tuple_indices[ff_flat_idx]] = ff_flat_idx;
	    FormFactors().flat_index_to_tuple_index_map.push_back(ff_tuple_indices[ff_flat_idx]);
    }

    MomentumGrid().form_factors_in_momentum_idx_space = FormFactors().container_ptr->make_form_factors_in_momentum_idx_space(*MomentumGrid().ptr);
    
    FineMomentumGrid().form_factors_in_momentum_idx_space = FormFactors().container_ptr->make_form_factors_in_momentum_idx_space(*FineMomentumGrid().ptr);
}

template<unsigned dim_>
double Hubbard<dim_>::InitCoarseToFineMomentumIdxMap()
{
    std::vector<unsigned> coarse_to_fine_mom_idx_map;

    if (MomentumGrid().ptr == nullptr || FineMomentumGrid().ptr == nullptr)
	throw "coarse and(or) fine grid were not initialised! They should be initialised first before the coarse to fine momentum map (and the map in the other direction) can be constructed.";

    auto &coarse_mesh = MomentumGrid().ptr->m_points;

    double max_error = -1;

    for (auto coarse_mom_coord : coarse_mesh){
	auto [fine_idx, error] = FineMomentumGrid().ptr->get_idx_from_coord(coarse_mom_coord);
	coarse_to_fine_mom_idx_map.push_back(fine_idx);
	max_error = std::max(error, max_error);
    }

    std::cout << "Max deviation in calculated coarse to fine momentum grid map: " << max_error << std::endl;

    CoarseToFineIdxMap() = coarse_to_fine_mom_idx_map;

    return max_error;
}

template<unsigned dim_>
double Hubbard<dim_>::InitFineToCoarseMomentumIdxMap()
{
    std::vector<unsigned> fine_to_coarse_mom_idx_map;

    if (MomentumGrid().ptr == nullptr || FineMomentumGrid().ptr == nullptr)
	throw "coarse and(or) fine grid were not initialised! They should be initialised first before the coarse to fine momentum map (and the map in the other direction) can be constructed.";

    auto &fine_mesh = FineMomentumGrid().ptr->m_points;

    double max_error = -1;

    for (auto fine_mom_coord : fine_mesh){
	auto [coarse_idx, error] = dynamic_cast<primitive_zone_t<dim_>*>(MomentumGrid().ptr)->get_idx_from_coord_unrefined(fine_mom_coord);
	fine_to_coarse_mom_idx_map.push_back(coarse_idx);
	max_error = std::max(error, max_error);
    }

    std::cout << "Max deviation in calculated fine to coarse momentum grid map: " << max_error << std::endl;

    FineToCoarseIdxMap() = fine_to_coarse_mom_idx_map;

    return max_error;
}

template<unsigned dim_>
coord_t<dim_> Hubbard<dim_>::GetMomentumCoordFromCoarseIdx(const int idx_k)
{   
    return MomentumGrid().ptr->m_points[idx_k];
}

template<unsigned dim_>
coord_t<dim_> Hubbard<dim_>::GetMomentumCoordFromFineIdx(const int idx_p)
{
    return FineMomentumGrid().ptr->m_points[idx_p];
}

template<unsigned dim_>
delta_function_vector_t<dim_> Hubbard<dim_>::GetFormFactor( const unsigned n )
{
    return FormFactors().container_ptr->get(FormFactors().flat_index_to_tuple_index_map[n]);
}

template<unsigned dim_>
ff_in_mom_idx_space_t &Hubbard<dim_>::GetFormFactorInMomentumIdxSpace(const unsigned n)
{
    auto &[shell_idx, ff_idx] = FormFactors().flat_index_to_tuple_index_map[n];
    return MomentumGrid().form_factors_in_momentum_idx_space[shell_idx][ff_idx];
}

template<unsigned dim_>
ff_in_mom_idx_space_t &Hubbard<dim_>::GetFormFactorInFineMomentumIdxSpace(const unsigned n)
{
    auto &[shell_idx, ff_idx] = FormFactors().flat_index_to_tuple_index_map[n];
    return FineMomentumGrid().form_factors_in_momentum_idx_space[shell_idx][ff_idx];
}

template<unsigned dim_>
unsigned Hubbard<dim_>::SumFineMomentaIdxes(const unsigned p1_idx, const unsigned p2_idx)
{
#ifndef PRECOMPUTE_MOMENTUM_ARITHMETIC
    const unsigned fine_dim = K_DIM * P_IN_K;

  unsigned p1_plus_p2_idx = 0;
  
  unsigned digit_value = 1;

  for (unsigned i = 0; i < dim_; i ++){
    const unsigned p1_coord_i = (p1_idx/digit_value);// % fine_dim;
    const unsigned p2_coord_i = (p2_idx/digit_value);// % fine_dim;

    const unsigned sum = (p1_coord_i + p2_coord_i) % fine_dim;

    p1_plus_p2_idx += sum * digit_value;
    
    digit_value *= fine_dim;
  }

  return p1_plus_p2_idx;
#else
  return FineMomentumGrid().momentum_plus_momentum_idx_map[p1_idx][p2_idx];
#endif
}

template<unsigned dim_>
unsigned Hubbard<dim_>::GetNegativeFineMomentumIdx(const unsigned p_idx)
{
#ifndef PRECOMPUTE_MOMENTUM_ARITHMETIC
    const unsigned fine_dim = K_DIM * P_IN_K;

  unsigned p_negative_idx = 0;
  
  unsigned digit_value = 1;

  for (unsigned i = 0; i < dim_; i ++){
    const unsigned p_coord_i = (p_idx/digit_value) % fine_dim;

    const unsigned negative = (fine_dim - p_coord_i) % fine_dim;

    p_negative_idx += negative * digit_value;
    
    digit_value *= fine_dim;
  }

  return p_negative_idx;
#else
  return FineMomentumGrid().negative_momentum_idx_map[p_idx];
#endif
}

template<unsigned dim_>
unsigned Hubbard<dim_>::GetNegativeCoarseMomentumIdx(const unsigned k_idx)
{
#ifndef PRECOMPUTE_MOMENTUM_ARITHMETIC
    /*const unsigned coarse_dim = K_DIM;

  unsigned k_negative_idx = 0;
  
  unsigned digit_value = 1;

  for (unsigned i = 0; i < dim_; i ++){
    const unsigned k_coord_i = (k_idx/digit_value) % coarse_dim;

    const unsigned negative = (coarse_dim - k_coord_i) % coarse_dim;

    k_negative_idx += negative * digit_value;
    
    digit_value *= coarse_dim;
  }

  return k_negative_idx;*/
  return GetCoarseMomentumIdxFromFine(GetNegativeFineMomentumIdx(GetFineMomentumIdxFromCoarse(k_idx)));
#else
  return MomentumGrid().negative_momentum_idx_map[k_idx];
#endif
}

template<unsigned dim_>
unsigned Hubbard<dim_>::GetFineMomentumIdxFromCoarse(const unsigned k_idx)
{
  /*#ifndef PRECOMPUTE_MOMENTUM_ARITHMETIC
// this part does not work with partial refinement

    const unsigned fine_dim = K_DIM * P_IN_K;

  unsigned p_idx = 0;
  
  unsigned p_digit_value = 1;
  unsigned k_digit_value = 1;
  
  for (unsigned i = 0; i < dim_; i ++){
    const unsigned k_coord_i = (k_idx/k_digit_value) % K_DIM;

    const unsigned p_coord_i = k_coord_i * P_IN_K;

    p_idx += p_coord_i * p_digit_value;

    k_digit_value *= K_DIM;
    p_digit_value *= fine_dim;
  }

  return p_idx;
  #else*/
  return CoarseToFineIdxMap()[k_idx];
  //#endif
}

template<unsigned dim_>
unsigned Hubbard<dim_>::GetCoarseMomentumIdxFromFine(const unsigned p_idx)
{
  /*#ifndef PRECOMPUTE_MOMENTUM_ARITHMETIC
    const unsigned fine_dim = K_DIM * P_IN_K;

  unsigned k_idx = 0;
  
  unsigned p_digit_value = 1;
  unsigned k_digit_value = 1;

  
  for (unsigned i = 0; i < dim_; i ++){
    const unsigned p_coord_i = (p_idx/p_digit_value) % fine_dim;

    const unsigned k_coord_i = p_coord_i / P_IN_K;

    k_idx += k_coord_i * k_digit_value;

    k_digit_value *= K_DIM;
    p_digit_value *= fine_dim;
  }
  return k_idx;
  #else*/
  return FineToCoarseIdxMap()[p_idx];
  //#endif
}



// instantiate the template for the following dimensions
template class Hubbard<0>;
template class Hubbard<1>;
template class Hubbard<2>;
template class Hubbard<3>;
