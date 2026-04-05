#include <ir_basis.h>


IR_2p_basis_data_t IR_basis_factory::make_2p_basis(int total_number_of_freqs, double bandwidth, Statistics zeta, bool augment_with_constant)
{
    initialise_julia_context_and_headers();

    jl_eval_string("beta = " + std::to_string(m_beta));
    jl_eval_string("omega_max = " + std::to_string(bandwidth));
    jl_eval_string("epsilon = " + std::to_string(m_error));
    
    jl_eval_string("basis_f, basis_b = SparseIR.finite_temp_bases(beta, omega_max, epsilon)");
    
    IR_2p_basis_data_t basis_data;

    std::string julia_basis_string;
    
    if (zeta == Statistics::FERMI){
	julia_basis_string = "basis_f";
    }else{
	julia_basis_string = "basis_b";
    }

    if (augment_with_constant)
	julia_basis_string = "SparseIR.MatsubaraSampling(AugmentedBasis(" + julia_basis_string+", MatsubaraConst))";

    jl_eval_string("matsubara_sampling = SparseIR.MatsubaraSampling("+julia_basis_string+")");
    jl_eval_string("everywhere_sampling = SparseIR.MatsubaraSampling("+julia_basis_string+", sampling_points="+make_julia_2p_freq_array(total_number_of_freqs)+")");
    get_julia_2Darray_complex("everywhere_sampling.matrix * pinv(matsubara_sampling.matrix)", basis_data.E);


    if (zeta == Statistics::FERMI)
	get_julia_array( "(Int.(matsubara_sampling.sampling_points) .- 1) ./ 2 ", basis_data.sampling_freqs);
    else
	get_julia_array( "Int.(matsubara_sampling.sampling_points) ./ 2 ", basis_data.sampling_freqs);
    
   
    destroy_julia_context();
    return basis_data;
}


IR_3p_basis_data_t IR_basis_factory::make_3p_basis(int total_number_of_freqs_bos, double bandwidth_bos, int total_number_of_freqs_ferm, double bandwidth_ferm, bool augment_with_constant)
{
    initialise_julia_context_and_headers();

    jl_eval_string("beta = " + std::to_string(m_beta));
    jl_eval_string("omega_max_b = " + std::to_string(bandwidth_bos));
    jl_eval_string("omega_max_f = " + std::to_string(bandwidth_ferm));
    jl_eval_string("epsilon = " + std::to_string(m_error));
    
    jl_eval_string("(basis_f, _) = SparseIR.finite_temp_bases(beta, omega_max_f, epsilon)");
    jl_eval_string("(_, basis_b) = SparseIR.finite_temp_bases(beta, omega_max_b, epsilon)");

    if (augment_with_constant)
	jl_eval_string("basis3 = OvercompleteBasis(OvercompleteIR.DEFAULT_THREEPOINT_SET,\
			       AugmentedBasis(basis_f, MatsubaraConst),	\
			       AugmentedBasis(basis_b, MatsubaraConst))");
    else
	jl_eval_string("basis3 = OvercompleteBasis(OvercompleteIR.DEFAULT_THREEPOINT_SET, basis_f, basis_b)");
    
    jl_eval_string("eval3_sampling = OvercompleteIR.MatsubaraEval(basis3)");
    get_julia_tuple4array<int>( "map(w -> Int.(Int(w[1])/2, (Int(w[2])-1)/2, (Int(w[3])-1)/2 ), eval3_sampling.wsample)", basis_data.sampling_freqs);
    
    {
	std::string pp_interpolating_freqs = make_julia_4p_freq_array([](int W, int w){ 
		int w2, w4;
		Symmetrised3P::pp_to_full_freqs(W, w, w2, w4);
		return std::make_tuple(2*W, 2*w2+1, 2*w4+1);
	    }, total_number_of_freqs_bos, total_number_of_freqs_ferm);
	std::string ph_interpolating_freqs = make_julia_4p_freq_array([](int W, int w){ 
		int w2, w3;
		Symmetrised3P::ph_to_full_freqs(W, w, w2, w3);
		return std::make_tuple(2*W, 2*w2+1, 2*w3+1);
	    }, total_number_of_freqs_bos, total_number_of_freqs_ferm);
	std::string xph_interpolating_freqs = make_julia_4p_freq_array([](int W, int w){ 
		int w1, w2;
		Symmetrised3P::xph_to_full_freqs(W, w, w1, w2);
		return std::make_tuple(2*W, 2*w1+1, 2*w2+1);
	    }, total_number_of_freqs_bos, total_number_of_freqs_ferm);
	
	jl_eval_string("eval3_pp = OvercompleteIR.MatsubaraEval(basis3, sampling_points="+pp_interpolating_freqs+")");
	jl_eval_string("eval3_ph = OvercompleteIR.MatsubaraEval(basis3, sampling_points="+ph_interpolating_freqs+")");
	jl_eval_string("eval3_xph = OvercompleteIR.MatsubaraEval(basis3, sampling_points="+xph_interpolating_freqs+")");
	
	get_julia_2Darray_complex("eval3_pp.trans * pinv(eval_3_sampling.trans)", basis_data.E_pp);
	get_julia_2Darray_complex("eval3_ph.trans * pinv(eval_3_sampling.trans)", basis_data.E_ph);
	get_julia_2Darray_complex("eval3_xph.trans * pinv(eval_3_sampling.trans)", basis_data.E_xph);
    }

    destroy_julia_context();
    return basis_data;
}


IR_4p_basis_data_t IR_basis_factory::make_4p_basis(int total_number_of_freqs_bos, double bandwidth_bos, int total_number_of_freqs_ferm, double bandwidth_ferm, bool augment_with_constant)
{
    initialise_julia_context_and_headers();

    jl_eval_string("beta = " + std::to_string(m_beta));
    jl_eval_string("omega_max_b = " + std::to_string(bandwidth_bos));
    jl_eval_string("omega_max_f = " + std::to_string(bandwidth_ferm));
    jl_eval_string("epsilon = " + std::to_string(m_error));
    
    jl_eval_string("(basis_f, _) = SparseIR.finite_temp_bases(beta, omega_max_f, epsilon)");
    jl_eval_string("(_, basis_b) = SparseIR.finite_temp_bases(beta, omega_max_b, epsilon)");

    if (augment_with_constant)
	jl_eval_string("basis4 = OvercompleteBasis(OvercompleteIR.DEFAULT_FOURPOINT_SET,\
			       AugmentedBasis(basis_f, MatsubaraConst),	\
			       AugmentedBasis(basis_b, MatsubaraConst))");
    else
	jl_eval_string("basis4 = OvercompleteBasis(OvercompleteIR.DEFAULT_FOURPOINT_SET, basis_f, basis_b)");
    
    jl_eval_string("eval4_sampling = OvercompleteIR.MatsubaraEval(basis4)");
    get_julia_tuple4array<int>( "map(w -> (Int.(w).- 1)./2, eval4_sampling.wsample)", basis_data.sampling_freqs);
    
    {
	std::string pp_interpolating_freqs = make_julia_4p_freq_array([](int W, int w, int wp){ 
		int w1, w2, w3, w4;
		Symmetrised4P::pp_to_full_freqs(W, w, wp, w1, w2, w3, w4);
		return std::make_tuple(2*w1+1, 2*w2+1, 2*w3+1, 2*w4+1);
	    }, total_number_of_freqs_bos, total_number_of_freqs_ferm, total_number_of_freqs_ferm);
	std::string ph_interpolating_freqs = make_julia_4p_freq_array([](int W, int w, int wp){ 
		int w1, w2, w3, w4;
		Symmetrised4P::ph_to_full_freqs(W, w, wp, w1, w2, w3, w4);
		return std::make_tuple(2*w1+1, 2*w2+1, 2*w3+1, 2*w4+1);
	    }, total_number_of_freqs_bos, total_number_of_freqs_ferm, total_number_of_freqs_ferm);
	std::string xph_interpolating_freqs = make_julia_4p_freq_array([](int W, int w, int wp){ 
		int w1, w2, w3, w4;
		Symmetrised4P::xph_to_full_freqs(W, w, wp, w1, w2, w3, w4);
		return std::make_tuple(2*w1+1, 2*w2+1, 2*w3+1, 2*w4+1);
	    }, total_number_of_freqs_bos, total_number_of_freqs_ferm, total_number_of_freqs_ferm);
	
	jl_eval_string("eval4_pp = OvercompleteIR.MatsubaraEval(basis4, sampling_points="+pp_interpolating_freqs+")");
	jl_eval_string("eval4_ph = OvercompleteIR.MatsubaraEval(basis4, sampling_points="+ph_interpolating_freqs+")");
	jl_eval_string("eval4_xph = OvercompleteIR.MatsubaraEval(basis4, sampling_points="+xph_interpolating_freqs+")");
	
	get_julia_2Darray_complex("eval4_pp.trans * pinv(eval4_sampling.trans)", basis_data.E_pp);
	get_julia_2Darray_complex("eval4_ph.trans * pinv(eval4_sampling.trans)", basis_data.E_ph);
	get_julia_2Darray_complex("eval4_xph.trans * pinv(eval4_sampling.trans)", basis_data.E_xph);
    }

    destroy_julia_context();
    return basis_data;
}


IR_2p_basis_data_t IR_basis_factory::make_2p_basis_with_max_sampling_frequency(int max_sampling_frequency, Statistics zeta, bool augment_with_constant)
{
    double bandwidth = obtain_bandwidth_from_max_sampling_frequency(max_sampling_frequency, zeta);
    return make_2p_basis(max_sampling_frequency, bandwidth, zeta, augment_with_constant);
}

IR_3p_basis_data_t IR_basis_factory::make_3p_basis_with_max_sampling_frequency(int max_sampling_frequency_bos, int max_sampling_frequency_ferm, bool augment_with_constant)
{
    double bandwidth_ferm = obtain_bandwidth_from_max_sampling_frequency(max_sampling_frequency, Statistics::FERMI);
    double bandwidth_bos = obtain_bandwidth_from_max_sampling_frequency(max_sampling_frequency, Statistics::BOSE);
    return make_3p_basis(max_sampling_frequency_bos, bandwidth_bos, max_sampling_frequency_ferm, bandwidth_ferm, augment_with_constant);
}

IR_4p_basis_data_t IR_basis_factory::make_3p_basis_with_max_sampling_frequency(int max_sampling_frequency_bos, int max_sampling_frequency_ferm, bool augment_with_constant)
{
    double bandwidth_ferm = obtain_bandwidth_from_max_sampling_frequency(max_sampling_frequency, Statistics::FERMI);
    double bandwidth_bos = obtain_bandwidth_from_max_sampling_frequency(max_sampling_frequency, Statistics::BOSE);
    return make_4p_basis(max_sampling_frequency_bos, bandwidth_bos, max_sampling_frequency_ferm, bandwidth_ferm, augment_with_constant);
}


double IR_basis_factory::obtain_bandwidth_from_max_sampling_frequency(unsigned max_sampling_frequency, Statistics zeta)
{
    initialise_julia_context_and_headers();
    
    jl_eval_string("beta = " + std::to_string(m_beta));
    jl_eval_string("epsilon = " + std::to_string(m_error));
    
    // note: function is monotonically increasing (step wise)
    auto delta_max_sampling_frequency = [max_sampling_frequency, zeta](double bandwidth){
	jl_eval_string("omega_max = " + std::to_string(bandwidth));
	jl_eval_string("basis_f, basis_b = SparseIR.finite_temp_bases(beta, omega_max, epsilon)");
	std::vector<int> sampling_points;  
	if (zeta == Statistics::FERMI)
	    get_julia_array( "(Int.(SparseIR.MatsubaraSampling(basis_f).sampling_points) .- 1) ./ 2 ", sampling_points);  
	else
	    get_julia_array( "Int.(SparseIR.MatsubaraSampling(basis_b).sampling_points) ./ 2 ", sampling_points);  
	return std::max(std::abs(sampling_points[0]), std::abs(sampling_points.back())) -  max_sampling_frequency;
    };
    
    // use bisection method
    double a = 0.1, b = 10000.0;
    int f_a = delta_max_sampling_frequency(a);
    int f_b = delta_max_sampling_frequency(b);
    assert((f_a * f_b < 0), "f_a and f_b must have opposite signs.");

    int f_c, c;
    do{
	c = b-a/2.0;
	f_c = delta_max_sampling_frequency(c);
	if (f_c > 0)
	    b = c;
	else
	    a = c;
		
    } while (f_c != 0);
    
    destroy_julia_context();

    int desired_bandwidth = c;
    return desired_bandwidth;
}

void IR_basis_factory::initialise_julia_context_and_headers()
{
    jl_init();

    jl_eval_string(R"(
using OvercompleteIR
using OvercompleteIR.SparseIR
using OvercompleteIR.LinearAlgebra
)");

}

void IR_basis_factory::destroy_julia_context()
{
    jl_atexit_hook(0);
}

std::string IR_basis_factory::make_julia_2p_freq_array(std::vector<int> freq_indices)
{
    std::string str = "SparseIR.Fermionic_Freq.([";

    str += std::to_string(freq_indices[0]);
    for (int i = 1; i < freq_indices.size(); i ++){
	str += "," + std::to_string(freq_indices[i] * 2 + 1);
    }
    str += "]);";

    return str;
}

std::string IR_basis_factory::make_julia_2p_freq_array(unsigned total_number_of_freqs)
{
    std::string str;

    if (total_number_of_freqs % 2 == 1)
	str = "SparseIR.Bosonic_Freq.([";
    else
	str = "SparseIR.Fermionic_Freq.([";
	
    int first_freq = -(total_number_of_freqs - 1);
    int final_freq = total_number_of_freqs - 1;
    

    str += std::to_string(first_freq);
    for (int i = first_freq + 2; i <= final_freq; i += 2){
	str += "," + std::to_string(i);
    } 

    str += "]);";

    return str;
}

std::string IR_basis_factory::make_julia_2p_freq_array(unsigned total_number_of_freqs)
{
    std::string str;
    if (total_number_of_freqs % 2 == 1)
	str = "SparseIR.Bosonic_Freq.([";
    else
	str = "SparseIR.Fermionic_Freq.([";
	
    int first_freq = -(total_number_of_freqs - 1);
    int final_freq = total_number_of_freqs - 1;
    
    str += std::to_string(first_freq);
    for (int i = first_freq + 2; i <= final_freq; i += 2)
	str += "," + std::to_string(i);
 
    str += "]);";

    return str;
}


std::string IR_basis_factory::make_julia_3p_freq_array(std::function<std::tuple<int, int, int>(int, int) > conv_freq_to_full_freq, unsigned total_number_of_freqs1, unsigned total_number_of_freqs2)
{
    std::string str = "[";

    auto make_3_freq = [](std::tuple<int, int, int> freqs){
	return "(Bosonic_Freq(" + std::to_string(std::get<0>(freqs)) + "),Fermionic_Freq(" + std::to_string(std::get<1>(freqs)) + "), Fermionic_Freq(" + std::to_string(std::get<2>(freqs)) + "))";
    };
    
    int first_freq1 = -(total_number_of_freqs1 - 1);
    int final_freq1 = total_number_of_freqs1 - 1;
    int first_freq2 = -(total_number_of_freqs2 - 1);
    int final_freq2 = total_number_of_freqs2 - 1;
    
    str += make_3_freq(std::tuple(first_freq1, first_freq2) );
    for (int omega = first_freq1 + 2; omega <= final_freq1; omega += 2)
	for (int nu = first_freq2 + 2; nu <= final_freq2; nu += 2)
	    str += "," + make_3_freq(std::tuple(omega, nu));
 
    str += "];";

    return str;
}

std::string IR_basis_factory::make_julia_4p_freq_array(std::function<std::tuple<int, int, int, int>(int, int, int) > conv_freq_to_full_freq, unsigned total_number_of_freqs1, unsigned total_number_of_freqs2, unsigned total_number_of_freqs3)
{
    std::string str = "[";

    auto make_4_freq = [](std::tuple<int, int, int, int> freqs){
	return "(Fermionic_Freq(" + std::to_string(std::get<0>(freqs)) + "),Fermionic_Freq(" + std::to_string(std::get<1>(freqs)) + "), Fermionic_Freq(" + std::to_string(std::get<2>(freqs)) + "), Fermionic_Freq(" + std::to_string(std::get<3>(freqs)) + "))";
    };
    
    int first_freq1 = -(total_number_of_freqs1 - 1);
    int final_freq1 = total_number_of_freqs1 - 1;
    int first_freq2 = -(total_number_of_freqs2 - 1);
    int final_freq2 = total_number_of_freqs2 - 1;
    int first_freq3 = -(total_number_of_freqs3 - 1);
    int final_freq3 = total_number_of_freqs3 - 1;
    
    str += make_4_freq(conv_freq_to_full_freq(first_freq1, first_freq2, first_freq3) );
    for (int omega = first_freq1 + 2; omega <= final_freq1; omega += 2)
	for (int nu = first_freq2 + 2; nu <= final_freq2; nu += 2)
	    for (int nup = first_freq3 + 2; nup <= final_freq3; nup += 2)
		str += "," + make_4_freq(conv_freq_to_full_freq(omega, nu, nup));
    
    str += "];";
    
    return str;
}
