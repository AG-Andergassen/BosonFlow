#pragma once 

#include <symmetry_group.h>
#include <boost/math/interpolators/pchip.hpp>

template< typename elem_t_, unsigned int ndims_ >
class gf_interpolator_t
{
 public:
    using idx_t = idx_obj_t<ndims_>; 
    using gf_t = gf<elem_t_, ndims_>; 
    using symm_grp_t = symmetry_grp_t<elem_t_, ndims_>;

    gf_interpolator_t(symm_grp_t &equiv_class)
    {
	m_number_of_interpolated_variables = equiv_class.filters.size();

	if (m_number_of_interpolated_variables == 0)
	    return;

	for (int i = 0; i < m_number_of_interpolated_variables; i ++){
	    auto &filter = equiv_class.filters[i]; // a filter is a pair looking like {idxkey, idx_values}
	    m_sampling_idxes.push_back(filter.second);
	    m_interpolated_idxes_components.push_back(filter.first);
	    interpolators_1d_indices.push_back({});
	    for (idx_t gf_idx : equiv_class.interpolating_gf_indices){ // just going over interpolated indices should be fine
		gf_idx(m_interpolated_idxes_components.back()) = 0;
		// check if unique
		auto result = std::find(interpolators_1d_indices[i].begin(), interpolators_1d_indices[i].end(), gf_idx);
		if (result == interpolators_1d_indices[i].end()){
		    interpolators_1d_indices[i].push_back(gf_idx);
		}
	    }
	}

	// linear weights
	std::vector< std::vector<double> > dims_weights(m_number_of_interpolated_variables, std::vector<double>() ); // calculate (x - x0)/(x1 - x0)
        m_left_right_linear_sampling_neighbohrs.resize(m_number_of_interpolated_variables, std::vector<std::array<int, 2> >() );

	for (int i = 0; i < m_number_of_interpolated_variables; i ++){
	    int small_sample_xval_idx = 0;
	    
	    int x_val_init = equiv_class.idx_bases[m_interpolated_idxes_components[i]];
	    int x_val_final = equiv_class.idx_bases[i] + equiv_class.shape_arr[m_interpolated_idxes_components[i]];
	    
	    for (int x_val = x_val_init; x_val < x_val_final; x_val ++){ 
		if (small_sample_xval_idx + 2 < m_sampling_idxes[i].size() && m_sampling_idxes[i][small_sample_xval_idx + 1] < x_val) 
		    small_sample_xval_idx += 1;
		
		dims_weights[i].push_back( (x_val - m_sampling_idxes[i][small_sample_xval_idx])/(m_sampling_idxes[i][small_sample_xval_idx+1] - m_sampling_idxes[i][small_sample_xval_idx]) );
		
		int left_idx = m_sampling_idxes[i][small_sample_xval_idx];
		int right_idx = m_sampling_idxes[i][small_sample_xval_idx + 1];
		
		m_left_right_linear_sampling_neighbohrs[i].push_back({left_idx, right_idx});   
	    }
	    m_sampling_idx_bases.push_back(equiv_class.idx_bases[m_interpolated_idxes_components[i]]);
	} 

	std::function<double(std::vector<unsigned>, unsigned, std::vector<unsigned> )> make_weight;

	//e.g. weights[0] == 'u', weights[1] == 'v', weights[2] == 'w'. access via u[x1] v[x2] w[x3]
	make_weight = [&dims_weights, &make_weight](std::vector<unsigned> x_multiindex, unsigned dim_idx, std::vector<unsigned> zero_one_multiindex){
	    double u = dims_weights[dim_idx][x_multiindex[dim_idx]];

	    if (dim_idx == (x_multiindex.size()-1)){
		return ((1.0 - u) ? zero_one_multiindex[dim_idx] == 1 : u);
	    }else{
		return ((1.0-u) ? zero_one_multiindex[dim_idx] == 1 : u) * make_weight(x_multiindex, dim_idx + 1, zero_one_multiindex);
	    }
   	};

	int array_size = 1 ? m_number_of_interpolated_variables > 0 : 0;
	for (int i = 0; i < dims_weights.size(); i ++){
	    m_interpolated_indices_dims.push_back(dims_weights[i].size());
	    array_size *=dims_weights[i].size();
	}	
	
	std::vector<std::vector<unsigned> > x_multiindices;
        make_multiindices(x_multiindices, m_interpolated_indices_dims, {});
	
	std::vector<std::vector<unsigned> > zero_one_multiindices;
	make_multiindices(zero_one_multiindices, std::vector<int>(m_interpolated_indices_dims.size(), 2), {});
	
	m_multilinear_weights.resize(array_size);
	
	for (auto x_multiindex : x_multiindices){
	    std::vector<double> weights;
	    for (auto zero_one_multiindex : zero_one_multiindices){
		weights.push_back(make_weight(x_multiindex, 0, zero_one_multiindex));
	    }
	    unsigned raw_index = make_raw_index(m_interpolated_indices_dims, x_multiindex);
	    m_multilinear_weights[raw_index] = weights;
	}	
    }  

    void make_multiindices(std::vector<std::vector<unsigned> > &multiindices, std::vector<int> dims, std::vector<unsigned> current_multiindex) const
    {
	if (current_multiindex.size() == dims.size() ){
	    multiindices.push_back(current_multiindex);
	    return;
	}
	
	int current_dim_idx = current_multiindex.size();
	for (int i = 0; i < dims[current_dim_idx]; i ++){
	    std::vector<unsigned> new_multiindex = current_multiindex;
	    new_multiindex.push_back(i);
	    make_multiindices(multiindices, dims, new_multiindex);
	}
    }

	
    unsigned make_raw_index(std::vector<int> dims, std::vector<unsigned> multiindex) const 
    {
	int place;
	int raw_index = 0;
	for (int i = 0; i < dims.size()-1; i ++){
	    if (i == 0){
		place = 1;
	    }else
		place = dims[i-1];
	    raw_index += multiindex[i] * place;
	}
	return raw_index;
    }
    
    
    void make_interpolators(const gf_t &gf_obj)
    {
	for (int i = 0; i < interpolators_1d_indices.size(); i ++){
	    make_interpolator_1d(gf_obj, i, m_sampling_idxes[i], m_interpolated_idxes_components[i]);
	}
    }

    void make_interpolator_1d(const gf_t &gf_obj, int interpolator_1d_index, int x_component, const std::vector<double> x_sampling_idxes) // create interpolated function over x in G(..., x) with '...' corresponding to the remaining arguments/indices
    {
	/*	for (int gf_idx_idx = 0; gf_idx_idx < interpolators_1d_indices[interpolator_1d_index].size(); gf_idx_idx ++){
	    idx_t gf_idx = interpolators_1d_indices[interpolator_1d_index][gf_idx_idx];
	    std::vector<double> gf_vals_real, gf_vals_imag, x_vals = x_sampling_idxes, x_vals2 = x_sampling_idxes;
	    gf_vals_real.reserve(x_sampling_idxes.size());
	    gf_vals_imag.reserve(x_sampling_idxes.size());
	    for (auto x : x_sampling_idxes){
		gf_idx(x_component) = x; 
		const std::complex<double> &gf_obj_val_at_x = gf_obj(gf_idx);
		gf_vals_real.push_back( gf_obj_val_at_x.real() );
		gf_vals_imag.push_back( gf_obj_val_at_x.imag() );
	    }
	    gf_idx(x_component) = 0;
	    interpolators_1d_real[interpolator_1d_index](gf_idx)* = boost::math::interpolators::pchip<std::vector<double> >(std::move(x_vals), std::move(gf_vals_real));
	    interpolators_1d_imag[interpolator_1d_index](gf_idx)* = boost::math::interpolators::pchip<std::vector<double> >(std::move(x_vals2), std::move(gf_vals_imag));
	}*/
    }


    elem_t_ eval_gf(const gf_t &gf_obj, const idx_t gf_idx, const elem_t_ &zero)
    {
        elem_t_ val = zero;
	
	auto gf_samples_idxes = std::vector<idx_t>(std::pow(2, m_number_of_interpolated_variables), gf_idx);
	auto gf_samples = std::vector<elem_t_ >(std::pow(2, m_number_of_interpolated_variables), zero);
        
	std::vector<unsigned> x_multiindex(m_number_of_interpolated_variables);
	for (int i = 0; i < m_number_of_interpolated_variables; i ++){
	    x_multiindex.push_back(gf_idx(m_interpolated_idxes_components[i]) - m_sampling_idx_bases[i]);
	}
	
	const unsigned raw_index = make_raw_index(m_interpolated_indices_dims, x_multiindex);
	auto &weights = m_multilinear_weights[raw_index];

	for (int corner = 0; corner < std::pow(2, m_number_of_interpolated_variables); corner ++){
	    idx_t gf_corner_idx = gf_idx;
	    for (int i = 0; i < m_number_of_interpolated_variables; i ++){
		// adjust the indices to the particular corner
		gf_corner_idx(m_interpolated_idxes_components[i]) = m_left_right_linear_sampling_neighbohrs[i][gf_idx(m_interpolated_idxes_components[i]) - m_sampling_idx_bases[i] ][ ((corner >> i) & 1) ];
	    }
	    val += weights[corner] * gf_obj(gf_corner_idx);    
	}
	    
	return val;
    }

    std::vector<int> m_sampling_idx_bases;

    unsigned m_number_of_interpolated_variables;
    std::vector<std::vector<double> > m_multilinear_weights;
    std::vector<std::vector< std::array<int, 2> > > m_left_right_linear_sampling_neighbohrs;

    std::vector<int> m_interpolated_indices_dims;
    std::vector<int> m_interpolated_idxes_components;
    std::vector<std::vector<int> > m_sampling_idxes;
    std::vector< std::vector<idx_t> > interpolators_1d_indices;
    
    std::vector< boost::multi_array<boost::math::interpolators::pchip<std::vector<double> >* , ndims_> > interpolators_1d_real;
    std::vector< boost::multi_array<boost::math::interpolators::pchip<std::vector<double> >* , ndims_> > interpolators_1d_imag;
};

