#include <form_factor_container.h>

template<unsigned dim>
form_factor_container_t<dim>::form_factor_container_t(unsigned shells_count) : m_shells_count(shells_count)
{
    for (int j = 0; j < shells_count; j++){
	m_form_factors_shells.push_back({});
    }
}

template <unsigned dim>
ff_in_mom_idx_space_t form_factor_container_t<dim>::make_form_factor_in_momentum_idx_space(const momentum_space_t<dim> &momentum_space, delta_function_vector_t<dim> form_factor)
{
    ff_in_mom_idx_space_t idx_form_factor;
    for (coord_t<dim> mom : momentum_space.m_points)
	idx_form_factor.push_back(form_factor.evaluate_fourier_transform(mom)/std::sqrt(momentum_space.get_volume()));
    return idx_form_factor;
}

template <unsigned dim>
std::vector<std::vector<ff_in_mom_idx_space_t> > form_factor_container_t<dim>::make_form_factors_in_momentum_idx_space(const momentum_space_t<dim> &momentum_space)
{
    std::vector<std::vector<ff_in_mom_idx_space_t> > FFs;

    for (unsigned shell_idx = 0; shell_idx < m_shells_count; shell_idx++){
	std::vector<ff_in_mom_idx_space_t> shells_FFs;
	for (unsigned ff_idx = 0; ff_idx < m_form_factors_shells[shell_idx].size(); ff_idx++){
	    shells_FFs.push_back(make_form_factor_in_momentum_idx_space(momentum_space, m_form_factors_shells[shell_idx][ff_idx]));
	}
	FFs.push_back(shells_FFs);
    }
    return FFs;
}



template <unsigned dim>
double form_factor_container_t<dim>::compare_idx_space_form_factors(ff_in_mom_idx_space_t ff1, ff_in_mom_idx_space_t ff2)
{
    double max_deviation = 0;
    bool is_first = true;
    
    if (ff1.size() != ff2.size())
	throw std::invalid_argument("Can't compare form factors defined on different domains! (of different sizes " + std::to_string(ff1.size()) + " and " + std::to_string(ff2.size()) + ")");

    for (unsigned idx = 0; idx < ff1.size(); idx++){
	double deviation = std::abs(ff1[idx] - ff2[idx]);
	if (is_first || deviation > max_deviation){
	    is_first = false;
	    max_deviation = deviation;
	}
    }

    return max_deviation;
}


template <unsigned dim>
std::vector<shell_ff_idx_t> form_factor_container_t<dim>::get_indices()
{
    std::vector<shell_ff_idx_t> indices;
    for (unsigned shell_idx = 0; shell_idx < m_shells_count; shell_idx++){
	for (unsigned ff_idx = 0; ff_idx < m_form_factors_shells[shell_idx].size(); ff_idx++){
	    indices.push_back(std::make_tuple(shell_idx, ff_idx));
	}
    }
    return indices;
}


template <unsigned dim>
delta_function_vector_t<dim> form_factor_container_t<dim>::get(unsigned shell_idx, unsigned ff_idx)
{
    return m_form_factors_shells[shell_idx][ff_idx];
}

template <unsigned dim>
delta_function_vector_t<dim> form_factor_container_t<dim>::get(shell_ff_idx_t shell_ff_idx)
{
    return m_form_factors_shells[std::get<0>(shell_ff_idx)][std::get<1>(shell_ff_idx)];
}

// todo:  code it so that it can be called from outside. Should also check that:
// 1- point group in mom idx space are initialised (in the point group class).
// 2-  

template <unsigned dim>
std::vector<std::map<shell_ff_idx_t, std::tuple<ff_idx_t, std::string> > > form_factor_container_t<dim>::make_point_group_in_shells_form_factor_idx_space(const std::vector<std::vector<unsigned> > &point_group_in_mom_idx_space, const std::vector<std::vector< ff_in_mom_idx_space_t > > &form_factors_shells_in_momentum_idx_space)
{
    std::vector<std::map<shell_ff_idx_t, std::tuple<ff_idx_t, std::string> > > point_group_in_shells_form_factor_idx_space;

   
    for (unsigned g_idx = 0; g_idx < point_group_in_mom_idx_space.size(); g_idx ++){
	std::map<shell_ff_idx_t, std::tuple<ff_idx_t, std::string> > g_in_shell_ff_idx_space;
	for (unsigned shell_idx = 0; shell_idx < m_shells_count; shell_idx++){
	    auto & shells_form_factors = form_factors_shells_in_momentum_idx_space[shell_idx];
	    auto group_in_shell = with_shell_make_point_group_in_form_factor_idx_space(shells_form_factors, point_group_in_mom_idx_space);
	
	    for (auto shell_iter = group_in_shell[g_idx].begin(); shell_iter != group_in_shell[g_idx].end(); shell_iter ++){
		
		shell_ff_idx_t preimage = std::make_tuple(shell_idx, shell_iter->first);
		// note that it's enough to specify the "ff index" of the image without the shell index, since the point group preserves lengths (i.e; keeps you at same shell)
	        g_in_shell_ff_idx_space[preimage] = shell_iter->second;
	    }
	}
	point_group_in_shells_form_factor_idx_space.push_back(g_in_shell_ff_idx_space);
    } 

    return point_group_in_shells_form_factor_idx_space;
    
}


template <unsigned dim>
std::vector<std::map<unsigned, std::tuple<unsigned, std::string> > > form_factor_container_t<dim>::with_shell_make_point_group_in_form_factor_idx_space(const std::vector<ff_in_mom_idx_space_t>  &form_factors_in_momentum_idx_space, const std::vector<std::vector<unsigned> > &point_group_in_mom_idx_space)
{
    std::vector<std::map<unsigned, std::tuple<unsigned, std::string> > > point_group_in_shell_form_factor_idx_space;

    for (auto &g_in_mom_idx_space : point_group_in_mom_idx_space){
	std::map<unsigned, std::tuple<unsigned, std::string> > g_in_form_factor_space;
	for (unsigned ff_idx = 0; ff_idx < form_factors_in_momentum_idx_space.size(); ff_idx ++){
	    const ff_in_mom_idx_space_t &ff = form_factors_in_momentum_idx_space[ff_idx];
	    ff_in_mom_idx_space_t ff_prime;
	    for (unsigned mom_idx = 0; mom_idx < ff.size(); mom_idx ++)
		ff_prime.push_back(ff[g_in_mom_idx_space[mom_idx]]);

	    bool found_image = false;

	    for (unsigned ff_idx2 = 0; ff_idx2 < form_factors_in_momentum_idx_space.size(); ff_idx2 ++){
	        const ff_in_mom_idx_space_t &ff2 = form_factors_in_momentum_idx_space[ff_idx2];
		double deviation = compare_idx_space_form_factors(ff2, ff_prime);
		if (deviation < epsilon){
		    g_in_form_factor_space[ff_idx] = std::make_tuple(ff_idx2, "identity");
		    found_image = true;
		    break;
		}

	        ff_in_mom_idx_space_t minus_ff2;
		for (unsigned mom_idx = 0; mom_idx < ff2.size(); mom_idx ++)
		    minus_ff2.push_back(-ff2[mom_idx]);
		deviation = compare_idx_space_form_factors(minus_ff2, ff_prime);
		if (deviation < epsilon){
		    g_in_form_factor_space[ff_idx] = std::make_tuple(ff_idx2, "minus");
		    found_image = true;
		    break;	   
		}
		// it turned out what follows is not really useful. We want to only consider symmetries of form factors that correspond to a change of the whole object it's in. A minus sign does that, since the form factors enter multiplicatively. But complex conjugation of the form factor doesn't. So keep the followign commented out (keeping it here in case the snippet gets a use in the future) 
		/*
	        ff_in_mom_idx_space_t cc_ff2;
		for (unsigned mom_idx = 0; mom_idx < ff2.size(); mom_idx ++)
		    cc_ff2.push_back(std::conj(ff2[mom_idx]));
		deviation = compare_idx_space_form_factors(cc_ff2, ff_prime);
		if (deviation < epsilon){
		    g_in_form_factor_space[ff_idx] = std::make_tuple(ff_idx2, "cc");
		    found_image = true;
		    break;	   
		}

	        ff_in_mom_idx_space_t cc_with_minus_ff2;
		for (unsigned mom_idx = 0; mom_idx < ff2.size(); mom_idx ++)
		    cc_with_minus_ff2.push_back(-std::conj(ff2[mom_idx]));
		deviation = compare_idx_space_form_factors(cc_with_minus_ff2, ff_prime);
		if (deviation < epsilon){
		    g_in_form_factor_space[ff_idx] = std::make_tuple(ff_idx2, "cc_with_minus");
		    found_image = true;
		    break;
		}
		*/
	    } 
	    if (!found_image)
		throw std::invalid_argument("This symmetry doesn't map the form-factor to any existing form factor. It means your set of form factors break the symmetry of the lattice. In case teh form factors were generated from the bonds or using the irreps method: This should not happen since the set of form factors that were projected (bond-form factors) map to one another under all the lattice symmetries! If it did, it's a bug.");
	}
	point_group_in_shell_form_factor_idx_space.push_back(g_in_form_factor_space);
    }
    return point_group_in_shell_form_factor_idx_space;
}



//
