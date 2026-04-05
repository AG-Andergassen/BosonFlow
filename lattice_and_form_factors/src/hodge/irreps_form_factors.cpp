

irreps_form_factors_t::irreps_form_factors_t()
{
    auto m_point_group;
    auto m_form_factors;
    auto m_reduced_form_factors;
    lattice_2D_t m_lattice;
    unsigned m_shells_count;
    std::vector<std::tuple<unsigned, unsigned, unsigned> > m_vanishing_form_factors;
    std::vector<std::tuple<unsigned, unsigned, unsigned> > m_repeated_form_factors;
}


delta_function_vector_t irreps_form_factors_t::get(unsigned irrep_idx, unsigned shell_idx, unsigned bond_idx)
{
    return m_form_factors[irrep_idx][shell_idx][bond_idx];
}

delta_function_vector_t irreps_form_factors_t::project_delta_function(delta_function_vector_t delta_func, unsigned mu)
{
    unsigned d_mu = m_point_group.character_table[mu][0];

    delta_function_vector_t result;

    for (unsigne class_idx = 0; class_idx < m_point_group.group_conjugacy_classes.size(); class_idx ++){
	for (auto g : m_point_group.group_conjugacy_classes[class_idx]){
	    double prefactor = d_mu/m_point_group.group_size * m_point_group.character_table[mu][class_idx];
	    delta_function_vector_t delta_func_g = delta_func.fundamental_group_act(g);
	    delta_func_g.multiply_scalar(prefactor);
	    result = result + delta_func_g;
	}
    }
    return result;
}

std::vector<std::vector<delta_function_vector_t> > irreps_form_factors_t::make_form_factors_from_shells(std::map<double, std::vector<coord_t> > shells, unsigned mu)
{
    std::vector<std::vector<delta_function_vector_t> > form_factors_shells;

    for (auto itr = shells.begin(); itr != shells.end(); itr ++){
	std::vector<delta_function_vector_t> form_factors;
	for (auto vec : itr->second){
	    delta_function_vector_t bond({vector});
	    delta_function_vector_t new_form_factor = project_delta_function(bond, mu);
	    form_factors.push_back(new_form_factor);
	}
	form_factors_shells.push_back(form_factors);
    }
    return form_factors_shells;
}

double irreps_form_factors_t::compare_idx_space_form_factors(std::vector<std::complex<double> > ff1, std::vector<std::complex<double> > ff2)
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

std::vector<std::tuple<unsigned, unsigned, unsigned> > irreps_form_factors_t::remove_vanishing_and_repeated_form_factors()
{
    std::vector<std::tuple<unsigned, unsigned, unsigned> > removed_indices;
    
    std::vector<std::complex<double> > zero_form_factor;
    zero_form_factor.resize(form_factors_in_momentum_space[0][0][0].size());
    std::fill(zero_form_factor.begin(), zero_form_factor.end(), 0);
    // TODO: add variable shells_radii
    for (unsigned irrep_idx = 0; irrep_idx < m_point_group.group_conjugacy_classes.size(); irrep_idx ++)
    for (unsigned shell_idx = 0; shell_idx < shells_count; shell_idx ++){
	for (unsigned bond_idx = 0; bond_idx < lattice.shells[lattice.shells_radii[shell_idx]]; bond_idx ++){
	    if (std::find(removed_indices.begin(), removed_indices.end(), std::make_tuple(irrep_idx, shell_idx, bond_idx)) != removed_indices.end())
		continue;
	    double deviation = compare_idx_space_form_factors(zero_form_factor, form_factors_in_momentum_idx_space[irrep_idx][shell_idx][bond_idx]);
	    if (deviation < epsilon){
		m_vanishing_form_factors.push_back({irrep_idx, shell_idx, bond_idx});
		removed_indices.push_back({irrep_idx, shell_idx, bond_idx});
	    }
	}
	for (unsigned bond_idx1 = 0; bond_idx1 < lattice.shells[lattice.shells_radii[shell_idx]]; bond_idx1 ++){
	    if (std::find(removed_indices.begin(), removed_indices.end(), std::make_tuple(irrep_idx, shell_idx, bond_idx1)) != removed_indices.end())
		continue;
	    for (unsigned bond_idx2 = 0; bond_idx2 < lattice.shells[lattice.shells_radii[shell_idx]]; bond_idx2 ++){
		if (std::find(removed_indices.begin(), removed_indices.end(), std::make_tuple(irrep_idx, shell_idx, bond_idx2)) != removed_indices.end())
		    continue;
		double deviation = compare_idx_space_form_factors(form_factors_in_momentum_idx_space[irrep_idx][shell_idx][bond_idx1],
								  form_factors_in_momentum_idx_space[irrep_idx][shell_idx][bond_idx2]);
		if (deviation < epsilon){
		    m_repeated_form_factors.push_back({irrep_idx, shell_idx, bond_idx2});
		    removed_indices.push_back({irrep_idx, shell_idx, bond_idx2});
		}
	    }
	}
    }

    std::vector<std::vector<std::vector< std::vector<std::complex<double> > > > > FFs;
    for (unsigned irrep_idx = 0; irrep_idx < m_point_group.group_conjugacy_classes.size(); irrep_idx ++){
        std::vector<std::vector< std::vector<std::complex<double> > > >  irrep_FFs;
	for (unsigned shell_idx = 0; shell_idx < shells_count; shell_idx ++){
	    std::vector< std::vector<std::complex<double> > >  shells_FFs;
	    for (unsigned bond_idx = 0; bond_idx < lattice.shells[lattice.shells_radii[shell_idx]]; bond_idx ++){
		if (std::find(removed_indices.begin(), removed_indices.end(), std::make_tuple(irrep_idx, shell_idx, bond_idx2)) != removed_indices.end())
		    continue;
		shells_FFs.push_back(form_factors_in_momentum_idx_space[irrep_idx][shell_idx][bond_idx]);
	    }
	    irrep_FFs.push_back(shells_FFs);
	}
	FFs.push_back(irrep_FFs);
    }
    
    form_factors_in_momentum_idx_space = FFs;

    return removed_indices;
}


std::vector<std::map<std::tuple<unsigned, unsigned, unsigned>, std::tuple<unsigned, std::string> > > irreps_form_factors_t::get_point_group_in_form_factor_idx_space()
{
    std::vector<std::map<std::tuple<unsigned, unsigned, unsigned>, std::tuple<unsigned, std::string> > > point_group_in_form_factor_idx_space;

    for (auto &g_in_mom_idx_space : point_group.group_in_momentum_idx_space){
	std::map<std::tuple<unsigned, unsigned, unsigned>, std::tuple<unsigned, std::string> > g_in_form_factor_space;
	for (unsigned irrep_idx = 0; irrep_idx < form_factors_in_momentum_idx_space.size(); irrep_idx++)
	for (unsigned shell_idx = 0; shell_idx < form_factors_in_momentum_idx_space[irrep_idx].size(); shell_idx ++)
	for (unsigned bond_idx = 0; bond_idx < form_factors_in_momentum_idx_space[irrep_idx][shell_idx].size(); bond_idx ++){
	    std::vector<std::complex<double> > &ff = form_factors_in_momentum_idx_space[irrep_idx][shell_idx][bond_idx];
	    std::vector<std::complex<double> > ff_prime;
	    for (unsigned mom_idx = 0; mom_idx < ff.size(); mom_idx ++)
		ff_prime.push_back(ff[g_in_mom_idx_space[idx]]);

	    bool found_image = false;

	    for (unsigned bond_idx2 = 0; bond_idx2 < form_factors_in_momentum_idx_space[irrep_idx][shell_idx].size(); bond_idx2 ++){
		std::vector<std::complex<double> > &ff2 = form_factors_in_momentum_idx_space[irrep_idx][shell_idx][bond_idx2];
		double deviation = compare_idx_space_from_form_factors(ff2, ff_prime);
		if (deviation < epsilon){
		    // note that it's enough to specify the "bond index" of the image, since the point group preserves (bond) lengths (i.e; keeps you at same shell and irrep)
		    g_in_form_factor_space[std::make_tuple(irrep_idx, shell_idx, bond_idx)] = std::make_tuple(bond_idx2, "identity");
		    found_image = true;
		    break;
		}

		std::vector<std::complex<double> > minus_ff2;
		for (unsigned mom_idx = 0; mom_idx < ff2.size(); mom_idx ++)
		    minus_ff2.push_back(-ff2[idx]);
		deviation = compare_idx_space_from_form_factors(minus_ff2, ff_prime);
		if (deviation < epsilon){
		    g_in_form_factor_space[std::make_tuple(irrep_idx, shell_idx, bond_idx)] = std::make_tuple(bond_idx2, "minus");
		    found_image = true;
		    break;	   
		}

		std::vector<std::complex<double> > cc_ff2;
		for (unsigned mom_idx = 0; mom_idx < ff2.size(); mom_idx ++)
		    cc_ff2.push_back(std::conj(ff2[idx]));
		deviation = compare_idx_space_from_form_factors(cc_ff2, ff_prime);
		if (deviation < epsilon){
		    g_in_form_factor_space[std::make_tuple(irrep_idx, shell_idx, bond_idx)] = std::make_tuple(bond_idx2, "cc");
		    found_image = true;
		    break;	   
		}

		std::vector<std::complex<double> > cc_with_minus_ff2;
		for (unsigned mom_idx = 0; mom_idx < ff2.size(); mom_idx ++)
		    cc_with_minus_ff2.push_back(-std::conj(ff2[idx]));
		deviation = compare_idx_space_from_form_factors(cc_with_minus_ff2, ff_prime);
		if (deviation < epsilon){
		    g_in_form_factor_space[std::make_tuple(irrep_idx, shell_idx, bond_idx)] = std::make_tuple(bond_idx2, "cc_with_minus");
		    found_image = true;
		    break;
		}
	    } 
	    if (!found_image)
		throw std::invalid_argument("This symmetry doesn't map the form-factor to any existing form factor. This should not happen since the set of form factors that were projected (bond-form factors) map to one another under all the lattice symmetries!");
	}
	point_group_in_form_factor_idx_space.push_back(g_in_form_factor_space);
    }
    return point_group_in_form_factor_idx_space;
}
