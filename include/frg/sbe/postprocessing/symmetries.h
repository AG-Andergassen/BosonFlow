
#pragma once

#include <frg/common_frg_gf_types.h>
#include <frg/symmetries_common.h>

template <typename Model>
class SymmetriesfRGSBEPostprocessing : public SymmetriesfRGCommon<Model>
{
 public:
    static void Init();
    
    static void CleanUp();


    // index equivalence classes for postproc susc sc/d/m
    static symmetry_grp_t<dcomplex,8>* &IdxEquivClasses_susc_sc_Ptr(){
	static symmetry_grp_t<dcomplex,8> *idx_equiv_classes_susc_sc_ptr(nullptr);
	return idx_equiv_classes_susc_sc_ptr;
    }
    static symmetry_grp_t<dcomplex,8>* &IdxEquivClasses_susc_d_Ptr(){
	static symmetry_grp_t<dcomplex,8> *idx_equiv_classes_susc_d_ptr(nullptr);
	return idx_equiv_classes_susc_d_ptr;
    }
    static symmetry_grp_t<dcomplex,8>* &IdxEquivClasses_susc_m_Ptr(){
	static symmetry_grp_t<dcomplex,8> *idx_equiv_classes_susc_m_ptr(nullptr);
	return idx_equiv_classes_susc_m_ptr;
    }

    // index equivalence classes for postproc polarisation sc/d/m
    static symmetry_grp_t<dcomplex,2>* &IdxEquivClasses_polarisation_sc_Ptr(){
	static symmetry_grp_t<dcomplex,2> *idx_equiv_classes_polarisation_sc_ptr(nullptr);
	return idx_equiv_classes_polarisation_sc_ptr;
    }
    static symmetry_grp_t<dcomplex,2>* &IdxEquivClasses_polarisation_d_Ptr(){
	static symmetry_grp_t<dcomplex,2> *idx_equiv_classes_polarisation_d_ptr(nullptr);
	return idx_equiv_classes_polarisation_d_ptr;
    }
    static symmetry_grp_t<dcomplex,2>* &IdxEquivClasses_polarisation_m_Ptr(){
	static symmetry_grp_t<dcomplex,2> *idx_equiv_classes_polarisation_m_ptr(nullptr);
	return idx_equiv_classes_polarisation_m_ptr;
    }


    // index equivalence classes for postproc lambda sc/d/m
    static symmetry_grp_t<dcomplex,2>* &IdxEquivClasses_static_lambda_sc_Ptr(){
	static symmetry_grp_t<dcomplex,2> *idx_equiv_classes_static_lambda_sc_ptr(nullptr);
	return idx_equiv_classes_static_lambda_sc_ptr;
    }
    static symmetry_grp_t<dcomplex,2>* &IdxEquivClasses_static_lambda_d_Ptr(){
	static symmetry_grp_t<dcomplex,2> *idx_equiv_classes_static_lambda_d_ptr(nullptr);
	return idx_equiv_classes_static_lambda_d_ptr;
    }
    static symmetry_grp_t<dcomplex,2>* &IdxEquivClasses_static_lambda_m_Ptr(){
	static symmetry_grp_t<dcomplex,2> *idx_equiv_classes_static_lambda_m_ptr(nullptr);
	return idx_equiv_classes_static_lambda_m_ptr;
    }


 protected:
    static void IncludeAlgebraicSymmetries();
    
    static void IncludeLatticeSymmetries();

    static void CalculateIndexEquivalenceClasses();

    static bool &AreLatticeSymmetriesIncluded()
    {
	static bool are_lattice_symmetries_included = false; // initialised value
	return are_lattice_symmetries_included;
    }

    static bool &AreAlgebraicSymmetriesIncluded()
    {
	static bool are_algebraic_symmetries_included = false;
	return are_algebraic_symmetries_included;
    }

    // Symmetry group actions for susc sc/d/m
    static std::vector<boost::function<operation (idx_susc_t<Model>&)> > &GroupActions_susc_sc(){
	static std::vector<boost::function<operation (idx_susc_t<Model>&)> > group_actions_susc_sc({});
	return group_actions_susc_sc;
    }
    static std::vector<boost::function<operation (idx_susc_t<Model>&)> > &GroupActions_susc_d(){
	static std::vector<boost::function<operation (idx_susc_t<Model>&)> > group_actions_susc_d({});
	return group_actions_susc_d;
    }
    static std::vector<boost::function<operation (idx_susc_t<Model>&)> > &GroupActions_susc_m(){
	static std::vector<boost::function<operation (idx_susc_t<Model>&)> > group_actions_susc_m({});
	return group_actions_susc_m;
    }

    // Symmetry group actions for polarisation sc/d/m
    static std::vector<boost::function<operation (idx_polarisation_t<Model>&)> > &GroupActions_polarisation_sc(){
	static std::vector<boost::function<operation (idx_polarisation_t<Model>&)> > group_actions_polarisation_sc({});
	return group_actions_polarisation_sc;
    }
    static std::vector<boost::function<operation (idx_polarisation_t<Model>&)> > &GroupActions_polarisation_d(){
	static std::vector<boost::function<operation (idx_polarisation_t<Model>&)> > group_actions_polarisation_d({});
	return group_actions_polarisation_d;
    }
    static std::vector<boost::function<operation (idx_polarisation_t<Model>&)> > &GroupActions_polarisation_m(){
	static std::vector<boost::function<operation (idx_polarisation_t<Model>&)> > group_actions_polarisation_m({});
	return group_actions_polarisation_m;
    }

    
    // Symmetry group actions for postproc lambda sc/d/m
    static std::vector<boost::function<operation (idx_static_lambda_t<Model>&)> > &GroupActions_static_lambda_sc(){
	static std::vector<boost::function<operation (idx_static_lambda_t<Model>&)> > group_actions_static_lambda_sc({});
	return group_actions_static_lambda_sc;
    }
    static std::vector<boost::function<operation (idx_static_lambda_t<Model>&)> > &GroupActions_static_lambda_d(){
	static std::vector<boost::function<operation (idx_static_lambda_t<Model>&)> > group_actions_static_lambda_d({});
	return group_actions_static_lambda_d;
    }
    static std::vector<boost::function<operation (idx_static_lambda_t<Model>&)> > &GroupActions_static_lambda_m(){
	static std::vector<boost::function<operation (idx_static_lambda_t<Model>&)> > group_actions_static_lambda_m({});
	return group_actions_static_lambda_m;
    }

    // Lattice Symmetries
    static std::function<operation(idx_susc_t<Model>& )> MakeLatticeSymmetryMap_susc(unsigned point_group_element_idx);

    
    static std::function<operation(idx_polarisation_t<Model>& )> MakeLatticeSymmetryMap_polarisation(unsigned point_group_element_idx);

    static std::function<operation(idx_static_lambda_t<Model>& )> MakeLatticeSymmetryMap_static_lambda(unsigned point_group_element_idx);
};


/**
 * \brief for postprocessing susceptibility
 */
template<typename Model> operation cconj_timerev_susc_sc( idx_susc_t<Model>& idx );	
template<typename Model> operation cconj_timerev_susc_d( idx_susc_t<Model>& idx );	
template<typename Model> operation cconj_timerev_susc_m( idx_susc_t<Model>& idx );

template<typename Model> operation timerev_susc_sc( idx_susc_t<Model>& idx);	
template<typename Model> operation hmirror_timerev_susc_d( idx_susc_t<Model>& idx);	
template<typename Model> operation timerev_susc_m( idx_susc_t<Model>& idx);	
