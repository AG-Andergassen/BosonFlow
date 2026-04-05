
#pragma once

#include <frg/sbe/state_gf_types.h>
#include <frg/symmetries_common.h>

template <typename Model>
class Symmetries1lfRGSBE : public SymmetriesfRGCommon<Model>
{
 public:
    static void Init();
    
    static void CleanUp();


    // index equivalence classes for w sc/d/m
    static symmetry_grp_t<dcomplex,2>* &IdxEquivClasses_w_sc_Ptr(){
	static symmetry_grp_t<dcomplex,2> *idx_equiv_classes_w_sc_ptr(nullptr);
	return idx_equiv_classes_w_sc_ptr;
    }
    static symmetry_grp_t<dcomplex,2>* &IdxEquivClasses_w_d_Ptr(){
	static symmetry_grp_t<dcomplex,2> *idx_equiv_classes_w_d_ptr(nullptr);
	return idx_equiv_classes_w_d_ptr;
    }
    static symmetry_grp_t<dcomplex,2>* &IdxEquivClasses_w_m_Ptr(){
	static symmetry_grp_t<dcomplex,2> *idx_equiv_classes_w_m_ptr(nullptr);
	return idx_equiv_classes_w_m_ptr;
    }

    // index equivalence classes for lambda sc/d/m
    static symmetry_grp_t<dcomplex,4>* &IdxEquivClasses_lambda_sc_Ptr(){
	static symmetry_grp_t<dcomplex,4> *idx_equiv_classes_lambda_sc_ptr(nullptr);
	return idx_equiv_classes_lambda_sc_ptr;
    }
    static symmetry_grp_t<dcomplex,4>* &IdxEquivClasses_lambda_d_Ptr(){
	static symmetry_grp_t<dcomplex,4> *idx_equiv_classes_lambda_d_ptr(nullptr);
	return idx_equiv_classes_lambda_d_ptr;
    }
    static symmetry_grp_t<dcomplex,4>* &IdxEquivClasses_lambda_m_Ptr(){
	static symmetry_grp_t<dcomplex,4> *idx_equiv_classes_lambda_m_ptr(nullptr);
	return idx_equiv_classes_lambda_m_ptr;
    }

    // index equivalence classes for M sc/d/m
    static symmetry_grp_t<dcomplex,6>* &IdxEquivClasses_M_sc_Ptr(){
	static symmetry_grp_t<dcomplex,6> *idx_equiv_classes_M_sc_ptr(nullptr);
	return idx_equiv_classes_M_sc_ptr;
    }
    static symmetry_grp_t<dcomplex,6>* &IdxEquivClasses_M_d_Ptr(){
	static symmetry_grp_t<dcomplex,6> *idx_equiv_classes_M_d_ptr(nullptr);
	return idx_equiv_classes_M_d_ptr;
    }
    static symmetry_grp_t<dcomplex,6>* &IdxEquivClasses_M_m_Ptr(){
	static symmetry_grp_t<dcomplex,6> *idx_equiv_classes_M_m_ptr(nullptr);
	return idx_equiv_classes_M_m_ptr;
    }
  
    // index equivalence classes for projected nabla sc/d/m
    static symmetry_grp_t<dcomplex,6>* &IdxEquivClasses_projected_nabla_pp_Ptr(){
	static symmetry_grp_t<dcomplex,6> *idx_equiv_classes_projected_nabla_pp_ptr(nullptr);
	return idx_equiv_classes_projected_nabla_pp_ptr;
    }
    static symmetry_grp_t<dcomplex,6>* &IdxEquivClasses_projected_nabla_ph_Ptr(){
	static symmetry_grp_t<dcomplex,6> *idx_equiv_classes_projected_nabla_ph_ptr(nullptr);
	return idx_equiv_classes_projected_nabla_ph_ptr;
    }
    static symmetry_grp_t<dcomplex,6>* &IdxEquivClasses_projected_nabla_xph_Ptr(){
	static symmetry_grp_t<dcomplex,6> *idx_equiv_classes_projected_nabla_xph_ptr(nullptr);
	return idx_equiv_classes_projected_nabla_xph_ptr;
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

    // Symmetry group actions for w sc/d/m
    static std::vector<boost::function<operation (idx_w_t<Model>&)> > &GroupActions_w_sc(){
	static std::vector<boost::function<operation (idx_w_t<Model>&)> > group_actions_w_sc({});
	return group_actions_w_sc;
    }
    static std::vector<boost::function<operation (idx_w_t<Model>&)> > &GroupActions_w_d(){
	static std::vector<boost::function<operation (idx_w_t<Model>&)> > group_actions_w_d({});
	return group_actions_w_d;
    }
    static std::vector<boost::function<operation (idx_w_t<Model>&)> > &GroupActions_w_m(){
	static std::vector<boost::function<operation (idx_w_t<Model>&)> > group_actions_w_m({});
	return group_actions_w_m;
    }

    // Symmetry group actions for lambda sc/d/m
    static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > &GroupActions_lambda_sc(){
	static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > group_actions_lambda_sc({});
	return group_actions_lambda_sc;
    }
    static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > &GroupActions_lambda_d(){
	static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > group_actions_lambda_d({});
	return group_actions_lambda_d;
    }
    static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > &GroupActions_lambda_m(){
	static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > group_actions_lambda_m({});
	return group_actions_lambda_m;
    }

    // Symmetry group actions for M sc/d/m
    static std::vector<boost::function<operation (idx_M_t<Model>&)> > &GroupActions_M_sc(){
	static std::vector<boost::function<operation (idx_M_t<Model>&)> > group_actions_M_sc({});
	return group_actions_M_sc;
    }
    static std::vector<boost::function<operation (idx_M_t<Model>&)> > &GroupActions_M_d(){
	static std::vector<boost::function<operation (idx_M_t<Model>&)> > group_actions_M_d({});
	return group_actions_M_d;
    }
    static std::vector<boost::function<operation (idx_M_t<Model>&)> > &GroupActions_M_m(){
	static std::vector<boost::function<operation (idx_M_t<Model>&)> > group_actions_M_m({});
	return group_actions_M_m;
    }

    // Lattice Symmetries
    static std::function<operation(idx_w_t<Model>& )> MakeLatticeSymmetryMap_w(unsigned point_group_element_idx);

    static std::function<operation(idx_lambda_t<Model>& )> MakeLatticeSymmetryMap_lambda(unsigned point_group_element_idx);

    static std::function<operation(idx_M_t<Model>& )> MakeLatticeSymmetryMap_M(unsigned point_group_element_idx);
};



/**
 * \brief for the hedin vertex lambda
 */
template<typename Model> operation cconj_timerev_lambda_sc( idx_lambda_t<Model>& idx );	
template<typename Model> operation cconj_timerev_lambda_d( idx_lambda_t<Model>& idx );	
template<typename Model> operation cconj_timerev_lambda_m( idx_lambda_t<Model>& idx );


/**
 * \brief for screened interaction w
 */
template<typename Model> operation hmirror_w_sc( idx_w_t<Model>& idx );	
template<typename Model> operation hmirror_w_d( idx_w_t<Model>& idx );	
template<typename Model> operation hmirror_w_m( idx_w_t<Model>& idx );
template<typename Model> operation cconj_w_sc( idx_w_t<Model>& idx );	
template<typename Model> operation cconj_w_ph( idx_w_t<Model>& idx );	
template<typename Model> operation cconj_w_m( idx_w_t<Model>& idx );	
template<typename Model> operation timerev_w_sc( idx_w_t<Model>& idx);	
template<typename Model> operation timerev_w_ph( idx_w_t<Model>& idx);	
template<typename Model> operation timerev_w_m( idx_w_t<Model>& idx);	


/**
 * \brief for U-irreducible but two-particle-reducible vertex M
 */
template<typename Model> operation cconj_timerev_M_sc( idx_M_t<Model>& idx );
template<typename Model> operation cconj_timerev_M_d( idx_M_t<Model>& idx );
template<typename Model> operation cconj_timerev_M_m( idx_M_t<Model>& idx );
template<typename Model> operation timerev_M_sc( idx_M_t<Model>& idx );
template<typename Model> operation timerev_M_m( idx_M_t<Model>& idx );
template<typename Model> operation hmirror_timerev_M_d( idx_M_t<Model>& idx );
