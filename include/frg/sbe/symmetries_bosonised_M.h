#pragma once

#include <frg/sbe/state_gf_types_bosonised_M.h>
#include <frg/sbe/symmetries.h>

template <typename Model>
class SymmetriesfRGSBEBosonised_M : public Symmetries1lfRGSBE<Model>
{
 public:
    static void Init();
    
    static void CleanUp();

    // index equivalence classes for w_M sc/d/m
    static symmetry_grp_t<dcomplex,4>* &IdxEquivClasses_w_M_sc_Ptr(){
	static symmetry_grp_t<dcomplex,4> *idx_equiv_classes_w_M_sc_ptr(nullptr);
	return idx_equiv_classes_w_M_sc_ptr;
    }
    static symmetry_grp_t<dcomplex,4>* &IdxEquivClasses_w_M_d_Ptr(){
	static symmetry_grp_t<dcomplex,4> *idx_equiv_classes_w_M_d_ptr(nullptr);
	return idx_equiv_classes_w_M_d_ptr;
    }
    static symmetry_grp_t<dcomplex,4>* &IdxEquivClasses_w_M_m_Ptr(){
	static symmetry_grp_t<dcomplex,4> *idx_equiv_classes_w_M_m_ptr(nullptr);
	return idx_equiv_classes_w_M_m_ptr;
    }
    
    // index equivalence classes for lambda_M sc/d/m
    static symmetry_grp_t<dcomplex,5>* &IdxEquivClasses_lambda_M_sc_Ptr(){
	static symmetry_grp_t<dcomplex,5> *idx_equiv_classes_lambda_M_sc_ptr(nullptr);
	return idx_equiv_classes_lambda_M_sc_ptr;
    }
    static symmetry_grp_t<dcomplex,5>* &IdxEquivClasses_lambda_M_d_Ptr(){
	static symmetry_grp_t<dcomplex,5> *idx_equiv_classes_lambda_M_d_ptr(nullptr);
	return idx_equiv_classes_lambda_M_d_ptr;
    }
    static symmetry_grp_t<dcomplex,5>* &IdxEquivClasses_lambda_M_m_Ptr(){
	static symmetry_grp_t<dcomplex,5> *idx_equiv_classes_lambda_M_m_ptr(nullptr);
	return idx_equiv_classes_lambda_M_m_ptr;
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

    // Symmetry group actions for w_M sc/d/m
    static std::vector<boost::function<operation (idx_w_M_t<Model>&)> > &GroupActions_w_M_sc(){
	static std::vector<boost::function<operation (idx_w_M_t<Model>&)> > group_actions_w_M_sc({});
	return group_actions_w_M_sc;
    }
    static std::vector<boost::function<operation (idx_w_M_t<Model>&)> > &GroupActions_w_M_d(){
	static std::vector<boost::function<operation (idx_w_M_t<Model>&)> > group_actions_w_M_d({});
	return group_actions_w_M_d;
    }
    static std::vector<boost::function<operation (idx_w_M_t<Model>&)> > &GroupActions_w_M_m(){
	static std::vector<boost::function<operation (idx_w_M_t<Model>&)> > group_actions_w_M_m({});
	return group_actions_w_M_m;
    }

    // Symmetry group actions for lambda_M sc/d/m
    static std::vector<boost::function<operation (idx_lambda_M_t<Model>&)> > &GroupActions_lambda_M_sc(){
	static std::vector<boost::function<operation (idx_lambda_M_t<Model>&)> > group_actions_lambda_M_sc({});
	return group_actions_lambda_M_sc;
    }
    static std::vector<boost::function<operation (idx_lambda_M_t<Model>&)> > &GroupActions_lambda_M_d(){
	static std::vector<boost::function<operation (idx_lambda_M_t<Model>&)> > group_actions_lambda_M_d({});
	return group_actions_lambda_M_d;
    }
    static std::vector<boost::function<operation (idx_lambda_M_t<Model>&)> > &GroupActions_lambda_M_m(){
	static std::vector<boost::function<operation (idx_lambda_M_t<Model>&)> > group_actions_lambda_M_m({});
	return group_actions_lambda_M_m;
    }

    // Lattice Symmetries
    static std::function<operation(idx_w_M_t<Model>& )> MakeLatticeSymmetryMap_w_M(unsigned point_group_element_idx);

    static std::function<operation(idx_lambda_M_t<Model>& )> MakeLatticeSymmetryMap_lambda_M(unsigned point_group_element_idx);
};
