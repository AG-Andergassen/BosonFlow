
#pragma once

#include <frg/sbe/state_gf_types.h>
//#include <symmetry_group.h>
#include <frg/sbe/symmetries.h>

template <typename Model>
class SymmetriesmfRGSBE : public Symmetries1lfRGSBE<Model>
{
 public:
    static void Init();
    
    static void CleanUp();

    // index equivalence classes for lambda left corrections sc/d/m
    static symmetry_grp_t<dcomplex,4>* &IdxEquivClasses_lambda_sc_LeftCorrections_Ptr(){
	static symmetry_grp_t<dcomplex,4> *idx_equiv_classes_lambda_sc_left_corrections_ptr(nullptr);
	return idx_equiv_classes_lambda_sc_left_corrections_ptr;
    }
    static symmetry_grp_t<dcomplex,4>* &IdxEquivClasses_lambda_d_LeftCorrections_Ptr(){
	static symmetry_grp_t<dcomplex,4> *idx_equiv_classes_lambda_d_left_corrections_ptr(nullptr);
	return idx_equiv_classes_lambda_d_left_corrections_ptr;
    }
    static symmetry_grp_t<dcomplex,4>* &IdxEquivClasses_lambda_m_LeftCorrections_Ptr(){
	static symmetry_grp_t<dcomplex,4> *idx_equiv_classes_lambda_m_left_corrections_ptr(nullptr);
	return idx_equiv_classes_lambda_m_left_corrections_ptr;
    }


    // index equivalence classes for M left correction sc/d/m
    static symmetry_grp_t<dcomplex,6>* &IdxEquivClasses_M_sc_LeftCorrections_Ptr(){
	static symmetry_grp_t<dcomplex,6> *idx_equiv_classes_M_sc_left_corrections_ptr(nullptr);
	return idx_equiv_classes_M_sc_left_corrections_ptr;
    }
    static symmetry_grp_t<dcomplex,6>* &IdxEquivClasses_M_d_LeftCorrections_Ptr(){
	static symmetry_grp_t<dcomplex,6> *idx_equiv_classes_M_d_left_corrections_ptr(nullptr);
	return idx_equiv_classes_M_d_left_corrections_ptr;
    }
    static symmetry_grp_t<dcomplex,6>* &IdxEquivClasses_M_m_LeftCorrections_Ptr(){
	static symmetry_grp_t<dcomplex,6> *idx_equiv_classes_M_m_left_corrections_ptr(nullptr);
	return idx_equiv_classes_M_m_left_corrections_ptr;
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

    // Symmetry group actions for lambda left corrections sc/d/m
    static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > &GroupActions_lambda_sc_LeftCorrections(){
	static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > group_actions_lambda_sc_left_corrections({});
	return group_actions_lambda_sc_left_corrections;
    }
    static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > &GroupActions_lambda_d_LeftCorrections(){
	static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > group_actions_lambda_d_left_corrections({});
	return group_actions_lambda_d_left_corrections;
    }
    static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > &GroupActions_lambda_m_LeftCorrections(){
	static std::vector<boost::function<operation (idx_lambda_t<Model>&)> > group_actions_lambda_m_left_corrections({});
	return group_actions_lambda_m_left_corrections;
    }


    // Symmetry group actions for M left corrections sc/d/m
    static std::vector<boost::function<operation (idx_M_t<Model>&)> > &GroupActions_M_sc_LeftCorrections(){
	static std::vector<boost::function<operation (idx_M_t<Model>&)> > group_actions_M_sc_left_corrections({});
	return group_actions_M_sc_left_corrections;
    }
    static std::vector<boost::function<operation (idx_M_t<Model>&)> > &GroupActions_M_d_LeftCorrections(){
	static std::vector<boost::function<operation (idx_M_t<Model>&)> > group_actions_M_d_left_corrections({});
	return group_actions_M_d_left_corrections;
    }
    static std::vector<boost::function<operation (idx_M_t<Model>&)> > &GroupActions_M_m_LeftCorrections(){
	static std::vector<boost::function<operation (idx_M_t<Model>&)> > group_actions_M_m_left_corrections({});
	return group_actions_M_m_left_corrections;
    }

    // Lattice Symmetries
    //static std::function<operation(idx_M_t<Model>& )> MakeLatticeSymmetryMap_M(unsigned point_group_element_idx);
};
