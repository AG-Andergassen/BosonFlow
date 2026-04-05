
#pragma once

#include <frg/sbe/state_gf_types.h>
#include <symmetry_group.h>
#include <models/concrete_available_models.h>

template <typename Model>
class SymmetriesfRGCommon
{
 public:
    static void Init();
    
    static void CleanUp();

    // index equivalence classes for bubbles
    static symmetry_grp_t<MatPatch, 8>* &IdxEquivClasses_bubblemat_pp_Ptr(){
	static symmetry_grp_t<MatPatch, 8> *idx_equiv_classes_bubblemat_pp_ptr(nullptr);
	return idx_equiv_classes_bubblemat_pp_ptr;
    }

    static symmetry_grp_t<MatPatch, 8>* &IdxEquivClasses_bubblemat_ph_Ptr(){
	static symmetry_grp_t<MatPatch, 8> *idx_equiv_classes_bubblemat_ph_ptr(nullptr);
	return idx_equiv_classes_bubblemat_ph_ptr;
    }

    // index equivalence classes for sig
    static symmetry_grp_t<dcomplex, 4>* &IdxEquivClasses_sig_Ptr(){
	static symmetry_grp_t<dcomplex, 4> *idx_equiv_classes_sig_ptr(nullptr);
	return idx_equiv_classes_sig_ptr;
    }

    // index equivalence classes for the MatPatch-valued sig
    static symmetry_grp_t<MatPatch,3>* &IdxEquivClasses_sigkmat_Ptr(){
	static symmetry_grp_t<MatPatch,3> *idx_equiv_classes_sigkmat_ptr(nullptr);
	return idx_equiv_classes_sigkmat_ptr;
    }

    // index equivalence classes for the Green function in real space
    static symmetry_grp_t<MatReal,3>* &IdxEquivClasses_Gvec_real_Ptr(){
	static symmetry_grp_t<MatReal,3> *idx_equiv_classes_Gvec_real_ptr(nullptr);
	return idx_equiv_classes_Gvec_real_ptr;
    }

    // index equivalence classes for the projection matrices
    static symmetry_grp_t<dcomplex, 6>* &IdxEquivClasses_projection_matrix_Ptr(){
	static symmetry_grp_t<dcomplex, 6> *idx_equiv_classes_projection_matrix_ptr(nullptr);
	return idx_equiv_classes_projection_matrix_ptr;
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


    // Symmetry group actions for the bubble
    static std::vector<boost::function<operation (idx_bubble_mat_t<Model>&)> > &GroupActions_bubblemat(){
	static std::vector<boost::function<operation (idx_bubble_mat_t<Model>&)> > group_actions_bubblemat({});
	return group_actions_bubblemat;
    }

    // Symmetry group actions for sig
    static std::vector<boost::function<operation (idx_1p_t<Model>&)> > &GroupActions_sig(){
	static std::vector<boost::function<operation (idx_1p_t<Model>&)> > group_actions_sig({});
	return group_actions_sig;
    }

    // Symmetry group actions for the MatPatch-valued sig
    static std::vector<boost::function<operation (idx_Sig_kMat_t<Model>&)> > &GroupActions_sigkmat(){
	static std::vector<boost::function<operation (idx_Sig_kMat_t<Model>&)> > group_actions_sigkmat({});
	return group_actions_sigkmat;
    }

    // Symmetry group actions for the Green function in real space
    static std::vector<boost::function<operation (idx_Sig_kMat_t<Model>&)> > &GroupActions_Gvec_real(){
	static std::vector<boost::function<operation (idx_Sig_kMat_t<Model>&)> > group_actions_Gvec_real({});
	return group_actions_Gvec_real;
    }

    // Symmetry group actions for the TU projection matrices
    static std::vector<boost::function<operation (idx_proj_matrix_t<Model>&)> > &GroupActions_projection_matrix(){
	static std::vector<boost::function<operation (idx_proj_matrix_t<Model>&)> > group_actions_projection_matrix({});
	return group_actions_projection_matrix;
    }

    // Lattice Symmetries
    static std::function<operation(idx_1p_t<Model>& )> MakeLatticeSymmetryMap_sig(unsigned point_group_element_idx);

    static std::function<operation(idx_proj_matrix_t<Model>& )> MakeLatticeSymmetryMap_projection_matrix(unsigned point_group_element_idx);
};


