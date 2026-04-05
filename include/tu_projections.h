

/*******************************************************************************************//** @file
 *  		
 * 	file: 		tu_projections.h
 * 	contents:  	Definition of projection-related quantities 
 * 
****************************************************************************************************/
#pragma once

#include <base_gf_types.h>


//using proj_matrix_t = gf_proj_matrix_t;
//using vert_bare_ff_t = gf_vert_bare_ff_t;


template <typename Model> 
class TUProjections
{
 public:
    static void Init();

    static gf_proj_matrix_t<Model> &Matrix_xph_to_ph()
    {
	static gf_proj_matrix_t<Model> matrix_xph_to_ph;
	return matrix_xph_to_ph;
    }

    static gf_proj_matrix_t<Model> &Matrix_xph_to_pp()
    {
	static gf_proj_matrix_t<Model> matrix_xph_to_pp;
	return matrix_xph_to_pp;
    }

    static gf_proj_matrix_t<Model> &Matrix_ph_to_xph()
    {
	static gf_proj_matrix_t<Model> matrix_ph_to_xph;
	return matrix_ph_to_xph;
    }

    static gf_proj_matrix_t<Model> &Matrix_ph_to_pp()
    {
	static gf_proj_matrix_t<Model> matrix_ph_to_pp;
	return matrix_ph_to_pp;
    }

    static gf_proj_matrix_t<Model> &Matrix_pp_to_ph()
    {
	static gf_proj_matrix_t<Model> matrix_pp_to_ph;
	return matrix_pp_to_ph;
    }

    static gf_proj_matrix_t<Model> &Matrix_pp_to_xph()
    {
	static gf_proj_matrix_t<Model> matrix_pp_to_xph;
	return matrix_pp_to_xph;
    }

    static gf_vert_bare_ff_t<Model> &VertexBare()
    {
	static gf_vert_bare_ff_t<Model> vertex_4pt_bare;
	return vertex_4pt_bare;
    }

 private:
    static void CalculateDiagrammaticProjections();
    
    static delta_function_vector_t<Model::dim> GetFormFactor(const unsigned m);

    static dcomplex CalculateProjectionElement_xph_to_pp( const idx_proj_matrix_t<Model>& idx );
    static dcomplex CalculateProjectionElement_ph_to_pp( const idx_proj_matrix_t<Model>& idx );
    static dcomplex CalculateProjectionElement_pp_to_xph( const idx_proj_matrix_t<Model>& idx );
    static dcomplex CalculateProjectionElement_pp_to_ph( const idx_proj_matrix_t<Model>& idx );
    static dcomplex CalculateProjectionElement_ph_to_xph( const idx_proj_matrix_t<Model>& idx );
    static dcomplex CalculateProjectionElement_xph_to_ph( const idx_proj_matrix_t<Model>& idx );

    static dcomplex CalculateProjectionVertexBare(const idx_vert_bare_ff_t<Model>& idx);

};
