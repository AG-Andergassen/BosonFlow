
/******************************************************************************************//** @file
 *  		
 * 	file: 		tu_projections.cpp
 * 	contents:  	See tu_projections.h
 * 
 ****************************************************************************************************/

#include <cmath>
#include <models/concrete_available_models.h>
#include <tu_projections.h>
#include <frg/symmetries_common.h>

template <typename Model> 
void TUProjections<Model>::Init( )
{
    std::cout << "Calculating TU projection matrices.." << std::endl;
    CalculateDiagrammaticProjections();
    std::cout << "Done calculating TU projection matrices." << std::endl;
}


template <typename Model>
void TUProjections<Model>::CalculateDiagrammaticProjections()
{

#pragma omp parallel default(shared)
    {
	// todo: split symmetries of projections away from fRG directory
	// initialise the main gf objects for projections
	SymmetriesfRGCommon<Model>::IdxEquivClasses_projection_matrix_Ptr()->init_batched(Matrix_xph_to_ph(), TUProjections<Model>::CalculateProjectionElement_xph_to_ph);
	SymmetriesfRGCommon<Model>::IdxEquivClasses_projection_matrix_Ptr()->init_batched(Matrix_xph_to_pp(), TUProjections<Model>::CalculateProjectionElement_xph_to_pp);
	SymmetriesfRGCommon<Model>::IdxEquivClasses_projection_matrix_Ptr()->init_batched(Matrix_ph_to_xph(), TUProjections<Model>::CalculateProjectionElement_ph_to_xph);
	SymmetriesfRGCommon<Model>::IdxEquivClasses_projection_matrix_Ptr()->init_batched(Matrix_ph_to_pp(), TUProjections<Model>::CalculateProjectionElement_ph_to_pp);
	SymmetriesfRGCommon<Model>::IdxEquivClasses_projection_matrix_Ptr()->init_batched(Matrix_pp_to_ph(), TUProjections<Model>::CalculateProjectionElement_pp_to_ph);
	SymmetriesfRGCommon<Model>::IdxEquivClasses_projection_matrix_Ptr()->init_batched(Matrix_pp_to_xph(), TUProjections<Model>::CalculateProjectionElement_pp_to_xph);
    }

    VertexBare().init(TUProjections<Model>::CalculateProjectionVertexBare);
}


template <typename Model> 
dcomplex TUProjections<Model>::CalculateProjectionElement_xph_to_pp( const idx_proj_matrix_t<Model>& idx )
{
    dcomplex val( 0.0, 0.0 );

    unsigned K_idx = idx(0), Q_idx = idx(1); // match momenta indices
    unsigned m = idx(2), n = idx(3), mp = idx(4), np = idx(5); // match form factor indices

    coord_t<Model::dim> K = Model::MomentumGrid().ptr->m_points.at(K_idx); // refined lattice
    coord_t<Model::dim> Q = Model::MomentumGrid().ptr->m_points.at(Q_idx); // refined lattice

    auto r1_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r2_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r3_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 2);

    for (coord_t<Model::dim> r1 : r1_points)
    for (coord_t<Model::dim> r2 : r2_points)
    for (coord_t<Model::dim> r3 : r3_points)
	val+= std::exp(I*((Q-K).dot(r3))) * 
	    std::conj(Model::GetFormFactor(m).evaluate(r1-r3))*
	    Model::GetFormFactor(n).evaluate(r2+r3)* 
	    Model::GetFormFactor(mp).evaluate(r1)*
	    std::conj(Model::GetFormFactor(np).evaluate(r2));	
    
    if(std::abs(val)<CHOP_ERR) 
	val=0;

    return val;

}

// A^{pp, ph}
template <typename Model> 
dcomplex TUProjections<Model>::CalculateProjectionElement_ph_to_pp( const idx_proj_matrix_t<Model>& idx )
{
    dcomplex val( 0.0, 0.0 );

    unsigned K_idx = idx(0), Q_idx = idx(1); // match momenta indices
    unsigned m = idx(2), n = idx(3), mp = idx(4), np = idx(5); // match form factor indices
   
    coord_t<Model::dim> K = Model::MomentumGrid().ptr->m_points.at(K_idx); // refined lattice
    coord_t<Model::dim> Q = Model::MomentumGrid().ptr->m_points.at(Q_idx); // refined lattice

    auto r1_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r2_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r3_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 2);


    for (coord_t<Model::dim> r1 : r1_points)
    for (coord_t<Model::dim> r2 : r2_points)
    for (coord_t<Model::dim> r3 : r3_points)
	val+= std::exp(I*(Q.dot(r3) + K.dot(r2))) * 
	    std::conj(Model::GetFormFactor(m).evaluate(r1-r3))*
	    Model::GetFormFactor(n).evaluate(-r2-r3)* 
	    Model::GetFormFactor(mp).evaluate(r1)*
	    std::conj(Model::GetFormFactor(np).evaluate(r2));	
    
    if(std::abs(val)<CHOP_ERR) 
	val=0;
    
    return val;
}

template <typename Model> 
dcomplex TUProjections<Model>::CalculateProjectionElement_pp_to_xph( const idx_proj_matrix_t<Model>& idx )
{
    dcomplex val( 0.0, 0.0 );
    
    unsigned K_idx = idx(0), Q_idx = idx(1); // match momenta indices
    unsigned m = idx(2), n = idx(3), mp = idx(4), np = idx(5); // match form factor indices
    
    coord_t<Model::dim> K = Model::MomentumGrid().ptr->m_points.at(K_idx); // refined lattice
    coord_t<Model::dim> Q = Model::MomentumGrid().ptr->m_points.at(Q_idx); // refined lattice

    
    auto r1_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r2_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r3_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 2);

    for (coord_t<Model::dim> r1 : r1_points)
    for (coord_t<Model::dim> r2 : r2_points)
    for (coord_t<Model::dim> r3 : r3_points)
	val+= std::exp(I*((Q-K).dot(r3))) * 
	    std::conj(Model::GetFormFactor(m).evaluate(r1+r3))*
	    Model::GetFormFactor(n).evaluate(r2-r3)* 
	    Model::GetFormFactor(mp).evaluate(r1)*
	    std::conj(Model::GetFormFactor(np).evaluate(r2));	
    
    if(std::abs(val)<CHOP_ERR) 
	val=0;

    return val;
}

template <typename Model> 
dcomplex TUProjections<Model>::CalculateProjectionElement_pp_to_ph( const idx_proj_matrix_t<Model>& idx )
{
    dcomplex val( 0.0, 0.0 );
    
    unsigned K_idx = idx(0), Q_idx = idx(1); // match momenta indices
    unsigned m = idx(2), n = idx(3), mp = idx(4), np = idx(5); // match form factor indices
    
    coord_t<Model::dim> K = Model::MomentumGrid().ptr->m_points.at(K_idx); // refined lattice
    coord_t<Model::dim> Q = Model::MomentumGrid().ptr->m_points.at(Q_idx); // refined lattice

    auto r1_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r2_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r3_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 2);

    for (coord_t<Model::dim> r1 : r1_points)
    for (coord_t<Model::dim> r2 : r2_points)
    for (coord_t<Model::dim> r3 : r3_points)
	val+= std::exp(I*(K.dot(r3) + Q.dot(r2 - r3))) * 
	    std::conj(Model::GetFormFactor(m).evaluate(r1-r3))*
	    Model::GetFormFactor(n).evaluate(r3-r2)* 
	    Model::GetFormFactor(mp).evaluate(r1)*
	    std::conj(Model::GetFormFactor(np).evaluate(r2));	
    
    if(std::abs(val)<CHOP_ERR) 
	val=0;

    return val;
}


template <typename Model> 
dcomplex TUProjections<Model>::CalculateProjectionElement_ph_to_xph( const idx_proj_matrix_t<Model>& idx )
{
    dcomplex val( 0.0, 0.0 );
    
    unsigned K_idx = idx(0), Q_idx = idx(1); // match momenta indices
    unsigned m = idx(2), n = idx(3), mp = idx(4), np = idx(5); // match form factor indices
    
    coord_t<Model::dim> K = Model::MomentumGrid().ptr->m_points.at(K_idx); // refined lattice
    coord_t<Model::dim> Q = Model::MomentumGrid().ptr->m_points.at(Q_idx); // refined lattice
    

    auto r1_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r2_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r3_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);

    for (coord_t<Model::dim> r1 : r1_points)
    for (coord_t<Model::dim> r2 : r2_points)
    for (coord_t<Model::dim> r3 : r3_points)
	val+= std::exp(I*(Q.dot(r3) + K.dot(r2))) * 
	    std::conj(Model::GetFormFactor(m).evaluate(r1-r2-r3))*
	    Model::GetFormFactor(n).evaluate(-r3)* 
	    Model::GetFormFactor(mp).evaluate(r1)*
	    std::conj(Model::GetFormFactor(np).evaluate(r2));	
    
    if(std::abs(val)<CHOP_ERR) 
	val=0;

    return val;
}

template <typename Model> 
dcomplex TUProjections<Model>::CalculateProjectionElement_xph_to_ph( const idx_proj_matrix_t<Model>& idx )
{
    dcomplex val( 0.0, 0.0 );
    
    unsigned K_idx = idx(0), Q_idx = idx(1); // match momenta indices
    unsigned m = idx(2), n = idx(3), mp = idx(4), np = idx(5); // match form factor indices
    
    coord_t<Model::dim> K = Model::MomentumGrid().ptr->m_points.at(K_idx); // refined lattice
    coord_t<Model::dim> Q = Model::MomentumGrid().ptr->m_points.at(Q_idx); // refined lattice
    
    
    auto r1_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r2_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);
    auto r3_points = Model::RealLattice().ptr->get_subset(Model::FormFactors().container_ptr->m_shells_count * 1);

    for (coord_t<Model::dim> r1 : r1_points)
    for (coord_t<Model::dim> r2 : r2_points)
    for (coord_t<Model::dim> r3 : r3_points)
	val+= std::exp(I*(Q.dot(r3) + K.dot(r2))) * 
	    std::conj(Model::GetFormFactor(m).evaluate(r1-r2-r3))*
	    Model::GetFormFactor(n).evaluate(-r3)* 
	    Model::GetFormFactor(mp).evaluate(r1)*
	    std::conj(Model::GetFormFactor(np).evaluate(r2));	
    
    if(std::abs(val)<CHOP_ERR) 
	val=0;

    return val;
}

//BARE VERTEX PROJECTION
template <typename Model> 
dcomplex TUProjections<Model>::CalculateProjectionVertexBare(const idx_vert_bare_ff_t<Model>& idx)
{
    unsigned n_in = idx(0), n_out = idx(1);
    unsigned s1_in = idx(2), s2_in = idx(3), s1_out = idx(4), s2_out = idx(5); // extract orbital indices

    // TODO: add proper momentum dependence (for models with extended interactions)
    if( n_in==0 && n_out==0)
	return Model::vertex_4pt_bare(-00, -00, -00, -00,
				  -00, -00, -00, -00,
                  s1_in,s2_in,s1_out,s2_out) * 4.0 * M_PI * M_PI;
    else
	return dcomplex(0.0,0.0);
}


// instantiate class for the particular models
#define Instantiate(MODEL) template class TUProjections<MODEL>;
WITH_MODELS(Instantiate)




