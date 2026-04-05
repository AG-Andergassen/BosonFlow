
/*******************************************************************************************//** @file
 *  		
 * 	file: 		symmetry_group.h
 * 	contents:  	Defines a class that allows to initialize elements of the gf container
 * 			using symmetry relations in an efficient way
 * 
 ****************************************************************************************************/

#pragma once

#include <gf.h>
#include <complex>
#include <iostream>
/**
 *	Set of possible operations after symmetry operation
 */
class operation : public std::pair<bool, bool>
{
   public:
      using base_t = std::pair<bool, bool>; 
      using base_t::first; 
      using base_t::second; 

      operation( bool first_, bool second_ ) ///< Constructor taking two bools as argument
      {
	 first = first_ ;	// Apply minus sign?
	 second = second_;	// Apply complex conjugation?
      }

      operation  operator*( const operation& b )	///< Overload multiplication operator for successive application of operations
      {
	 return operation(  (*this ).first xor b.first, (*this ).second xor b.second ) ;
      }

      std::complex<double> operator()( const std::complex<double> val )		///< Apply operation
      {
	 if( first )
	 {
	    if( second )
	       return -conj( val ); 
	    return -val; 
	 }

	 if( second )
	    return conj( val ); 

	 return val; 
      }
      
      Eigen::Matrix<dcomplex, Eigen::Dynamic, 1> operator()( const Eigen::Matrix<dcomplex, Eigen::Dynamic, 1> val )		///< Apply operation

      {
	 if( first )
	 {
	    if( second )
	       return -val.conjugate(); 
	    return -val; 
	 }

	 if( second )
	    return val.conjugate() ; 

	 return val; 
      }

   private:
};

/**
 *	Elements of vertex tensor. States index of independent coupling array and possible operations on it
 */
struct symm_idx_t
{
   public:
      unsigned long int idx;	///< Index in the raw data array of container
      operation oper; 	///< Possible operations that relate two tensor elements. ( sign change, complex conjugation )

      symm_idx_t():
	 idx( 0 ), oper( false, false )
   {}

      symm_idx_t( int idx_, bool first_ = false, bool second_ = false ):
	 idx( idx_ ), oper( first_, second_ )
   {}

      symm_idx_t( int idx_, operation oper_ ):
	 idx( idx_ ), oper( oper_ )
   {}
};

template< typename elem_t_, unsigned int ndims_ >
class symmetry_grp_t
{
   public:
      using dcomplex = std::complex<double>; 

      using idx_t = idx_obj_t<ndims_>; 
      using type = symmetry_grp_t<elem_t_, ndims_>; 
      using gf_t = gf<elem_t_, ndims_>; 
      using symm_func_t = boost::function<operation ( idx_t& idx )>; 
      using init_func_t = boost::function<elem_t_ ( const idx_t& idx )>; 
      using symm_class_t = std::vector< symm_idx_t >;

      static constexpr unsigned int ndims = ndims_;                     ///< The number of dimensions        

      void init( gf<elem_t_,ndims>& gf_obj, init_func_t init_func )
      {
#pragma omp parallel for schedule( static )
	 for( unsigned int i = 0; i < symm_classes.size(); ++i )
	 {
	    elem_t_ val = init_func( gf_obj.get_idx( symm_classes[i][0].idx ) ); 
	    for( auto symm_idx : symm_classes[i] )
	       gf_obj( symm_idx.idx ) = symm_idx.oper( val ); 
	 }

      }
      
      void init_at_sampling_indices( gf<elem_t_,ndims>& gf_obj, init_func_t init_func )
      {
#pragma omp parallel for schedule( static )
	  for( unsigned int i = 0; i <  sampling_indices.size(); i ++)
	      {
		  elem_t_ val = init_func( gf_obj.get_idx( symm_classes[sampling_indices[i] ][0].idx ) ); 
		  for( auto symm_idx : symm_classes[sampling_indices[i] ] )
		      gf_obj( symm_idx.idx ) = symm_idx.oper( val ); 
	      }
      }
      

      void init_batched( gf<elem_t_,ndims>& gf_obj, init_func_t init_func )
      {
#pragma omp for schedule(nonmonotonic:dynamic) nowait
          for( unsigned int i = 0; i < symm_classes.size(); ++i )
          {
              elem_t_ val = init_func( gf_obj.get_idx( symm_classes[i][0].idx ) );
              for( auto symm_idx : symm_classes[i] )
              gf_obj( symm_idx.idx ) = symm_idx.oper( val );
          }
       }


      void init_batched_at_sampling_indices(gf<elem_t_,ndims>& gf_obj, init_func_t init_func )
      {
#pragma omp for schedule(nonmonotonic:dynamic) nowait
          for( unsigned int i = 0; i < sampling_indices.size(); i++)
          {
              elem_t_ val = init_func( gf_obj.get_idx( symm_classes[sampling_indices[i] ][0].idx ) );
              for( auto symm_idx : symm_classes[sampling_indices[i] ] )
		  gf_obj( symm_idx.idx ) = symm_idx.oper( val );
          }
       }

      void init_batched_at_interpolating_indices(gf<elem_t_,ndims>& gf_obj, init_func_t init_func )
      {
#pragma omp for schedule(nonmonotonic:dynamic) nowait
          for( unsigned int i = 0; i < interpolating_indices.size(); i++)
          {
              elem_t_ val = init_func( gf_obj.get_idx( symm_classes[interpolating_indices[i] ][0].idx ) );
              for( auto symm_idx : symm_classes[interpolating_indices[i] ] )
		  gf_obj( symm_idx.idx ) = symm_idx.oper( val );
          }
       }


      void filter_classes( const gf_t& gf_obj ) // each element in the argument looks like {idxkey, idx_values}
      {
	  sampling_indices = {};
	  interpolating_indices = {};
	  interpolating_gf_indices = {};
	  //#pragma omp parallel for schedule( static ) // push_back is not thread safe
	  for( unsigned int i = 0; i < symm_classes.size(); ++i ) {
	      bool is_filtered = false;
	      for (auto symm_idx : symm_classes[i]){
		  for (unsigned j = 0; j < filters.size(); j++){
		      const unsigned idx_key = filters[j].first;
		      const std::vector<int> kept_idxes_components = filters[j].second;
		      auto result = std::find(kept_idxes_components.begin(), kept_idxes_components.end(), gf_obj.get_idx( symm_idx.idx )(idx_key));
		      if (result == kept_idxes_components.end())
			  break;
		      if (j+1 == filters.size()){ 
			  sampling_indices.push_back(i);
			  is_filtered = true;
			  goto exit_loop;
		      }
		  }
	      }
	  exit_loop:
	      if (! is_filtered){
		  interpolating_indices.push_back(i);
		  interpolating_gf_indices.push_back( gf_obj.get_idx(symm_classes[i][0].idx) );
	      }
	  }
      }

      symmetry_grp_t( const type& symm_grp ):
      symm_lst( symm_grp.symm_lst ), symm_classes( symm_grp.symm_classes ), sampling_indices( symm_grp.sampling_indices ), interpolating_indices( symm_grp.interpolating_indices ), interpolating_gf_indices( symm_grp.interpolating_gf_indices )
   {}

      symmetry_grp_t( type&& symm_grp ) noexcept:
      symm_lst( symm_grp.symm_lst ), symm_classes( std::move( symm_grp.symm_classes ) ), sampling_indices( std::move( symm_grp.sampling_indices ) ), interpolating_indices( std::move( symm_grp.interpolating_indices ) ), interpolating_gf_indices( std::move( symm_grp.interpolating_gf_indices ) )
   {}

      type& operator=( const type& symm_grp )
      {
	 symm_lst = symm_grp.symm_lst; 
	 symm_classes = symm_grp.symm_classes;
	 sampling_indices = symm_grp.sampling_indices;
	 interpolating_indices = symm_grp.interpolating_indices;
	 interpolating_gf_indices = symm_grp.interpolating_gf_indices;
	 return *this; 
      }

      type& operator=( type&& symm_grp ) noexcept
	  {
	      symm_lst = symm_grp.symm_lst; 
	      symm_classes.operator=( std::move( symm_grp.symm_classes ) );
	      sampling_indices.operator=( std::move( symm_grp.sampling_indices ) );  
	      interpolating_indices.operator=( std::move( symm_grp.interpolating_indices ) );  
	      interpolating_gf_indices.operator=( std::move( symm_grp.interpolating_gf_indices ) );  

	      return *this; 
	  }

 symmetry_grp_t( const gf_t& gf_obj, const std::vector< symm_func_t >& symm_lst_,  std::vector< std::pair<int, std::vector<int> > > filters_) : symmetry_grp_t( gf_obj, symm_lst_ )
	  {
	      filters = filters_;
	      // todo: make child class: sampled_symmetry_grp_t.init_at_sampled/interpolating should be virtual
	      filter_classes(gf_obj);
	  }

 symmetry_grp_t( const gf_t& gf_obj, const std::vector< symm_func_t >& symm_lst_ ):
      symm_lst(symm_lst_),
	  shape_arr( reinterpret_cast<const boost::array< boost::multi_array_types::size_type, ndims >& >( *(gf_obj.shape_arr) ) ), 
	  idx_bases( reinterpret_cast<const boost::array< boost::multi_array_types::index, ndims >& >( *(gf_obj.idx_bases) ) )
	      {
		  std::cout << " Initializing symmetry group for container of Length : " << ndims << std::endl; 
		  //std::cout << " Class (internal) name : " << typeid(this).name() << std::endl; 
		  std::cout << "   Shape: "; 
		  for( unsigned int i = 0; i < ndims; ++i )
		      std::cout << shape_arr[i] << " "; 
		  std::cout << "   Idx_bases: "; 
		  for( unsigned int i = 0; i < ndims; ++i )
		      std::cout << idx_bases[i] << " ";
		  
		  // For frequencies, consider going out of range with symmetry functions
		  
		  auto extended_idx_bases = idx_bases;
		  auto extended_shape_arr = shape_arr;

		  for( unsigned int i = 0; i < ndims; ++i ){
		      if( extended_idx_bases[i] != 0 ){
			  extended_idx_bases[i] *= 2; 
			  extended_shape_arr[i] *= 2; 
		      }
		  } 

		  // Create bool array to track checked elements
		  gf< bool, ndims > checked( extended_shape_arr );
		  checked.reindex( extended_idx_bases ); 
		  checked.init( []( const idx_t& idx ){ return false; } ); 

		  for( long unsigned int iter = 0; iter < gf_obj.num_elements(); ++iter )
		      {
			  idx_t idx( gf_obj.get_idx( iter ) ); 

			  if ( !(checked( idx )) )  								// if tensor object not yet related to any other
			      {
				  checked( idx ) = true;
				  std::vector< symm_idx_t > current_class { symm_idx_t(iter) }; 	// initialize new symmetry class
				  operation track_op( false, false );
				  iterate( idx, gf_obj, track_op, checked, current_class ); 		// start iterating on index 
				  symm_classes.push_back( current_class );				// Add current symmetry class to list
			      }
			  //	else
			  //	   std::cout << std::endl << " Indices " << idx(0) << idx(1) << idx(2) << std::endl; 
		      }
		  
		  //std::cout << std::endl << " Symmetry class size " <<  symm_classes.size() << "gf_obj num elements" << gf_obj.num_elements() << std::endl; 
		  std::cout << std::endl << " Symmetries reduction " << 1.0 * symm_classes.size() / gf_obj.num_elements() << std::endl; 
	      }
      
      std::vector<unsigned> interpolating_indices;
      std::vector<idx_t> interpolating_gf_indices;
      std::vector< std::pair<int, std::vector<int> > > filters;

      boost::array<long, ndims> idx_bases; 
      boost::array<size_t, ndims> shape_arr; 
      
   private:
      std::vector< symm_func_t > symm_lst; 
      std::vector< symm_class_t > symm_classes; 
      std::vector<unsigned> sampling_indices;


      void iterate( const idx_t& idx_it, const gf<elem_t_,ndims>& gf_obj,  const operation& track_op, gf<bool,ndims>& checked, std::vector< symm_idx_t >& current_class  )
      {
	 for( const auto& symm: symm_lst ) 				// iterate over list of all symmetries specified
	 {
	    idx_t idx = idx_it;					// copy idx
	    operation curr_op = symm( idx ) * track_op;		// apply symmetry operation and track operations applied

	    if( !checked( idx ) )				// if index not yet checked
	    { 
	       checked( idx ) = true;

	       if( gf_obj.is_valid( idx ) ) 
		  current_class.push_back( symm_idx_t( gf_obj.get_pos_1d( idx ), curr_op ) ); // if valid index, add to current symmetry class

	       iterate( idx, gf_obj, curr_op, checked, current_class );	// iterate further	
	    }
	 }
      }
}; 
