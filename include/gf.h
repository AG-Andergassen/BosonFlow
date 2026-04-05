
/******************************************************************************************//** @file
 *  		
 * 	file: 		gf.h
 * 	contents:	Definition of Green function container class based on boost::multi_array
 * 			and associated index objects
 * 
 ****************************************************************************************************/


#pragma once

#include <iostream>
#include <vector>
#include <complex>
#include <boost/multi_array.hpp>
#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/bool.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/and.hpp>
#include <boost/function.hpp>
#include <boost/bind.hpp>

using boost::multi_array_types::extent_range;
using boost::multi_array_types::extent_gen;
using boost::multi_array; 

   template <class T, std::size_t N>
std::ostream& operator<<( std::ostream& os, const std::array<T, N>& arr )
{
   std::copy(arr.cbegin(), arr.cend(), std::ostream_iterator<T>(os, " "));
   return os;
}
/********************* Index type ********************/

template< unsigned int ndims_ >
class idx_obj_t							///< Index type that will be associated to each container
{
   public:
      static constexpr unsigned int ndims = ndims_; 				///< The number of dimensions

      using arr_t = std::array<int, ndims>; 	///< Array type

      arr_t idx_arr; 				///< Data array

      template< typename integer_t >
	 int& operator()( integer_t val )		///< Accesss by name, e.g. idx(name_t::w)
	 {
	    return idx_arr[int(val)]; 
	 }

      template< typename integer_t >
	 int operator()( integer_t val) const	///< Accesss by name, e.g. idx(name_t::w)
	 {
	    return idx_arr[int(val)]; 
	 }

      bool operator==(const idx_obj_t &idx2) const
      {
	  return idx_arr == idx2.idx_arr;
      }

      inline int* data()			///< Pointer to first data element
      {
	 return idx_arr.data(); 
      }

      idx_obj_t( const arr_t& idx_arr_ ):
	 idx_arr( idx_arr_ )	 
   {}

      idx_obj_t( arr_t&& idx_arr_ ):
	 idx_arr( std::move(idx_arr_) )	 
   {}

      idx_obj_t():
	 idx_arr()	 
   {}

      friend std::ostream& operator<<( std::ostream& os, idx_obj_t idx )	///< Output idx_t object to ostream
      {
	 return os << idx.idx_arr;
      }

}; 

/********************* Container ********************/

template< typename value_t_, unsigned int ndims_>
class gf: public boost::multi_array<value_t_, ndims_>
{
   public:

      // -------  Typedefs
      using value_t = value_t_; 						///< Type of container values
      using base_t = boost::multi_array<value_t, ndims_>; 			///< Type of base class
      using extents_t = boost::detail::multi_array::extent_gen<ndims_>; 	///< Type of extents object passed at contruction
      using idx_t = idx_obj_t<ndims_>; 
      using type = gf< value_t_, ndims_ >;
      using init_func_t = boost::function<value_t ( const idx_t& idx )>; 	///< Initalization function type that returns value_t for a given idx_t object
      
      // -------- Member variables
      const typename base_t::size_type* shape_arr; 			///< Array containing the extent for each dimension
      const typename base_t::index* stride_arr; 			///< Array containing the stride for each dimension
      const typename base_t::index* idx_bases; 				///< Array containing the base index for each dimension
      static constexpr unsigned int ndims = ndims_; 			///< The number of dimensions


      inline int get_pos_1d( const idx_t& idx ) const			///< For a given idx_t object, returns the corresponding position in a 1d array
      {
	 return get_pos_1d( idx.idx_arr ); 
      }

      int get_pos_1d( const std::array<int, ndims>& idx_arr ) const	///< For a given index array, returns the corresponding position in a 1d array
      {
	 int val = 0; 
	 for( unsigned int i = 0; i < ndims; ++i )
	 {
	    val += stride_arr[i]*( idx_arr[i] - idx_bases[i] ); 
	 }
	 return val; 
      }

      idx_t get_idx( int pos_1d ) const					///< For a given index in a 1d array, returns the corresponding idx_t object
      {
	 assert( pos_1d < base_t::num_elements() ); 
	 std::array<int, ndims> idx_arr; 

	 for( unsigned int i = 0; i < ndims; ++i )
	 {
	    idx_arr[i] = idx_bases[i] + pos_1d / stride_arr[i];
	    pos_1d -= ( idx_arr[i] - idx_bases[i] ) * stride_arr[i]; 
	 }

	 return idx_t( idx_arr ); 
      }

      void fill_idx_lst( std::vector<idx_t>& idx_lst )	const		///< Fills a std::vector<idx_t> with all possible sets of indeces
      {
	 idx_lst.resize( base_t::num_elements() ); 

	 for( unsigned int i = 0; i < base_t::num_elements(); ++i )
	    idx_lst[i] = get_idx( i ); 
      }

      void init( init_func_t init_func )				///< Initializes values with a given initialization function
      {
#pragma omp parallel for schedule( static )
	 for( unsigned int i = 0; i < base_t::num_elements(); ++i )
	    operator()( i ) = init_func( get_idx( i ) ); 
      }

      void init_at_sampling_indices(init_func_t init_func, const std::vector<unsigned> &sampling_indices)
      {
#pragma omp parallel for schedule( static )
	  for( unsigned i = 0; i <  sampling_indices.size(); i ++ )
	      operator()( sampling_indices[i] ) = init_func( get_idx( sampling_indices[i] ) ); 	  
      }

      std::pair<std::vector<unsigned>, std::vector<unsigned> > filter_indices( std::vector< std::pair<unsigned, std::vector<int> > > filters) // each element looks like {idxkey, idx_values}
      {
	  std::vector<unsigned> filtered_indices;
	  std::vector<unsigned> residual_indices;
#pragma omp parallel for schedule( static )
	  for( unsigned int i = 0; i < base_t::num_elements(); ++i ) {
	      bool is_filtered = false;
	      for (auto f : filters){
		  const unsigned idx_key = f.first;
		  const std::vector<int> kept_idxes_components = f.second;
		  auto result = std::find(kept_idxes_components.begin(), kept_idxes_components.end(), get_idx(i)(idx_key));
		  if (result != kept_idxes_components.end()){
		      filtered_indices.push_back(i);
		      is_filtered = true;
		      break;
		  } 
	      }
	      if (! is_filtered)
		  residual_indices.push_back(i);
	  }

	  return std::make_pair(filtered_indices, residual_indices);
      }

      bool is_valid( const idx_t& idx ) const				///< Return value of container for given index object
      {
	 for( unsigned int i = 0; i < ndims; ++i )
	    if( idx(i) >= idx_bases[i] + static_cast<int>(shape_arr[i]) || idx(i) < idx_bases[i] ) // fails without cast as shape_arr contains unsigned
	       return false;
	 return true; 
      }

      value_t& operator()( const idx_t& idx )				///< Return value of container for given index object
      {
	 return data_ptr[ get_pos_1d( idx.idx_arr ) ]; 
      }

      const value_t& operator()( const idx_t& idx ) const		///< Return value of container for given index object
      {
	 return data_ptr[ get_pos_1d( idx.idx_arr ) ]; 
      }

      value_t& operator()( int i )					///< Return value of container for given index object
      {
	 return data_ptr[ i ]; 
      }

      const value_t& operator()( int i ) const				///< Return value of container for given index object
      {
	 return data_ptr[ i ]; 
      }

      gf( const type& gf_obj ):
	 base_t( static_cast<const base_t&>(gf_obj) ), shape_arr( base_t::shape() ), stride_arr( base_t::strides() ), idx_bases( base_t::index_bases() ), data_ptr( base_t::data() ) 
   {}; 

      gf( type&& gf_obj ):
	 base_t( std::move( static_cast<base_t&>(gf_obj) ) ), shape_arr( base_t::shape() ), stride_arr( base_t::strides() ), idx_bases( base_t::index_bases() ), data_ptr( base_t::data() ) 
   {};

      ~gf(){
      }

      type& operator=( const type& gf_obj )
      {
	 base_t::operator=( static_cast<const base_t&>(gf_obj) ); 
	 shape_arr = base_t::shape(); 
	 stride_arr = base_t::strides(); 
	 idx_bases = base_t::index_bases(); 
	 data_ptr = base_t::data(); 
	 return *this; 
      } 

      type& operator=( type&& gf_obj )
      {
	 base_t::operator=( std::move( static_cast< base_t& >(gf_obj) ) ); 
	 shape_arr = base_t::shape(); 
	 stride_arr = base_t::strides(); 
	 idx_bases = base_t::index_bases(); 
	 data_ptr = base_t::data(); 
	 return *this; 
      } 

      gf( extents_t idx_ranges ):
	 base_t( idx_ranges ), shape_arr( base_t::shape() ), stride_arr( base_t::strides() ), idx_bases( base_t::index_bases() ), data_ptr( base_t::data() ) 
   {};

      gf( const boost::array<size_t, ndims>& idx_range_arr ):
	 base_t( idx_range_arr ), shape_arr( base_t::shape() ), stride_arr( base_t::strides() ), idx_bases( base_t::index_bases() ), data_ptr( base_t::data() ) 
   {};

      gf( extents_t idx_ranges, init_func_t init_func ):
	 base_t( idx_ranges ), shape_arr( base_t::shape() ), stride_arr( base_t::strides() ), idx_bases( base_t::index_bases() ), data_ptr( base_t::data() ) 
   {
      (*this).init( init_func ); 
   };

   private:
      value_t* data_ptr; 
};

struct gf_conversion_tester
{
   template <typename value_t_, unsigned int ndims_>
      gf_conversion_tester(const gf< value_t_, ndims_ > &);
};

template <class test_class>
struct is_instance_of_gf: boost::is_convertible<test_class,gf_conversion_tester>
{}; 

//template< typename gf_t_ >
//struct is_gf : std::false_type { };  

//template< typename value_t_, unsigned int ndims_ >
//struct is_gf< gf< value_t_, ndims_ > > : std::true_type { };  

template <class test_class>
struct is_scalar : std::false_type
{}; 

template<> struct is_scalar<double> : std::true_type { };  
template<> struct is_scalar<int> : std::true_type { };  
template<> struct is_scalar<std::complex<double>> : std::true_type { };  

// --------------- Abs and Norm

   template< typename gf_t_ >
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_ >::type abs( const gf_t_& lhs )
{
   using std::abs; 

   gf_t_ res( lhs ); 
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      res(i) = abs(res(i)); 
   return res; 
}

   template< typename gf_t_ >
typename boost::enable_if< is_instance_of_gf< gf_t_ >, double >::type norm( const gf_t_& lhs )
{
   gf_t_ gf_abs = abs( lhs ); 
   using value_t = typename gf_t_::value_t; 
   return std::real( *( std::max_element( gf_abs.data(), gf_abs.data() + gf_abs.num_elements(), []( value_t& a, value_t& b )->bool{ return std::abs(a) < std::abs(b); } ) ) ); 
}

// -------------- OPERATORS 

// ------ Single gf operations

/// -gf
   template< typename gf_t_ >
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_ >::type operator-( const gf_t_& lhs )
{
   gf_t_ res( lhs ); 
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      res(i) = -res(i); 
   return res; 
}

// ------ Two gf operations

/// gf1 += gf2
   template< typename gf_t_>
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_& >::type operator+=( gf_t_& lhs, const gf_t_& rhs )
{
   assert( lhs.num_elements() == rhs.num_elements() && " Adding gf's of different size " ); 
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      lhs(i) += rhs(i); 
   return lhs; 
}

/// gf1 + gf2
   template< typename gf_t_>
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_ >::type operator+( const gf_t_& lhs, const gf_t_& rhs )
{
   gf_t_ res( lhs ); 
   res += rhs; 
   return res; 
}

/// gf1 -= gf2
   template< typename gf_t_>
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_ >::type operator-=( gf_t_& lhs, const gf_t_& rhs )
{
   assert( lhs.num_elements() == rhs.num_elements() && " Adding gf's of different size " ); 
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      lhs(i) -= rhs(i); 
   return lhs; 
}

/// gf1 - gf2
   template< typename gf_t_>
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_ >::type operator-( const gf_t_& lhs, const gf_t_& rhs )
{
   gf_t_ res( lhs ); 
   res -= rhs; 
   return res; 
}

/// gf1 *= gf2
   template< typename gf_t_>
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_ >::type operator*=( gf_t_& lhs, const gf_t_& rhs )
{
   assert( lhs.num_elements() == rhs.num_elements() && " Adding gf's of different size " ); 
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      lhs(i) *= rhs(i); 
   return lhs; 
}

/// gf1 * gf2
   template< typename gf_t_>
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_ >::type operator*( const gf_t_& lhs, const gf_t_& rhs )
{
   gf_t_ res( lhs ); 
   res *= rhs;  
   return res; 
}

/// gf1 /= gf2
   template< typename gf_t_>
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_ >::type operator/=( gf_t_& lhs, const gf_t_& rhs )
{
   assert( lhs.num_elements() == rhs.num_elements() && " Adding gf's of different size " ); 
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      lhs(i) /= rhs(i); 
   return lhs; 
}

/// gf1 / gf2
   template< typename gf_t_>
typename boost::enable_if< is_instance_of_gf< gf_t_ >, gf_t_ >::type operator/( const gf_t_& lhs, const gf_t_& rhs )
{
   gf_t_ res( lhs ); 
   res /= rhs; 
   return res; 
}


// ------------ Scalar operations

/// gf += s
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_& >::type 
operator+=( gf_t_& lhs, const scalar_t_& rhs )
{
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      lhs(i) += rhs; 
   return lhs; 
}

/// gf + s
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_ >::type 
operator+( const gf_t_& lhs, const scalar_t_& rhs )
{
   gf_t_ res( lhs ); 
   res += rhs; 
   return res; 
}

/// s + gf
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_ >::type 
operator+( const scalar_t_& lhs, const gf_t_& rhs )
{
   gf_t_ res( rhs ); 
   res += lhs; 
   return res; 
}

/// gf -= s
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_& >::type 
operator-=( gf_t_& lhs, const scalar_t_& rhs )
{
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      lhs(i) -= rhs; 
   return lhs; 
}

/// gf - s
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_ >::type 
operator-( const gf_t_& lhs, const scalar_t_& rhs )
{
   gf_t_ res( lhs ); 
   res -= rhs; 
   return res; 
}

/// s - gf
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_ >::type 
operator-( const scalar_t_& lhs,  const gf_t_& rhs )
{
   gf_t_ res( -rhs ); 
   res += lhs; 
   return res; 
}

/// gf *= s
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_& >::type 
operator*=( gf_t_& lhs, const scalar_t_& rhs )
{
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      lhs(i) *= rhs; 
   return lhs; 
}

/// gf * s
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_ >::type 
operator*( const gf_t_& lhs, const scalar_t_& rhs )
{
   gf_t_ res( lhs ); 
   res *= rhs; 
   return res; 
}

/// s * gf
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_ >::type 
operator*( const scalar_t_& lhs, const gf_t_& rhs )
{
   gf_t_ res( rhs ); 
   res *= lhs; 
   return res; 
}

/// gf /= s
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_& >::type 
operator/=( gf_t_& lhs, const scalar_t_& rhs )
{
   for( unsigned int i = 0; i < lhs.num_elements(); ++i )
      lhs(i) /= rhs; 
   return lhs; 
}

/// gf / s
   template< typename gf_t_, typename scalar_t_ >
typename boost::enable_if< boost::mpl::and_< is_instance_of_gf< gf_t_ >, is_scalar< scalar_t_> >, gf_t_ >::type 
operator/( const gf_t_& lhs, const scalar_t_& rhs )
{
   gf_t_ res( lhs ); 
   res /= rhs; 
   return res; 
}

/// ostream << gf
   template< typename gf_t_ >
typename boost::enable_if< is_instance_of_gf< gf_t_ >, std::ostream& >::type 
operator<<( std::ostream& lhs, const gf_t_& rhs )
{
   for( unsigned int i = 0; i < rhs.num_elements(); ++i )
      lhs << rhs(i) << " \t "; //std::endl; 
      //lhs << rhs.get_idx(i) << " : " << rhs(i) << " ;\t "; //std::endl; 
   return lhs; 
}
