
/************************************************************************************************//**
 *  		
 * 	file: 		arithmetic_tuple.h
 * 	contents:   	Wrapper around std::tuple to allow for element-wise arithmetic operations
 * 			Based on ReaK Library
 * 
 ****************************************************************************************************/


#pragma once

#include <tuple>

#include <boost/type_traits.hpp>
#include <boost/utility/enable_if.hpp>
#include <boost/mpl/size_t.hpp>
#include <boost/mpl/and.hpp>
#include <boost/mpl/or.hpp>
#include <boost/mpl/less.hpp>
#include <boost/mpl/greater.hpp>
#include <boost/mpl/equal_to.hpp>
#include <boost/mpl/prior.hpp>

namespace ReaK {

   using std::get;
   using std::abs; 
   using std::true_type;  
   using std::false_type;  

   /**
    * This meta-function computes a bool integral constant if the given type is an arithmetic-tuple.
    * \tparam Tuple The type to be tested as being an arithmetic-tuple or not.
    */
   template <typename Tuple>
      struct is_arithmetic_tuple : 
	 false_type { };

   /**
    * This meta-function computes an integral constant describing the size (or length) of an arithmetic-tuple.
    * \tparam Tuple The arithmetic-tuple type.
    */
   template <typename Tuple, bool val>
      struct arithmetic_tuple_size_impl
      {}; 

   template <typename Tuple>
      struct arithmetic_tuple_size_impl<Tuple, false> : 
      boost::mpl::size_t<0>
      {}; 

   template <typename Tuple>
      struct arithmetic_tuple_size_impl<Tuple, true> : 
      Tuple::tuple_size_t
      {}; 

   //template <typename Tuple, class Enable = void>
   //struct arithmetic_tuple_size : 
   //boost::mpl::size_t< 0 > { };

   /**
    * This class template is a simple wrapper of a tuple with the addition of arithmetic operators. 
    * This class is basically just a wrapper of the std::tuple class, and it provides 
    * a meta-programming interface that is equivalent to std::tuple and with the addition 
    * of the support for all the basic arithmetic operators, requiring, of course, that these 
    * arithmetic operators are also available on all the types contained in the tuple.
    * \tparam T The types contained in the arithmetic-tuple.
    */
   template <typename... T>
      class arithmetic_tuple : public std::tuple< T... > {
	 public:
	    typedef std::tuple< T... > base_t;

	    typedef boost::mpl::size_t< sizeof... (T) > tuple_size_t; 

	    constexpr arithmetic_tuple() : base_t() { };

	    explicit arithmetic_tuple( const double val ) : base_t() { };

	    template <typename... U>
	       explicit arithmetic_tuple(U&&... u) : base_t(std::forward<U>(u)...) { };


	    arithmetic_tuple( const arithmetic_tuple< T...>& Tuple ):    						
	       base_t( Tuple )							
	 {}

	    arithmetic_tuple( arithmetic_tuple< T...>&& Tuple ):    						
	       base_t( std::move(Tuple) )							
	 {}
        ~arithmetic_tuple(){
        }

	    arithmetic_tuple< T...>& operator=( const arithmetic_tuple< T... >& Tuple )					
	    {
            if( this == &Tuple)
                    return *this;
	       base_t::operator=( Tuple ); 	
	       return *this; 
	    } 							

	    arithmetic_tuple< T...>& operator=( arithmetic_tuple< T... >&& Tuple )						
	    {
           if( this == &Tuple )
                   return *this;
	       base_t::operator=( std::move( Tuple ) ); 
	       return *this; 
	    } 

      };

   /* Specialization, see general template docs. */
   template <typename... T>
      struct is_arithmetic_tuple< arithmetic_tuple< T... > > : 
      true_type {};

   struct arithmetic_tuple_conversion_tester
   {
      template <typename... T>
	 arithmetic_tuple_conversion_tester(const arithmetic_tuple<T...> &);
   };

   template <class test_class>
      struct is_instance_of_arithmetic_tuple: public boost::is_convertible<test_class,arithmetic_tuple_conversion_tester>
   {}; 

   template <typename Tuple>
      struct arithmetic_tuple_size: arithmetic_tuple_size_impl< Tuple, is_instance_of_arithmetic_tuple< Tuple >::value >
   {}; 

   template <class test_class>
      struct is_scalar : std::false_type
   {}; 

   template<> struct is_scalar<double> : std::true_type { };  
   template<> struct is_scalar<int> : std::true_type { };  
   template<> struct is_scalar<std::complex<double>> : std::true_type { };  

   /**
    * This function template can be used to create an arithmetic-tuple.
    * \tparam T The types contained in the arithmetic-tuple.
    * \param t The values that make up the arithmetic-tuple.
    * \return An arithmetic-tuple.
    */
   template <typename... T>
      inline 
      arithmetic_tuple< typename std::remove_reference<T>::type... > make_arithmetic_tuple(T&&... t) {
	 return arithmetic_tuple< typename std::remove_reference<T>::type... >(std::forward<T>(t)...);
      };

   namespace detail {

      /*****************************************************************************************
	Implementation details
       *****************************************************************************************/

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_addassign_impl( Tuple& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_addassign_impl( Tuple& lhs, const Tuple& rhs) {
	    tuple_addassign_impl< typename boost::mpl::prior<Idx>::type,Tuple >(lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(lhs) += get<boost::mpl::prior<Idx>::type::value>(rhs);
	 };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx,
	 boost::mpl::size_t<0>
	    >,
	 void >::type tuple_subassign_impl( Tuple& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx,
	 boost::mpl::size_t<0>
	    >,
	 void >::type tuple_subassign_impl( Tuple& lhs, const Tuple& rhs) {
	    tuple_subassign_impl< typename boost::mpl::prior<Idx>::type,Tuple>(lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(lhs) -= get<boost::mpl::prior<Idx>::type::value>(rhs);
	 };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_mulassign_impl( Tuple& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_mulassign_impl( Tuple& lhs, const Tuple& rhs) {
	    tuple_mulassign_impl<typename boost::mpl::prior<Idx>::type,Tuple>(lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(lhs) *= get<boost::mpl::prior<Idx>::type::value>(rhs);
	 };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_divassign_impl( Tuple& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_divassign_impl( Tuple& lhs, const Tuple& rhs) {
	    tuple_divassign_impl<typename boost::mpl::prior<Idx>::type,Tuple>(lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(lhs) /= get<boost::mpl::prior<Idx>::type::value>(rhs);
	 };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_smulassign_impl( Tuple& lhs, const Scalar& rhs) { };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_smulassign_impl( Tuple& lhs, const Scalar& rhs) {
	    tuple_smulassign_impl<typename boost::mpl::prior<Idx>::type,Tuple,Scalar>(lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(lhs) *= rhs;
	 };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_sdivassign_impl( Tuple& lhs, const Scalar& rhs) { };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_sdivassign_impl( Tuple& lhs, const Scalar& rhs) {
	    tuple_sdivassign_impl<typename boost::mpl::prior<Idx>::type,Tuple,Scalar>(lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(lhs) /= rhs;
	 };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_saddassign_impl( Tuple& lhs, const Scalar& rhs) { };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_saddassign_impl( Tuple& lhs, const Scalar& rhs) {
	    tuple_saddassign_impl<typename boost::mpl::prior<Idx>::type,Tuple,Scalar>(lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(lhs) += rhs;
	 };



      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_add_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_add_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) {
	    tuple_add_impl<typename boost::mpl::prior<Idx>::type,Tuple>(result,lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) + get<boost::mpl::prior<Idx>::type::value>(rhs);

	    //std::cout << " performing addition at level " << boost::mpl::prior<Idx>::type::value << std::endl; 

	    //if( boost::mpl::prior<Idx>::type::value == 2 )
	    //{
	       //std::cout << " adding lhs " << get<boost::mpl::prior<Idx>::type::value>(lhs)(2639) << std::endl; 
	       //std::cout << " with rhs " << get<boost::mpl::prior<Idx>::type::value>(rhs)(2639) << std::endl; 
	       //std::cout << " result " << get<boost::mpl::prior<Idx>::type::value>(rhs)(2639) << std::endl; 
	    //}
	 };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_sub_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_sub_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) {
	    tuple_sub_impl<typename boost::mpl::prior<Idx>::type,Tuple>(result,lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) - get<boost::mpl::prior<Idx>::type::value>(rhs);
	 };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_neg_impl( Tuple& result, const Tuple& lhs) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_neg_impl( Tuple& result, const Tuple& lhs) {
	    tuple_neg_impl<typename boost::mpl::prior<Idx>::type,Tuple>(result,lhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = -get<boost::mpl::prior<Idx>::type::value>(lhs);
	 };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_mul_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_mul_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) {
	    tuple_mul_impl<typename boost::mpl::prior<Idx>::type,Tuple>(result,lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) * get<boost::mpl::prior<Idx>::type::value>(rhs);
	 };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_div_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_div_impl( Tuple& result, const Tuple& lhs, const Tuple& rhs) {
	    tuple_div_impl<typename boost::mpl::prior<Idx>::type,Tuple>(result,lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) / get<boost::mpl::prior<Idx>::type::value>(rhs);
	 };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_muls_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) { };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_muls_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) {
	    tuple_muls_impl<typename boost::mpl::prior<Idx>::type,Tuple,Scalar>(result,lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) * rhs;
	 };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_smul_impl( Tuple& result, const Scalar& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_smul_impl( Tuple& result, const Scalar& lhs, const Tuple& rhs) {
	    tuple_smul_impl<typename boost::mpl::prior<Idx>::type,Tuple,Scalar>(result,lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = lhs * get<boost::mpl::prior<Idx>::type::value>(rhs);
	 };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_divs_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) { };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_divs_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) {
	    tuple_divs_impl<typename boost::mpl::prior<Idx>::type,Tuple,Scalar>(result,lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) / rhs;
	 };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_adds_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) { };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_adds_impl( Tuple& result, const Tuple& lhs, const Scalar& rhs) {
	    tuple_adds_impl<typename boost::mpl::prior<Idx>::type,Tuple,Scalar>(result,lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = get<boost::mpl::prior<Idx>::type::value>(lhs) + rhs;
	 };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_sadd_impl( Tuple& result, const Scalar& lhs, const Tuple& rhs) { };

      template <typename Idx, typename Tuple, typename Scalar>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_sadd_impl( Tuple& result, const Scalar& lhs, const Tuple& rhs) {
	    tuple_sadd_impl<typename boost::mpl::prior<Idx>::type,Tuple,Scalar>(result,lhs,rhs);
	    get<boost::mpl::prior<Idx>::type::value>(result) = lhs + get<boost::mpl::prior<Idx>::type::value>(rhs);
	 };


      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_abs_impl( Tuple& result, const Tuple& tpl ) { };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_abs_impl( Tuple& result, const Tuple& tpl ) {
	    tuple_abs_impl< typename boost::mpl::prior<Idx>::type,Tuple >(result,tpl);
	    get<boost::mpl::prior<Idx>::type::value>(result) = abs(get<boost::mpl::prior<Idx>::type::value>(tpl));
	 };


      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 double >::type tuple_norm_impl( const Tuple& tpl ) { return 0.0; };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 double >::type tuple_norm_impl( const Tuple& tpl ) {
	    return std::max( norm( get<boost::mpl::prior<Idx>::type::value>(tpl) ), tuple_norm_impl< typename boost::mpl::prior<Idx>::type,Tuple >(tpl) );
	 };


      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::equal_to< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_std_output_impl( std::ostream& lhs, const Tuple& rhs) { 
	    lhs << get<Idx::type::value>(rhs) << std::endl << std::endl;
	 };

      template <typename Idx, typename Tuple>
	 inline 
	 typename boost::enable_if< 
	 boost::mpl::greater< 
	 Idx, 
	 boost::mpl::size_t<0> 
	    >,
	 void >::type tuple_std_output_impl( std::ostream& lhs, const Tuple& rhs) {
	    tuple_std_output_impl<typename boost::mpl::prior<Idx>::type,Tuple>(lhs,rhs);
	    lhs << get<Idx::type::value>(rhs) << std::endl << std::endl;
	 };

   }; // detail

   // tuple + tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple >::type operator +(const Tuple& lhs, const Tuple& rhs) {
		  Tuple result;
		  detail::tuple_add_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs, rhs);
		  return result;
	       };

   // tuple - tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple >::type operator -(const Tuple& lhs, const Tuple& rhs) {
		  Tuple result;
		  detail::tuple_sub_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs, rhs);
		  return result;
	       };

   // -tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple >::type operator -(const Tuple& lhs) {
		  Tuple result;
		  detail::tuple_neg_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs);
		  return result;
	       };

   // tuple += tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple& >::type operator +=(Tuple& lhs, const Tuple& rhs) {
		  detail::tuple_addassign_impl<arithmetic_tuple_size<Tuple>,Tuple>(lhs, rhs);
		  return lhs;
	       };

   // tuple -= tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple& >::type operator -=(Tuple& lhs, const Tuple& rhs) {
		  detail::tuple_subassign_impl<arithmetic_tuple_size<Tuple>,Tuple>(lhs, rhs);
		  return lhs;
	       };

   // tuple * tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple >::type operator *(const Tuple& lhs, const Tuple& rhs) {
		  Tuple result;
		  detail::tuple_mul_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs, rhs);
		  return result;
	       };

   // tuple / tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple >::type operator /(const Tuple& lhs, const Tuple& rhs) {
		  Tuple result;
		  detail::tuple_div_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, lhs, rhs);
		  return result;
	       };

   // tuple * scalar
   template <typename Tuple, typename Scalar>
      typename boost::enable_if< boost::mpl::and_< is_instance_of_arithmetic_tuple< Tuple >, is_scalar< Scalar > >,
	       Tuple >::type operator *(const Tuple& lhs, const Scalar& rhs) {
		  Tuple result;
		  detail::tuple_muls_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(result, lhs, rhs);
		  return result;
	       };

   // scalar * tuple
   template <typename Tuple, typename Scalar>
      typename boost::enable_if< boost::mpl::and_< is_instance_of_arithmetic_tuple< Tuple >, is_scalar< Scalar > >,
	       Tuple >::type operator *(const Scalar& lhs, const Tuple& rhs) {
		  Tuple result;
		  detail::tuple_smul_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(result, lhs, rhs);
		  return result;
	       };

   // tuple * scalar
   template <typename Tuple, typename Scalar>
      typename boost::enable_if< boost::mpl::and_< is_instance_of_arithmetic_tuple< Tuple >, is_scalar< Scalar > >,
	       Tuple >::type operator /(const Tuple& lhs, const Scalar& rhs) {
		  Tuple result;
		  detail::tuple_divs_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(result, lhs, rhs);
		  return result;
	       };

   // tuple + scalar
   template <typename Tuple, typename Scalar>
      typename boost::enable_if< boost::mpl::and_< is_instance_of_arithmetic_tuple< Tuple >, is_scalar< Scalar > >,
	       Tuple >::type operator +(const Tuple& lhs, const Scalar& rhs) {
		  Tuple result;
		  detail::tuple_adds_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(result, lhs, rhs);
		  return result;
	       };

   // scalar + tuple
   template <typename Tuple, typename Scalar>
      typename boost::enable_if< boost::mpl::and_< is_instance_of_arithmetic_tuple< Tuple >, is_scalar< Scalar > >,
	       Tuple >::type operator +(const Scalar& lhs, const Tuple& rhs) {
		  Tuple result;
		  detail::tuple_sadd_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(result, lhs, rhs);
		  return result;
	       };

   // abs( tuple )
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple >::type abs(const Tuple& tpl) {
		  Tuple result;
		  detail::tuple_abs_impl<arithmetic_tuple_size<Tuple>,Tuple>(result, tpl);
		  return result;
	       };

   //norm( tuple )
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       double >::type norm(const Tuple& tpl) {
		  return detail::tuple_norm_impl<arithmetic_tuple_size<Tuple>,Tuple>(tpl);
	       };

   // tuple *= tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple& >::type operator *=(Tuple& lhs, const Tuple& rhs) {
		  detail::tuple_mulassign_impl<arithmetic_tuple_size<Tuple>,Tuple>(lhs, rhs);
		  return lhs;
	       };

   // tuple /= tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       Tuple& >::type operator /=(Tuple& lhs, const Tuple& rhs) {
		  detail::tuple_divassign_impl<arithmetic_tuple_size<Tuple>,Tuple>(lhs, rhs);
		  return lhs;
	       };

   // tuple *= scalar
   template <typename Tuple, typename Scalar>
      typename boost::enable_if< boost::mpl::and_< is_instance_of_arithmetic_tuple< Tuple >, is_scalar< Scalar > >,
	       Tuple& >::type operator *=(Tuple& lhs, const Scalar& rhs) {
		  detail::tuple_smulassign_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(lhs, rhs);
		  return lhs;
	       };

   // tuple /= scalar
   template <typename Tuple, typename Scalar>
      typename boost::enable_if< boost::mpl::and_< is_instance_of_arithmetic_tuple< Tuple >, is_scalar< Scalar > >,
	       Tuple& >::type operator /=(Tuple& lhs, const Scalar& rhs) {
		  detail::tuple_sdivassign_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(lhs, rhs);
		  return lhs;
	       };

   // tuple += scalar
   template <typename Tuple, typename Scalar>
      typename boost::enable_if< boost::mpl::and_< is_instance_of_arithmetic_tuple< Tuple >, is_scalar< Scalar > >,
	       Tuple& >::type operator +=(Tuple& lhs, const Scalar& rhs) {
		  detail::tuple_saddassign_impl<arithmetic_tuple_size<Tuple>,Tuple,Scalar>(lhs, rhs);
		  return lhs;
	       };

   // ostream << tuple
   template <typename Tuple>
      typename boost::enable_if< is_instance_of_arithmetic_tuple<Tuple>,
	       std::ostream& >::type operator <<(std::ostream& out, const Tuple& rhs) {
		  detail::tuple_std_output_impl<typename boost::mpl::prior< arithmetic_tuple_size<Tuple> >::type, Tuple>(out,rhs);
	       };


   template <int Idx, typename Tuple>
      struct arithmetic_tuple_element {
	 typedef typename std::tuple_element< Idx, Tuple >::type type;
      };

   /* Specialization, see general template docs. */
   template <int Idx, typename... T>
      struct arithmetic_tuple_element< Idx, arithmetic_tuple<T...> > {
	 typedef typename std::tuple_element< Idx, typename arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
      };

   /* Specialization, see general template docs. */
   template <int Idx, typename... T>
      struct arithmetic_tuple_element< Idx, const arithmetic_tuple<T...> > {
	 typedef typename std::tuple_element< Idx, const typename arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
      };

   /* Specialization, see general template docs. */
   template <int Idx, typename... T>
      struct arithmetic_tuple_element< Idx, volatile arithmetic_tuple<T...> > {
	 typedef typename std::tuple_element< Idx, volatile typename arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
      };

   /* Specialization, see general template docs. */
   template <int Idx, typename... T>
      struct arithmetic_tuple_element< Idx, const volatile arithmetic_tuple<T...> > {
	 typedef typename std::tuple_element< Idx, const volatile typename arithmetic_tuple<T...>::arithmetic_tuple_base_class >::type type;
      };


   namespace detail {

      template <typename Idx, typename T, typename Tuple>
	 struct arithmetic_tuple_index_of_impl; // forward-decl

      template <typename Idx, typename T, typename Tuple>
	 struct arithmetic_tuple_index_of_dispatch {
	    typedef boost::mpl::prior<Idx> prior_Idx;
	    typedef typename boost::mpl::if_<
	       boost::is_same< T, typename arithmetic_tuple_element<prior_Idx::type::value, Tuple>::type >,
	       prior_Idx,
	       arithmetic_tuple_index_of_impl< typename prior_Idx::type, T, Tuple > >::type::type type;
	 };

      template <typename T, typename Tuple>
	 struct arithmetic_tuple_index_of_failure { };

      template <typename Idx, typename T, typename Tuple>
	 struct arithmetic_tuple_index_of_impl {
	    typedef typename boost::mpl::if_<
	       boost::mpl::greater< Idx, boost::mpl::size_t<0> >,
	       arithmetic_tuple_index_of_dispatch< Idx, T, Tuple >,
	       arithmetic_tuple_index_of_failure< T, Tuple > >::type::type type;
	 };


      template <typename Idx, typename T, typename Tuple>
	 struct arithmetic_tuple_has_type_impl; // forward-decl

      template <typename Idx, typename T, typename Tuple>
	 struct arithmetic_tuple_has_type_dispatch {
	    typedef boost::mpl::prior<Idx> prior_Idx;
	    typedef typename boost::mpl::if_<
	       boost::is_same< T, typename arithmetic_tuple_element<prior_Idx::type::value, Tuple>::type >,
	       true_type,
	       arithmetic_tuple_has_type_impl< typename prior_Idx::type, T, Tuple > >::type::type type;
	 };

      template <typename Idx, typename T, typename Tuple>
	 struct arithmetic_tuple_has_type_impl {
	    typedef typename boost::mpl::if_<
	       boost::mpl::greater< Idx, boost::mpl::size_t<0> >,
	       arithmetic_tuple_has_type_dispatch< Idx, T, Tuple >,
	       false_type >::type::type type;
	 };



   }; // Detail


   template <typename T, typename Tuple>
      struct arithmetic_tuple_index_of {
	 typedef typename detail::arithmetic_tuple_index_of_impl< arithmetic_tuple_size<Tuple>, T, Tuple >::type type;
      };


   template <typename T, typename Tuple>
      struct arithmetic_tuple_has_type {
	 typedef typename detail::arithmetic_tuple_has_type_impl< arithmetic_tuple_size<Tuple>, T, Tuple >::type type;
	 BOOST_STATIC_CONSTANT(bool, value = type::value);
      };


}; // Reak
