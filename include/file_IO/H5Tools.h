
/*******************************************************************************************//** @file
 *  		
 * 	file: 		H5Tools.h
 * 	contents:  	Helper Functions to write containers to hdf5 file
 * 
 ****************************************************************************************************/


#pragma once

#include <string>
#include <frequencies/matsubara_space.h>
#include <grid.h>
#include <vector>
#include <boost/multi_array.hpp>
#include <file_IO/H5Compat.h>

typedef std::complex<double> dcomplex; 

// -- Convenient write functions to HDF5 group
void write_fermionic_matsubara_space( const fermionic_matsubara_space_t& fgrid, H5::Location& group, const std::string& dataset_name = std::string("fgrid") ); 	///<   	Write fermionic grid to given group of hdf5 file
void write_bosonic_matsubara_space( const bosonic_matsubara_space_t& bgrid, H5::Location& group, const std::string& dataset_name = std::string("bgrid") ); 	///<   	Write bosonic grid to given group of hdf5 file
//void write( const grid_t<2>& momgrid, H5::Location& group, const std::string& dataset_name = std::string("momgrid") ); 	///<   	Write bosonic momentum grid to given group of hdf5 file
void write( const std::vector<double >& points, unsigned dim, unsigned points_count, H5::Location& group, const std::string& dataset_name  = std::string("momgrid") );

// deprecated
void write( double scalar, H5::Location& group, const std::string& dataset_name ); 				///<   	Write a scalar value to given group of hdf5 file


void write_double( double scalar, H5::Location& group, const std::string& dataset_name ); 				///<   	Write a scalar value to given group of hdf5 file

void write( int integer, H5::Location& group, const std::string& dataset_name ); 					///<   	Write a integer value to given group of hdf5 file
void write( std::string text, H5::Location& group, const std::string& dataset_name ); 				///<   	Write a scalar value to given group of hdf5 file
void write( std::vector<double> vec_scalar, H5::Location& group, const std::string& dataset_name ); 		///<   	Write a scalar valued vector to given group of hdf5 file

void write( std::vector<std::string > vec_text, H5::Location& group, const std::string& dataset_name ); 		///<   	Write a scalar valued vector to given group of hdf5 file

void write_vec_int( std::vector<int> vec_scalar_int, H5::Location& group, const std::string& dataset_name ); 		///<   	Write a scalar valued vector to given group of hdf5 file

void write( std::vector< std::vector<int> > vec_scalar_indices, H5::Location& group, const std::string& dataset_name ); 		///<   	Write a scalar valued vector to given group of hdf5 file


void write( std::vector<std::complex<double> > vec_scalar_complex, H5::Location& group, const std::string& dataset_name ); 		///<   	Write a scalar valued vector to given group of hdf5 file

void write( std::complex<double> scalar_complex, H5::Location& group, const std::string& dataset_name ); 		///<   	Write a complex scalar to given group of hdf5 file


void write( std::vector<unsigned> vec_scalar, H5::Location& group, const std::string& dataset_name ); 		///<   	Write an unsigned valued vector to given group of hdf5 file

// -- Convenient write functions for Boost Multi-Array and Grids
template< size_t ndims > void write_real_valued_gf_object( const boost::multi_array<double, ndims>& my_arr, H5::Location& group, const std::string& dataset_name );  ///< Writes double valued boost multi_array to given group of hdf5 file
template< size_t ndims > void write_complex_valued_gf_object( const boost::multi_array<std::complex<double>, ndims>& my_arr, H5::Location& group, const std::string& dataset_name = std::string("") );  ///< Writes complex<double> valued boost multi_array to given group of hdf5 file

template< size_t ndims > void read( boost::multi_array<double, ndims>& my_arr, H5::Location& group, const std::string& dataset_name );  ///< Reads double valued boost multi_array to given group from hdf5 file
template< size_t ndims > void read( boost::multi_array<std::complex<double>, ndims>& my_arr, H5::Location& group, const std::string& dataset_name = std::string("") );  ///< Reads complex<double> valued boost multi_array from given group of hdf5 file

void read_double_attribute(double& attributeValue, H5::Location& group, const std::string& attributeName);
void read_int_attribute(int& attributeValue, H5::Location& group, const std::string& attributeName);

template< size_t ndims > void read_complex_alternative( boost::multi_array<dcomplex, ndims>& my_arr, H5::Location& group, const std::string& dataset_name );

#include <file_IO/H5Tools_impl.h> 	// contains implementations of template functions

