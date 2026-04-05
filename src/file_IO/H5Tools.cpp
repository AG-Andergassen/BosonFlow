#include <file_IO/H5Tools.h>
#include <iostream>

using namespace H5; 
using std::string; 
using std::vector; 


void write_fermionic_matsubara_space( const fermionic_matsubara_space_t& fgrid, Location& group, const string& dataset_name )
{
   
   hsize_t grid_dim[1] = { 2* fgrid.get_pos_freq_count() };
   DataSpace grid_dataspace( 1, grid_dim );

   DataSet grid_dataset = group.createDataSet( dataset_name, PredType::NATIVE_DOUBLE, grid_dataspace );

   grid_dataset.write( fgrid.data(), PredType::NATIVE_DOUBLE );

   grid_dataset.close();
   grid_dataspace.close();
}

void write_bosonic_matsubara_space( const bosonic_matsubara_space_t& bgrid, Location& group, const string& dataset_name )
{
   
   hsize_t grid_dim[1] = { 2* bgrid.get_pos_freq_count() + 1 };
   DataSpace grid_dataspace( 1, grid_dim );

   DataSet grid_dataset = group.createDataSet( dataset_name, PredType::NATIVE_DOUBLE, grid_dataspace );

   grid_dataset.write( bgrid.data(), PredType::NATIVE_DOUBLE );

   grid_dataset.close();
   grid_dataspace.close();

}

/*
void write( const grid_t<2>& momgrid, Group& group, const string& dataset_name )
{
   
   hsize_t grid_dim[2] = { momgrid.m_points.size(), 2 };
   hsize_t hndims(2); 		
   DataSpace grid_dataspace( hndims, grid_dim );

   DataSet grid_dataset = group.createDataSet( dataset_name, PredType::NATIVE_DOUBLE, grid_dataspace );
 
   grid_dataset.write( momgrid.m_points.data(), PredType::NATIVE_DOUBLE );

}
*/


void write( const std::vector<double>& points, unsigned dim, unsigned points_count, Location& group, const string& dataset_name )
{
    std::cout << dataset_name << std::endl;
   
    hsize_t grid_dim[2] = { points_count, dim };
    hsize_t hndims(2); 		
    DataSpace grid_dataspace( hndims, grid_dim );
    
    DataSet grid_dataset = group.createDataSet( dataset_name, PredType::NATIVE_DOUBLE, grid_dataspace );
    
    grid_dataset.write( points.data(), PredType::NATIVE_DOUBLE );

    grid_dataset.close();
    grid_dataspace.close();
}


void write( double scalar, Location& group, const string& attr_name )  //deprecated
{

    DataSpace attr_dsp = DataSpace ( H5S_SCALAR );
    Attribute attr = group.createAttribute( attr_name, PredType::NATIVE_DOUBLE, attr_dsp );
    attr.write( PredType::NATIVE_DOUBLE, &scalar );

    attr.close();
    attr_dsp.close();
}

void write_double( double scalar, Location& group, const string& attr_name )
{
    DataSpace attr_dsp = DataSpace ( H5S_SCALAR );
    Attribute attr = group.createAttribute( attr_name, PredType::NATIVE_DOUBLE, attr_dsp );
    attr.write( PredType::NATIVE_DOUBLE, &scalar );

    attr.close();
    attr_dsp.close();
}


void write( int integer, Location& group, const string& attr_name )
{

    DataSpace attr_dsp = DataSpace ( H5S_SCALAR );
    Attribute attr = group.createAttribute( attr_name, PredType::NATIVE_INT, attr_dsp );
    attr.write( PredType::NATIVE_INT, &integer );
   
    attr.close();
    attr_dsp.close();
}

void write( string text, H5::Location& group, const string& dataset_name )
{
   
    DataSpace attr_dsp = DataSpace ( H5S_SCALAR );
    StrType strdatatype( PredType::C_S1, H5T_VARIABLE ); // String type with variable length
    Attribute attr = group.createAttribute( dataset_name, strdatatype, attr_dsp );
    const char *text_cstr = text.c_str();
    attr.write( strdatatype, &text_cstr );

    attr.close();
    attr_dsp.close();
}

void write( vector<double> vec_scalar, H5::Location& group, const string& dataset_name )
{

    hsize_t dims[1] = { vec_scalar.size() }; 
    DataSpace dataspace( 1, dims );
    DataSet dataset = group.createDataSet( dataset_name, PredType::NATIVE_DOUBLE, dataspace );
    dataset.write( vec_scalar.data(), PredType::NATIVE_DOUBLE );
   
    dataset.close();
    dataspace.close();
}

void write( vector<unsigned> vec_scalar, H5::Location& group, const string& dataset_name )
{

    hsize_t dims[1] = { vec_scalar.size() }; 
    DataSpace dataspace( 1, dims );
    DataSet dataset = group.createDataSet( dataset_name, PredType::NATIVE_UINT, dataspace );
    dataset.write( vec_scalar.data(), PredType::NATIVE_UINT );
   
    dataset.close();
    dataspace.close();
}

void write_vec_int( vector<int> vec_scalar, H5::Location& group, const string& dataset_name )
{

    hsize_t dims[1] = { vec_scalar.size() }; 
    DataSpace dataspace( 1, dims );
    DataSet dataset = group.createDataSet( dataset_name, PredType::NATIVE_INT, dataspace );
    dataset.write( vec_scalar.data(), PredType::NATIVE_INT );
   
    dataset.close();
    dataspace.close();
}


void write( vector< vector<int> > vec_indices, H5::Location& group, const string& dataset_name )
{
    unsigned idx_dim = 1;
    if (vec_indices.size() > 0){
	idx_dim = vec_indices[0].size();
    }

    std::vector<int> points;
    for (auto idx : vec_indices){
	for (auto i : idx)
	    points.push_back(i);
    }

    hsize_t grid_dim[2] = { vec_indices.size()*idx_dim, idx_dim };
    hsize_t hndims(2); 		
    DataSpace grid_dataspace( hndims, grid_dim );
    
    DataSet grid_dataset = group.createDataSet( dataset_name, PredType::NATIVE_INT, grid_dataspace );
    

    grid_dataset.write( points.data(), PredType::NATIVE_INT );
    
    grid_dataset.close();
    grid_dataspace.close();
}


void write( vector<string> vec_text, H5::Location& group, const string& dataset_name )
{
    vector<const char *> vec_text_cstr;
    for (int ii = 0 ; ii < vec_text.size(); ii ++){
	vec_text_cstr.push_back(vec_text[ii].c_str());
    }

    StrType strdatatype( PredType::C_S1, H5T_VARIABLE ); // String type with variable length

    hsize_t dims[1] = { vec_text.size() }; 
    DataSpace dataspace( 1, dims );
    DataSet dataset = group.createDataSet( dataset_name, strdatatype, dataspace );

    dataset.write( vec_text_cstr.data(), strdatatype );

    dataset.close();
    dataspace.close();
}

void write( std::vector<std::complex<double> > vec_scalar, H5::Location& group, const string& dataset_name )
{
    std::cout << dataset_name << std::endl;

    std::vector<double> vec_scalar_real, vec_scalar_imaginary;

    vec_scalar_real.resize(vec_scalar.size());
    vec_scalar_imaginary.resize(vec_scalar.size());

    std::transform( vec_scalar.data(), vec_scalar.data() + vec_scalar.size(), vec_scalar_real.data(), [](std::complex<double> a){ return std::real(a); } ); 
    std::transform( vec_scalar.data(), vec_scalar.data() + vec_scalar.size(), vec_scalar_imaginary.data(), [](std::complex<double> a){ return std::imag(a); } ); 
   

    hsize_t dims[1] = { vec_scalar.size() }; 
    DataSpace dataspace( 1, dims );
    DataSet dataset_real = group.createDataSet( string("RE") + dataset_name, PredType::NATIVE_DOUBLE, dataspace );
    DataSet dataset_imaginary = group.createDataSet( string("IM") + dataset_name, PredType::NATIVE_DOUBLE, dataspace );
    dataset_real.write( vec_scalar_real.data(), PredType::NATIVE_DOUBLE );
    dataset_imaginary.write( vec_scalar_imaginary.data(), PredType::NATIVE_DOUBLE );  

    dataset_real.close();
    dataset_imaginary.close();
    dataspace.close();
}


void write( std::complex<double> scalar, H5::Location& group, const string& attr_name )
{
    double real_part = std::real(scalar);
    double imag_part = std::imag(scalar);

    DataSpace attr_dsp_real = DataSpace ( H5S_SCALAR );
    Attribute attr_real = group.createAttribute( string("RE") + attr_name, PredType::NATIVE_DOUBLE, attr_dsp_real );
   
    DataSpace attr_dsp_imag = DataSpace ( H5S_SCALAR );
    Attribute attr_imag = group.createAttribute( string("IM") + attr_name, PredType::NATIVE_DOUBLE, attr_dsp_imag );
   
    attr_real.write( PredType::NATIVE_DOUBLE, &real_part );
    attr_imag.write( PredType::NATIVE_DOUBLE, &imag_part );

    attr_real.close();
    attr_imag.close();
    attr_dsp_real.close();
    attr_dsp_imag.close();
}



void read_double_attribute(double& attributeValue, H5::Location& group, const std::string& attributeName)
{
    try {
	H5::Attribute attribute = group.openAttribute(attributeName);
	
	H5::DataType attributeType = attribute.getDataType();
        if (attributeType.getSize() == sizeof(double)) {
            attribute.read(attributeType, &attributeValue);    
	} else {
	    std::cerr << "Attribute is not of type double." << std::endl;
	}
        attribute.close();
    } catch (H5::Exception& e) {
        // Handle exceptions here
        e.printErrorStack();
    }
}

void read_int_attribute(int& attributeValue, H5::Location& group, const std::string& attributeName)
{
    try {
	H5::Attribute attribute = group.openAttribute(attributeName);
	
	H5::DataType attributeType = attribute.getDataType();
        if (attributeType.getSize() == sizeof(int)) {
            attribute.read(attributeType, &attributeValue);
	} else {
	    std::cerr << "Attribute is not of type double." << std::endl;
	}
	attribute.close();
    } catch (H5::Exception& e) {
        // Handle exceptions here
        e.printErrorStack();
    }
}
