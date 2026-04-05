#include <frequencies/matsubara_space.h>
#include <fstream>
#include <string>
#include <tuple>

matsubara_space_t::matsubara_space_t( unsigned int pos_freq_count_, double step_size_, Statistics zeta )
    :pos_freq_count( pos_freq_count_ ), step_size( step_size_ ){
#ifndef STATIC_CALCULATION 
    for( int i = -pos_freq_count; i < pos_freq_count + (1 ? zeta == Statistics::BOSE : 0); ++i )
	grid_points.push_back( step_size * i + static_cast<unsigned>(zeta)*step_size / 2.0 );
#else
    grid_points.push_back(0);
#endif
}


void matsubara_space_t::print( const std::string fname ){
    std::fstream data;
    data.open("dat/" + fname, std::ios::out | std::ios::trunc );
    
    for( double freq : grid_points )
        data << freq << std::endl;
    
    data.flush();
    data.close();
}

unsigned int matsubara_space_t::get_pos_freq_count() const{
    return pos_freq_count;
}


double matsubara_space_t::get_step_size() const{
    return step_size;
}

