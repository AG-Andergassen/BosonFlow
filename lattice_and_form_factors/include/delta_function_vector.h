#pragma once 

#include <tuple>
#include <vector>
#include <string>
#include <complex>
#include <Eigen/Core>


template <unsigned dim>
class delta_function_vector_t
{
 public:
    delta_function_vector_t(std::vector<coord_t<dim> > delta_peaks_, std::vector<std::complex<double> > coeffs_);

    delta_function_vector_t<dim> operator+(const delta_function_vector_t<dim> &other);

    void multiply_scalar(std::complex<double> number);

    delta_function_vector_t<dim> group_fundamental_act(matrix_t<dim> g);

    std::complex<double> evaluate(coord_t<dim> x);

    std::complex<double> evaluate_fourier_transform(coord_t<dim> p);

    std::string get_pretty_string();

    std::vector<coord_t<dim> > m_delta_peaks;
    std::vector<std::complex<double> > m_coeffs;
    
    double epsilon = 1e-10;
};


#include "../src/delta_function_vector.tpp"
