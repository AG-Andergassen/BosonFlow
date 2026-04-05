#include <delta_function_vector.h>

template <unsigned dim>
delta_function_vector_t<dim>::delta_function_vector_t(std::vector<coord_t<dim> > delta_peaks, std::vector<std::complex<double> > coeffs): m_delta_peaks(delta_peaks), m_coeffs(coeffs)
{
}


template <unsigned dim>
delta_function_vector_t<dim> delta_function_vector_t<dim>::operator+(const delta_function_vector_t<dim> &other)
{
    std::vector<std::complex<double> > new_coeffs;
    new_coeffs.insert(new_coeffs.end(), m_coeffs.begin(), m_coeffs.end());
    new_coeffs.insert(new_coeffs.end(), other.m_coeffs.begin(), other.m_coeffs.end());

    std::vector<coord_t<dim> > new_delta_peaks;
    new_delta_peaks.insert(new_delta_peaks.end(), m_delta_peaks.begin(), m_delta_peaks.end());
    new_delta_peaks.insert(new_delta_peaks.end(), other.m_delta_peaks.begin(), other.m_delta_peaks.end());
    
    return delta_function_vector_t<dim>(new_delta_peaks, new_coeffs);
}


template <unsigned dim>
void delta_function_vector_t<dim>::multiply_scalar(std::complex<double> number)
{
    for (unsigned i = 0; i < m_coeffs.size(); i ++)
	m_coeffs[i] *= number;
}


template <unsigned dim>
delta_function_vector_t<dim> delta_function_vector_t<dim>::group_fundamental_act(matrix_t<dim> g)
{
    std::vector<coord_t<dim> > new_delta_peaks;
    for (unsigned i = 0; i < m_delta_peaks.size(); i ++){
	vector_t<dim> v;
	v = g * m_delta_peaks[i];
	new_delta_peaks.push_back(v);
    }
    return delta_function_vector_t<dim>(new_delta_peaks, m_coeffs);
}


template <unsigned dim>
std::complex<double> delta_function_vector_t<dim>::evaluate(coord_t<dim> x)
{
    std::complex<double> result = 0;
    for (unsigned i = 0; i < m_delta_peaks.size(); i ++){
	double d = (x - m_delta_peaks[i]).dot(x - m_delta_peaks[i]);
	if (d < epsilon)
	    result += m_coeffs[i];
    }
    return result;
}


template <unsigned dim>
std::complex<double> delta_function_vector_t<dim>::evaluate_fourier_transform(coord_t<dim> p)
{
    std::complex<double> result = 0;
    for (unsigned j = 0; j < m_delta_peaks.size(); j ++){
	result += m_coeffs[j] * std::exp(std::complex<double>(0,-1) *(p.dot(m_delta_peaks[j])) );
    }
    return result;
}


template <unsigned dim>
std::string delta_function_vector_t<dim>::get_pretty_string()
{
    std::string out = "";
    for (unsigned i = 0; i < m_delta_peaks.size(); i ++){
	if (out != "")
	    out += "+ (";
	out += std::to_string(std::real(m_coeffs[i]));
	out += "+";
	out += std::to_string(std::imag(m_coeffs[i]));
	out += "j)*Delta(" + std::to_string(m_delta_peaks[i](0)) + "," + std::to_string(m_delta_peaks[i](1))+")";
    }
    return out;
}
