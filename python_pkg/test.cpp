#include <iostream>
#include <complex>
#include <vector>
#include <boost/multiprecision/cpp_complex.hpp>
#include <boost/multiprecision/cpp_int.hpp>
#include <boost/math/special_functions/trunc.hpp>
#include <boost/math/special_functions/round.hpp>

using namespace boost::multiprecision;

int main()
{
    std::cout << std::setprecision(20);
    // cpp_int x = pow(cpp_int(2), 1024);
    // std::cout << x << std::endl;
    
    cpp_int x = pow(cpp_int(2), 100);
    cpp_complex<100> y(3.231647865345, 4.47389678979);
    auto real = y.real();
    real = round(real*1e80) / 1e80;
    std::cout << std::pow(2, -1) << std::endl;
    return 0;
}