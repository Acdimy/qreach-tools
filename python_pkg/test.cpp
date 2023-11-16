#include <iostream>
#include <complex>
#include <vector>
// #include <boost/multiprecision/cpp_complex.hpp>

// namespace mp = boost::multiprecision;

// class QuantumGate {
//     /* data */
//     public:
//     std::string name;
//     std::vector<unsigned int> index;
//     std::vector<double> vars;
//     public:
//     QuantumGate(std::string nam, std::vector<unsigned int> idx, std::vector<double> pars) {
//         this->name = nam;
//         this->index = idx;
//         this->vars = pars;
//     }
//     ~QuantumGate() {
//         // necessary?
//         this->index.clear();
//         this->vars.clear();
//     }
// };

class QuantumOperator
{
private:
    int numChannel;
    std::vector<int> operators;
public:
    QuantumOperator() {
        numChannel = 0;
    }
    // void appendChannel() {
    //     numChannel++;
    //     std::vector<QuantumGate> v;
    //     operator.push_back(v);
    // }
    // void appendGate(std::string name, std::vector<unsigned int> idx, std::vector<double> par) {
    //     operator[numChannel-1].push_back(QuantumGate(name, idx, par));
    // }
    ~QuantumOperator(){};
};

int main()
{
    // std::cout << "First, some operations we usually perform with std::complex:\n";
    // complex_number_examples<std::complex<double>>();
    std::cout << "\nNow the same operations performed using quad precision complex numbers:\n";
    // complex_number_examples<boost::multiprecision::cpp_complex_100>();
    QuantumOperator c;

    return 0;
}