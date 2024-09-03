#ifndef _QREACH_CIRCUIT
#define _QREACH_CIRCUIT

#include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include <random>
#include <queue>
#include <vector>

namespace QREACH
{
    class quantumOpe {
        /* data */
        public:
        std::string name;
        std::vector<unsigned int> index;
        std::vector<double> vars;
        public:
        quantumOpe(std::string nam, std::vector<unsigned int> idx, std::vector<double> pars) {
            this->name = nam;
            this->index = idx;
            this->vars = pars;
        }
        ~quantumOpe() {
            // necessary?
            this->index.clear();
            this->vars.clear();
        }
    };
    class quantumCircuit
    {
    public:
        quantumCircuit(unsigned int numQubits, int seed);
        quantumCircuit();
        virtual ~quantumCircuit();
        virtual void setNumQubits(unsigned int numQubits) = 0;
        virtual void addGate(quantumOpe gate) = 0;
    protected:
        unsigned int numQubits;
        unsigned int realQubits;
        std::mt19937 mt;
        std::vector<quantumOpe> gateList;
    };
    
    quantumCircuit::quantumCircuit(unsigned int numQubits, int seed) :  numQubits(numQubits) 
    {
        mt.seed(seed);
        srand(seed);
        realQubits = 0;
    }
    quantumCircuit::quantumCircuit() : numQubits(0), realQubits(0) {}
    
    quantumCircuit::~quantumCircuit() {}
    
}

using namespace QREACH;
class CFLOBDDQuantumCircuit : public quantumCircuit {
    public:
    CFLOBDDQuantumCircuit(unsigned int numQubits,  int seed);
    CFLOBDDQuantumCircuit();
    ~CFLOBDDQuantumCircuit();
    void setNumQubits(unsigned int numQubits);
    void addGate(quantumOpe gate);
};

CFLOBDDQuantumCircuit::CFLOBDDQuantumCircuit(unsigned int numQubits, int seed) : QuantumCircuit(numQubits, seed)
{
    // Initialize
        CFLOBDDNodeHandle::InitNoDistinctionTable();
        CFLOBDDNodeHandle::InitAdditionInterleavedTable();
        CFLOBDDNodeHandle::InitReduceCache();
        InitPairProductCache();
        InitTripleProductCache();
        Matrix1234ComplexFloatBoost::Matrix1234Initializer();
        VectorComplexFloatBoost::VectorInitializer();
    //
        unsigned int level = ceil(log2(numQubits)) + 1;
}

CFLOBDDQuantumCircuit::CFLOBDDQuantumCircuit()
{
    // Initialize
        CFLOBDDNodeHandle::InitNoDistinctionTable();
        CFLOBDDNodeHandle::InitAdditionInterleavedTable();
        CFLOBDDNodeHandle::InitReduceCache();
        InitPairProductCache();
        InitTripleProductCache();
        Matrix1234ComplexFloatBoost::Matrix1234Initializer();
        VectorComplexFloatBoost::VectorInitializer();
        numQubits = 0;
    //
}

CFLOBDDQuantumCircuit::~CFLOBDDQuantumCircuit()
{
    DisposeOfTripleProductCache();
	DisposeOfPairProductCache();
	CFLOBDDNodeHandle::DisposeOfReduceCache();
}

#endif