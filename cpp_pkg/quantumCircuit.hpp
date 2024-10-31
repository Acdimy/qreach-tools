#ifndef _QREACH_CIRCUIT
#define _QREACH_CIRCUIT

#include "../cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include <random>
#include <queue>
#include <vector>

namespace QREACH
{
    using namespace CFL_OBDD;
    

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

    class pureState {
        public:
        // root->level
        unsigned int numQubits = 0;
        CFLOBDD_COMPLEX_BIG state;
        pureState(){}
        pureState(CFLOBDD_COMPLEX_BIG s, int n) {
            state = s;
            numQubits = n;
        }
        // Is it necessary?
        pureState(pureState ps) {
            state = ps.state;
            numQubits = ps.numQubits;
        }
        void applyQuantumOpe(quantumOpe e) {}
        void applyGate(std::string name) {}

        BIG_COMPLEX_FLOAT applyProjector(std::vector<pureState> eigens) {
            assert(eigens.empty() == false);
            unsigned int level = ceil(log2(numQubits));
            CFLOBDD_COMPLEX_BIG res = VectorComplexFloatBoost::NoDistinctionNode(level+1, 0);
            auto tmpVec = Matrix1234ComplexFloatBoost::MatrixTranspose(state);
            tmpVec = Matrix1234ComplexFloatBoost::MatrixConjugate(tmpVec);
            BIG_COMPLEX_FLOAT total = 0;
            for(int i = 0; i < eigens.size(); i++) {
                auto tmp = Matrix1234ComplexFloatBoost::MatrixMultiplyV4(tmpVec, eigens[i].state);
                // Same as norm
                assert(tmp.root->rootConnection.returnMapHandle.Size() <= 2);
                auto resMap = tmp.root->rootConnection.returnMapHandle;
                BIG_COMPLEX_FLOAT amp;
                if(resMap.Size() == 2)
                    amp = (resMap[0] != 0) ? resMap[0] : resMap[1];
                else
                    amp = resMap[0];
                total += (amp * conj(amp));
                amp = conj(amp);
                // std::cout << "ApplyProjectorToState: " << amp << std::endl;
                res = res + (amp * eigens[i].state);
            }
            state = res;
            return total;
        }

        void applyInnerProduct(pureState s) {
            auto tmpVec = Matrix1234ComplexFloatBoost::MatrixTranspose(s.state);
            tmpVec = Matrix1234ComplexFloatBoost::MatrixConjugate(tmpVec);
            state = Matrix1234ComplexFloatBoost::MatrixMultiplyV4(tmpVec, state);
        }
        BIG_COMPLEX_FLOAT norm() {
            auto tmpVec = Matrix1234ComplexFloatBoost::MatrixTranspose(state);
            tmpVec = Matrix1234ComplexFloatBoost::MatrixConjugate(tmpVec);
            auto tmp = Matrix1234ComplexFloatBoost::MatrixMultiplyV4(tmpVec, state);
            // Same as applyProjector
            auto resMap = tmp.root->rootConnection.returnMapHandle;
            BIG_COMPLEX_FLOAT amp;
            if(resMap.Size() == 2)
                amp = (resMap[0] != 0) ? resMap[0] : resMap[1];
            else
                amp = resMap[0];
            // amp should be real
            // amp = conj(amp);
            return amp;
        }
        void applyDistance(pureState s) {}

    };

    class atomicProp {
        public:
        std::vector<pureState> space;
        int l;
        int r;
        atomicProp() {}
        bool equal(atomicProp p) {return false;}
        bool satisfy(pureState s) {
            return false;
        }
        bool satisfy(std::vector<pureState> rho) {return false;}
    };

    class CTransitionSys {

    };

    class QTransitionSys {
        public:
        std::vector<atomicProp> ap;
        std::vector<pureState> initState;
        void setAP();
        void setInitialState();
        void setTransition();
        void prepareCTL();
        void preparePCTL();

    };

    class quantumCircuit {
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