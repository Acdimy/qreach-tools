#ifndef _QUANTUM_CIRCUIT
#define _QUANTUM_CIRCUIT

#include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include <random>
#include <queue>
#include <vector>


class QuantumGate {
    /* data */
    public:
    std::string name;
    std::vector<unsigned int> index;
    std::vector<double> vars;
    public:
    QuantumGate(std::string nam, std::vector<unsigned int> idx, std::vector<double> pars) {
        this->name = nam;
        this->index = idx;
        this->vars = pars;
    }
    // QuantumGate(std::string name, unsigned int i) {
    //     this->name = name;
    //     this->index.push_back(i);
    // }
    // QuantumGate(std::string name, int idx1, int idx2) {
    //     this->name = name;
    //     this->index.push_back((unsigned int) idx1);
    //     this->index.push_back((unsigned int) idx2);
    // }
    // QuantumGate(std::string name, int idx1, int idx2, int idx3) {
    //     this->name = name;
    //     this->index.push_back((unsigned int) idx1);
    //     this->index.push_back((unsigned int) idx2);
    //     this->index.push_back((unsigned int) idx3);
    // }
    // QuantumGate(std::string name, double f) {
    //     // Global phase
    //     this->name = name;
    //     this->vars.push_back(f);
    // }
    // QuantumGate(std::string name, int idx, double f) {
    //     // Global phase
    //     this->name = name;
    //     this->index.push_back((unsigned int) idx);
    //     this->vars.push_back(f);
    // }
    // QuantumGate(std::string name, int idx1, int idx2, double f) {
    //     // Global phase
    //     this->name = name;
    //     this->index.push_back((unsigned int) idx1);
    //     this->index.push_back((unsigned int) idx2);
    //     this->vars.push_back(f);
    // }
    // QuantumGate(std::string name, int idx, double f1, double f2, double f3) {
    //     // Global phase
    //     this->name = name;
    //     this->index.push_back((unsigned int) idx);
    //     this->vars.push_back(f1);
    //     this->vars.push_back(f2);
    //     this->vars.push_back(f3);
    // }
    ~QuantumGate() {
        // necessary?
        this->index.clear();
        this->vars.clear();
    }
};

class QuantumOperator
{
public:
    int numChannel = 0;
    std::vector<std::vector<QuantumGate>> operations;
public:
    QuantumOperator() {
        numChannel = 0;
    }
    void appendChannel() {
        numChannel++;
        std::vector<QuantumGate> v;
        operations.push_back(v);
    }
    void appendGate(std::string name, std::vector<unsigned int> idx, std::vector<double> par) {
        operations[numChannel-1].push_back(QuantumGate(name, idx, par));
    }
    ~QuantumOperator(){};
};



class QuantumCircuit {
    public:
        // Constructor
        QuantumCircuit(unsigned int numQubits, int seed);
        // Constructor
        QuantumCircuit();
        // Destructor
        virtual ~QuantumCircuit();
        // set qubit count;
        virtual void setNumQubits(unsigned int numQubits) = 0;
        // I = [[1 0] [0 1]]
        // For no-op
        virtual void ApplyIdentityGate(unsigned int index) = 0;
        // H = [[1 1] [1 -1]]
        // Also called Walsh Gate
        virtual void ApplyHadamardGate(unsigned int index) = 0;
        // X = [[0 1] [1 0]]
        // Also called Pauli-X or bit-flip
        virtual void ApplyNOTGate(unsigned int index) = 0;
        // Y = [[0 -i] [i 0]]
        virtual void ApplyPauliYGate(unsigned int index) = 0;
        // Z = [[1 0] [0 -1]]
        // Also called phase-flip
        virtual void ApplyPauliZGate(unsigned int index) = 0;
        // S = [[1 0] [0 i]]
        // Also called √(Z) Gate
        virtual void ApplySGate(unsigned int index) = 0;
        // CNOT = [[1 0 0 0] [0 1 0 0] [0 0 0 1] [0 0 1 0]]
        virtual void ApplyCNOTGate(long int controller, long int controlled) = 0;
        // Ph = e^{i * π * phase} I
        virtual void ApplyGlobalPhase(double phase) = 0;
        // SWAP = [[1 0 0 0] [0 0 1 0] [0 1 0 0] [0 0 0 1]]
        virtual void ApplySwapGate(long int index1, long int index2) = 0;
        // iSWAP = [[1 0 0 0] [0 0 i 0] [0 i 0 0] [0 0 0 1]]
        virtual void ApplyiSwapGate(long int index1, long int index2) = 0;
        // CZ = [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 -1]]
        // Also called Controlled Phase Flip Gate
        virtual void ApplyCZGate(long int controller, long int controlled) = 0;
        // CPhase = [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 e^(i * π * θ)]]
        virtual void ApplyCPGate(long int controller, long int controlled, double theta) = 0;
        // P = [[1 0] [0 e^{i * π * θ}]]
        virtual void ApplyPhaseShiftGate(unsigned int index, double theta) = 0;
        // T = P(1/4) = [[1 0] [0 e^{i * π * 1/4}]]
        virtual void ApplyTGate(unsigned int index) = 0;
        // CS = [[1 0 0 0] [0 1 0 0] [0 0 1 0] [0 0 0 e^(i * π * 1/2)]]
        virtual void ApplyCSGate(long int controller, long int controlled) = 0;
        // CCNOT or Toffoli gate
        virtual void ApplyCCNOTGate(long int controller1, long int controller2, long int controlled) = 0;
        // CSWAP or Fredkin gate
        virtual void ApplyCSwapGate(long int controller, long int index1, long int index2) = 0; 
        // Obtain Probability
        virtual long double GetProbability(std::map<unsigned int, int>& qubit_vals) = 0;
        // Measure
        virtual std::string Measure() = 0;
        // Get Path Counts
        virtual unsigned long long int GetPathCount(long double prob) = 0;
    
    protected:
        unsigned int numQubits;
        unsigned int realQubits;
        unsigned int hadamard_count;
        std::mt19937 mt;
};

using namespace CFL_OBDD;

CFLOBDD_COMPLEX_BIG ApplyGateF(unsigned int n, unsigned int i, CFLOBDD_COMPLEX_BIG(*f)(unsigned int));
CFLOBDD_COMPLEX_BIG ApplyGateFWithParam(unsigned int n, unsigned int i, CFLOBDD_COMPLEX_BIG(*f)(unsigned int, double), double theta);
CFLOBDD_COMPLEX_BIG ApplyGateFWithParamVec(unsigned int n, unsigned int i, CFLOBDD_COMPLEX_BIG(*f)(unsigned int, std::vector<double>), std::vector<double> v);


class CFLOBDDQuantumCircuit : public QuantumCircuit {
    public:
        CFLOBDDQuantumCircuit(unsigned int numQubits,  int seed);
        CFLOBDDQuantumCircuit();
        ~CFLOBDDQuantumCircuit();
        void setNumQubits(unsigned int numQubits);
        void ApplyIdentityGate(unsigned int index);
        void ApplyHadamardGate(unsigned int index);
        void ApplyNOTGate(unsigned int index);
        void ApplyPauliYGate(unsigned int index);
        void ApplyPauliZGate(unsigned int index);
        void ApplySGate(unsigned int index);
        void ApplyCNOTGate(long int controller, long int controlled);
        void ApplyGlobalPhase(double phase);
        void ApplySwapGate(long int index1, long int index2);
        void ApplyiSwapGate(long int index1, long int index2);
        void ApplyCZGate(long int controller, long int controlled);
        void ApplyCPGate(long int controller, long int controlled, double theta);
        void ApplyPhaseShiftGate(unsigned int index, double theta);
        void ApplyTGate(unsigned int index);
        void ApplyCSGate(long int controller, long int controlled);
        void ApplyCCNOTGate(long int controller1, long int controller2, long int controlled);
        void ApplyCSwapGate(long int controller, long int index1, long int index2); 
        void ApplyU3Gate(unsigned int index, double theta, double phi, double lambda);
        void ApplyArbitraryGate(unsigned int index, std::vector<double> v);
        long double GetProbability(std::map<unsigned int, int>& qubit_vals);
        std::string Measure();
        unsigned long long int GetPathCount(long double prob);


        // Equivalence checking used
        void setBasicStateVector(std::string s);
        void setRealQubits(unsigned int q);
        void setInitGate();
        void setMaximalEntangledVec();
        void setProjector(std::vector<double> real, std::vector<double> imag, std::vector<int> index);
        void setProjectorFromS();
        void ApplyProjectorToEntangle();
        void ApplyProjectorToState();
        void pushStateToCache();
        void swapStateAndCache(unsigned int index = 0);
        void addCacheToState();

        void appendGateSeries(std::string name, std::vector<unsigned int> index, std::vector<double> vars, bool newChannel=false);
        void clearGateSeries();
        int ApplyGateSeries(int channelIdx);
        unsigned int reachability();

        void printStateColMajor();
        void printStateColHead();
        void printYield();
        void print();
        void printRV(std::string type = "state");
        unsigned int printSize(std::string type = "state");
        void printProjector();
        int getRealQubits();
    private:
        CFLOBDD_COMPLEX_BIG stateVector;
        std::vector<CFLOBDD_COMPLEX_BIG> stateVectorCache;
        std::vector<CFLOBDD_COMPLEX_BIG> stateProjector;
        std::deque<unsigned int> state_queue;
        QuantumOperator circuitGates;
        // std::vector<CFLOBDD_COMPLEX_BIG> kraus_operators;
        //unsigned int numQubits;
};


#include "cflobdd/CFLOBDD/wmatrix1234_complex_fb_mul.h"

class WeightedBDDQuantumCircuit : public QuantumCircuit
{
    public:
        // Constructor
        WeightedBDDQuantumCircuit(unsigned int numQubits, int seed);
        // Constructor
        WeightedBDDQuantumCircuit();
        // Destructor
        ~WeightedBDDQuantumCircuit();
        // set qubit count;
        void setNumQubits(unsigned int numQubits);
        void ApplyIdentityGate(unsigned int index);
        void ApplyHadamardGate(unsigned int index);
        void ApplyNOTGate(unsigned int index);
        void ApplyPauliYGate(unsigned int index);
        void ApplyPauliZGate(unsigned int index);
        void ApplySGate(unsigned int index);
        void ApplyCNOTGate(long int controller, long int controlled);
        void ApplyGlobalPhase(double phase);
        void ApplySwapGate(long int index1, long int index2);
        void ApplyiSwapGate(long int index1, long int index2);
        void ApplyCZGate(long int controller, long int controlled);
        void ApplyCPGate(long int controller, long int controlled, double theta);
        void ApplyPhaseShiftGate(unsigned int index, double theta);
        void ApplyTGate(unsigned int index);
        void ApplyCSGate(long int controller, long int controlled);
        void ApplyCCNOTGate(long int controller1, long int controller2, long int controlled);
        void ApplyCSwapGate(long int controller, long int index1, long int index2); 
        long double GetProbability(std::map<unsigned int, int>& qubit_vals);
        std::string Measure();
        unsigned long long int GetPathCount(long double prob);
    private:
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL stateVector;
};

#endif
