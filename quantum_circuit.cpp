#include "quantum_circuit.h"
#include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include "cflobdd/CFLOBDD/vector_complex_float_boost.h"
#include <random>

QuantumCircuit::QuantumCircuit(unsigned int numQubits, int seed) :  numQubits (numQubits) 
{
    mt.seed(seed);
    srand(seed);
    hadamard_count = 0;
    realQubits = 0;
}
QuantumCircuit::QuantumCircuit() :  numQubits (0), hadamard_count (0), realQubits (0) {}
QuantumCircuit::~QuantumCircuit() {}

// using namespace CFL_OBDD;

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
    stateVector = VectorComplexFloatBoost::MkBasisVector(level, 0);
    // auto H = ApplyGateF(std::pow(2, level-1), 0, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
    // kraus_operators.push_back(H);
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

void CFLOBDDQuantumCircuit::setNumQubits(unsigned int num)
{
    numQubits = num;
    unsigned int level = ceil(log2(numQubits)) + 1;
    stateVector = VectorComplexFloatBoost::MkBasisVector(level, 0); 
}

CFLOBDD_COMPLEX_BIG ApplyGateF(unsigned int n, unsigned int i, CFLOBDD_COMPLEX_BIG(*f)(unsigned int))
{
    // i is the index of the applied qubit
    if (n == 1)
    {
        // but here i is the level
        return f(1);
    }
    else {
        int level = ceil(log2(n/2));
        if (i < n/2)
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            CFLOBDD_COMPLEX_BIG H = ApplyGateF(n/2, i, f);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(H, T);
        }
        else
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(T, ApplyGateF(n/2, i - n/2, f)); 
        }
    }
}

CFLOBDD_COMPLEX_BIG ApplyGateFWithParam(unsigned int n, unsigned int i, CFLOBDD_COMPLEX_BIG(*f)(unsigned int, double), double theta)
{
    if (n == 1)
    {
        return f(1, theta);
    }
    else {
        int level = ceil(log2(n/2));
        if (i < n/2)
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            CFLOBDD_COMPLEX_BIG H = ApplyGateFWithParam(n/2, i, f, theta);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(H, T);
        }
        else
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(T, ApplyGateFWithParam(n/2, i - n/2, f, theta)); 
        }
    }
}

CFLOBDD_COMPLEX_BIG ApplyGateFWithParamVec(unsigned int n, unsigned int i, CFLOBDD_COMPLEX_BIG(*f)(unsigned int, std::vector<double>), std::vector<double> v)
{
    if (n == 1)
    {
        return f(1, v);
    }
    else {
        int level = ceil(log2(n/2));
        if (i < n/2)
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            CFLOBDD_COMPLEX_BIG H = ApplyGateFWithParamVec(n/2, i, f, v);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(H, T);
        }
        else
        {
            CFLOBDD_COMPLEX_BIG T = Matrix1234ComplexFloatBoost::MkIdRelationInterleaved(level + 1);
            return Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(T, ApplyGateFWithParamVec(n/2, i - n/2, f, v)); 
        }
    }
}

bool checkForInit(unsigned int numQubits)
{
    return numQubits != 0;
}

void CFLOBDDQuantumCircuit::ApplyIdentityGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    // std::cout << stateVector.root->level << std::endl;
    auto H = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
    // std::cout << H.root->level << std::endl;
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyHadamardGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto H = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkWalshInterleaved);
    // kraus_operators[0] = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(kraus_operators[0], H);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, stateVector);
    hadamard_count++;
}

void CFLOBDDQuantumCircuit::ApplyNOTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto X = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkNegationMatrixInterleaved);
    // kraus_operators[0] = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(kraus_operators[0], X);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(X, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyPauliYGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto Y = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkPauliYMatrixInterleaved);
    int size = Y.root->rootConnection.returnMapHandle.Size();
    // std::cout << "check Y entries" << std::endl;
    // for(int i=0; i<size; i++) {
    //     std::cout << "(" << Y.root->rootConnection.returnMapHandle[i].real() << ", " 
    //     << Y.root->rootConnection.returnMapHandle[i].imag() << ") ";
    // } std::cout << std::endl;
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(Y, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyPauliZGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto Z = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkPauliZMatrixInterleaved);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(Z, stateVector);
}

void CFLOBDDQuantumCircuit::ApplySGate(unsigned int index)
{
    auto S = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkSGateInterleaved);
    // kraus_operators[0] = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(kraus_operators[0], S);
    // int size = S.root->rootConnection.returnMapHandle.Size();
    // std::cout << "check S entries" << std::endl;
    // for(int i=0; i<size; i++) {
    //     std::cout << "(" << S.root->rootConnection.returnMapHandle[i].real() << ", " 
    //     << S.root->rootConnection.returnMapHandle[i].imag() << ") ";
    // } std::cout << std::endl;
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyU3Gate(unsigned int index, double theta, double phi, double lambda)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    std::vector<double> v;
    v.push_back(theta); v.push_back(phi); v.push_back(lambda);
    auto U = ApplyGateFWithParamVec(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkU3GateInterleaved, v);
    // std::cout << U << std::endl;
    // std::cout << U.root->rootConnection.returnMapHandle << std::endl;
    // kraus_operators[0] = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(kraus_operators[0], U);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(U, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyArbitraryGate(unsigned int index, std::vector<double> v)
{
    //double v1_r, double v1_i, double v2_r, double v2_i, double v3_r, double v3_i, double v4_r, double v4_i
    assert(v.size() == 8);
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    // std::vector<double> v;
    // v.push_back(v1_r); v.push_back(v1_i); v.push_back(v2_r); v.push_back(v2_i); v.push_back(v3_r); v.push_back(v3_i); v.push_back(v4_r); v.push_back(v4_i);
    auto U = ApplyGateFWithParamVec(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkArbitraryGateInterleaved, v);
    // auto U = ApplyGateF(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
    // VectorComplexFloatBoost::VectorPrintColumnMajor(U, std::cout);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(U, stateVector);
}

void CFLOBDDQuantumCircuit::ApplyCNOTGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = Matrix1234ComplexFloatBoost::MkCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controller, controlled);
        // Test, the total scale of a circuit!
        // kraus_operators[0] = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(kraus_operators[0], C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controlled, controller);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        // kraus_operators[0] = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(kraus_operators[0], C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplySwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(index1 != index2);
    
    if (index1 < index2)
    {
        auto C = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index1, index2);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto C = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index2, index1);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyiSwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(index1 != index2);
    
    if (index1 < index2)
    {
        auto C = Matrix1234ComplexFloatBoost::MkiSwapGate(stateVector.root->level, index1, index2);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto C = Matrix1234ComplexFloatBoost::MkiSwapGate(stateVector.root->level, index2, index1);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyCZGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = Matrix1234ComplexFloatBoost::MkCPGate(stateVector.root->level, controller, controlled, 1.0);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCPGate(stateVector.root->level, controlled, controller, 1.0);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyCPGate(long int controller, long int controlled, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = Matrix1234ComplexFloatBoost::MkCPGate(stateVector.root->level, controller, controlled, theta);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else
    {
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCPGate(stateVector.root->level, controlled, controller, theta);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    } 
}

void CFLOBDDQuantumCircuit::ApplyCSGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ApplyCPGate(controller, controlled, 0.5);
}

void CFLOBDDQuantumCircuit::ApplyPhaseShiftGate(unsigned int index, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto S = ApplyGateFWithParam(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, theta);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, stateVector); 
}

void CFLOBDDQuantumCircuit::ApplyTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto S = ApplyGateFWithParam(std::pow(2, stateVector.root->level-1), index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, 0.25);
    stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, stateVector); 
}

void CFLOBDDQuantumCircuit::ApplyCCNOTGate(long int controller1, long int controller2, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller1 != controlled);
    assert(controller2 != controlled);
    assert(controller1 != controller2);
    if (controller1 < controller2 && controller2 < controlled)
    {
        // a b c
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controller1, controller2, controlled);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (controller1 < controlled && controlled < controller2)
    {
        // a c b   
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller2);
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controller1, controlled, controller2);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (controller2 < controller1 && controller1 < controlled)
    {
        // b a c
        ApplyCCNOTGate(controller2, controller1, controlled);
    }
    else if (controller2 < controlled && controlled < controller1)
    {
        // b c a
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller1);
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controller2, controlled, controller1);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector); 
    }
    else if (controlled < controller1 && controller1 < controller2)
    {
        // c a b
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller2);
        // b a c
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controlled, controller1, controller2);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (controlled < controller2 && controller2 < controller1)
    {
        // c b a
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, controlled, controller1);
        // a b c
        auto C = Matrix1234ComplexFloatBoost::MkCCNOT(stateVector.root->level, std::pow(2, stateVector.root->level - 1), controlled, controller2, controller1);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyCSwapGate(long int controller, long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != index1);
    assert(controller != index2);
    assert(index1 != index2);
    
    if (controller < index1 && index1 < index2)
    {
        // a b c
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, controller, index1, index2);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (controller < index2 && index2 < index1)
    {
        // a c b   
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, controller, index2, index1);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (index1 < controller && controller < index2)
    {
        // b a c
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index1, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, index1, controller, index2);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (index1 < index2 && index2 < controller)
    {
        // b c a
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index1, controller);
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, index1, index2, controller);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (index2 < controller && controller < index1)
    {
        // c a b
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index2, controller);
        // b a c
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, index2, controller, index1);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
    else if (index2 < index1 && index1 < controller)
    {
        // c b a
        auto S = Matrix1234ComplexFloatBoost::MkSwapGate(stateVector.root->level, index2, controller);
        // a b c
        auto C = Matrix1234ComplexFloatBoost::MkCSwapGate(stateVector.root->level, index2, index1, controller);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
        C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
        stateVector = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, stateVector);
    }
}

void CFLOBDDQuantumCircuit::ApplyGlobalPhase(double phase)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto cos_v = boost::math::cos_pi(phase);
    auto sin_v = boost::math::sin_pi(phase);
    BIG_COMPLEX_FLOAT phase_complex(cos_v, sin_v);
    stateVector = phase_complex * stateVector;
}

long double CFLOBDDQuantumCircuit::GetProbability(std::map<unsigned int, int>& qubit_vals)
{
    // stateVector.print(std::cout);
    auto tmp = VectorComplexFloatBoost::VectorWithAmplitude(stateVector);
    std::string s(std::pow(2, tmp.root->level-1), 'X');
    for (unsigned int i = 0; i < numQubits; i++)
    {
        if (qubit_vals.find(i) != qubit_vals.end())
        {
            if (qubit_vals[i] == 0)
                s[i] = '0';
            else if (qubit_vals[i] == 1)
                s[i] = '1';   
        }
    }
    auto restricted = Matrix1234ComplexFloatBoost::MkRestrictMatrix(tmp.root->level, s);
    tmp = tmp * restricted;
    tmp.CountPaths();
    return VectorComplexFloatBoost::getNonZeroProbability(tmp);
}

std::string CFLOBDDQuantumCircuit::Measure() 
{
    auto tmp = VectorComplexFloatBoost::VectorWithAmplitude(stateVector);
    tmp.CountPaths();
    return VectorComplexFloatBoost::Sampling(tmp, true).substr(0, numQubits); 
}

unsigned long long int CFLOBDDQuantumCircuit::GetPathCount(long double prob)
{
    auto tmp = VectorComplexFloatBoost::VectorWithAmplitude(stateVector);
    tmp.CountPaths();
    return VectorComplexFloatBoost::GetPathCount(tmp, prob);  
}



void CFLOBDDQuantumCircuit::appendGateSeries(std::string name, std::vector<unsigned int> index, std::vector<double> vars, bool newChannel)
{
    // need to add some defaults
    if (newChannel) {
        circuitGates.appendChannel();
    }
    circuitGates.appendGate(name, index, vars);
}

CFLOBDD_COMPLEX_BIG normalize(CFLOBDD_COMPLEX_BIG c) {
    auto H = ApplyGateF(std::pow(2, c.root->level-1), 0, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
    CFLOBDD_COMPLEX_BIG c1 = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, c);
    CFLOBDD_COMPLEX_BIG c1_conj = Matrix1234ComplexFloatBoost::MatrixConjugate(c1);
    // VectorComplexFloatBoost::VectorPrintColumnHead(c1_conj, std::cout);
    c1_conj = Matrix1234ComplexFloatBoost::MatrixTranspose(c1_conj);
    auto mulres = Matrix1234ComplexFloatBoost::MatrixMultiplyV4(c1_conj, c1);
    auto resMap = mulres.root->rootConnection.returnMapHandle;
    // Maybe #BUGS here!
    double dimfactor = std::pow(double(2), double(std::pow(2, c.root->level-1)-1));
    // double dimfactor = 1;
    assert(resMap.Size() <= 2);
    BIG_COMPLEX_FLOAT amp;
    if(resMap.Size() == 2) {
        amp = (resMap[0] != 0) ? resMap[0] : resMap[1];
        // std::cout << amp.imag() << " " << amp.real() << std::endl;
        assert(abs(amp.imag()*dimfactor) < 1e-7 && amp.real() > 0);
        double factor = double(sqrt(amp.real()));
        // std::cout << "Normalize: " << factor << std::endl;
        c1 = (1/factor) * c1;
    } else {
        assert(abs(amp.imag()*dimfactor) < 1e-7 && abs(amp.real()*dimfactor) < 1e-7);
        // std::cout << "Nodistinction!!!" << std::endl;
        c1 = VectorComplexFloatBoost::NoDistinctionNode(c.root->level, 0);
    }
    return c1;
}

bool checkifzero(CFLOBDD_COMPLEX_BIG c) {
    // Maybe #BUGS here!
    double dimfactor = std::pow(double(2), double(std::pow(2, c.root->level-1)-1));
    // double dimfactor = 1;
    double threshold = 1e-7;
    auto resMap = c.root->rootConnection.returnMapHandle;
    for(int i = 0; i < resMap.Size(); i++) {
        if(abs(resMap[i].real()*dimfactor) + abs(resMap[i].imag()*dimfactor) > threshold) {
            // std::cout << abs(resMap[i].real()*dimfactor) + abs(resMap[i].imag()*dimfactor) << " ";
            return false;
        }
    }
    return true;
}

int CFLOBDDQuantumCircuit::ApplyGateSeries(int channelIdx)
{
    assert(channelIdx < circuitGates.numChannel);
    for(int i = 0; i < circuitGates.operations[channelIdx].size(); i++) {
        QuantumGate g = circuitGates.operations[channelIdx][i];
        if(g.name == "h") {
            ApplyHadamardGate(g.index[0]);
        } else if(g.name == "x") {
            ApplyNOTGate(g.index[0]);
        } else if(g.name == "y") {
            ApplyPauliYGate(g.index[0]);
        } else if(g.name == "z") {
            ApplyPauliZGate(g.index[0]);
        } else if(g.name == "s") {
            ApplySGate(g.index[0]);
        } else if(g.name == "t") {
            ApplyTGate(g.index[0]);
        } else if(g.name == "tdg") {
            ApplyU3Gate(g.index[0], 0, 0, -0.25);
        } else if(g.name == "u1") {
            // ApplyPhaseShiftGate
            ApplyU3Gate(g.index[0], 0, 0, g.vars[0]);
        } else if(g.name == "u2") {
            ApplyU3Gate(g.index[0], 0.5, g.vars[0], g.vars[1]);
        } else if(g.name == "u3") {
            ApplyU3Gate(g.index[0], g.vars[0], g.vars[1], g.vars[2]);
        } else if(g.name == "cx") {
            ApplyCNOTGate(g.index[0], g.index[1]);
        } else if(g.name == "cz") {
            ApplyCZGate(g.index[0], g.index[1]);
        } else if(g.name == "ccx") {
            ApplyCCNOTGate(g.index[0], g.index[1], g.index[2]);
        } else if(g.name == "ad") {
            // unnormalized
            if(g.index[1] == 1) {
                std::vector<double> v{1,0,0,0,0,0,sqrt(1-g.vars[0]),0};
                ApplyArbitraryGate(g.index[0], v);
                // ApplyNOTGate(g.index[0]);
            } else if(g.index[1] == 2) {
                std::vector<double> v{0,0,sqrt(g.vars[0]),0,0,0,0,0};
                ApplyArbitraryGate(g.index[0], v);
                // ApplyNOTGate(g.index[0]);
            } else {
                // CP CSWAP CS
                assert(1 == 0);
            }
        } else if(g.name == "measure") {
            // normalized
            if(g.index[1] == 1) {
                std::vector<double> v{1,0,0,0,0,0, 0,0};
                // VectorComplexFloatBoost::VectorPrintColumnHead(stateVector, std::cout);
                ApplyArbitraryGate(g.index[0], v);
                // VectorComplexFloatBoost::VectorPrintColumnHead(stateVector, std::cout);
                stateVector = normalize(stateVector);
                if(checkifzero(stateVector)) {
                    return -1;
                }
            } else if(g.index[1] == 2) {
                std::vector<double> v{0,0,0,0,0,0,1,0};
                ApplyArbitraryGate(g.index[0], v);
                stateVector = normalize(stateVector);
                if(checkifzero(stateVector)) {
                    return -1;
                }
            } else {
                assert(2 == 0);
            }
        } else if(g.name == "partial") {
            // reset = partial trace + append
            
        } else {
            assert(3 == 0);
        }
    }
    return 0;
}

unsigned int CFLOBDDQuantumCircuit::reachability()
{
    // Initialize queue, projector and cnt
    state_queue.clear();
    for(int i = 0; i < stateProjector.size(); i++) {
        state_queue.push_back(i);
    }
    unsigned int cnt = stateProjector.size();
    assert(cnt != 0 && state_queue.size() != 0);
    unsigned int d;
    if(realQubits == 0) {
        d = 1 << numQubits;
    } else {
        d = 1 << realQubits;
        // std::cout << d << std::endl;
    }
    /// main loop: if queue is not empty and cnt < d
    while(state_queue.empty() == false && cnt <= d) {
        // if(cnt % 10 == 0) {
        //     std::cout << cnt << std::endl;
        // }
        // #DEBUG: Avoid long looping!!
        // if(cnt >= 30) break;
        auto cur_state_idx = state_queue.front();
        state_queue.pop_front();
        std::vector<CFLOBDD_COMPLEX_BIG> extended_states;
        for(int j = 0; j < circuitGates.numChannel; j++) {
            //// make sure transpose value!!
            stateVector = stateProjector[cur_state_idx];
            //// Here, may introduce much complexity!!
            int ck = ApplyGateSeries(j);
            // VectorComplexFloatBoost::VectorPrintColumnHead(stateVector, std::cout);
            if(ck == 0) {
                extended_states.push_back(stateVector);
            }
        }
        //// arbitrary clear
        stateVector = VectorComplexFloatBoost::MkBasisVector(1, 0);
        for(int j = 0; j < extended_states.size(); j++) {
            CFLOBDD_COMPLEX_BIG e_state = extended_states[j];
            stateVector = e_state;
            ApplyProjectorToState();
            CFLOBDD_COMPLEX_BIG tmp_state = e_state + (-1)*stateVector;
            // Normalization is necessary!!!
            // Here, put normalization into checkifzero block!
            
            // Here, introduce some epsilon
            if(checkifzero(tmp_state) == false) {
                // Need some benchmarks!!
                tmp_state = normalize(tmp_state);
                state_queue.push_back(cnt);
                cnt++;
                stateProjector.push_back(tmp_state);
            }
        }
    }
    return cnt;
}

std::string getIndexStringFromBitString(std::string s){
			std::string ans(s.length() * 2, '0');
			unsigned int j = 0;
			for (int i = 0; i < s.length(); i++){
				ans[j] = s[i];
				ans[j + 1] = s[i];
				j+=2;
			}
			return ans;
		}

void CFLOBDDQuantumCircuit::setBasicStateVector(std::string s)
{
    unsigned int level = ceil(log2(numQubits));
    stateVector = VectorComplexFloatBoost::MkBasisVector(level, s);
    stateVector = VectorComplexFloatBoost::VectorToMatrixInterleaved(stateVector);
    return;
}

void CFLOBDDQuantumCircuit::setRealQubits(unsigned int q)
{
    realQubits = q;
    return;
}

void CFLOBDDQuantumCircuit::setInitGate()
{
    unsigned int level = ceil(log2(numQubits));
    auto H = ApplyGateF(std::pow(2, level), 0, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
    stateVector = H;
}


void CFLOBDDQuantumCircuit::setProjector(std::vector<double> real, std::vector<double> imag, std::vector<int> index)
{
    assert(real.size() == imag.size() && real.size() == index.size());
    assert(real.empty() == false);
    unsigned int level = ceil(log2(numQubits));
    auto ProjV = VectorComplexFloatBoost::MkBasisVector(level, index[0]);
    BIG_COMPLEX_FLOAT amp(real[0], imag[0]);
    ProjV = amp * ProjV;
    for(int i = 1; i < real.size(); i++) {
        ProjV = ProjV + BIG_COMPLEX_FLOAT(real[i], imag[i]) * VectorComplexFloatBoost::MkBasisVector(level, index[i]);
    }
    ProjV = VectorComplexFloatBoost::VectorToMatrixInterleaved(ProjV);
    stateProjector.push_back(ProjV);
    // auto tmp = Matrix1234ComplexFloatBoost::MatrixTranspose(ProjV);
    // stateProjector.push_back(tmp);
}

void CFLOBDDQuantumCircuit::setProjectorFromS()
{
    // set projector from statevector
    stateProjector.push_back(stateVector);
    // auto tmp = Matrix1234ComplexFloatBoost::MatrixTranspose(stateVector);
    // stateProjector.push_back(tmp);
}

void CFLOBDDQuantumCircuit::ApplyProjectorToEntangle()
{
    assert(stateProjector.empty() == false);
    unsigned int level = ceil(log2(numQubits));
    stateVector = VectorComplexFloatBoost::NoDistinctionNode(level+1, 0);
    for(int i = 0; i <= ((1 << numQubits)-1); i++) {
    // for(int i = 0; i <= 3; i++) {
        auto indexVec = VectorComplexFloatBoost::MkBasisVector(level, i);
        indexVec = VectorComplexFloatBoost::VectorToMatrixInterleaved(indexVec);
        auto tmpVec = Matrix1234ComplexFloatBoost::MatrixTranspose(indexVec);
        CFLOBDD_COMPLEX_BIG ithColP = VectorComplexFloatBoost::NoDistinctionNode(level+1, 0);
        //// make an all-zero vector
        for(int j = 0; j < stateProjector.size(); j++) {
            // calculate: \ket{j}\bra{j}\ket{i}
            // std::cout << tmpVec.root->level << " " << stateProjector[j].root->level << std::endl;
            auto tmp = Matrix1234ComplexFloatBoost::MatrixMultiplyV4(tmpVec, stateProjector[j]);
            // return value of tmp should be only one non-zero.
            // std::cout << tmp.root->rootConnection.returnMapHandle.Size() << std::endl;
            assert(tmp.root->rootConnection.returnMapHandle.Size() <= 2);
            auto resMap = tmp.root->rootConnection.returnMapHandle;
            BIG_COMPLEX_FLOAT amp;
            if(resMap.Size() == 2)
                amp = (resMap[0] != 0) ? resMap[0] : resMap[1];
            else
                amp = resMap[0];
            amp = conj(amp);
            ithColP = ithColP + (amp * stateProjector[j]);
        }
        CFLOBDD_COMPLEX_BIG ChoiUi = Matrix1234ComplexFloatBoost::KroneckerProduct2Vocs(ithColP, indexVec);
        stateVector = stateVector + ChoiUi;
    }
}


void CFLOBDDQuantumCircuit::ApplyProjectorToState()
{
    assert(stateProjector.empty() == false);
    unsigned int level = ceil(log2(numQubits));
    CFLOBDD_COMPLEX_BIG res = VectorComplexFloatBoost::NoDistinctionNode(level+1, 0);
    auto tmpVec = Matrix1234ComplexFloatBoost::MatrixTranspose(stateVector);
    tmpVec = Matrix1234ComplexFloatBoost::MatrixConjugate(tmpVec);
    for(int i = 0; i < stateProjector.size(); i++) {
        auto tmp = Matrix1234ComplexFloatBoost::MatrixMultiplyV4(tmpVec, stateProjector[i]);
        assert(tmp.root->rootConnection.returnMapHandle.Size() <= 2);
        auto resMap = tmp.root->rootConnection.returnMapHandle;
        BIG_COMPLEX_FLOAT amp;
        if(resMap.Size() == 2)
            amp = (resMap[0] != 0) ? resMap[0] : resMap[1];
        else
            amp = resMap[0];
        amp = conj(amp);
        // std::cout << "ApplyProjectorToState: " << amp << std::endl;
        res = res + (amp * stateProjector[i]);
    }
    stateVector = res;
    return;
}
// void CFLOBDDQuantumCircuit::setSuperPosVector(std::string s)
// {
//     return;
// }

void CFLOBDDQuantumCircuit::pushStateToCache()
{
    stateVectorCache.push_back(stateVector);
}

void CFLOBDDQuantumCircuit::addCacheToState()
{
    stateVector = stateVector + stateVectorCache.back();
    stateVectorCache.pop_back();
}

void CFLOBDDQuantumCircuit::printStateColMajor()
{
    VectorComplexFloatBoost::VectorPrintColumnMajor(stateVector, std::cout);
}

void CFLOBDDQuantumCircuit::printStateColHead()
{
    VectorComplexFloatBoost::VectorPrintColumnHead(stateVector, std::cout);
}

void CFLOBDDQuantumCircuit::printYield()
{
    stateVector.PrintYield(&std::cout);
    std::cout << std::endl;
}

void CFLOBDDQuantumCircuit::print()
{
    std::cout << stateVector << std::endl;
}

void CFLOBDDQuantumCircuit::printRV(std::string type)
{
    if (type == "state") {
        std::cout << stateVector.root->rootConnection.returnMapHandle << std::endl;
        std::cout << "Number of values: " << 
                    stateVector.root->rootConnection.returnMapHandle.mapContents->mapArray.size() << std::endl;
    } else {
        std::cout << "test" << std::endl;
        // std::cout << kraus_operators[0].root->rootConnection.returnMapHandle << std::endl;
        // std::cout << "Number of values: " << 
        //             kraus_operators[0].root->rootConnection.returnMapHandle.mapContents->mapArray.size() << std::endl;
    }
}

unsigned int CFLOBDDQuantumCircuit::printSize(std::string type)
{
    if (type == "state") {
        unsigned int nodeCount, edgeCount;
        unsigned int returnEdgesCount, returnEdgesObjCount;
        stateVector.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
        std::cout << "Number of values: " << 
                    stateVector.root->rootConnection.returnMapHandle.mapContents->mapArray.size() << std::endl;
        std::cout << "Nodes & Edges in stateVector: " << nodeCount << ", " << edgeCount << std::endl;
    } else if (type == "projector") {
        unsigned int nodeCount, edgeCount;
        unsigned int returnEdgesCount, returnEdgesObjCount;
        unsigned int nodeCountTotal = 0, edgeCountTotal = 0;
        for(int i = 0; i < stateProjector.size(); i++) {
            stateProjector[i].CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
            nodeCountTotal += nodeCount;
            edgeCountTotal += edgeCount;
            // std::cout << i << " " << nodeCount << " " << edgeCount << std::endl;
            // std::cout << "Number of values " << i << ": "
            //         stateProjector[i].root->rootConnection.returnMapHandle.mapContents->mapArray.size() << std::endl;
        }
        std::cout << "Nodes & Edges in projectors: " << nodeCountTotal << ", " << edgeCountTotal << std::endl;
    }
    
    return 0;
}

void CFLOBDDQuantumCircuit::printProjector()
{
    if(stateProjector.size() == 0) {
        std::cout << "Empty space!\n";
        return;
    }
    for(int i = 0; i < stateProjector.size(); i++) {
        VectorComplexFloatBoost::VectorPrintColumnHead(stateProjector[i], std::cout);
    }
    return;
}

void CFLOBDDQuantumCircuit::test()
{
    stateVector = VectorComplexFloatBoost::MkBasisVector(2, 5);
    std::cout << stateVector << std::endl;
    CFLOBDD_COMPLEX_BIG tmpVector = VectorComplexFloatBoost::MkBasisVector(2, 1);
    stateVector = stateVector + tmpVector;
    std::cout << stateVector << std::endl;
    stateVector = stateVector + tmpVector;
    std::cout << stateVector << std::endl;
    // std::vector<double> vec;
    // vec.push_back(0);
    // vec.push_back(0);
    // vec.push_back(0.25);
    // auto U = Matrix1234ComplexFloatBoost::MkU3GateInterleaved(1, vec);
    // std::cout << U << std::endl;

    // double cos_theta = boost::math::cos_pi(0/2);
    // double sin_theta = boost::math::sin_pi(0/2);
    // double cos_lambda = boost::math::cos_pi(0.25);
    // double sin_lambda = boost::math::sin_pi(0.25);
    // double cos_phi = boost::math::cos_pi(0);
    // double sin_phi = boost::math::sin_pi(0);
    // BIG_COMPLEX_FLOAT exp_lambda(cos_lambda, sin_lambda);
    // BIG_COMPLEX_FLOAT exp_phi(cos_phi, sin_phi);
    // auto a = -exp_lambda*sin_theta;
    // auto b = exp_phi*sin_theta;
    // std::cout << a << " " << b << std::endl;
    // std::cout << (a==b) << std::endl;
    // std::cout << typeid(a).name() << " " << typeid(b).name() << std::endl;
    // int size = stateVector.root->rootConnection.returnMapHandle.Size();
    // std::cout << size << std::endl;
    // for(int i=0; i<size; i++) {
    //     std::cout << "(" << stateVector.root->rootConnection.returnMapHandle[i].real() << ", " 
    //     << stateVector.root->rootConnection.returnMapHandle[i].imag() << ") ";
    // }
    // std::cout << std::endl;

    // auto tmp = VectorComplexFloatBoost::VectorWithAmplitude(stateVector);
    // for(int i=0; i<size; i++) {
    //     std::cout << "(" << tmp.root->rootConnection.returnMapHandle[i].real() << ", " 
    //     << tmp.root->rootConnection.returnMapHandle[i].imag() << ") ";
    // }
    // std::cout << std::endl;

    // Matrix1234ComplexFloatBoost::MatrixConjugate(tmp);
    // for(int i=0; i<size; i++) {
    //     std::cout << "(" << tmp.root->rootConnection.returnMapHandle[i].real() << ", " 
    //     << tmp.root->rootConnection.returnMapHandle[i].imag() << ") ";
    // }
    // std::cout << std::endl;

    // BIG_COMPLEX_FLOAT c{0.2,1.3};
    // std::cout << c << std::endl;
    return;
}

