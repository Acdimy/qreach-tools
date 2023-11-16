#include "quantum_circuit.h"
#include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include "cflobdd/CFLOBDD/vector_complex_float_boost.h"
#include <random>

QuantumCircuit::QuantumCircuit(unsigned int numQubits, int seed) :  numQubits (numQubits) 
{
    mt.seed(seed);
    srand(seed);
    hadamard_count = 0;
}
QuantumCircuit::QuantumCircuit() :  numQubits (0), hadamard_count (0) {}
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

void CFLOBDDQuantumCircuit::ApplyGateSeries(int channelIdx)
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
        } else if(g.name == "u1") {
            ApplyU3Gate(g.index[0], 0, 0, g.vars[0]);
        } else if(g.name == "u2") {
            ApplyU3Gate(g.index[0], 0.5, g.vars[0], g.vars[1]);
        } else if(g.name == "u3") {
            ApplyU3Gate(g.index[0], g.vars[0], g.vars[1], g.vars[2]);
        } else if(g.name == "cx") {
            ApplyCNOTGate(g.index[0], g.index[1]);
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
            }
        } else if(g.name == "tr") {
            
        } else if(g.name == "dp") {
            
        } else {
            assert(1 == 0);
        }
    }
}

CFLOBDD_COMPLEX_BIG nomalize(CFLOBDD_COMPLEX_BIG c) {
    auto H = ApplyGateF(std::pow(2, c.root->level-1), 0, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
    CFLOBDD_COMPLEX_BIG c1 = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, c);
    CFLOBDD_COMPLEX_BIG c1_conj = Matrix1234ComplexFloatBoost::MatrixTranspose(c1);
    c1_conj = Matrix1234ComplexFloatBoost::MatrixConjugate(c1_conj);
    auto mulres = Matrix1234ComplexFloatBoost::MatrixMultiplyV4(c1_conj, c1);
    auto resMap = mulres.root->rootConnection.returnMapHandle;
    assert(resMap.Size() <= 2);
    BIG_COMPLEX_FLOAT amp;
    if(resMap.Size() == 2) {
        amp = (resMap[0] != 0) ? resMap[0] : resMap[1];
        assert(amp.imag() == 0 && amp.real() > 0);
        double factor = double(sqrt(amp.real()));
        c1 = (1/factor) * c1;
    }
    return c1;
}

bool checkifzero(CFLOBDD_COMPLEX_BIG c) {
    double threshold = 1e-9;
    auto resMap = c.root->rootConnection.returnMapHandle;
    for(int i = 0; i < resMap.Size(); i++) {
        if(abs(resMap[i].real()) + abs(resMap[i].imag()) > threshold) {
            return false;
        }
    }
    return true;
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
    unsigned int d = 2 << numQubits;
    /// main loop: if queue is not empty or cnt < d
    while(state_queue.empty() == false && cnt <= d) {
        // std::cout << cnt << std::endl;
        auto cur_state_idx = state_queue.front();
        state_queue.pop_front();
        std::vector<CFLOBDD_COMPLEX_BIG> extended_states;
        for(int j = 0; j < circuitGates.numChannel; j++) {
            //// make sure transpose value!!
            stateVector = stateProjector[cur_state_idx];
            //// Here, may introduce much complexity!!
            ApplyGateSeries(j);
            extended_states.push_back(stateVector);
        }
        //// arbitrary clear
        stateVector = VectorComplexFloatBoost::MkBasisVector(1, 0);
        for(int j = 0; j < extended_states.size(); j++) {
            CFLOBDD_COMPLEX_BIG e_state = extended_states[j];
            stateVector = e_state;
            ApplyProjectorToState();
            CFLOBDD_COMPLEX_BIG tmp_state = e_state + (-1)*stateVector;
            tmp_state = nomalize(tmp_state);
            // Here, introduce some epsilon
            if(checkifzero(tmp_state) == false) {
                // Need some benchmarks!!
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

void CFLOBDDQuantumCircuit::setInitGate()
{
    unsigned int level = ceil(log2(numQubits));
    auto H = ApplyGateF(std::pow(2, level), 0, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
    stateVector = H;
}

// void CFLOBDDQuantumCircuit::setMaximalEntangledVec()
// {
//     // Make stateVector maximally entangled in space of 2*numQubits.
//     // A uniform phase is omitted
//     unsigned int level = ceil(log2(numQubits)) + 1;
//     stateVector = VectorComplexFloatBoost::MkBasisVector(level, 0);
//     for(int j = 1; j <= ((1 << numQubits)-1); j++) {
//         int index = j + (j << numQubits);
//         // std::cout << index << " ";
//         CFLOBDD_COMPLEX_BIG tmpVector = VectorComplexFloatBoost::MkBasisVector(level, index);
//         stateVector = stateVector + tmpVector;
//     }
//     // printSize("state");
//     stateVector = VectorComplexFloatBoost::VectorToMatrixInterleaved(stateVector);
//     // printSize("state");
//     // std::cout << std::endl;
// }

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
    unsigned int nodeCount, edgeCount;
    unsigned int returnEdgesCount, returnEdgesObjCount;
    if (type == "state") {
        stateVector.CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
        std::cout << "Number of values: " << 
                    stateVector.root->rootConnection.returnMapHandle.mapContents->mapArray.size() << std::endl;
    } else {
        std::cout << "test" << std::endl;
        // kraus_operators[0].CountNodesAndEdges(nodeCount, edgeCount, returnEdgesCount, returnEdgesObjCount);
    }
    std::cout << nodeCount << ", " << edgeCount << std::endl;
    return edgeCount;
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


// *******************
// Weighted BDD
// *******************

#include "cflobdd/CFLOBDD/wvector_complex_fb_mul.h"
#include "cflobdd/CFLOBDD/weighted_cross_product.h"
#include "cflobdd/CFLOBDD/weighted_cross_product_bdd.h"

WeightedBDDQuantumCircuit::WeightedBDDQuantumCircuit(unsigned int numQubits, int seed) : QuantumCircuit(numQubits, seed)
{
    // Initialize
    WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable();
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable_Ann();
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitIdentityNodeTable();	
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitReduceCache();
	WeightedMatrix1234ComplexFloatBoostMul::Matrix1234Initializer();
	WeightedVectorComplexFloatBoostMul::VectorInitializer();
	InitWeightedPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();

    WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitLeafNodes();
	InitWeightedBDDPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
    //
    std::string s(2 * numQubits, '0');
    stateVector = WeightedVectorComplexFloatBoostMul::MkBasisVector(2 * numQubits, s, 0);
}

WeightedBDDQuantumCircuit::WeightedBDDQuantumCircuit()
{
    // Initialize
    WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable();
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitNoDistinctionTable_Ann();
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitIdentityNodeTable();	
	WeightedCFLOBDDNodeHandleT<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitReduceCache();
	WeightedMatrix1234ComplexFloatBoostMul::Matrix1234Initializer();
	WeightedVectorComplexFloatBoostMul::VectorInitializer();
	InitWeightedPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();

    WeightedBDDNodeHandle<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>::InitLeafNodes();
	InitWeightedBDDPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
    numQubits = 0;
    //
}

WeightedBDDQuantumCircuit::~WeightedBDDQuantumCircuit()
{
    DisposeOfWeightedPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
    DisposeOfWeightedBDDPairProductCache<BIG_COMPLEX_FLOAT, std::multiplies<BIG_COMPLEX_FLOAT>>();
}

void WeightedBDDQuantumCircuit::setNumQubits(unsigned int num)
{
    numQubits = num;
    stateVector = WeightedVectorComplexFloatBoostMul::MkBasisVector(numQubits, 0, 0);
}

WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ApplyGateF(unsigned int n, unsigned int i, WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(*f)(unsigned int, int))
{
    WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL H = f(2, 0);
    if (i == 0)
    {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n-1), 0);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(H, I); 
    }
    else if (i == n - 1)
    {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n-1), 0);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(I, H);
    }
    else {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I1 = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * i, 0);
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I2 = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n - i - 1), 0);
        auto T = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(H, I2);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(I1, T);
    }
}

WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL ApplyGateFWithParam(unsigned int n, unsigned int i, WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL(*f)(unsigned int, double, int), double theta)
{
    WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL H = f(2, theta, 0);
    if (i == 0)
    {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n-1), 0);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(H, I); 
    }
    else if (i == n - 1)
    {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n-1), 0);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(I, H);
    }
    else {
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I1 = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * i, 0);
        WEIGHTED_CFLOBDD_COMPLEX_FLOAT_BOOST_MUL I2 = WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved(2 * (n - i - 1), 0);
        auto T = WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(H, I2);
        return WeightedMatrix1234ComplexFloatBoostMul::KroneckerProduct2Vocs(I1, T);
    }
}


void WeightedBDDQuantumCircuit::ApplyIdentityGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto H = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkIdRelationInterleaved);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(H, stateVector);
}

void WeightedBDDQuantumCircuit::ApplyHadamardGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto H = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkWalshInterleaved);
    // H.print(std::cout);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(H, stateVector);
    // stateVector.print(std::cout);
    // std::cout << std::endl;
    hadamard_count++;
}

void WeightedBDDQuantumCircuit::ApplyNOTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    // stateVector.print(std::cout);
    auto X = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkNegationMatrixInterleaved);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(X, stateVector);
    // stateVector.print(std::cout);
}

void WeightedBDDQuantumCircuit::ApplyPauliYGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto Y = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkPauliYGate);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(Y, stateVector);
}

void WeightedBDDQuantumCircuit::ApplyPauliZGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto Z = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkPauliZGate);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(Z, stateVector);
}

void WeightedBDDQuantumCircuit::ApplySGate(unsigned int index)
{
    auto S = ApplyGateF(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkSGate);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, stateVector);
}

void WeightedBDDQuantumCircuit::ApplyCNOTGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCNOT(2 * numQubits, numQubits, controller, controlled, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
        // stateVector.print(std::cout);
    }
    else
    {
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCNOT(2 * numQubits, numQubits, controlled, controller, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplySwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(index1 != index2);
    
    if (index1 < index2)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index1, index2, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index2, index1, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyiSwapGate(long int index1, long int index2)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(index1 != index2);
    
    if (index1 < index2)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkiSwapGate(2 * numQubits, index1, index2, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkiSwapGate(2 * numQubits, index2, index1, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyCZGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(2 * numQubits, controller, controlled, 1.0, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else
    {
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(2 * numQubits, controlled, controller, 1.0, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyCPGate(long int controller, long int controlled, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    assert(controller != controlled);
    
    if (controller < controlled)
    {
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(2 * numQubits, controller, controlled, theta, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else
    {
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCPGate(2 * numQubits, controlled, controller, theta, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    } 
}

void WeightedBDDQuantumCircuit::ApplyCSGate(long int controller, long int controlled)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    ApplyCPGate(controller, controlled, 0.5);
}

void WeightedBDDQuantumCircuit::ApplyPhaseShiftGate(unsigned int index, double theta)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto S = ApplyGateFWithParam(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkPhaseShiftGate, theta);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, stateVector); 
}

void WeightedBDDQuantumCircuit::ApplyTGate(unsigned int index)
{
    if (checkForInit(numQubits) == false)
    {
        std::cout << "Number of Qubits is unset" << std::endl;
        abort();   
    }
    auto S = ApplyGateFWithParam(numQubits, index, WeightedMatrix1234ComplexFloatBoostMul::MkPhaseShiftGate, 0.25);
    stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, stateVector); 
}

void WeightedBDDQuantumCircuit::ApplyCCNOTGate(long int controller1, long int controller2, long int controlled)
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
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controller1, controller2, controlled, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (controller1 < controlled && controlled < controller2)
    {
        // a c b   
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller2, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controller1, controlled, controller2, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (controller2 < controller1 && controller1 < controlled)
    {
        // b a c
        ApplyCCNOTGate(controller2, controller1, controlled);
    }
    else if (controller2 < controlled && controlled < controller1)
    {
        // b c a
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller1, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controller2, controlled, controller1, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector); 
    }
    else if (controlled < controller1 && controller1 < controller2)
    {
        // c a b
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller2, 0);
        // b a c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controlled, controller1, controller2, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (controlled < controller2 && controller2 < controller1)
    {
        // c b a
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, controlled, controller1, 0);
        // a b c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCCNOT(2 * numQubits, controlled, controller2, controller1, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyCSwapGate(long int controller, long int index1, long int index2)
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
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, controller, index1, index2, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (controller < index2 && index2 < index1)
    {
        // a c b   
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, controller, index2, index1, 0);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (index1 < controller && controller < index2)
    {
        // b a c
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index1, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, index1, controller, index2, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (index1 < index2 && index2 < controller)
    {
        // b c a
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index1, controller, 0);
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, index1, index2, controller, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (index2 < controller && controller < index1)
    {
        // c a b
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index2, controller, 0);
        // b a c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, index2, controller, index1, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
    else if (index2 < index1 && index1 < controller)
    {
        // c b a
        auto S = WeightedMatrix1234ComplexFloatBoostMul::MkSwapGate(2 * numQubits, index2, controller, 0);
        // a b c
        auto C = WeightedMatrix1234ComplexFloatBoostMul::MkCSwapGate(2 * numQubits, index2, index1, controller, 0);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, S);
        C = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(S, C);
        stateVector = WeightedMatrix1234ComplexFloatBoostMul::MatrixMultiplyV4(C, stateVector);
    }
}

void WeightedBDDQuantumCircuit::ApplyGlobalPhase(double phase)
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

long double WeightedBDDQuantumCircuit::GetProbability(std::map<unsigned int, int>& qubit_vals)
{
    auto tmp = stateVector;
    std::string s(numQubits, 'X');
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
    auto restricted = WeightedMatrix1234ComplexFloatBoostMul::MkRestrictMatrix(2 * numQubits, s, 0);
    tmp = tmp * restricted;
    tmp.ComputeWeightOfPathsAsAmpsToExits();
    return WeightedVectorComplexFloatBoostMul::getNonZeroProbability(tmp);
}

std::string WeightedBDDQuantumCircuit::Measure() 
{
    auto tmp = stateVector;
    // tmp.print(std::cout);
    tmp.ComputeWeightOfPathsAsAmpsToExits();
    std::uniform_real_distribution<double> dis(0.0, 1.0);
    return WeightedVectorComplexFloatBoostMul::Sampling(tmp, true, mt, dis).substr(0, numQubits); 
}

unsigned long long int WeightedBDDQuantumCircuit::GetPathCount(long double prob)
{
    std::cout << "Error! Operation not supported in WBDDs" << std::endl;
    abort();
} 






