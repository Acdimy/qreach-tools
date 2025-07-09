#ifndef _QUANTUM_OPERATION
#define _QUANTUM_OPERATION

#include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include "cflobdd/CFLOBDD/vector_complex_float_boost.h"
#include <random>
#include <queue>
#include <vector>
#include <stdexcept>
#include <memory>
#include <algorithm>
#include <cctype>
#include <optional>
#include <cmath>
#include <iomanip>
#include <sstream>

using namespace CFL_OBDD;

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

std::string stringPadding(std::string str, unsigned int length) {
    if (str.length() >= length) {
        return str;
    }
    return str + std::string(length - str.length(), '0');
}

// 1-norm
bool checkifzero(CFLOBDD_COMPLEX_BIG c) {
    double threshold = 1e-8;
    auto resMap = c.root->rootConnection.returnMapHandle;
    if(resMap.Size() == 0) {
        return true;
    }
    auto sum = abs(resMap[0].real()) + abs(resMap[0].imag());
    for(int i = 1; i < resMap.Size(); i++) {
        // Hide an inequality!
        sum += (abs(resMap[i].real()) + abs(resMap[i].imag()));
        // std::cout << "checkifzero: " << sum << std::endl;
        if(sum > threshold) {
            return false;
        }
    }
    return true;
}

std::string toLower(const std::string& input) {
    std::string result = input;
    std::transform(result.begin(), result.end(), result.begin(),
                   [](unsigned char c){ return std::tolower(c); });
    return result;
}

class PauliString {
    unsigned int length;
    bool sign;
    std::vector<unsigned int> pauliList; // 0: I, 1: X, 2: Y, 3: Z
};

/*
Quantum term: two types of atomic items.
1. A quantum gate, the same as QReach
2. A CFLOBDD item.
*/
class QuantumTerm {
    /* data */
    public:
    // bool type;
    unsigned int qNum;
    
    CFLOBDD_COMPLEX_BIG content;
    public:
    virtual ~QuantumTerm() {}
    virtual bool getType() const = 0;
    virtual std::unique_ptr<QuantumTerm> clone() const = 0;
};

class QuantumGateTerm : public QuantumTerm {
    public:
    std::string name;
    std::vector<unsigned int> index;
    std::vector<double> vars;
    unsigned int level;
    bool isConcret = false;
    bool zeroOperator = false;
    bool ideOperator = false;
    // construct as a gate
    QuantumGateTerm() {}
    // QuantumGateTerm(unsigned int qubit) {this->qNum = qubit;}
    QuantumGateTerm(std::string nam, std::vector<unsigned int> idx, std::vector<double> pars, unsigned int qNum)  {
        assert((qNum & (qNum - 1)) == 0 && qNum != 0);
        // type = true;
        this->name = nam;
        this->index = idx;
        this->vars = pars;
        // Make sure the qNum is exp.
        this->level = ceil(log2(qNum)) + 1;
        this->qNum = std::pow(2, this->level-1);
        content = VectorComplexFloatBoost::NoDistinctionNode(1, 0);
    }
    QuantumGateTerm(bool setconstant) {
        if (setconstant) {
            this->ideOperator = 1;
        } else {
            this->zeroOperator = 1;
        }
    }
    bool getType() const override {return true;}
    std::unique_ptr<QuantumTerm> clone() const override {
        return std::make_unique<QuantumGateTerm>(*this);
    }

    CFLOBDD_COMPLEX_BIG concretize() const {
        std::string name = toLower(this->name);
        unsigned index = this->index[0];
        CFLOBDD_COMPLEX_BIG res;
        // std::pow(2, this->content.root->level-1);
        // unsigned int level = ceil(log2(numQubits)) + 1;
        if (name == "x") {
            auto X = ApplyGateF(this->qNum, index, Matrix1234ComplexFloatBoost::MkNegationMatrixInterleaved);
            res = X;
        } else if (name == "y") {
            auto Y = ApplyGateF(this->qNum, index, Matrix1234ComplexFloatBoost::MkPauliYMatrixInterleaved);
            res = Y;
        } else if (name == "z") {
            auto Z = ApplyGateF(this->qNum, index, Matrix1234ComplexFloatBoost::MkPauliZMatrixInterleaved);
            res = Z;
        } else if (name == "h") {
            auto H = ApplyGateF(this->qNum, index, Matrix1234ComplexFloatBoost::MkWalshInterleaved);
            res = H;
        } else if (name == "i") {
            auto H = ApplyGateF(this->qNum, index, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
            res = H;
        } else if (name == "s") {
            auto S = ApplyGateF(this->qNum, index, Matrix1234ComplexFloatBoost::MkSGateInterleaved);
            res = S;
        } else if (name == "t") {
            auto S = ApplyGateFWithParam(this->qNum, index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, 0.25);
            res = S; 
        } else if (name == "p") {
            double theta = this->vars[0];
            auto S = ApplyGateFWithParam(this->qNum, index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, theta);
            res = S;
        } else if (name == "cx") {
            assert(this->qNum && (this->qNum & (this->qNum - 1)) == 0);
            unsigned int controller = this->index[0];
            unsigned int controlled = this->index[1];
            unsigned int state_level = ceil(log2(this->qNum)) + 1;
            assert(controller != controlled);
    
            if (controller < controlled)
            {
                auto C = Matrix1234ComplexFloatBoost::MkCNOT(state_level, std::pow(2, state_level - 1), controller, controlled);
                res = C;
            }
            else
            {
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, controlled, controller);
                auto C = Matrix1234ComplexFloatBoost::MkCNOT(state_level, std::pow(2, state_level - 1), controlled, controller);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = C;
            }
        } else if (name == "u3") {
            assert(this->vars.size() == 3);
            std::vector<double> v;
            double theta = this->vars[0];
            double phi = this->vars[1];
            double lambda = this->vars[2];
            v.push_back(theta); v.push_back(phi); v.push_back(lambda);
            auto U = ApplyGateFWithParamVec(this->qNum, index, Matrix1234ComplexFloatBoost::MkU3GateInterleaved, v);
            res = U;
        } else if (name == "meas0") {
            std::vector<double> v{1,0,0,0,0,0,0,0};
            auto U = ApplyGateFWithParamVec(this->qNum, index, Matrix1234ComplexFloatBoost::MkArbitraryGateInterleaved, v);
            res = U;
        } else if (name == "meas1") {
            std::vector<double> v{0,0,0,0,0,0,1,0};
            auto U = ApplyGateFWithParamVec(this->qNum, index, Matrix1234ComplexFloatBoost::MkArbitraryGateInterleaved, v);
            res = U;
        } else if (name == "reset0") {
            // Reset the index qubit to |0>
            std::vector<double> v{0,0,1,0,0,0,0,0};
            auto U = ApplyGateFWithParamVec(this->qNum, index, Matrix1234ComplexFloatBoost::MkArbitraryGateInterleaved, v);
            res = U;
        } else if (name == "resetAll") {
            CFLOBDD_COMPLEX_BIG U = VectorComplexFloatBoost::NoDistinctionNode(ceil(log2(this->qNum)), 1);
            U = VectorComplexFloatBoost::VectorToMatrixInterleaved(U);
            U = Matrix1234ComplexFloatBoost::MatrixConjugate(U);
            U = Matrix1234ComplexFloatBoost::MatrixTranspose(U);
            res = U;
        } else if (name == "swap") {
            
        } else if (name == "iswap") {
            
        } else if (name == "cz") {
            // Assert qNum is a power of 2
            assert(this->qNum && (this->qNum & (this->qNum - 1)) == 0);
            unsigned int controller = this->index[0];
            unsigned int controlled = this->index[1];
            unsigned int state_level = ceil(log2(this->qNum)) + 1;
            assert(controller != controlled);
    
            if (controller < controlled)
            {
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(state_level, controller, controlled, 1.0);
                res = C;
            }
            else
            {
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, controlled, controller);
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(state_level, controlled, controller, 1.0);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = C;
            }
        } else if (name == "cp") {

        } else if (name == "cs") {
            
        } /* CCNOT, CSWAP, GlobalPhase */
        else {
            std::cout << "Unknown quantum gate: " << name << std::endl;
            throw std::runtime_error("Unknown quantum gate.");
        }
        return res;
    }

    void concretizeInline() {
        this->isConcret = true;
        if (this->content.root->level != 1) {
            this->content = this->concretize();
        }
    }
    QuantumGateTerm cascade(const QuantumGateTerm& other) const {
        QuantumGateTerm res;
        if (this->zeroOperator || other.zeroOperator) {
            return QuantumGateTerm(false);
        } else if (this->ideOperator) {
            res = other;
        } else if (other.ideOperator) {
            res = *this;
        } else {
            throw std::runtime_error("Not supported gate cascade.");
        }
        return res;
    }
};

class SingleVecTerm : public QuantumTerm {
    public:
    // construct as a projector or state
    SingleVecTerm() {}
    SingleVecTerm(unsigned int qubits) {this->qNum = qubits;}
    SingleVecTerm(std::string s, unsigned int qubits) {
        // auto tmp = CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
        // s.size() is the realQubits, qNum is the total physical qubits.
        unsigned int level = ceil(log2(qubits));
        this->qNum = std::pow(2, level);
        CFLOBDD_COMPLEX_BIG stateVector = VectorComplexFloatBoost::MkBasisVector(level, s);
        stateVector = VectorComplexFloatBoost::VectorToMatrixInterleaved(stateVector);
        this->content = stateVector;
    }
    SingleVecTerm(CFLOBDD_COMPLEX_BIG x) {
        // Need copy?
        // type = false;
        this->content = x;
        this->qNum = std::pow(2, x.root->level-1);
    }
    bool getType() const override {return false;}
    std::unique_ptr<QuantumTerm> clone() const override {
        return std::make_unique<SingleVecTerm>(*this);
    }

    BIG_COMPLEX_FLOAT dot(const SingleVecTerm& other) const {
        /* bra(other) * ket(this) Note the complex number conjugation */
        // assert(this->type == false && other.type == false);
        unsigned int level = ceil(log2(this->qNum));
        auto tmpVec = Matrix1234ComplexFloatBoost::MatrixTranspose(other.content);
        tmpVec = Matrix1234ComplexFloatBoost::MatrixConjugate(tmpVec);
        auto tmp = Matrix1234ComplexFloatBoost::MatrixMultiplyV4(tmpVec, this->content);
        assert(tmp.root->rootConnection.returnMapHandle.Size() <= 2);
        auto resMap = tmp.root->rootConnection.returnMapHandle;
        BIG_COMPLEX_FLOAT amp;
        if(resMap.Size() == 2)
            amp = (resMap[0] != 0) ? resMap[0] : resMap[1];
        else
            amp = resMap[0];
        // amp = conj(amp);
        return amp;
    }
    CFLOBDD_COMPLEX_BIG normalize() const {
        // assert(this->type == false);
        auto H = ApplyGateF(std::pow(2, content.root->level-1), 0, Matrix1234ComplexFloatBoost::MkIdRelationInterleaved);
        CFLOBDD_COMPLEX_BIG c1 = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, content);
        CFLOBDD_COMPLEX_BIG c1_conj = Matrix1234ComplexFloatBoost::MatrixConjugate(c1);
        // VectorComplexFloatBoost::VectorPrintColumnHead(c1_conj, std::cout);
        c1_conj = Matrix1234ComplexFloatBoost::MatrixTranspose(c1_conj);
        auto mulres = Matrix1234ComplexFloatBoost::MatrixMultiplyV4(c1_conj, c1);
        auto resMap = mulres.root->rootConnection.returnMapHandle;
        // Maybe #BUGS here!
        double dimfactor = std::pow(double(2), double(std::pow(2, content.root->level-1)-1));
        // double dimfactor = 1;
        assert(resMap.Size() <= 2);
        BIG_COMPLEX_FLOAT amp;
        if(resMap.Size() == 2) {
            amp = (resMap[0] != 0) ? resMap[0] : resMap[1];
            // std::cout << "SingleVecTerm::normalize() amp = " << amp << std::endl;
            assert(abs(amp.imag()*dimfactor) < 1e-8 && amp.real() > 0);
            double factor = double(sqrt(amp.real()));
            // std::cout << "SingleVecTerm::normalize() factor = " << factor << std::endl;
            c1 = (1/factor) * c1;
        } else {
            std::cout << "Warning: SingleVecTerm::normalize() has only one factor!" << std::endl;
            // Here, assump the only factor is zero!
            amp = resMap[0];
            assert(abs(amp.imag()*dimfactor) < 1e-8 && abs(amp.real()*dimfactor) < 1e-8);
            c1 = VectorComplexFloatBoost::NoDistinctionNode(content.root->level, 0);
        }
        return c1;
    }
    void normalizeInline() {
        this->content = this->normalize();
    }
    CFLOBDD_COMPLEX_BIG projectOnto(const SingleVecTerm& other) const {
        /* project this onto other, \ket(other)\bra(other)\ket(this) == \bra(other)\ket(this)\ket(other) */
        // The result is not a normalized vector.
        // assert(this->type == false);
        // No need for other.dot(other)!!
        BIG_COMPLEX_FLOAT norm = this->dot(other) / other.dot(other);
        return norm * other.content;
    }
    bool isZero() const {
        return sqrt(this->dot(*this).real()) < 1e-8;
    }
    CFLOBDD_COMPLEX_BIG applyGate(const QuantumGateTerm& other, bool direction) const {
        /* If direction is false, apply an inverse gate. */
        /* This function is INPLACE! */
        CFLOBDD_COMPLEX_BIG operand = other.concretize();
        CFLOBDD_COMPLEX_BIG res ;
        if (direction) {
            // Forward induction
            // Temporarily use the content as the operand.

            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(operand, this->content);
        } else {
            // Backward induction
            operand = Matrix1234ComplexFloatBoost::MatrixConjugate(operand);
            operand = Matrix1234ComplexFloatBoost::MatrixTranspose(operand);
            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(operand, this->content);
        }
        return res;
    }
    void applyGateInline(const QuantumGateTerm& other, bool direction) {
        this->content = this->applyGate(other, direction);
    }
    double getMod() {
        auto amp = this->dot(*this);
        double res = double(sqrt(amp.real()));
        return res;
    }
};

class CliffordTerm : public QuantumTerm {
    // Represent a Clifford gate, which is a special type of quantum gate.
    bool getType() const override {return true;}
};

class TableauTerm : public QuantumTerm {
    // Represent a stabilizer subspace, in the special a stabilizer state.
    bool getType() const override {return false;}
};

/*
Quantum Opertion: The integration of two types of quantum items:
1. A quantum operator: list of {gate(name), unitary(CFLOBDD), vector(CFLOBDD)}
2. A projection operator: list of orthogonal vectors(CFLOBDD).
*/
class QOperation {
    public:
    /* type == true means the QOperation is a gate-type operation, the QuantumTerm in oplist are quantum gates or SignleVecTerm
    * type == false means the QOperation is a projective operation, the QuantumTerm in oplist are CFLOBDD items (SingleVecTerm)
    */
    bool type;
    std::vector<std::unique_ptr<QuantumTerm>> oplist;
    bool normalized = 0;
    unsigned int qNum = 0;
    unsigned int realqNum = 0;
    bool isIdentity = false;
    int isProj = -1;
    // std::unique_ptr<Node> ast = nullptr;
    public:
    /* The type must be specified */
    QOperation() : type(0) {}
    QOperation(std::vector<std::string> strings) {
        // Construct a QOperation of basic vectors.
        assert(strings.size() > 0);
        for (const auto& str : strings) {
            assert(str.size() > 0);
            assert(str.size() == strings[0].size());
        }
        this->type = false;
        this->realqNum = strings[0].size();
        this->qNum = std::pow(2, ceil(log2(strings[0].size())));
        for (const auto& str : strings) {
            SingleVecTerm term(stringPadding(str, this->qNum), std::pow(2, ceil(log2(str.size()))));
            oplist.push_back(std::make_unique<SingleVecTerm>(term));
        }
        this->normalized = true;
    }
    QOperation(std::string nam, unsigned int qNum, std::vector<unsigned int> idx, std::vector<double> pars=std::vector<double>{}) {
        // Construct a QOperation of quantum gate.
        // Could merge
        this->type = true;
        this->realqNum = qNum;
        unsigned int logicqNum = std::pow(2, ceil(log2(qNum)));
        this->qNum = logicqNum;
        // Just a copy of CreateProjectiveMeasQO
        if (nam == "meas0") {
            this->isProj = idx[0];
            QuantumGateTerm tmp("meas0", std::vector<unsigned int>{idx[0]}, std::vector<double>{}, logicqNum);
            // Whether it is necessary to concretize?
            tmp.concretizeInline();
            // Attention!
            this->oplist.push_back(tmp.clone());
        } else if (nam == "meas1") {
            this->isProj = idx[0];
            QuantumGateTerm tmp("meas1", std::vector<unsigned int>{idx[0]}, std::vector<double>{}, logicqNum);
            tmp.concretizeInline();
            // Attention!
            this->oplist.push_back(tmp.clone());
        } else if (nam == "reset") {
            // The reset operator creates a mixed state
            this->isProj = idx[0];
            QuantumGateTerm cond0("meas0", std::vector<unsigned int>{idx[0]}, std::vector<double>{}, logicqNum);
            // Whether it is necessary to concretize?
            cond0.concretizeInline();
            // Attention!
            this->oplist.push_back(cond0.clone());
            QuantumGateTerm cond1("reset0", std::vector<unsigned int>{idx[0]}, std::vector<double>{}, logicqNum);
            cond1.concretizeInline();
            this->oplist.push_back(cond1.clone());
        }
        else {
            oplist.push_back(std::make_unique<QuantumGateTerm>(nam, idx, pars, logicqNum));
            // Here, isIdentity is used in type == true case.
            if (toLower(nam) == "i") {
                isIdentity = true;
            } else {
                isIdentity = false;
            }
        }
    }
    QOperation(const QOperation& other) : type(other.type), normalized(other.normalized), qNum(other.qNum), isIdentity(other.isIdentity), isProj(other.isProj) {
        // 深拷贝 oplist
        oplist.reserve(other.oplist.size()); // 预分配空间以提高效率
        for (const auto& term : other.oplist) {
            if (term) { // 检查指针是否为空
                oplist.push_back(term->clone()); // 使用 clone 方法进行深拷贝
            } else {
                oplist.push_back(nullptr);
            }
        }
        // this->ast = std::make_unique<Node>(*other.ast);
    }
    QOperation(const QOperation& other1, const QOperation& other2) {
        assert(other1.type == other2.type);
        if (other1.qNum != other2.qNum) {
            std::cout << other1.qNum << " != " << other2.qNum << std::endl;
        }
        assert(other1.qNum == other2.qNum);
        this->qNum = other1.qNum;
        this->type = other1.type;
        this->normalized = false;
        oplist.reserve(other1.oplist.size() + other2.oplist.size());
        for (const auto& term : other1.oplist) {
            if (term) { // 检查指针是否为空
                oplist.push_back(term->clone()); // 使用 clone 方法进行深拷贝
            } else {
                oplist.push_back(nullptr);
            }
        }
        for (const auto& term : other2.oplist) {
            if (term) { // 检查指针是否为空
                oplist.push_back(term->clone()); // 使用 clone 方法进行深拷贝
            } else {
                oplist.push_back(nullptr);
            }
        }
    }
    QOperation(QOperation&& other) noexcept : type(other.type), oplist(std::move(other.oplist)), normalized(other.normalized), qNum(other.qNum), isIdentity(other.isIdentity), isProj(other.isProj) {
        // Move constructor
        // this->ast = std::make_unique<Node>(*other.ast);
        // this->ast = std::move(other.ast);
    }
    QOperation(bool t) : type(t) {}
    QOperation(bool t, std::vector<std::unique_ptr<QuantumTerm>>&& c) : type(t), oplist(std::move(c)), normalized(false) {
        // Move constructor
        if (!c.empty() && c[0]) {
            qNum = c[0]->qNum;
        } else {
            qNum = 0;
        }
    }
    QOperation(bool t, std::vector<std::unique_ptr<QuantumTerm>>&& c, bool n) : type(t), oplist(std::move(c)), normalized(n) {
        // Move constructor with normalization
        if (!c.empty() && c[0]) {
            qNum = c[0]->qNum;
        } else {
            qNum = 0;
        }
    }

    QOperation& operator=(const QOperation& other) {
        if (this != &other) {
            type = other.type;
            normalized = other.normalized;
            isIdentity = other.isIdentity;
            isProj = other.isProj;
            qNum = other.qNum;
            oplist.clear();
            oplist.reserve(other.oplist.size());
            for (const auto& term : other.oplist) {
                if (term) {
                    oplist.push_back(term->clone());
                } else {
                    oplist.push_back(nullptr);
                }
            }
        }
        return *this;
    }
    
    std::string getName() const {
        if (this->type == false) {
            return "Projective Operation";
        } else {
            std::string res;
            for (const auto& term : this->oplist) {
                if (term && term->getType()) {
                    auto* gateTerm = dynamic_cast<QuantumGateTerm*>(term.get());
                    if (gateTerm) {
                        res += gateTerm->name + " ";
                    } else {
                        std::cout << "Unknown quantum term type in getName()." << std::endl;
                    }
                } else {
                    std::cout << "Null quantum term in getName()." << std::endl;
                }
            }
            return res.empty() ? "Empty Quantum Operation" : res;
        }
    }

    void append(std::unique_ptr<QuantumTerm> qt) {
        if (qt->getType() != this->type) {
            std::cout << "Append a wrong type of quantum term.";
            return;
        }
        oplist.push_back(std::move(qt));
    }
    std::vector<std::unique_ptr<QuantumTerm>> fetch(unsigned int begin,unsigned int end) {
        // fetch the oplist from begin to end, using move semantics
        assert(end < static_cast<int>(this->oplist.size()));
        std::vector<std::unique_ptr<QuantumTerm>> res;
        res.reserve(end - begin + 1);
        for (size_t i = begin; i <= end; i++) {
            if (this->oplist[i] == nullptr) {
                res.push_back(nullptr);
            } else {
                res.push_back(this->oplist[i]->clone()); // Use clone to ensure deep copy
            }
        }
        return res;
    }

    int findVectorContent(const SingleVecTerm& vec) const {
        // Check if the vector is in the oplist
        assert(this->type == false);
        unsigned int i = 0;
        for (const auto& term : this->oplist) {
            if (term && term->getType() == false) {
                auto* ivec = dynamic_cast<SingleVecTerm*>(term.get());
                if (ivec && ivec->content == vec.content) {
                    return i;
                }
            }
            i++;
        }
        return -1;
    }

    void genProjMeasSpace() {
        assert(this->isProj >= 0);
        if (this->oplist.size() > 1) {
            std::cout << "Already generated measure support vectors." << std::endl;
            return;
        }
        assert(this->oplist.size() == 1);
        // Check the equivalence of std::string
        bool measureType = (dynamic_cast<QuantumGateTerm*>(this->oplist[0].get())->name == "meas1");
        int dim_num = std::pow(2, this->qNum-1);
        // Clear the oplist first
        this->oplist.clear();
        this->oplist.reserve(dim_num); // Reserve space for the basis vectors
        // This is a very time-consuming operation!!!
        for (int i = 0; i <= dim_num - 1; i++) {
            // Generate the ith bit string of the basis. The value of index of isProj is measureType, iterator all other indexes except isProj.
            std::string basis_str;
            for (int j = 0; j < this->qNum; j++) {
                if (j == this->isProj) {
                    basis_str += measureType ? "1" : "0"; // Measure type is 1, otherwise is 0.
                } else if (j < this->isProj) {
                    basis_str += (i & (1 << j)) ? "1" : "0"; // The jth bit of i.
                } else if (j > this->isProj) {
                    basis_str += (i & (1 << (j - 1))) ? "1" : "0"; // The j-1th bit of i.
                }
            }
            // Create a SingleVecTerm with the basis string and qNum
            SingleVecTerm basis(basis_str, std::pow(2, ceil(log2(this->qNum))));
            // Append the basis to the oplist
            this->oplist.push_back(std::make_unique<SingleVecTerm>(basis));
        }
        this->type = false; // Change the type to support vectors of subspaces.
        this->isProj = -1; // Reset isProj to indicate this is a projective operation.
        this->normalized = true; // Mark as normalized since we have generated the basis vectors.
    }

    SingleVecTerm projectIn(const SingleVecTerm& vec) const {
        // project the vector onto the QOperation
        if (this->oplist[0]->getType() == false) {
            assert(vec.getType() == false);
            assert(this->oplist[0]->qNum == vec.qNum);
            unsigned int level = ceil(log2(vec.qNum));
            CFLOBDD_COMPLEX_BIG content = VectorComplexFloatBoost::NoDistinctionNode(level+1, 0); // Check: Initialization
            for (size_t i = 0; i < this->oplist.size(); i++) {
                auto* ivec = dynamic_cast<SingleVecTerm*>(this->oplist[i].get());
                if (!ivec) continue;
                content = content + vec.projectOnto(*ivec);
            }
            // std::cout << "SingleVecTerm projectIn: content is zero? " << checkifzero(content) << std::endl;
            SingleVecTerm res(content);
            res.qNum = this->qNum;
            return res;
        } else {
            // this is a projective operation, this->oplist[0] is a projective operator.
            assert(vec.getType() == false);
            assert(this->oplist[0]->qNum == vec.qNum);
            // TODO: Do the projection.
            CFLOBDD_COMPLEX_BIG content; // Check: Initialization
            // Direction: true.
            content = vec.applyGate(*dynamic_cast<QuantumGateTerm*>(this->oplist[0].get()), true);
            SingleVecTerm res(content);
            res.qNum = this->qNum;
            return res;
        }
    }

    void GramSchmidt(int begin, int end, std::tuple<unsigned int, unsigned int> existedOrthogonal = std::make_tuple(0, 0)) {
        // Here need to optimize the orthogonalBasis
        // TODO: If the size is 1
        assert(this->type == false);
        assert(end < static_cast<int>(this->oplist.size()));
        if (this->oplist.size() == 0) {
            return;
        }
        std::vector<std::unique_ptr<QuantumTerm>> orthogonalBasis;
        for (size_t i = begin; i <= end; i++) {
            auto* ivec = dynamic_cast<SingleVecTerm*>(this->oplist[i].get());
            if (!ivec) {std::cout << "strange i nullptr";continue;}
            for (const auto& basis : orthogonalBasis) {
                auto* jvec = dynamic_cast<SingleVecTerm*>(basis.get());
                if (!jvec) {std::cout << "strange j nullptr";continue;}
                ivec->content = ivec->content + (-1) * ivec->projectOnto(*jvec);
            }
            // if (!ivec->isZero()) {
            if (!checkifzero(ivec->content)) {
                ivec->normalizeInline();
                // the orthogonal basis must be normalized! make sure the other.dot(other) == 1 in the project onto function.
                orthogonalBasis.push_back(std::move(this->oplist[i]));
            }
        }
        for (size_t i = 0; i < orthogonalBasis.size(); i++) {
            this->oplist[begin+i] = std::move(orthogonalBasis[i]);
        }
        this->oplist.erase(this->oplist.begin() + begin + orthogonalBasis.size(), this->oplist.begin() + end + 1);
        // this->oplist = std::move(orthogonalBasis);
        this->normalized = true;
    }

    QOperation minus(const QOperation& other) const {
        // Return a set of basis vectors that are in this but not in other.
        if (other.isProj < 0) {
            int dimOther = -1;
            std::optional<QOperation> totalDisj;
            if (other.normalized) {
                dimOther = other.oplist.size();
                totalDisj.emplace(other,*this);
            } else {
                std::cout << "Warning: An operator is not normalized, will normalize it." << std::endl;
                QOperation otherCopy(other);
                otherCopy.GramSchmidt(0, static_cast<int>(otherCopy.oplist.size())-1);
                dimOther = other.oplist.size();
                totalDisj.emplace(otherCopy,*this);
            }
            // Optimize: add some already-orthogonal basis.
            totalDisj->GramSchmidt(0, static_cast<int>(totalDisj->oplist.size())-1);
            std::vector<std::unique_ptr<QuantumTerm>> orthogonalBasis = totalDisj->fetch(dimOther, static_cast<int>(totalDisj->oplist.size())-1);
            QOperation res(false);
            for (int i = 0; i < orthogonalBasis.size(); i++) {
                SingleVecTerm projVec = this->projectIn(*dynamic_cast<SingleVecTerm*>(orthogonalBasis[i].get()));
                // if (!projVec.isZero()) {
                if (!checkifzero(projVec.content)) {
                    res.append(projVec.clone()); // Use clone to ensure deep copy
                }
            }
            res.GramSchmidt(0, static_cast<int>(res.oplist.size())-1);
            res.qNum = this->qNum;
            return res;
        }
        assert(other.isProj >= 0);
        assert(this->isProj < 0);
        assert(this->type == false);
        assert(this->oplist[0]->getType() == false);
        // not (not A or not B) == A - (A - B) == B - (B - A)? Yes, but A - B cannot be computed by GramSchmidt directly.
        QOperation temp = *this; // Make a copy of this QOperation
        std::vector<std::unique_ptr<QuantumTerm>> orthogonalBasis;
        // for loop: find a set of orthogonal basis that is orthogonal to other (A binary projection)
        for (size_t i = 0; i < temp.oplist.size(); i++) {
            auto* ivec = dynamic_cast<SingleVecTerm*>(temp.oplist[i].get());
            if (!ivec) continue;
            for (const auto& basis : orthogonalBasis) {
                auto* jvec = dynamic_cast<SingleVecTerm*>(basis.get());
                if (!jvec) continue;
                ivec->content = ivec->content + (-1) * ivec->projectOnto(*jvec);
            }
            // TODO: Check the projectIn function
            ivec->content = ivec->content + (-1) * other.projectIn(*ivec).content;
            // if (!ivec->isZero()) {
            if (!checkifzero(ivec->content)) {
                ivec->normalizeInline();
                orthogonalBasis.push_back(std::move(temp.oplist[i]));
            }
        }
        // Here, it seems no need to copy the orthogonalBasis to temp, just use the orthogonalBasis directly.
        for (size_t i = 0; i < orthogonalBasis.size(); i++) {
            temp.oplist[i] = std::move(orthogonalBasis[i]);
        }
        temp.oplist.erase(temp.oplist.begin() + orthogonalBasis.size(), temp.oplist.end());
        // this->oplist = std::move(orthogonalBasis);
        QOperation res;
        for (size_t i = 0; i < temp.oplist.size(); i++) {
            SingleVecTerm projVec = this->projectIn(*dynamic_cast<SingleVecTerm*>(temp.oplist[i].get()));
            // if (!projVec.isZero()) {
            if (!checkifzero(projVec.content)) {
                res.append(projVec.clone()); // Use clone to ensure deep copy
            }
        }
        res.GramSchmidt(0, static_cast<int>(res.oplist.size())-1);
        res.qNum = this->qNum;
        return res;
    }

    QOperation negation() const {
        // compute the negation of some of operators. For example, a projective operator, the negation is the set of all other orthogonal vectors.
        // For now only projective operators are supported.
        assert(this->type == true);
        assert(this->isProj >= 0);
        
        std::string name = dynamic_cast<QuantumGateTerm*>(this->oplist[0].get())->name;
        unsigned int idx = dynamic_cast<QuantumGateTerm*>(this->oplist[0].get())->index[0];
        assert(idx == this->isProj);
        if (name == "meas0") {
            QOperation res(std::string("meas1"), this->qNum, std::vector<unsigned int>{idx});
            return res;
        } else if (name == "meas1") {
            QOperation res(std::string("meas0"), this->qNum, std::vector<unsigned int>{idx});
            return res;
        } else {
            throw std::runtime_error("Unknown negation operation.");
        }
        return *this;
    }
    QOperation conjunction(const QOperation& other) const {
        /*** To do the conjunction:
         * 1. preserve the abstract semantic tree;
         * 2. widening function;
         * 3. hard calculate \neg(\neg A \lor \neg B).
         ***/ 
        // Make Sure! At least one operand is a projective operator!
        assert(this->isProj < 0);
        assert(!this->type);
        if (other.isIdentity) {
            return *this;
        }
        if (this->isIdentity) {
            return other;
        }
        // TODO: What if this->oplist.size() == 0?
        if (other.oplist.size() == 0) {
            return QOperation();
        }
        assert(!other.oplist[0]->getType()); // Is this right?
        if (!this->normalized) {
            throw std::runtime_error("The QOperation is not normalized.");
            // this->GramSchmidt(0, this->oplist.size()-1);
            // this->normalized = true;
        }
        assert(other.normalized);
        // if (!other.normalized) {
        //     other.GramSchmidt(0, other.oplist.size()-1);
        //     other.normalized = true;
        // }
        QOperation op1(*this, other);
        // Modify: optimize the GramSchmidt length incremental
        op1.GramSchmidt(0, static_cast<int>(op1.oplist.size())-1);
        std::vector<std::unique_ptr<QuantumTerm>> negthis = op1.fetch(this->oplist.size(), static_cast<int>(op1.oplist.size())-1);
        QOperation op2(other, *this);
        op2.GramSchmidt(0, static_cast<int>(op2.oplist.size())-1);
        std::vector<std::unique_ptr<QuantumTerm>> negother = op2.fetch(other.oplist.size(), static_cast<int>(op2.oplist.size())-1);
        std::vector<std::unique_ptr<QuantumTerm>> remaining = op1.fetch(0, static_cast<int>(this->oplist.size())-1);
        // create a new QOperation with oplist negthis and negother
        QOperation op3(false, std::move(negthis));
        op3.oplist.insert(op3.oplist.end(), std::make_move_iterator(negother.begin()), std::make_move_iterator(negother.end()));
        op3.GramSchmidt(0, static_cast<int>(op3.oplist.size())-1);
        int dimNeg1OrNeg2 = op3.oplist.size();
        op3.oplist.insert(op3.oplist.end(), std::make_move_iterator(remaining.begin()), std::make_move_iterator(remaining.end()));
        op3.GramSchmidt(0, static_cast<int>(op3.oplist.size())-1);
        std::vector<std::unique_ptr<QuantumTerm>> resvec = op3.fetch(dimNeg1OrNeg2, static_cast<int>(op3.oplist.size())-1);
        QOperation res(false, std::move(resvec));
        res.normalized = true;
        res.qNum = this->qNum;
        return res;
    }

    QOperation conjunction_simp(const QOperation& other) const {
        /*
        * In two cases this function is invoked:
        * 1. when computing the preImage of a projective operator, using to compute the Sasaki-hook;
        * 2. when merging different pre-conditions.
        */
        // We must make sure at least one operand is not a projective operator!
        // Already Schmidted.
        // this: a subspace, other: a subspace (type == false) or a projective measurement (type == true).
        assert(this->isProj < 0);
        if (this->type || !this->normalized) {
            std::cout << this->type << " " << this->normalized << std::endl;
        }
        assert(!this->type && this->normalized);
        if (other.isIdentity) {
            return *this;
        }
        if (this->isIdentity && other.isProj < 0) {
            return other;
        }
        if (this->isIdentity && other.isProj >= 0) {
            QOperation othercopy = other;
            othercopy.genProjMeasSpace();
            return othercopy;
        }
        if (other.oplist.size() == 0 || this->oplist.size() == 0) {
            // If one of the operands is empty, return an empty QOperation
            QOperation res(false);
            res.qNum = this->qNum;
            res.normalized = true;
            return res;
        }
        if (other.oplist[0]->getType() == true) {
            // Case 1: other.oplist[0] is a quantumGate type: a projective measurement
            // std::cout << "Conjunction of a supp subspace and a projective operator." << std::endl;
            // TODO: This is the most time-consuming part!
            assert(other.isProj >= 0);
            assert(other.oplist.size() == 1);
            // auto *jgate = dynamic_cast<QuantumGateTerm*>(other.oplist[0].get());
            // TODO: Specialize for the case when dimension 1
            if (this->oplist.size() == 1) {}
            QOperation minusQO = this->minus(other);
            int dimThisMinusOther = minusQO.oplist.size();
            minusQO = QOperation(minusQO, *this);
            minusQO.GramSchmidt(0, static_cast<int>(minusQO.oplist.size())-1);
            std::vector<std::unique_ptr<QuantumTerm>> negthis = minusQO.fetch(dimThisMinusOther, static_cast<int>(minusQO.oplist.size())-1);
            QOperation res(false, std::move(negthis));
            res.normalized = true;
            res.qNum = this->qNum;
            return res;
        } else {
            // TODO
            if (this->oplist.size() == 1) {}
            if (other.oplist.size() == 1) {}
            // A normal conjunction between two sets of subspaces.
            assert(other.normalized);
            QOperation thisMinusOther = this->minus(other);
            int dimThisMinusOther = thisMinusOther.oplist.size();
            QOperation op1(thisMinusOther, *this);
            // Modify: optimize the GramSchmidt length incremental
            op1.GramSchmidt(0, static_cast<int>(op1.oplist.size())-1);
            std::vector<std::unique_ptr<QuantumTerm>> resvec = op1.fetch(dimThisMinusOther, static_cast<int>(op1.oplist.size())-1);
            QOperation res(false, std::move(resvec));
            res.normalized = true;
            res.qNum = this->qNum;
            return res;
        }
    }

    QOperation disjunction(const QOperation& other) const {
        /*** To do the disjunction:
         * 1. preserve the abstract semantic tree;
         * 2. widening function;
         * 3. hard calculate \neg(\neg A \land \neg B).
         * Already Schmidted.
         ***/
        assert(!this->type && !other.type);
        if (other.oplist.size() == 0) {
            // The vector remains unnormalized here!
            return *this;
        }
        if (this->oplist.size() == 0) {
            // The vector remains unnormalized here!
            return other;
        }
        if (this->isIdentity || other.isIdentity || this->oplist.size() == std::pow(2, this->qNum) || other.oplist.size() == std::pow(2, other.qNum)) {
            QOperation res(false);
            res.normalized = true;
            res.qNum = this->qNum;
            res.isIdentity = true;
            return res;
        }
        assert(!other.oplist[0]->getType());
        QOperation res(*this, other);
        // std::cout << "Disjunction: " << res.oplist.size() << std::endl;
        res.GramSchmidt(0, static_cast<int>(res.oplist.size())-1);
        // std::cout << "Disjunction after GramSchmidt: " << res.oplist.size() << std::endl;
        return res;
    }

    QOperation preImage(const QOperation& other) {
        /* The pre-image of a quantum operator */
        assert(this->type == false && other.type == true);
        QOperation res(false);
        /* Two cases here: other.oplist contains only QuantumGateTerm, or only SingleVecTerm. 
                In the later case, we should use jvec->projectOnto() */
        if (this->oplist.size() == 0) {
            res.normalized = true;
            res.qNum = this->qNum;
            return res; // If this is an empty operator, return an empty operator.
        }
        if (other.oplist[0]->getType() == false) {
            // Case 0.5: other.oplist[0] is a SingleVecTerm, other is a set of orthogonal vectors.
            // std::cout << "Case 0.5: other.oplist[0] is a SingleVecTerm, other is a set of orthogonal vectors." << std::endl;
            res = this->conjunction_simp(other); // This is not a const operator!
        } else if (other.oplist[0]->getType() == true && other.isProj >= 0) {
            // Case 1: other is a binary projection operator
            // other.genProjMeasSpace();
            // std::cout << "Case 1: other is a binary projection operator." << std::endl;
            res = this->conjunction_simp(other);
            // TODO: Should negOther be saved? Save: reduce the time complexity; Don't save: reduce the space complexity.
            // I need to fix the following codes to reduce the time complexity!!!
            QOperation negOther = other.negation();
            negOther.genProjMeasSpace();
            res = res.disjunction(negOther);
        } else {
            // Case 2: other.oplist[0] is a QuantumGateTerm
            // std::cout << "Case 2: other.oplist[0] is a QuantumGateTerm. " << other.oplist[0]->getType() << " " << other.isProj << std::endl;
            for (size_t i = 0; i < this->oplist.size(); i++) {
                auto* ivec = dynamic_cast<SingleVecTerm*>(this->oplist[i].get());
                if (!ivec) {std::cout << "Strange nullptr" << std::endl; continue;}
                
                // if (other.oplist.size() == 0) continue;
                for (size_t j = 0; j < other.oplist.size(); j++) {
                    auto* jgate = dynamic_cast<QuantumGateTerm*>(other.oplist[j].get());
                    if (!jgate) continue;
                    assert(jgate->getType() == true);
                    // Backward induction
                    auto tmp = ivec->applyGate(*jgate, false);
                    res.append(std::make_unique<SingleVecTerm>(SingleVecTerm(tmp)));
                    /* We have GramSchmidt in preImage, no need to handle the res uniqueness */
                    // SingleVecTerm preVec(tmp);
                    // if (res.oplist.size() == 0) {
                    //     res.append(std::make_unique<SingleVecTerm>(preVec));
                    // } else if (res.findVectorContent(preVec) == -1) {
                    //     // If the preVec is not in the res, append it.
                    //     res.append(std::make_unique<SingleVecTerm>(preVec));
                    // } else if (res.findVectorContent(preVec) >= 0) {
                    //     // If the preVec is in the res, add the amplitude to the existing vector.
                    // }
                }
            }
        }
        // TODO: Here need a judgement: if all operators are unitary, and the space before preimage is orthogonal, then no need to GramSchmidt.
        // if ((other.oplist[0]->getType() == true && other.isProj < 0 && this->normalized)) {
        //     // If the res.oplist.size() <= 1, no need to GramSchmidt.
        //     // If other is a projective operator, no need to GramSchmidt.
        //     // Normalize the results.
        //     res.normalized = true;
        //     res.qNum = this->qNum;
        //     return res;
        // }
        res.GramSchmidt(0, static_cast<int>(res.oplist.size())-1);
        res.qNum = this->qNum;
        return res;
    }
    
    QOperation postImage(const QOperation& other) const {
        /* The post-image of a quantum operator */
        assert(this->type == false && other.type == true);
        QOperation res;
        for (size_t i = 0; i < this->oplist.size(); i++) {
            auto* ivec = dynamic_cast<SingleVecTerm*>(this->oplist[i].get());
            if (!ivec) continue;
            for (size_t j = 0; j < other.oplist.size(); j++) {
                auto* jgate = dynamic_cast<QuantumGateTerm*>(other.oplist[j].get());
                if (!jgate) continue;
                auto tmp = ivec->applyGate(*jgate, true);
                res.append(std::make_unique<SingleVecTerm>(SingleVecTerm(tmp)));
                /* We don't have GramSchmidt in postImage (for some interface of probability), we need to handle the uniqueness */
                // SingleVecTerm postVec(tmp);
                // if (res.oplist.size() == 0) {
                //     res.append(std::make_unique<SingleVecTerm>(postVec));
                // } else if (res.findVectorContent(postVec) == -1) {
                //     // If the postVec is not in the res, append it.
                //     res.append(std::make_unique<SingleVecTerm>(postVec));
                // } else if (res.findVectorContent(postVec) >= 0) {
                //     // If the postVec is in the res, add the amplitude to the existing vector.
                // }
            }
        }
        res.normalized = false;
        res.qNum = this->qNum;
        // Use GramSchmidt to handle the uniqueness of the vectors.
        if (res.oplist.size() > 1) {
            res.GramSchmidt(0, static_cast<int>(res.oplist.size())-1);
            res.normalized = true;
        }
        return res;
    }

    // QOperation operator+(const QOperation& other) const {
    //     // Still unclear: when is the plus used?
    //     if (type == true) {
    //         if (other.type == true) {
    //             /* gate plus gate */
    //             /* append gates to the end */
    //             if (this->oplist.size() > 0 && other.oplist.size() > 0) {
    //                 assert(this->oplist[0]->getType() == other.oplist[0]->getType());
    //             }
    //             std::vector<std::unique_ptr<QuantumTerm>> newlist = this->oplist;
    //             newlist.insert(newlist.end(), other.oplist.begin(), other.oplist.end());
    //             QOperation res(true, newlist);
    //             return res;
    //         } else {
    //             /* gate plus cflobdd */
    //             throw std::runtime_error("Cannot support plus between gates and operators");
    //         }
    //     } else {
    //         if (other.type == true) {
    //             /* cflobdd plus cflobdd */
    //             /* disjunction */
    //             std::vector<std::unique_ptr<QuantumTerm>> newlist = this->oplist;
    //             newlist.insert(newlist.end(), other.oplist.begin(), other.oplist.end());
    //             QOperation res(true, newlist, false);
    //             return res;
    //         } else {
    //             /* cflobdd plus gate */
    //             throw std::runtime_error("Cannot support plus between gates and operators");
    //         }
    //     }
    // }
    // QOperation operator*(const QOperation& other) const {
    //     if (type == true) {
    //         if (other.type == true) {
    //             /* gate times gate */
    //             /* cascade */
    //             // throw std::runtime_error("Logically, it should be supported");
    //             for (size_t i = 0; i < this->oplist.size(); i++) {
    //                 auto* gate1 = dynamic_cast<QuantumGateTerm*>(this->oplist[i].get());
    //                 for (size_t j = 0; i < other.oplist.size(); j++) {
    //                     CFLOBDD_COMPLEX_BIG newItem;
    //                     auto* gate2 = dynamic_cast<QuantumGateTerm*>(other.oplist[j].get());
    //                     if (this->oplist[0]->getType() == true) {
    //                         // Assert all type in a operation are consistant
    //                         // all operators are gates
    //                         auto res = gate1->cascade(*gate2);
    //                         if (res.isConcret) {
    //                             newItem = res.concretize();
    //                         } else {
    //                             newItem = res.content;
    //                         }
    //                     } else {
    //                         throw std::runtime_error("Cannot support the cascade of two operators");
    //                     }
    //                 }
                    
    //             }
                
    //         } else {
    //             /* gate times cflobdd */
    //             QOperation res;
    //             for (size_t i = 0; i < other.oplist.size(); i++) {
    //                 CFLOBDD_COMPLEX_BIG newItem;
    //                 for (size_t j = 0; j < this->oplist.size(); j++) {
    //                     if (this->oplist[0]->getType() == true) {
    //                         // Assert all type in a operation are consistant
    //                         // all operators are gates
    //                         auto* gate1 = dynamic_cast<QuantumGateTerm*>(this->oplist[j].get());
    //                     } else {
    //                         // all operators are supporting vector of a space
    //                     }
                        
    //                 }
    //             }
                
    //         }
    //     } else {
    //         if (other.type == true) {
    //             /* cflobdd times cflobdd */
    //             /* conjunction */
    //         } else {
    //             /* cflobdd times gate*/
    //         }
    //     }
    // }
    
    
    int compare(const QOperation& other) const {
        /* Compare the relation of this and other */
        /* 0: this is included in other; 1: other is included in this; 2: exclude; 3: intersect but not include; 4: equality*/
        assert(this->normalized && other.normalized);
        if (this->oplist.size() > 0) {
            assert(this->oplist[0]->getType() == false);
        }
        if (other.oplist.size() > 0) {
            assert(other.oplist[0]->getType() == false);
        }
        if (this->oplist.size() == 0 && other.oplist.size() == 0) {
            return -1;
        }
        int inside_cnt = 0;
        for (size_t i = 0; i < this->oplist.size(); i++) {
            auto* ivec = dynamic_cast<SingleVecTerm*>(this->oplist[i].get());
            if (!ivec) continue;
            // Check if every single vector in this can be totally projected onto other
            SingleVecTerm tmp = other.projectIn(*ivec);
            auto tmpdot = tmp.dot(tmp); // TODO: A more efficient way to justify the difference?
            if (abs(tmpdot.real()-1) < 1e-8 && abs(tmpdot.imag()) < 1e-8) {
                // this is a subspace of other
                if (i == static_cast<int>(this->oplist.size())-1) {
                    inside_cnt++;
                }
            }
        }
        if (inside_cnt == this->oplist.size()) {
            if (this->oplist.size() == other.oplist.size()) {
                return 4;
            } else {
                return 0;
            }
        }
        if (inside_cnt == 0) {
            return 2;
        }
        if (inside_cnt == other.oplist.size()) {
            return 1;
        }
        return 3;
    }
    bool operator==(const QOperation& other) const {
        /* If the support space of this is the subspace of other */
        assert(this->normalized && other.normalized);
        if (this->oplist.size() > 0) {
            assert(this->oplist[0]->getType() == false);
        }
        if (other.oplist.size() > 0) {
            assert(other.oplist[0]->getType() == false);
        }
        if (this->oplist.size() != other.oplist.size()) {
            return false;
        }
        if (this->compare(other) == 4) {
            return true;
        }
        return false;
    }
    bool operator<=(const QOperation& other) const {
        /* If the support space of this is the subspace of other */
        return this->compare(other) == 0;
    }
    
    // bool satisfy(const QOperation& other) const {
    //     /* If the support space of this is the subspace of other */
    //     return false;
    // }
    void print() const {
        std::cout << "Support vectors:" << std::endl;
        if (this->isIdentity) {
            std::cout << "Identity operator." << std::endl;
            return;
        }
        for (size_t i = 0; i < this->oplist.size(); i++) {
            VectorComplexFloatBoost::VectorPrintColumnHead(this->oplist[i]->content, std::cout);
            std::cout << std::endl;
        }
    }

    void printFormal() const {
        std::cout << "Support vectors:" << std::endl;
        if (this->isIdentity) {
            std::cout << "Identity operator." << std::endl;
            return;
        }
        for (size_t i = 0; i < this->oplist.size(); i++) {
            std::ostringstream oss;
            VectorComplexFloatBoost::VectorPrintColumnHead(this->oplist[i]->content, oss);
            std::string vecStr = oss.str();
            // Get a list of complex numbers from the string, the numbers are separated by spaces.
            std::istringstream iss(vecStr);
            std::vector<std::complex<double>> vec;
            std::string token;
            std::complex<double> cplxNum;
            while (iss >> token) {
                if (token.front() == '(' && token.back() == ')') {
                    double real, imag;
                    sscanf(token.c_str(), "(%lf,%lf)", &real, &imag);
                    cplxNum = std::complex<double>(real, imag);
                } else {
                    cplxNum = std::complex<double>(std::stod(token), 0.0);
                }
                vec.push_back(cplxNum);
            }
            std::ostringstream result;
            bool first = true;
            for (int k = 0; k < vec.size(); k++) {
                if (std::abs(vec[k]) > 1e-8) {
                    if (!first) {
                        result << " + ";
                    }
                    first = false;
                    const auto& amp = vec[k];
                    std::string idxStr;
                    for (int j = this->qNum - 1; j >= 0; --j) {
                        idxStr += ((k >> j) & 1) ? '1' : '0';
                    }
                    result << "(" << amp.real() << "," << amp.imag() << ")";
                    result << "|" << idxStr << ">";
                }
            }
            std::cout << result.str() << std::endl;
        }
        
    }
};

/* Const QOperation:
I: the identity operator as the top of the QOperation lattice
  --> properties: for all QOperation A, conjunction(A, I) = A
  --> properties: for all QOperation A, disjunction(A, I) = I
*/

QOperation CreateIdentityQO(unsigned int qNum) {
    unsigned int logicqNum = std::pow(2, ceil(log2(qNum)));
    QOperation res(false);
    res.realqNum = qNum;
    res.qNum = logicqNum;
    res.isIdentity = true;
    res.normalized = true;
    return res;
}

QOperation CreateZeroQO(unsigned int qNum) {
    unsigned int logicqNum = std::pow(2, ceil(log2(qNum)));
    QOperation res(false);
    res.realqNum = qNum;
    res.qNum = logicqNum;
    res.normalized = true;
    return res;
}

// Should be mearged into the constructor of QOperation
QOperation CreateProjectiveMeasQO(unsigned int qNum, unsigned int i, bool val) {
    unsigned int logicqNum = std::pow(2, ceil(log2(qNum)));
    QOperation res(true);
    res.qNum = logicqNum;
    res.realqNum = qNum;
    res.isProj = i;
    if (!val) {
        QuantumGateTerm tmp("meas0", std::vector<unsigned int>{i}, std::vector<double>{}, qNum);
        tmp.concretizeInline();
        // Attention!
        res.oplist.push_back(tmp.clone());
        // res.oplist.push_back(std::make_unique<QuantumGateTerm>(tmp));
    } else {
        QuantumGateTerm tmp("meas1", std::vector<unsigned int>{i}, std::vector<double>{}, qNum);
        tmp.concretizeInline();
        res.oplist.push_back(tmp.clone());
        // res.oplist.push_back(std::make_unique<QuantumGateTerm>(tmp));
    }
    return res;
}

// QOperation CreatePartialProjectiveQO(std::vector<std::string> parBasis, std::vector<unsigned int> parIndex, unsigned int qNum) {
//     unsigned int logicqNum = std::pow(2, ceil(log2(qNum)));
//     QOperation res(true);
//     res.qNum = logicqNum;
//     res.realqNum = qNum;
//     res.isProj = parIndex[0]; // Assume the first index is the projective index


// }

#endif
