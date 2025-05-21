#ifndef _QUANTUM_OPERATION
#define _QUANTUM_OPERATION

#include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include "cflobdd/CFLOBDD/vector_complex_float_boost.h"
#include "AST.hpp"
#include <random>
#include <queue>
#include <vector>
#include <stdexcept>
#include <memory>

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

class PauliString {
    unsigned int length;
    
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
    CFLOBDD_COMPLEX_BIG concretize() const {
        std::string name = this->name;
        unsigned index = this->index[0];
        CFLOBDD_COMPLEX_BIG res;
        std::pow(2, this->content.root->level-1);
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
            unsigned int controller = this->index[0];
            unsigned int controlled = this->index[1];
            assert(controller != controlled);
    
            if (controller < controlled)
            {
                auto C = Matrix1234ComplexFloatBoost::MkCNOT(this->content.root->level, std::pow(2, this->content.root->level - 1), controller, controlled);
                res = C;
            }
            else
            {
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(this->content.root->level, controlled, controller);
                auto C = Matrix1234ComplexFloatBoost::MkCNOT(this->content.root->level, std::pow(2, this->content.root->level - 1), controlled, controller);
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
        } else if (name == "arbitrary") {
            
        } else if (name == "swap") {
            
        } else if (name == "iswap") {
            
        } else if (name == "cz") {
            unsigned int controller = this->index[0];
            unsigned int controlled = this->index[1];
            assert(controller != controlled);
    
            if (controller < controlled)
            {
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(this->content.root->level, controller, controlled, 1.0);
                res = C;
            }
            else
            {
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(this->content.root->level, controlled, controller);
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(this->content.root->level, controlled, controller, 1.0);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = C;
            }
        } else if (name == "cp") {

        } else if (name == "cs") {
            
        } /* CCNOT, CSWAP, GlobalPhase */
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
            return QuantumGateTerm(0);
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
    SingleVecTerm(CFLOBDD_COMPLEX_BIG x) {
        // Need copy?
        // type = false;
        this->content = x;
        this->qNum = std::pow(2, x.root->level-1);
    }
    bool getType() const override {return false;}
    BIG_COMPLEX_FLOAT dot(const SingleVecTerm& other) const {
        /* bra(other) * ket(this) Note the complex number conjugation */
        // assert(this->type == false && other.type == false);
        unsigned int level = ceil(log2(this->qNum));
        auto tmpVec = Matrix1234ComplexFloatBoost::MatrixTranspose(other.content);
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
            assert(abs(amp.imag()*dimfactor) < 1e-8 && amp.real() > 0);
            double factor = double(sqrt(amp.real()));
            c1 = (1/factor) * c1;
        } else {
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
        // assert(this->type == false);
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
};

class CliffordTerm : public QuantumTerm {
    bool getType() const override {return true;}
};

class TableauTerm : public QuantumTerm {
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
    bool isIdentity = 0;
    // std::unique_ptr<Node> ast = nullptr;
    public:
    /* The type must be specified */
    QOperation() : type(0) {}
    QOperation(const QOperation& other) : type(other.type), normalized(other.normalized), qNum(other.qNum) {
        // 深拷贝 oplist
        oplist.reserve(other.oplist.size()); // 预分配空间以提高效率
        for (const auto& term : other.oplist) {
            if (term) { // 检查指针是否为空
                oplist.push_back(std::make_unique<QuantumTerm>(*term));
            } else {
                oplist.push_back(nullptr);
            }
        }
        // this->ast = std::make_unique<Node>(*other.ast);
    }
    QOperation(const QOperation& other1, const QOperation& other2) {
        assert(other1.type == other2.type);
        assert(other1.qNum == other2.qNum);
        this->qNum = other1.qNum;
        this->type = other1.type;
        this->normalized = false;
        oplist.reserve(other1.oplist.size() + other2.oplist.size());
        for (const auto& term : other1.oplist) {
            if (term) { // 检查指针是否为空
                oplist.push_back(std::make_unique<QuantumTerm>(*term));
            } else {
                oplist.push_back(nullptr);
            }
        }
        for (const auto& term : other2.oplist) {
            if (term) { // 检查指针是否为空
                oplist.push_back(std::make_unique<QuantumTerm>(*term));
            } else {
                oplist.push_back(nullptr);
            }
        }
    }
    QOperation(QOperation&& other) : type(other.type), oplist(std::move(other.oplist)), normalized(other.normalized), qNum(other.qNum) {
        // this->ast = std::make_unique<Node>(*other.ast);
        // this->ast = std::move(other.ast);
    }
    QOperation(bool t) : type(t) {}
    QOperation(bool t, std::vector<std::unique_ptr<QuantumTerm>> c) : type(t), oplist(std::move(c)), normalized(false) {
        if (!c.empty() && c[0]) {
            qNum = c[0]->qNum;
        } else {
            qNum = 0;
        }
    }
    QOperation(bool t, std::vector<std::unique_ptr<QuantumTerm>> c, bool n) : type(t), oplist(std::move(c)), normalized(n) {
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
            qNum = other.qNum;
            oplist.clear();
            oplist.reserve(other.oplist.size());
            for (const auto& term : other.oplist) {
                if (term) {
                    oplist.push_back(std::make_unique<QuantumTerm>(*term));
                } else {
                    oplist.push_back(nullptr);
                }
            }
        }
        return *this;
    }
    // void setASTdefaultADD() {
    //     // set the default AST as a ADD operation
    //     if (this->ast == nullptr) {
    //         this->ast = std::make_unique<LeafNode>(0, this->oplist.size()-1);
    //     }
    // }
    void append(std::unique_ptr<QuantumTerm> qt) {
        if (qt->getType() != this->type) {
            std::cout << "Append a wrong type of quantum term.";
            return;
        }
        oplist.push_back(std::move(qt));
    }
    std::vector<std::unique_ptr<QuantumTerm>> fetch(unsigned int begin,unsigned int end) {
        // fetch the oplist from begin to end
        assert(end < this->oplist.size());
        std::vector<std::unique_ptr<QuantumTerm>> res;
        res.reserve(end - begin + 1);
        for (size_t i = begin; i <= end; i++) {
            if (this->oplist[i] == nullptr) {
                res.push_back(nullptr);
            } else {
                res.push_back(std::make_unique<QuantumTerm>(*this->oplist[i]));
            }
        }
        return res;
    }
    SingleVecTerm projectIn(const SingleVecTerm& vec) const {
        // project the vector onto the QOperation
        assert(this->oplist[0]->getType() == false);
        assert(vec.getType() == false);
        assert(this->oplist[0]->qNum == vec.qNum);
        CFLOBDD_COMPLEX_BIG content; // Check: Initialization
        for (size_t i = 0; i < this->oplist.size(); i++) {
            auto* ivec = dynamic_cast<SingleVecTerm*>(this->oplist[i].get());
            if (!ivec) continue;
            content = content + vec.projectOnto(*ivec);
        }
        SingleVecTerm res(content);
        res.qNum = this->qNum;
        return res;
    }
    void GramSchmidt(unsigned int begin, unsigned int end) {
        // Here need to optimize the orthogonalBasis
        assert(this->type == false);
        assert(end < this->oplist.size());
        std::vector<std::unique_ptr<QuantumTerm>> orthogonalBasis;
        for (size_t i = begin; i <= end; i++) {
            auto* ivec = dynamic_cast<SingleVecTerm*>(this->oplist[i].get());
            if (!ivec) continue;
            for (const auto& basis : orthogonalBasis) {
                auto* jvec = dynamic_cast<SingleVecTerm*>(basis.get());
                if (!jvec) continue;
                ivec->content = ivec->content + (-1) * jvec->projectOnto(*ivec);
            }
            ivec->normalizeInline();
            if (!ivec->isZero()) {
                orthogonalBasis.push_back(std::move(this->oplist[i]));
            }
        }
        for (size_t i = 0; i <= (end-begin); i++) {
            this->oplist[begin+i] = std::move(orthogonalBasis[i]);
        }
        // this->oplist = std::move(orthogonalBasis);
        this->normalized = true;
    }
    QOperation negation() {
        // compute the negation of some of operators (as a disjunction)
        assert(!this->type);
        return *this;
    }
    QOperation conjunction(QOperation& other) const {
        /*** To do the conjunction:
         * 1. preserve the abstract semantic tree;
         * 2. widening function;
         * 3. hard calculate \neg(\neg A \lor \neg B).
         ***/ 
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
        assert(!other.oplist[0]->getType());
        if (!this->normalized) {
            throw std::runtime_error("The QOperation is not normalized.");
            // this->GramSchmidt(0, this->oplist.size()-1);
            // this->normalized = true;
        }
        if (!other.normalized) {
            other.GramSchmidt(0, other.oplist.size()-1);
            other.normalized = true;
        }
        QOperation op1(*this, other);
        // Modify: optimize the GramSchmidt length incremental
        op1.GramSchmidt(0, op1.oplist.size()-1);
        std::vector<std::unique_ptr<QuantumTerm>> negthis = op1.fetch(this->oplist.size(), op1.oplist.size()-1);
        QOperation op2(other, *this);
        op2.GramSchmidt(0, op2.oplist.size()-1);
        std::vector<std::unique_ptr<QuantumTerm>> negother = op2.fetch(other.oplist.size(), op2.oplist.size()-1);
        std::vector<std::unique_ptr<QuantumTerm>> remaining = op1.fetch(0, this->oplist.size()-1);
        // create a new QOperation with oplist negthis and negother
        QOperation op3(false, negthis);
        op3.oplist.insert(op3.oplist.end(), std::make_move_iterator(negother.begin()), std::make_move_iterator(negother.end()));
        op3.GramSchmidt(0, op3.oplist.size()-1);
        int dimNeg1AndNeg2 = op3.oplist.size();
        op3.oplist.insert(op3.oplist.end(), std::make_move_iterator(remaining.begin()), std::make_move_iterator(remaining.end()));
        op3.GramSchmidt(0, op3.oplist.size()-1);
        std::vector<std::unique_ptr<QuantumTerm>> resvec = op3.fetch(dimNeg1AndNeg2, op3.oplist.size()-1);
        QOperation res(false, resvec);
        res.normalized = true;
        res.qNum = this->qNum;
        return res;
    }

    QOperation disjunction(const QOperation& other) const {
        /*** To do the disjunction:
         * 1. preserve the abstract semantic tree;
         * 2. widening function;
         * 3. hard calculate \neg(\neg A \land \neg B).
         ***/
        assert(!this->type);
        if (other.oplist.size() == 0) {
            return *this;
        }
        assert(!other.oplist[0]->getType());
        QOperation res(*this, other);
        res.GramSchmidt(0, res.oplist.size()-1);
        return res;
    }

    QOperation preImage(QOperation& other) {
        /* The pre-image of a quantum operator */
        assert(this->type == false && other.type == true);
        QOperation res;
        /* Two cases here: other.oplist contains only QuantumGateTerm, or only SingleVecTerm. 
                In the later case, we should use jvec->projectOnto() */
        if (other.oplist[0]->getType() == false) {
            // Case 1: other.oplist[0] is a SingleVecTerm
            res = this->conjunction(other); // This is not a const operator!
        } else {
            // Case 2: other.oplist[0] is a QuantumGateTerm
            for (size_t i = 0; i < this->oplist.size(); i++) {
                auto* ivec = dynamic_cast<SingleVecTerm*>(this->oplist[i].get());
                if (!ivec) continue;
                
                // if (other.oplist.size() == 0) continue;
                for (size_t j = 0; j < other.oplist.size(); j++) {
                    auto* jgate = dynamic_cast<QuantumGateTerm*>(other.oplist[j].get());
                    if (!jgate) continue;
                    assert(jgate->getType() == true);
                    // Backward induction
                    auto tmp = ivec->applyGate(*jgate, false);
                    res.append(std::make_unique<SingleVecTerm>(SingleVecTerm(tmp)));
                }
            }
        }
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
            }
        }
        return res;
    }
    QOperation operator+(const QOperation& other) const {
        // Still unclear: when is the plus used?
        if (type == true) {
            if (other.type == true) {
                /* gate plus gate */
                /* append gates to the end */
                if (this->oplist.size() > 0 && other.oplist.size() > 0) {
                    assert(this->oplist[0]->getType() == other.oplist[0]->getType());
                }
                std::vector<std::unique_ptr<QuantumTerm>> newlist = this->oplist;
                newlist.insert(newlist.end(), other.oplist.begin(), other.oplist.end());
                QOperation res(true, newlist);
                return res;
            } else {
                /* gate plus cflobdd */
                throw std::runtime_error("Cannot support plus between gates and operators");
            }
        } else {
            if (other.type == true) {
                /* cflobdd plus cflobdd */
                /* disjunction */
                std::vector<std::unique_ptr<QuantumTerm>> newlist = this->oplist;
                newlist.insert(newlist.end(), other.oplist.begin(), other.oplist.end());
                QOperation res(true, newlist, false);
                return res;
            } else {
                /* cflobdd plus gate */
                throw std::runtime_error("Cannot support plus between gates and operators");
            }
        }
    }
    QOperation operator*(const QOperation& other) const {
        if (type == true) {
            if (other.type == true) {
                /* gate times gate */
                /* cascade */
                // throw std::runtime_error("Logically, it should be supported");
                for (size_t i = 0; i < this->oplist.size(); i++) {
                    auto* gate1 = dynamic_cast<QuantumGateTerm*>(this->oplist[i].get());
                    for (size_t j = 0; i < other.oplist.size(); j++) {
                        CFLOBDD_COMPLEX_BIG newItem;
                        auto* gate2 = dynamic_cast<QuantumGateTerm*>(other.oplist[j].get());
                        if (this->oplist[0]->getType() == true) {
                            // Assert all type in a operation are consistant
                            // all operators are gates
                            auto res = gate1->cascade(*gate2);
                            if (res.isConcret) {
                                newItem = res.concretize();
                            } else {
                                newItem = res.content;
                            }
                        } else {
                            throw std::runtime_error("Cannot support the cascade of two operators");
                        }
                    }
                    
                }
                
            } else {
                /* gate times cflobdd */
                QOperation res;
                for (size_t i = 0; i < other.oplist.size(); i++) {
                    CFLOBDD_COMPLEX_BIG newItem;
                    for (size_t j = 0; j < this->oplist.size(); j++) {
                        if (this->oplist[0]->getType() == true) {
                            // Assert all type in a operation are consistant
                            // all operators are gates
                            auto* gate1 = dynamic_cast<QuantumGateTerm*>(this->oplist[j].get());
                        } else {
                            // all operators are supporting vector of a space
                        }
                        
                    }
                }
                
            }
        } else {
            if (other.type == true) {
                /* cflobdd times cflobdd */
                /* conjunction */
            } else {
                /* cflobdd times gate*/
            }
        }
    }
    
    
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
                if (i == this->oplist.size()-1) {
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
};

/* Const QOperation:
I: the identity operator as the top of the QOperation lattice
  --> properties: for all QOperation A, conjunction(A, I) = A
  --> properties: for all QOperation A, disjunction(A, I) = I
*/

QOperation CreateIdentityQO(unsigned int qNum) {
    QOperation res(true);
    res.qNum = qNum;
    res.isIdentity = true;
    return res;
}

QOperation CreateZeroQO(unsigned int qNum) {
    QOperation res(true);
    res.qNum = qNum;
    return res;
}

#endif