#ifndef _QUANTUM_OPERATION
#define _QUANTUM_OPERATION

#include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include "cflobdd/CFLOBDD/vector_complex_float_boost.h"
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
    // construct as a gate
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
    QuantumGateTerm cascade(QuantumGateTerm other) const {
        return;
    }
};

class SingleVecTerm : public QuantumTerm {
    public:
    // construct as a projector or state
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
    CFLOBDD_COMPLEX_BIG applyGate(const QuantumGateTerm& other) const {
        /* This function is INPLACE! */
        CFLOBDD_COMPLEX_BIG operand = other.concretize();
        std::string name = other.name;
        unsigned index = other.index[0];
        CFLOBDD_COMPLEX_BIG res;
        if (name == "x" || name == "y" || name == "z" || name == "h" || name == "s" || name == "t" || name == "p") {
            auto X = ApplyGateF(std::pow(2, this->content.root->level-1), index, Matrix1234ComplexFloatBoost::MkNegationMatrixInterleaved);
            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(operand, this->content);
        } else if (name == "y") {
            auto Y = ApplyGateF(std::pow(2, this->content.root->level-1), index, Matrix1234ComplexFloatBoost::MkPauliYMatrixInterleaved);
            int size = Y.root->rootConnection.returnMapHandle.Size();
            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(Y, this->content);
        } else if (name == "z") {
            auto Z = ApplyGateF(std::pow(2, this->content.root->level-1), index, Matrix1234ComplexFloatBoost::MkPauliZMatrixInterleaved);
            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(Z, this->content);
        } else if (name == "h") {
            auto H = ApplyGateF(std::pow(2, this->content.root->level-1), index, Matrix1234ComplexFloatBoost::MkWalshInterleaved);
            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, this->content);
        } else if (name == "s") {
            auto S = ApplyGateF(std::pow(2, this->content.root->level-1), index, Matrix1234ComplexFloatBoost::MkSGateInterleaved);
            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, this->content);
        } else if (name == "t") {
            auto S = ApplyGateFWithParam(std::pow(2, this->content.root->level-1), index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, 0.25);
            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, this->content); 
        } else if (name == "p") {
            double theta = other.vars[0];
            auto S = ApplyGateFWithParam(std::pow(2, this->content.root->level-1), index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, theta);
            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, this->content);
        } else if (name == "cx") {
            unsigned int controller = other.index[0];
            unsigned int controlled = other.index[1];
            assert(controller != controlled);
    
            if (controller < controlled)
            {
                auto C = Matrix1234ComplexFloatBoost::MkCNOT(this->content.root->level, std::pow(2, this->content.root->level - 1), controller, controlled);
                res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, this->content);
            }
            else
            {
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(this->content.root->level, controlled, controller);
                auto C = Matrix1234ComplexFloatBoost::MkCNOT(this->content.root->level, std::pow(2, this->content.root->level - 1), controlled, controller);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, this->content);
            }
        } else if (name == "u3") {
            assert(other.vars.size() == 3);
            std::vector<double> v;
            double theta = other.vars[0];
            double phi = other.vars[1];
            double lambda = other.vars[2];
            v.push_back(theta); v.push_back(phi); v.push_back(lambda);
            auto U = ApplyGateFWithParamVec(std::pow(2, this->content.root->level-1), index, Matrix1234ComplexFloatBoost::MkU3GateInterleaved, v);
            res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(U, this->content);
        } else if (name == "arbitrary") {
            
        } else if (name == "swap") {
            
        } else if (name == "iswap") {
            
        } else if (name == "cz") {
            unsigned int controller = other.index[0];
            unsigned int controlled = other.index[1];
            assert(controller != controlled);
    
            if (controller < controlled)
            {
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(this->content.root->level, controller, controlled, 1.0);
                res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, this->content);
            }
            else
            {
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(this->content.root->level, controlled, controller);
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(this->content.root->level, controlled, controller, 1.0);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, this->content);
            }
        } else if (name == "cp") {

        } else if (name == "cs") {
            
        } /* CCNOT, CSWAP, GlobalPhase */
        return res;
    }
    void applyGateInline(const QuantumGateTerm& other) {
        this->content = this->applyGate(other);
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
    bool type;
    std::vector<std::unique_ptr<QuantumTerm>> oplist;
    bool normalized = 0;
    unsigned int qNum = 0;
    public:
    /* The type must be specified */
    // QOperation() {}
    QOperation() : type(0) {}
    QOperation(bool t) : type(t) {}
    QOperation(bool t, std::vector<std::unique_ptr<QuantumTerm>> c) : type(t), oplist(c) {
        this->qNum = c[0]->qNum;
    }
    QOperation(bool t, std::vector<std::unique_ptr<QuantumTerm>> c, bool n) : type(t), oplist(c), normalized(n) {
        this->qNum = c[0]->qNum;
    }
    void append(std::unique_ptr<QuantumTerm> qt) {
        if (qt->getType() != this->type) {
            std::cout << "Append a wrong type of quantum term.";
            return;
        }
        oplist.push_back(qt);
    }
    void GramSchmidt(unsigned int begin, unsigned int end) {
        // A better sulotion is to record every index of operands. Because it may not be continuous.
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
    }
    QOperation negation(unsigned int begin, unsigned int end) {
        // compute the negation of some of operators (as a disjunction)
        
        return;
    }
    QOperation conjunction(unsigned int begin, unsigned int end) {
        /*** To do the conjunction:
         * 1. preserve the abstract semantic tree;
         * 2. widening function;
         * 3. hard calculate \neg(\neg A \lor \neg B).
         ***/ 
    }
    QOperation operator+(const QOperation& other) const {
        if (type == true) {
            if (other.type == true) {
                /* gate plus gate */
                /* append gates to the end */
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
                    for (size_t j = 0; i < other.oplist.size(); j++) {
                        /* code */
                    }
                    
                }
                
            } else {
                /* gate times cflobdd */
                QOperation res;
                for (size_t i = 0; i < other.oplist.size(); i++) {
                    for (size_t j = 0; j < this->oplist.size(); j++) {
                        CFLOBDD_COMPLEX_BIG newItem;
                        if (this->oplist[0]->getType() == true) {
                            // Assert all type in a operation are consistant
                        } else {
                            
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
    bool satisfy(const QOperation& other) const {
        /* If the support space of this is the subspace of other */
        return false;
    }
};

#endif