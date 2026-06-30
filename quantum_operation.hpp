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
#include <cstdlib>
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

CFLOBDD_COMPLEX_BIG InitializeWithVector(unsigned int qnum, std::vector<double> vec_raw) {
    // The length of vec must be a power of 2
    // TODO: the qubits of the vector may not be a power of 2!!!
    unsigned int n = vec_raw.size();
    assert((n & (n - 1)) == 0 && n != 0);
    // assert(n == 4);
    assert(2 * (1 << qnum) == n);
    std::vector<std::complex<double>> vec(n/2);
    for (unsigned int i = 0; i < n/2; i++) {
        vec[i] = std::complex<double>(vec_raw[i], vec_raw[i+n/2]);
    }
    // assert vec is normalized
    double norm = 0;
    for (unsigned int i = 0; i < n/2; i++) {
        norm += std::norm(vec[i]);
    }
    assert(abs(norm - 1.0) < 1e-8);
    unsigned int level = ceil(log2(qnum));
    CFLOBDD_COMPLEX_BIG res = VectorComplexFloatBoost::NoDistinctionNode(level, 0);
    for (unsigned int i = 0; i < n/2; i++) {
        if (abs(vec[i]) > 1e-10) {
            auto basisVec = VectorComplexFloatBoost::MkBasisVector(level, i);
            auto scaledVec = BIG_COMPLEX_FLOAT(vec[i]) * basisVec;
            res = res + scaledVec;
        }
    }
    res = VectorComplexFloatBoost::VectorToMatrixInterleaved(res);
    // VectorComplexFloatBoost::VectorPrintColumnHead(res, std::cout);
    return res;
}

std::string stringPadding(std::string str, unsigned int length) {
    if (str.length() >= length) {
        return str;
    }
    return str + std::string(length - str.length(), '0');
}

bool isSimpleProductStateString(const std::string& str) {
    if (str.empty()) {
        return false;
    }
    for (char c : str) {
        if (c != '0' && c != '1' && c != '+' && c != '-') {
            return false;
        }
    }
    return true;
}

bool hasHadamardBasisSymbol(const std::string& str) {
    return str.find('+') != std::string::npos || str.find('-') != std::string::npos;
}

std::vector<double> simpleProductStateAmplitudes(const std::string& str, unsigned int physicalQubits) {
    assert(isSimpleProductStateString(str));
    assert(physicalQubits >= str.size());
    std::string padded = stringPadding(str, physicalQubits);
    unsigned int basisSize = 1 << physicalQubits;
    std::vector<double> realAmps(basisSize, 0.0);

    unsigned int hadamardCount = 0;
    for (char c : padded) {
        if (c == '+' || c == '-') {
            ++hadamardCount;
        }
    }
    double baseAmp = 1.0 / std::sqrt(static_cast<double>(1 << hadamardCount));

    for (unsigned int basis = 0; basis < basisSize; ++basis) {
        double amp = baseAmp;
        bool compatible = true;
        for (unsigned int pos = 0; pos < physicalQubits; ++pos) {
            unsigned int bit = (basis >> (physicalQubits - pos - 1)) & 1U;
            char symbol = padded[pos];
            if ((symbol == '0' && bit != 0) || (symbol == '1' && bit != 1)) {
                compatible = false;
                break;
            }
            if (symbol == '-' && bit == 1) {
                amp = -amp;
            }
        }
        if (compatible) {
            realAmps[basis] = amp;
        }
    }

    std::vector<double> amps;
    amps.reserve(2 * basisSize);
    amps.insert(amps.end(), realAmps.begin(), realAmps.end());
    amps.resize(2 * basisSize, 0.0);
    return amps;
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

namespace qoprof {

struct Stats {
    size_t qoperation_default_ctor = 0;
    size_t qoperation_basis_ctor = 0;
    size_t qoperation_amplitude_ctor = 0;
    size_t qoperation_gate_ctor = 0;
    size_t qoperation_copy_ctor = 0;
    size_t qoperation_merge_ctor = 0;
    size_t qoperation_move_ctor = 0;
    size_t qoperation_bool_ctor = 0;
    size_t qoperation_vector_ctor = 0;
    size_t qoperation_vector_norm_ctor = 0;
    size_t qoperation_copy_assign = 0;
    size_t qoperation_dtor = 0;
    size_t qoperation_alive = 0;
    size_t qoperation_peak_alive = 0;
    size_t qoperation_live_terms = 0;
    size_t qoperation_peak_live_terms = 0;

    size_t quantum_gate_term_clones = 0;
    size_t single_vec_term_clones = 0;

    size_t append_calls = 0;
    size_t append_terms = 0;
    size_t add_calls = 0;
    size_t add_terms_cloned = 0;
    size_t fetch_calls = 0;
    size_t fetch_terms_cloned = 0;
    size_t gen_proj_meas_space_calls = 0;
    size_t gen_proj_generated_terms = 0;

    size_t gram_schmidt_calls = 0;
    size_t gram_schmidt_total_input_terms = 0;
    size_t gram_schmidt_total_output_terms = 0;
    size_t gram_schmidt_peak_input_terms = 0;
    size_t gram_schmidt_peak_output_terms = 0;

    size_t minus_calls = 0;
    size_t conjunction_calls = 0;
    size_t conjunction_simp_calls = 0;
    size_t disjunction_calls = 0;
    size_t preimage_calls = 0;
    size_t postimage_calls = 0;
    size_t postimage_total_input_terms = 0;
    size_t postimage_total_gate_terms = 0;
    size_t postimage_total_raw_terms = 0;
    size_t postimage_total_final_terms = 0;
    size_t postimage_peak_input_terms = 0;
    size_t postimage_peak_gate_terms = 0;
    size_t postimage_peak_raw_terms = 0;
    size_t postimage_peak_final_terms = 0;
};

inline bool enabled() {
    static const bool enabled_flag = []() {
        const char* value = std::getenv("QOP_PROFILE");
        return value && std::string(value) != "0";
    }();
    return enabled_flag;
}

inline bool trace_enabled() {
    static const bool trace_flag = []() {
        const char* value = std::getenv("QOP_PROFILE_TRACE");
        return value && std::string(value) != "0";
    }();
    return trace_flag;
}

inline Stats& stats() {
    static Stats value;
    return value;
}

inline void on_qoperation_create(size_t terms, size_t* counter) {
    if (!enabled()) return;
    Stats& s = stats();
    ++(*counter);
    ++s.qoperation_alive;
    s.qoperation_peak_alive = std::max(s.qoperation_peak_alive, s.qoperation_alive);
    s.qoperation_live_terms += terms;
    s.qoperation_peak_live_terms = std::max(s.qoperation_peak_live_terms, s.qoperation_live_terms);
}

inline void on_qoperation_move_create(size_t terms) {
    if (!enabled()) return;
    Stats& s = stats();
    ++s.qoperation_move_ctor;
    ++s.qoperation_alive;
    s.qoperation_peak_alive = std::max(s.qoperation_peak_alive, s.qoperation_alive);
    s.qoperation_peak_live_terms = std::max(s.qoperation_peak_live_terms, s.qoperation_live_terms);
    (void)terms;
}

inline void on_qoperation_destroy(size_t terms) {
    if (!enabled()) return;
    Stats& s = stats();
    ++s.qoperation_dtor;
    if (s.qoperation_alive > 0) {
        --s.qoperation_alive;
    }
    if (s.qoperation_live_terms >= terms) {
        s.qoperation_live_terms -= terms;
    } else {
        s.qoperation_live_terms = 0;
    }
}

inline void on_qoperation_copy_assign(size_t old_terms, size_t new_terms) {
    if (!enabled()) return;
    Stats& s = stats();
    ++s.qoperation_copy_assign;
    if (new_terms >= old_terms) {
        s.qoperation_live_terms += (new_terms - old_terms);
    } else {
        s.qoperation_live_terms -= (old_terms - new_terms);
    }
    s.qoperation_peak_live_terms = std::max(s.qoperation_peak_live_terms, s.qoperation_live_terms);
}

inline void on_quantum_gate_clone() {
    if (!enabled()) return;
    ++stats().quantum_gate_term_clones;
}

inline void on_single_vec_clone() {
    if (!enabled()) return;
    ++stats().single_vec_term_clones;
}

inline void on_append() {
    if (!enabled()) return;
    ++stats().append_calls;
    ++stats().append_terms;
}

inline void on_add(size_t cloned_terms) {
    if (!enabled()) return;
    ++stats().add_calls;
    stats().add_terms_cloned += cloned_terms;
}

inline void on_fetch(size_t cloned_terms) {
    if (!enabled()) return;
    ++stats().fetch_calls;
    stats().fetch_terms_cloned += cloned_terms;
}

inline void on_gen_proj_meas_space(size_t generated_terms) {
    if (!enabled()) return;
    ++stats().gen_proj_meas_space_calls;
    stats().gen_proj_generated_terms += generated_terms;
}

inline void on_gram_schmidt(size_t input_terms, size_t output_terms) {
    if (!enabled()) return;
    Stats& s = stats();
    ++s.gram_schmidt_calls;
    s.gram_schmidt_total_input_terms += input_terms;
    s.gram_schmidt_total_output_terms += output_terms;
    s.gram_schmidt_peak_input_terms = std::max(s.gram_schmidt_peak_input_terms, input_terms);
    s.gram_schmidt_peak_output_terms = std::max(s.gram_schmidt_peak_output_terms, output_terms);
    if (trace_enabled()) {
        std::cout << "[qoprof-trace] gramschmidt input=" << input_terms
                  << " output=" << output_terms
                  << " peak_input=" << s.gram_schmidt_peak_input_terms
                  << " peak_output=" << s.gram_schmidt_peak_output_terms
                  << std::endl;
    }
}

inline void on_minus_call() {
    if (!enabled()) return;
    ++stats().minus_calls;
}

inline void on_conjunction_call() {
    if (!enabled()) return;
    ++stats().conjunction_calls;
}

inline void on_conjunction_simp_call() {
    if (!enabled()) return;
    ++stats().conjunction_simp_calls;
}

inline void on_disjunction_call() {
    if (!enabled()) return;
    ++stats().disjunction_calls;
}

inline void on_preimage_call() {
    if (!enabled()) return;
    ++stats().preimage_calls;
}

inline void on_postimage_call(size_t input_terms, size_t gate_terms, size_t raw_terms, size_t final_terms) {
    if (!enabled()) return;
    Stats& s = stats();
    ++s.postimage_calls;
    s.postimage_total_input_terms += input_terms;
    s.postimage_total_gate_terms += gate_terms;
    s.postimage_total_raw_terms += raw_terms;
    s.postimage_total_final_terms += final_terms;
    s.postimage_peak_input_terms = std::max(s.postimage_peak_input_terms, input_terms);
    s.postimage_peak_gate_terms = std::max(s.postimage_peak_gate_terms, gate_terms);
    s.postimage_peak_raw_terms = std::max(s.postimage_peak_raw_terms, raw_terms);
    s.postimage_peak_final_terms = std::max(s.postimage_peak_final_terms, final_terms);
    if (trace_enabled()) {
        std::cout << "[qoprof-trace] postimage input_terms=" << input_terms
                  << " gate_terms=" << gate_terms
                  << " raw_terms=" << raw_terms
                  << " final_terms=" << final_terms
                  << " peak_raw=" << s.postimage_peak_raw_terms
                  << " peak_final=" << s.postimage_peak_final_terms
                  << " calls=" << s.postimage_calls
                  << std::endl;
    }
}

inline void print_summary(std::ostream& os = std::cout) {
    const Stats& s = stats();
    os << "[qoprof] qoperation_default_ctor=" << s.qoperation_default_ctor
       << " basis_ctor=" << s.qoperation_basis_ctor
       << " amplitude_ctor=" << s.qoperation_amplitude_ctor
       << " gate_ctor=" << s.qoperation_gate_ctor
       << " copy_ctor=" << s.qoperation_copy_ctor
       << " merge_ctor=" << s.qoperation_merge_ctor
       << " move_ctor=" << s.qoperation_move_ctor
       << " bool_ctor=" << s.qoperation_bool_ctor
       << " vector_ctor=" << s.qoperation_vector_ctor
       << " vector_norm_ctor=" << s.qoperation_vector_norm_ctor
       << " copy_assign=" << s.qoperation_copy_assign
       << " dtor=" << s.qoperation_dtor
       << " alive_peak=" << s.qoperation_peak_alive
       << " live_terms_peak=" << s.qoperation_peak_live_terms
       << std::endl;

    os << "[qoprof] gate_clones=" << s.quantum_gate_term_clones
       << " vec_clones=" << s.single_vec_term_clones
       << " append_calls=" << s.append_calls
       << " append_terms=" << s.append_terms
       << " add_calls=" << s.add_calls
       << " add_terms_cloned=" << s.add_terms_cloned
       << " fetch_calls=" << s.fetch_calls
       << " fetch_terms_cloned=" << s.fetch_terms_cloned
       << " gen_proj_calls=" << s.gen_proj_meas_space_calls
       << " gen_proj_terms=" << s.gen_proj_generated_terms
       << std::endl;

    os << "[qoprof] gramschmidt_calls=" << s.gram_schmidt_calls
       << " gramschmidt_input_total=" << s.gram_schmidt_total_input_terms
       << " gramschmidt_output_total=" << s.gram_schmidt_total_output_terms
       << " gramschmidt_input_peak=" << s.gram_schmidt_peak_input_terms
       << " gramschmidt_output_peak=" << s.gram_schmidt_peak_output_terms
       << std::endl;

    os << "[qoprof] minus_calls=" << s.minus_calls
       << " conjunction_calls=" << s.conjunction_calls
       << " conjunction_simp_calls=" << s.conjunction_simp_calls
       << " disjunction_calls=" << s.disjunction_calls
       << " preimage_calls=" << s.preimage_calls
       << " postimage_calls=" << s.postimage_calls
       << " postimage_input_total=" << s.postimage_total_input_terms
       << " postimage_gate_total=" << s.postimage_total_gate_terms
       << " postimage_raw_total=" << s.postimage_total_raw_terms
       << " postimage_final_total=" << s.postimage_total_final_terms
       << " postimage_input_peak=" << s.postimage_peak_input_terms
       << " postimage_gate_peak=" << s.postimage_peak_gate_terms
       << " postimage_raw_peak=" << s.postimage_peak_raw_terms
       << " postimage_final_peak=" << s.postimage_peak_final_terms
       << std::endl;
}

struct Reporter {
    ~Reporter() {
        if (enabled()) {
            print_summary();
        }
    }
};

inline Reporter reporter;

} // namespace qoprof

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
    bool isEqual(const QuantumGateTerm& other) const {
        return this->name == other.name && this->index == other.index && this->vars == other.vars && this->qNum == other.qNum;
    }
    std::unique_ptr<QuantumTerm> clone() const override {
        qoprof::on_quantum_gate_clone();
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
        } else if (name == "sdg") {
            auto Sdg = ApplyGateFWithParam(this->qNum, index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, -0.5);
            res = Sdg;
        } else if (name == "t") {
            auto S = ApplyGateFWithParam(this->qNum, index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, 0.25);
            res = S; 
        } else if (name == "p") {
            double theta = this->vars[0];
            auto S = ApplyGateFWithParam(this->qNum, index, Matrix1234ComplexFloatBoost::MkPhaseShiftGateInterleaved, theta);
            res = S;
        } else if (name == "sx") {
            auto H = ApplyGateF(this->qNum, index, Matrix1234ComplexFloatBoost::MkWalshInterleaved);
            auto S = ApplyGateF(this->qNum, index, Matrix1234ComplexFloatBoost::MkSGateInterleaved);
            S = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, S);
            S = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, H);
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
        } else if (name == "csx") {
            // Use the identity CSX = H(target) CP[pi/2](control, target) H(target)
            assert(this->qNum && (this->qNum & (this->qNum - 1)) == 0);
            unsigned int controller = this->index[0];
            unsigned int controlled = this->index[1];
            unsigned int state_level = ceil(log2(this->qNum)) + 1;
            assert(controller != controlled);

            auto H = ApplyGateF(this->qNum, controlled, Matrix1234ComplexFloatBoost::MkWalshInterleaved);
            if (controller < controlled)
            {
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(state_level, controller, controlled, 1/2);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, C);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, H);
                res = C;
            }
            else
            {
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, controlled, controller);
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(state_level, controlled, controller, 1/2);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(H, C);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, H);
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
        } else if (name == "arb") {
            assert(this->vars.size() == 8);
            auto U = ApplyGateFWithParamVec(this->qNum, index, Matrix1234ComplexFloatBoost::MkArbitraryGateInterleaved, this->vars);
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
        } else if (name == "resetall") {
            CFLOBDD_COMPLEX_BIG U = VectorComplexFloatBoost::NoDistinctionNode(ceil(log2(this->qNum)), 1);
            U = VectorComplexFloatBoost::VectorToMatrixInterleaved(U);
            U = Matrix1234ComplexFloatBoost::MatrixConjugate(U);
            U = Matrix1234ComplexFloatBoost::MatrixTranspose(U);
            res = U;
        } else if (name == "init") {
            // Assume the state is in a tensor state, all indexed qubits are |0>
            // Assume the indexes are sequentially ordered
            // Prepare the state from this->vars
            // TODO: Out of compromise, we make the following assumption: only one qubit is initialized.
            unsigned int numVars = this->vars.size();
            assert(numVars && (numVars & (numVars - 1)) == 0);
            // stateVec: [a+bi, c+di]
            // auto stateVec = InitializeWithVector(1, this->vars);
            // Prepare the projecor |init><0|
            // stateVec = VectorComplexFloatBoost::VectorToMatrixInterleaved(stateVec);
            // Append this->vars with zeros, as the imaginary part, double the size to vec_raw.
            auto vec_raw = this->vars;
            vec_raw.resize(2 * vec_raw.size(), 0.0);
            auto U = ApplyGateFWithParamVec(this->qNum, index, InitializeWithVector, vec_raw);
            // Check the indexes that before and after the applied indexes, padding them (through tensor) with identity
            res = U;
            
        } else if (name == "swap") {
            assert(this->qNum && (this->qNum & (this->qNum - 1)) == 0);
            unsigned int index1 = this->index[0];
            unsigned int index2 = this->index[1];
            unsigned int state_level = ceil(log2(this->qNum)) + 1;
            assert(index1 != index2);
            if (index1 < index2)
            {
                auto C = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, index1, index2);
                res = C;
            }
            else
            {
                auto C = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, index2, index1);
                res = C;
            }
        } else if (name == "iswap") {
            assert(this->qNum && (this->qNum & (this->qNum - 1)) == 0);
            unsigned int index1 = this->index[0];
            unsigned int index2 = this->index[1];
            unsigned int state_level = ceil(log2(this->qNum)) + 1;
            assert(index1 != index2);
            if (index1 < index2)
            {
                auto C = Matrix1234ComplexFloatBoost::MkiSwapGate(state_level, index1, index2);
                res = C;
            }
            else
            {
                auto C = Matrix1234ComplexFloatBoost::MkiSwapGate(state_level, index2, index1);
                res = C;
            }
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
            assert(this->qNum && (this->qNum & (this->qNum - 1)) == 0);
            unsigned int controller = this->index[0];
            unsigned int controlled = this->index[1];
            unsigned int state_level = ceil(log2(this->qNum)) + 1;
            double theta = this->vars[0];
            assert(controller != controlled);

            if (controller < controlled)
            {
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(state_level, controller, controlled, theta);
                res = C;
            }
            else
            {
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, controlled, controller);
                auto C = Matrix1234ComplexFloatBoost::MkCPGate(state_level, controlled, controller, theta);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = C;
            }
        } else if (name == "cs") {
            
        } /* CCNOT, CSWAP */ else if (name == "ccx") {
            assert(this->qNum && (this->qNum & (this->qNum - 1)) == 0);
            unsigned int controller1 = this->index[0];
            unsigned int controller2 = this->index[1];
            unsigned int controlled = this->index[2];
            unsigned int state_level = ceil(log2(this->qNum)) + 1;
            assert(controller1 != controlled);
            assert(controller2 != controlled);
            assert(controller1 != controller2);
            if (controller1 < controller2 && controller2 < controlled)
            {
                // a b c
                auto C = Matrix1234ComplexFloatBoost::MkCCNOT(state_level, std::pow(2, state_level - 1), controller1, controller2, controlled);
                res = C;
            }
            else if (controller1 < controlled && controlled < controller2)
            {
                // a c b   
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, controlled, controller2);
                auto C = Matrix1234ComplexFloatBoost::MkCCNOT(state_level, std::pow(2, state_level - 1), controller1, controlled, controller2);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = C;
            }
            else if (controller2 < controller1 && controller1 < controlled)
            {
                // b a c
                auto C = Matrix1234ComplexFloatBoost::MkCCNOT(state_level, std::pow(2, state_level - 1), controller2, controller1, controlled);
                res = C;
            }
            else if (controller2 < controlled && controlled < controller1)
            {
                // b c a
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, controlled, controller1);
                auto C = Matrix1234ComplexFloatBoost::MkCCNOT(state_level, std::pow(2, state_level - 1), controller2, controlled, controller1);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = C;
            }
            else if (controlled < controller1 && controller1 < controller2)
            {
                // c a b
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, controlled, controller2);
                // b a c
                auto C = Matrix1234ComplexFloatBoost::MkCCNOT(state_level, std::pow(2, state_level - 1), controlled, controller1, controller2);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = C;
            }
            else if (controlled < controller2 && controller2 < controller1)
            {
                // c b a
                auto S = Matrix1234ComplexFloatBoost::MkSwapGate(state_level, controlled, controller1);
                // a b c
                auto C = Matrix1234ComplexFloatBoost::MkCCNOT(state_level, std::pow(2, state_level - 1), controlled, controller2, controller1);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(C, S);
                C = Matrix1234ComplexFloatBoost::MatrixMultiplyV4WithInfo(S, C);
                res = C;
            }
        }
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
    SingleVecTerm(std::vector<double> amp, unsigned int qubits) {
        // qubits is the realQubits, qNum is the total physical qubits.
        unsigned level = ceil(log2(qubits));
        this->qNum = std::pow(2, level);
        assert(amp.size() == 2 * (1 << qubits));
        // Padding the amp to 2^(this->qNum) with 0s
        while(amp.size() < (1 << this->qNum)) {
            amp.push_back(0);
        }
        CFLOBDD_COMPLEX_BIG stateVector = InitializeWithVector(this->qNum, amp);
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
        qoprof::on_single_vec_clone();
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
    // realqNum may not be an exp of 2, it is the real number of qubits in the ideal algorithm.
    unsigned int realqNum = 0;
    bool isIdentity = false;
    int isProj = -1;
    // std::unique_ptr<Node> ast = nullptr;
    public:
    /* The type must be specified */
    QOperation() : type(0) { qoprof::on_qoperation_create(oplist.size(), &qoprof::stats().qoperation_default_ctor); }
    QOperation(std::vector<std::string> strings) {
        // Construct a QOperation of simple product states.
        // A string may contain 0/1 computational-basis symbols and +/- Hadamard-basis symbols:
        // '+' means H|0> = (|0> + |1>)/sqrt(2), '-' means H|1> = (|0> - |1>)/sqrt(2).
        assert(strings.size() > 0);
        for (const auto& str : strings) {
            assert(str.size() > 0);
            assert(str.size() == strings[0].size());
            assert(isSimpleProductStateString(str));
        }
        this->type = false;
        this->realqNum = strings[0].size();
        this->qNum = std::pow(2, ceil(log2(strings[0].size())));
        for (const auto& str : strings) {
            if (hasHadamardBasisSymbol(str)) {
                SingleVecTerm term(simpleProductStateAmplitudes(str, this->qNum), this->qNum);
                oplist.push_back(std::make_unique<SingleVecTerm>(term));
            } else {
                SingleVecTerm term(stringPadding(str, this->qNum), std::pow(2, ceil(log2(str.size()))));
                oplist.push_back(std::make_unique<SingleVecTerm>(term));
            }
        }
        this->normalized = true;
        qoprof::on_qoperation_create(oplist.size(), &qoprof::stats().qoperation_basis_ctor);
    }
    QOperation(std::vector<double> amps, unsigned int qubits) {
        // Construct a QOperation of arbitrary amplitude vectors.
        // qubits is the real qubits, qNum is the CLFOBDD qubits.
        // assert it is normalized.
        assert(amps.size() > 0);
        assert(amps.size() == 2 * (1 << qubits));
        this->type = false;
        this->realqNum = qubits;
        this->qNum = std::pow(2, ceil(log2(qubits)));
        // Tensor the amps to the next power of 2 with |0>s, that is, padding 0s after every elements. For example, qubits = 3, amps = [1,1,0,0,1,1,0,0], 
        // then we tensor it to qubits = 4, amps = [1,0,1,0,0,0,0,0,1,0,1,0,0,0,0,0]
        std::vector<double> newamps;
        unsigned int targetSize = 2 * (1 << this->qNum);
        for (unsigned int i = 0; i < amps.size(); i++) {
            newamps.push_back(amps[i]);
            // padding
            for (unsigned int j = 0; j < (targetSize / amps.size() - 1); j++) {
                newamps.push_back(0);
            }
        }
        assert(newamps.size() == targetSize);
        SingleVecTerm term(newamps, this->qNum);
        oplist.push_back(std::make_unique<SingleVecTerm>(term));
        this->normalized = true;
        qoprof::on_qoperation_create(oplist.size(), &qoprof::stats().qoperation_amplitude_ctor);
    }
    QOperation(std::string nam, unsigned int qNum, std::vector<unsigned int> idx, std::vector<double> pars=std::vector<double>{}) {
        /*
        Construct a QOperation of quantum gate.
        Params:
            nam: the name of the gate, such as "X", "H", "CX", "U3", "meas0", "reset", etc.
            qNum: the real qubits of the operation.
            idx: the indexes of the qubits that the gate acts on.
            pars: the parameters of the gate, such as the angles of the U3 gate.
        */
        this->type = true;
        this->realqNum = qNum;
        unsigned int logicqNum = std::pow(2, ceil(log2(qNum)));
        this->qNum = logicqNum;
        // Just a copy of CreateProjectiveMeasQO
        if (nam == "meas0") {
            this->isProj = idx[0];
            this->normalized = true;
            QuantumGateTerm tmp("meas0", std::vector<unsigned int>{idx[0]}, std::vector<double>{}, logicqNum);
            // Whether it is necessary to concretize?
            tmp.concretizeInline();
            // Attention!
            this->oplist.push_back(tmp.clone());
        } else if (nam == "meas1") {
            this->isProj = idx[0];
            this->normalized = true;
            QuantumGateTerm tmp("meas1", std::vector<unsigned int>{idx[0]}, std::vector<double>{}, logicqNum);
            tmp.concretizeInline();
            // Attention!
            this->oplist.push_back(tmp.clone());
        } else if (nam == "reset") {
            // The reset operator creates a mixed state
            // This is another projective operation!!!
            this->isProj = idx[0];
            this->normalized = true;
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
        qoprof::on_qoperation_create(oplist.size(), &qoprof::stats().qoperation_gate_ctor);
    }
    QOperation(const QOperation& other) : type(other.type), normalized(other.normalized), qNum(other.qNum), isIdentity(other.isIdentity), isProj(other.isProj), realqNum(other.realqNum) {
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
        qoprof::on_qoperation_create(oplist.size(), &qoprof::stats().qoperation_copy_ctor);
    }
    QOperation(const QOperation& other1, const QOperation& other2) {
        assert(other1.type == other2.type);
        if (other1.qNum != other2.qNum) {
            std::cout << other1.qNum << " != " << other2.qNum << std::endl;
        }
        assert(other1.qNum == other2.qNum);
        this->qNum = other1.qNum;
        this->realqNum = other1.realqNum;
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
        qoprof::on_qoperation_create(oplist.size(), &qoprof::stats().qoperation_merge_ctor);
    }
    QOperation(QOperation&& other) noexcept : type(other.type), oplist(std::move(other.oplist)), normalized(other.normalized), qNum(other.qNum), isIdentity(other.isIdentity), isProj(other.isProj), realqNum(other.realqNum) {
        // Move constructor
        // this->ast = std::make_unique<Node>(*other.ast);
        // this->ast = std::move(other.ast);
        qoprof::on_qoperation_move_create(oplist.size());
    }
    QOperation(bool t) : type(t) { qoprof::on_qoperation_create(oplist.size(), &qoprof::stats().qoperation_bool_ctor); }
    QOperation(bool t, std::vector<std::unique_ptr<QuantumTerm>>&& c) : type(t), oplist(std::move(c)), normalized(false) {
        // Move constructor
        if (!c.empty() && c[0]) {
            qNum = c[0]->qNum;
            realqNum = c[0]->qNum;
        } else {
            qNum = 0;
            realqNum = 0;
        }
        qoprof::on_qoperation_create(oplist.size(), &qoprof::stats().qoperation_vector_ctor);
    }
    QOperation(bool t, std::vector<std::unique_ptr<QuantumTerm>>&& c, bool n) : type(t), oplist(std::move(c)), normalized(n) {
        // Move constructor with normalization
        if (!c.empty() && c[0]) {
            qNum = c[0]->qNum;
            realqNum = c[0]->qNum;
        } else {
            qNum = 0;
            realqNum = 0;
        }
        qoprof::on_qoperation_create(oplist.size(), &qoprof::stats().qoperation_vector_norm_ctor);
    }
    ~QOperation() {
        qoprof::on_qoperation_destroy(oplist.size());
    }

    bool isZeroValue() const {
        return !isIdentity && oplist.empty();
    }

    bool isZeroSubspace() const {
        return !type && isZeroValue();
    }

    bool isZeroOperator() const {
        return type && isZeroValue();
    }

    const SingleVecTerm* singletonVectorTerm() const {
        if (type || isIdentity || !normalized || oplist.size() != 1) {
            return nullptr;
        }
        return dynamic_cast<const SingleVecTerm*>(oplist[0].get());
    }

    bool singletonSameSpan(const QOperation& other) const {
        const SingleVecTerm* lhs = singletonVectorTerm();
        const SingleVecTerm* rhs = other.singletonVectorTerm();
        if (!lhs || !rhs) {
            return false;
        }
        BIG_COMPLEX_FLOAT overlap = lhs->dot(*rhs);
        double overlap_norm = std::sqrt(double(overlap.real() * overlap.real() + overlap.imag() * overlap.imag()));
        return std::abs(overlap_norm - 1.0) < 1e-8;
    }

    bool singletonOrthogonalTo(const QOperation& other) const {
        const SingleVecTerm* lhs = singletonVectorTerm();
        const SingleVecTerm* rhs = other.singletonVectorTerm();
        if (!lhs || !rhs) {
            return false;
        }
        BIG_COMPLEX_FLOAT overlap = lhs->dot(*rhs);
        double overlap_norm = std::sqrt(double(overlap.real() * overlap.real() + overlap.imag() * overlap.imag()));
        return overlap_norm < 1e-8;
    }

    bool isOperation() const {
        // type == true and isProj == -1
        return type == true && isProj < 0;
    }
    bool isProjection() const {
        // type == false or isProj >= 0
        return type == false || isProj >= 0;
    }

    QOperation& operator=(const QOperation& other) {
        if (this != &other) {
            const size_t old_terms = oplist.size();
            type = other.type;
            normalized = other.normalized;
            isIdentity = other.isIdentity;
            isProj = other.isProj;
            qNum = other.qNum;
            realqNum = other.realqNum;
            oplist.clear();
            oplist.reserve(other.oplist.size());
            for (const auto& term : other.oplist) {
                if (term) {
                    oplist.push_back(term->clone());
                } else {
                    oplist.push_back(nullptr);
                }
            }
            qoprof::on_qoperation_copy_assign(old_terms, oplist.size());
        }
        return *this;
    }
    
    std::string getName() const {
        if(this->oplist.empty()) {
            return "0";
        }
        if(this->isIdentity) {
            return "I";
        }
        if (this->type == false) {
            return this->printFormal(false);
        } else {
            std::string res;
            for (const auto& term : this->oplist) {
                if (term && term->getType()) {
                    auto* gateTerm = dynamic_cast<QuantumGateTerm*>(term.get());
                    if (gateTerm) {
                        // Append the qubit indexes to the name of the gate, for example, "CX(0,1)", "U3(0,1,2)", etc.
                        res += gateTerm->name;
                        // Append gateTerm->index
                        res += "(";
                        for (size_t i = 0; i < gateTerm->index.size(); i++) {
                            res += std::to_string(gateTerm->index[i]);
                            if (i != gateTerm->index.size() - 1) {
                                res += ",";
                            }
                        }
                        res += ")";
                        // Append gateTerm->vars if it is not empty, for example, "U3(0,1,2)[theta,phi,lambda]", "CP(0,1)[theta]", etc.
                        res += "[";
                        for (size_t i = 0; i < gateTerm->vars.size(); i++) {
                            res += std::to_string(gateTerm->vars[i]);
                            if (i != gateTerm->vars.size() - 1) {
                                res += ",";
                            }
                        }
                        res += "]";
                    } else {
                        std::cout << "Unknown quantum term type in getName()." << std::endl;
                    }
                } else {
                    std::cout << "Null quantum term in getName()." << std::endl;
                }
            }
            return res.empty() ? "Emp" : res;
        }
    }

    void append(std::unique_ptr<QuantumTerm> qt) {
        if (qt->getType() != this->type) {
            std::cout << "Append a wrong type of quantum term.";
            return;
        }
        qoprof::on_append();
        oplist.push_back(std::move(qt));
    }
    QOperation add(const QOperation& other) const {
        // std::cout << this->getName() << this->qNum << " + " << other.getName() << other.qNum << std::endl;
        assert(this->qNum == other.qNum);
        // If one of the QOperations is empty, return the other one.
        if (this->oplist.empty()) {
            return other;
        }
        if (other.oplist.empty()) {
            return *this;
        }
        assert(this->type == other.type && this->type == true);
        // Just concatenate the oplist, and mark as not normalized.
        QOperation res = *this;
        size_t cloned_terms = 0;
        for (const auto& term : other.oplist) {
            if (term) {
                res.oplist.push_back(term->clone());
                ++cloned_terms;
            } else {
                res.oplist.push_back(nullptr);
            }
        }
        qoprof::on_add(cloned_terms);
        return res;
    }
    std::vector<std::unique_ptr<QuantumTerm>> fetch(unsigned int begin,unsigned int end) {
        // fetch the oplist from begin to end, using move semantics
        assert(end < static_cast<int>(this->oplist.size()));
        std::vector<std::unique_ptr<QuantumTerm>> res;
        res.reserve(end - begin + 1);
        size_t cloned_terms = 0;
        for (size_t i = begin; i <= end; i++) {
            if (this->oplist[i] == nullptr) {
                res.push_back(nullptr);
            } else {
                res.push_back(this->oplist[i]->clone()); // Use clone to ensure deep copy
                ++cloned_terms;
            }
        }
        qoprof::on_fetch(cloned_terms);
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
        qoprof::on_gen_proj_meas_space(this->qNum > 0 ? std::pow(2, this->qNum - 1) : 0);
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
        const size_t input_terms = this->oplist.size();
        if (this->oplist.size() == 0) {
            qoprof::on_gram_schmidt(input_terms, input_terms);
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
        qoprof::on_gram_schmidt(input_terms, this->oplist.size());
    }

    QOperation minus(const QOperation& other) const {
        qoprof::on_minus_call();
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
        qoprof::on_conjunction_call();
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
        qoprof::on_conjunction_simp_call();
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
        qoprof::on_disjunction_call();
        /*** To do the disjunction:
         * 1. preserve the abstract semantic tree;
         * 2. widening function;
         * 3. hard calculate \neg(\neg A \land \neg B).
         * Already Schmidted.
         ***/
        assert(!this->type && !other.type);
        if (other.oplist.size() == 0) {
            // The vector remains unnormalized here! And without judgement of zero vector!!!
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
        if (this->singletonSameSpan(other)) {
            return *this;
        }
        if (this->singletonOrthogonalTo(other)) {
            QOperation res(*this, other);
            res.normalized = true;
            res.qNum = this->qNum;
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
        qoprof::on_preimage_call();
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
        if (other.isIdentity) {
            res = *this;
            return res;
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
        const size_t input_terms = this->oplist.size();
        const size_t gate_terms = other.oplist.size();
        /* The post-image of a quantum operator */
        assert(this->type == false && other.type == true);
        QOperation res;
        for (size_t i = 0; i < this->oplist.size(); i++) {
            auto* ivec = dynamic_cast<SingleVecTerm*>(this->oplist[i].get());
            if (!ivec) continue;
            if (checkifzero(ivec->content)) continue; // A key optimization!
            for (size_t j = 0; j < other.oplist.size(); j++) {
                auto* jgate = dynamic_cast<QuantumGateTerm*>(other.oplist[j].get());
                if (!jgate) continue;
                auto tmp = ivec->applyGate(*jgate, true);
                // If tmp is non-trivial, we insert it to the res.
                if (!checkifzero(tmp)) {
                    res.append(std::make_unique<SingleVecTerm>(SingleVecTerm(tmp)));
                } else {
                    // std::cout << "PostImage: a zero vector is generated." << jgate->name << std::endl;
                    // this->printFormal();
                }
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
        // res.normalized = true;
        res.qNum = this->qNum;
        const size_t raw_terms = res.oplist.size();
        // Use GramSchmidt to handle the uniqueness of the vectors. And remove zero vectors.
        if (res.oplist.size() > 1) {
            res.GramSchmidt(0, static_cast<int>(res.oplist.size())-1);
            // res.normalized = true;
        } else {
            // normalize the single vector if there is only one vector.
            if (res.oplist.size() == 1) {
                auto* ivec = dynamic_cast<SingleVecTerm*>(res.oplist[0].get());
                if (!ivec) {
                    std::cout << "Strange nullptr" << std::endl;
                } else {
                    ivec->normalizeInline();
                }
            }
            res.normalized = true;
        }
        qoprof::on_postimage_call(input_terms, gate_terms, raw_terms, res.oplist.size());
        return res;
    }
    
    int compare(const QOperation& other) const {
        /* Compare the relation of this and other */
        /* 0: this is included in other; 1: other is included in this; 2: exclude; 3: intersect but not include; 4: equality*/
        // std::cout << this->normalized << " " << other.normalized << std::endl;
        assert(this->normalized && other.normalized);
        assert(this->qNum == other.qNum);
        if (this->oplist.size() == 0 && other.oplist.size() == 0) {
            return 4;
        }
        if (this->isIdentity && other.isIdentity) {
            return 4; // Both are identity operators
        }
        if (this->isIdentity) {
            return 1;
        }
        if (other.isIdentity) {
            return 0;
        }
        if (this->isProj >=0 && other.isProj >=0) {
            
        } else if (this->isProj >= 0 && other.isProj < 0) {
            // this is a projective operator, other is a subspace
        } else if (this->isProj < 0 && other.isProj >= 0) {
            // this is a subspace, other is a projective operator
        }
        if (this->oplist.size() > 0) {
            assert(this->oplist[0]->getType() == false);
        } else {
            // If other is not empty, return 0; else return 4
            if (other.oplist.size() > 0) {
                return 0; // The empty operator is a subspace of any operator.
            } else {
                return 4; // The empty operator is equal to the empty operator.
            }
        }
        if (other.oplist.size() > 0) {
            assert(other.oplist[0]->getType() == false);
        } else {
            return 1; // Any non zero operator is a super-space of the empty operator.
        }
        if (this->singletonSameSpan(other)) {
            return 4;
        }
        if (this->singletonOrthogonalTo(other)) {
            return 2;
        }
        // In case this is support-vector like subspace
        // Compute the disjunction of this and other, if the result dimension is:
        // 1. the same as both, then this is equal to other, return 4;
        // 2. the same as other, then this is included in other, return 0;
        // 3. the same as this, then other is included in this, return 1;
        // 4. the same as the sum of the dimensions of this and other, then this is disjoint with other, return 2;
        // 5. otherwise, this is intersecting with other but not included, return 3.
        QOperation disj = this->disjunction(other);
        if (disj.oplist.size() == this->oplist.size() && disj.oplist.size() == other.oplist.size()) {
            return 4;
        }
        if (disj.oplist.size() == other.oplist.size()) {
            return 0;
        }
        if (disj.oplist.size() == this->oplist.size()) {
            return 1;
        }
        if (disj.oplist.size() == this->oplist.size() + other.oplist.size()) {
            return 2;
        }
        return 3;
    }
    bool operator==(const QOperation& other) const {
        if (this->qNum != other.qNum ||
            this->realqNum != other.realqNum ||
            this->isIdentity != other.isIdentity ||
            this->isProj != other.isProj) {
            return false;
        }
        if(this->isOperation() != other.isOperation()) {
            return false;
        }
        // Either both are projective measurements or both are gates.
        if(this->type && other.type) {
            // We assume the order of the gates in oplist is consistent, which is guaranteed by the constructor of QOperation.
            if (this->oplist.size() != other.oplist.size()) {
                return false;
            }
            if (this->oplist.size() == 0) {
                return true;
            }
            for (size_t i = 0; i < this->oplist.size(); i++) {
                auto* gate1 = dynamic_cast<QuantumGateTerm*>(this->oplist[i].get());
                auto* gate2 = dynamic_cast<QuantumGateTerm*>(other.oplist[i].get());
                if (!gate1 || !gate2) {
                    return false;
                }
                if (!gate1->isEqual(*gate2)) {
                    return false;
                }
            }
            return true;
        }
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
        if (this->singletonSameSpan(other)) {
            return true;
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

    std::string printFormal(bool print=true) const {
        std::string finalresult;
        if(print) {
            std::cout << "Support vectors:" << std::endl;
        }
        if (this->isIdentity) {
            if(print) {
                std::cout << "Identity operator." << std::endl;
            }
            return "I";
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
            if (print) {
                std::cout << result.str() << std::endl;
            }
            finalresult += result.str() + "\n";
        }
        // If finalresult is empty, return "0".
        if (finalresult.empty()) {
            finalresult = "0";
        }
        return finalresult;
    }
};

inline QOperation SpanQOperations(const std::vector<QOperation>& ops) {
    if (ops.empty()) {
        throw std::runtime_error("Cannot construct the span of an empty QOperation list.");
    }

    QOperation res(false);
    bool initialized = false;

    for (const auto& op : ops) {
        if (op.isIdentity) {
            QOperation identity(false);
            identity.normalized = true;
            identity.qNum = op.qNum;
            identity.realqNum = op.realqNum;
            identity.isIdentity = true;
            return identity;
        }
        if (op.type) {
            throw std::runtime_error("SpanQOperations expects subspace/state QOperations, not gate-type QOperations.");
        }
        if (op.oplist.empty()) {
            continue;
        }
        if (!initialized) {
            res.qNum = op.qNum;
            res.realqNum = op.realqNum;
            initialized = true;
        } else if (res.qNum != op.qNum) {
            throw std::runtime_error("Cannot span QOperations with different qNum values.");
        }

        QOperation normalizedOp(op);
        if (!normalizedOp.normalized && !normalizedOp.oplist.empty()) {
            normalizedOp.GramSchmidt(0, static_cast<int>(normalizedOp.oplist.size()) - 1);
        }
        for (const auto& term : normalizedOp.oplist) {
            if (term) {
                res.oplist.push_back(term->clone());
            }
        }
    }

    if (!initialized || res.oplist.empty()) {
        res.normalized = true;
        return res;
    }
    if (res.qNum == 1 && res.oplist.size() > 1) {
        // The current CFLOBDD matrix-multiply path used by GramSchmidt/dot asserts
        // on the one-qubit level-1 case.  Any two distinct one-qubit states span
        // the full one-qubit space; callers that build from strings should dedup
        // identical inputs before reaching this lower-level helper.
        return QOperation(std::vector<std::string>{"0", "1"});
    }
    if (res.oplist.size() > 1) {
        res.GramSchmidt(0, static_cast<int>(res.oplist.size()) - 1);
    } else {
        // Single-vector inputs from QOperation string/amplitude constructors are already
        // normalized.  Avoid normalizeInline() here because the current CFLOBDD
        // MatrixMultiplyV4WithInfo path asserts for the one-qubit level-1 case.
        res.normalized = true;
    }
    return res;
}

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

QOperation CreateZeroQO(unsigned int qNum, bool optype = false) {
    unsigned int logicqNum = std::pow(2, ceil(log2(qNum)));
    QOperation res(optype);
    res.realqNum = qNum;
    res.qNum = logicqNum;
    res.normalized = true;
    // res.isProj = -1; // No need to set isProj for zero operator.
    return res;
}

// Should be mearged into the constructor of QOperation
QOperation CreateProjectiveMeasQO(unsigned int qNum, unsigned int i, bool val) {
    unsigned int logicqNum = std::pow(2, ceil(log2(qNum)));
    QOperation res(true);
    res.qNum = logicqNum;
    res.realqNum = qNum;
    res.isProj = i;
    res.normalized = true;
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
