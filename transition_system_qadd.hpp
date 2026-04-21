#ifndef TRANSITION_SYSTEM_QADD_HPP
#define TRANSITION_SYSTEM_QADD_HPP

#include <vector>
#include <cmath>
#include <cassert>

#include "qadd.hpp"
#include "quantum_operation.hpp"
#include "cl_proposition.hpp"

namespace qts {

using namespace qadd;

// ========================================
// TransitionSystem (QADD-based)
// ========================================

int MAX_NUM_VARS = 20; // Maximum number of bits for location encoding

void initializeTransitionSystem() {
    CFLOBDDNodeHandle::InitNoDistinctionTable();
    CFLOBDDNodeHandle::InitAdditionInterleavedTable();
    CFLOBDDNodeHandle::InitReduceCache();
    InitPairProductCache();
    InitTripleProductCache();
    Matrix1234ComplexFloatBoost::Matrix1234Initializer();
    VectorComplexFloatBoost::VectorInitializer();
}

class TransitionSystem {
public:
    TransitionSystem(int num_qubits, int max_locations = 0) : num_locations(0), num_qubits(num_qubits) {
        if(max_locations > 0) {
            num_vars = std::ceil(std::log2(max_locations));
        } else {
            num_vars = MAX_NUM_VARS;
        }
        annotation = make_terminal(CreateZeroQO(num_qubits, false));
        relation   = make_terminal(CreateZeroQO(num_qubits, true));
    }

    // =============================
    // Location Encoding
    // =============================

    int addLocation() {
        int id = num_locations++;
        update_num_vars();

        encodings.push_back(encode(id));
        return id;
    }

    // =============================
    // Annotation
    // =============================

    void setAnnotation(int loc, const QOperation& val) {
        assert(loc < num_locations);

        QADDNode* indicator = build_state_indicator(0, encodings[loc], val);

        annotation = Apply(ApplyOp::JOIN, annotation, indicator);
    }

    // =============================
    // Relation
    // =============================

    void addRelation(int src, int dst, const QOperation& op) {
        assert(src < num_locations && dst < num_locations);

        QADDNode* indicator =
            build_relation_indicator(0,
                                     encodings[src],
                                     encodings[dst],
                                     op);

        relation = Apply(ApplyOp::ADD, relation, indicator);
    }

    // =============================
    // Satisfaction of classical propositions
    // =============================
    std::vector<std::string> satisfyBit(int loc, std::vector<unsigned int> indexs, std::vector<bool> values) {
        assert(loc < num_locations);
        return cp[loc].satisfyBit(indexs, values);
    }

    std::vector<std::string> unsatisfyBit(int loc, std::vector<unsigned int> indexs, std::vector<bool> values) {
        assert(loc < num_locations);
        return cp[loc].unsatisfyBit(indexs, values);
    }

    // =============================
    // Model Checking
    // =============================
    void postConditions()
    {
        // Step 1: Apply (mixed)
        QADDNode* tmp = Apply(ApplyOp::APPLY, relation, annotation);

        // Step 2: existential abstraction over x (0..n-1)
        QADDNode* eliminated = exists_vars(tmp);

        // Step 3: rename x' -> x
        QADDNode* next = rename_vars(eliminated);

        annotation = Apply(ApplyOp::JOIN, annotation, next);
    }

    // =============================
    // Getters
    // =============================

    QADDNode* getAnnotation() const { return annotation; }
    QADDNode* getRelation() const { return relation; }

    int getNumVars() const { return num_vars; }
    int getNumLocations() const { return num_locations; }
    int getNumQubits() const { return num_qubits; }

    // =============================
    // Visualization
    // =============================
    void printAnnotation(const std::string& filename="annotation.dot") const {
        std::cout << "Annotation:" << std::endl;
        printQADD(annotation, filename);
    }

    void printRelation(const std::string& filename="relation.dot") const {
        std::cout << "Relation:" << std::endl;
        printQADD(relation, filename);
    }

public:

    // =============================
    // Internal State
    // =============================

    int num_vars;          // number of bits
    int num_locations;
    int num_qubits;        // number of qubits in the quantum operations

    std::vector<std::vector<bool>> encodings;

    std::vector<ClassicalProposition> cp;

    QADDNode* annotation;
    QADDNode* relation;

    // =============================
    // Encoding Helpers
    // =============================

    void update_num_vars() {
        int needed = std::ceil(std::log2(std::max(1, num_locations)));
        if (needed > num_vars) {
            num_vars = needed;
        }
    }

    std::vector<bool> encode(int id) {
        std::vector<bool> bits(num_vars, false);
        for (int i = 0; i < num_vars; ++i) {
            bits[i] = (id >> i) & 1;
        }
        return bits;
    }

    // =============================
    // Indicator Builders
    // =============================

    // build: χ_loc(x)
    QADDNode* build_state_indicator(
        int i,
        const std::vector<bool>& enc,
        const QOperation& val)
    {
        if (i == num_vars) {
            return make_terminal(val);
        }

        QADDNode* zero = make_terminal(CreateZeroQO(num_qubits));

        if (enc[i]) {
            return make_node(i,
                zero,
                build_state_indicator(i+1, enc, val));
        } else {
            return make_node(i,
                build_state_indicator(i+1, enc, val),
                zero);
        }
    }

    // build: χ_{(src,dst)}(x,x')
    QADDNode* build_relation_indicator(
        int i,
        const std::vector<bool>& src,
        const std::vector<bool>& dst,
        const QOperation& val)
    {
        if (i == 2 * num_vars) {
            return make_terminal(val);
        }

        QADDNode* zero = make_terminal(CreateZeroQO(num_qubits, true));

        bool bit;
        if (i < num_vars) {
            bit = src[i];
        } else {
            bit = dst[i - num_vars];
        }

        if (bit) {
            return make_node(i,
                zero,
                build_relation_indicator(i+1, src, dst, val));
        } else {
            return make_node(i,
                build_relation_indicator(i+1, src, dst, val),
                zero);
        }
    }

    // =============================
    // Model Checking Helpers
    // =============================
        QADDNode* restrict_var(QADDNode* node, int var, bool value) {
        if (is_terminal(node)) return node;

        if (node->var == var) {
            return value ? node->high : node->low;
        }

        QADDNode* low = restrict_var(node->low, var, value);
        QADDNode* high = restrict_var(node->high, var, value);

        return make_node(node->var, low, high);
    }

    QADDNode* exists_var(QADDNode* node, int var) {
        QADDNode* f0 = restrict_var(node, var, false);
        QADDNode* f1 = restrict_var(node, var, true);

        return Apply(ApplyOp::JOIN, f0, f1);
    }

    QADDNode* exists_vars(QADDNode* node) {
        for (int i = 0; i < num_vars; ++i) {
            node = exists_var(node, i);
        }
        return node;
    }

    QADDNode* rename_vars(QADDNode* node) {
        if (is_terminal(node)) return node;

        QADDNode* low = rename_vars(node->low);
        QADDNode* high = rename_vars(node->high);

        return make_node(node->var - num_vars, low, high);
    }
};

} // namespace qts

#endif