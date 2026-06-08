#ifndef TRANSITION_SYSTEM_QADD_HPP
#define TRANSITION_SYSTEM_QADD_HPP

#include <vector>
#include <cmath>
#include <cassert>
#include <algorithm>
#include <unordered_map>
#include <unordered_set>
#include <iostream>
#include <cstdint>
#include <chrono>
#include <cstdlib>

#include "qadd.hpp"
#include "quantum_operation.hpp"
#include "cl_proposition.hpp"

namespace qts {

using namespace qadd;

namespace tsprof {

struct Stats {
    size_t post_conditions_calls = 0;
    size_t post_iterations = 0;
    size_t max_living_count = 0;
    size_t max_post_count = 0;
    double extract_delta_ms = 0.0;
    double apply_ms = 0.0;
    double exists_ms = 0.0;
    double rename_ms = 0.0;
    double join_ms = 0.0;
    double dim_before_ms = 0.0;
    double dim_after_ms = 0.0;
    double get_post_ids_ms = 0.0;
};

inline bool enabled() {
    static const bool enabled_flag = []() {
        const char* value = std::getenv("TS_PROFILE");
        return value && std::string(value) != "0";
    }();
    return enabled_flag;
}

inline Stats& stats() {
    static Stats value;
    return value;
}

inline double elapsed_ms(std::chrono::steady_clock::time_point start,
                         std::chrono::steady_clock::time_point end) {
    return std::chrono::duration_cast<std::chrono::duration<double, std::milli>>(end - start).count();
}

inline void print_summary(std::ostream& os = std::cout) {
    if (!enabled()) return;
    const Stats& s = stats();
    os << "[tsprof] post_conditions_calls=" << s.post_conditions_calls
       << " post_iterations=" << s.post_iterations
       << " max_living=" << s.max_living_count
       << " max_postids=" << s.max_post_count
       << std::endl;
    os << "[tsprof] get_post_ids_ms=" << s.get_post_ids_ms
       << " dim_before_ms=" << s.dim_before_ms
       << " extract_delta_ms=" << s.extract_delta_ms
       << " apply_ms=" << s.apply_ms
       << " exists_ms=" << s.exists_ms
       << " rename_ms=" << s.rename_ms
       << " join_ms=" << s.join_ms
       << " dim_after_ms=" << s.dim_after_ms
       << std::endl;
}

inline int iteration_limit_override() {
    static const int value = []() {
        const char* raw = std::getenv("TS_MAX_POST_ITER");
        if (!raw || std::string(raw).empty()) {
            return -1;
        }
        return std::atoi(raw);
    }();
    return value;
}

} // namespace tsprof

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
    struct SymbolicPostPairKey {
        QADDNode* relation_node;
        QADDNode* delta_node;

        bool operator==(const SymbolicPostPairKey& other) const {
            return relation_node == other.relation_node && delta_node == other.delta_node;
        }
    };

    struct SymbolicPostPairKeyHash {
        size_t operator()(const SymbolicPostPairKey& key) const {
            return (std::hash<QADDNode*>()(key.relation_node) << 1) ^
                   std::hash<QADDNode*>()(key.delta_node);
        }
    };

    TransitionSystem() : num_locations(0), num_qubits(0), num_vars(MAX_NUM_VARS) {
        annotation = make_terminal(CreateZeroQO(num_qubits, false));
        relation   = make_terminal(CreateZeroQO(num_qubits, true));
    }
    
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
        return addLocation(ClassicalProposition(), "");
    }

    int addLocation(const ClassicalProposition& proposition, const std::string& identifier = "") {
        int id = num_locations++;

        bool grew = update_num_vars();
        if (grew) {
            extend_existing_encodings();
        }

        encodings.push_back(encode(id));
        cp.push_back(proposition);
        identifiers.push_back(identifier);
        labels.emplace_back();
        post_locations.emplace_back();
        post_relation_terminals.emplace_back();
        post_marks.push_back(0);
        return id;
    }

    void setInitLocation(int loc) {
        assert(loc >= 0 && loc < num_locations);
        init_location = loc;
    }

    int getInitLocation() const {
        return init_location;
    }

    std::vector<int> getLocationIDs() const {
        std::vector<int> ids(num_locations);
        for (int id = 0; id < num_locations; ++id) {
            ids[id] = id;
        }
        return ids;
    }

    // =============================
    // Annotation
    // =============================

    void setAnnotation(int loc, const QOperation& val) {
        assert(loc < num_locations);

        QADDNode* indicator = build_state_indicator(0, encodings[loc], val);

        annotation = Apply(ApplyOp::JOIN, annotation, indicator);
        append_unique(initAnnotationID, loc);
        append_unique(livingAnnotationID, loc);
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

        auto& outgoing_locations = post_locations[src];
        auto& outgoing_terms = post_relation_terminals[src];
        auto existing = std::find(outgoing_locations.begin(), outgoing_locations.end(), dst);
        if (existing == outgoing_locations.end()) {
            outgoing_locations.push_back(dst);
            outgoing_terms.push_back(make_terminal(op));
        } else {
            size_t index = static_cast<size_t>(std::distance(outgoing_locations.begin(), existing));
            outgoing_terms[index] = apply_terminal(ApplyOp::ADD, outgoing_terms[index], make_terminal(op));
        }
    }

    std::vector<int> getPostIDs(const std::vector<int>& ids) {
        std::vector<int> postids;
        postids.reserve(ids.size() * 2);

        ++post_mark_epoch;
        if (post_mark_epoch == 0) {
            std::fill(post_marks.begin(), post_marks.end(), 0);
            post_mark_epoch = 1;
        }

        for (int id : ids) {
            assert(id >= 0 && id < num_locations);
            for (int dst : post_locations[id]) {
                if (post_marks[dst] != post_mark_epoch) {
                    post_marks[dst] = post_mark_epoch;
                    postids.push_back(dst);
                }
            }
        }
        return postids;
    }

    // Extract a delta P-DD containing only the annotations for the given locations.
    QADDNode* extract_delta(const std::vector<int>& ids) {
        QADDNode* delta = make_terminal(CreateZeroQO(num_qubits, false));
        for (int id : ids) {
            assert(id >= 0 && id < num_locations);
            QADDNode* term = get_location_terminal(annotation, id);
            QOperation val = term->val; // copy terminal QOperation
            QADDNode* indicator = build_state_indicator(0, encodings[id], val);
            delta = Apply(ApplyOp::JOIN, delta, indicator);
        }
        return delta;
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

    void appendClassicalAP(int loc, const std::string& ap) {
        assert(loc >= 0 && loc < num_locations);
        cp[loc].addTerm(ap);
    }

    void setClassicalValue(int loc, unsigned int index, bool value) {
        assert(loc >= 0 && loc < num_locations);
        cp[loc].setValue(index, value);
    }

    bool find(int loc, const std::string& ap) const {
        assert(loc >= 0 && loc < num_locations);
        return cp[loc].find(ap);
    }

    int termNum(int loc) const {
        assert(loc >= 0 && loc < num_locations);
        return static_cast<int>(cp[loc].terms.size());
    }

    ClassicalProposition getClassicalProposition(int loc) const {
        assert(loc >= 0 && loc < num_locations);
        return cp[loc];
    }

    void setClassicalProposition(int loc, const ClassicalProposition& proposition) {
        assert(loc >= 0 && loc < num_locations);
        cp[loc] = proposition;
    }

    void setIdentifier(int loc, const std::string& identifier) {
        assert(loc >= 0 && loc < num_locations);
        identifiers[loc] = identifier;
    }

    std::string getIdentifier(int loc) const {
        assert(loc >= 0 && loc < num_locations);
        return identifiers[loc];
    }

    void setLabel(int loc, const std::string& label) {
        assert(loc >= 0 && loc < num_locations);
        auto& loc_labels = labels[loc];
        if (std::find(loc_labels.begin(), loc_labels.end(), label) == loc_labels.end()) {
            loc_labels.push_back(label);
        }
    }

    std::vector<std::string> getLabels(int loc) const {
        assert(loc >= 0 && loc < num_locations);
        return labels[loc];
    }

    // =============================
    // Model Checking
    // =============================
    void postOneStep() {
        // Step 1: Apply (mixed)
        QADDNode* tmp = Apply(ApplyOp::APPLY, relation, annotation);

        // APPLY cache can be dropped before existential elimination.
        clear_compute_table();

        // Step 2: existential abstraction over source vars x_i (2*i)
        QADDNode* eliminated = exists_vars(tmp);

        // Step 3: rename target vars x_i' (2*i+1) -> source vars x_i (2*i)
        QADDNode* next = rename_vars(eliminated);

        clear_compute_table();
        annotation = Apply(ApplyOp::JOIN, annotation, next);

        // Apply cache is temporary; clearing it avoids unbounded growth across iterations.
        clear_compute_table();
    }

    // Perform one step but only for a delta P-DD built from `ids`.
    void postOneStepDelta(const std::vector<int>& ids) {
        if (ids.empty()) return;

        const bool profile = tsprof::enabled();
        auto step_start = std::chrono::steady_clock::now();
        QADDNode* delta = extract_delta(ids);
        if (profile) {
            auto after_extract = std::chrono::steady_clock::now();
            tsprof::stats().extract_delta_ms += tsprof::elapsed_ms(step_start, after_extract);
            step_start = after_extract;
        }

        // Step 1: Apply relation to delta
        QADDNode* tmp = Apply(ApplyOp::APPLY, relation, delta);
        if (profile) {
            auto after_apply = std::chrono::steady_clock::now();
            tsprof::stats().apply_ms += tsprof::elapsed_ms(step_start, after_apply);
            step_start = after_apply;
        }
        clear_compute_table();

        // Step 2: existential abstraction over source vars x_i (2*i)
        QADDNode* eliminated = exists_vars(tmp);
        if (profile) {
            auto after_exists = std::chrono::steady_clock::now();
            tsprof::stats().exists_ms += tsprof::elapsed_ms(step_start, after_exists);
            step_start = after_exists;
        }
        clear_compute_table();

        // Step 3: rename target vars x_i' (2*i+1) -> source vars x_i (2*i)
        QADDNode* next = rename_vars(eliminated);
        if (profile) {
            auto after_rename = std::chrono::steady_clock::now();
            tsprof::stats().rename_ms += tsprof::elapsed_ms(step_start, after_rename);
            step_start = after_rename;
        }
        clear_compute_table();

        // Join incremental next into global annotation
        annotation = Apply(ApplyOp::JOIN, annotation, next);
        if (profile) {
            auto after_join = std::chrono::steady_clock::now();
            tsprof::stats().join_ms += tsprof::elapsed_ms(step_start, after_join);
        }
        clear_compute_table();
    }

    QADDNode* compute_symbolic_post(QADDNode* delta) {
        const bool profile = tsprof::enabled();
        auto step_start = std::chrono::steady_clock::now();

        std::vector<QOperation> next_values(
            static_cast<size_t>(num_locations),
            CreateZeroQO(num_qubits, false));
        std::vector<unsigned int> touched_marks(static_cast<size_t>(num_locations), 0);
        std::vector<int> touched_locations;
        unsigned int touched_epoch = 1;

        for (int src = 0; src < num_locations; ++src) {
            const QOperation& source_value = get_location_terminal(delta, src)->val;
            if (source_value.isZeroSubspace()) {
                continue;
            }

            const auto& outgoing_locations = post_locations[src];
            const auto& outgoing_terms = post_relation_terminals[src];
            for (size_t edge_index = 0; edge_index < outgoing_locations.size(); ++edge_index) {
                int dst = outgoing_locations[edge_index];
                const QOperation& relation_value = outgoing_terms[edge_index]->val;
                QOperation image = source_value.postImage(relation_value);
                if (image.isZeroSubspace()) {
                    continue;
                }

                QOperation& slot = next_values[static_cast<size_t>(dst)];
                if (slot.isZeroSubspace()) {
                    slot = image;
                } else {
                    slot = slot.disjunction(image);
                }

                if (touched_marks[static_cast<size_t>(dst)] != touched_epoch) {
                    touched_marks[static_cast<size_t>(dst)] = touched_epoch;
                    touched_locations.push_back(dst);
                }
            }
        }

        QADDNode* next = make_terminal(CreateZeroQO(num_qubits, false));
        for (int dst : touched_locations) {
            const QOperation& value = next_values[static_cast<size_t>(dst)];
            if (value.isZeroSubspace()) {
                continue;
            }
            QADDNode* indicator = build_state_indicator(0, encodings[dst], value);
            next = Apply(ApplyOp::JOIN, next, indicator);
        }

        if (profile) {
            auto after_apply = std::chrono::steady_clock::now();
            tsprof::stats().apply_ms += tsprof::elapsed_ms(step_start, after_apply);
            step_start = after_apply;
        }
        return next;
    }

    bool is_zero_delta(QADDNode* node) const {
        return is_terminal(node) && !node->val.isIdentity && get_dimension(node->val) == 0;
    }

    void postConditions()
    {
        int MAX_ITER = (1 << num_qubits) * 2; // upper bound of iterations to prevent infinite loop
        int iter = 0;
        std::cout << MAX_ITER << " maximum iterations allowed. " << num_locations << " " << num_qubits << std::endl;
        if (tsprof::enabled()) {
            ++tsprof::stats().post_conditions_calls;
        }

        QADDNode* delta = annotation;

        while (true) {
            if (iter++ > MAX_ITER) {
                std::cerr << "Warning: postConditions reached maximum iterations. Possible non-convergence." << std::endl;
                break;
            }
            if (tsprof::iteration_limit_override() >= 0 && iter > tsprof::iteration_limit_override()) {
                std::cerr << "Warning: postConditions stopped early due to TS_MAX_POST_ITER="
                          << tsprof::iteration_limit_override() << std::endl;
                break;
            }

            if (is_zero_delta(delta)) {
                break;
            }

            if (tsprof::enabled()) {
                ++tsprof::stats().post_iterations;
            }

            QADDNode* old_annotation = annotation;
            QADDNode* next = compute_symbolic_post(delta);

            auto stage_start = std::chrono::steady_clock::now();
            QADDNode* joined = Apply(ApplyOp::JOIN, annotation, next);
            if (tsprof::enabled()) {
                auto stage_end = std::chrono::steady_clock::now();
                tsprof::stats().join_ms += tsprof::elapsed_ms(stage_start, stage_end);
                stage_start = stage_end;
            }

            QADDNode* new_delta = Apply(ApplyOp::DIFF, next, old_annotation);
            if (tsprof::enabled()) {
                auto stage_end = std::chrono::steady_clock::now();
                tsprof::stats().extract_delta_ms += tsprof::elapsed_ms(stage_start, stage_end);
            }
            clear_compute_table();

            annotation = joined;
            delta = new_delta;
        }

        tsprof::print_summary();
    }

    // =============================
    // Getters
    // =============================

    QADDNode* getAnnotation() const { return annotation; }
    QADDNode* getRelation() const { return relation; }

    QOperation getLocationAnnotation(int loc) const {
        return get_location_terminal(annotation, loc)->val;
    }

    int getLocationDimension(int loc) const {
        return get_dimension(getLocationAnnotation(loc));
    }

    bool locationHasNonZeroAnnotation(int loc) const {
        QOperation value = getLocationAnnotation(loc);
        if (value.isIdentity) {
            return true;
        }
        return get_dimension(value) > 0;
    }

    std::pair<int, int> printDims(int loc) const {
        assert(loc >= 0 && loc < num_locations);
        int dim = getLocationDimension(loc);
        return std::make_pair(dim, dim);
    }

    std::vector<int> getPostLocations(int loc) const {
        assert(loc >= 0 && loc < num_locations);
        return post_locations[loc];
    }

    std::string getRelationName(int src, int dst) const {
        assert(src >= 0 && src < num_locations);
        assert(dst >= 0 && dst < num_locations);
        return get_relation_terminal(src, dst)->val.getName();
    }

    std::vector<int> filterReachableLocations() const {
        std::vector<int> reachable;
        for (int loc = 0; loc < num_locations; ++loc) {
            if (locationHasNonZeroAnnotation(loc)) {
                reachable.push_back(loc);
            }
        }
        return reachable;
    }

    bool satisfy(int loc, const QOperation& spec) const {
        assert(loc >= 0 && loc < num_locations);
        const QOperation actual = getLocationAnnotation(loc);
        const int relation = actual.compare(spec);
        return relation == 0 || relation == 4;
    }

    bool isLeafLoc(int loc) const {
        assert(loc >= 0 && loc < num_locations);
        return post_locations[loc].empty() ||
               (post_locations[loc].size() == 1 && post_locations[loc][0] == loc);
    }

    int getNumVars() const { return num_vars; }
    int getNumLocations() const { return num_locations; }
    int getNumQubits() const { return num_qubits; }

    size_t getAnnotationNodeCount() const { return count_nodes(annotation); }
    size_t getRelationNodeCount() const { return count_nodes(relation); }
    size_t getTotalUniqueNodeCount() const { return uniqueTable.size() + terminalTable.size(); }

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
    int init_location = -1;

    std::vector<std::vector<bool>> encodings;

    std::vector<ClassicalProposition> cp;
    std::vector<std::string> identifiers;
    std::vector<std::vector<std::string>> labels;

    QADDNode* annotation;
    QADDNode* relation;

    std::vector<int> initAnnotationID;
    std::vector<int> livingAnnotationID;
    std::vector<std::vector<int>> post_locations;
    std::vector<std::vector<QADDNode*>> post_relation_terminals;
    std::vector<unsigned int> post_marks;
    unsigned int post_mark_epoch = 0;

    // =============================
    // Encoding Helpers
    // =============================

    bool update_num_vars() {
        int needed = std::ceil(std::log2(std::max(1, num_locations)));
        if (needed > num_vars) {
            num_vars = needed;
            return true;
        }
        return false;
    }

    void extend_existing_encodings() {
        for (auto& bits : encodings) {
            bits.resize(num_vars, false);
        }
    }

    std::vector<bool> encode(int id) {
        std::vector<bool> bits(num_vars, false);
        for (int i = 0; i < num_vars; ++i) {
            bits[i] = (id >> i) & 1;
        }
        return bits;
    }

    int src_var(int bit) const { return 2 * bit; }
    int dst_var(int bit) const { return 2 * bit + 1; }

    void append_unique(std::vector<int>& vec, int id) {
        if (std::find(vec.begin(), vec.end(), id) == vec.end()) {
            vec.push_back(id);
        }
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
        int var = src_var(i);

        if (enc[i]) {
            return make_node(var,
                zero,
                build_state_indicator(i+1, enc, val));
        } else {
            return make_node(var,
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

        int bit_idx = i / 2;
        bool bit;
        if ((i % 2) == 0) {
            bit = src[bit_idx];
        } else {
            bit = dst[bit_idx];
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
    QADDNode* get_location_terminal(QADDNode* node, int loc) const {
        assert(loc >= 0 && loc < num_locations);

        QADDNode* cur = node;
        const std::vector<bool>& enc = encodings[loc];

        while (!is_terminal(cur)) {
            int bit_idx = cur->var / 2;
            assert(bit_idx >= 0 && bit_idx < num_vars);
            bool bit = enc[bit_idx];
            cur = bit ? cur->high : cur->low;
        }
        return cur;
    }

    QADDNode* get_relation_terminal(int src, int dst) const {
        assert(src >= 0 && src < num_locations);
        assert(dst >= 0 && dst < num_locations);

        QADDNode* cur = relation;
        const std::vector<bool>& src_enc = encodings[src];
        const std::vector<bool>& dst_enc = encodings[dst];

        while (!is_terminal(cur)) {
            int bit_idx = cur->var / 2;
            assert(bit_idx >= 0 && bit_idx < num_vars);
            bool bit = (cur->var % 2 == 0) ? src_enc[bit_idx] : dst_enc[bit_idx];
            cur = bit ? cur->high : cur->low;
        }
        return cur;
    }

    size_t count_nodes(QADDNode* root) const {
        if (!root) return 0;

        std::unordered_set<QADDNode*> visited;
        std::vector<QADDNode*> stack{root};
        while (!stack.empty()) {
            QADDNode* node = stack.back();
            stack.pop_back();
            if (!visited.insert(node).second || is_terminal(node)) {
                continue;
            }
            stack.push_back(node->low);
            stack.push_back(node->high);
        }
        return visited.size();
    }

    std::vector<int> collect_dimensions(QADDNode* node, const std::vector<int>& ids) const {
        std::vector<int> dims;
        dims.reserve(ids.size());

        for (int id : ids) {
            QADDNode* terminal = get_location_terminal(node, id);
            dims.push_back(get_dimension(terminal->val));
        }
        return dims;
    }

    QADDNode* restrict_var_impl(
        QADDNode* node,
        int var,
        bool value,
        std::unordered_map<QADDNode*, QADDNode*>& memo)
    {
        auto it = memo.find(node);
        if (it != memo.end()) {
            return it->second;
        }

        if (is_terminal(node)) return node;

        // Variable ordering property: if current variable already exceeds target,
        // the target variable cannot appear in this subtree.
        if (node->var > var) {
            memo[node] = node;
            return node;
        }

        if (node->var == var) {
            QADDNode* res = value ? node->high : node->low;
            memo[node] = res;
            return res;
        }

        QADDNode* low = restrict_var_impl(node->low, var, value, memo);
        QADDNode* high = restrict_var_impl(node->high, var, value, memo);

        QADDNode* res = make_node(node->var, low, high);
        memo[node] = res;
        return res;
    }

    QADDNode* restrict_var(QADDNode* node, int var, bool value) {
        std::unordered_map<QADDNode*, QADDNode*> memo;
        return restrict_var_impl(node, var, value, memo);
    }

    QADDNode* exists_var_impl(
        QADDNode* node,
        int var,
        std::unordered_map<QADDNode*, QADDNode*>& memo)
    {
        auto it = memo.find(node);
        if (it != memo.end()) {
            return it->second;
        }

        if (is_terminal(node)) {
            memo[node] = node;
            return node;
        }

        if (node->var > var) {
            memo[node] = node;
            return node;
        }

        QADDNode* res = nullptr;
        if (node->var == var) {
            // Ordered DD: target variable cannot reappear below this node.
            res = Apply(ApplyOp::JOIN, node->low, node->high);
        } else {
            QADDNode* low = exists_var_impl(node->low, var, memo);
            QADDNode* high = exists_var_impl(node->high, var, memo);
            res = make_node(node->var, low, high);
        }

        memo[node] = res;
        return res;
    }

    QADDNode* exists_var(QADDNode* node, int var) {
        std::unordered_map<QADDNode*, QADDNode*> memo;
        return exists_var_impl(node, var, memo);
    }

    QADDNode* exists_vars(QADDNode* node) {
        for (int i = 0; i < num_vars; ++i) {
            node = exists_var(node, src_var(i));
        }
        return node;
    }

    QADDNode* rename_vars_impl(
        QADDNode* node,
        std::unordered_map<QADDNode*, QADDNode*>& memo)
    {
        auto it = memo.find(node);
        if (it != memo.end()) {
            return it->second;
        }

        if (is_terminal(node)) {
            memo[node] = node;
            return node;
        }

        QADDNode* low = rename_vars_impl(node->low, memo);
        QADDNode* high = rename_vars_impl(node->high, memo);

        int var = node->var;
        if ((var % 2) == 1) {
            var = var - 1;
        }

        QADDNode* res = make_node(var, low, high);
        memo[node] = res;
        return res;
    }

    QADDNode* rename_vars(QADDNode* node) {
        std::unordered_map<QADDNode*, QADDNode*> memo;
        return rename_vars_impl(node, memo);
    }
};

} // namespace qts

#endif