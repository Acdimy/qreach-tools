#pragma once

#include <unordered_map>
#include <stdexcept>
#include <functional>
#include <fstream>
#include <iostream>
#include <unordered_set>
#include <vector>

#include "quantum_operation.hpp"

namespace qadd {

// =======================
// Node Definition
// =======================

enum class NodeType {
    INTERNAL,
    TERMINAL
};

struct QADDNode {
    NodeType type;

    // internal node
    int var = -1;
    QADDNode* low = nullptr;
    QADDNode* high = nullptr;

    // terminal (unified)
    QOperation val;

    // constructors
    QADDNode(const QOperation& v)
        : type(NodeType::TERMINAL), val(v) {}

    QADDNode(int v, QADDNode* l, QADDNode* h)
        : type(NodeType::INTERNAL), var(v), low(l), high(h) {}
};

// =======================
// Utility
// =======================

inline bool is_terminal(QADDNode* n) {
    return n->type == NodeType::TERMINAL;
}

inline int var(QADDNode* n) {
    return is_terminal(n) ? INT32_MAX : n->var;
}

// =======================
// Unique Table
// =======================

struct NodeKey {
    int var;
    QADDNode* low;
    QADDNode* high;

    bool operator==(const NodeKey& other) const {
        return var == other.var &&
               low == other.low &&
               high == other.high;
    }
};

struct NodeKeyHash {
    size_t operator()(const NodeKey& k) const {
        return std::hash<int>()(k.var) ^
               (std::hash<QADDNode*>()(k.low) << 1) ^
               (std::hash<QADDNode*>()(k.high) << 2);
    }
};

static std::unordered_map<NodeKey, QADDNode*, NodeKeyHash> uniqueTable;

// =======================
// Compute Table
// =======================

enum class ApplyOp {
    ADD,
    COMPOSE,
    JOIN,
    MEET,
    APPLY
};

struct ApplyKey {
    ApplyOp op;
    QADDNode* a;
    QADDNode* b;

    bool operator==(const ApplyKey& other) const {
        return op == other.op && a == other.a && b == other.b;
    }
};

struct ApplyKeyHash {
    size_t operator()(const ApplyKey& k) const {
        return std::hash<int>()((int)k.op) ^
               (std::hash<QADDNode*>()(k.a) << 1) ^
               (std::hash<QADDNode*>()(k.b) << 2);
    }
};

static std::unordered_map<ApplyKey, QADDNode*, ApplyKeyHash> computeTable;

inline void clear_compute_table() {
    computeTable.clear();
}

inline size_t compute_table_size() {
    return computeTable.size();
}

// =======================
// Terminal Unique Table
// =======================

struct TerminalKey {
    QOperation val;
    bool operator==(const TerminalKey& other) const {
        return val == other.val;
    }
};

struct TerminalKeyHash {
    size_t operator()(const TerminalKey& k) const {
        return std::hash<std::string>()(k.val.getName());
    }
};

static std::unordered_map<TerminalKey, QADDNode*, TerminalKeyHash> terminalTable;

// =======================
// Terminal Constructor
// =======================

inline QADDNode* make_terminal(const QOperation& val) {
    TerminalKey key{val};
    auto it = terminalTable.find(key);
    if (it != terminalTable.end()) {
        return it->second;
    }
    QADDNode* node = new QADDNode(val);
    terminalTable.emplace(key, node);
    return node;
}

// =======================
// Node Constructor
// =======================

inline QADDNode* make_node(int var, QADDNode* low, QADDNode* high) {
    if (low == high) return low;

    NodeKey key{var, low, high};
    auto it = uniqueTable.find(key);
    if (it != uniqueTable.end()) return it->second;

    QADDNode* node = new QADDNode(var, low, high);
    uniqueTable[key] = node;
    return node;
}

// =======================
// Terminal Apply
// =======================

inline QADDNode* apply_terminal(ApplyOp op, QADDNode* a, QADDNode* b) {
    const QOperation& A = a->val;
    const QOperation& B = b->val;

    switch (op) {
        case ApplyOp::ADD:
            return make_terminal(A.add(B));

        // case ApplyOp::COMPOSE:
        //     return make_terminal(A.compose(B));

        case ApplyOp::JOIN:
            return make_terminal(A.disjunction(B));

        case ApplyOp::MEET:
            return make_terminal(A.conjunction_simp(B));

        case ApplyOp::APPLY:
            // 语义：B.postImage(A)
            return make_terminal(B.postImage(A));

        default:
            throw std::runtime_error("Unknown ApplyOp");
    }
}

// =======================
// Apply Algorithm
// =======================

inline QADDNode* Apply(ApplyOp op, QADDNode* u1, QADDNode* u2) {
    ApplyKey key{op, u1, u2};

    auto it = computeTable.find(key);
    if (it != computeTable.end()) return it->second;

    if (is_terminal(u1) && is_terminal(u2)) {
        QADDNode* res = apply_terminal(op, u1, u2);
        computeTable[key] = res;
        return res;
    }

    int top = std::min(var(u1), var(u2));

    QADDNode* u1_low  = (var(u1) == top) ? u1->low  : u1;
    QADDNode* u1_high = (var(u1) == top) ? u1->high : u1;

    QADDNode* u2_low  = (var(u2) == top) ? u2->low  : u2;
    QADDNode* u2_high = (var(u2) == top) ? u2->high : u2;

    QADDNode* low  = Apply(op, u1_low,  u2_low);
    QADDNode* high = Apply(op, u1_high, u2_high);

    QADDNode* res = make_node(top, low, high);
    computeTable[key] = res;
    return res;
}

// =======================
// Utility for QOperation
// =======================

inline int get_dimension(const QOperation& op) {
    assert(op.type == false);
    if (op.isIdentity) {
        return (1 << op.realqNum);
    }
    // Assume that a projection is represented by orthogonal vectors.
    return op.oplist.size();
}

inline size_t count_nodes(QADDNode* root) {
    if (!root) return 0;

    std::unordered_set<QADDNode*> visited;
    std::vector<QADDNode*> worklist{root};

    while (!worklist.empty()) {
        QADDNode* node = worklist.back();
        worklist.pop_back();

        if (!node || visited.find(node) != visited.end()) {
            continue;
        }

        visited.insert(node);
        if (!is_terminal(node)) {
            worklist.push_back(node->low);
            worklist.push_back(node->high);
        }
    }

    return visited.size();
}

inline size_t total_unique_node_count() {
    return uniqueTable.size() + terminalTable.size();
}


void printQADD(QADDNode* node, const std::string& filename) {
    // Generate a DOT file for visualization
    std::ofstream outFile(filename);
    if (!outFile.is_open()) {
        std::cerr << "File cannot open: " << filename << std::endl;
        return;
    }

    std::unordered_map<QADDNode*, int> nodeIds;
    int idCounter = 0;

    std::function<void(QADDNode*)> dfs = [&](QADDNode* n) {
        if (nodeIds.count(n)) return;
        int id = idCounter++;
        nodeIds[n] = id;

        if (is_terminal(n)) {
            outFile << "  node" << id << " [label=\"" << n->val.getName() << "\", shape=ellipse];" << std::endl;
        } else {
            outFile << "  node" << id << " [label=\"X" << n->var << "\", shape=box];" << std::endl;
            dfs(n->low);
            dfs(n->high);
            outFile << "  node" << id << " -> node" << nodeIds[n->low] << " [label=\"0\"];" << std::endl;
            outFile << "  node" << id << " -> node" << nodeIds[n->high] << " [label=\"1\"];" << std::endl;
        }
    };

    outFile << "digraph QADD {" << std::endl;
    dfs(node);
    outFile << "}" << std::endl;

    outFile.close();
    std::cout << "DOT file complete: " << filename << std::endl;
}

// Debug helper: dump terminal QOperation fields grouped by getName().
// If firstLocationMap is provided, it also prints the first location id that
// reaches each terminal node.
inline void debug_print_terminals(
    QADDNode* root,
    const std::string& outFilename = "terminals.txt",
    const std::unordered_map<QADDNode*, int>* firstLocationMap = nullptr)
{
    if (!root) return;

    std::unordered_map<std::string, std::vector<QADDNode*>> termMap;
    std::unordered_set<QADDNode*> visited;

    std::function<void(QADDNode*)> dfs = [&](QADDNode* n) {
        if (!n) return;
        if (visited.find(n) != visited.end()) return;
        visited.insert(n);

        if (is_terminal(n)) {
            termMap[n->val.getName()].push_back(n);
            return;
        }

        dfs(n->low);
        dfs(n->high);
    };

    dfs(root);

    std::ofstream ofs(outFilename);
    if (!ofs.is_open()) {
        std::cerr << "Cannot open terminal dump file: " << outFilename << std::endl;
        return;
    }

    auto dump_qoperation = [&](std::ostream& os, const QOperation& op) {
        os << "    getName: " << op.getName() << "\n";
        os << "    type: " << (op.type ? "gate/super-op" : "projection") << "\n";
        os << "    normalized: " << (op.normalized ? "true" : "false") << "\n";
        os << "    isIdentity: " << (op.isIdentity ? "true" : "false") << "\n";
        os << "    isProj: " << op.isProj << "\n";
        os << "    qNum: " << op.qNum << "\n";
        os << "    realqNum: " << op.realqNum << "\n";
        os << "    oplist.size: " << op.oplist.size() << "\n";
    };

    for (const auto& p : termMap) {
        const std::string& name = p.first;
        const auto& nodes = p.second;

        ofs << "Name group: \"" << name << "\"  count=" << nodes.size() << "\n";
        QADDNode* ref = nodes.empty() ? nullptr : nodes.front();
        for (size_t i = 0; i < nodes.size(); ++i) {
            QADDNode* np = nodes[i];
            ofs << "  terminal[" << i << "] ptr=" << np << "\n";
            if (firstLocationMap) {
                auto lit = firstLocationMap->find(np);
                if (lit != firstLocationMap->end()) {
                    ofs << "    first_location_id: " << lit->second << "\n";
                } else {
                    ofs << "    first_location_id: -1\n";
                }
            }
            dump_qoperation(ofs, np->val);
            if (ref && np != ref) {
                ofs << "    same_as_first_by_operator_eq: "
                    << (np->val.compare(ref->val)) << "\n";
            } else {
                ofs << "    same_as_first_by_operator_eq: true\n";
            }
            ofs << "\n";
        }
    }

    ofs.close();
    std::cout << "Terminal dump written to " << outFilename << std::endl;
}

} // namespace qadd