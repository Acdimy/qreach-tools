#pragma once

#include <unordered_map>
#include <stdexcept>
#include <functional>
#include <fstream>

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

} // namespace qadd