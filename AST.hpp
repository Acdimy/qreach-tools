#include <iostream>
#include <memory>

// 抽象基类 Node
class Node {
public:
    virtual ~Node() = default;
    virtual void print() const = 0; // 用于调试，输出树结构
};

// 叶子节点，表示 vector 的一个区间 [left, right]
class LeafNode : public Node {
public:
    LeafNode(size_t left, size_t right) : left_(left), right_(right) {
        if (left_ > right_) {
            throw std::invalid_argument("Invalid range: left > right");
        }
    }
    void print() const override {
        std::cout << "[" << left_ << ", " << right_ << "]";
    }
    size_t getLeft() const { return left_; }
    size_t getRight() const { return right_; }
private:
    size_t left_;
    size_t right_;
};

// 操作节点，表示加法或乘法
class OpNode : public Node {
public:
    enum OpType { ADD, MUL };
    OpNode(OpType type, std::unique_ptr<Node> left, std::unique_ptr<Node> right)
        : type_(type), left_(std::move(left)), right_(std::move(right)) {}

    void print() const override {
        std::cout << "(";
        left_->print();
        std::cout << (type_ == ADD ? " + " : " * ");
        right_->print();
        std::cout << ")";
    }
private:
    OpType type_;
    std::unique_ptr<Node> left_;
    std::unique_ptr<Node> right_;
};

// 连接两个 AST 的函数
std::unique_ptr<Node> combineTrees(std::unique_ptr<Node> tree1, std::unique_ptr<Node> tree2, OpNode::OpType op) {
    return std::make_unique<OpNode>(op, std::move(tree1), std::move(tree2));
}

// // 测试代码
// int main() {
//     // 创建第一个 AST: 表示区间 [0, 2]
//     auto tree1 = std::make_unique<LeafNode>(0, 2);

//     // 创建第二个 AST: 表示区间 [3, 5] + [6, 7]
//     auto leaf1 = std::make_unique<LeafNode>(3, 5);
//     auto leaf2 = std::make_unique<LeafNode>(6, 7);
//     auto tree2 = std::make_unique<OpNode>(OpNode::ADD, std::move(leaf1), std::move(leaf2));

//     // 将 tree1 和 tree2 用乘法连接成新树
//     auto newTree = combineTrees(std::move(tree1), std::move(tree2), OpNode::MUL);

//     // 输出新树的结构
//     std::cout << "New AST: ";
//     newTree->print();
//     std::cout << std::endl; // 预期输出: ([0, 2] * ([3, 5] + [6, 7]))

//     return 0;
// }
