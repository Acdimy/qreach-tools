#ifndef _QASM_PARSER
#define _QASM_PARSER

#include <iostream>
#include <fstream>
#include <sstream>
#include <string>
#include <vector>
#include <regex>
#include <stdexcept>
#include <optional>

struct ClassicalCondition {
    std::string creg;
    int index;
    int value;
};

struct QASMOperation {
    std::string gate;
    std::vector<std::pair<std::string, int>> args;
    std::optional<ClassicalCondition> condition; // 可选控制条件
};

struct QASMProgram {
    std::string qreg_name;
    int qreg_size;
    std::string creg_name;
    int creg_size;
    std::vector<QASMOperation> operations;
};

QASMProgram parse_qasm_file(const std::string& filename) {
    std::ifstream infile(filename);
    if (!infile) {
        throw std::runtime_error("Failed to open file.");
    }

    std::string line;
    QASMProgram program;

    // 跳过前两行（版本声明与 include）
    std::getline(infile, line);
    std::getline(infile, line);

    // 解析量子寄存器
    std::getline(infile, line);
    std::smatch match;
    if (std::regex_match(line, match, std::regex(R"(qreg\s+(\w+)\[(\d+)\];)"))) {
        program.qreg_name = match[1];
        program.qreg_size = std::stoi(match[2]);
    } else {
        throw std::runtime_error("Invalid qreg line.");
    }

    // 解析经典寄存器
    std::getline(infile, line);
    if (std::regex_match(line, match, std::regex(R"(creg\s+(\w+)\[(\d+)\];)"))) {
        program.creg_name = match[1];
        program.creg_size = std::stoi(match[2]);
    } else {
        throw std::runtime_error("Invalid creg line.");
    }

    // 解析操作指令
    while (std::getline(infile, line)) {
        if (line.empty()) continue;
        QASMOperation op;

        // measure q[x] -> c[y];
        if (std::regex_match(line, match, std::regex(R"(measure\s+(\w+)\[(\d+)\]\s*->\s*(\w+)\[(\d+)\];)"))) {
            op.gate = "measure";
            op.args.emplace_back(match[1], std::stoi(match[2]));
            op.args.emplace_back(match[3], std::stoi(match[4]));
        }
        // 两比特门：cx q[a], q[b];
        else if (std::regex_match(line, match, std::regex(R"((\w+)\s+(\w+)\[(\d+)\],\s*(\w+)\[(\d+)\];)"))) {
            op.gate = match[1];
            op.args.emplace_back(match[2], std::stoi(match[3]));
            op.args.emplace_back(match[4], std::stoi(match[5]));
        }
        // 单比特门：x q[a];
        else if (std::regex_match(line, match, std::regex(R"((\w+)\s+(\w+)\[(\d+)\];)"))) {
            op.gate = match[1];
            op.args.emplace_back(match[2], std::stoi(match[3]));
        } else if (std::regex_match(line, match, std::regex(R"(^if\s*\(\s*(\w+)\[(\d+)\]\s*==\s*(\d+)\s*\)\s*(.+);$)"))) {
            ClassicalCondition cond { match[1], std::stoi(match[2]), std::stoi(match[3]) };
            std::string op_line = match[4];

            // 递归处理 op_line：
            QASMOperation op;

            if (std::regex_match(op_line, match, std::regex(R"((\w+)\s+(\w+)\[(\d+)\],\s*(\w+)\[(\d+)\])"))) {
                op.gate = match[1];
                op.args = {{match[2], std::stoi(match[3])}, {match[4], std::stoi(match[5])}};
            } else if (std::regex_match(op_line, match, std::regex(R"((\w+)\s+(\w+)\[(\d+)\])"))) {
                op.gate = match[1];
                op.args = {{match[2], std::stoi(match[3])}};
            } else if (std::regex_match(op_line, match, std::regex(R"(measure\s+(\w+)\[(\d+)\]\s*->\s*(\w+)\[(\d+)\])"))) {
                op.gate = "measure";
                op.args = {{match[1], std::stoi(match[2])}, {match[3], std::stoi(match[4])}};
            } else {
                throw std::runtime_error("Invalid quantum operation in if-statement: " + op_line);
            }

            op.condition = cond;
            program.operations.push_back(op);
        } else {
            throw std::runtime_error("Invalid operation line: " + line);
        }

        program.operations.push_back(op);
    }

    return program;
}

// 示例打印函数
void print_qasm_program(const QASMProgram& program) {
    std::cout << "qreg " << program.qreg_name << "[" << program.qreg_size << "];\n";
    std::cout << "creg " << program.creg_name << "[" << program.creg_size << "];\n";
    for (const auto& op : program.operations) {
        if (op.condition) {
            std::cout << "if (" << op.condition->creg << "[" << op.condition->index << "] == "
                      << op.condition->value << ") ";
        }
        std::cout << op.gate << " ";
        for (size_t i = 0; i < op.args.size(); ++i) {
            std::cout << op.args[i].first << "[" << op.args[i].second << "]";
            if (op.gate == "measure" && i == 0) std::cout << " -> ";
            else if (i < op.args.size() - 1) std::cout << ", ";
        }
        std::cout << ";\n";
    }
}

#endif