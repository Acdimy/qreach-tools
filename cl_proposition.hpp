#ifndef CLASSICAL_PROPOSITION
#define CLASSICAL_PROPOSITION

#include <vector>
#include <string>
#include <iostream>
#include <algorithm>
#include <cassert>

/* Classical proposition: 
    An integer n: the number of variables (Boolean type)
    A vector of bit strings: each of length n, represent the disjunction of atomic terms. */

class ClassicalProposition {
private:
    bool normalized = false;
public:
    unsigned int n; // Number of variables
    std::vector<std::string> terms;
    ClassicalProposition() : n(0) {}
    ClassicalProposition(unsigned int n) : n(n) {
        terms.resize(1 << n, std::string(n, '0'));
    }
    ClassicalProposition(unsigned int n, std::vector<std::string> terms) : n(n), terms(terms) {
        assert(terms.size() == (1 << n));
    }
    ClassicalProposition(const ClassicalProposition& other) : n(other.n), terms(other.terms) {}
    void normalize() {
        // Normalize the terms by removing duplicates and sorting
        std::sort(terms.begin(), terms.end());
        terms.erase(std::unique(terms.begin(), terms.end()), terms.end());
        normalized = true;
    }
    void setValue(unsigned int index, const bool value) {
        assert(index < n);
        for (auto& term : terms) {
            term[index] = value ? '1' : '0';
        }
    }
    ClassicalProposition& operator=(const ClassicalProposition& other) {
        if (this != &other) {
            n = other.n;
            terms = other.terms;
        }
        return *this;
    }
    bool operator==(const ClassicalProposition& other) const {
        // Two propositions are equal if they have the same number of variables and the same normalized terms.
        ClassicalProposition normThis = *this;
        ClassicalProposition normOther = other;
        // normThis.normalize();
        // normOther.normalize();
        return (n == other.n) && (normThis.terms == normOther.terms);
    }
    bool satisfy(const std::string& assignment) const {
        assert(assignment.size() == n);
        for (const auto& term : terms) {
            bool satisfied = true;
            for (unsigned int i = 0; i < n; ++i) {
                if (term[i] == '1' && assignment[i] != '1') {
                    satisfied = false;
                    break;
                }
            }
            if (satisfied) {
                return true;
            }
        }
        return false;
    }
    std::vector<std::string> satisfyBit(std::vector<unsigned int> indexs, std::vector<bool> values) const {
        // Returns a list of terms that satisfy the proposition with all the given bit assignment.
        assert(indexs.size() == values.size());
        std::vector<std::string> satisfiedTerms;
        for (const auto& term : terms) {
            bool satisfied = true;
            for (size_t i = 0; i < indexs.size(); ++i) {
                unsigned int index = indexs[i];
                if (index >= n || term[index] != (values[i] ? '1' : '0')) {
                    satisfied = false;
                    break;
                }
            }
            if (satisfied) {
                satisfiedTerms.push_back(term);
            }
        }
        return satisfiedTerms;
    }
    std::vector<std::string> unsatisfyBit(std::vector<unsigned int> indexs, std::vector<bool> values) const {
        // Returns a list of terms that do not satisfy the proposition with all the given bit assignment.
        assert(indexs.size() == values.size());
        std::vector<std::string> unsatisfiedTerms;
        for (const auto& term : terms) {
            bool unsatisfied = false;
            for (size_t i = 0; i < indexs.size(); ++i) {
                unsigned int index = indexs[i];
                if (index >= n || term[index] != (values[i] ? '1' : '0')) {
                    unsatisfied = true;
                    break;
                }
            }
            if (unsatisfied) {
                unsatisfiedTerms.push_back(term);
            }
        }
        return unsatisfiedTerms;
    }
    void addTerm(const std::string& term) {
        if (n == 0) {
            assert(terms.empty());
            n = term.size();
            terms.push_back(term);
            return;
        }
        auto it = std::lower_bound(terms.begin(), terms.end(), term);
        if (it == terms.end() || *it != term) {
            terms.insert(it, term);
        }
    }
    void removeTerm(const std::string& term) {
        assert(term.size() == n);
        auto it = std::find(terms.begin(), terms.end(), term);
        if (it != terms.end()) {
            terms.erase(it);
        }
    }
    std::string toString() const {
        std::string result;
        for (const auto& term : terms) {
            if (!result.empty()) {
                result += " + ";
            }
            result += term;
        }
        return result.empty() ? "0" : result;
    }
    void print() const {
        std::cout << "Classical Proposition with " << n << " variables:" << std::endl;
        for (const auto& term : terms) {
            std::cout << term << std::endl;
        }
    }
    // Additional methods can be added for more functionality, such as conjunction, disjunction, etc.
};

#endif