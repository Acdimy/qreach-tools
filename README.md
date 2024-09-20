# QReach: A Reachability Analysis Tool for Quantum Markov Chains

### Introduction
This tool aims to providing reachability analysis for quantum Markov chains, which is a basic step of further quantum model checking tasks.

Quantum Markov chains (QMCs for short) have been adopted as a fundamental model of many quantum information processing systems (e.g. quantum communication protocols, semantics of quantum programs, etc). So, we choose to use QMCs as our system model. As is well-known, reachability analysis is a core task in classical model checking algorithms. In the quantum case, indeed, reachability analysis has been applied in quantum communication, quantum control and termination analysis of quantum programs among many others. Therefore, we focus on the issue of reachability analysis of quantum Markov chains.

In this tool, we incorporate quantum DDs into quantum model checking for the first time. The backend we used in this tool is Context-Free-Language Ordered Binary Decision Diagrams (CFLOBDD for short), one of the most efficient quantum DDs as the backend of our tool to provide support for functionalities. It is flexible for changing other backends. Supported by the efficiency of the DD representation, our tool is well scalable and has the potential to be expanded into large-scale quantum circuit model checkers in the future.

Some previous works are referred in this project: [Quasimodo simulator](https://github.com/trishullab/Quasimodo), [Original CFLOBDD](https://github.com/trishullab/cflobdd). We build some components of QReach based on these programs, but did some important modifications.

### Usage and Installation

### Demo

下一阶段工作内容：
标准化的输入，包括接受更多量子转移系统作为输入，以及搭建一个顶层的框架，希望能够进行partition与circuit knitting等预处理
抽象化的中间表示，包括尽量使用形式化、符号化的逻辑语言
拓展更多数据结构接口

最近看protocal的一个启发：分支不仅仅可以是measurement或者non-unitary，而可能是协议，一些non deterministic的操作。我的模型似乎处理不了这种情况，而SVMC可以。
写一个Algorithm出来！
