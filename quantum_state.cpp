#include "quantum_state.h"

QuantumState::QuantumState() {}
QuantumState::~QuantumState() {}

// ********* CFLOBDDQuantumState *********
CFLOBDDQuantumState::CFLOBDDQuantumState()
{
}

CFLOBDDQuantumState::CFLOBDDQuantumState(CFLOBDD_COMPLEX_BIG g)
{
    state = g;   
}

CFLOBDDQuantumState::~CFLOBDDQuantumState()
{
}

void CFLOBDDQuantumState::Print()
{
    state.print(std::cout);
}