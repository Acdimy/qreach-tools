#ifndef _QUANTUM_STATE
#define _QUANTUM_STATE

#include "cflobdd/CFLOBDD/matrix1234_complex_float_boost.h"
#include <random>


class QuantumState {
    public:
        // Constructor
        QuantumState();
        // Destructor
        virtual ~QuantumState();
        virtual void Print() = 0;
};

using namespace CFL_OBDD;

class CFLOBDDQuantumState : public QuantumState {
    public:
        CFLOBDDQuantumState(CFLOBDD_COMPLEX_BIG state);
        CFLOBDDQuantumState();
        ~CFLOBDDQuantumState();
        void Print();
        CFLOBDDQuantumState GetState() {return state;}

    private:
        CFLOBDD_COMPLEX_BIG state;
};

#endif