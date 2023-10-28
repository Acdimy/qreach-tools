#ifndef CFLOBDD_FLOAT_BOOST_FLOAT_BOOST_GUARD
#define CFLOBDD_FLOAT_BOOST_FLOAT_BOOST_GUARD

//
//    Copyright (c) 1999, 2017 Thomas W. Reps
//    All Rights Reserved.
//
//    This software is furnished under a license and may be used and
//    copied only in accordance with the terms of such license and the
//    inclusion of the above copyright notice.  This software or any
//    other copies thereof or any derivative works may not be provided
//    or otherwise made available to any other person.  Title to and
//    ownership of the software and any derivative works is retained
//    by Thomas W. Reps.
//
//    THIS IMPLEMENTATION MAY HAVE BUGS, SOME OF WHICH MAY HAVE SERIOUS
//    CONSEQUENCES.  THOMAS W. REPS PROVIDES THIS SOFTWARE IN ITS "AS IS"
//    CONDITION, AND ANY EXPRESS OR IMPLIED WARRANTIES, INCLUDING,
//    BUT NOT LIMITED TO, THE IMPLIED WARRANTIES OF MERCHANTABILITY
//    AND FITNESS FOR A PARTICULAR PURPOSE ARE DISCLAIMED.  IN NO EVENT SHALL
//    THOMAS W. REPS BE LIABLE FOR ANY DIRECT, INDIRECT, INCIDENTAL,
//    SPECIAL, EXEMPLARY, OR CONSEQUENTIAL DAMAGES (INCLUDING, BUT NOT LIMITED
//    TO, PROCUREMENT OF SUBSTITUTE GOODS OR SERVICES; LOSS OF USE, DATA, OR
//    PROFITS; OR BUSINESS INTERRUPTION) HOWEVER CAUSED AND ON ANY THEORY OF
//    LIABILITY, WHETHER IN CONTRACT, STRICT LIABILITY, OR TORT (INCLUDING
//    NEGLIGENCE OR OTHERWISE) ARISING IN ANY WAY OUT OF THE USE OF THIS
//    SOFTWARE, EVEN IF ADVISED OF THE POSSIBILITY OF SUCH DAMAGE.
//

#include <iostream>
#include <fstream>
#include <boost/multiprecision/cpp_dec_float.hpp>
#include "cflobdd_t.h"

namespace mp = boost::multiprecision;

namespace CFL_OBDD {

	typedef mp::cpp_dec_float_100 BIG_FLOAT;
	//typedef mp::number<mp::cpp_dec_float<1000> > BIG_FLOAT;
	typedef CFLOBDD_T<BIG_FLOAT> CFLOBDD_FLOAT_BOOST;
	// typedef CFLOBDD_T<int> CFLOBDD;

	// CFLOBDD_FLOAT_BOOST-creation operations --------------------------------------
// 	extern CFLOBDD_FLOAT_BOOST MkTrue();                             // Representation of \x.true
// 	extern CFLOBDD_FLOAT_BOOST MkFalse();                            // Representation of \x.false
// 	extern CFLOBDD_FLOAT_BOOST MkProjection(unsigned int i);         // Representation of \x.x_i
// 	extern CFLOBDD_FLOAT_BOOST MkAdditionInterleavedRecursive();     // Representation of addition relation { (xi yi zi _)* | vec{x} + vec{y} = vec{z} }, created recursively (expensive)
// 	extern CFLOBDD_FLOAT_BOOST MkAdditionInterleavedBruteForce();    // Representation of addition relation { (xi yi zi _)* | vec{x} + vec{y} = vec{z} }, created by brute force (expensive)
// 	extern CFLOBDD_FLOAT_BOOST MkAdditionInterleaved();              // Representation of addition relation { (xi yi zi _)* | vec{x} + vec{y} = vec{z} }, created by tabulation
// 	extern CFLOBDD_FLOAT_BOOST MkX(unsigned int i);                  // Bits of i -> vec{x} for addition relation
// 	extern CFLOBDD_FLOAT_BOOST MkY(unsigned int i);                  // Bits of i -> vec{y} for addition relation
// 	extern CFLOBDD_FLOAT_BOOST MkZ(unsigned int i);                  // Bits of i -> vec{z} for addition relation
// 	extern CFLOBDD_FLOAT_BOOST MkXBit(unsigned int j);
// 	extern CFLOBDD_FLOAT_BOOST MkYBit(unsigned int j);
// 	extern CFLOBDD_FLOAT_BOOST MkZBit(unsigned int j);
// 	extern CFLOBDD_FLOAT_BOOST MkMultiplicationInterleavedBruteForce(); // Representation of multiplication relation { (xi yi zi _)* | vec{x} * vec{y} = vec{z} }
// 	extern bool FactorZ(CFLOBDD_FLOAT_BOOST R, unsigned int z, CFLOBDD_FLOAT_BOOST &f1, CFLOBDD_FLOAT_BOOST &f2, unsigned int &v1, unsigned int &v2);   // Return true if z has a non-trivial factorization
// 	extern CFLOBDD_FLOAT_BOOST MkParity();                           // Representation of parity function
// 	extern CFLOBDD_FLOAT_BOOST MkStepUpOneFourth();                  // Representation of step function
// 	extern CFLOBDD_FLOAT_BOOST MkStepDownOneFourth();                // Representation of step function
// #ifdef ARBITRARY_STEP_FUNCTIONS
// 	extern CFLOBDD_FLOAT_BOOST MkStepUp(unsigned int i);           // Representation of step function
// 	extern CFLOBDD_FLOAT_BOOST MkStepDown(unsigned int i);         // Representation of step function
// 	extern CFLOBDD_FLOAT_BOOST MkPulseUp(unsigned int i, unsigned int j);
// 	// Representation of pulse function
// 	extern CFLOBDD_FLOAT_BOOST MkPulseDown(unsigned int i, unsigned int j);
// 	// Representation of step function
// #endif

// 	extern CFLOBDD_FLOAT_BOOST shiftAtoB(const CFLOBDD_FLOAT_BOOST f, const unsigned int levelAtWhichToShift = CFLOBDD_MAX_LEVEL);
// 	// Special operation used for Karatsuba multiplication
// 	// with ASCENDING_INTERLEAVED variable ordering
// 	extern CFLOBDD_FLOAT_BOOST shiftBtoA(const CFLOBDD_FLOAT_BOOST f, const unsigned int levelAtWhichToShift = CFLOBDD_MAX_LEVEL);
// 	// Special operation used for Karatsuba multiplication
// 	// with DESCENDING_INTERLEAVED variable ordering

// 	extern CFLOBDD_FLOAT_BOOST shiftAtoBAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, const CFLOBDD_FLOAT_BOOST f);
// 	extern CFLOBDD_FLOAT_BOOST shiftBtoAAtLevelOne(Hashset<CFLOBDDNode> *visitedNodes, unsigned int &totalVisitCount, unsigned int &redundantVisitCount, const CFLOBDD_FLOAT_BOOST f);
// 	extern CFLOBDD_FLOAT_BOOST duplicateAinBAtLevelOne(const CFLOBDD_FLOAT_BOOST f);

// 	// Unary operations on CFLOBDD_FLOAT_BOOSTs --------------------------------------
// 	extern CFLOBDD_FLOAT_BOOST MkNot(CFLOBDD_FLOAT_BOOST f);                     // \f.!f

// 	// Binary operations on CFLOBDD_FLOAT_BOOSTs --------------------------------------
// 	extern CFLOBDD_FLOAT_BOOST MkAnd(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);          // \f.\g.(f && g)
// 	extern CFLOBDD_FLOAT_BOOST MkNand(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);         // \f.\g.!(f && g)
// 	extern CFLOBDD_FLOAT_BOOST MkOr(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);           // \f.\g.(f || g)
// 	extern CFLOBDD_FLOAT_BOOST MkNor(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);          // \f.\g.!(f || g)
// 	extern CFLOBDD_FLOAT_BOOST MkIff(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);          // \f.\g.(f == g)
// 	extern CFLOBDD_FLOAT_BOOST MkExclusiveOr(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);  // \f.\g.(f != g)
// 	extern CFLOBDD_FLOAT_BOOST MkImplies(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);      // \f.\g.(!f || g)
// 	extern CFLOBDD_FLOAT_BOOST MkMinus(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);        // \f.\g.(f && !g)
// 	extern CFLOBDD_FLOAT_BOOST MkQuotient(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);     // \f.\g.(!g || f)
// 	extern CFLOBDD_FLOAT_BOOST MkNotQuotient(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);  // \f.\g.(g && !f)
// 	extern CFLOBDD_FLOAT_BOOST MkFirst(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);        // \f.\g.f
// 	extern CFLOBDD_FLOAT_BOOST MkNotFirst(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);     // \f.\g.!f
// 	extern CFLOBDD_FLOAT_BOOST MkSecond(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);       // \f.\g.g
// 	extern CFLOBDD_FLOAT_BOOST MkNotSecond(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);    // \f.\g.!g

// 	// extern CFLOBDD_FLOAT_BOOST MkPlus(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);      // \f.\g.(f + g)
// 	// extern CFLOBDD_FLOAT_BOOST MkTimes(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g);     // \f.\g.(f * g)

// 	// N-ary operations on CFLOBDD_FLOAT_BOOSTs --------------------------------------
// 	extern CFLOBDD_FLOAT_BOOST MkAnd(int N, ...);                    // \f1. ... \fk.(f1 && ... && fk)
// 	extern CFLOBDD_FLOAT_BOOST MkNand(int N, ...);                   // \f1. ... \fk.!(f1 && ... && fk)
// 	extern CFLOBDD_FLOAT_BOOST MkOr(int N, ...);                     // \f1. ... \fk.(f1 || ... || fk)
// 	extern CFLOBDD_FLOAT_BOOST MkNor(int N, ...);                    // \f1. ... \fk.!(f1 || ... || fk)

// 	// Ternary operations on CFLOBDD_FLOAT_BOOSTs --------------------------------------
// 	extern CFLOBDD_FLOAT_BOOST MkIfThenElse(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g, CFLOBDD_FLOAT_BOOST h);  // \a.\b.\c.(a && b) || (!a && c)
// 	extern CFLOBDD_FLOAT_BOOST MkNegMajority(CFLOBDD_FLOAT_BOOST f, CFLOBDD_FLOAT_BOOST g, CFLOBDD_FLOAT_BOOST h); // \a.\b.\c.(b && !a) || (c && !a) || (b && c)

// 	extern CFLOBDD_FLOAT_BOOST MkRestrict(CFLOBDD_FLOAT_BOOST f, unsigned int i, bool val);  // \f. f | (x_i = val)
// 	extern CFLOBDD_FLOAT_BOOST MkExists(CFLOBDD_FLOAT_BOOST f, unsigned int i);              // \f. exists x_i : f
// 	extern CFLOBDD_FLOAT_BOOST MkForall(CFLOBDD_FLOAT_BOOST f, unsigned int i);              // \f. forall x_i : f
// 	extern CFLOBDD_FLOAT_BOOST MkCompose(CFLOBDD_FLOAT_BOOST f, int i, CFLOBDD_FLOAT_BOOST g);  // \f. f | (x_i = g)

// 	// Other operations on CFLOBDD_FLOAT_BOOSTs ---------------------------------------------------
// 	extern bool DependsOn(CFLOBDD_FLOAT_BOOST f, int i);
// 	extern bool IsPositiveCube(CFLOBDD_FLOAT_BOOST f);
// 	extern bool IsPositiveCubeInt(CFLOBDD_FLOAT_BOOST f, int least);
// 	extern bool SupportSetIs(CFLOBDD_FLOAT_BOOST f, const apvector<int> &ss);
// 	extern apvector<int> *GetSupportSet(CFLOBDD_FLOAT_BOOST f);

} // namespace CFL_OBDD

#endif
