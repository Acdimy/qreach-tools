//
//    Copyright (c) 2017, 2018 Thomas W. Reps
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

#include <cassert>
#include <cstdio>
#include <iostream>
#include <fstream>
#include <map>
#include <unordered_map>
#include <tuple>
#include <algorithm>
#include "matrix1234_node.h"
#include "cflobdd_node.h"
#include "reduction_map.h"
#include "hash.h"
#include "hashset.h"
#include "cflobdd_top_node_matmult_map.h"

namespace CFL_OBDD {
	//CFLOBDDLinearMapMemoTableRefPtr projectMemoTable;
	std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash> matmult_hashMap;
	std::unordered_map<MatMultPairWithInfo, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPairWithInfo::MatMultPairWithInfoHash> matmult_hashMap_info;
	std::unordered_map<std::string, CFLOBDDNodeHandle> cnot_hashMap;
	std::unordered_map<std::string, CFLOBDDNodeHandle> cp_hashMap;
	std::unordered_map<std::string, CFLOBDDNodeHandle> cswap_hashMap;
	std::unordered_map<std::string, CFLOBDDNodeHandle> swap_hashMap;
	std::unordered_map<std::string, CFLOBDDNodeHandle> iswap_hashMap;
	std::unordered_map<std::string, CFLOBDDNodeHandle> ccnot_hashMap;
	std::vector<ReturnMapHandle<int>> commonly_used_return_maps;// m0, m1, m01, m10

	void InitReturnMapHandles(){
		ReturnMapHandle<int> m0, m1, m01, m10;
		m0.AddToEnd(0);
		m0.Canonicalize();
		m1.AddToEnd(1);
		m1.Canonicalize();
		m01.AddToEnd(0);
		m01.AddToEnd(1);
		m01.Canonicalize();
		m10.AddToEnd(1);
		m10.AddToEnd(0);
		m10.Canonicalize();
		commonly_used_return_maps.push_back(m0);
		commonly_used_return_maps.push_back(m1);
		commonly_used_return_maps.push_back(m01);
		commonly_used_return_maps.push_back(m10);
	}

	void Matrix1234InitializerNode()
	{
		// Empty for now
		//projectMemoTable = new CFLOBDDLinearMapMemoTable;
		InitReturnMapHandles();
		return;
	}
	
	void clearMultMap(){
		// std::cout << "mapSize: " << matmult_hashMap.size() << std::endl;
		matmult_hashMap.clear();
	}

	// Recursive routine to create a representation of the IdRelation where bits are interleaved
	// 
	// Note: with CFLOBDDMaxLevel == 1024, this routine crashed. (Stack too deep?)
	CFLOBDDNodeHandle MkIdRelationInterleavedNode(unsigned int level)
	{
		CFLOBDDInternalNode *n;

		if (level == 0) {
			return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
		}
		else if (level == 1) {
			n = new CFLOBDDInternalNode(level);

			CFLOBDDNodeHandle temp = CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
			n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			n->BConnection[1] = Connection(temp, commonly_used_return_maps[3]);//m10
		}
		else {  // Create an appropriate CFLOBDDInternalNode
			n = new CFLOBDDInternalNode(level);

			CFLOBDDNodeHandle temp = MkIdRelationInterleavedNode(level - 1);
			n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[1]);//m1
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	} // MkIdRelationInterleavedNode

	CFLOBDDNodeHandle MkNegationMatrixInterleavedNode(unsigned int level)
	{
		CFLOBDDInternalNode *n;

		if (level == 0) {
			return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
		}
		else if (level == 1) {
			n = new CFLOBDDInternalNode(level);

			CFLOBDDNodeHandle temp = CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
			n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			n->BConnection[1] = Connection(temp, commonly_used_return_maps[3]);//m10
		}
		else {  // Create an appropriate CFLOBDDInternalNode
			n = new CFLOBDDInternalNode(level);

			CFLOBDDNodeHandle temp = MkNegationMatrixInterleavedNode(level - 1);
			n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[1] = Connection(temp, commonly_used_return_maps[2]);//m01
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[1]);//m1
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkWalshInterleavedNode(unsigned int i)
	{
		assert(i >= 1);
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 1) {  // Base case
			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
		}
		else {
			CFLOBDDNodeHandle temp = MkWalshInterleavedNode(i - 1);
			n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			n->BConnection[1] = Connection(temp, commonly_used_return_maps[3]);//m10
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	} // MkWalshInterleavedNode

	CFLOBDDNodeHandle MkPauliYInterleavedNode(unsigned int i)
	{
		assert(i >= 1);
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 1) {  // Base case
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			CFLOBDDReturnMapHandle m20;
			m20.AddToEnd(2); m20.AddToEnd(0); m20.Canonicalize();
			// HERE!!! m20
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m20);//m20
		}
		else {
			CFLOBDDNodeHandle temp = MkPauliYInterleavedNode(i - 1);
			CFLOBDDReturnMapHandle m012;
			m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
			n->AConnection = Connection(temp, m012);
			n->numBConnections = 3;
			n->BConnection = new Connection[n->numBConnections];
			CFLOBDDReturnMapHandle m0, m021;
			m0.AddToEnd(0); m0.Canonicalize();
			m021.AddToEnd(0); m021.AddToEnd(2); m021.AddToEnd(1); m021.Canonicalize();
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i - 1], m0);
			n->BConnection[1] = Connection(temp, m012);
			n->BConnection[2] = Connection(temp, m021);
		}
		n->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkPauliZInterleavedNode(unsigned int i)
	{
		assert(i >= 1);
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 1) {  // Base case
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			CFLOBDDReturnMapHandle m12;
			m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);//m12
		}
		else {
			CFLOBDDNodeHandle temp = MkPauliZInterleavedNode(i - 1);
			CFLOBDDReturnMapHandle m012;
			m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
			n->AConnection = Connection(temp, m012);
			n->numBConnections = 3;
			n->BConnection = new Connection[n->numBConnections];
			CFLOBDDReturnMapHandle m1, m210;
			m1.AddToEnd(1); m1.Canonicalize();
			m210.AddToEnd(2); m210.AddToEnd(1); m210.AddToEnd(0); m210.Canonicalize();
			n->BConnection[0] = Connection(temp, m012);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i-1], m1);
			n->BConnection[2] = Connection(temp, m210);
		}
		n->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkSGateInterleavedNode(unsigned int i)
	{	// USUALLY OCCURS I=1 CASE ONLY!!! It is different from the conjugatation.
		assert(i >= 1);
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 1) {  // Base case
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			CFLOBDDReturnMapHandle m12;
			m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);//m12
			n->numExits = 3;
		}
		else if (i == 2) {
			CFLOBDDNodeHandle temp = MkSGateInterleavedNode(i - 1);
			CFLOBDDReturnMapHandle m012;
			m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
			n->AConnection = Connection(temp, m012);
			n->numBConnections = 3;
			n->BConnection = new Connection[n->numBConnections];
			CFLOBDDReturnMapHandle m1, m213;
			m1.AddToEnd(1); m1.Canonicalize();
			m213.AddToEnd(2); m213.AddToEnd(1); m213.AddToEnd(3); m213.Canonicalize();
			n->BConnection[0] = Connection(temp, m012);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i-1], m1);
			n->BConnection[2] = Connection(temp, m213);
			n->numExits = 4;
		}
		else if (i == 3) {
			CFLOBDDNodeHandle temp = MkSGateInterleavedNode(i - 1);
			CFLOBDDReturnMapHandle m0123;
			m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
			n->AConnection = Connection(temp, m0123);
			n->numBConnections = 4;
			n->BConnection = new Connection[n->numBConnections];
			CFLOBDDReturnMapHandle m1, m2134, m3140;
			m1.AddToEnd(1); m1.Canonicalize();
			m2134.AddToEnd(2); m2134.AddToEnd(1); m2134.AddToEnd(3); m2134.AddToEnd(4); m2134.Canonicalize();
			m3140.AddToEnd(3); m3140.AddToEnd(1); m3140.AddToEnd(4); m3140.AddToEnd(0); m3140.Canonicalize();
			n->BConnection[0] = Connection(temp, m0123);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i-1], m1);
			n->BConnection[2] = Connection(temp, m2134);
			n->BConnection[3] = Connection(temp, m3140);
			n->numExits = 5;	
		}
		else if (i == 4) {
			CFLOBDDNodeHandle temp = MkSGateInterleavedNode(i - 1);
			CFLOBDDReturnMapHandle m01234;
			m01234.AddToEnd(0); m01234.AddToEnd(1); m01234.AddToEnd(2); m01234.AddToEnd(3); m01234.AddToEnd(4); m01234.Canonicalize();
			n->AConnection = Connection(temp, m01234);
			n->numBConnections = 5;
			n->BConnection = new Connection[n->numBConnections];
			CFLOBDDReturnMapHandle m1, m21340, m31402, m41023;
			m1.AddToEnd(1); m1.Canonicalize();
			m21340.AddToEnd(2); m21340.AddToEnd(1); m21340.AddToEnd(3); m21340.AddToEnd(4); m21340.AddToEnd(0); m21340.Canonicalize();
			m31402.AddToEnd(3); m31402.AddToEnd(1); m31402.AddToEnd(4); m31402.AddToEnd(0); m31402.AddToEnd(2); m31402.Canonicalize();
			m41023.AddToEnd(4); m41023.AddToEnd(1); m41023.AddToEnd(0); m41023.AddToEnd(2); m41023.AddToEnd(3); m41023.Canonicalize();
			n->BConnection[0] = Connection(temp, m01234);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i-1], m1);
			n->BConnection[2] = Connection(temp, m21340);
			n->BConnection[3] = Connection(temp, m31402);
			n->BConnection[4] = Connection(temp, m41023);
			n->numExits = 5;	
		}
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);	
	}

	CFLOBDDNodeHandle MkU3GateInterleavedNode(unsigned int i, unsigned int eq_mark)
	{
		// Same as MkArbitraryGateInterleavedNode
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);

		CFLOBDDReturnMapHandle m01;
		m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
		CFLOBDDReturnMapHandle m23;
		m23.AddToEnd(2); m23.AddToEnd(3); m23.Canonicalize();
		CFLOBDDReturnMapHandle m20;
		m20.AddToEnd(2); m20.AddToEnd(0); m20.Canonicalize();
		CFLOBDDReturnMapHandle m02;
		m02.AddToEnd(0); m02.AddToEnd(2); m02.Canonicalize();
		CFLOBDDReturnMapHandle m0;
		m0.AddToEnd(0); m0.Canonicalize();
		CFLOBDDReturnMapHandle m1;
		m1.AddToEnd(1); m1.Canonicalize();
		CFLOBDDReturnMapHandle m2;
		m2.AddToEnd(2); m2.Canonicalize();
		CFLOBDDReturnMapHandle m12;
		m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
		CFLOBDDReturnMapHandle m21;
		m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
		CFLOBDDReturnMapHandle m10;
		m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
		

		n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
		n->numBConnections = 2;
		n->BConnection = new Connection[n->numBConnections];

		switch (eq_mark)
		{
		case 0b000000:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m23);//m23
			n->numExits = 4;
			break;
		case 0b001000:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m20);//m20
			n->numExits = 3;
			break;
		case 0b010000:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m02);//m02
			n->numExits = 3;
			break;
		case 0b011001:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);//m0
			n->numExits = 2;
			break;
		case 0b100000:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);//m0
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);//m12
			n->numExits = 3;
			break;
		case 0b101010:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);//m0
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m10);//m10
			n->numExits = 2;
			break;
		case 0b110100:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);//m0
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->numExits = 2;
			break;
		case 0b001100:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m10);//m10
			n->numExits = 2;
			break;
		case 0b000100:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);//m12
			n->numExits = 3;
			break;
		case 0b000111:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);//m1
			n->numExits = 2;
			break;
		case 0b000010:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m21);//m21
			n->numExits = 3;
			break;
		case 0b000001:
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m2);//m2
			n->numExits = 3;
			break;
		default:
			assert(false);
			break;
		}
		
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkArbitraryGateInterleavedNode(unsigned int i)
	{
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		CFLOBDDReturnMapHandle m01;
		m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
		n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01

		n->numBConnections = 2;
		n->BConnection = new Connection[n->numBConnections];
		n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
		CFLOBDDReturnMapHandle m23;
		m23.AddToEnd(2); m23.AddToEnd(3); m23.Canonicalize();
		n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m23);//m23
		n->numExits = 4;
		return CFLOBDDNodeHandle(n);
	}
	
	CFLOBDDNodeHandle MkCNOTInterleavedNode(unsigned int i)
	{
		assert(i >= 2);
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 2) {  // Base case

			CFLOBDDInternalNode *nh = new CFLOBDDInternalNode(i - 1);
			CFLOBDDReturnMapHandle m12;

			m12.AddToEnd(1);
			m12.AddToEnd(2);
			m12.Canonicalize();
			
			CFLOBDDNodeHandle temp = CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
			nh->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
			nh->numBConnections = 2;
			nh->BConnection = new Connection[nh->numBConnections];
			nh->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			nh->BConnection[1] = Connection(temp, m12);
			nh->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			nh->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nhHandle(nh);  // Create handle and canonicalize

			CFLOBDDReturnMapHandle m012;
			m012.AddToEnd(0);
			m012.AddToEnd(1);
			m012.AddToEnd(2);
			m012.Canonicalize();


			CFLOBDDInternalNode *nB = new CFLOBDDInternalNode(i - 1);
			nB->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
			nB->numBConnections = 2;
			nB->BConnection = new Connection[nB->numBConnections];
			nB->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			nB->BConnection[1] = Connection(temp, commonly_used_return_maps[3]);//m10
			nB->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
			nB->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nBHandle(nB);  // Create handle and canonicalize

			
			n->AConnection = Connection(nhHandle, m012);
			n->numBConnections = 3;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(nBHandle, commonly_used_return_maps[2]);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i - 1], commonly_used_return_maps[1]);//m1
			n->BConnection[2] = Connection(nBHandle, commonly_used_return_maps[3]);//m10
		}
		else {
			CFLOBDDNodeHandle temp = MkCNOTInterleavedNode(i - 1);
			n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i - 1], commonly_used_return_maps[1]);//m1
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkExchangeInterleavedNode(unsigned int level)
	{
		CFLOBDDInternalNode *n;

		if (level == 0) {
			return CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
		}
		else if (level == 1) {
			n = new CFLOBDDInternalNode(level);

			CFLOBDDNodeHandle temp = CFLOBDDNodeHandle::CFLOBDDForkNodeHandle;
			n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			n->BConnection[1] = Connection(temp, commonly_used_return_maps[3]);//m10
		}
		else {  // Create an appropriate CFLOBDDInternalNode
			n = new CFLOBDDInternalNode(level);

			CFLOBDDNodeHandle temp = MkExchangeInterleavedNode(level - 1);
			n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[0]);//m0
			n->BConnection[1] = Connection(temp, commonly_used_return_maps[2]);//m01
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}	// MkExchangeInterleavedNode

	CFLOBDDNodeHandle MkCNOTNode(unsigned int level, unsigned int n, long int controller, long int controlled)
	{
		assert(level >= 2);
		std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(controlled);
		if (cnot_hashMap.find(p) != cnot_hashMap.end()){
			return cnot_hashMap[p];
		}
		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(level);

		if (level == 2)
		{
			if (controller == -1 && controlled == -1)
			{
				return MkIdRelationInterleavedNode(level);
			}
			else if (controller != -1 && controlled != -1)
			{
				assert(controller == 0 && controlled == 1);
				return MkCNOTInterleavedNode(level);
			}
			else if (controller != -1)
			{
				CFLOBDDInternalNode *tmp = new CFLOBDDInternalNode(level - 1);
				tmp->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
				tmp->numBConnections = 2;
				tmp->BConnection = new Connection[2];
				tmp->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
				CFLOBDDReturnMapHandle m12;
				m12.AddToEnd(1);
				m12.AddToEnd(2);
				m12.Canonicalize();
				tmp->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);
				tmp->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
				tmp->InstallPathCounts();
#endif
				if (controller == 0)
				{
					CFLOBDDReturnMapHandle m012;
					m012.AddToEnd(0);
					m012.AddToEnd(1);
					m012.AddToEnd(2);
					m012.Canonicalize();
					CFLOBDDNodeHandle tmp_c = CFLOBDDNodeHandle(tmp);
					g->AConnection = Connection(tmp_c, m012);
					g->numBConnections = 3;
					g->BConnection = new Connection[3];
					CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
					g->BConnection[0] = Connection(Id, commonly_used_return_maps[2]);//m01
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[1]);//m1
					CFLOBDDReturnMapHandle m21;
					m21.AddToEnd(2);
					m21.AddToEnd(1);
					m21.Canonicalize();
					g->BConnection[2] = Connection(Id, m21);
					g->numExits = 3;
				}
				else
				{
					CFLOBDDNodeHandle Id =  MkIdRelationInterleavedNode(level - 1);
					g->AConnection = Connection(Id, commonly_used_return_maps[2]);//m01
					g->numBConnections = 2;
					g->BConnection = new Connection[2];
					CFLOBDDReturnMapHandle m012;
					m012.AddToEnd(0);
					m012.AddToEnd(1);
					m012.AddToEnd(2);
					m012.Canonicalize();
					g->BConnection[0] = Connection(tmp, m012);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[1]);//m1
					g->numExits = 3;
				}
			}
			else if (controlled != -1)
			{
				if (controlled == 0)
				{
					CFLOBDDNodeHandle Id =  MkIdRelationInterleavedNode(level - 1);
					g->AConnection = Connection(Id, commonly_used_return_maps[2]);//m01
					g->numBConnections = 2;
					g->BConnection = new Connection[2];
					g->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[0]);//m0
					g->BConnection[1] = Connection(Id, commonly_used_return_maps[3]);//m10
					g->numExits = 2;
				}
				else
				{
					CFLOBDDNodeHandle Id =  MkIdRelationInterleavedNode(level - 1);
					g->AConnection = Connection(Id, commonly_used_return_maps[2]);//m01
					g->numBConnections = 2;
					g->BConnection = new Connection[2];
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[0]);//m0
					g->BConnection[0] = Connection(Id, commonly_used_return_maps[2]);//m01
					g->numExits = 2;
				}
			}
		}
		else
		{
			long int numVarsOfLowerLevel = (1 << (level - 2));
			long int controller_A = (controller >= numVarsOfLowerLevel) ? -1 : controller;
			long int controlled_A = (controlled >= numVarsOfLowerLevel) ? -1 : controlled;
			CFLOBDDNodeHandle ANodeHandle;
			if (controller_A == -1 && controlled_A == -1)
				ANodeHandle = MkIdRelationInterleavedNode(level - 1);
			else
				ANodeHandle = MkCNOTNode(level - 1, n, controller_A, controlled_A);
			g->numBConnections = ANodeHandle.handleContents->numExits;
			CFLOBDDReturnMapHandle mI;
			for (int i = 0; i < g->numBConnections; i++)
				mI.AddToEnd(i);
			mI.Canonicalize();
			g->AConnection = Connection(ANodeHandle, mI);
			g->BConnection = new Connection[g->numBConnections];

			long int controller_B = (controller >= numVarsOfLowerLevel) ? (controller - numVarsOfLowerLevel) : -1;
			long int controlled_B = (controlled >= numVarsOfLowerLevel) ? (controlled - numVarsOfLowerLevel) : -1;
			if (g->numBConnections == 2)
			{
				CFLOBDDNodeHandle B0;
				if (controller_B == -1 && controlled_B == -1){
					if (controlled_A != -1)
						B0 = CFLOBDDNodeHandle::NoDistinctionNode[level - 1];
					else
						B0 = MkIdRelationInterleavedNode(level - 1);
				}
				else
					B0 = MkCNOTNode(level - 1, n, controller_B, controlled_B);
				CFLOBDDReturnMapHandle mI_B;
				for (int i = 0; i < B0.handleContents->numExits; i++)
					mI_B.AddToEnd(i);
				mI_B.Canonicalize();
				g->BConnection[0] = Connection(B0, mI_B);
				// NoDistinct Node
				if (controlled_A != -1){
					CFLOBDDNodeHandle Id =  MkIdRelationInterleavedNode(level - 1);
					g->BConnection[1] = Connection(Id, commonly_used_return_maps[3]);//m10
				}
				else{
					if (controlled_B == -1){
						g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[1]);//m1
					}
					else{
						g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[0]);//m0
					}
				}
				g->numExits = mI_B.Size() > 2 ? mI_B.Size() : 2;
			}
			else
			{
				// ID Connection
				CFLOBDDNodeHandle B0 = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(B0, commonly_used_return_maps[2]);//m01
				// NoDistinct Node
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], commonly_used_return_maps[1]);//m1

				CFLOBDDNodeHandle B2;
				if (controlled == -1){
					B2 = MkIdRelationInterleavedNode(level - 1);
					CFLOBDDReturnMapHandle m21;
					m21.AddToEnd(2);
					m21.AddToEnd(1);
					m21.Canonicalize();
					g->BConnection[2] = Connection(B2, m21);
					g->numExits = 3;
				}
				else{
					B2 = MkCNOTNode(level - 1, n, controller_B, controlled_B);
					g->BConnection[2] = Connection(B2, commonly_used_return_maps[3]);//m10
					g->numExits = 2;
				}
			}
		}

#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle = CFLOBDDNodeHandle(g);
		cnot_hashMap.insert(std::make_pair(p, gHandle));
		return gHandle;
	}


	CFLOBDDNodeHandle MkCNOT2Node(unsigned int level, unsigned int n, long int controller, long int controlled)
	{
		std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(controlled);
		if (cnot_hashMap.find(p) != cnot_hashMap.end()){
			return cnot_hashMap[p];
		}	

		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(level);

		if (level == 1)
		{
			if (controller == 0)
			{
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[g->numBConnections];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				CFLOBDDReturnMapHandle m12;
				m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);
				g->numExits = 3;
			}
			else if (controlled == 0)
			{
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[g->numBConnections];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				CFLOBDDReturnMapHandle m10;
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m10);
				g->numExits = 2;
			}
		}
		else
		{
			if (controller < n/2 && controlled < n/2 && controlled >= 0 && controller >= 0)
			{
				// Case 1: Both in A Connection
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto aa = MkCNOT2Node(level-1, n/2, controller, controlled);
				g->AConnection = Connection(aa, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller >= n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 2: Both in B Connection
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				auto bb = MkCNOT2Node(level-1, n/2, controller - n/2, controlled - n/2);
				g->BConnection[0] = Connection(bb, m01);
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller < n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 3: controller in A and controlled in B
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				auto aa = MkCNOT2Node(level-1, n/2, controller, -1);
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDReturnMapHandle m01, m1, m10;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				auto bb = MkCNOT2Node(level-1, n/2, -1, controlled - n/2);
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(bb, m10);
				g->numExits = 2;
			}
			else if (controlled == -1 && controller < n/2 && controller >= 0)
			{
				// Case 4: controller in A and controlled == -1
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				auto aa = MkCNOT2Node(level-1, n/2, controller, -1);
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[g->numBConnections];
				CFLOBDDReturnMapHandle m01, m1, m21;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else if (controlled == -1 && controller >= n/2)
			{
				// Case 5: controller in B and controlled == -1
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				auto bb = MkCNOT2Node(level-1, n/2, controller - n/2, -1);
				CFLOBDDReturnMapHandle m012, m1;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[0] = Connection(bb, m012);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 3;
			}
			else if (controller == -1 && controlled >= 0 && controlled < n/2)
			{
				// Case 6: controller == -1 and controlled in A
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto aa = MkCNOT2Node(level-1, n/2, -1, controlled);
				g->AConnection = Connection(aa, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m0, m10;
				m0.AddToEnd(0); m0.Canonicalize();
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[1] = Connection(Id, m10);
				g->numExits = 2;
			}
			else if (controller == -1 && controlled >= n/2)
			{
				// Case 7: controller == -1 and controlled in B
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				auto bb = MkCNOT2Node(level-1, n/2, -1, controlled - n/2);
				CFLOBDDReturnMapHandle m0;
				m0.AddToEnd(0); m0.Canonicalize();
				g->BConnection[0] = Connection(bb, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
				g->numExits = 2;
			}
		}

#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle = CFLOBDDNodeHandle(g);
		cnot_hashMap.insert(std::make_pair(p, gHandle));
		return gHandle;
	}

	CFLOBDDNodeHandle MkCCNOTNode(unsigned int level, unsigned int n, long int controller1, long int controller2, long int controlled)
	{
		std::string p = std::to_string(level) + ";" + std::to_string(controller1) + ";" + std::to_string(controller2) + ";" + std::to_string(controlled);
		if (ccnot_hashMap.find(p) != ccnot_hashMap.end()){
			return ccnot_hashMap[p];
		}	

		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(level);

		if (level == 1)
		{
			if (controller1 == 0 || controller2 == 0)
			{
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[g->numBConnections];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				CFLOBDDReturnMapHandle m12;
				m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);
				g->numExits = 3;
			}
			else if (controlled == 0)
			{
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[g->numBConnections];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				CFLOBDDReturnMapHandle m10;
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m10);
				g->numExits = 2;
			}
		}
		else
		{
			if (controller1 < n/2 && controller2 < n/2 && controlled < n/2 && controlled >= 0 && controller1 >= 0 && controller2 >= 0)
			{
				// Case 1: CR1, CR2 and CD in A Connection
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto aa = MkCCNOTNode(level-1, n/2, controller1, controller2, controlled);
				g->AConnection = Connection(aa, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller1 >= n/2 && controller2 >= n/2 && controlled >= n/2 && controller1 >= 0 && controller2 >= 0 && controlled >= 0)
			{
				// Case 2: CR1, CR2 and CD in B Connection
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				auto bb = MkCCNOTNode(level-1, n/2, controller1 - n/2, controller2 - n/2, controlled - n/2);
				g->BConnection[0] = Connection(bb, m01);
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller1 < n/2 && controller2 >= n/2 && controlled >= n/2 && controller1 >= 0 && controlled >= 0 && controller2 >= 0)
			{
				// Case 3: CR1 in A, CR2 and CD in B
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				auto aa = MkCCNOTNode(level-1, n/2, controller1, -1, -1);
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDReturnMapHandle m01, m1, m10;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				auto bb = MkCCNOTNode(level-1, n/2, -1, controller2 - n/2, controlled - n/2);
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(bb, m01);
				g->numExits = 2;
			}
			else if (controlled == -1 && controller1 < n/2 && controller1 >= 0 && controller2 == -1)
			{
				// Case 4: CR1 in A and CR2 == -1 and CD == -1
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				auto aa = MkCCNOTNode(level-1, n/2, controller1, -1, -1);
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[g->numBConnections];
				CFLOBDDReturnMapHandle m01, m1, m21;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else if (controlled == -1 && controller1 >= n/2 && controller2 == -1)
			{
				// Case 5: CR1 in B and CR2 == -1 and CD == -1
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				auto bb = MkCCNOTNode(level-1, n/2, controller1 - n/2, -1, -1);
				CFLOBDDReturnMapHandle m012, m1;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[0] = Connection(bb, m012);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 3;
			}
			else if (controller1 == -1 && controlled >= 0 && controlled < n/2 && controller2 == -1)
			{
				// Case 6: controller == -1 and controlled in A
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto aa = MkCCNOTNode(level-1, n/2, -1, -1, controlled);
				g->AConnection = Connection(aa, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m0, m10;
				m0.AddToEnd(0); m0.Canonicalize();
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[1] = Connection(Id, m10);
				g->numExits = 2;
			}
			else if (controller1 == -1 && controlled >= n/2 && controller2 == -1)
			{
				// Case 7: controller == -1 and controlled in B
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				auto bb = MkCCNOTNode(level-1, n/2, -1, -1, controlled - n/2);
				CFLOBDDReturnMapHandle m0;
				m0.AddToEnd(0); m0.Canonicalize();
				g->BConnection[0] = Connection(bb, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
				g->numExits = 2;
			}
			else if (controller1 < n/2 && controller2 < n/2 && controlled >= n/2 && controller1 >= 0 && controlled >= 0 && controller2 >= 0)
			{
				// Case 8: CR1 and CR2 in A, CD in B
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				auto aa = MkCCNOTNode(level-1, n/2, controller1, controller2, -1);
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDReturnMapHandle m01, m1, m10;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				auto bb = MkCCNOTNode(level-1, n/2, -1, -1, controlled - n/2);
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(bb, m10);
				g->numExits = 2;
			}
			else if (controller1 < n/2 && controller2 < n/2 && controlled == -1 && controller1 >= 0 && controller2 >= 0)
			{
				// Case 9: CR1 and CR2 in A, CD == -1
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				auto aa = MkCCNOTNode(level-1, n/2, controller1, controller2, -1);
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDReturnMapHandle m01, m1, m21;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else if (controller1 >= n/2 && controller2 >= n/2 && controlled == -1 && controller1 >= 0 && controller2 >= 0)
			{
				// Case 10: CR1 and CR2 in B, CD == -1
				auto Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m012, m1;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				auto bb = MkCCNOTNode(level-1, n/2, controller1 - n/2, controller2 - n/2, -1);
				g->BConnection[0] = Connection(bb, m012);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 3;
			}
			else if (controller1 < n/2 && controller2 >= n/2 && controlled == -1 && controller1 >= 0 && controller2 >= 0)
			{
				// Case 11: CR1 in A, CR2 in B, CD == -1
				auto Id = MkIdRelationInterleavedNode(level - 1);
				auto aa = MkCCNOTNode(level-1, n/2, controller1, -1, -1);
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDReturnMapHandle m01, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				auto bb = MkCCNOTNode(level-1, n/2, -1, controller2 - n/2, -1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(bb, m012);
				g->numExits = 3;
			}
			else if (controller1 == -1 && controller2 < n/2 && controlled == -1 && controller2 >= 0)
			{
				// Case 12: CR1 == -1, CR2 in A, CD == -1
				auto Id = MkIdRelationInterleavedNode(level - 1);
				auto aa = MkCCNOTNode(level-1, n/2, -1, controller2, -1);
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDReturnMapHandle m01, m1, m21;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else if (controller1 == -1 && controller2 >= n/2 && controlled == -1 && controller2 >= 0)
			{
				// Case 12: CR1 == -1, CR2 in B, CD == -1
				auto Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m012, m1;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				auto bb = MkCCNOTNode(level-1, n/2, -1, controller2 - n/2, -1);
				g->BConnection[0] = Connection(bb, m012);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 3;
			}
			else if (controller1 == -1 && controller2 < n/2 && controlled < n/2 && controller2 >= 0 && controlled >= 0)
			{
				// Case 13: CR1 == -1, CR2 and CD in A
				auto Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto aa = MkCCNOTNode(level-1, n/2, -1, controller2, controlled);
				g->AConnection = Connection(aa, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller1 == -1 && controller2 >= n/2 && controlled >= n/2 && controller2 >= 0 && controlled >= 0)
			{
				// Case 14: CR1 == -1, CR2 and CD in B
				auto Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				auto bb = MkCCNOTNode(level-1, n/2, -1, controller2 - n/2, controlled - n/2);
				g->BConnection[0] = Connection(bb, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller1 == -1 && controller2 < n/2 && controlled >= n/2 && controller2 >= 0 && controlled >= 0)
			{
				// Case 15: CR1 == -1, CR2 in A and CD in B
				auto Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				auto aa = MkCCNOTNode(level-1, n/2, -1, controller2, -1);
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDReturnMapHandle m01, m1, m10;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				auto bb = MkCCNOTNode(level-1, n/2, -1, -1, controlled - n/2);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(bb, m10);
				g->numExits = 2;
			}
		}

#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle = CFLOBDDNodeHandle(g);
		ccnot_hashMap.insert(std::make_pair(p, gHandle));
		return gHandle;
	}

	CFLOBDDNodeHandle MkCPGateNode(unsigned int level, long int controller, long int controlled)
	{
		std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(controlled);
		if (cp_hashMap.find(p) != cp_hashMap.end()){
			return cp_hashMap[p];
		}
		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(level);
		if (level == 1)
		{
			// std::cout << "CaseBase" << std::endl;
			CFLOBDDReturnMapHandle m01, m12;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
			g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
			g->numBConnections = 2;
			g->BConnection = new Connection[2];
			g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
			g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);
			g->numExits = 3;
		}
		else 
		{
			unsigned int n = pow(2, level - 1);
			if (controller < n/2 && controlled < n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 1: Both fall in A
				// std::cout << "Case1" << std::endl;
				CFLOBDDNodeHandle aTmp = MkCPGateNode(level-1, controller, controlled);
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0);
				m012.AddToEnd(1);
				m012.AddToEnd(2);
				m012.Canonicalize();
				g->AConnection = Connection(aTmp, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDReturnMapHandle m01, m1, m21;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else if (controller >= n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 2: Both fall in B region
				// std::cout << "Case2" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				CFLOBDDNodeHandle btmp = MkCPGateNode(level-1, controller - n/2, controlled - n/2);
				g->BConnection[0] = Connection(btmp, m012);
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 3;
			}
			else if (controller < n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 3: controller in A and controlled in B
				// std::cout << "Case3" << std::endl;
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0);
				m012.AddToEnd(1);
				m012.AddToEnd(2);
				m012.Canonicalize();
				CFLOBDDNodeHandle atmp = MkCPGateNode(level-1, controller, -1);
				g->AConnection = Connection(atmp, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDReturnMapHandle m01, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				CFLOBDDNodeHandle btmp = MkCPGateNode(level-1, -1, controlled - n/2);
				g->BConnection[2] = Connection(btmp, m012);
				g->numExits = 3;
			}
			else if (controller < n/2 && controlled == -1 && controller >= 0)
			{
				// Case 4: controller in A and controlled == -1
				// std::cout << "Case4" << std::endl;
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0);
				m012.AddToEnd(1);
				m012.AddToEnd(2);
				m012.Canonicalize();
				CFLOBDDNodeHandle atmp = MkCPGateNode(level-1, controller, -1);
				g->AConnection = Connection(atmp, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				CFLOBDDReturnMapHandle m01, m21, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else if (controller >= n/2 && controlled == -1 && controller >= 0)
			{
				// Case 5: controller in B and controlled == -1
				// std::cout << "Case5" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				CFLOBDDNodeHandle btmp = MkCPGateNode(level-1, controller - n/2, -1);
				CFLOBDDReturnMapHandle m012, m1;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(btmp, m012);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 3;
			} 
			else if (controller == -1 && controlled < n/2 && controlled >= 0)
			{
				// Case 6: controller == -1 && controlled in A
				// std::cout << "Case6" << std::endl;
				CFLOBDDNodeHandle atmp = MkCPGateNode(level-1, -1, controlled);
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				g->AConnection = Connection(atmp, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01, m21, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else if (controller == -1 && controlled >= n/2 && controlled >= 0)
			{
				// Case 7: controller == -1 && controlled in B
				// std::cout << "Case7" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m012, m1;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				CFLOBDDNodeHandle btmp = MkCPGateNode(level-1, -1, controlled - n/2);
				g->BConnection[0] = Connection(btmp, m012);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 3;
			}
			
		}
		#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle = CFLOBDDNodeHandle(g);
		cp_hashMap.insert(std::make_pair(p, gHandle));
		return gHandle;
	}

	CFLOBDDNodeHandle MkSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num)
	{
		std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(controlled) + ";" + std::to_string(case_num);
		if (swap_hashMap.find(p) != swap_hashMap.end()){
			return swap_hashMap[p];
		}
		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(level);
		if (level == 1)
		{
			// std::cout << "CaseBase" << std::endl;
			if (case_num == -1)
			{
				assert(controlled == -1);
				CFLOBDDReturnMapHandle m01, m23;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m23.AddToEnd(2); m23.AddToEnd(3); m23.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m23);
				g->numExits = 4;
			}
			else if (case_num == 0)
			{
				// [[1 0] [0 0]]
				assert(controller == -1);
				CFLOBDDReturnMapHandle m01, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);
				g->numExits = 2;
			}
			else if (case_num == 1)
			{
				// [[0 1] [0 0]]
				assert(controller == -1);
				CFLOBDDReturnMapHandle m01, m0;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				g->numExits = 2;	
			}
			else if (case_num == 2)
			{
				// [[0 0] [1 0]]
				assert(controller == -1);
				CFLOBDDReturnMapHandle m01, m0, m10;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m10);
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				g->numExits = 2;	
			}
			else if (case_num == 3)
			{
				// [[0 0] [0 1]]
				assert(controller == -1);
				CFLOBDDReturnMapHandle m01, m0;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				g->numExits = 2;	
			}
		}
		else if (level == 2 && controller == 0 && controlled == 1)
		{
			CFLOBDDNodeHandle atmp = MkSwapGateNode(level - 1, controller, -1, case_num);
			CFLOBDDReturnMapHandle m0123;
			m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
			g->AConnection = Connection(atmp, m0123);
			CFLOBDDNodeHandle b0 = MkSwapGateNode(level-1, -1, controlled, 0);
			CFLOBDDNodeHandle b1 = MkSwapGateNode(level-1, -1, controlled, 1);
			CFLOBDDNodeHandle b2 = MkSwapGateNode(level-1, -1, controlled, 2);
			CFLOBDDNodeHandle b3 = MkSwapGateNode(level-1, -1, controlled, 3);
			CFLOBDDReturnMapHandle m01, m10;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
			g->numBConnections = 4;
			g->BConnection = new Connection[4];
			g->BConnection[0] = Connection(b0, m01);
			g->BConnection[1] = Connection(b2, m10);
			g->BConnection[2] = Connection(b1, m10);
			g->BConnection[3] = Connection(b3, m10);
			g->numExits = 2;
		}
		else if (level == 2 && controller == 0 && controlled == -1)
		{
			CFLOBDDNodeHandle atmp = MkSwapGateNode(level - 1, controller, -1, case_num);
			CFLOBDDReturnMapHandle m0123;
			m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
			g->AConnection = Connection(atmp, m0123);
			CFLOBDDReturnMapHandle m01, m21, m31, m41;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
			m31.AddToEnd(3); m31.AddToEnd(1); m31.Canonicalize();
			m41.AddToEnd(4); m41.AddToEnd(1); m41.Canonicalize();
			g->numBConnections = 4;
			g->BConnection = new Connection[4];
			CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
			g->BConnection[0] = Connection(Id, m01);
			g->BConnection[1] = Connection(Id, m21);
			g->BConnection[2] = Connection(Id, m31);
			g->BConnection[3] = Connection(Id, m41);
			g->numExits = 5;	
		}
		else if (level == 2 && controller == 1 && controlled == -1)
		{
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
			g->AConnection = Connection(Id, m01);
			g->numBConnections = 2;
			CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, 0, -1, case_num);
			CFLOBDDReturnMapHandle m0123, m4;
			m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
			m4.AddToEnd(4); m4.Canonicalize();
			g->BConnection = new Connection[2];
			g->BConnection[0] = Connection(btmp, m0123);
			g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m4);
			g->numExits = 5;

		}
		else 
		{
			unsigned int n = pow(2, level - 1);
			if (controller < n/2 && controlled < n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 1: Both fall in A
				// std::cout << "Case1" << std::endl;
				CFLOBDDNodeHandle aTmp = MkSwapGateNode(level-1, controller, controlled, case_num);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m01.Canonicalize();
				g->AConnection = Connection(aTmp, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller >= n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 2: Both fall in B region
				// std::cout << "Case2" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, controller - n/2, controlled - n/2, case_num);
				g->BConnection[0] = Connection(btmp, m01);
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller < n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 3: controller in A and controlled in B
				// std::cout << "Case3" << std::endl;
				CFLOBDDReturnMapHandle m01234;
				m01234.AddToEnd(0);
				m01234.AddToEnd(1);
				m01234.AddToEnd(2);
				m01234.AddToEnd(3);
				m01234.AddToEnd(4);
				m01234.Canonicalize();
				CFLOBDDNodeHandle atmp = MkSwapGateNode(level-1, controller, -1, case_num);
				g->AConnection = Connection(atmp, m01234);
				g->numBConnections = 5;
				g->BConnection = new Connection[5];
				CFLOBDDNodeHandle b0 = MkSwapGateNode(level-1, -1, controlled - n/2, 0);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->BConnection[0] = Connection(b0, m01);
				CFLOBDDNodeHandle b1 = MkSwapGateNode(level-1, -1, controlled - n/2, 2);
				CFLOBDDNodeHandle b2 = MkSwapGateNode(level-1, -1, controlled - n/2, 1);
				CFLOBDDNodeHandle b3 = MkSwapGateNode(level-1, -1, controlled - n/2, 3);
				CFLOBDDReturnMapHandle m10, m1;
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				if (controller == n/2 - 1)
				{
					// std::cout << b1 << std::endl;
					// std::cout << b2 << std::endl;
					// std::cout << b3 << std::endl;
					g->BConnection[1] = Connection(b1, m10);
					g->BConnection[2] = Connection(b2, m10);
					g->BConnection[3] = Connection(b3, m10);
					g->BConnection[4] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else
				{
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
					g->BConnection[2] = Connection(b1, m10);
					g->BConnection[3] = Connection(b2, m10);
					g->BConnection[4] = Connection(b3, m10);
				}
				g->numExits = 2;
			}
			else if (controller < n/2 && controlled == -1 && controller >= 0)
			{
				// Case 4: controller in A and controlled == -1
				// std::cout << "Case4" << std::endl;
				CFLOBDDReturnMapHandle m01234;
				m01234.AddToEnd(0);
				m01234.AddToEnd(1);
				m01234.AddToEnd(2);
				m01234.AddToEnd(3);
				m01234.AddToEnd(4);
				m01234.Canonicalize();
				CFLOBDDNodeHandle atmp = MkSwapGateNode(level-1, controller, -1, case_num);
				g->AConnection = Connection(atmp, m01234);
				g->numBConnections = 5;
				g->BConnection = new Connection[5];
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				CFLOBDDReturnMapHandle m01, m21, m31, m41, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				m31.AddToEnd(3); m31.AddToEnd(1); m31.Canonicalize();
				m41.AddToEnd(4); m41.AddToEnd(1); m41.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				if (controller == n/2 - 1)
				{
					g->BConnection[0] = Connection(Id, m01);
					g->BConnection[1] = Connection(Id, m21);
					g->BConnection[2] = Connection(Id, m31);
					g->BConnection[3] = Connection(Id, m41);
					g->BConnection[4] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else
				{
					g->BConnection[0] = Connection(Id, m01);
					g->BConnection[2] = Connection(Id, m21);
					g->BConnection[3] = Connection(Id, m31);
					g->BConnection[4] = Connection(Id, m41);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);	
				}
				g->numExits = 5;
			}
			else if (controller >= n/2 && controlled == -1 && controller >= 0)
			{
				// Case 5: controller in B and controlled == -1
				// std::cout << "Case5" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, controller - n/2, -1, case_num);
				CFLOBDDReturnMapHandle m0123, m1, m2, m4, m01234;
				m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
				m01234.AddToEnd(0); m01234.AddToEnd(1); m01234.AddToEnd(2); m01234.AddToEnd(3); m01234.AddToEnd(4); m01234.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m2.AddToEnd(2); m2.Canonicalize();
				m4.AddToEnd(4); m4.Canonicalize();
				g->BConnection = new Connection[2];
				// if (controller == n - 1)
				// {
				// 	g->BConnection[0] = Connection(btmp, m0123);
				// 	g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m4);
				// }
				// else
				// {
					g->BConnection[0] = Connection(btmp, m01234);
					if (controller == n - 1)
						g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m4);
					else
						g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);

				// }
				g->numExits = 5;
			} 
			else if (controller == -1 && controlled < n/2 && controlled >= 0)
			{
				// Case 6: controller == -1 && controlled in A
				// std::cout << "Case6" << std::endl;
				CFLOBDDNodeHandle atmp = MkSwapGateNode(level-1, -1, controlled, case_num);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(atmp, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m1, m0;
				m1.AddToEnd(1); m1.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				if (case_num == 0)
				{
					g->BConnection[0] = Connection(Id, m01);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else{
					CFLOBDDReturnMapHandle m10;
					m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
					g->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
					g->BConnection[1] = Connection(Id, m10);
				}
				g->numExits = 2;
			}
			else if (controller == -1 && controlled >= n/2 && controlled >= 0)
			{
				// Case 7: controller == -1 && controlled in B
				// std::cout << "Case7" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m1, m0;
				m1.AddToEnd(1); m1.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, -1, controlled - n/2, case_num);
				if (case_num == 0)
				{
					g->BConnection[0] = Connection(btmp, m01);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else
				{
					g->BConnection[0] = Connection(btmp, m01);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
				}
				g->numExits = 2;
			}
			
		}
		#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle = CFLOBDDNodeHandle(g);
		swap_hashMap.insert(std::make_pair(p, gHandle));
		return gHandle;
	}

	
	CFLOBDDNodeHandle MkiSwapGateNode(unsigned int level, long int controller, long int controlled, int case_num)
	{
		std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(controlled) + ";" + std::to_string(case_num);
		if (iswap_hashMap.find(p) != iswap_hashMap.end()){
			return iswap_hashMap[p];
		}
		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(level);
		if (level == 1)
		{
			// std::cout << "CaseBase" << std::endl;
			if (case_num == -1)
			{
				assert(controlled == -1);
				CFLOBDDReturnMapHandle m01, m23;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m23.AddToEnd(2); m23.AddToEnd(3); m23.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m23);
				g->numExits = 4;
			}
			else if (case_num == 0)
			{
				// [[1 0] [0 0]]
				assert(controller == -1);
				CFLOBDDReturnMapHandle m01, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);
				g->numExits = 2;
			}
			else if (case_num == 1)
			{
				// [[0 i] [0 0]]
				assert(controller == -1);
				CFLOBDDReturnMapHandle m01, m0;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				g->numExits = 2;	
			}
			else if (case_num == 2)
			{
				// [[0 0] [i 0]]
				assert(controller == -1);
				CFLOBDDReturnMapHandle m01, m0, m10;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m10);
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				g->numExits = 2;	
			}
			else if (case_num == 3)
			{
				// [[0 0] [0 1]]
				assert(controller == -1);
				CFLOBDDReturnMapHandle m01, m0;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				g->numExits = 2;	
			}
		}
		else if (level == 2 && controller == 0 && controlled == 1)
		{
			CFLOBDDNodeHandle atmp = MkSwapGateNode(level - 1, controller, -1, case_num);
			CFLOBDDReturnMapHandle m0123;
			m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
			g->AConnection = Connection(atmp, m0123);
			CFLOBDDNodeHandle b0 = MkSwapGateNode(level-1, -1, controlled, 0);
			CFLOBDDNodeHandle b1 = MkSwapGateNode(level-1, -1, controlled, 1);
			CFLOBDDNodeHandle b2 = MkSwapGateNode(level-1, -1, controlled, 2);
			CFLOBDDNodeHandle b3 = MkSwapGateNode(level-1, -1, controlled, 3);
			CFLOBDDReturnMapHandle m01, m10, m12;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
			m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
			g->numBConnections = 4;
			g->BConnection = new Connection[4];
			g->BConnection[0] = Connection(b0, m01);
			g->BConnection[1] = Connection(b2, m12);
			g->BConnection[2] = Connection(b1, m12);
			g->BConnection[3] = Connection(b3, m10);
			g->numExits = 2;
		}
		else if (level == 2 && controller == 0 && controlled == -1)
		{
			CFLOBDDNodeHandle atmp = MkSwapGateNode(level - 1, controller, -1, case_num);
			CFLOBDDReturnMapHandle m0123;
			m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
			g->AConnection = Connection(atmp, m0123);
			CFLOBDDReturnMapHandle m01, m21, m31, m41;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
			m31.AddToEnd(3); m31.AddToEnd(1); m31.Canonicalize();
			m41.AddToEnd(4); m41.AddToEnd(1); m41.Canonicalize();
			g->numBConnections = 4;
			g->BConnection = new Connection[4];
			CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
			g->BConnection[0] = Connection(Id, m01);
			g->BConnection[1] = Connection(Id, m21);
			g->BConnection[2] = Connection(Id, m31);
			g->BConnection[3] = Connection(Id, m41);
			g->numExits = 5;	
		}
		else if (level == 2 && controller == 1 && controlled == -1)
		{
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
			g->AConnection = Connection(Id, m01);
			g->numBConnections = 2;
			CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, 0, -1, case_num);
			CFLOBDDReturnMapHandle m0123, m4;
			m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
			m4.AddToEnd(4); m4.Canonicalize();
			g->BConnection = new Connection[2];
			g->BConnection[0] = Connection(btmp, m0123);
			g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m4);
			g->numExits = 5;

		}
		else 
		{
			unsigned int n = pow(2, level - 1);
			if (controller < n/2 && controlled < n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 1: Both fall in A
				// std::cout << "Case1" << std::endl;
				CFLOBDDNodeHandle aTmp = MkSwapGateNode(level-1, controller, controlled, case_num);
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0);
				m012.AddToEnd(1);
				m012.AddToEnd(2);
				m012.Canonicalize();
				g->AConnection = Connection(aTmp, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01, m21;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else if (controller >= n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 2: Both fall in B region
				// std::cout << "Case2" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, controller - n/2, controlled - n/2, case_num);
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				g->BConnection[0] = Connection(btmp, m012);
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 3;
			}
			else if (controller < n/2 && controlled >= n/2 && controller >= 0 && controlled >= 0)
			{
				// Case 3: controller in A and controlled in B
				// std::cout << "Case3" << std::endl;
				CFLOBDDReturnMapHandle m01234;
				m01234.AddToEnd(0);
				m01234.AddToEnd(1);
				m01234.AddToEnd(2);
				m01234.AddToEnd(3);
				m01234.AddToEnd(4);
				m01234.Canonicalize();
				CFLOBDDNodeHandle atmp = MkSwapGateNode(level-1, controller, -1, case_num);
				g->AConnection = Connection(atmp, m01234);
				g->numBConnections = 5;
				g->BConnection = new Connection[5];
				CFLOBDDNodeHandle b0 = MkSwapGateNode(level-1, -1, controlled - n/2, 0);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->BConnection[0] = Connection(b0, m01);
				CFLOBDDNodeHandle b1 = MkSwapGateNode(level-1, -1, controlled - n/2, 2);
				CFLOBDDNodeHandle b2 = MkSwapGateNode(level-1, -1, controlled - n/2, 1);
				CFLOBDDNodeHandle b3 = MkSwapGateNode(level-1, -1, controlled - n/2, 3);
				CFLOBDDReturnMapHandle m10, m1, m12;
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
				if (controller == n/2 - 1)
				{
					// std::cout << b1 << std::endl;
					// std::cout << b2 << std::endl;
					// std::cout << b3 << std::endl;
					g->BConnection[1] = Connection(b1, m12);
					g->BConnection[2] = Connection(b2, m12);
					g->BConnection[3] = Connection(b3, m10);
					g->BConnection[4] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else
				{
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
					g->BConnection[2] = Connection(b1, m12);
					g->BConnection[3] = Connection(b2, m12);
					g->BConnection[4] = Connection(b3, m10);
				}
				g->numExits = 3;
			}
			else if (controller < n/2 && controlled == -1 && controller >= 0)
			{
				// Case 4: controller in A and controlled == -1
				// std::cout << "Case4" << std::endl;
				CFLOBDDReturnMapHandle m01234;
				m01234.AddToEnd(0);
				m01234.AddToEnd(1);
				m01234.AddToEnd(2);
				m01234.AddToEnd(3);
				m01234.AddToEnd(4);
				m01234.Canonicalize();
				CFLOBDDNodeHandle atmp = MkSwapGateNode(level-1, controller, -1, case_num);
				g->AConnection = Connection(atmp, m01234);
				g->numBConnections = 5;
				g->BConnection = new Connection[5];
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				CFLOBDDReturnMapHandle m01, m21, m31, m41, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				m31.AddToEnd(3); m31.AddToEnd(1); m31.Canonicalize();
				m41.AddToEnd(4); m41.AddToEnd(1); m41.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				if (controller == n/2 - 1)
				{
					g->BConnection[0] = Connection(Id, m01);
					g->BConnection[1] = Connection(Id, m21);
					g->BConnection[2] = Connection(Id, m31);
					g->BConnection[3] = Connection(Id, m41);
					g->BConnection[4] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else
				{
					g->BConnection[0] = Connection(Id, m01);
					g->BConnection[2] = Connection(Id, m21);
					g->BConnection[3] = Connection(Id, m31);
					g->BConnection[4] = Connection(Id, m41);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);	
				}
				g->numExits = 5;
			}
			else if (controller >= n/2 && controlled == -1 && controller >= 0)
			{
				// Case 5: controller in B and controlled == -1
				// std::cout << "Case5" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, controller - n/2, -1, case_num);
				CFLOBDDReturnMapHandle m0123, m1, m2, m4, m01234;
				m0123.AddToEnd(0); m0123.AddToEnd(1); m0123.AddToEnd(2); m0123.AddToEnd(3); m0123.Canonicalize();
				m01234.AddToEnd(0); m01234.AddToEnd(1); m01234.AddToEnd(2); m01234.AddToEnd(3); m01234.AddToEnd(4); m01234.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m2.AddToEnd(2); m2.Canonicalize();
				m4.AddToEnd(4); m4.Canonicalize();
				g->BConnection = new Connection[2];
				// if (controller == n - 1)
				// {
				// 	g->BConnection[0] = Connection(btmp, m0123);
				// 	g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m4);
				// }
				// else
				// {
					g->BConnection[0] = Connection(btmp, m01234);
					if (controller == n - 1)
						g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m4);
					else
						g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);

				// }
				g->numExits = 5;
			} 
			else if (controller == -1 && controlled < n/2 && controlled >= 0)
			{
				// Case 6: controller == -1 && controlled in A
				// std::cout << "Case6" << std::endl;
				CFLOBDDNodeHandle atmp = MkSwapGateNode(level-1, -1, controlled, case_num);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(atmp, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m1, m0;
				m1.AddToEnd(1); m1.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				if (case_num == 0)
				{
					g->BConnection[0] = Connection(Id, m01);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else{
					CFLOBDDReturnMapHandle m10;
					m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
					g->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
					g->BConnection[1] = Connection(Id, m10);
				}
				g->numExits = 2;
			}
			else if (controller == -1 && controlled >= n/2 && controlled >= 0)
			{
				// Case 7: controller == -1 && controlled in B
				// std::cout << "Case7" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m1, m0;
				m1.AddToEnd(1); m1.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, -1, controlled - n/2, case_num);
				if (case_num == 0)
				{
					g->BConnection[0] = Connection(btmp, m01);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else
				{
					g->BConnection[0] = Connection(btmp, m01);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
				}
				g->numExits = 2;
			}
			
		}
		#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle = CFLOBDDNodeHandle(g);
		iswap_hashMap.insert(std::make_pair(p, gHandle));
		return gHandle;
	}
	

	CFLOBDDNodeHandle MkCSwapGate2Node(unsigned int level, long int controller, long int index1, long int index2, int case_num)
	{
		std::string p = std::to_string(level) + ";" + std::to_string(controller) + ";" + std::to_string(index1) + ";" + std::to_string(index2) + ";" + std::to_string(case_num);
		if (cswap_hashMap.find(p) != cswap_hashMap.end()){
			return cswap_hashMap[p];
		}

		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(level);
		if (level == 1)
		{
			// std::cout << "CaseBase" << std::endl;
			if (controller == 0 && index1 == -1 && index2 == -1)
			{
				CFLOBDDReturnMapHandle m01, m12;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);
				g->numExits = 3;
			}
		}
		else 
		{
			unsigned int n = pow(2, level - 1);
			if (controller < n/2 && index1 < n/2 && index2 < n/2 && controller >= 0 && index1 >= 0 && index2 >= 0)
			{
				// Case 1: All fall in A
				// std::cout << "Case1" << std::endl;
				CFLOBDDNodeHandle aTmp = MkCSwapGate2Node(level-1, controller, index1, index2, case_num);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m01.Canonicalize();
				g->AConnection = Connection(aTmp, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller >= n/2 && index1 >= n/2 && index2 >= n/2 && controller >= 0 && index1 >= 0 && index2 >= 0)
			{
				// Case 2: All fall in B region
				// std::cout << "Case2" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDNodeHandle btmp = MkCSwapGate2Node(level-1, controller - n/2, index1 - n/2, index2 - n/2, case_num);
				g->BConnection[0] = Connection(btmp, m01);
				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller < n/2 && index1 >= n/2 && index2 >= n/2 && controller >= 0 && index1 >= 0 && index2 >= 0)
			{
				// Case 3: controller in A and index1 and index2 in B
				// std::cout << "Case3" << std::endl;
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				CFLOBDDNodeHandle atmp = MkCSwapGate2Node(level-1, controller, -1, -1, case_num);
				g->AConnection = Connection(atmp, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDNodeHandle bb = MkSwapGateNode(level-1, index1 - n/2, index2 - n/2, -1);
				CFLOBDDReturnMapHandle m01, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				auto Id = MkIdRelationInterleavedNode(level - 1);
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(bb, m01);
				g->numExits = 2;
			}
			else if (controller < n/2 && index1 == -1 && index2 == -1 && controller >= 0)
			{
				// Case 4: controller in A and index1 == index2 == -1
				// std::cout << "Case4" << std::endl;
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0);
				m012.AddToEnd(1);
				m012.AddToEnd(2);
				m012.Canonicalize();
				CFLOBDDNodeHandle atmp = MkCSwapGate2Node(level-1, controller, -1, -1, case_num);
				g->AConnection = Connection(atmp, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				CFLOBDDReturnMapHandle m01, m21, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else if (controller >= n/2 && index1 == -1 && index2 == -1 && controller >= 0)
			{
				// Case 5: controller in B and index1 == index2 == -1
				// std::cout << "Case5" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				CFLOBDDNodeHandle btmp = MkCSwapGate2Node(level-1, controller - n/2, -1, -1, case_num);
				CFLOBDDReturnMapHandle m012, m1;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(btmp, m012);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level - 1], m1);
				g->numExits = 3;
			} 
			else if (controller == -1 && index2 < n/2 && index2 >= 0 && index1 == -1)
			{
				// Case 6: controller == -1 && index2 in A && index1 == -1
				// std::cout << "Case6" << std::endl;
				CFLOBDDNodeHandle atmp = MkSwapGateNode(level-1, -1, index2, case_num);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(atmp, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m1, m0;
				m1.AddToEnd(1); m1.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				if (case_num == 0)
				{
					g->BConnection[0] = Connection(Id, m01);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else{
					CFLOBDDReturnMapHandle m10;
					m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
					g->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
					g->BConnection[1] = Connection(Id, m10);
				}
				g->numExits = 2;
			}
			else if (controller == -1 && index2 >= n/2 && index2 >= 0 && index1 == -1)
			{
				// Case 7: controller == -1 && index2 in B && index1 == -1
				// std::cout << "Case7" << std::endl;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level - 1);
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDReturnMapHandle m1, m0;
				m1.AddToEnd(1); m1.Canonicalize();
				m0.AddToEnd(0); m0.Canonicalize();
				CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, -1, index2 - n/2, case_num);
				if (case_num == 0)
				{
					g->BConnection[0] = Connection(btmp, m01);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else
				{
					g->BConnection[0] = Connection(btmp, m01);
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
				}
				g->numExits = 2;
			}
			else if (controller < n/2 && index1 < n/2 && index2 >= n/2 && controller >= 0 && index1 >= 0 && index2 >= 0)
			{
				// Case 8: CR in A, index1 in A, index2 in B
				CFLOBDDReturnMapHandle m012345;
				m012345.AddToEnd(0);
				m012345.AddToEnd(1);
				m012345.AddToEnd(2);
				m012345.AddToEnd(3);
				m012345.AddToEnd(4);
				m012345.AddToEnd(5);
				m012345.Canonicalize();

				auto aa = MkCSwapGate2Node(level - 1, controller, index1, -1, case_num);
				g->AConnection = Connection(aa, m012345);
				g->numBConnections = 6;
				g->BConnection = new Connection[6];
				auto Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01, m1, m10;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				auto b0 = MkSwapGateNode(level-1, -1, index2 - n/2, 0);
				auto b1 = MkSwapGateNode(level-1, -1, index2 - n/2, 2);
				auto b2 = MkSwapGateNode(level-1, -1, index2 - n/2, 1);
				auto b3 = MkSwapGateNode(level-1, -1, index2 - n/2, 3);
				g->BConnection[2] = Connection(b0, m01);
				g->BConnection[3] = Connection(b1, m10);
				g->BConnection[4] = Connection(b2, m10);
				g->BConnection[5] = Connection(b3, m10);
				g->numExits = 2;
			}
			else if (controller < n/2 && index1 < n/2 && index2 == -1 && controller >= 0 && index1 >= 0)
			{
				// Case 9: CR and index1 in A, index2 == -1
				CFLOBDDReturnMapHandle m012345;
				m012345.AddToEnd(0);
				m012345.AddToEnd(1);
				m012345.AddToEnd(2);
				m012345.AddToEnd(3);
				m012345.AddToEnd(4);
				m012345.AddToEnd(5);
				m012345.Canonicalize();

				auto aa = MkCSwapGate2Node(level - 1, controller, index1, -1, case_num);
				g->AConnection = Connection(aa, m012345);
				g->numBConnections = 6;
				g->BConnection = new Connection[6];
				auto Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01, m1, m21, m31, m41, m51;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				m31.AddToEnd(2); m31.AddToEnd(1); m31.Canonicalize();
				m41.AddToEnd(2); m41.AddToEnd(1); m41.Canonicalize();
				m51.AddToEnd(2); m51.AddToEnd(1); m51.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->BConnection[2] = Connection(Id, m21);
				g->BConnection[3] = Connection(Id, m31);
				g->BConnection[4] = Connection(Id, m41);
				g->BConnection[5] = Connection(Id, m51);
				g->numExits = 6;
			}
			else if (controller >= n/2 && index1 >= n/2 && index2 == -1 && controller >= 0 && index1 >= 0)
			{
				// Case 9: CR and index1 in B, index2 == -1
				auto Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				CFLOBDDReturnMapHandle m012345;
				m012345.AddToEnd(0);
				m012345.AddToEnd(1);
				m012345.AddToEnd(2);
				m012345.AddToEnd(3);
				m012345.AddToEnd(4);
				m012345.AddToEnd(5);
				m012345.Canonicalize();

				CFLOBDDReturnMapHandle m1;
				m1.AddToEnd(1); m1.Canonicalize();

				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				auto aa = MkCSwapGate2Node(level - 1, controller, index1, -1, case_num);
				g->BConnection[0] = Connection(aa, m012345);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 2;
			}
			else if (controller < n/2 && index1 >= n/2 && index2 == -1 && controller >= 0)
			{
				// Case 10: CR in A, index1 in B and index2 == -1
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				auto aa = MkCSwapGate2Node(level - 1, controller, -1, -1, case_num);
				g->AConnection = Connection(aa, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				auto Id = MkIdRelationInterleavedNode(level - 1);
				CFLOBDDReturnMapHandle m01, m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				auto bb = MkSwapGateNode(level-1, index1 - n/2, -1, case_num);
				CFLOBDDReturnMapHandle m21345;
				m21345.AddToEnd(2);
				m21345.AddToEnd(1);
				m21345.AddToEnd(3);
				m21345.AddToEnd(4);
				m21345.AddToEnd(5);
				m21345.Canonicalize();
				g->BConnection[2] = Connection(bb, m21345);
				g->numExits = 6;
			}
			
		}
		#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle = CFLOBDDNodeHandle(g);
		swap_hashMap.insert(std::make_pair(p, gHandle));
		return gHandle;
	}

	// i and j are always in the second half of the variables at top node
	CFLOBDDNodeHandle MkCSwapGateNode(unsigned int level, long int controller, long int i, long int j, int flag)
	{
		// std::cout << level << " " << controller << " " << i << " " << j << " " << flag << std::endl;
		// top level
		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(level);
		if (flag == 1)
		{
			assert(level >= 3);
			CFLOBDDNodeHandle atmp = MkCSwapGateNode(level-1, controller, i, j, 0);
			CFLOBDDReturnMapHandle m012;
			m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
			g->AConnection = Connection(atmp, m012);
			g->numBConnections = 3;
			g->BConnection = new Connection[3];
			CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			g->BConnection[0] = Connection(Id, m01);
			CFLOBDDReturnMapHandle m1;
			m1.AddToEnd(1); m1.Canonicalize();
			g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
			unsigned int n = pow(2, level - 1);
			CFLOBDDNodeHandle btmp = MkSwapGateNode(level-1, i-n/2,j-n/2, -1);
			g->BConnection[2] = Connection(btmp, m01);	
			g->numExits = 2;

		}
		else
		{
			unsigned int n = pow(2, level - 1);
			if (level == 1)
			{
				CFLOBDDReturnMapHandle m01, m12;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m12.AddToEnd(1); m12.AddToEnd(2); m12.Canonicalize();
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);
				g->numExits = 3;
			}
			else if (controller < n/2)
			{
				CFLOBDDNodeHandle atmp = MkCSwapGateNode(level-1,controller,i,j,flag);
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				g->AConnection = Connection(atmp, m012);
				g->numBConnections = 3;
				g->BConnection = new Connection[3];
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				CFLOBDDReturnMapHandle m01,m1,m21;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				m21.AddToEnd(2); m21.AddToEnd(1); m21.Canonicalize();
				g->BConnection[0] = Connection(Id, m01);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1],m1);
				g->BConnection[2] = Connection(Id, m21);
				g->numExits = 3;
			}
			else
			{
				CFLOBDDNodeHandle Id = MkIdRelationInterleavedNode(level-1);
				CFLOBDDReturnMapHandle m01,m1;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				g->AConnection = Connection(Id, m01);
				g->numBConnections = 2;
				g->BConnection = new Connection[2];
				CFLOBDDNodeHandle btmp = MkCSwapGateNode(level-1, controller - n/2, i,j,flag);
				CFLOBDDReturnMapHandle m012;
				m012.AddToEnd(0); m012.AddToEnd(1); m012.AddToEnd(2); m012.Canonicalize();
				m1.AddToEnd(1); m1.Canonicalize();
				g->BConnection[0] = Connection(btmp, m012);
				g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				g->numExits = 3;
			}
			
		}
#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle = CFLOBDDNodeHandle(g);
		return gHandle;
	}

	CFLOBDDNodeHandle ReverseColumnsNode(CFLOBDDNodeHandle nh)
	{
		CFLOBDDInternalNode *nhNode = (CFLOBDDInternalNode *)nh.handleContents;
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(nhNode->level);
		if (nhNode->level == 1)
		{
			n->AConnection = nhNode->AConnection;
			n->numBConnections = nhNode->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int i = 0; i < n->numBConnections; i++)
			{
				if (*(nhNode->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle){
					n->BConnection[i] == nhNode->BConnection[i];
				}
				else
				{
					CFLOBDDReturnMapHandle m;
					for (unsigned int j = nhNode->BConnection[i].returnMapHandle.Size() - 1; j >= 0; j--)  // TWR: Infinite loop -- unsigned int always >=0
					{
						m.AddToEnd(nhNode->BConnection[i].returnMapHandle[i]);   // TWR: should j be used?  E.g., nhNode->BConnection[i].returnMapHandle[j]
					}
					m.Canonicalize();
					n->BConnection[i] = Connection(*(nhNode->BConnection[i].entryPointHandle), m);
				}
			}
		}
		else
		{
			CFLOBDDNodeHandle reverseTmp = ReverseColumnsNode(*(nhNode->AConnection.entryPointHandle));
			n->AConnection = Connection(reverseTmp,nhNode->AConnection.returnMapHandle);
			n->numBConnections = nhNode->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int i = 0; i < n->numBConnections; i++)
			{
				CFLOBDDNodeHandle Btmp = ReverseColumnsNode(*(nhNode->BConnection[i].entryPointHandle));
				n->BConnection[i] = Connection(Btmp, nhNode->BConnection[i].returnMapHandle);
			}
		}
		n->numExits = nhNode->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkInverseReedMullerInterleavedNode(unsigned int i)
	{
		assert(i >= 1);

		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 1) {  // Base case
			CFLOBDDReturnMapHandle m1, m2;

			m1.AddToEnd(0);
			m1.AddToEnd(1);
			m1.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m1);

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			m2.AddToEnd(2);
			m2.AddToEnd(0);
			m2.Canonicalize();
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m1);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m2);
		}
		else {
			CFLOBDDReturnMapHandle m1, m2, m3;

			m1.AddToEnd(0);
			m1.AddToEnd(1);
			m1.AddToEnd(2);
			m1.Canonicalize();
			CFLOBDDNodeHandle temp = MkInverseReedMullerInterleavedNode(i - 1);
			n->AConnection = Connection(temp, m1);

			n->numBConnections = 3;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, m1);

			m2.AddToEnd(1);
			m2.Canonicalize();
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i - 1], m2);

			m3.AddToEnd(2);
			m3.AddToEnd(1);
			m3.AddToEnd(0);
			m3.Canonicalize();
			n->BConnection[2] = Connection(temp, m3);
		}
		n->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	} // MkInverseReedMullerInterleavedNode
	
	CFLOBDDNodeHandle MkWalshVoc13Node(unsigned int i)
	{
		assert(i >= 2);
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 2) {  // Base case
			CFLOBDDReturnMapHandle m0, m1, m01;

			m0.AddToEnd(0);
			m0.Canonicalize();
			m1.AddToEnd(1);
			m1.Canonicalize();
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();

			CFLOBDDInternalNode *nA = new CFLOBDDInternalNode(1);
			nA->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
			nA->numBConnections = 2;
			nA->BConnection = new Connection[nA->numBConnections];
			nA->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
			nA->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);
			nA->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
			nA->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nAHandle(nA);  // Create handle and canonicalize

			n->AConnection = Connection(nAHandle, m01);
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i - 1], m0);
			n->BConnection[1] = Connection(nAHandle, m01);
		}
		else {
			CFLOBDDReturnMapHandle m01, m10;

			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			CFLOBDDNodeHandle temp = MkWalshVoc13Node(i - 1);
			n->AConnection = Connection(temp, m01);

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, m01);
			m10.AddToEnd(1);
			m10.AddToEnd(0);
			m10.Canonicalize();
			n->BConnection[1] = Connection(temp, m10);
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	} // MkWalshVoc13Node


	CFLOBDDNodeHandle MkWalshVoc12Node(unsigned int i)
	{
		assert(i >= 2);
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 2) {  // Base case
			CFLOBDDReturnMapHandle m0, m1, m01;

			m0.AddToEnd(0);
			m0.Canonicalize();
			m1.AddToEnd(1);
			m1.Canonicalize();
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();

			CFLOBDDInternalNode *nA = new CFLOBDDInternalNode(1);
			nA->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
			nA->numBConnections = 2;
			nA->BConnection = new Connection[nA->numBConnections];
			nA->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
			nA->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
			nA->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
			nA->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nAHandle(nA);  // Create handle and canonicalize

			n->AConnection = Connection(nAHandle, m01);
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[1], m0);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[1], m1);
		}
		else {
			CFLOBDDReturnMapHandle m01, m10;

			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			CFLOBDDNodeHandle temp = MkWalshVoc12Node(i - 1);
			n->AConnection = Connection(temp, m01);

			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, m01);
			m10.AddToEnd(1);
			m10.AddToEnd(0);
			m10.Canonicalize();
			n->BConnection[1] = Connection(temp, m10);
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	} // MkWalshVoc12Node

	/*CFLOBDDNodeHandle MatrixConvertVocs1234To12Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh)) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}
		CFLOBDDInternalNode *nhNode = (CFLOBDDInternalNode *)nh.handleContents;
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(nhNode->level + 1);
		if (nhNode->level == 1)
		{
			CFLOBDDNodeHandle tempHandle(nh);
			CFLOBDDReturnMapHandle m1;
			for (int i = 0; i < tempHandle.handleContents->numExits; i++)
			{
				m1.AddToEnd(i);
			}
			m1.Canonicalize();
			n->AConnection = Connection(nh, m1);
			n->numBConnections = tempHandle.handleContents->numExits;
			n->BConnection = new Connection[n->numBConnections];
			for (int i = 0; i < n->numBConnections; i++)
			{
				CFLOBDDReturnMapHandle m2;
				m2.AddToEnd(i);
				m2.Canonicalize();
				n->BConnection[i] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[1], m2);
			}
		}
		else
		{
			CFLOBDDNodeHandle temp = MatrixConvertVocs1234To12Node(memoTable, nhNode->AConnection.entryPointHandle);
			n->AConnection = Connection(temp, nhNode->AConnection.returnMapHandle);
			n->numBConnections = nhNode->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			for (int i = 0; i < n->numBConnections; i++){
				n->BConnection[i] = Connection(MatrixConvertVocs1234To12Node(memoTable, nhNode->BConnection[i].entryPointHandle), nhNode->BConnection[i].returnMapHandle);
			}
		}
		n->numExits = nhNode->numExits;
		CFLOBDDNodeHandle nHandle(n);
		memoTable->Insert(nh, nHandle);
		return nHandle;
	}
	*/



	CFLOBDDNodeHandle MatrixConjugateNode(CFLOBDDNodeHandle nh)
	{
		CFLOBDDInternalNode *nhNode = (CFLOBDDInternalNode *)nh.handleContents;
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(nhNode->level);
		n->AConnection = nhNode->AConnection;
		n->numBConnections = nhNode->numBConnections;
		n->BConnection = new Connection[n->numBConnections];
		for (unsigned int i = 0; i < n->numBConnections; i++)
		{
			n->BConnection[i] == nhNode->BConnection[i];
		}
		n->numExits = nhNode->numExits;
		return CFLOBDDNodeHandle(n);
	}

	// Vocabulary shift in a matrix
	CFLOBDDNodeHandle MatrixShiftVocs13To24Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && (((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2)) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		assert(nh.handleContents->level >= 1);
		assert(nh.handleContents->NodeKind() == CFLOBDD_INTERNAL);

		CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)nh.handleContents;

		if (n->level == 1) {
			if (nh == CFLOBDDNodeHandle::NoDistinctionNode[1]) {
				return nh;
			}
			else {
				assert(*(n->AConnection.entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle);
				assert(*(n->BConnection[0].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle);
				assert(*(n->BConnection[1].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle);

				CFLOBDDInternalNode *g = new CFLOBDDInternalNode(1);
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
				g->numBConnections = 1;
				g->BConnection = new Connection[g->numBConnections];
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
				g->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
				g->InstallPathCounts();
#endif
				CFLOBDDNodeHandle gHandle(g);  // Create handle and canonicalize
				memoTable->Insert(nh, gHandle);  // Memoize result
				return gHandle;
			}
		}
		else { // Invoke the shift recursively on the AConnection and BConnection
			CFLOBDDInternalNode *g = new CFLOBDDInternalNode(n->level);
			CFLOBDDNodeHandle atmp = MatrixShiftVocs13To24Node(memoTable, *(n->AConnection.entryPointHandle));
			g->AConnection = Connection(atmp, n->AConnection.returnMapHandle);
			g->numBConnections = n->numBConnections;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int i = 0; i < g->numBConnections; i++) {
				CFLOBDDNodeHandle temp = MatrixShiftVocs13To24Node(memoTable, *(n->BConnection[i].entryPointHandle));
				g->BConnection[i] = Connection(temp, n->BConnection[i].returnMapHandle);
			}
			g->numExits = n->numExits;
#ifdef PATH_COUNTING_ENABLED
			g->InstallPathCounts();
#endif
			CFLOBDDNodeHandle gHandle(g);  // Create handle and canonicalize
			memoTable->Insert(nh, gHandle);  // Memoize result
			return gHandle;
		}
	} // MatrixShiftVocs13To24Node

	std::pair<CFLOBDDNodeHandle, CFLOBDDReturnMapHandle>
		MatrixTransposeNode(std::unordered_map<CFLOBDDNodeHandle, std::pair<CFLOBDDNodeHandle, CFLOBDDReturnMapHandle>, 
		CFLOBDDNodeHandle::CFLOBDDNodeHandle_Hash>& hashMap,
		CFLOBDDNodeHandle nh)
	{
		
		if (hashMap.find(nh) != hashMap.end() && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			return hashMap[nh];
		}
		
		CFLOBDDInternalNode *nhNode = (CFLOBDDInternalNode *)nh.handleContents;
		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(nhNode->level);
		CFLOBDDReturnMapHandle return_handle;
		if (nhNode->level == 1)
		{
			if (nh == CFLOBDDNodeHandle::NoDistinctionNode[1]){
				n = nhNode;
				return_handle.AddToEnd(0);
			}
			else if (*(nhNode->AConnection.entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle)
			{
				CFLOBDDReturnMapHandle m01, m0, m1;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m01.Canonicalize();

				m0.AddToEnd(0);
				m0.Canonicalize();

				m1.AddToEnd(1);
				m1.Canonicalize();

				n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				n->numBConnections = 2;
				n->BConnection = new Connection[2];
				assert(nhNode->BConnection[0].returnMapHandle.Lookup(0) == 0);
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);
				return_handle.AddToEnd(0);
				return_handle.AddToEnd(1);
			}
			else
			{
				bool isFork = false;
				for (unsigned int i = 0; i < nhNode->numBConnections; i++)
				{
					if (*(nhNode->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle)
						isFork = true;
				}
				if (isFork)
				{
					n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
					n->numBConnections = 2;
					CFLOBDDReturnMapHandle m0;
					if (nhNode->BConnection[0].returnMapHandle.Lookup(0) == nhNode->BConnection[1].returnMapHandle.Lookup(0)){
						m0.AddToEnd(0);
						return_handle.AddToEnd(nhNode->BConnection[0].returnMapHandle.Lookup(0));
					}
					else
					{
						m0.AddToEnd(0); m0.AddToEnd(1);
						return_handle.AddToEnd(nhNode->BConnection[0].returnMapHandle.Lookup(0));
						return_handle.AddToEnd(nhNode->BConnection[1].returnMapHandle.Lookup(0));
					}
					m0.Canonicalize();

					CFLOBDDReturnMapHandle m1;
					if (nhNode->BConnection[0].returnMapHandle.Lookup(nhNode->BConnection[0].returnMapHandle.Size() - 1)
						== nhNode->BConnection[1].returnMapHandle.Lookup(nhNode->BConnection[1].returnMapHandle.Size() - 1)){
						int v = nhNode->BConnection[0].returnMapHandle.Lookup(nhNode->BConnection[0].returnMapHandle.Size() - 1);
						if (return_handle.LookupInv(v) == -1){
							return_handle.AddToEnd(v);
							m1.AddToEnd(m0.Size());
						}
						else{
							m1.AddToEnd(v);
						}
					}
					else
					{
						int v1 = nhNode->BConnection[0].returnMapHandle.Lookup(nhNode->BConnection[0].returnMapHandle.Size() - 1);
						int v2 = nhNode->BConnection[1].returnMapHandle.Lookup(nhNode->BConnection[1].returnMapHandle.Size() - 1);
						if (return_handle.LookupInv(v1) == -1){
							return_handle.AddToEnd(v1);
							m1.AddToEnd(m0.Size());
						}
						else{
							m1.AddToEnd(return_handle.LookupInv(v1));
						}
						if (return_handle.LookupInv(v2) == -1){
							return_handle.AddToEnd(v2);
							m1.AddToEnd(v2);
						}
						else{
							m1.AddToEnd(return_handle.LookupInv(v2));
						}
					}
					m1.Canonicalize();

					n->BConnection = new Connection[2];
					if (m0.Size() == 1)
						n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
					else
						n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m0);

					if (m1.Size() == 1)
						n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);
					else
						n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m1);
				}
				else
				{
					n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
					n->numBConnections = 1;
					n->BConnection = new Connection[1];
					n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, nhNode->AConnection.returnMapHandle);
					return_handle.AddToEnd(nhNode->AConnection.returnMapHandle[0]);
					return_handle.AddToEnd(nhNode->AConnection.returnMapHandle[1]);
				}

			}
		}
		else
		{
			auto a_handle = MatrixTransposeNode(hashMap, *(nhNode->AConnection.entryPointHandle));
			n->AConnection = Connection(a_handle.first, nhNode->AConnection.returnMapHandle);
			n->numBConnections = nhNode->numBConnections;
			n->BConnection = new Connection[n->numBConnections];
			n->numExits = 0;
			std::unordered_map<unsigned int, unsigned int> reduction_map;
			for (unsigned int i = 0; i < n->numBConnections; i++)
			{
				CFLOBDDNodeHandle btmp = *(nhNode->BConnection[a_handle.second[i]].entryPointHandle);
				auto b_handle = MatrixTransposeNode(hashMap, btmp);
				CFLOBDDReturnMapHandle m;
				for (unsigned int j = 0; j < b_handle.second.Size(); j++){
					CFLOBDDReturnMapHandle b_return_map_handle = nhNode->BConnection[a_handle.second[i]].returnMapHandle;
					unsigned int val = b_return_map_handle[b_handle.second[j]];
					if (reduction_map.find(val) == reduction_map.end()){
						m.AddToEnd(n->numExits++);
						reduction_map.emplace(val, n->numExits - 1);
						return_handle.AddToEnd(val);
					}
					else{
						m.AddToEnd(reduction_map[val]);
					}
				}
				m.Canonicalize();
				n->BConnection[i] = Connection(b_handle.first, m);
			}
		}

		return_handle.Canonicalize();
		n->numExits = nhNode->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif

		CFLOBDDNodeHandle gHandle(n);  // Create handle and canonicalize
		hashMap.emplace(nh, std::make_pair(gHandle, return_handle));
		return std::make_pair(gHandle, return_handle);

	}

	// Vocabulary shift in a matrix
	// Requires that Vocs 3 and 4 are not in use
	CFLOBDDNodeHandle MatrixShiftVocs12To34Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		//std::cout << "node" << std::endl;
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			//std::cout << "node lookup" << std::endl;
			return lookupResult;
		}

		assert(nh.handleContents->level >= 2);
		assert(nh.handleContents->NodeKind() == CFLOBDD_INTERNAL);

		CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)nh.handleContents;

		if (n->level == 2) {
			if (nh == CFLOBDDNodeHandle::NoDistinctionNode[2]) {
				return nh;
			}
			else {  // Shift the AConnection to BConnection[0]
				for (unsigned int i = 0; i < n->numBConnections; i++) {  // Make sure that Vocs 3 and 4 are not in use
					assert(*(n->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[1]);
				}


				CFLOBDDInternalNode *g = new CFLOBDDInternalNode(2);
				g->AConnection = Connection(CFLOBDDNodeHandle::NoDistinctionNode[1], commonly_used_return_maps[0]);//m0
				g->numBConnections = 1;
				g->BConnection = new Connection[g->numBConnections];
				g->BConnection[0] = n->AConnection;
				g->numExits = n->numBConnections;
#ifdef PATH_COUNTING_ENABLED
				g->InstallPathCounts();
#endif
				CFLOBDDNodeHandle gHandle(g);  // Create handle and canonicalize
				memoTable->Insert(nh, gHandle);  // Memoize result
				return gHandle;
			}
		}
		else { // Invoke the shift recursively on the AConnection and BConnection
			CFLOBDDInternalNode *g = new CFLOBDDInternalNode(n->level);
			CFLOBDDNodeHandle atmp = MatrixShiftVocs12To34Node(memoTable, *(n->AConnection.entryPointHandle));
			g->AConnection = Connection(atmp, n->AConnection.returnMapHandle);
			g->numBConnections = n->numBConnections;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int i = 0; i < g->numBConnections; i++) {
				CFLOBDDNodeHandle temp = MatrixShiftVocs12To34Node(memoTable, *(n->BConnection[i].entryPointHandle));
				g->BConnection[i] = Connection(temp, n->BConnection[i].returnMapHandle);
			}
			g->numExits = n->numExits;
#ifdef PATH_COUNTING_ENABLED
			g->InstallPathCounts();
#endif
			CFLOBDDNodeHandle gHandle(g);  // Create handle and canonicalize
			memoTable->Insert(nh, gHandle);  // Memoize result
			return gHandle;
		}
	} // MatrixShiftVocs12To34Node

	// Vocabulary shift in a matrix
	// Assumes that Voc 3 is "don't care" everywhere
	CFLOBDDNodeHandle MatrixShiftVoc43Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		assert(nh.handleContents->level >= 2);
		assert(nh.handleContents->NodeKind() == CFLOBDD_INTERNAL);

		CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)nh.handleContents;

		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(n->level);
		if (n->level == 2) {
			g->AConnection = Connection(*(n->AConnection.entryPointHandle), n->AConnection.returnMapHandle); // Vocabularies 1 and 2 preserved
			g->numBConnections = n->numBConnections;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int i = 0; i < g->numBConnections; i++) {
				if (*(n->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[1]) {  // Nothing to do
					g->BConnection[i] = n->BConnection[i];
				}
				else { // Shift n->BConnection[i]'s BConnections to g->BConnection[i]'s AConnection
					assert(n->BConnection[i].entryPointHandle->handleContents->NodeKind() == CFLOBDD_INTERNAL);
					CFLOBDDInternalNode *nb_i = (CFLOBDDInternalNode *)n->BConnection[i].entryPointHandle->handleContents;  // A node at level 1
					assert(*(nb_i->AConnection.entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle);
					assert(*(nb_i->BConnection[0].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle);


					CFLOBDDInternalNode *gb = new CFLOBDDInternalNode(1);
					gb->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
					gb->numBConnections = 2;
					gb->BConnection = new Connection[gb->numBConnections];
					gb->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
					gb->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[1]);//m1
					gb->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
					gb->InstallPathCounts();
#endif
					CFLOBDDNodeHandle gbHandle(gb);  // Create handle and canonicalize

					g->BConnection[i] = Connection(gbHandle, n->BConnection[i].returnMapHandle);
				}
			}
		}
		else { // Invoke the shift recursively on the AConnection and BConnection
			CFLOBDDNodeHandle atmp = MatrixShiftVoc43Node(memoTable, *(n->AConnection.entryPointHandle));
			g->AConnection = Connection(atmp, n->AConnection.returnMapHandle);
			g->numBConnections = n->numBConnections;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int i = 0; i < g->numBConnections; i++) {
				CFLOBDDNodeHandle temp = MatrixShiftVoc43Node(memoTable, *(n->BConnection[i].entryPointHandle));
				g->BConnection[i] = Connection(temp, n->BConnection[i].returnMapHandle);
			}
		}
		g->numExits = n->numExits;
#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif

		CFLOBDDNodeHandle gHandle(g);  // Create handle and canonicalize
		memoTable->Insert(nh, gHandle);  // Memoize result
		return gHandle;
	} // MatrixShiftVoc43Node

	// Vocabulary shift in a matrix
	// Requires that Vocs 2 and 3 are not in use (i.e., have been projected out)
	CFLOBDDNodeHandle MatrixShiftVoc42Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		assert(nh.handleContents->level >= 2);
		assert(nh.handleContents->NodeKind() == CFLOBDD_INTERNAL);

		CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)nh.handleContents;

		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(n->level);
		if (n->level == 2) {
			CFLOBDDInternalNode *gA = new CFLOBDDInternalNode(1);
			CFLOBDDInternalNode *nA = (CFLOBDDInternalNode *)n->AConnection.entryPointHandle->handleContents;
			gA->AConnection = nA->AConnection;
			gA->numBConnections = nA->numBConnections;
			gA->BConnection = new Connection[gA->numBConnections];
			for (unsigned int i = 0; i < nA->numBConnections; i++) {  // Make sure that Vocs 2 and 3 not in use
				assert(*(nA->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle);  // Check Voc 2
				CFLOBDDInternalNode *nB = (CFLOBDDInternalNode *)n->BConnection[i].entryPointHandle->handleContents;
				assert(*(nB->AConnection.entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle);  // Check Voc 3
				CFLOBDDNodeHandle temp = *(nB->BConnection[0].entryPointHandle);
				gA->BConnection[i] = Connection(temp, n->BConnection[i].returnMapHandle);  // Node n holds the appropriate return map
			}
			gA->numExits = n->numExits;
#ifdef PATH_COUNTING_ENABLED
			gA->InstallPathCounts();
#endif
			CFLOBDDNodeHandle gAHandle(gA);  // Create handle and canonicalize
			CFLOBDDReturnMapHandle mId;
			for (unsigned int j = 0; j < n->numExits; j++) {
				mId.AddToEnd(j);
			}
			mId.Canonicalize();

			g->AConnection = Connection(gAHandle, mId);
			g->numBConnections = n->numExits;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int k = 0; k < g->numBConnections; k++) {
				CFLOBDDReturnMapHandle m;
				m.AddToEnd(k);
				m.Canonicalize();
				g->BConnection[k] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[1], m);
			}
		}
		else { // Invoke the shift recursively on the AConnection and BConnection
			CFLOBDDNodeHandle atmp = MatrixShiftVoc42Node(memoTable, *(n->AConnection.entryPointHandle));
			g->AConnection = Connection(atmp, n->AConnection.returnMapHandle);
			g->numBConnections = n->numBConnections;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int i = 0; i < g->numBConnections; i++) {
				CFLOBDDNodeHandle temp = MatrixShiftVoc42Node(memoTable, *(n->BConnection[i].entryPointHandle));
				g->BConnection[i] = Connection(temp, n->BConnection[i].returnMapHandle);
			}
		}
		g->numExits = n->numExits;
#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif

		CFLOBDDNodeHandle gHandle(g);  // Create handle and canonicalize
		memoTable->Insert(nh, gHandle);  // Memoize result
		return gHandle;
	} // MatrixShiftVoc42Node

	// Create representation of a matrix in which vocabularies 2 and 3 are constrained to be equal:
	// (W,X,Y,Z) s.t. X==Y with interleaved variables
	CFLOBDDNodeHandle MkDetensorConstraintInterleavedNode(unsigned int i)
	{
		assert(i >= 2);

		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 2) {  // Base case

			CFLOBDDInternalNode *nA = new CFLOBDDInternalNode(1);
			CFLOBDDReturnMapHandle m0, m01;
			m0.AddToEnd(0);
			m0.Canonicalize();
			nA->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);//m0
			nA->numBConnections = 1;
			nA->BConnection = new Connection[nA->numBConnections];
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			nA->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			nA->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
			nA->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nAHandle(nA);  // Create handle and canonicalize

			CFLOBDDInternalNode *nB = new CFLOBDDInternalNode(1);
			CFLOBDDReturnMapHandle m1;
			nB->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);//m01
			nB->numBConnections = 2;
			nB->BConnection = new Connection[nB->numBConnections];
			nB->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);//m0
			m1.AddToEnd(1);
			m1.Canonicalize();
			nB->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);//m1
			nB->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
			nB->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nBHandle(nB);  // Create handle and canonicalize

			n->AConnection = Connection(nAHandle, m01);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(nBHandle, m01);//m01
			CFLOBDDReturnMapHandle m10;
			m10.AddToEnd(1);
			m10.AddToEnd(0);
			m10.Canonicalize();
			n->BConnection[1] = Connection(nBHandle, m10);//m10
		}
		else {

			CFLOBDDNodeHandle temp = MkDetensorConstraintInterleavedNode(i - 1);
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			n->AConnection = Connection(temp, m01);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, m01);//m01
			CFLOBDDReturnMapHandle m1;
			m1.AddToEnd(1);
			m1.Canonicalize();
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i - 1], m1);//m1
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	} // MkDetensorConstraintInterleavedNode

	// Create representation of a matrix in which vocabularies 1 and 4 are constrained to be equal:
	// (W,X,Y,Z) s.t. W==Z with interleaved variables
	CFLOBDDNodeHandle MkCFLOBDDMatrixEqVoc14Node(unsigned int i)
	{
		assert(i >= 2);

		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 2) {  // Base case

			CFLOBDDInternalNode *nA = new CFLOBDDInternalNode(1);
			nA->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
			nA->numBConnections = 2;
			nA->BConnection = new Connection[nA->numBConnections];
			nA->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
			nA->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[1]);//m1
			nA->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
			nA->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nAHandle(nA);  // Create handle and canonicalize

			CFLOBDDInternalNode *nB = new CFLOBDDInternalNode(1);
			nB->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
			nB->numBConnections = 1;
			nB->BConnection = new Connection[nB->numBConnections];
			nB->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
			nB->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
			nB->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nBHandle(nB);  // Create handle and canonicalize

			n->AConnection = Connection(nAHandle, commonly_used_return_maps[2]);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(nBHandle, commonly_used_return_maps[2]);//m01
			n->BConnection[1] = Connection(nBHandle, commonly_used_return_maps[3]);//m10
		}
		else {
			
			CFLOBDDNodeHandle temp = MkCFLOBDDMatrixEqVoc14Node(i - 1);
			n->AConnection = Connection(temp, commonly_used_return_maps[2]);//m01
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i - 1], commonly_used_return_maps[1]);//m1
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	} // MkCFLOBDDMatrixEqVoc14Node

	CFLOBDDNodeHandle MkFourierDiagonalComponentNode(unsigned int i)
	{
		assert(i >= 2);

		CFLOBDDInternalNode *n = new CFLOBDDInternalNode(i);
		if (i == 2) {  // Base case
			CFLOBDDReturnMapHandle m12, m012, m314;
			m12.AddToEnd(1);
			m12.AddToEnd(2);
			m12.Canonicalize();
			m012.AddToEnd(0);
			m012.AddToEnd(1);
			m012.AddToEnd(2);
			m012.Canonicalize();
			m314.AddToEnd(3);
			m314.AddToEnd(1);
			m314.AddToEnd(4);
			m314.Canonicalize();

			CFLOBDDInternalNode *nA = new CFLOBDDInternalNode(1);
			nA->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
			nA->numBConnections = 2;
			nA->BConnection = new Connection[nA->numBConnections];
			nA->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
			nA->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m12);
			nA->numExits = 3;
#ifdef PATH_COUNTING_ENABLED
			nA->InstallPathCounts();
#endif
			CFLOBDDNodeHandle nAHandle(nA);  // Create handle and canonicalize

			n->AConnection = Connection(nAHandle, m012);
			n->numBConnections = 3;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(nAHandle, m012);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i - 1], commonly_used_return_maps[1]);//m1
			n->BConnection[2] = Connection(nAHandle, m314);
			n->numExits = 5;
		}
		else {

			unsigned long int M = 1UL + (1UL << (1 << (i - 2)));
			CFLOBDDReturnMapHandle mA;
			for (unsigned long int j = 0; j < M; j++) {
				mA.AddToEnd(j);
			}
			mA.Canonicalize();

			CFLOBDDNodeHandle temp = MkFourierDiagonalComponentNode(i - 1);
			n->AConnection = Connection(temp, mA);
			n->numBConnections = M;
			n->BConnection = new Connection[n->numBConnections];
			n->BConnection[0] = Connection(temp, mA);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[i - 1], commonly_used_return_maps[1]);//m1
			unsigned long int k = M;
			for (unsigned long int j = 2; j < M; j++) {
				CFLOBDDReturnMapHandle mB;
				mB.AddToEnd(k);
				k++;
				mB.AddToEnd(1);
				for (unsigned long int p = 2; p < M; p++) {
					mB.AddToEnd(k);
					k++;
				}
				mB.Canonicalize();
				n->BConnection[j] = Connection(temp, mB);
			}
			n->numExits = 1UL + (1UL << (1 << (i - 1)));
		}
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	} // MkFourierDiagonalComponentNode

	CFLOBDDNodeHandle PromoteInterleavedTo12Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		assert(nh.handleContents->level >= 1);
		assert(nh.handleContents->NodeKind() == CFLOBDD_INTERNAL);

		CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)nh.handleContents;
		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(n->level + 1);

		if (n->level == 1) {
			CFLOBDDReturnMapHandle mId;
			for (unsigned int k = 0; k < n->numExits; k++) {
				mId.AddToEnd(k);
			}
			mId.Canonicalize();
			g->AConnection = Connection(nh, mId);
			g->numBConnections = n->numExits;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int k = 0; k < g->numBConnections; k++) {
				CFLOBDDReturnMapHandle mk;
				mk.AddToEnd(k);
				mk.Canonicalize();
				g->BConnection[k] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[1],mk);
			}
		}
		else { // Invoke PromoteInterleavedTo12Node recursively on the AConnection and BConnection
			CFLOBDDNodeHandle atmp = PromoteInterleavedTo12Node(memoTable, *(n->AConnection.entryPointHandle));
			g->AConnection = Connection(atmp, n->AConnection.returnMapHandle);
			g->numBConnections = n->numBConnections;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int i = 0; i < g->numBConnections; i++) {
				CFLOBDDNodeHandle temp = PromoteInterleavedTo12Node(memoTable, *(n->BConnection[i].entryPointHandle));
				g->BConnection[i] = Connection(temp, n->BConnection[i].returnMapHandle);
			}
		}
		g->numExits = n->numExits;
#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle(g);  // Create handle and canonicalize
		memoTable->Insert(nh, gHandle);  // Memoize result
		return gHandle;
	} // PromoteInterleavedTo12Node	
	
	CFLOBDDNodeHandle Demote12ToInterleavedNode(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		assert(nh.handleContents->level >= 2);
		assert(nh.handleContents->NodeKind() == CFLOBDD_INTERNAL);

		CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)nh.handleContents;

		if (n->level == 2) {
			if (nh == CFLOBDDNodeHandle::NoDistinctionNode[2]) {
				return CFLOBDDNodeHandle::NoDistinctionNode[1];
			}
			else {
				// Make sure that Vocs 3 and 4 are not in use
				for (unsigned int i = 0; i < n->numBConnections; i++) {
					assert(*(n->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::NoDistinctionNode[1]);
				}
				// Return the AConnection
				return *(n->AConnection.entryPointHandle);
			}
		}
		else { // Invoke Demote12ToInterleavedNode recursively on the AConnection and BConnection
			CFLOBDDInternalNode *g = new CFLOBDDInternalNode(n->level - 1);
			CFLOBDDNodeHandle atmp = Demote12ToInterleavedNode(memoTable, *(n->AConnection.entryPointHandle));
			g->AConnection = Connection(atmp, n->AConnection.returnMapHandle);
			g->numBConnections = n->numBConnections;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int i = 0; i < g->numBConnections; i++) {
				CFLOBDDNodeHandle temp = Demote12ToInterleavedNode(memoTable, *(n->BConnection[i].entryPointHandle));
				g->BConnection[i] = Connection(temp, n->BConnection[i].returnMapHandle);
			}
			g->numExits = n->numExits;
#ifdef PATH_COUNTING_ENABLED
			g->InstallPathCounts();
#endif
			CFLOBDDNodeHandle gHandle(g);    // Create handle and canonicalize
			memoTable->Insert(nh, gHandle);  // Memoize result
			return gHandle;
		}
	} // Demote12ToInterleavedNode

	CFLOBDDNodeHandle PromoteInterleavedTo13Node(CFLOBDDNodeMemoTableRefPtr memoTable, CFLOBDDNodeHandle nh)
	{
		if (memoTable->Lookup(nh) && ((((CFLOBDDInternalNode *)(nh.handleContents))->GetRefCount() >= 2))) {  // Check if already processed
			CFLOBDDNodeHandle lookupResult;
			memoTable->Fetch(nh, lookupResult);
			return lookupResult;
		}

		CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)nh.handleContents;
		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(n->level + 1);

		if (n->level == 0) {
			CFLOBDDReturnMapHandle mId;
			for (unsigned int k = 0; k < n->numExits; k++) {
				mId.AddToEnd(k);
			}
			mId.Canonicalize();
			g->AConnection = Connection(nh, mId);
			g->numBConnections = n->numExits;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int k = 0; k < g->numBConnections; k++) {
				CFLOBDDReturnMapHandle mk;
				mk.AddToEnd(k);
				mk.Canonicalize();
				g->BConnection[k] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[0], mk);
			}
		}
		else { // Invoke PromoteInterleavedTo13Node recursively on the AConnection and BConnection
			CFLOBDDNodeHandle atmp = PromoteInterleavedTo13Node(memoTable, *(n->AConnection.entryPointHandle));
			g->AConnection = Connection(atmp, n->AConnection.returnMapHandle);
			g->numBConnections = n->numBConnections;
			g->BConnection = new Connection[g->numBConnections];
			for (unsigned int i = 0; i < g->numBConnections; i++) {
				CFLOBDDNodeHandle temp = PromoteInterleavedTo13Node(memoTable, *(n->BConnection[i].entryPointHandle));
				g->BConnection[i] = Connection(temp, n->BConnection[i].returnMapHandle);
			}
		}
		g->numExits = n->numExits;
#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif
		CFLOBDDNodeHandle gHandle(g);  // Create handle and canonicalize
		memoTable->Insert(nh, gHandle);  // Memoize result
		return gHandle;
	} // PromoteInterleavedTo13Node	

	CFLOBDDNodeHandle SMatrixNode(std::string s)
	{
		CFLOBDDInternalNode* n = new CFLOBDDInternalNode(log2(s.length()) + 1);

		if (s.find('1') == std::string::npos)
		{
			return CFLOBDDNodeHandle::NoDistinctionNode[n->level];
		}
		else if (s.length() == 1)
		{
			if (s[0] == '1')
			{
				n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
				n->numBConnections = 2;
				n->BConnection = new Connection[n->numBConnections];
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
				n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[1]);//m1
				n->numExits = 2;
			}
			else
			{
				return CFLOBDDNodeHandle::NoDistinctionNode[n->level];
			}
		}
		else
		{
			CFLOBDDNodeHandle tempHandle = SMatrixNode(s.substr(0, s.length() / 2));
			CFLOBDDReturnMapHandle m;
			for (unsigned int i = 0; i < tempHandle.handleContents->numExits; i++)
			{
				m.AddToEnd(i);
			}
			m.Canonicalize();
			n->AConnection = Connection(tempHandle, m);
			n->numBConnections = tempHandle.handleContents->numExits;
			n->BConnection = new Connection[n->numBConnections];
			CFLOBDDNodeHandle temp = SMatrixNode(s.substr(s.length() / 2));
			assert(n->numBConnections == 1 || n->numBConnections == 2);
			if (n->numBConnections == 1)
			{
				assert(temp.handleContents->numExits == 2);
				n->BConnection[0] = Connection(temp, commonly_used_return_maps[2]);//m01
			}
			else
			{
				assert(temp.handleContents->numExits == 1);
				n->BConnection[0] = Connection(temp, commonly_used_return_maps[0]);//m0
				n->BConnection[1] = Connection(temp, commonly_used_return_maps[1]);//m1
			}
			n->numExits = (s.find('1') == std::string::npos) ? 1 : 2;
		}
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		CFLOBDDNodeHandle nHandle(n);  // Create handle and canonicalize
		return nHandle;
	}


	CFLOBDDNodeHandle MatrixShiftToAConnectionNode(CFLOBDDNodeHandle c)
	{
		CFLOBDDInternalNode* n = new CFLOBDDInternalNode(c.handleContents->level + 1);
		CFLOBDDReturnMapHandle m;
		for (unsigned int i = 0; i < c.handleContents->numExits; i++)
			m.AddToEnd(i);
		m.Canonicalize();
		n->AConnection = Connection(c, m);
		n->numBConnections = c.handleContents->numExits;
		n->BConnection = new Connection[n->numBConnections];
		for (unsigned int i = 0; i < n->numBConnections; i++)
		{
			CFLOBDDReturnMapHandle mi;
			mi.AddToEnd(i);
			mi.Canonicalize();
			n->BConnection[i] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[c.handleContents->level], mi);
		}
		n->numExits = n->numBConnections;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MatrixShiftToBConnectionNode(CFLOBDDNodeHandle c)
	{
		CFLOBDDInternalNode* n = new CFLOBDDInternalNode(c.handleContents->level + 1);
		n->AConnection = Connection(CFLOBDDNodeHandle::NoDistinctionNode[c.handleContents->level], commonly_used_return_maps[0]);//m0
		n->numBConnections = 1;
		n->BConnection = new Connection[n->numBConnections];
		
		CFLOBDDReturnMapHandle m;
		for (unsigned int i = 0; i < c.handleContents->numExits; i++)
			m.AddToEnd(i);
		m.Canonicalize();

		n->BConnection[0] = Connection(c, m);
		n->numExits = c.handleContents->numExits;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}



	CFLOBDDNodeHandle MkDistinctionTwoVarsNode(int x, int y, unsigned int var_level, unsigned int matrix_level)
	{
		CFLOBDDInternalNode* n = new CFLOBDDInternalNode(matrix_level);

		if (matrix_level == 1)
		{
			assert(var_level == 0);
			if (x == 0 && y == 0)
			{ // 1, 0
				CFLOBDDReturnMapHandle m01, m1;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m1.AddToEnd(1);
				m01.Canonicalize();
				m1.Canonicalize();
				n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				n->numBConnections = 2;
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);
				n->numExits = 2;
			}
			else if (x == 0 && y == 1)
			{ // 0,1
				CFLOBDDReturnMapHandle m01, m0;
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m0.AddToEnd(0);
				m01.Canonicalize();
				m0.Canonicalize();
				n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				n->numBConnections = 2;
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				n->numExits = 2;
			}
			else if (x == 1 && y == 0)
			{ // 0, 1
				CFLOBDDReturnMapHandle m10, m0, m01;
				m10.AddToEnd(1);
				m10.AddToEnd(0);
				m0.AddToEnd(0);
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m10.Canonicalize();
				m0.Canonicalize();
				m01.Canonicalize();
				n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				n->numBConnections = 2;
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m10);
				n->numExits = 2;
			}
			else
			{ // 0,1
				CFLOBDDReturnMapHandle m0, m01;
				m0.AddToEnd(0);
				m01.AddToEnd(0);
				m01.AddToEnd(1);
				m0.Canonicalize();
				m01.Canonicalize();
				n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				n->numBConnections = 2;
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
				n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
				n->numExits = 2;
			}
		}
		else if (var_level < (unsigned int)(1 << (matrix_level - 1)))
		{
			CFLOBDDNodeHandle tmp = MkDistinctionTwoVarsNode(x, y, var_level, matrix_level - 1);
			n->numBConnections = tmp.handleContents->numExits;
			CFLOBDDReturnMapHandle mi;
			for (int i = 0; i < n->numBConnections; i++)
				mi.AddToEnd(i);
			mi.Canonicalize();
			n->AConnection = Connection(tmp, mi);
			n->BConnection = new Connection[n->numBConnections];
			for (int i = 0; i < n->numBConnections; i++)
			{
				CFLOBDDReturnMapHandle m;
				m.AddToEnd(i);
				m.Canonicalize();
				n->BConnection[i] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[matrix_level - 1], m);
			}
			n->numExits = n->numBConnections;
		}
		else
		{
			CFLOBDDNodeHandle tmp = MkDistinctionTwoVarsNode(x, y, var_level / ((unsigned int)(1 << (matrix_level - 1))), matrix_level - 1);
			n->numBConnections = 1;
			CFLOBDDReturnMapHandle m0;
			m0.AddToEnd(0);
			m0.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::NoDistinctionNode[matrix_level-1], m0);
			n->BConnection = new Connection[n->numBConnections];
			CFLOBDDReturnMapHandle mi;
			for (int i = 0; i < tmp.handleContents->numExits; i++)
				mi.AddToEnd(i);
			mi.Canonicalize();
			n->BConnection[0] = Connection(tmp, mi);
			n->numExits = tmp.handleContents->numExits;
		}

#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	std::pair<CFLOBDDNodeHandle, std::vector<std::vector<std::pair<int, int>>>> MultiplyOperationNode(CFLOBDDNodeHandle c)
	{
		CFLOBDDReturnMapHandle m;
		for (unsigned int i = 0; i < c.handleContents->numExits; i++)
			m.AddToEnd(i);
		m.Canonicalize();
		return MultiplyOperationNodeInternal(c, 'A', c.handleContents->level, m);
	}

	std::pair<CFLOBDDNodeHandle, std::vector<std::vector<std::pair<int, int>>>> MultiplyOperationNodeInternal(CFLOBDDNodeHandle c, char position, unsigned int maxLevel, CFLOBDDReturnMapHandle cReturnMapHandle)
	{
		CFLOBDDInternalNode *n = (CFLOBDDInternalNode *)c.handleContents;
		CFLOBDDInternalNode *g = new CFLOBDDInternalNode(n->level - 1);
		if (n->level == 1)
		{
				std::vector<std::vector<std::pair<int, int>>> returnGeneralMap;
				if (position == 'A'){
					if (!(n->numBConnections == 2 && *(n->BConnection[0].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle &&
					*(n->BConnection[1].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle && n->BConnection[0].returnMapHandle.Lookup(0) == n->BConnection[1].returnMapHandle.Lookup(1) && 
					cReturnMapHandle[n->BConnection[0].returnMapHandle[1]] == cReturnMapHandle[n->BConnection[1].returnMapHandle[0]])){
						for (int i = 0; i < n->numBConnections; i++)
						{
							std::vector<std::pair<int, int>> returnMapI;
							if (*(n->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle)
							{
								returnMapI.push_back(std::make_pair(2, cReturnMapHandle[n->BConnection[i].returnMapHandle[0]]));
							}
							else
							{
								returnMapI.push_back(std::make_pair(1, std::min(cReturnMapHandle[n->BConnection[i].returnMapHandle[0]], cReturnMapHandle[n->BConnection[i].returnMapHandle[1]])));
								returnMapI.push_back(std::make_pair(1, std::max(cReturnMapHandle[n->BConnection[i].returnMapHandle[1]], cReturnMapHandle[n->BConnection[i].returnMapHandle[0]])));
							}
							returnGeneralMap.push_back(returnMapI);
						}
						return std::make_pair(*(n->AConnection.entryPointHandle), returnGeneralMap);
					}
					else{
						std::vector<std::pair<int, int>> returnMap;
						returnMap.push_back(std::make_pair(1, std::min(cReturnMapHandle[n->BConnection[0].returnMapHandle[0]], cReturnMapHandle[n->BConnection[0].returnMapHandle[1]])));
						returnMap.push_back(std::make_pair(1, std::max(cReturnMapHandle[n->BConnection[0].returnMapHandle[1]], cReturnMapHandle[n->BConnection[0].returnMapHandle[0]])));
						returnGeneralMap.push_back(returnMap);
						return std::make_pair(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, returnGeneralMap);
					}
				}
				else{
					if (*(n->AConnection.entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle)
					{
						std::vector<std::pair<int, int>> returnMapI;
						if (*(n->BConnection[0].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle && 
						*(n->BConnection[1].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle)
						{
							returnMapI.push_back(std::make_pair(1, std::min(cReturnMapHandle[n->BConnection[0].returnMapHandle[0]], cReturnMapHandle[n->BConnection[1].returnMapHandle[0]])));
							returnMapI.push_back(std::make_pair(1, std::max(cReturnMapHandle[n->BConnection[1].returnMapHandle[0]], cReturnMapHandle[n->BConnection[0].returnMapHandle[0]])));
							returnGeneralMap.push_back(returnMapI);
							return std::make_pair(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, returnGeneralMap);
						}
						else if (*(n->BConnection[0].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle && 
						*(n->BConnection[1].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle)
						{
							for (int i = 0; i < 2; i++){
								if (cReturnMapHandle[n->BConnection[0].returnMapHandle[0]] == cReturnMapHandle[n->BConnection[1].returnMapHandle[i]]){
									returnMapI.push_back(std::make_pair(2, cReturnMapHandle[n->BConnection[0].returnMapHandle[0]]));
								}
								else {
									returnMapI.push_back(std::make_pair(1, std::min(cReturnMapHandle[n->BConnection[0].returnMapHandle[0]], cReturnMapHandle[n->BConnection[1].returnMapHandle[i]])));
									returnMapI.push_back(std::make_pair(1, std::max(cReturnMapHandle[n->BConnection[0].returnMapHandle[0]], cReturnMapHandle[n->BConnection[1].returnMapHandle[i]])));
								}
								returnGeneralMap.push_back(returnMapI);
								returnMapI.clear();
							}
							return std::make_pair(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, returnGeneralMap);
						}
						else if (*(n->BConnection[1].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle && 
						*(n->BConnection[0].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle)
						{
							for (int i = 0; i < 2; i++){
								if (cReturnMapHandle[n->BConnection[1].returnMapHandle[0]] != cReturnMapHandle[n->BConnection[0].returnMapHandle[i]]){
									returnMapI.push_back(std::make_pair(1, std::min(cReturnMapHandle[n->BConnection[1].returnMapHandle[0]], cReturnMapHandle[n->BConnection[0].returnMapHandle[i]])));
									returnMapI.push_back(std::make_pair(1, std::max(cReturnMapHandle[n->BConnection[1].returnMapHandle[0]], cReturnMapHandle[n->BConnection[0].returnMapHandle[i]])));
								}
								else{
									returnMapI.push_back(std::make_pair(2, cReturnMapHandle[n->BConnection[0].returnMapHandle[i]]));
								}
								returnGeneralMap.push_back(returnMapI);
								returnMapI.clear();
							}
							return std::make_pair(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, returnGeneralMap);
						}
						else if (*(n->BConnection[0].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle && 
						*(n->BConnection[1].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDForkNodeHandle)
						{
							for (int i = 0; i < 2; i++){
								if (cReturnMapHandle[n->BConnection[0].returnMapHandle[i]] != cReturnMapHandle[n->BConnection[1].returnMapHandle[i]]){
									returnMapI.push_back(std::make_pair(1, std::min(cReturnMapHandle[n->BConnection[0].returnMapHandle[i]], cReturnMapHandle[n->BConnection[1].returnMapHandle[i]])));
									returnMapI.push_back(std::make_pair(1, std::max(cReturnMapHandle[n->BConnection[0].returnMapHandle[i]], cReturnMapHandle[n->BConnection[1].returnMapHandle[i]])));
								}
								else{
									returnMapI.push_back(std::make_pair(2, cReturnMapHandle[n->BConnection[1].returnMapHandle[i]]));
								}
								returnGeneralMap.push_back(returnMapI);
								returnMapI.clear();
							}
							return std::make_pair(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, returnGeneralMap);
						}
					}
					else{
						std::vector<std::pair<int, int>> returnMapI;
						if (*(n->BConnection[0].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle)
						{
							returnMapI.push_back(std::make_pair(2, cReturnMapHandle[n->BConnection[0].returnMapHandle[0]]));
							returnGeneralMap.push_back(returnMapI);
							return std::make_pair(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, returnGeneralMap);
						}
						else{
							returnMapI.push_back(std::make_pair(2, cReturnMapHandle[n->BConnection[0].returnMapHandle[0]]));
							returnGeneralMap.push_back(returnMapI);
							std::vector<std::pair<int, int>> tmpReturnMapI;
							tmpReturnMapI.push_back(std::make_pair(2, cReturnMapHandle[n->BConnection[0].returnMapHandle[1]]));
							returnGeneralMap.push_back(tmpReturnMapI);
							return std::make_pair(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, returnGeneralMap);
						}
					}
				}
		}
		else
		{
				auto aConnection = MultiplyOperationNodeInternal(*(n->AConnection.entryPointHandle),
				 					maxLevel == n->level ? 'A' : position, maxLevel, n->AConnection.returnMapHandle);
				CFLOBDDReturnMapHandle mi;
				g->numBConnections = aConnection.first.handleContents->numExits;
				g->BConnection = new Connection[g->numBConnections];
				for (unsigned int i = 0; i < g->numBConnections; i++)
					mi.AddToEnd(i);
				mi.Canonicalize();
				g->AConnection = Connection(aConnection.first, mi);
				CFLOBDDNodeHandle* tmpBConnections = new CFLOBDDNodeHandle[n->numBConnections];
				std::vector<std::vector<std::vector<std::pair<int, int>>>> tmpBConnectionReturnMaps;

				for (unsigned int i = 0; i < n->numBConnections; i++)
				{
					CFLOBDDReturnMapHandle m;
					for (unsigned int j = 0; j < n->BConnection[i].returnMapHandle.Size(); j++)
						m.AddToEnd(j);
					m.Canonicalize();
					auto tmpReturn = MultiplyOperationNodeInternal(*(n->BConnection[i].entryPointHandle), maxLevel == n->level ? 'B' : position, maxLevel, m);
					tmpBConnections[i] = tmpReturn.first;
					//if (n->level > 1)
						changeGeneralMap(tmpReturn.second, n->BConnection[i].returnMapHandle);
					tmpBConnectionReturnMaps.push_back(tmpReturn.second);
				}

				std::map<std::vector<std::pair<int, int>>, int> mapOfGeneralToReturn;
				std::vector<std::vector<std::pair<int, int>>> returnGeneralMap;

				for (unsigned int i = 0; i < aConnection.second.size(); i++){
					CFLOBDDNodeHandle tmpBConnectionI = tmpBConnections[aConnection.second[i][0].second];
					auto tmpBConnectionIGeneralMap = multiplyGeneralMap(aConnection.second[i][0].first, tmpBConnectionReturnMaps[aConnection.second[i][0].second]);
					for (unsigned int j = 1; j < aConnection.second[i].size(); j++){
						auto tmpBConnectionJGeneralMap = multiplyGeneralMap(aConnection.second[i][j].first, tmpBConnectionReturnMaps[aConnection.second[i][j].second]);
						auto nodeAndGeneralMap = addNodes(tmpBConnectionI, tmpBConnections[aConnection.second[i][j].second], tmpBConnectionIGeneralMap, tmpBConnectionJGeneralMap);
						tmpBConnectionIGeneralMap = nodeAndGeneralMap.second;
						tmpBConnectionI = nodeAndGeneralMap.first;
					}
					CFLOBDDReturnMapHandle returnMapHandle;
					for (unsigned int k = 0; k < tmpBConnectionIGeneralMap.size(); k++)
					{
						if (mapOfGeneralToReturn.find(tmpBConnectionIGeneralMap[k]) == mapOfGeneralToReturn.end())
						{
							mapOfGeneralToReturn[tmpBConnectionIGeneralMap[k]] = returnGeneralMap.size();
							returnMapHandle.AddToEnd(returnGeneralMap.size());
							returnGeneralMap.push_back(tmpBConnectionIGeneralMap[k]);
						}
						else{
							returnMapHandle.AddToEnd(mapOfGeneralToReturn[tmpBConnectionIGeneralMap[k]]);
						}
					}
					returnMapHandle.Canonicalize();
					g->BConnection[i] = Connection(tmpBConnectionI, returnMapHandle);
				}

				g->numExits = returnGeneralMap.size();
#ifdef PATH_COUNTING_ENABLED
				g->InstallPathCounts();
#endif
				return std::pair<CFLOBDDNodeHandle, std::vector<std::vector<std::pair<int, int>>>>(CFLOBDDNodeHandle(g), returnGeneralMap);
		}
	}

	void changeGeneralMap(std::vector<std::vector<std::pair<int, int>>> &generalMap, CFLOBDDReturnMapHandle returnMapHandle){
		for (unsigned int i = 0; i < generalMap.size(); i++){
			for (unsigned int j = 0; j < generalMap[i].size(); j++){
				generalMap[i][j].second = returnMapHandle[generalMap[i][j].second];
			}
			std::sort(generalMap[i].begin(), generalMap[i].end(), [](std::pair<int,int> &left, std::pair<int,int> &right){
				return left.second < right.second;
			});
		}
	}

	std::vector<std::vector<std::pair<int, int>>> multiplyGeneralMap(int multiple, std::vector<std::vector<std::pair<int, int>>> &tmpReturnMap)
	{
		auto returnMap = tmpReturnMap;
		for (int i = 0; i < tmpReturnMap.size(); i++)
		{
			for (int j = 0; j < tmpReturnMap[i].size(); j++)
			{
				returnMap[i][j].first *= multiple;
			}
		}
		return returnMap;
	}

	std::vector<std::pair<int, int>> addGeneralMap(std::vector<std::pair<int, int>> &n1, std::vector<std::pair<int, int>>& n2)
	{
		std::vector<std::pair<int, int>> ans;
		unsigned int p1 = 0, p2 = 0;

		while (p1 < n1.size() && p2 < n2.size())
		{
			if (n1[p1].second < n2[p2].second){
				ans.push_back(n1[p1]);
				p1++;
			}
			else if (n1[p1].second > n2[p2].second){
				ans.push_back(n2[p2]);
				p2++;
			}else{
				ans.push_back(std::make_pair(n1[p1].first + n2[p2].first, n1[p1].second));
				p1++; p2++;
			}
		}

		while (p1 < n1.size())
		{
			ans.push_back(n1[p1]);
			p1++;
		}

		while (p2 < n2.size())
		{
			ans.push_back(n2[p2]);
			p2++;
		}

		return ans;
	}

	std::pair<CFLOBDDNodeHandle, std::vector<std::vector<std::pair<int, int>>>> addNodes(CFLOBDDNodeHandle n1, CFLOBDDNodeHandle n2, std::vector<std::vector<std::pair<int, int>>> n1GeneralMap, std::vector<std::vector<std::pair<int, int>>> n2GeneralMap)
	{
		// Perform 2-way cross product of n1 and n2
		PairProductMapHandle MapHandle;
		CFLOBDDNodeHandle n = PairProduct(n1, n2, MapHandle);

		// Create returnMapHandle from MapHandle: Fold the pairs in MapHandle by applying
		// [n1->rootConnection.returnMapHandle, n2->rootConnection.returnMapHandle]
		// (component-wise) to each pair.

		CFLOBDDReturnMapHandle returnMapHandle;
		//PairProductMapBodyIterator MapIterator(*MapHandle.mapContents);
		//MapIterator.Reset();

		std::map<std::vector<std::pair<int, int>>, int> mapOfGeneralToReturn;
		std::vector<std::vector<std::pair<int, int>>> returnGeneralMap;
		unsigned int iterator = 0;
		//while (!MapIterator.AtEnd()) {
		while (iterator < MapHandle.Size()){
			std::vector<std::pair<int, int>> c1, c2;
			int first, second;
			//first = MapIterator.Current().First();
			//second = MapIterator.Current().Second();
			first = MapHandle[iterator].First();
			second = MapHandle[iterator].Second();
			c1 = n1GeneralMap[first];
			c2 = n2GeneralMap[second];
			auto ans = addGeneralMap(c1, c2);
			if (mapOfGeneralToReturn.find(ans) == mapOfGeneralToReturn.end())
			{
				mapOfGeneralToReturn[ans] = returnGeneralMap.size();
				returnMapHandle.AddToEnd(returnGeneralMap.size());
				returnGeneralMap.push_back(ans);
			}
			else{
				returnMapHandle.AddToEnd(mapOfGeneralToReturn[ans]);
			}
			
			//MapIterator.Next();
			iterator++;
		}
		
		returnMapHandle.Canonicalize();

		// Perform reduction on n, with respect to the common elements that returnMapHandle maps together
		ReductionMapHandle inducedReductionMapHandle;
		CFLOBDDReturnMapHandle inducedReturnMap;
		returnMapHandle.InducedReductionAndReturnMap(inducedReductionMapHandle, inducedReturnMap);
		//CFLOBDDNodeHandle::InitReduceCache();
		CFLOBDDNodeHandle reduced_n = n.Reduce(inducedReductionMapHandle, inducedReturnMap.Size());
		//CFLOBDDNodeHandle::DisposeOfReduceCache();
		// Create and return CFLOBDDTopNode
		return std::pair<CFLOBDDNodeHandle, std::vector<std::vector<std::pair<int, int>>>>(reduced_n, returnGeneralMap);
	}

	CFLOBDDNodeHandle MkMatrixMultiplyConstraintNode(unsigned int level)
	{
		CFLOBDDInternalNode* n = new CFLOBDDInternalNode(level);
		CFLOBDDReturnMapHandle mI;
		for (unsigned int i = 0; i < pow(2, pow(2, level - 2)); i++)
			mI.AddToEnd(i);
		mI.Canonicalize();
		CFLOBDDNodeHandle atmp = MkMatrixMultiplyConstraintNodeInternal(level - 1);
		n->AConnection = Connection(atmp, mI);
		n->numBConnections = pow(2, pow(2, level - 2));
		n->BConnection = new Connection[n->numBConnections];
		std::map<std::pair<int, int>, CFLOBDDNodeHandle>  memoTable;
		for (unsigned int i = 0; i < n->numBConnections; i++){
			CFLOBDDReturnMapHandle m;
			if (i == 0){
				m.AddToEnd(0);
				m.AddToEnd(1);
			}
			else{
				m.AddToEnd(1);
				m.AddToEnd(0);
			}
			m.Canonicalize();
			CFLOBDDNodeHandle btmp = MkRowWithOne(pow(2, level - 2), i, memoTable);
			n->BConnection[i] = Connection(btmp, m);
		}
		n->numExits = 2;
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkMatrixMultiplyConstraintNodeInternal(unsigned int level)
	{
		CFLOBDDInternalNode* n = new CFLOBDDInternalNode(level);
		if (level == 1)
		{
			CFLOBDDReturnMapHandle m0;
			m0.AddToEnd(0);
			m0.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
			n->numBConnections = 1;
			n->BConnection = new Connection[n->numBConnections];
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
			n->numExits = 2;
		}
		else
		{
			CFLOBDDNodeHandle tmpHandle = MkMatrixMultiplyConstraintNodeInternal(level - 1);
			CFLOBDDReturnMapHandle mI;
			for (unsigned int i = 0; i < pow(2, pow(2, level - 2)); i++)
				mI.AddToEnd(i);
			mI.Canonicalize();
			n->AConnection = Connection(tmpHandle, mI);
			n->numBConnections = pow(2, pow(2, level - 2));
			n->BConnection = new Connection[n->numBConnections];
			for (unsigned int i = 0; i < n->numBConnections; i++)
			{
				CFLOBDDReturnMapHandle mi;
				for (unsigned int j = 0; j < n->numBConnections; j++)
					mi.AddToEnd(i*n->numBConnections + j);
				mi.Canonicalize();
				n->BConnection[i] = Connection(tmpHandle, mi);
			}
			n->numExits = pow(2, pow(2, level - 1));
		}
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		return CFLOBDDNodeHandle(n);
	}

	CFLOBDDNodeHandle MkRowWithOne(unsigned int bits, unsigned int rowNo, std::map<std::pair<int, int>, CFLOBDDNodeHandle>&  memoTable)
	{
		if (memoTable.find(std::make_pair(bits, rowNo)) != memoTable.end())
			return memoTable[std::make_pair(bits, rowNo)];
		CFLOBDDInternalNode* n = new CFLOBDDInternalNode(log2(bits) + 1);
		if (bits == 1)
		{
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			n->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
			n->numBConnections = 2;
			n->BConnection = new Connection[2];
			CFLOBDDReturnMapHandle m0, m1;
			m0.AddToEnd(0);
			m0.Canonicalize();
			m1.AddToEnd(1);
			m1.Canonicalize();
			n->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
			n->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);
			n->numExits = 2;
		}
		else
		{
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0);
			m01.AddToEnd(1);
			m01.Canonicalize();
			CFLOBDDNodeHandle atmp = MkRowWithOne(bits / 2, rowNo / bits, memoTable);
			n->AConnection = Connection(atmp, m01);
			n->numBConnections = 2;
			n->BConnection = new Connection[n->numBConnections];
			if (rowNo < bits){
				CFLOBDDNodeHandle btmp = MkRowWithOne(bits / 2, rowNo % bits, memoTable);
				n->BConnection[0] = Connection(btmp, m01);
				CFLOBDDReturnMapHandle m;
				if (rowNo % bits == 0)
					m.AddToEnd(1);
				else
					m.AddToEnd(0);
				m.Canonicalize();
				n->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[(int)log2(bits)], m);
				n->numExits = 2;
			}
			else{
				CFLOBDDReturnMapHandle m;
				if (rowNo % bits == 0){
					m.AddToEnd(1);
					m.AddToEnd(0);
				}
				else{
					m.AddToEnd(0);
					m.AddToEnd(1);
				}
				m.Canonicalize();
				CFLOBDDNodeHandle btmp = MkRowWithOne(bits / 2, rowNo % bits, memoTable);
				n->BConnection[1] = Connection(btmp, m);
				CFLOBDDReturnMapHandle m0;
				m0.AddToEnd(0);
				m0.Canonicalize();
				n->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[(int)log2(bits)], m0);
				n->numExits = 2;
			}
		}
#ifdef PATH_COUNTING_ENABLED
		n->InstallPathCounts();
#endif
		CFLOBDDNodeHandle nhHandle(n);
		memoTable[std::make_pair(bits, rowNo)] = nhHandle;
		return nhHandle;
	}

	// Naive Matrix Multiplication
	
	CFLOBDDTopNodeMatMultMapRefPtr
		MatrixMultiplyV4Node(std::unordered_map<MatMultPair, CFLOBDDTopNodeMatMultMapRefPtr, MatMultPair::MatMultPairHash>& hashMap,
		CFLOBDDNodeHandle c1, CFLOBDDNodeHandle c2)
	{

		MatMultPair mmp(c1, c2);
		/*if (hashMap.find(mmp) != hashMap.end()){
			return hashMap[mmp];
		}*/
		
		if (matmult_hashMap.find(mmp) != matmult_hashMap.end() && 
			((((CFLOBDDInternalNode *)(c1.handleContents))->GetRefCount() >= 2) || 
			(((CFLOBDDInternalNode *)(c2.handleContents))->GetRefCount() >= 2))){
			return matmult_hashMap[mmp];
		}

		//std::cout << "node: " << c1.handleContents->level << std::endl;

		CFLOBDDMatMultMapHandle g_return_map;
		CFLOBDDInternalNode* g = new CFLOBDDInternalNode(c1.handleContents->level);
		ReductionMapHandle reductionMapHandle;
		
		if (c1.handleContents->level == 1){
			CFLOBDDInternalNode* c1_internal = (CFLOBDDInternalNode *)c1.handleContents;
			CFLOBDDInternalNode* c2_internal = (CFLOBDDInternalNode *)c2.handleContents;
			std::vector<unsigned int> m1_indices{ 0, 0, 0, 0 }, m2_indices{ 0, 0, 0, 0 };
			for (unsigned int i = 0; i < c1_internal->numBConnections; i++){
				if (*(c1_internal->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle){
					m1_indices[2 * i] = c1_internal->BConnection[i].returnMapHandle[0];
					m1_indices[2 * i + 1] = c1_internal->BConnection[i].returnMapHandle[0];
				}
				else{
					m1_indices[2 * i] = c1_internal->BConnection[i].returnMapHandle[0];
					m1_indices[2 * i + 1] = c1_internal->BConnection[i].returnMapHandle[1];
				}
			}

			if (c1_internal->numBConnections == 1){
				m1_indices[2] = m1_indices[0];
				m1_indices[3] = m1_indices[1];
			}

			for (unsigned int i = 0; i < c2_internal->numBConnections; i++){
				if (*(c2_internal->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle){
					m2_indices[2 * i] = c2_internal->BConnection[i].returnMapHandle[0];
					m2_indices[2 * i + 1] = c2_internal->BConnection[i].returnMapHandle[0];
				}
				else{
					m2_indices[2 * i] = c2_internal->BConnection[i].returnMapHandle[0];
					m2_indices[2 * i + 1] = c2_internal->BConnection[i].returnMapHandle[1];
				}
			}

			if (c2_internal->numBConnections == 1){
				m2_indices[2] = m2_indices[0];
				m2_indices[3] = m2_indices[1];
			}

			MatMultMapHandle m0, m1, m2, m3;
			VAL_TYPE one = 1;
			m0.Add(std::make_pair(m1_indices[0], m2_indices[0]), one);
			m0.Add(std::make_pair(m1_indices[1], m2_indices[2]), one);
			m0.Canonicalize();

			m1.Add(std::make_pair(m1_indices[0], m2_indices[1]), one);
			m1.Add(std::make_pair(m1_indices[1], m2_indices[3]), one);
			m1.Canonicalize();

			m2.Add(std::make_pair(m1_indices[2], m2_indices[0]), one);
			m2.Add(std::make_pair(m1_indices[3], m2_indices[2]), one);
			m2.Canonicalize();

			m3.Add(std::make_pair(m1_indices[2], m2_indices[1]), one);
			m3.Add(std::make_pair(m1_indices[3], m2_indices[3]), one);
			m3.Canonicalize();
			if ((m0 == m2) && (m1 == m3)){
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
				g->numBConnections = 1;
				g->BConnection = new Connection[g->numBConnections];

				if (m0 == m1){
					g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
					g->numExits = 1;
					g_return_map.AddToEnd(m0);
					reductionMapHandle.AddToEnd(0);
				}
				else{
					g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
					g->numExits = 2;
					g_return_map.AddToEnd(m0);
					g_return_map.AddToEnd(m1);
					reductionMapHandle.AddToEnd(0);
					reductionMapHandle.AddToEnd(1);
				}
			}
			else{
				g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
				g->numBConnections = 2;
				g->BConnection = new Connection[g->numBConnections];


				if (m0 == m1){
					g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
					g_return_map.AddToEnd(m0);
					reductionMapHandle.AddToEnd(0);
				}
				else{
					g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
					g_return_map.AddToEnd(m0);
					g_return_map.AddToEnd(m1);
					reductionMapHandle.AddToEnd(0);
					reductionMapHandle.AddToEnd(1);
				}

				if (m2 == m3){
					unsigned int k;
					for (k = 0; k < g_return_map.Size(); k++){
						if (g_return_map[k] == m2)
							break;
					}
					CFLOBDDReturnMapHandle mtmp;
					mtmp.AddToEnd(k);
					mtmp.Canonicalize();
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, mtmp);
					if (k == g_return_map.Size()){
						g_return_map.AddToEnd(m2);
						reductionMapHandle.AddToEnd(k);
					}
				}
				else{
					unsigned int k1 = g_return_map.Size(), k2 = g_return_map.Size();
					for (unsigned int k = 0; k < g_return_map.Size(); k++){
						if (g_return_map[k] == m2){
							k1 = k;
						}
						if (g_return_map[k] == m3){
							k2 = k;
						}
					}
					CFLOBDDReturnMapHandle mtmp;
					mtmp.AddToEnd(k1);
					if (k1 == g_return_map.Size() && k2 == g_return_map.Size())
						k2++;
					mtmp.AddToEnd(k2);
					mtmp.Canonicalize();
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, mtmp);
					if (k1 == g_return_map.Size()){
						g_return_map.AddToEnd(m2);
						reductionMapHandle.AddToEnd(k1);
					}
					if (k2 == g_return_map.Size()){
						g_return_map.AddToEnd(m3);
						reductionMapHandle.AddToEnd(k2);
					}
				}
			}

		}
		else{
			
			CFLOBDDInternalNode* c1_internal = (CFLOBDDInternalNode *)c1.handleContents;
			CFLOBDDInternalNode* c2_internal = (CFLOBDDInternalNode *)c2.handleContents;

			CFLOBDDTopNodeMatMultMapRefPtr aa = MatrixMultiplyV4Node(hashMap, *(c1_internal->AConnection.entryPointHandle),
				*(c2_internal->AConnection.entryPointHandle));
			CFLOBDDReturnMapHandle mI;
			for (unsigned int i = 0; i < aa->rootConnection.returnMapHandle.Size(); i++)
				mI.AddToEnd(i);
			mI.Canonicalize();

			g->AConnection = Connection(*(aa->rootConnection.entryPointHandle), mI);
			g->numBConnections = mI.Size();
			g->BConnection = new Connection[g->numBConnections];
			g->numExits = 0;
			//std::unordered_map<std::string, unsigned int> mapFromHandleToIndex;
			std::unordered_map<unsigned int, unsigned int> mapFromHandleToIndex;
			for (unsigned int i = 0; i < g->numBConnections; i++){
				MatMultMapHandle matmult_returnmap = aa->rootConnection.returnMapHandle[i];
				CFLOBDDTopNodeMatMultMapRefPtr ans;
				bool first = true;
				// Consider Multiplication of M1 and M2
				for (auto &v : matmult_returnmap.mapContents->map){
					unsigned int M1_index = v.first.first;
					unsigned int M2_index = v.first.second;
					//VAL_TYPE factor = v.second;
					CFLOBDDTopNodeMatMultMapRefPtr bb_old = 
						MatrixMultiplyV4Node(hashMap, *(c1_internal->BConnection[M1_index].entryPointHandle),
						*(c2_internal->BConnection[M2_index].entryPointHandle));
					CFLOBDDMatMultMapHandle new_bb_return;
					for (unsigned int j = 0; j < bb_old->rootConnection.returnMapHandle.Size(); j++)
					{
						MatMultMapHandle tmp;
						for (auto& it : bb_old->rootConnection.returnMapHandle[j].mapContents->map){
							tmp.Add(std::make_pair(c1_internal->BConnection[M1_index].returnMapHandle[it.first.first],
								c2_internal->BConnection[M2_index].returnMapHandle[it.first.second]), it.second);
						}
						tmp.Canonicalize();
						new_bb_return.AddToEnd(tmp);
					}
					new_bb_return.Canonicalize();
					CFLOBDDTopNodeMatMultMapRefPtr bb =
						new CFLOBDDTopNodeMatMultMap(*(bb_old->rootConnection.entryPointHandle), new_bb_return);
					bb = MkLeftScalarTimesTopNode<MatMultMapHandle, VAL_TYPE>(v.second, bb);
					// bb = v.second * bb;
					if (first){
						ans = bb;
						first = false;
					}
					else{
						ans = MkPlusTopNode<MatMultMapHandle>(ans, bb);
						// ans = ans + bb;
					}
				}

				CFLOBDDReturnMapHandle ans_return_map;
				for (unsigned int j = 0; j < ans->rootConnection.returnMapHandle.Size(); j++){
					//std::string map_as_string = ans->rootConnection.returnMapHandle[j].ToString();
					unsigned int map_hash_check = ans->rootConnection.returnMapHandle[j].mapContents->hashCheck;
					if (mapFromHandleToIndex.find(map_hash_check) == mapFromHandleToIndex.end()){
						ans_return_map.AddToEnd(g->numExits++);
						g_return_map.AddToEnd(ans->rootConnection.returnMapHandle[j]);
						reductionMapHandle.AddToEnd(g->numExits - 1);
						mapFromHandleToIndex[map_hash_check] = g->numExits - 1;
					}
					else{
						unsigned int index = mapFromHandleToIndex[map_hash_check];
						ans_return_map.AddToEnd(index);
					}
				}
				ans_return_map.Canonicalize();
				g->BConnection[i] = Connection(*(ans->rootConnection.entryPointHandle), ans_return_map);
			}

		}
		g_return_map.Canonicalize();
		reductionMapHandle.Canonicalize();
		g->numExits = g_return_map.Size();

#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif

		CFLOBDDNodeHandle gHandle(g);
		gHandle = gHandle.Reduce(reductionMapHandle, g_return_map.Size(), true);
		CFLOBDDTopNodeMatMultMapRefPtr return_ans = new CFLOBDDTopNodeMatMultMap(gHandle, g_return_map);
		//hashMap[mmp] = return_ans;
		matmult_hashMap[mmp] = return_ans;
		return return_ans;
	}

	CFLOBDDTopNodeMatMultMapRefPtr MatrixMultiplyV4WithInfoNode(
		std::unordered_map<ZeroValNodeInfo, ZeroIndicesMapHandle, ZeroValNodeInfo::ZeroValNodeInfoHash>& hashMap,
		CFLOBDDNodeHandle c1, CFLOBDDNodeHandle c2, int c1_zero_index, int c2_zero_index)
	{
		MatMultPairWithInfo mmp(c1, c2, c1_zero_index, c2_zero_index);

		if (matmult_hashMap_info.find(mmp) != matmult_hashMap_info.end() 
			&& ((((CFLOBDDInternalNode *)(c1.handleContents))->GetRefCount() >= 2) ||
			(((CFLOBDDInternalNode *)(c2.handleContents))->GetRefCount() >= 2))
			){
			return matmult_hashMap_info[mmp];
		}

		CFLOBDDMatMultMapHandle g_return_map;
		CFLOBDDInternalNode* g = new CFLOBDDInternalNode(c1.handleContents->level);

		CFLOBDDInternalNode* c1_internal = (CFLOBDDInternalNode *)c1.handleContents;
		CFLOBDDInternalNode* c2_internal = (CFLOBDDInternalNode *)c2.handleContents;

		ReductionMapHandle reductionMapHandle;

		if ((c1_internal->numExits == 1 && c1_zero_index == 0) || (c2_internal->numExits == 1 && c2_zero_index == 0))
		{
			g = (CFLOBDDInternalNode *)CFLOBDDNodeHandle::NoDistinctionNode[c1.handleContents->level].handleContents;
			MatMultMapHandle m;
			VAL_TYPE one = 1;
			m.ForceAdd(std::make_pair(-1, -1), one);
			m.Canonicalize();
			g_return_map.AddToEnd(m);
			reductionMapHandle.AddToEnd(0);
		}
		else{
			if (c1.handleContents->level == 1){
				std::vector<int> m1_indices{ 0, 0, 0, 0 }, m2_indices{ 0, 0, 0, 0 };
				for (unsigned int i = 0; i < c1_internal->numBConnections; i++){
					if (*(c1_internal->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle){
						m1_indices[2 * i] = c1_internal->BConnection[i].returnMapHandle[0];
						m1_indices[2 * i + 1] = c1_internal->BConnection[i].returnMapHandle[0];
					}
					else{
						m1_indices[2 * i] = c1_internal->BConnection[i].returnMapHandle[0];
						m1_indices[2 * i + 1] = c1_internal->BConnection[i].returnMapHandle[1];
					}
				}

				if (c1_internal->numBConnections == 1){
					m1_indices[2] = m1_indices[0];
					m1_indices[3] = m1_indices[1];
				}

				for (unsigned int i = 0; i < c2_internal->numBConnections; i++){
					if (*(c2_internal->BConnection[i].entryPointHandle) == CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle){
						m2_indices[2 * i] = c2_internal->BConnection[i].returnMapHandle[0];
						m2_indices[2 * i + 1] = c2_internal->BConnection[i].returnMapHandle[0];
					}
					else{
						m2_indices[2 * i] = c2_internal->BConnection[i].returnMapHandle[0];
						m2_indices[2 * i + 1] = c2_internal->BConnection[i].returnMapHandle[1];
					}
				}

				if (c2_internal->numBConnections == 1){
					m2_indices[2] = m2_indices[0];
					m2_indices[3] = m2_indices[1];
				}

				for (unsigned int i = 0; i < m1_indices.size(); i++){
					if (m1_indices[i] == c1_zero_index)
						m1_indices[i] = -1;
				}

				for (unsigned int i = 0; i < m2_indices.size(); i++){
					if (m2_indices[i] == c2_zero_index)
						m2_indices[i] = -1;
				}

				MatMultMapHandle m0, m1, m2, m3;
				VAL_TYPE one = 1;
				m0.Add(std::make_pair(m1_indices[0], m2_indices[0]), one);
				m0.Add(std::make_pair(m1_indices[1], m2_indices[2]), one);
				if (m0.Size() == 0)
					m0.ForceAdd(std::make_pair(-1, -1), one);
				m0.Canonicalize();

				m1.Add(std::make_pair(m1_indices[0], m2_indices[1]), one);
				m1.Add(std::make_pair(m1_indices[1], m2_indices[3]), one);
				if (m1.Size() == 0)
					m1.ForceAdd(std::make_pair(-1, -1), one);
				m1.Canonicalize();

				m2.Add(std::make_pair(m1_indices[2], m2_indices[0]), one);
				m2.Add(std::make_pair(m1_indices[3], m2_indices[2]), one);
				if (m2.Size() == 0)
					m2.ForceAdd(std::make_pair(-1, -1), one);
				m2.Canonicalize();

				m3.Add(std::make_pair(m1_indices[2], m2_indices[1]), one);
				m3.Add(std::make_pair(m1_indices[3], m2_indices[3]), one);
				if (m3.Size() == 0)
					m3.ForceAdd(std::make_pair(-1, -1), one);
				m3.Canonicalize();
				if ((m0 == m2) && (m1 == m3)){
					g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
					g->numBConnections = 1;
					g->BConnection = new Connection[g->numBConnections];

					if (m0 == m1){
						g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
						g->numExits = 1;
						g_return_map.AddToEnd(m0);
						reductionMapHandle.AddToEnd(0);
					}
					else{
						g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
						g->numExits = 2;
						g_return_map.AddToEnd(m0);
						g_return_map.AddToEnd(m1);
						reductionMapHandle.AddToEnd(0);
						reductionMapHandle.AddToEnd(1);
					}
				}
				else{
					g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
					g->numBConnections = 2;
					g->BConnection = new Connection[g->numBConnections];


					if (m0 == m1){
						g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, commonly_used_return_maps[0]);//m0
						g_return_map.AddToEnd(m0);
						reductionMapHandle.AddToEnd(0);
					}
					else{
						g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, commonly_used_return_maps[2]);//m01
						g_return_map.AddToEnd(m0);
						g_return_map.AddToEnd(m1);
						reductionMapHandle.AddToEnd(0);
						reductionMapHandle.AddToEnd(1);
					}

					if (m2 == m3){
						unsigned int k;
						for (k = 0; k < g_return_map.Size(); k++){
							if (g_return_map[k] == m2)
								break;
						}
						CFLOBDDReturnMapHandle mtmp;
						mtmp.AddToEnd(k);
						mtmp.Canonicalize();
						g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, mtmp);
						if (k == g_return_map.Size()){
							g_return_map.AddToEnd(m2);
							reductionMapHandle.AddToEnd(k);
						}
					}
					else{
						unsigned int k1 = g_return_map.Size(), k2 = g_return_map.Size();
						for (unsigned int k = 0; k < g_return_map.Size(); k++){
							if (g_return_map[k] == m2){
								k1 = k;
							}
							if (g_return_map[k] == m3){
								k2 = k;
							}
						}
						CFLOBDDReturnMapHandle mtmp;
						mtmp.AddToEnd(k1);
						if (k1 == g_return_map.Size() && k2 == g_return_map.Size())
							k2++;
						mtmp.AddToEnd(k2);
						mtmp.Canonicalize();
						g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, mtmp);
						if (k1 == g_return_map.Size()){
							g_return_map.AddToEnd(m2);
							reductionMapHandle.AddToEnd(k1);
						}
						if (k2 == g_return_map.Size()){
							g_return_map.AddToEnd(m3);
							reductionMapHandle.AddToEnd(k2);
						}
					}
				}

			}
			else{

				std::vector<int> b_zero_indices_c1(c1_internal->numBConnections, -1), b_zero_indices_c2(c2_internal->numBConnections, -1);
				int a_zero_index_c1 = -1, a_zero_index_c2 = -1;
				if (c1_zero_index != -1){
					ZeroValNodeInfo c1_zero_val_node(c1, c1_zero_index);
					auto it = hashMap.find(c1_zero_val_node);
					if (it != hashMap.end()){
						b_zero_indices_c1 = it->second.mapContents->b_indices;
						a_zero_index_c1 = it->second.mapContents->a_index;
					}
					else{
						for (unsigned int i = 0; i < b_zero_indices_c1.size(); i++){
							b_zero_indices_c1[i] = c1_internal->BConnection[i].returnMapHandle.LookupInv(c1_zero_index);
							if (c1_internal->BConnection[i].returnMapHandle.Size() == 1 && c1_internal->BConnection[i].returnMapHandle[0] == c1_zero_index)
							{
								a_zero_index_c1 = i;
							}
						}
						ZeroIndicesMapHandle zm;
						zm.mapContents->b_indices = b_zero_indices_c1;
						zm.Set_AIndex(a_zero_index_c1);
						zm.Canonicalize();
						hashMap.emplace(c1_zero_val_node, zm);
					}
				}

				if (c2_zero_index != -1){
					ZeroValNodeInfo c2_zero_val_node(c2, c2_zero_index);
					auto it = hashMap.find(c2_zero_val_node);
					if (it != hashMap.end()){
						b_zero_indices_c2 = it->second.mapContents->b_indices;
						a_zero_index_c2 = it->second.mapContents->a_index;
					}
					else{
						for (unsigned int i = 0; i < b_zero_indices_c2.size(); i++){
							b_zero_indices_c2[i] = c2_internal->BConnection[i].returnMapHandle.LookupInv(c2_zero_index);
							if (c2_internal->BConnection[i].returnMapHandle.Size() == 1 && c2_internal->BConnection[i].returnMapHandle[0] == c2_zero_index)
							{
								a_zero_index_c2 = i;
							}
						}
						ZeroIndicesMapHandle zm;
						zm.mapContents->b_indices = b_zero_indices_c2;
						zm.Set_AIndex(a_zero_index_c2);
						zm.Canonicalize();
						hashMap.emplace(c2_zero_val_node, zm);
					}
				}

				CFLOBDDTopNodeMatMultMapRefPtr aa = MatrixMultiplyV4WithInfoNode(hashMap, *(c1_internal->AConnection.entryPointHandle),
					*(c2_internal->AConnection.entryPointHandle), a_zero_index_c1, a_zero_index_c2);
				CFLOBDDReturnMapHandle mI;
				for (unsigned int i = 0; i < aa->rootConnection.returnMapHandle.Size(); i++)
					mI.AddToEnd(i);
				mI.Canonicalize();

				g->AConnection = Connection(*(aa->rootConnection.entryPointHandle), mI);
				g->numBConnections = mI.Size();
				g->BConnection = new Connection[g->numBConnections];
				g->numExits = 0;
				//std::unordered_map<std::string, unsigned int> mapFromHandleToIndex;
				std::unordered_map<unsigned int, unsigned int> mapFromHandleToIndex;
				for (unsigned int i = 0; i < g->numBConnections; i++){
					MatMultMapHandle matmult_returnmap = aa->rootConnection.returnMapHandle[i];
					CFLOBDDTopNodeMatMultMapRefPtr ans;
					bool first = true;
					if (matmult_returnmap.Size() == 1 &&
						(matmult_returnmap.mapContents->map.find(std::make_pair(-1, -1)) != matmult_returnmap.mapContents->map.end())){
						CFLOBDDMatMultMapHandle tmp_return_map;
						tmp_return_map.AddToEnd(matmult_returnmap);
						tmp_return_map.Canonicalize();
						ans = new CFLOBDDTopNodeMatMultMap(CFLOBDDNodeHandle::NoDistinctionNode[c1.handleContents->level - 1], tmp_return_map);
					}
					else{
						// Consider Multiplication of M1 and M2
						for (auto &v : matmult_returnmap.mapContents->map){
							unsigned int M1_index = v.first.first;
							unsigned int M2_index = v.first.second;
							//VAL_TYPE factor = v.second;
							CFLOBDDTopNodeMatMultMapRefPtr bb_old =
								MatrixMultiplyV4WithInfoNode(hashMap, *(c1_internal->BConnection[M1_index].entryPointHandle),
								*(c2_internal->BConnection[M2_index].entryPointHandle), b_zero_indices_c1[M1_index], b_zero_indices_c2[M2_index]);
							CFLOBDDMatMultMapHandle new_bb_return;
							for (unsigned int j = 0; j < bb_old->rootConnection.returnMapHandle.Size(); j++)
							{
								MatMultMapHandle tmp;
								for (auto& it : bb_old->rootConnection.returnMapHandle[j].mapContents->map){
									if (it.first.first != -1 && it.first.second != -1)
										tmp.Add(std::make_pair(c1_internal->BConnection[M1_index].returnMapHandle[it.first.first],
											c2_internal->BConnection[M2_index].returnMapHandle[it.first.second]), it.second);
								}
								if (tmp.Size() == 0){
									VAL_TYPE one = 1;
									tmp.ForceAdd(std::make_pair(-1, -1), one);
								}
								tmp.Canonicalize();
								new_bb_return.AddToEnd(tmp);
							}
							new_bb_return.Canonicalize();
							CFLOBDDTopNodeMatMultMapRefPtr bb =
								new CFLOBDDTopNodeMatMultMap(*(bb_old->rootConnection.entryPointHandle), new_bb_return);
							if (!(new_bb_return.Size() == 1 &&
								new_bb_return[0].mapContents->contains_zero_val))
								bb = MkLeftScalarTimesTopNode<MatMultMapHandle, VAL_TYPE>(v.second, bb);
								// bb = v.second * bb;
							if (first){
								ans = bb;
								first = false;
							}
							else{
								if (!(ans->rootConnection.returnMapHandle.Size() == 1 &&
									ans->rootConnection.returnMapHandle[0].mapContents->contains_zero_val))
									// ans = ans + bb;
									ans = MkPlusTopNode<MatMultMapHandle>(ans, bb);
								else
									ans = bb;
							}
						}
					}
					CFLOBDDReturnMapHandle ans_return_map;
					for (unsigned int j = 0; j < ans->rootConnection.returnMapHandle.Size(); j++){
						//std::string map_as_string = ans->rootConnection.returnMapHandle[j].ToString();
						unsigned int map_hash_check = ans->rootConnection.returnMapHandle[j].mapContents->hashCheck;
						if (mapFromHandleToIndex.find(map_hash_check) == mapFromHandleToIndex.end()){
							ans_return_map.AddToEnd(g->numExits++);
							g_return_map.AddToEnd(ans->rootConnection.returnMapHandle[j]);
							reductionMapHandle.AddToEnd(g->numExits - 1);
							mapFromHandleToIndex[map_hash_check] = g->numExits - 1;
						}
						else{
							unsigned int index = mapFromHandleToIndex[map_hash_check];
							ans_return_map.AddToEnd(index);
						}
					}
					ans_return_map.Canonicalize();
					g->BConnection[i] = Connection(*(ans->rootConnection.entryPointHandle), ans_return_map);
				}

			}
		}
		g_return_map.Canonicalize();
		g->numExits = g_return_map.Size();
		reductionMapHandle.Canonicalize();
#ifdef PATH_COUNTING_ENABLED
		g->InstallPathCounts();
#endif

		CFLOBDDNodeHandle gHandle(g);
		//gHandle = gHandle.Reduce(reductionMapHandle, g_return_map.Size(), true);
		CFLOBDDTopNodeMatMultMapRefPtr return_ans = new CFLOBDDTopNodeMatMultMap(gHandle, g_return_map);
		//hashMap[mmp] = return_ans;
		matmult_hashMap_info[mmp] = return_ans;
		return return_ans;
	}

	std::pair<CFLOBDDNodeHandle, int> MkRestrictNode(unsigned int level, std::string s)
	{
		if (s.find('0') == std::string::npos && s.find('1') == std::string::npos)
			return std::make_pair(CFLOBDDNodeHandle::NoDistinctionNode[level], -1);
		CFLOBDDInternalNode* g = new CFLOBDDInternalNode(level);
		if (level == 1)
		{
			assert(s.length() == 1);
			if (s[0] == 'X')
			{
				return std::make_pair(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, -1);
			}
			CFLOBDDReturnMapHandle m01;
			m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
			g->AConnection = Connection(CFLOBDDNodeHandle::CFLOBDDForkNodeHandle, m01);
			g->numBConnections = 2;
			g->BConnection = new Connection[g->numBConnections];
			CFLOBDDReturnMapHandle m0, m1;
			m0.AddToEnd(0); m0.Canonicalize();
			m1.AddToEnd(1); m1.Canonicalize();
			g->BConnection[0] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m0);
			g->BConnection[1] = Connection(CFLOBDDNodeHandle::CFLOBDDDontCareNodeHandle, m1);
			g->numExits = 2;
			return std::make_pair(CFLOBDDNodeHandle(g), (s[0] == '0') ? 1 : 0);
		}
		else
		{
			auto aa = MkRestrictNode(level-1, s.substr(0, s.length()/2));
			CFLOBDDReturnMapHandle mI;
			for (unsigned int i = 0; i < aa.first.handleContents->numExits; i++)
				mI.AddToEnd(i);
			mI.Canonicalize();
			g->AConnection = Connection(aa.first, mI);
			g->numBConnections = mI.Size();
			g->BConnection = new Connection[g->numBConnections];
			if (aa.second == -1){
				auto bb = MkRestrictNode(level-1, s.substr(s.length()/2));
				CFLOBDDReturnMapHandle m;
				for (unsigned int i = 0; i < bb.first.handleContents->numExits; i++)
					m.AddToEnd(i);
				m.Canonicalize();
				g->BConnection[0] = Connection(bb.first, m);
				g->numExits = m.Size();
				return std::make_pair(CFLOBDDNodeHandle(g), bb.second);
			}
			else if (aa.second == 0)
			{
				CFLOBDDReturnMapHandle m0;
				m0.AddToEnd(0); m0.Canonicalize();
				g->BConnection[0] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
				auto bb = MkRestrictNode(level-1, s.substr(s.length()/2));
				if (bb.second == 0){
					CFLOBDDReturnMapHandle m01;
					m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
					g->BConnection[1] = Connection(bb.first, m01);
				}
				else if (bb.second == 1)
				{
					CFLOBDDReturnMapHandle m10;
					m10.AddToEnd(1); m10.AddToEnd(0); m10.Canonicalize();
					g->BConnection[1] = Connection(bb.first, m10);
				}
				else if (bb.second == -1)
				{
					CFLOBDDReturnMapHandle m1;
					m1.AddToEnd(1);
					m1.Canonicalize();
					g->BConnection[1] = Connection(bb.first, m1);
				}
				g->numExits = 2;
				return std::make_pair(CFLOBDDNodeHandle(g), 0);
			}
			else
			{
				int index = -1;
				CFLOBDDReturnMapHandle m01;
				m01.AddToEnd(0); m01.AddToEnd(1); m01.Canonicalize();
				auto bb = MkRestrictNode(level-1, s.substr(s.length()/2));
				if (bb.second != -1){
					g->BConnection[0] = Connection(bb.first, m01);
					index = bb.second;
				}
				else
				{
					CFLOBDDReturnMapHandle m0;
					m0.AddToEnd(0); m0.Canonicalize();
					g->BConnection[0] = Connection(bb.first, m0);
					index = 1;
				}
				if (bb.second == 0){	
					CFLOBDDReturnMapHandle m0;
					m0.AddToEnd(0); m0.Canonicalize();
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m0);
				}
				else if (bb.second == 1)
				{
					CFLOBDDReturnMapHandle m1;
					m1.AddToEnd(1); m1.Canonicalize();
					g->BConnection[1] = Connection(CFLOBDDNodeHandle::NoDistinctionNode[level-1], m1);
				}
				else if (bb.second == -1)
				{
					CFLOBDDReturnMapHandle m1;
					m1.AddToEnd(1);
					m1.Canonicalize();
					g->BConnection[1] = Connection(bb.first, m1);
				}
				g->numExits = 2;
				return std::make_pair(CFLOBDDNodeHandle(g), index);
			}
		}
	}

}  // namespace CFL_OBDD
