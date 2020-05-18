#include <iostream>
#include <cstdint>
#include <string>
#include <vector>
#include <map>
#include <set>
#include <tuple>
#include <algorithm>
#include <utility>
#include <stdlib.h>     /* atoi */

#include "BCData.hpp"
#include "configBC.hpp"

using namespace std;
typedef unsigned int uint;

//Here are exemple functions on how to reproduce results
//For most of these, one first needs to generate a BCData object which contains information about the block cipher
//See the configBC.hpp and configBC.cpp files for examples

void computeNumberTrails(BCData & BCD,
						 vector<uint8_t> const & input,
						 vector<uint8_t> const & output,
						 vector<vector<uint8_t>> const & keyval){
	/*
		Compute the number of trails from #input to #output using the key #keyval
		the vectors should be of the appropriate sizes (i.e. #BCD.blockSize)
	*/

	auto nbTrail_fullOpti = BCD.countTrailsFullCipher(input,output,keyval);
	if(!nbTrail_fullOpti.second){
		cout << "Could not compute all trails, but found " << nbTrail_fullOpti.first << " trails with the given limits" << endl;
	}
	else{
		cout << nbTrail_fullOpti.first << " trails" << endl;
	}
}

void proveFullDegreeAllOutputBit(BCData & BCD){
	/*
		Try to prove full degree for each output bit of the block cipher represented by #BCD
		The number of rounds is fixed on the creation of the BCD object
	*/

	for(uint i = 0; i < BCD.blockSize; i++){
		//Search for an input/key that allows to get full degree
		auto iok = BCD.improvedDynamicSearch(i);

		//If we were able to do so, we get non-empty vectors
		if(get<0>(iok).size() > 0){
			cout << "Output bit " << i << " is of full degree using " << endl;
			cout << "Input : "; BCD.printParitySetHex(get<0>(iok)); cout << endl;
			cout << "Keys : " << endl;
			auto const & allKeys = get<2>(iok);
			for(auto const & key : allKeys){
				BCD.printParitySetHex(key); cout << endl;
			}
			//Output bit is fixed here so it does not need to be printed
		}
		else{
			cout << "Could not prove full degree for output bit " << i << endl;
		}
	}
}

void computeLowerBoundsAllOutputBitAllRounds(BCData & BCD){
	/*
		Try to compute lower bounds on the degree for each number of rounds between 1 and #BCD.rMax
	*/

	vector<vector<uint>> allLB(BCD.blockSize, vector<uint>(BCD.rMax+1,0));
	for(uint i = 0; i < BCD.blockSize; i++){
		//Getting the lower bounds is based on first proving the full degree for the max number of rounds (#BCD.rMax)
		auto iok = BCD.improvedDynamicSearch(i);

		if(get<0>(iok).size() > 0){
			//Once we found an input/key proving full degree for output bit i, we can comput the lower bounds for a smaller numbre of rounds
			auto lb = BCD.getLowerBounds(get<0>(iok),get<1>(iok),get<2>(iok));
			lb.emplace_back(BCD.blockSize-1);
			cout << "For output bit " << i << " : " << endl;
			for(uint r = 1; r <= BCD.rMax; r++)
				cout << r << " rounds, lower bound : " << lb[r] << endl;

			allLB[i] = std::move(lb);
		}
		else{
			cout << "Could not prove full degree for output bit " << i << endl;
		}
	}

	for(uint i = 0; i < BCD.blockSize; i++){
		cout << "For output bit " << i << " : " << endl;
		cout << "[" << allLB[i][1];
		for(uint j = 2; j <= BCD.rMax; j++)
			cout << "," << allLB[i][j];
		cout << "]" << endl;
	}
}

void proveFullMinDegreeAllBlocks(BCData & BCD){
	/*
		Try to prove the full min-degree for each block (Sbox or Super Sbox depending on where the key is added)
	*/
	BCD.searchMinDegree();
}

void proveFullMinDegreeAllBlocksAllInput(BCData & BCD){
	/*
		Try to prove the "Full monomial property" for each block (Sbox or Super Sbox depending on where the key is added)
		Essentially, it tries to prove full min-degree using each possible input monomial of degree #BCD.blockSize-1
	*/
	BCD.searchMinDegree_allInput();
}

void computeUpperBoundsAllOutput(BCData & BCD){
	/*
		Get the upper bounds on #BCD.rMax rounds for all output bits
		If one needs the upper bounds for different rounds, we need to recreate the BCD object, see next function
	*/
	BCD.getUpperBoundAllOutput();
}

void computeUpperBoundsAllOutputAllRounds(uint rMax){
	/*
		Get the upper bound for all output bits and for all number of rounds between 1 and #rMax (included)
	*/

	for(uint r = 1; r <= rMax; r++){
		auto BCD = genDataPRESENT(r); //Change the BCD generator function for another block cipher
		BCD.getUpperBoundAllOutput();
	}
}

string permuteNibbleSkinny(string const & s){
	string tmp(s);
	vector<uint> SR({0,5,10,15,4,9,14,3,8,13,2,7,12,1,6,11});
	for(uint i = 0; i < 16; i++)
		tmp[SR[i]] = s[i];
	return tmp;
}


int main(){
	
	
	//Generate the data for the block cipher (see configBC files)
	uint rMax = 10;
	auto BCD = genDataGIFT(rMax);
	// auto BCD = genDataPRESENT(rMax,true); 
	//For PRESENT our results on the min-degree were obtained by slightly changing the last sbox layer (while still maintaining correctness), hence the additional boolean

	//Do some stuff (just examples)

	//Check the number of trails
	vector<uint8_t>  input(hexToVec("FFFFFFEFFFFFFFFF"));
	vector<uint8_t> output(hexToVec("1000000000000000"));
	vector<vector<uint8_t>> keyval({
							hexToVec("0000000000000000"),
							hexToVec("0000000000000000"),
							hexToVec("0000000081004254"),
							hexToVec("2A82401244260E0D"),
							hexToVec("52A201812E100100"),
							hexToVec("0130100006100000"),
							hexToVec("0000000000000000"),
							hexToVec("0000000000000000"),
							hexToVec("0000000000000000")});
	computeNumberTrails(BCD,input,output,keyval); //With these values for Gift, should be 1 trail

	// //Prove the full degree for all output bits
	proveFullDegreeAllOutputBit(BCD);
	
	


	
	//Another exemple for Skinny to check the number of trails
	//Because the results for Skinny were obtained with a different implementation where the key was added before the SR operation (while here it is added after SR)
	//We need to permute the nibbles to get the same results, hence the use of the permuteNibbleSkinny function

	// uint rMax = 11;
	// auto BCD = genDataSkinny64(rMax);
	// vector<uint8_t>  input(hexToVec("FFFFFFFEFFFFFFFF"));
	// vector<uint8_t> output(hexToVec("1000000000000000"));
	// vector<vector<uint8_t>> keyval({
	// 						hexToVec(permuteNibbleSkinny("0000000800000000")),
	// 						hexToVec(permuteNibbleSkinny("0000000000000000")),
	// 						hexToVec(permuteNibbleSkinny("0000000000008004")),
	// 						hexToVec(permuteNibbleSkinny("0000000480040000")),
	// 						hexToVec(permuteNibbleSkinny("8084800004080400")),
	// 						hexToVec(permuteNibbleSkinny("84C04C4088044028")),
	// 						hexToVec(permuteNibbleSkinny("0400C880C0040000")),
	// 						hexToVec(permuteNibbleSkinny("0000C00000040000")),
	// 						hexToVec(permuteNibbleSkinny("0000400000000000")),
	// 						hexToVec(permuteNibbleSkinny("0000000000400000"))});
	// computeNumberTrails(BCD,input,output,keyval); //Should be 5211
	
	
}
