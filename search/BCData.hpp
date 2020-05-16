#ifndef H_BCDATA
#define H_BCDATA
#include <cstdint>
#include <vector>
#include <set>
#include <iostream>
#include <fstream>
#include <string>
#include <utility>
#include <cmath>
#include <algorithm>
#include <omp.h>
#include <bitset>
#include <random>
#include <tuple>
#include <map>

#include <boost/dynamic_bitset.hpp>
#include <gmpxx.h>

#include "gurobi_c++.h"
#include "configBC.hpp"
#include "customCallback.hpp"
#include "SmallMatrix.hpp"

/*
	Most of the stuff is happening here.
	The BCData class contains informations about the block cipher
	Functions to do the different searches are thus defined as members of this class
*/


/*
	A lot of bitvectors are implemented using a vector<uint8_t> instead of bitset/dynamic_bitset
	While it does use more memory to use vector<uint8_t>, I believe access-time is faster
	bitsets are stil used when the number of vector<uint8_t> could grow too high though
	Overall, it might be less pretty, but it still works nicely and does not really hinder performances, as the main bottleneck is the MILP solving
*/

class BCData{
	//Some data about the block cipher using a specific key #k for each SSB in the middle layer

	public:

		std::string name; 								//name of the block cipher
		uint rMax;										//Number of rounds
		unsigned int blockSize; 						//block size in bits
		unsigned int sboxSize; 							//sbox size in bits
		unsigned int SSBSize; 							//SSB size in bits
		unsigned int nbSSB;								//Number of SSB
		unsigned int nbSbox;							//Number of sboxes
		bool linAsSbox;									//Indicate whether the linear layer should be modelized with Lboxes
		bool keyAfterMC;								//Indicate whether the key is added after or before MC
		std::vector<unsigned int> S; 					//Sbox
		std::vector<unsigned int> P; 					//Permutation layer, given over the words (NOT bits)
		std::vector<unsigned int> invP;					//Inverse of P
		std::vector<unsigned int> invPBitLevel;			//Inverse of P at a bit level
		std::vector<std::vector<uint8_t>> Mmodel;		//MixColumns matrix to generate the model
		GRBEnv gurobiEnv;								//Unique Gurobi environment
		std::set<std::pair<std::vector<uint8_t>, std::vector<std::vector<uint8_t>>>> forbiddenInputKey; //list of input/key that we don't want to use. Useful in the lin indep search
		bool dedicatedPRESENTlastLayer; //Indicate if we are using the tweaked last Sbox layer in Present (see configBC files)

		BCData(std::string const & xname,
			   uint const xrMax,
			   unsigned int const xblockSize,
			   unsigned int const xsboxSize,
			   unsigned int const xSSBSize,
			   bool const xlinAsSbox,
			   bool const xkeyAfterMC,
			   std::vector<unsigned int> const & xS,
			   std::vector<unsigned int> const & xP,
			   std::vector<std::vector<uint8_t>> const & xMmodel,
			   bool const xdedicatedPRESENTlastLayer = false);

		uint64_t countSSBTrail(std::vector<uint8_t> & u,
							   std::vector<uint8_t> & k,
							   std::vector<uint8_t> & v) const;
		uint64_t countSSBTrail_dedicatedPresent(std::vector<uint8_t> & u,
							   std::vector<uint8_t> & k,
							   std::vector<uint8_t> & v) const;
		/*
			Return the number of trails over the SSB with
			#u the input division property
			#k the key division property
			#v the output division property
			For u,k and v, the division property is encoded over an uint32_t
			bit i is obtained with u & (1 << i)
			The number of relevant bits in the uint32_t is given by #this.SSBSize
		*/

		std::pair<uint64_t,bool> countTrailsFullCipher(std::vector<uint8_t> const & input,
							   		   std::vector<uint8_t> const & output, 
							   		   std::vector<std::vector<uint8_t>> const & keyval);
		/*
			Count the number of trails from #input to #output using #keyval for the keys' division property
			Return a pair (nbTrail, check) where :
			- #check indicates if the full number of trails was able to be computed in the given limits (true if full number, false if a limit was reached)
			- nbTrail indicates the number of trails found
			Thus if check=true, nbTrail is the total nuber of trails
		*/

		std::pair<std::vector<std::vector<std::vector<uint8_t>>>, 
				  std::vector<std::vector<boost::dynamic_bitset<>>>>
		searchKeyPattern(std::string const & modelFile,
						 int const nbRounds,
						 std::vector<uint8_t> const & input,
						 std::vector<uint8_t> const & output,
						 bool const startAfterSSB,
						 std::vector<std::vector<boost::dynamic_bitset<>>> const & allTriedKeyPattern = std::vector<std::vector<boost::dynamic_bitset<>>>());
		/*
			Search for a key pattern with an odd number of trails from input to output over nbRounds rounds
			Return a pair (Lsol, L) where Lsol contains a list of possible key patterns and L is the list of all key pattern tested during the search (to avoid trying them later on if needed)
			- #modelFile is the name of the base .mps file for the model (without the .mps extension)
			- #nbRounds is the number of rounds covered by the model
			- #input is the division property at the input
			- #output is the division property at the output
			- #startAfterSSB indicates whether the model start before or after the middle layer SSB
			- #allTriedKeyPattern contains the previously tested key pattern
			  these are removed from the solution pool to find new ones
		*/

		void searchMinDegree();
		/*
			Try to prove the full min-degree for each block (Sbox or Super Sbox depending on where the key is added)
		*/

		void searchMinDegree_allInput();
		/*
			Try to prove the "Full monomial property" for each block (Sbox or Super Sbox depending on where the key is added)
			Essentially, it tries to prove full min-degree using each possible input monomial of degree #BCD.blockSize-1
		*/

		std::tuple<std::vector<uint8_t>, std::vector<uint8_t>, std::vector<std::vector<uint8_t>>>
		improvedDynamicSearch(uint const indexOutput,
							  int const indexInput = -1);
		/*
			Try to prove full degree for the ouptut bit of index #indexOutput
			#indexInput is optional and allows to also fix the input bit set to 0 (thus defining the input monomial of degree 63 as the one without x[#indexInput])
			By default, #indexInput = -1, meaning that any input bit can be used
			If #indexInput > -1, then it fix the input as described above
			Returns a tuple (input,output,keyval) where :
			- #input is the input parity set
			- #output is the output parity set
			- #keyval contains the parity setof each key
			If the search fails for any reason (essentially, infeasible model), then the tuple contains empty vectors
		*/

		std::vector<uint> getLowerBounds(std::vector<uint8_t> const & input,
							std::vector<uint8_t> const & output,
							std::vector<std::vector<uint8_t>> const & keyVal);
		/*
			Compute lower bounds for 1 to #this.rMax-1 rounds using #input, #output and #keyval as helpers
			As described in the paper, the idea is to use the input/key/output obtained from the full degree proof to get lower bounds on a lower number of rounds
			Thus #input,#output,#keyval should describe the parity set needed to get an odd number of trails over #this.rMax rounds from #input to #output using #keyval for the keys
			Especially, #keyval should contain #this.rMax-1 vectors

			Note that technically, #input does not necessarily need to be of weight #this.blockSize-1 for this function to work, it only needs the total number of trails over #this.rMax rounds to be odd.
			Return a vector LB of length #this.rMax such that LB[r] is a lower bound on the degree for #r rounds
		*/

		std::vector<uint> getUpperBoundAllOutput();
		/*
			Get an upper bounds on the degree of each output bit for #this.rMax rounds
			Returns a vector UB of size #this.blockSize such that UB[i] is the upper bound for bit #i
		*/

		void printParitySetHex(std::vector<uint8_t> const & x) const;
		/*
			Hex printer for more compact prints
			LSB of x is x[0]
			An hex character thus represent (x[i],...,x[i+3]) with x[i] LSB
			e.g.
			x = {1,0,0,0} is printed as 1
			x = {0,1,0,0,1,1,0,1} is printed as 2B
		*/
};

uint selectNextOutputBit(SmallMatrix const & M);
/*
	Help selecting the next output bit to consider
	If the matrix is empty, return a random column
	Else, return a random column that does not have a unit vector
*/

uint selectNextOutputBit_dedicatedPresent(SmallMatrix const & M);

void printBits(uint64_t const x, uint nbBits);
//print the first #nbBits of x in cout

void printVec(std::vector<uint8_t> const & x);
//print the elements of x, wihtout separators

std::string hex_char_to_bin(char c);
//hex to bin conversion, LSB on the left in the binary representation

std::vector<uint8_t> hexToVec(std::string const & s);
/*
	From an hex string get the corresponding vector of bits
	LSB of each hex char is put on the left
	e.g. 
	s = "1"  returns {1,0,0,0}
	s = "2B" returns {0,1,0,0,1,1,0,1}
*/

#endif