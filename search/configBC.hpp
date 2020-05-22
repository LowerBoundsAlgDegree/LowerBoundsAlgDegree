#ifndef H_CONGIFBC
#define H_CONGIFBC
#include <vector>
#include <cstdint>
#include <fstream>

#include "BCData.hpp"

/*
	Contains the function used to create the BCData object for different block ciphers
	As well as the function to interface this code with the SageMath code to generate the MILP models

	Note that the MixColumns Matrix has to be provided according to the linAsSbox parameter
	If linAsSbox = true, the MC operation will be modelized using L-boxes
	In this case, the matrix needs to be the matrix of the L-box
	E.g. for Skinny, a 4x4 matrix 
	1011
	1000
	0110
	1010
	This is the better way if applicable (the matrix of MC needs to be binary in GF(2^n))

	If linAsSbox = false, MC is modelized using copy+XOR
	In that case, the matrix needs to be the full matrix over GF2
	E.g for Skinny, it would need to be 
	1000000010001000
	0100000001000100
	0010000000100010
	0001000000010001
	1000000000000000
	0100000000000000
	0010000000000000
	0001000000000000
	0000100010000000
	0000010001000000
	0000001000100000
	0000000100010000
	1000000010000000
	0100000001000000
	0010000000100000
	0001000000010000
*/

class BCData;


BCData genDataSkinny64(uint const rMax);
//Generate the data required for Skinny64 over #rMax rounds

BCData genDataPRESENT(uint const rMax, bool const dedicatedPRESENTlastLayer = false);
/*
	Generate the data required for PRESENT over #rMax rounds
	The #dedicatedPRESENTlastLayer allows to slightly change the last sbox layer
	If we denote by (y0,y1,y2,y3) = S(x0,x1,x2,x3) the output of the original Sbox, the last sbox layer use (y0,y1+y3,y2,y3) instead if #dedicatedPRESENTlastLayer = true
	This allows to prove the min-degree much more easily, while maintaining correctness (since it essentially applies a invertible linear mapping on the output of the last sbox layer)
*/


BCData genDataGIFT(uint const rMax);
//Generate the data required for Craft over #rMax rounds

void generateSSBModel(std::string const & name,
					  std::vector<uint> const & S,
					  std::vector<std::vector<uint8_t>> const & M,
					  uint const sboxSize,
					  uint const SSBSize,
					  bool linAsSbox,
					  bool keyAfterMC,
					  bool dedicatedPRESENTlastLayer = false);
/*
	Generate a model for a Super sbox using genSSBModel.sage
	- #name is the name of the model, the model file will be called #name_SSB.mps
	- #S is the Sbox
	- #M is the MixColumn matrix (careful, the format of this matrix must fit with the linAsSbox variable)
	- #sboxSize is the number of bits of the Sbox
	- #SSBSize is the number of bits of the SSB
	- #linAsSbox is true is MC is modelized as parallel Lboxes, false for COPY+XOR modelization
	  Lboxes modelization is better, but requires the matrix to be binary over GF(2^n) (with a few exceptions like PRINCE, but it's not handled for now)
	- #keyAfterMC gives the form of the SSB:
		- true if the SSB is S -> MC -> ARK -> S
		- false if the SSB is S -> ARK -> MC -> S
*/

void generateBaseModel(std::string const & name,
					   uint const rMax,
					   std::vector<uint> const & S,
					   std::vector<uint> const & P,
					   std::vector<std::vector<uint8_t>> const & M,
					   bool linAsSbox,
					   bool keyAfterMC,
					   bool startAfterSSB,
					   bool noLastMC,
					   bool dedicatedPRESENTlastLayer = false);
/*
	Generate a model for the block cipher using aes_like.sage
	- #name is the name of the model, the model file will be called #name.mps
	- #rMax is the number of rounds
	- #S is the Sbox
	- #P is the permutation, given on WORD level
	- #M is the MixColumn matrix (careful, the format of this matrix must fit with the linAsSbox variable)
	- #linAsSbox is true is MC is modelized as parallel Lboxes, false for COPY+XOR modelization
	  Lboxes modelization is better, but requires the matrix to be binary over GF(2^n) (with a few exceptions like PRINCE, but it's not handled for now)
	- #keyAfterMC gives the form of the round function:
		- true if the round function is S -> P -> MC -> ARK
		- false if the round function is S -> P -> ARK -> MC
	- #startAfterSSB indicates whether the model should start after the Sbox layer
	  if so, the model adds another half round at the beginning
	- #noLastMC indicates if the lqst layer should contains the linear layer
	- #dedicatedPRESENTlastLayer indicates whether we should tweak the last sbox layer for Present (see the comment in genDataPRESENT)
*/

std::string exec(const char* cmd);
//Exec cmd and grab the stdout

#endif