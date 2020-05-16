#include "configBC.hpp"

#define DISPLAY_GEN_OUTPUT false

using namespace std;
typedef unsigned int uint;

BCData genDataMidori64(uint const rMax){
	//Generate the data required for Midori64 over #rMax rounds

	string name = "Midori";
	uint blockSize = 64;
	uint sboxSize = 4;
	uint SSBSize = 16;
	bool keyAfterMC = true;
	bool linAsSbox = true;
	vector<uint> S({0xc,0xa,0xd,0x3,0xe,0xb,0xf,0x7,0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6});
	vector<uint> P({0,7,14,9,5,2,11,12,15,8,1,6,10,13,4,3});
	vector<vector<uint8_t>> M_GF4({{0,1,1,1},{1,0,1,1},{1,1,0,1},{1,1,1,0}});

	//Generate the SSB model
	generateSSBModel(name+"_SSB",S,M_GF4,sboxSize,SSBSize,linAsSbox,keyAfterMC);

	return BCData(name,rMax,blockSize,sboxSize,SSBSize,linAsSbox,keyAfterMC,S,P,M_GF4);
}

BCData genDataSkinny64(uint const rMax){
	//Generate the data required for Skinny64 over #rMax rounds

	string name = "Skinny64";
	uint blockSize = 64;
	uint sboxSize = 4;
	uint SSBSize = 16;
	bool keyAfterMC = false;
	bool linAsSbox = true;
	vector<uint> S({12,6,9,0,1,10,2,11,3,8,5,13,4,14,7,15});
	vector<uint> P({0,5,10,15,4,9,14,3,8,13,2,7,12,1,6,11});
	vector<vector<uint8_t>> M_GF4({{1,0,1,1},{1,0,0,0},{0,1,1,0},{1,0,1,0}});

	//Generate the SSB model
	generateSSBModel(name+"_SSB",S,M_GF4,sboxSize,SSBSize,linAsSbox,keyAfterMC);

	return BCData(name,rMax,blockSize,sboxSize,SSBSize,linAsSbox,keyAfterMC,S,P,M_GF4);
}

BCData genDataCRAFT(uint const rMax){
	//Generate the data required for Craft over #rMax rounds

	string name = "CRAFT";
	uint blockSize = 64;
	uint sboxSize = 4;
	uint SSBSize = 16;
	bool keyAfterMC = false;
	bool linAsSbox = true;
	vector<uint> S({0xc,0xa,0xd,0x3,0xe,0xb,0xf,0x7,0x8,0x9,0x1,0x5,0x0,0x2,0x4,0x6});
	vector<uint> P({15,10,9,4,3,6,5,8,7,2,1,12,11,14,13,0});
	vector<vector<uint8_t>> M_GF4({{1,0,1,1},{0,1,0,1},{0,0,1,0},{0,0,0,1}});

	//Generate the SSB model
	generateSSBModel(name+"_SSB",S,M_GF4,sboxSize,SSBSize,linAsSbox,keyAfterMC);

	return BCData(name,rMax,blockSize,sboxSize,SSBSize,linAsSbox,keyAfterMC,S,P,M_GF4);
}

BCData genDataPRESENT(uint const rMax, bool const dedicatedPRESENTlastLayer){
	/*
		Generate the data required for PRESENT over #rMax rounds
		The #dedicatedPRESENTlastLayer allows to slightly change the last sbox layer
		If we denote by (y0,y1,y2,y3) = S(x0,x1,x2,x3) the output of the original Sbox, the last sbox layer use (y0,y1+y3,y2,y3) instead if #dedicatedPRESENTlastLayer = true
		This allows to prove the min-degree much more easily, while maintaining correctness (since it essentially applies a invertible linear mapping on the output of the last sbox layer)
	*/

	string name = "PRESENT";
	uint blockSize = 64;
	uint sboxSize = 4;
	uint SSBSize = 16;
	bool keyAfterMC = true;
	bool linAsSbox = false;
	vector<uint> S({0xc,0x5,0x6,0xb,0x9,0x0,0xa,0xd,0x3,0xe,0xf,0x8,0x4,0x7,0x1,0x2});
	vector<uint> P({0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15});

	//Original bit permutation
	vector<uint> bp(64,0);
	for(uint i = 0; i < 63; i++)
		bp[i] = (16*i)%63;
	bp[63] = 63;
	
	//SSB Permutation (= MC)
	vector<uint> ssbp(16, 0);
	for(uint i = 0; i < 16; i++)
		ssbp[i] = bp[i]/4 + bp[i]%4;

	//Create the MC matrix
	vector<vector<uint8_t>> M_GF2(16, vector<uint8_t>(16,0));
	for(uint i = 0; i < 16; i++)
		M_GF2[ssbp[i]][i] = 1;

	//Generate the SSB model
	generateSSBModel(name+"_SSB",S,M_GF2,sboxSize,SSBSize,linAsSbox,keyAfterMC);
	if(dedicatedPRESENTlastLayer)
		generateSSBModel(name+"_SSB_dedicated",S,M_GF2,sboxSize,SSBSize,linAsSbox,keyAfterMC,true);

	return BCData(name,rMax,blockSize,sboxSize,SSBSize,linAsSbox,keyAfterMC,S,P,M_GF2,dedicatedPRESENTlastLayer);
}

BCData genDataGIFT(uint const rMax){
	//Generate the data required for Craft over #rMax rounds

	string name = "GIFT";
	uint blockSize = 64;
	uint sboxSize = 4;
	uint SSBSize = 16;
	bool keyAfterMC = true;
	bool linAsSbox = false;
	vector<uint> S({0x1,0xa,0x4,0xc,0x6,0xf,0x3,0x9,0x2,0xd,0xb,0x7,0x5,0x0,0x8,0xe});
	vector<uint> P({0,4,8,12,1,5,9,13,2,6,10,14,3,7,11,15});

	//SSB Permutation (= MC)
	vector<uint> ssbp({0,5,10,15,12,1,6,11,8,13,2,7,4,9,14,3});

	//Create the MC matrix
	vector<vector<uint8_t>> M_GF2(16, vector<uint8_t>(16,0));
	for(uint i = 0; i < 16; i++)
		M_GF2[ssbp[i]][i] = 1;

	//Generate the SSB model
	generateSSBModel(name+"_SSB",S,M_GF2,sboxSize,SSBSize,linAsSbox,keyAfterMC);

	return BCData(name,rMax,blockSize,sboxSize,SSBSize,linAsSbox,keyAfterMC,S,P,M_GF2);
}

void generateSSBModel(string const & name,
					  vector<uint> const & S,
					  vector<vector<uint8_t>> const & M,
					  uint const sboxSize,
					  uint const SSBSize,
					  bool linAsSbox,
					  bool keyAfterMC,
					  bool dedicatedPRESENTlastLayer){
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
	
	ofstream paramfile;
	paramfile.open("SSBModel_parameters.sage", std::ofstream::out | std::ofstream::trunc);

	//Name
	paramfile << "name = \"" << name << "\"" << endl;
	//Sbox
	paramfile << "S = [";
	for(uint i = 0; i < S.size(); i++)
		paramfile << S[i] << ",";
	paramfile << "]" << endl;
	//Matrix
	paramfile << "M = matrix(GF(2),[" << endl;
	for(uint i = 0; i < M.size(); i++){
		auto const & M_i = M[i];
		paramfile << "    [";
		for(uint j = 0; j < M_i.size(); j++)
			paramfile << int(M_i[j]) << ",";
		paramfile << "]," << endl;
	}
	paramfile << "    ])" << endl;
	//sboxSize
	paramfile << "sboxSize = " << sboxSize << endl;
	//blockSize
	paramfile << "blockSize = " << SSBSize << endl;
	//linAsSbox
	paramfile << "linAsSbox = ";
	if(linAsSbox) paramfile << "True" << endl;
	else paramfile << "False" << endl;
	//keyAfterMC
	paramfile << "keyAfterMC = ";
	if(keyAfterMC) paramfile << "True" << endl;
	else paramfile << "False" << endl;
	//dedicatedPRESENTlastLayer
	paramfile << "dedicatedPRESENTlastLayer = ";
	if(dedicatedPRESENTlastLayer) paramfile << "True" << endl;
	else paramfile << "False" << endl;

	paramfile.close();

	auto sysstdout = exec("sage genSSBModel.sage");
	if(DISPLAY_GEN_OUTPUT)
		cout << sysstdout << endl;
}

void generateBaseModel(string const & name,
					   uint const rMax,
					   vector<uint> const & S,
					   vector<uint> const & P,
					   vector<vector<uint8_t>> const & M,
					   bool linAsSbox,
					   bool keyAfterMC,
					   bool startAfterSSB,
					   bool noLastMC,
					   bool dedicatedPRESENTlastLayer){
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
	- #noLastMC indicates if the last layer should contains the linear layer-
	- #dedicatedPRESENTlastLayer indicates whether we should tweak the last sbox layer for Present (see the comment in genDataPRESENT)
*/
	
	ofstream paramfile;
	paramfile.open("aes_like_parameters.sage", std::ofstream::out | std::ofstream::trunc);

	//Name
	paramfile << "name = \"" << name << "\"" << endl;
	//rMax
	paramfile << "rMax = " << rMax << endl;
	//Sbox
	paramfile << "S = [";
	for(uint i = 0; i < S.size(); i++)
		paramfile << S[i] << ",";
	paramfile << "]" << endl;
	//Permutation
	paramfile << "P = [";
	for(uint i = 0; i < P.size(); i++)
		paramfile << P[i] << ",";
	paramfile << "]" << endl;
	//Matrix
	paramfile << "M = matrix(GF(2),[" << endl;
	for(uint i = 0; i < M.size(); i++){
		auto const & M_i = M[i];
		paramfile << "    [";
		for(uint j = 0; j < M_i.size(); j++)
			paramfile << int(M_i[j]) << ",";
		paramfile << "]," << endl;
	}
	paramfile << "    ])" << endl;
	//linAsSbox
	paramfile << "linAsSbox = ";
	if(linAsSbox) paramfile << "True" << endl;
	else paramfile << "False" << endl;
	//keyAfterMC
	paramfile << "keyAfterMC = ";
	if(keyAfterMC) paramfile << "True" << endl;
	else paramfile << "False" << endl;
	//startAfterSSB
	paramfile << "startAfterSSB = ";
	if(startAfterSSB) paramfile << "True" << endl;
	else paramfile << "False" << endl;
	//noLastMC
	paramfile << "noLastMC = ";
	if(noLastMC) paramfile << "True" << endl;
	else paramfile << "False" << endl;
	//dedicatedPRESENTlastLayer
	paramfile << "dedicatedPRESENTlastLayer = ";
	if(dedicatedPRESENTlastLayer) paramfile << "True" << endl;
	else paramfile << "False" << endl;

	paramfile.close();

	auto sysstdout = exec("sage aes_like.sage");
	if(DISPLAY_GEN_OUTPUT)
		cout << sysstdout << endl;
}

string exec(const char* cmd){
//Exec cmd and grab the stdout
    array<char, 128> buffer;
    string result;
    std::shared_ptr<FILE> pipe(popen(cmd, "r"), pclose);
    if (!pipe) throw std::runtime_error("popen() failed!");
    while (!feof(pipe.get())) {
        if (fgets(buffer.data(), 128, pipe.get()) != nullptr)
            result += buffer.data();
    }
    return result;
}