#include "BCData.hpp"

#define MAX_GRBMIPGapAbs 8						/*Accept non-optimal solutions in searchKeyPattern if the gap between the objective and the best bound on the objective is less than MAX_GRBMIPGapAbs*/
#define MAX_NUMBER_SEARCHKEY 10					/*Maximum number of key patterns to search in searchKeyPattern*/
#define TIME_LIMIT_SEARCHKEY 300				/*Time limit (in seconds) in searchKeyPattern*/
#define MAX_GRB_THREAD 8						/*Max number of thread to use by Gurobi. Setting it to 0 let Gurobi decide it itself*/
#define ENFORCE_NONZERO_LASTKEY true			/*If true, enforce the last key to be non-zero. NECESSARY for the min-degree search, sometimes improve the speed for the degree search. Most results on the algebraic degree were obtained with this parameter set to false*/

#define WEIGHT_THRESHOLD_STEP1 40				/*Threshold on the wheight of the input of a round to go to the next step in the dynamic search algorithm*/
#define TIME_LIMIT_ONESTEPDYNAMIC 30			/*Time limit (in seconds) in improvedDynamicSearch to find one iteration in the first step*/
#define DYNAMIC_STEP1_MIPFOCUS 1				/*Hint for Gurobi to focus more on finding solution than proving optimality. Seems to improve the speed, see the Gurobi doc for more details*/

#define MAX_NUMBER_SOLCOUNT_FULLCIPHER 100000	/*Solution limit when counting the total number of trails over the whole cipher in countTrailsFullCipher*/
#define TIME_LIMIT_COUNTTRAILS_FULLCIPER 240	/*Time limit (in seconds) when counting the total number of trails over the whole cipher in countTrailsFullCipher. Set to GRB_INFINITY for unlimited time*/

/*The values of the last two limits are the one used to prove the min-degree.
For the algebraic degree for PRESENT, higher limits were used but were not saved.
It should be something like 1 000 000 max solutions and unlimited time
*/

/*
	Macros are used in the code to define different limits on the search (e.g. time limits)
	It's a bit dirty but easier to change and than handling a ton of parameters everywhere
	Our results were obtained with the current values
*/

using namespace std;
using namespace boost;

typedef unsigned int uint;
typedef pair<uint32_t,uint32_t> pairuint32;
typedef tuple<uint32_t,uint32_t,uint64_t> tup3uint;

std::random_device rd;  //Will be used to obtain a seed for the random number engine
std::mt19937 prng(rd()); //Standard mersenne_twister_engine seeded with rd()
std::uniform_int_distribution<uint> randintSeedGurobi(0, 1999999999);
std::uniform_int_distribution<uint> randintColumn;

//Some data about the block cipher using a specific key #k
BCData::BCData(string const & xname,
			   uint const xrMax,
			   unsigned int const xblockSize,
			   unsigned int const xsboxSize,
			   unsigned int const xSSBSize,
			   bool const xlinAsSbox,
			   bool const xkeyAfterMC,
			   vector<unsigned int> const & xS,
			   vector<unsigned int> const & xP,
			   vector<vector<uint8_t>> const & xMmodel,
			   bool const xdedicatedPRESENTlastLayer) :
			   //Initialization list
			   name(xname),
			   rMax(xrMax),
			   blockSize(xblockSize),
			   sboxSize(xsboxSize),
			   SSBSize(xSSBSize),
			   nbSSB(xblockSize/xSSBSize),
			   nbSbox(xP.size()),
			   linAsSbox(xlinAsSbox),
			   keyAfterMC(xkeyAfterMC),
			   S(xS),
			   P(xP),
			   invP(xP.size()),
			   invPBitLevel(blockSize),
			   Mmodel(xMmodel),
			   gurobiEnv(),
			   forbiddenInputKey(),
			   dedicatedPRESENTlastLayer(xdedicatedPRESENTlastLayer)
{	
	if(rMax < 2)
		cout << "WARNING : rMax < 2, this can break a lot of stuff as it was not designed to handle this edge case" << endl;

	//Compute the inverse of P
	for(uint i = 0; i < P.size(); i++){
		invP[P[i]] = i;
	}

	//Compute the inverse of P at a bit level
	for(uint isbox = 0; isbox < P.size(); isbox++){
		for(uint j = 0; j < sboxSize; j++){
			invPBitLevel[isbox*sboxSize+j] = invP[isbox]*sboxSize + j;
		}
	}

	//Create a single Gurobi environment for the search
	gurobiEnv.set(GRB_IntParam_OutputFlag,0);
}

uint64_t BCData::countSSBTrail(vector<uint8_t> & u,
							   vector<uint8_t> & k,
							   vector<uint8_t> & v) const{
/*
	Return the number of trails over the SSB with
	#u the input division property
	#k the key division property
	#v the output division property
	For u,k and v, the division property is encoded over an uint32_t
	bit i is obtained with u & (1 << i)
	The number of relevant bits in the uint32_t is given by #this.SSBSize
*/

	//Read the model
	GRBModel m(gurobiEnv, name+"_SSB.mps");
	m.set(GRB_IntParam_Threads,1);

	//Get u,k,v variables
	vector<GRBVar> uvar(SSBSize);
	vector<GRBVar> kvar(SSBSize);
	vector<GRBVar> vvar(SSBSize);
	for(uint i = 0; i < SSBSize; i++){
		uvar[i] = m.getVarByName("u"+to_string(i));
		kvar[i] = m.getVarByName("k"+to_string(i));
		vvar[i] = m.getVarByName("v"+to_string(i));
	}

	//Set the value of u,k,v
	for(uint i = 0; i < SSBSize; i++){
		m.addConstr(uvar[i] == u[i]);
		m.addConstr(kvar[i] == k[i]);
		m.addConstr(vvar[i] == v[i]);
	}

	//Setup params to enumerate solutions
	m.set(GRB_IntParam_PoolSearchMode, 2);
	m.set(GRB_IntParam_PoolSolutions, 2000000000); //maximum number of solution allowed by Gurobi

	m.optimize();

	return m.get(GRB_IntAttr_SolCount);
}

pair<uint64_t,bool> BCData::countTrailsFullCipher(vector<uint8_t> const & input,
							   vector<uint8_t> const & output, 
							   vector<vector<uint8_t>> const & keyval){
	/*
		Count the number of trails from #input to #output using #keyval for the keys' division property
		Return a pair (nbTrail, check) where :
		- #check indicates if the full number of trails was able to be computed in the given limits (true if full number, false if a limit was reached)
		- nbTrail indicates the number of trails found
		Thus if check=true, nbTrail is the total nuber of trails
	*/

	//Count the overall number of trails if we remove the middle layer constraints
	generateBaseModel(name,rMax,S,P,Mmodel,linAsSbox,keyAfterMC,false,true,dedicatedPRESENTlastLayer);
	GRBModel m(gurobiEnv, name+".mps");

	int rMax_minus_1 = (int)(rMax) - 1;
	//Read variables
	//input/output variables
	vector<GRBVar> inputVar(blockSize);
	vector<GRBVar> outputVar(blockSize);
	//input vars are x_0_i
	//output vars are x_rMax_i
	for(uint i = 0; i < blockSize; i++){
		inputVar[i] = m.getVarByName("x_0_"+to_string(i));
		outputVar[i] = m.getVarByName("y_"+to_string(rMax_minus_1)+"_"+to_string(i));
	}
	//Key variables
	vector<vector<GRBVar>> allKeyVar(rMax_minus_1, vector<GRBVar>(blockSize));
	for(int r = 0; r < rMax_minus_1; r++){
		auto & allKeyVar_r = allKeyVar[r];
		for(uint i = 0; i < blockSize; i++){
			allKeyVar_r[i] = m.getVarByName("k_"+to_string(r)+"_"+to_string(i));
		}
	}

	// cout << "Input  : "; printVec(input); cout << endl;
	// cout << "Output : "; printVec(output); cout << endl;
	// cout << "Keys   : " << endl;
	// for(auto const & keypattern : keyval){
	// 	printVec(keypattern);
	// 	cout << endl;
	// }
	// cout << endl;

	//Value constraints
	//input/output
	for(uint i = 0; i < blockSize; i++){
		m.addConstr(inputVar[i] == input[i]);
		m.addConstr(outputVar[i] == output[i]);
	}
	//keys
	for(int r = 0; r < rMax_minus_1; r++){
		auto & allKeyVar_r = allKeyVar[r];
		auto const & keyval_r = keyval[r];
		for(uint i = 0; i < blockSize; i++){
			m.addConstr(allKeyVar_r[i] == keyval_r[i]);
		}
	}

	//Count trails
	m.set(GRB_IntParam_PoolSearchMode, 2);
	m.set(GRB_IntParam_PoolSolutions, 2000000000); //maximum number of solution allowed by Gurobi
	auto seedGurobi = randintSeedGurobi(prng);
	m.set(GRB_IntParam_Seed, seedGurobi);
	//Dummy objective
	GRBLinExpr objExpr;
	for(uint i = 0; i < blockSize; i++){
		objExpr += m.getVarByName("x_"+to_string(rMax/2)+"_"+to_string(i));
	}
	m.setObjective(objExpr);

	cout << "Counting the total number of trails" << endl;
	callbackCount cb = callbackCount(MAX_NUMBER_SOLCOUNT_FULLCIPHER);
	m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_COUNTTRAILS_FULLCIPER);
	m.set(GRB_IntParam_Threads,MAX_GRB_THREAD);
	m.setCallback(&cb);	
	m.update();
	m.optimize();

	if(m.get(GRB_IntAttr_Status) == GRB_INFEASIBLE)
		return make_pair(0,true);

	uint totalNumberTrails = m.get(GRB_IntAttr_SolCount);
	bool fullOptimized = true;
	if(!cb.wasAborted && m.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
		cout << endl << totalNumberTrails << " total trails" << endl;
	}
	else if(m.get(GRB_IntAttr_Status) == GRB_TIME_LIMIT){
		cout << endl << "counting solutions aborted, too much time spent" << endl;
		cout << totalNumberTrails << " trails found" << endl;
		fullOptimized = false;
	}
	else{
		cout << endl << "counting solutions aborted, too many solutions or too much time spent" << endl;
		cout << totalNumberTrails << " trails found" << endl;
		fullOptimized = false;
	}

	return make_pair(totalNumberTrails,fullOptimized);
}


pair<vector<vector<vector<uint8_t>>>, vector<vector<dynamic_bitset<>>>>
BCData::searchKeyPattern(string const & modelFile,
						   int const nbRounds,
						   vector<uint8_t> const & input,
						   vector<uint8_t> const & output,
						   bool const startAfterSSB,
						   vector<vector<dynamic_bitset<>>> const & allTriedKeyPattern){
/*
	Search for a key pattern with an odd number of trails from input to output over nbRounds rounds
	Return a pair (kp, L) where kp is the key pattern and L is the list of all key pattern tested during the search (to avoid trying them later on if needed)
	- #modelFile is the name of the base .mps file for the model (without the .mps extension)
	- #nbRounds is the number of rounds covered by the model
	- #input is the division property at the input
	- #output is the division property at the output
	- #startAfterSSB indicates whether the model start before or after the middle layer SSB 
	- #allTriedKeyPattern contains the previously tested key pattern (optional)
	  these are removed from the solution pool to find new ones
*/
	//Generate the model
	if(startAfterSSB)
		generateBaseModel(modelFile,nbRounds,S,P,Mmodel,linAsSbox,keyAfterMC,startAfterSSB,true,dedicatedPRESENTlastLayer);
	else
		generateBaseModel(modelFile,nbRounds,S,P,Mmodel,linAsSbox,keyAfterMC,startAfterSSB,false,dedicatedPRESENTlastLayer);

	//Read the base model
	GRBModel m(gurobiEnv, modelFile+".mps");
	m.set(GRB_IntParam_LazyConstraints, 1);

	uint offset = 0;
	if(startAfterSSB) offset = 1;

	//Read variables
	//input/output variables
	vector<GRBVar> inputVar(blockSize);
	vector<GRBVar> outputVar(blockSize);
	if(!startAfterSSB){
		//input vars are x_0_i
		//output vars are x_nbRounds_i
		for(uint i = 0; i < blockSize; i++){
			inputVar[i] = m.getVarByName("x_0_"+to_string(i));
			outputVar[i] = m.getVarByName("x_"+to_string(nbRounds)+"_"+to_string(i));
		}
	}
	else{
		//input vars are y_0_i
		//output vars are x_(nbRounds+1)_i
		for(uint i = 0; i < blockSize; i++){
			inputVar[i] = m.getVarByName("y_0_"+to_string(i));
			outputVar[i] = m.getVarByName("y_"+to_string(nbRounds)+"_"+to_string(i));
		}
	}
	//All key variables (one key is not involved in an SSB)
	uint const nbroundspoffset = nbRounds;
	vector<vector<GRBVar>> allKeyVar(nbroundspoffset, vector<GRBVar>(blockSize));
	for(uint r = 0; r < nbroundspoffset; r++){
		auto & allKeyVar_r = allKeyVar[r];
		for(uint i = 0; i < blockSize; i++){
			allKeyVar_r[i] = m.getVarByName("k_"+to_string(r)+"_"+to_string(i));
		}
	}


	//SSB variables
	//index [r][i][j], r-th rounds, i-th SSB, j-th variable
	vector<vector<vector<GRBVar>>> inSSBVar(nbRounds-1, 
											vector<vector<GRBVar>>(nbSSB,
											vector<GRBVar>(SSBSize)));
	vector<vector<vector<GRBVar>>> outSSBVar(nbRounds-1, 
											 vector<vector<GRBVar>>(nbSSB,
											 vector<GRBVar>(SSBSize)));
	vector<vector<vector<GRBVar>>> keySSBVar(nbRounds-1, 
											 vector<vector<GRBVar>>(nbSSB,
											 vector<GRBVar>(SSBSize)));

	for(int r = 0; r < nbRounds-1; r++){
		auto & inSSBVar_r = inSSBVar[r];
		auto & outSSBVar_r = outSSBVar[r];
		auto & keySSBVar_r = keySSBVar[r];
		uint rpoffset = r + offset;
		for(uint i = 0; i < nbSSB; i++){
			auto & inSSBVar_r_i = inSSBVar_r[i];
			auto & outSSBVar_r_i = outSSBVar_r[i];
			auto & keySSBVar_r_i = keySSBVar_r[i];
			for(uint j = 0; j < SSBSize; j++){
				uint bit = SSBSize*i + j;
				//input of the SSB
				inSSBVar_r_i[j] = m.getVarByName("x_"+to_string(rpoffset)+"_"+to_string(invPBitLevel[bit]));
				outSSBVar_r_i[j] = m.getVarByName("y_"+to_string(rpoffset+1)+"_"+to_string(bit));
				keySSBVar_r_i[j] = m.getVarByName("k_"+to_string(rpoffset)+"_"+to_string(bit));
			}
		}
	}

	//Add input/output constraints
	for(uint i = 0; i < blockSize; i++){
		m.addConstr(inputVar[i] == input[i]);
		m.addConstr(outputVar[i] == output[i]);
	}

	//Remove known key pattern
	uint const allKeyVar_size = allKeyVar.size();
	for(auto const & keypattern : allTriedKeyPattern){
		GRBLinExpr cutExpr(0);
		for(uint r = 0; r < allKeyVar_size; r++){
			auto & keypattern_r = keypattern[r];
			auto & allKeyVar_r = allKeyVar[r];
			for(uint i = 0; i < blockSize; i++){
				if(keypattern_r[i] == 0) cutExpr += allKeyVar_r[i];
				else cutExpr += (1 - allKeyVar_r[i]);
			}
		}
		m.addConstr(cutExpr >= 1);
	}

	//If needed, enforce the last key to be non zero
	if(ENFORCE_NONZERO_LASTKEY && startAfterSSB){
		GRBLinExpr sumLastKey(0);
		for(uint i = 0; i < blockSize; i++)
			sumLastKey += m.getVarByName("k_"+to_string(nbRounds-1)+"_"+to_string(i));
		m.addConstr(sumLastKey >= 1);
	}



	//Dummy objective
	GRBLinExpr objExpr(0);
	for(uint r = 0; r < allKeyVar_size; r++){
		auto & allKeyVar_r = allKeyVar[r];
		for(uint i = 0; i < blockSize; i++){
			objExpr += allKeyVar_r[i];
		}
	}
	m.setObjective(objExpr, GRB_MAXIMIZE);

	//setup params for enumerating solutions
	m.set(GRB_IntParam_PoolSearchMode, 2);
	m.set(GRB_IntParam_PoolSolutions, MAX_NUMBER_SEARCHKEY);
	m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_SEARCHKEY);
	m.set(GRB_DoubleParam_MIPGapAbs,MAX_GRBMIPGapAbs);
	auto seedGurobi = randintSeedGurobi(prng);
	m.set(GRB_IntParam_Seed, seedGurobi);

	//Search for key patterns so that we have an odd number of trail from input to output
	//Initialize Callback
	callbackSearchKeyPattern cb(*this,nbRounds,inputVar,outputVar,allKeyVar,inSSBVar,outSSBVar,keySSBVar);
	m.setCallback(&cb);
	m.set(GRB_IntParam_Threads,MAX_GRB_THREAD);
	
	m.update();
	m.optimize();
	
	uint const nbSol = m.get(GRB_IntAttr_SolCount);
	if(nbSol > 0){
		cout << "Number of optimal solutions : " << nbSol << endl;
		//We found a key pattern
		vector<vector<vector<uint8_t>>> Allkeyval(nbSol, vector<vector<uint8_t>>(nbroundspoffset, vector<uint8_t>(blockSize)));
		for(uint indexSol = 0; indexSol < nbSol; indexSol++){
			m.set(GRB_IntParam_SolutionNumber, indexSol);
			auto & keyval = Allkeyval[indexSol];
			for(uint r = 0; r < nbroundspoffset; r++){
				auto & keyval_r = keyval[r];
				auto & allKeyVar_r = allKeyVar[r];
				for(uint i = 0; i < blockSize; i++)
					keyval_r[i] = uint8_t(round(allKeyVar_r[i].get(GRB_DoubleAttr_Xn)));
			}
		}
		return make_pair(Allkeyval, cb.allSolutions);
	}
	else{
		//No key pattern found
		return make_pair(vector<vector<vector<uint8_t>>>(), vector<vector<dynamic_bitset<>>>());
	}
}

void BCData::searchMinDegree(){
	/*
		Try to prove the full min-degree for each block (Sbox or Super Sbox depending on where the key is added)
	*/

	if(!ENFORCE_NONZERO_LASTKEY){
		cout << "WARNING : Asking for min-degree search but the last key is not enforced to be non-zero" << endl;
		cout << "Fix by setting ENFORCE_NONZERO_LASTKEY to true" << endl;
		cout << "Press enter to continue without it, but the results might be wrong" << endl;
		getchar();
	}
	uint lenBlock; //Size of one block diagonal block in the matrix
	uint nbBlock;  //Number of diagonal blocks
	if(keyAfterMC){
		lenBlock = sboxSize;
		nbBlock = nbSbox;
	}
	else{
		lenBlock = SSBSize;
		nbBlock = nbSSB;
	}
	randintColumn = std::uniform_int_distribution<uint>(0,lenBlock-1);

	for(uint iblock = 0; iblock < nbBlock; iblock++){

		forbiddenInputKey = set<pair<vector<uint8_t>, vector<vector<uint8_t>>>>();
		SmallMatrix M(0,lenBlock);
		map<uint,uint> listPivots;

		//Keep track of the partial row we found to use them later on
		map<vector<uint8_t>, pair<vector<uint8_t>, vector<vector<uint8_t>>>> knownPartialRows;
		//knownPartialRows[r] = {i,k}, where #r is a partial row, #i and #k the corresponding input/key

		vector<tuple<vector<uint8_t>, vector<vector<uint8_t>>, uint64_t>> goodInputKeyRow;
		while(M.nrows < lenBlock){
			//Get an output bit to get an input/key
			uint index = selectNextOutputBit(M);
			uint indexOutput = iblock*lenBlock + index;

			//Get an input/key
			auto const iok = improvedDynamicSearch(indexOutput);
			auto const & input = get<0>(iok);
			auto const & keyval = get<2> (iok);

			if(input.size() > 0){

				forbiddenInputKey.emplace(make_pair(input, keyval));

				cout << "Found an input and key" << endl;

				//We know that for this output the number of trails is odd
				uint64_t newRow = uint64_t(1) << index;
				vector<bool> countedBit(lenBlock,false);
				countedBit[index] = true;

				//Count the number of trails for the remaining output bits
				//Don't count if we already have a unit vector for this bit
				for(uint ioutput = 0; ioutput < lenBlock; ioutput++){
					if(ioutput != index){
						//Check if we already have a unit vector for this bit
						bool hasUnitVector = false;
						for(uint i = 0; i < M.nrows; i++){
							auto row = M.get(i);
							if(__builtin_popcountll(row) == 1 && (row & (uint64_t(1) << ioutput))){
								//row is a unit vector and has 1 in the correct column
								hasUnitVector = true;
								break;
							}
						}

						if(!hasUnitVector){

							vector<uint8_t> output(blockSize, 0);
							output[iblock*lenBlock + ioutput] = 1;

							cout << "Input  : "; printParitySetHex(input); cout << endl;
							cout << "Output : "; printParitySetHex(output); cout << endl;
							cout << "Keys   : " << endl;
							for(auto const & keypattern : keyval){
								printParitySetHex(keypattern);
								cout << endl;
							}
							cout << "Counting trails..." << endl;

							auto nbTrail_fullOpti = countTrailsFullCipher(input,output,keyval);
							auto nbTrail = nbTrail_fullOpti.first;
							cout << nbTrail << " trails for output " << ioutput << endl;
							if(nbTrail_fullOpti.second){
								if(nbTrail % 2 == 1) newRow |= (uint64_t(1) << ioutput);
								countedBit[ioutput] = true;
							}
							else{
								cout << "could not compute all trails" << endl;
							}
						}
					}
				}
				bool countedAll = true;
				for(auto const & countedBiti : countedBit) countedAll &= countedBiti;

				if(!countedAll){
					cout << "Partial row : ";
					for(uint i = 0; i < lenBlock; i++){
						if(countedBit[i]){
							if(newRow & (uint64_t(1) << i)) cout <<  1;
							else cout << 0;
						}
						else cout << "*";
					}
					cout << endl;
					//We have a partial row, check if we can reduce it
					//Get the partial row pattern
					vector<uint8_t> partialRow(lenBlock,0);
					for(uint i = 0; i < lenBlock; i++){
						if(countedBit[i]){
							if(newRow & (uint64_t(1) << i)) partialRow[i] = 1;
							else partialRow[i] = 0;
						}
						else partialRow[i] = 2;
					}

					//Try to reduce it with unit vectors from the matrix
					for(uint i = 0; i < lenBlock; i++){
						if(partialRow[i] == 2){
							//Check if the i-th unit vector belongs to the matrix
							for(uint j = 0; j < M.nrows; j++){
								uint64_t row = M.get(j);
								if(__builtin_popcountll(row) == 1 && (row & (uint64_t(1) << i))){
									//i-th unit vector, reduce the partial row
									partialRow[i] = 0;
								}
							}
						}
					}

					cout << "Reduced Partial row : ";
					for(uint i = 0; i < lenBlock; i++){
						if(partialRow[i] != 2){
							if(newRow & (uint64_t(1) << i)) cout <<  1;
							else cout << 0;
						}
						else cout << "*";
					}
					cout << endl;

					//Check if the partial row is actually still partial after reduction
					countedAll = true;
					for(uint i = 0; i < lenBlock && countedAll; i++){
						if(partialRow[i] == 2) countedAll = false;
					}

					//If countedAll = true at this point, then it will behave like we find an actual full row
					//Otherwise it will go for a partial row
					if(!countedAll) knownPartialRows[partialRow] = make_pair(input,keyval);
				}

				if(countedAll){
					//we might need to update the partial rows so we need a loop
					map<uint64_t, pair<vector<uint8_t>, vector<vector<uint8_t>>>> newFullRows;
					newFullRows[newRow] = make_pair(input,keyval);
					while(newFullRows.size() > 0){

						//Pick one full row, update everything with it
						auto newRow_ik_it = newFullRows.begin();
						auto newRow_ik = *(newRow_ik_it);
						newFullRows.erase(newRow_ik_it);

						uint64_t currentNewRow = newRow_ik.first;
						auto const & currentInput = newRow_ik.second.first;
						auto const & currentKeyVal = newRow_ik.second.second;

						//Check if the new row is linearly independent from the others
						uint prev_nrows = M.nrows;
						checkLinIndep(M,currentNewRow,listPivots);

						cout << "For block " << iblock << endl;
						cout << "Input  : "; printParitySetHex(currentInput); cout << endl;
						cout << "Keys   : " << endl;
						for(auto const & keypattern : currentKeyVal){
							printParitySetHex(keypattern);
							cout << endl;
						}
						cout << "newRow : "; printBits(currentNewRow, lenBlock); cout << endl;
						cout << "Current Matrix :" << endl;
						M.print();
						cout << endl;

						ofstream outfile;
						outfile.open (name+"_"+to_string(rMax)+"r_linIndep.txt", ios::out | ios::app);

						outfile << "For block " << iblock << endl;
						outfile << "Input  : "; 
						for(auto const & xx : currentInput) outfile << int(xx);
						outfile << endl;

						outfile << "Keys   : " << endl;
						for(auto const & keypattern : currentKeyVal){
							for(auto const & xx : keypattern) outfile << int(xx); 
							outfile << endl;
						}
						outfile << endl;

						outfile << "newRow : "; 
						for(uint i = 0; i < lenBlock; i++){
							if(currentNewRow & (uint64_t(1) << i)) outfile << 1;
							else outfile << 0;
						}
						outfile << endl;

						outfile << "Current Matrix :" << endl;
						for(uint i = 0; i < M.nrows; i++){
							for(uint j = 0; j < M.ncols; j++)
								outfile << M.get(i,j);
							outfile << endl;
						}
						outfile << endl;

						if(M.nrows > prev_nrows){
							//We found a new linearly independent row
							cout << "Found new independent row" << endl;
							goodInputKeyRow.emplace_back(make_tuple(std::move(currentInput), std::move(currentKeyVal), currentNewRow));

							for(uint inewRow = 0; inewRow < M.nrows; inewRow++){
								uint64_t reducedNewRow = M.get(inewRow);
								if(__builtin_popcountll(reducedNewRow) == 1){
									//We found a new unit row, so we might be able to use previously found partial rows by reducing them

									//Get the index of the 1 in the row
									uint indexUnit = 0;
									while(!(reducedNewRow & (uint64_t(1) << indexUnit))) indexUnit++;

									//Update all known Partial rows
									map<vector<uint8_t>, pair<vector<uint8_t>, vector<vector<uint8_t>>>> updatedKnownPartialRows;
									for(auto const & pr_ik : knownPartialRows){
										auto partialRow = pr_ik.first;

										//If the coefficient in the partial row is undetermined and can be reduced, reduce it
										if(partialRow[indexUnit] == 2) partialRow[indexUnit] = 0;

										//Check if now the partial row if fully determined
										bool fullyDetermined = true;
										for(uint i = 0; i < lenBlock && fullyDetermined; i++){
											if(partialRow[i] == 2) fullyDetermined = false;
										}

										if(fullyDetermined){
											//Compute the row
											uint64_t newFullRow = 0;
											for(uint i = 0; i < lenBlock; i++)
												newFullRow |= (uint64_t(partialRow[i]) << i);
											//Add it to the set of known full rows
											newFullRows[newFullRow] = pr_ik.second;
										}
										else{
											updatedKnownPartialRows[partialRow] = pr_ik.second;
										}
									}
									knownPartialRows = std::move(updatedKnownPartialRows);
								}
							}
						}
					}
				}
				else{
					cout << "Found a new partial row" << endl;
					cout << "For block " << iblock << endl;
					cout << "Input  : "; printParitySetHex(input); cout << endl;
					cout << "Keys   : " << endl;
					for(auto const & keypattern : keyval){
						printParitySetHex(keypattern);
						cout << endl;
					}
					cout << "newRow : ";
					for(uint i = 0; i < lenBlock; i++){
						if(countedBit[i]){
							if(newRow & (uint64_t(1) << i)) cout << 1;
							else cout << 0;
						}
						else cout << "*";
					}
					cout << endl;

					cout << "Current Matrix :" << endl;
					M.print();
					cout << endl;

					ofstream outfile;
					outfile.open (name+"_"+to_string(rMax)+"r_linIndep.txt", ios::out | ios::app);

					outfile << "Found a new partial row" << endl;
					outfile << "For block " << iblock << endl;
					outfile << "Input  : "; 
					for(auto const & xx : input) outfile << int(xx);
					outfile << endl;

					outfile << "Keys   : " << endl;
					for(auto const & keypattern : keyval){
						for(auto const & xx : keypattern) outfile << int(xx); 
						outfile << endl;
					}
					outfile << endl;

					outfile << "newRow : "; 
					for(uint i = 0; i < lenBlock; i++){
						if(countedBit[i]){
							if(newRow & (uint64_t(1) << i)) outfile << 1;
							else outfile << 0;
						}
						else outfile << "*";
					}
					outfile << endl;

					outfile << "Current Matrix :" << endl;
					for(uint i = 0; i < M.nrows; i++){
						for(uint j = 0; j < M.ncols; j++)
							outfile << M.get(i,j);
						outfile << endl;
					}
					outfile << endl;
				}
			}
			else{
				//Could not find an input/key
				//For now, just exit
				cout << "Could not find input/key for block " << iblock << endl;
				exit(0);
			}
		}

		//We proved linear independence for block #iblock
		cout << "Block " << iblock << " is of max degree and linearly independent" << endl;
		cout << "Input/key used : " << endl;
		for(auto const & ikr : goodInputKeyRow){
			auto const & input = get<0>(ikr);
			auto const & keyval = get<1>(ikr);
			cout << "Input  : "; printParitySetHex(input); cout << endl;
			cout << "Keys   : " << endl;
			for(auto const & keypattern : keyval){
				printParitySetHex(keypattern);
				cout << endl;
			}
			cout << endl;
		}
		cout << "Final matrix : " << endl;
		for(auto const & ikr : goodInputKeyRow){
			printBits(get<2>(ikr), lenBlock);
			cout << endl;
		}

		ofstream outfile;
  		outfile.open (name+"_"+to_string(rMax)+"r_linIndep.txt", ios::out | ios::app);

  		outfile << "Block " << iblock << " is of max degree and linearly independent" << endl;
		outfile << "Input/key used : " << endl;
		for(auto const & ikr : goodInputKeyRow){
			auto const & input = get<0>(ikr);
			auto const & keyval = get<1>(ikr);

			outfile << "Input  : "; 
			for(auto const & xx : input) outfile << int(xx);
			outfile << endl;

			outfile << "Keys   : " << endl;
			for(auto const & keypattern : keyval){
				for(auto const & xx : keypattern) outfile << int(xx); 
				outfile << endl;
			}
			outfile << endl;
		}
		outfile << "Final matrix : " << endl;
		for(auto const & ikr : goodInputKeyRow){
			uint64_t y = get<2>(ikr);
			for(uint i = 0; i < lenBlock; i++){
				outfile << (y & 1);
				y >>= 1;
			}
			outfile << endl;
		}
		outfile << endl;
		outfile.close();

	}
}

void BCData::searchMinDegree_allInput(){
	/*
		Try to prove the "Full monomial property" for each block (Sbox or Super Sbox depending on where the key is added)
		Essentially, it tries to prove full min-degree using each possible input monomial of degree #BCD.blockSize-1
	*/

	if(!ENFORCE_NONZERO_LASTKEY){
		cout << "WARNING : Asking for linear independance search but the last key is not enforced to be non-zero" << endl;
		cout << "Fix by setting ENFORCE_NONZERO_LASTKEY to true" << endl;
		cout << "Press enter to continue without it, but the results might be wrong" << endl;
		getchar();
	}
	uint lenBlock; //Size of one block diagonal block in the matrix
	uint nbBlock;  //Number of diagonal blocks
	if(keyAfterMC){
		lenBlock = sboxSize;
		nbBlock = nbSbox;
	}
	else{
		lenBlock = SSBSize;
		nbBlock = nbSSB;
	}
	randintColumn = std::uniform_int_distribution<uint>(0,lenBlock-1);

	for(uint iblock = 0; iblock < nbBlock; iblock++){

		for(uint indexInput = 0; indexInput < blockSize; indexInput++){
			forbiddenInputKey = set<pair<vector<uint8_t>, vector<vector<uint8_t>>>>();
			SmallMatrix M(0,lenBlock);
			map<uint,uint> listPivots;

			//Keep track of the partial row we found to use them later on
			map<vector<uint8_t>, pair<vector<uint8_t>, vector<vector<uint8_t>>>> knownPartialRows;
			//knownPartialRows[r] = {i,k}, where #r is a partial row, #i and #k the corresponding input/key

			vector<tuple<vector<uint8_t>, vector<vector<uint8_t>>, uint64_t>> goodInputKeyRow;
			while(M.nrows < lenBlock){
				//Get an output bit to get an input/key
				uint index = selectNextOutputBit(M);
				uint indexOutput = iblock*lenBlock + index;

				//Get an input/key
				auto const iok = improvedDynamicSearch(indexOutput,indexInput);
				auto const & input = get<0>(iok);
				auto const & keyval = get<2> (iok);

				if(input.size() > 0){

					forbiddenInputKey.emplace(make_pair(input, keyval));

					cout << "Found a key" << " for block " << iblock << " and input " << indexInput << endl;

					//We know that for this output the number of trails is odd
					uint64_t newRow = uint64_t(1) << index;
					vector<bool> countedBit(lenBlock,false);
					countedBit[index] = true;

					//Count the number of trails for the remaining output bits
					//Don't count if we already have a unit vector for this bit
					for(uint ioutput = 0; ioutput < lenBlock; ioutput++){
						if(ioutput != index){
							//Check if we already have a unit vector for this bit
							bool hasUnitVector = false;
							for(uint i = 0; i < M.nrows; i++){
								auto row = M.get(i);
								if(__builtin_popcountll(row) == 1 && (row & (uint64_t(1) << ioutput))){
									//row is a unit vector and has 1 in the correct column
									hasUnitVector = true;
									break;
								}
							}

							if(!hasUnitVector){

								vector<uint8_t> output(blockSize, 0);
								output[iblock*lenBlock + ioutput] = 1;

								cout << "Input  : "; printParitySetHex(input); cout << endl;
								cout << "Output : "; printParitySetHex(output); cout << endl;
								cout << "Keys   : " << endl;
								for(auto const & keypattern : keyval){
									printParitySetHex(keypattern);
									cout << endl;
								}
								cout << "Counting trails..." << endl;

								auto nbTrail_fullOpti = countTrailsFullCipher(input,output,keyval);
								auto nbTrail = nbTrail_fullOpti.first;
								cout << nbTrail << " trails for output " << ioutput << endl;
								if(nbTrail_fullOpti.second){
									if(nbTrail % 2 == 1) newRow |= (uint64_t(1) << ioutput);
									countedBit[ioutput] = true;
								}
								else{
									cout << "could not compute all trails" << endl;
								}
							}
						}
					}
					bool countedAll = true;
					for(auto const & countedBiti : countedBit) countedAll &= countedBiti;

					if(!countedAll){
						cout << "Partial row : ";
						for(uint i = 0; i < lenBlock; i++){
							if(countedBit[i]){
								if(newRow & (uint64_t(1) << i)) cout <<  1;
								else cout << 0;
							}
							else cout << "*";
						}
						cout << endl;
						//We have a partial row, check if we can reduce it
						//Get the partial row pattern
						vector<uint8_t> partialRow(lenBlock,0);
						for(uint i = 0; i < lenBlock; i++){
							if(countedBit[i]){
								if(newRow & (uint64_t(1) << i)) partialRow[i] = 1;
								else partialRow[i] = 0;
							}
							else partialRow[i] = 2;
						}

						//Try to reduce it with unit vectors from the matrix
						for(uint i = 0; i < lenBlock; i++){
							if(partialRow[i] == 2){
								//Check if the i-th unit vector belongs to the matrix
								for(uint j = 0; j < M.nrows; j++){
									uint64_t row = M.get(j);
									if(__builtin_popcountll(row) == 1 && (row & (uint64_t(1) << i))){
										//i-th unit vector, reduce the partial row
										partialRow[i] = 0;
									}
								}
							}
						}

						cout << "Reduced Partial row : ";
						for(uint i = 0; i < lenBlock; i++){
							if(partialRow[i] != 2){
								if(newRow & (uint64_t(1) << i)) cout <<  1;
								else cout << 0;
							}
							else cout << "*";
						}
						cout << endl;

						//Check if the partial row is actually still partial after reduction
						countedAll = true;
						for(uint i = 0; i < lenBlock && countedAll; i++){
							if(partialRow[i] == 2) countedAll = false;
						}

						//If countedAll = true at this point, then it will behave like we find an actual full row
						//Otherwise it will go for a partial row
						if(!countedAll) knownPartialRows[partialRow] = make_pair(input,keyval);
					}

					if(countedAll){
						//we might need to update the partial rows so we need a loop
						map<uint64_t, pair<vector<uint8_t>, vector<vector<uint8_t>>>> newFullRows;
						newFullRows[newRow] = make_pair(input,keyval);
						while(newFullRows.size() > 0){

							//Pick one full row, update everything with it
							auto newRow_ik_it = newFullRows.begin();
							auto newRow_ik = *(newRow_ik_it);
							newFullRows.erase(newRow_ik_it);

							uint64_t currentNewRow = newRow_ik.first;
							auto const & currentInput = newRow_ik.second.first;
							auto const & currentKeyVal = newRow_ik.second.second;

							//Check if the new row is linearly independent from the others
							uint prev_nrows = M.nrows;
							checkLinIndep(M,currentNewRow,listPivots);

							cout << "For block " << iblock << " and input " << indexInput << endl;
							cout << "Input  : "; printParitySetHex(currentInput); cout << endl;
							cout << "Keys   : " << endl;
							for(auto const & keypattern : currentKeyVal){
								printParitySetHex(keypattern);
								cout << endl;
							}
							cout << "newRow : "; printBits(currentNewRow, lenBlock); cout << endl;
							cout << "Current Matrix :" << endl;
							M.print();
							cout << endl;

							ofstream outfile;
							outfile.open (name+"_"+to_string(rMax)+"r_linIndep_allInput.txt", ios::out | ios::app);

							outfile << "For block " << iblock << " and input " << indexInput << endl;
							outfile << "Input  : "; 
							for(auto const & xx : currentInput) outfile << int(xx);
							outfile << endl;

							outfile << "Keys   : " << endl;
							for(auto const & keypattern : currentKeyVal){
								for(auto const & xx : keypattern) outfile << int(xx); 
								outfile << endl;
							}
							outfile << endl;

							outfile << "newRow : "; 
							for(uint i = 0; i < lenBlock; i++){
								if(currentNewRow & (uint64_t(1) << i)) outfile << 1;
								else outfile << 0;
							}
							outfile << endl;

							outfile << "Current Matrix :" << endl;
							for(uint i = 0; i < M.nrows; i++){
								for(uint j = 0; j < M.ncols; j++)
									outfile << M.get(i,j);
								outfile << endl;
							}
							outfile << endl;

							if(M.nrows > prev_nrows){
								//We found a new linearly independent row
								cout << "Found new independent row" << endl;
								goodInputKeyRow.emplace_back(make_tuple(std::move(currentInput), std::move(currentKeyVal), currentNewRow));

								for(uint inewRow = 0; inewRow < M.nrows; inewRow++){
									uint64_t reducedNewRow = M.get(inewRow);
									if(__builtin_popcountll(reducedNewRow) == 1){
										//We found a new unit row, so we might be able to use previously found partial rows by reducing them

										//Get the index of the 1 in the row
										uint indexUnit = 0;
										while(!(reducedNewRow & (uint64_t(1) << indexUnit))) indexUnit++;

										//Update all known Partial rows
										map<vector<uint8_t>, pair<vector<uint8_t>, vector<vector<uint8_t>>>> updatedKnownPartialRows;
										for(auto const & pr_ik : knownPartialRows){
											auto partialRow = pr_ik.first;

											//If the coefficient in the partial row is undetermined and can be reduced, reduce it
											if(partialRow[indexUnit] == 2) partialRow[indexUnit] = 0;

											//Check if now the partial row if fully determined
											bool fullyDetermined = true;
											for(uint i = 0; i < lenBlock && fullyDetermined; i++){
												if(partialRow[i] == 2) fullyDetermined = false;
											}

											if(fullyDetermined){
												//Compute the row
												uint64_t newFullRow = 0;
												for(uint i = 0; i < lenBlock; i++)
													newFullRow |= (uint64_t(partialRow[i]) << i);
												//Add it to the set of known full rows
												newFullRows[newFullRow] = pr_ik.second;
											}
											else{
												updatedKnownPartialRows[partialRow] = pr_ik.second;
											}
										}
										knownPartialRows = std::move(updatedKnownPartialRows);
									}
								}
							}
						}
					}
					else{
						cout << "Found a new partial row" << endl;
						cout << "For block " << iblock << " and input " << indexInput << endl;
						cout << "Input  : "; printParitySetHex(input); cout << endl;
						cout << "Keys   : " << endl;
						for(auto const & keypattern : keyval){
							printParitySetHex(keypattern);
							cout << endl;
						}
						cout << "newRow : ";
						for(uint i = 0; i < lenBlock; i++){
							if(countedBit[i]){
								if(newRow & (uint64_t(1) << i)) cout << 1;
								else cout << 0;
							}
							else cout << "*";
						}
						cout << endl;

						cout << "Current Matrix :" << endl;
						M.print();
						cout << endl;

						ofstream outfile;
						outfile.open (name+"_"+to_string(rMax)+"r_linIndep_allInput.txt", ios::out | ios::app);

						outfile << "Found a new partial row" << endl;
						outfile << "For block " << iblock << " and input " << indexInput << endl;
						outfile << "Input  : "; 
						for(auto const & xx : input) outfile << int(xx);
						outfile << endl;

						outfile << "Keys   : " << endl;
						for(auto const & keypattern : keyval){
							for(auto const & xx : keypattern) outfile << int(xx); 
							outfile << endl;
						}
						outfile << endl;

						outfile << "newRow : "; 
						for(uint i = 0; i < lenBlock; i++){
							if(countedBit[i]){
								if(newRow & (uint64_t(1) << i)) outfile << 1;
								else outfile << 0;
							}
							else outfile << "*";
						}
						outfile << endl;

						outfile << "Current Matrix :" << endl;
						for(uint i = 0; i < M.nrows; i++){
							for(uint j = 0; j < M.ncols; j++)
								outfile << M.get(i,j);
							outfile << endl;
						}
						outfile << endl;
					}
				}
				else{
					//Could not find an input/key
					//For now, just exit
					cout << "Could not find input/key for block " << iblock << " with input " << indexInput << endl;
					ofstream outfile;
		  			outfile.open (name+"_"+to_string(rMax)+"r_linIndep_allInput.txt", ios::out | ios::app);
		  			outfile << "Could not find input/key for block " << iblock << " with input " << indexInput << endl;
		  			outfile.close();
					break;
				}
			}

			if(M.nrows == lenBlock){
				//We proved linear independence for block #iblock
				cout << "Block " << iblock << " with  input " << indexInput << " is of max degree and linearly independent" << endl;
				cout << "Input/key used : " << endl;
				for(auto const & ikr : goodInputKeyRow){
					auto const & input = get<0>(ikr);
					auto const & keyval = get<1>(ikr);
					cout << "Input  : "; printParitySetHex(input); cout << endl;
					cout << "Keys   : " << endl;
					for(auto const & keypattern : keyval){
						printParitySetHex(keypattern);
						cout << endl;
					}
					cout << endl;
				}
				cout << "Final matrix : " << endl;
				for(auto const & ikr : goodInputKeyRow){
					printBits(get<2>(ikr), lenBlock);
					cout << endl;
				}

				ofstream outfile;
		  		outfile.open (name+"_"+to_string(rMax)+"r_linIndep_allInput.txt", ios::out | ios::app);

		  		outfile << "Block " << iblock << " with input " << indexInput << " is of max degree and linearly independent" << endl;
				outfile << "Input/key used : " << endl;
				for(auto const & ikr : goodInputKeyRow){
					auto const & input = get<0>(ikr);
					auto const & keyval = get<1>(ikr);

					outfile << "Input  : "; 
					for(auto const & xx : input) outfile << int(xx);
					outfile << endl;

					outfile << "Keys   : " << endl;
					for(auto const & keypattern : keyval){
						for(auto const & xx : keypattern) outfile << int(xx); 
						outfile << endl;
					}
					outfile << endl;
				}
				outfile << "Final matrix : " << endl;
				for(auto const & ikr : goodInputKeyRow){
					uint64_t y = get<2>(ikr);
					for(uint i = 0; i < lenBlock; i++){
						outfile << (y & 1);
						y >>= 1;
					}
					outfile << endl;
				}
				outfile << endl;
				outfile.close();
			}
		}
	}
}

tuple<vector<uint8_t>, vector<uint8_t>, vector<vector<uint8_t>>>
BCData::improvedDynamicSearch(uint const indexOutput,
							  int const indexInput){
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

	vector<uint8_t> output(blockSize, 0);
	output[indexOutput] = 1;

	int rMax_minus_1 = (int)(rMax) - 1;
	int rMax_minus_2 = (int)(rMax) - 2;
	//Generate the model
	generateBaseModel(name+"_ImpDynamic",rMax,S,P,Mmodel,linAsSbox,keyAfterMC,false,true,dedicatedPRESENTlastLayer);

	vector<uint> possibleInput;
	if(indexInput > -1)
		possibleInput = vector<uint>({(uint)(indexInput)});
	else{
		possibleInput = vector<uint>(blockSize);
		for(uint i = 0; i < blockSize; i++) possibleInput[i] = i;
	}
	while(possibleInput.size() > 0){
		//Prepare the vectors to receive the iteratively built values
		vector<vector<uint8_t>> valX(rMax, vector<uint8_t>(blockSize, 0));
		vector<vector<uint8_t>> valK(rMax_minus_1, vector<uint8_t>(blockSize, 0));

		// //If the inputIndex is fixed, use it
		// //Otherwise use a random one
		// uint fixedInput = 0;
		// if(indexInput > -1)
		// 	fixedInput = indexInput;
		// else{
		// 	std::uniform_int_distribution<uint> randintBlockSize(0,blockSize-1);
		// 	fixedInput = randintBlockSize(prng);
		// }

		//Select an input
		std::uniform_int_distribution<uint> randintSelectInput(0,possibleInput.size()-1);
		uint indexSelectedInput = randintSelectInput(prng);
		uint fixedInput = possibleInput[indexSelectedInput];



		uint goodrMiddle = 0;
		bool foundMiddle = false;
		vector<uint8_t> input(blockSize);
		for(int rMiddle = rMax_minus_2; rMiddle >= 0; rMiddle--){
			//Iteratively build the last rounds
			//We are looking to compute X[rMiddle] and K[rMiddle] at this step

			//Read the model
	        GRBModel m(gurobiEnv,name+"_ImpDynamic.mps");
	        m.set(GRB_IntParam_LazyConstraints, 1);

	        if(dedicatedPRESENTlastLayer){
	        	//Get the index of the output bit in the sbox
				uint indexInSbox = indexOutput%sboxSize;
				uint indexActiveSbox = indexOutput/sboxSize;
		        //Force the input of the active sbox in the last round
		        if(indexInSbox == 3){
		        	//Input forced to 1110 = 0x7
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 0)) == 1);
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 1)) == 1);
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 2)) == 1);
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 3)) == 0);
		        }
		        else if(indexInSbox == 2){
		        	//Input forced to 1101 = 0xB
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 0)) == 1);
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 1)) == 1);
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 2)) == 0);
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 3)) == 1);
		        }
		        else if(indexInSbox == 1){
		        	//Input forced to 0xA = 0101 or 0xC = 0011
		        	//These should be the only ones of weight 2 so just add that
		        	GRBLinExpr sumExpr(0);
		        	for(uint i = 0; i < sboxSize; i++)
		        		sumExpr += m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + i));
		        	m.addConstr(sumExpr == 2);
		        }
		        else if(indexInSbox == 0){
		        	//Input forced to 0110 = 0x6
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 0)) == 0);
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 1)) == 1);
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 2)) == 1);
		        	m.addConstr(m.getVarByName("x_"+to_string(rMax_minus_1)+"_"+to_string(sboxSize*indexActiveSbox + 3)) == 0);
		        }
	        }

			//If needed, enforce the last key to be non zero
			if(ENFORCE_NONZERO_LASTKEY){
				GRBLinExpr sumLastKey(0);
				for(uint i = 0; i < blockSize; i++)
					sumLastKey += m.getVarByName("k_"+to_string(rMax_minus_2)+"_"+to_string(i));
				m.addConstr(sumLastKey >= 1);
			}

			//Input of weight blockSize-1
			GRBLinExpr weightInput(0);
			for(uint i = 0; i < blockSize; i++)
				weightInput += m.getVarByName("x_0_"+to_string(i));
			m.addConstr(weightInput == blockSize-1);

			//Fix the input
			m.addConstr(m.getVarByName("x_0_"+to_string(fixedInput)) == 0);


			//Fix the output
			for(uint i = 0; i < blockSize; i++)
				m.addConstr(m.getVarByName("y_"+to_string(rMax_minus_1)+"_"+to_string(i)) == output[i]);

			//Fix the known keys
			for(int r = rMiddle+1; r < rMax_minus_1; r++){
				auto const & valK_r = valK[r];
				for(uint i = 0; i < blockSize; i++)
					m.addConstr(m.getVarByName("k_"+to_string(r)+"_"+to_string(i)) == valK_r[i]);
			}

			//Fix the known round input
			for(int r = rMiddle+1; r < rMax_minus_1; r++){
				auto const & valX_r = valX[r];
				for(uint i = 0; i < blockSize; i++)
					m.addConstr(m.getVarByName("x_"+to_string(r)+"_"+to_string(i)) == valX_r[i]);
			}

			//Objective variables
			vector<GRBVar> XmidVar(blockSize);
			vector<GRBVar> KmidVar(blockSize);
			for(uint i = 0; i < blockSize; i++){
				XmidVar[i] = m.getVarByName("x_"+to_string(rMiddle)+"_"+to_string(i));
				KmidVar[i] = m.getVarByName("k_"+to_string(rMiddle)+"_"+to_string(i));
			}

			//Set the objective
			GRBLinExpr objExpr(0);
			for(uint i = 0; i < blockSize; i++){
				objExpr	+= XmidVar[i];
				objExpr -= m.getVarByName("y_"+to_string(rMiddle+1)+"_"+to_string(i));
				objExpr += 0.01*KmidVar[i];
			}
			m.setObjective(objExpr, GRB_MAXIMIZE);

			//Seed
			auto seedGurobi = randintSeedGurobi(prng);
			m.set(GRB_IntParam_Seed, seedGurobi);

			//Prepare the callback variables
			//SSB variables for the second half
			//index [r][i][j], r-th rounds (from rMiddle, e.g. 0 is for round rMiddle), i-th SSB, j-th variable
			vector<vector<vector<GRBVar>>> inSSBVar(rMax-rMiddle-1, 
													vector<vector<GRBVar>>(nbSSB,
													vector<GRBVar>(SSBSize)));
			vector<vector<vector<GRBVar>>> outSSBVar(rMax-rMiddle-1, 
													 vector<vector<GRBVar>>(nbSSB,
													 vector<GRBVar>(SSBSize)));
			vector<vector<vector<GRBVar>>> keySSBVar(rMax-rMiddle-1, 
													 vector<vector<GRBVar>>(nbSSB,
													 vector<GRBVar>(SSBSize)));

			for(uint r = 0; r < rMax-rMiddle-1; r++){
				auto & inSSBVar_r = inSSBVar[r];
				auto & outSSBVar_r = outSSBVar[r];
				auto & keySSBVar_r = keySSBVar[r];
				uint rpoffset = r + rMiddle;
				for(uint i = 0; i < nbSSB; i++){
					auto & inSSBVar_r_i = inSSBVar_r[i];
					auto & outSSBVar_r_i = outSSBVar_r[i];
					auto & keySSBVar_r_i = keySSBVar_r[i];
					for(uint j = 0; j < SSBSize; j++){
						uint bit = SSBSize*i + j;
						//input of the SSB
						inSSBVar_r_i[j] = m.getVarByName("x_"+to_string(rpoffset)+"_"+to_string(invPBitLevel[bit]));
						outSSBVar_r_i[j] = m.getVarByName("y_"+to_string(rpoffset+1)+"_"+to_string(bit));
						keySSBVar_r_i[j] = m.getVarByName("k_"+to_string(rpoffset)+"_"+to_string(bit));
					}
				}
			}

			//Time limit
			m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_ONESTEPDYNAMIC);

			//MIPFocus
			m.set(GRB_IntParam_MIPFocus, DYNAMIC_STEP1_MIPFOCUS);

			callbackDynamic cb(*this, rMiddle, output, valK, XmidVar, KmidVar, inSSBVar, outSSBVar, keySSBVar);
			m.setCallback(&cb);	
			//Optimize
			m.optimize();

			if(m.get(GRB_IntAttr_Status) == GRB_INFEASIBLE){
				cout << "Infeasible for rMiddle = " << rMiddle << ", starting over (if possible inputs remain...)" << endl;
				// if(indexInput != -1) //Infeasible with this input
				// 	return make_tuple(vector<uint8_t>(), vector<uint8_t>(), vector<vector<uint8_t>());
				// else
				// 	break;
				if(rMiddle == rMax_minus_2)
					possibleInput.erase(possibleInput.begin()+indexSelectedInput);
				break;
			}

			if(m.get(GRB_IntAttr_SolCount) == 0){
				//Didn't found any solution, start over
				cout << "Could not find a solution for rMiddle = " << rMiddle << ", starting over" << endl;
				break;
			}

			//Else, we have at least one solution, so extract it
			//The first solution in the solution pool is the best one
			auto & valX_rMiddle = valX[rMiddle];
			auto & valK_rMiddle = valK[rMiddle];
			for(uint i = 0; i < blockSize; i++){
				valX_rMiddle[i] = uint8_t(round(XmidVar[i].get(GRB_DoubleAttr_X)));
				valK_rMiddle[i] = uint8_t(round(KmidVar[i].get(GRB_DoubleAttr_X)));
			}

			cout << "Found a solution for rMiddle = " << rMiddle  << " in " << m.get(GRB_DoubleAttr_Runtime) << " seconds" << endl;
			cout << "x" << rMiddle << " : "; printParitySetHex(valX[rMiddle]); cout << endl;
			cout << "k" << rMiddle << " : "; printParitySetHex(valK[rMiddle]); cout << endl;

			//Check the hamming weight of Xmiddle to see if it exceed the threshold
			uint hw = 0;
			for(uint i = 0; i < blockSize; i++)
				hw += valX_rMiddle[i];
			cout << "Weight of x" << rMiddle << " : " << hw << endl;

			if(hw > WEIGHT_THRESHOLD_STEP1){
				//We found a middle round
				foundMiddle = true;
				goodrMiddle = rMiddle;
				//Also extract the input
				for(uint i = 0; i < blockSize; i++)
					input[i] = uint8_t(round(m.getVarByName("x_0_"+to_string(i)).get(GRB_DoubleAttr_X)));
				break;
			}
			//Else, keep going
		}

		//We either found a good middle round, or need to restart
		if(foundMiddle){
			cout << "Found a good middle round, starting from rMiddle = " << goodrMiddle << endl;
			cout << "Input : "; printParitySetHex(input); cout << endl;
			cout << "Keys : " << endl;
			for(int r = goodrMiddle; r < rMax_minus_1; r++){
				printParitySetHex(valK[r]); cout << endl;
			}
			cout << endl;
			cout << "Searching for key pattern from input to middle" << endl;
			//We now have found a possibly good middle round and key pattern for the second half
			//Now look for a key pattern in the first half
			auto Lsol_allsol = searchKeyPattern(name,goodrMiddle,input,valX[goodrMiddle],false);
			auto const & Lsol = Lsol_allsol.first;
			if(Lsol.size() > 0){
				//We found at least one solution for the first half
				for(auto const & firstKeys : Lsol){

					//Complete the current key pattern
					for(uint r = 0; r < goodrMiddle; r++){
						auto & valK_r = valK[r];
						auto const & firstKeys_r = firstKeys[r];
						for(uint i = 0; i < blockSize; i++)
							valK_r[i] = firstKeys_r[i];
					}

					cout << "Input  : "; printParitySetHex(input); cout << endl;
					cout << "Output : "; printParitySetHex(output); cout << endl;
					cout << "Keys   : " << endl;
					for(auto const & keypattern : valK){
						printParitySetHex(keypattern);
						cout << endl;
					}
					cout << "Counting trails..." << endl;
					//Now count the actual number of trails
					auto nbTrail_fullOpti = countTrailsFullCipher(input,output,valK);
					auto nbTrail = nbTrail_fullOpti.first;
					cout << nbTrail << " trails" << endl;
					if(nbTrail_fullOpti.second && nbTrail%2 == 1)
						return make_tuple(input,output,valK);
					else if(nbTrail_fullOpti.second && nbTrail%2 == 0)
						cout << "But even..." << endl;
					else
						cout << "But was not able to compute everything..." << endl;
					//Else try next key
				}
			}
			//Else restart
		}
		//Else restart
	}

	//If we get here, then no input can lead to at least one trail
	return make_tuple(vector<uint8_t>(), vector<uint8_t>(), vector<vector<uint8_t>>());
}

vector<uint> BCData::getLowerBounds(vector<uint8_t> const & input,
							vector<uint8_t> const & output,
							vector<vector<uint8_t>> const & keyVal){
	/*
		Compute lower bounds for 1 to #this.rMax-1 rounds using #input, #output and #keyval as helpers
		As described in the paper, the idea is to use the input/key/output obtained from the full degree proof to get lower bounds on a lower number of rounds
		Thus #input,#output,#keyval should describe the parity set needed to get an odd number of trails over #this.rMax rounds from #input to #output using #keyval for the keys
		Especially, #keyval should contain #this.rMax-1 vectors

		Note that technically, #input does not necessarily need to be of weight #this.blockSize-1 for this function to work, it only needs the total number of trails over #this.rMax rounds to be odd.
		Return a vector LB of length #this.rMax such that LB[r] is a lower bound on the degree for #r rounds
	*/

	//Generate the model
	generateBaseModel(name+"_lowerBounds",rMax,S,P,Mmodel,linAsSbox,keyAfterMC,false,true,dedicatedPRESENTlastLayer);

	int rMax_minus_1 = (int)(rMax) - 1;
	vector<uint> lowerBound(rMax,0);
	//For each round between rMax-1 and 1
	for(int rMiddle = rMax_minus_1; rMiddle > 0; rMiddle--){
		uint nbRounds = rMax - rMiddle;

		//Extract the keys
		vector<vector<uint8_t>> tmpKey(nbRounds, vector<uint8_t>(blockSize));
		for(uint r = 0; r < nbRounds-1; r++)
        	tmpKey[r] = keyVal[rMiddle+r];

        //First step, find the maximum input weight from rMiddle
		
		//Read the model
        GRBModel m(gurobiEnv,name+"_lowerBounds.mps");

        //Fix the input/output
        for(uint i = 0; i < blockSize; i++){
        	m.addConstr(m.getVarByName("x_"+to_string(0)+"_"+to_string(i)) == input[i]);
			m.addConstr(m.getVarByName("y_"+to_string(rMax_minus_1)+"_"+to_string(i)) == output[i]);
        }

        //Fix the keys
        for(int r = 0; r < rMax_minus_1; r++){
        	auto const & keyVal_r = keyVal[r];
        	for(uint i = 0; i < blockSize; i++)
        		m.addConstr(m.getVarByName("k_"+to_string(r)+"_"+to_string(i)) == keyVal_r[i]);
        }

        //Objective : maximize the weight of x[rMiddle]
        GRBLinExpr objExpr(0);
        for(uint i = 0; i < blockSize; i++)
        	objExpr += m.getVarByName("x_"+to_string(rMiddle)+"_"+to_string(i));
       	m.setObjective(objExpr, GRB_MAXIMIZE);

       	m.optimize();


       	uint maxWeight = uint(round(m.get(GRB_DoubleAttr_ObjVal)));

       	//We have the maximum weight, search for an input of this weight with an odd number of trails
       	//If none, decrease the weight
       	for(uint weight = maxWeight; weight > 0; weight--){

       		cout << nbRounds << " rounds, trying with weight " << weight << endl;
       		//The model is already created
       		//Read the model
        	GRBModel m(gurobiEnv,name+"_lowerBounds.mps");

        	//Fix the input/output
	        for(uint i = 0; i < blockSize; i++){
	        	m.addConstr(m.getVarByName("x_"+to_string(0)+"_"+to_string(i)) == input[i]);
				m.addConstr(m.getVarByName("y_"+to_string(rMax_minus_1)+"_"+to_string(i)) == output[i]);
	        }

        	//Fix the keys
	        for(int r = 0; r < rMax_minus_1; r++){
	        	auto const & keyVal_r = keyVal[r];
	        	for(uint i = 0; i < blockSize; i++)
	        		m.addConstr(m.getVarByName("k_"+to_string(r)+"_"+to_string(i)) == keyVal_r[i]);
	        }

	        //Input middle variables
	        vector<GRBVar> inputMiddleVar(blockSize);
	        for(uint i = 0; i < blockSize; i++)
	        	inputMiddleVar[i] = m.getVarByName("x_"+to_string(rMiddle)+"_"+to_string(i));

	        //Fix the weight of the input middle
	        GRBLinExpr inputExpr(0);
	       	for(uint i = 0; i < blockSize; i++)
	       		inputExpr += inputMiddleVar[i];
	       	m.addConstr(inputExpr == weight);

	       	//Dummy objective
	       	GRBLinExpr objExpr(0);
	       	for(uint i = 0; i < blockSize; i++)
	       		objExpr += m.getVarByName("x_"+to_string(nbRounds/2)+"_"+to_string(i));
	       	m.setObjective(objExpr, GRB_MAXIMIZE);

	       	//Setup the callback
	       	m.set(GRB_IntParam_LazyConstraints, 1);
	       	m.set(GRB_IntParam_SolutionLimit, 1);
	       	callbackLowerBound cb(*this,nbRounds,inputMiddleVar,output,tmpKey);
	       	m.setCallback(&cb);

	       	m.optimize();

	       	if(m.get(GRB_IntAttr_SolCount) > 0){
	       		cout << "Found good input of weight " << weight << " for " << nbRounds << " rounds : ";
	       		vector<uint8_t> inputVal(blockSize);
	       		for(uint i = 0; i < blockSize; i++)
	       			inputVal[i] = int(round(inputMiddleVar[i].get(GRB_DoubleAttr_X)));
	       		printParitySetHex(inputVal);
	       		cout << endl;
	       		lowerBound[nbRounds] = weight;
	       		break;
	       	}
	       	else
	       		cout << "No good input of weight " << weight << " for " << nbRounds << " rounds" << endl;
       	}
	}

	return lowerBound;
}

vector<uint> BCData::getUpperBoundAllOutput(){

	/*
		Get an upper bounds on the degree of each output bit for #this.rMax rounds
		Returns a vector UB of size #this.blockSize such that UB[i] is the upper bound for bit #i
	*/

	//Generate the model
	generateBaseModel(name+"_upperBounds",rMax,S,P,Mmodel,linAsSbox,keyAfterMC,false,true,dedicatedPRESENTlastLayer);

	int rMax_minus_1 = (int)(rMax) - 1;
	vector<uint> UB(blockSize,0);
	for(uint indexOutput = 0; indexOutput < blockSize; indexOutput++){
		//Read the model
	    GRBModel m(gurobiEnv,name+"_upperBounds.mps");

	    vector<uint8_t> output(blockSize, 0);
		output[indexOutput] = 1;

		//Fix the output
		for(uint i = 0; i < blockSize; i++)
			m.addConstr(m.getVarByName("y_"+to_string(rMax_minus_1)+"_"+to_string(i)) == output[i]);

		//Objective : maximize the weight of the input
		GRBLinExpr objExpr(0);
		for(uint i = 0; i < blockSize; i++)
			objExpr += m.getVarByName("x_"+to_string(0)+"_"+to_string(i));
		m.setObjective(objExpr, GRB_MAXIMIZE);

		//optimize
		m.optimize();

		UB[indexOutput] = uint(round(m.get(GRB_DoubleAttr_ObjVal)));
		cout << rMax << " rounds bit " << indexOutput << " : " << UB[indexOutput] << endl;
	}

	return UB;
}

void BCData::printParitySetHex(vector<uint8_t> const & x) const{
	/*
		Hex printer for more compact prints
		LSB of x is x[0]
		An hex character thus represent (x[i],...,x[i+3]) with x[i] LSB
		e.g.
		x = {1,0,0,0} is printed as 1
		x = {0,1,0,0,1,1,0,1} is printed as 2B
	*/
	for(uint s = 0; s < nbSbox; s++){
		uint val = 0;
		for(uint i = 0; i < sboxSize; i++)
			val |= (x[sboxSize*s + i] << i);
		cout << hex << uppercase << val << dec;
	}
}



//Non-member functions
uint selectNextOutputBit(SmallMatrix const & M){
	/*
	Help selecting the next output bit to consider
	If the matrix is empty, return a random column
	Else, return a random column that does not have a unit vector
	*/
	if(M.nrows == 0) return randintColumn(prng);
	else{
		//Check which columns have a unit vector
		vector<bool> unitCol(M.ncols, false);
		for(uint i = 0; i < M.nrows; i++){
			auto const row = M.get(i);
			if(__builtin_popcountll(row) == 1){
				uint indexUnit = 0;
				while(!(row & (uint64_t(1) << indexUnit))) indexUnit++;
				unitCol[indexUnit] = true;
			}
		}
		vector<uint> nonUnitCol;
		for(uint i = 0; i < M.ncols; i++){
			if(!unitCol[i]) nonUnitCol.emplace_back(i);
		}
		std::uniform_int_distribution<uint> tmprand(0,nonUnitCol.size()-1);

		return nonUnitCol[tmprand(prng)];
	}
}

void printBits(uint64_t const x, uint nbBits){
	//print the first #nbBits of x in cout
	uint64_t y = x;
	for(uint i = 0; i < nbBits; i++){
		cout << (y & 1);
		y >>= 1;
	}
}

void printVec(std::vector<uint8_t> const & x){
	//print the elements of x, wihtout separators
	for(auto const & xx : x) cout << int(xx);
}

string hex_char_to_bin(char c)
{
	//hex to bin conversion, LSB on the left in the binary representation
    switch(toupper(c))
    {
        case '0': return "0000";
        case '1': return "1000";
        case '2': return "0100";
        case '3': return "1100";
        case '4': return "0010";
        case '5': return "1010";
        case '6': return "0110";
        case '7': return "1110";
        case '8': return "0001";
        case '9': return "1001";
        case 'A': return "0101";
        case 'B': return "1101";
        case 'C': return "0011";
        case 'D': return "1011";
        case 'E': return "0111";
        case 'F': return "1111";
        default : return ""+c;
    }
}

vector<uint8_t> hexToVec(string const & s){
	/*
		From an hex string get the corresponding vector of bits
		LSB of each hex char is put on the left
		e.g. 
		s = "1"  returns {1,0,0,0}
		s = "2B" returns {0,1,0,0,1,1,0,1}
	*/

    string bin;
    for(unsigned i = 0; i != s.size(); ++i)
       bin += hex_char_to_bin(s[i]);
    
    vector<uint8_t> v(bin.size());
    for(uint i = 0; i < v.size(); i++){
    	if(bin[i] == '0') v[i] = 0;
    	else v[i] = 1;
    }

    return v;
}

