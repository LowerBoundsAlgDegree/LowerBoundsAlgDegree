#include "customCallback.hpp"

#define TIME_LIMIT_CBSEARCHKEY 60		/*Time limit (in seconds) when counting the actual number of solutions in callbackSearchKeyPattern*/
#define MAX_NUMBER_CBSEARCHKEY 100000	/*Max number of solution allowed when counting the actual number of solutions in callbackSearchKeyPattern*/
#define MAX_GRB_THREAD 8 				/*Max number of thread to use by Gurobi. Setting it to 0 let Gurobi decide it itself*/

#define TIME_LIMIT_CALLBACKDYNAMIC 60	/*Time limit (in seconds) when counting the actual number of solutions in callbackDynamic*/
#define THRESHOLD_NBSOL_DYNAMIC 1		/*Max number of solution allowed when counting the actual number of solutions in callbackDynamic*/
#define THRESHOLD_NBSOL_SSB_DYNAMIC 1	/*Max number of solution allowed in one SSB when counting the actual number of solutions in callbackDynamic*/

//Note that if THRESHOLD_NBSOL_DYNAMIC > 1, the actual threshold is set to the number of rounds times THRESHOLD_NBSOL_DYNAMIC to manage the induced growth in the number of trails

/*
	Macros are used in the code to define different limits on the search (e.g. time limits)
	It's a bit dirty but easier than handling a ton of parameters everywhere
	Our results were obtained with the current values
*/

using namespace std;
using namespace boost;

typedef unsigned int uint;

callbackSearchKeyPattern::callbackSearchKeyPattern(BCData const & xBCD,
						uint const xnbRounds,
						vector<GRBVar> const & xinputVar,
						vector<GRBVar> const & xoutputVar,
						vector<vector<GRBVar>> const & xallKeyVar,
                        vector<vector<vector<GRBVar>>> const & xinSSBVar,
                        vector<vector<vector<GRBVar>>> const & xoutSSBVar,
                        vector<vector<vector<GRBVar>>> const & xkeySSBVar,
                        bool const xstartAfterSSB) :
	BCD(xBCD),
	nbRounds(xnbRounds),
	inputVar(xinputVar),
	outputVar(xoutputVar),
	allKeyVar(xallKeyVar),
	inSSBVar(xinSSBVar),
	outSSBVar(xoutSSBVar),
	keySSBVar(xkeySSBVar),
	ctrsolutions(0),
	allSolutions(),
	startAfterSSB(xstartAfterSSB)
{}

void callbackSearchKeyPattern::callback(){
	try {
        if (where == GRB_CB_MIPSOL) { //If the solver found a solution

        	//Increment the solution counter
        	ctrsolutions++;
        	//cout << ctrsolutions << " solutions up to now\r" << flush;

        	//Get the solution for each SSB
        	vector<vector<vector<uint8_t>>> valInSSB(nbRounds-1,
        											vector<vector<uint8_t>>(BCD.nbSSB,
        											vector<uint8_t>(BCD.SSBSize)));
        	vector<vector<vector<uint8_t>>> valOutSSB(nbRounds-1,
        											vector<vector<uint8_t>>(BCD.nbSSB,
        											vector<uint8_t>(BCD.SSBSize)));
        	vector<vector<vector<uint8_t>>> valKeySSB(nbRounds-1,
        											vector<vector<uint8_t>>(BCD.nbSSB,
        											vector<uint8_t>(BCD.SSBSize)));
        	for(uint r = 0; r < nbRounds-1; r++){
        		auto & valInSSB_r = valInSSB[r];
        		auto & valOutSSB_r = valOutSSB[r];
        		auto & valKeySSB_r = valKeySSB[r];

        		auto & inSSBVar_r = inSSBVar[r];
        		auto & outSSBVar_r = outSSBVar[r];
        		auto & keySSBVar_r = keySSBVar[r];
        		for(uint i = 0; i < BCD.nbSSB; i++){
        			auto & valInSSB_r_i = valInSSB_r[i];
	        		auto & valOutSSB_r_i = valOutSSB_r[i];
	        		auto & valKeySSB_r_i = valKeySSB_r[i];

	        		auto & inSSBVar_r_i = inSSBVar_r[i];
	        		auto & outSSBVar_r_i = outSSBVar_r[i];
	        		auto & keySSBVar_r_i = keySSBVar_r[i];
        			for(uint j = 0; j < BCD.SSBSize; j++){
        				valInSSB_r_i[j] = int(round(getSolution(inSSBVar_r_i[j])));
        				valOutSSB_r_i[j] = int(round(getSolution(outSSBVar_r_i[j])));
        				valKeySSB_r_i[j] = int(round(getSolution(keySSBVar_r_i[j])));
        			}
        		}
        	}

        	//Get the full key solution
        	uint const allKeyVar_size = allKeyVar.size();
        	vector<dynamic_bitset<>> keySolution(allKeyVar_size, 
        										 dynamic_bitset<>(BCD.blockSize,0));
			for(uint r = 0; r < allKeyVar_size; r++){
        		auto & allKeyVar_r = allKeyVar[r];
        		auto & keySolution_r = keySolution[r];
        		for(uint i = 0; i < BCD.blockSize; i++){
        			keySolution_r[i] = uint(round(getSolution(allKeyVar_r[i])));
        		}
        	}
        	//save the solution
	        allSolutions.emplace_back(keySolution);

        	//Check if this key pattern can lead to an odd number of trails by checking all SSB
        	bool oddSSB = true;
        	for(uint r = 0; r < nbRounds-1; r++){
        		auto & valInSSB_r = valInSSB[r];
        		auto & valOutSSB_r = valOutSSB[r];
        		auto & valKeySSB_r = valKeySSB[r];
        		for(uint i = 0; i < BCD.nbSSB; i++){
        			auto & valInSSB_r_i = valInSSB_r[i];
        			auto & valKeySSB_r_i = valKeySSB_r[i];
        			auto & valOutSSB_r_i = valOutSSB_r[i];
        			uint64_t nbTrail;
        			if(r == nbRounds-2 && BCD.dedicatedPRESENTlastLayer && startAfterSSB){
        				//Dedicated SSB for the last SSB if we modified the last layer for Present
        				nbTrail = BCD.countSSBTrail_dedicatedPresent(valInSSB_r_i, valKeySSB_r_i, valOutSSB_r_i);
        			}
        			else
        				nbTrail = BCD.countSSBTrail(valInSSB_r_i, valKeySSB_r_i, valOutSSB_r_i);
        			if(nbTrail%2 == 0){ //Even number of trails in this SSB, remove it
        				oddSSB = false;
        				auto & inSSBVar_r_i = inSSBVar[r][i];
		        		auto & outSSBVar_r_i = outSSBVar[r][i];
		        		auto & keySSBVar_r_i = keySSBVar[r][i];
        				GRBLinExpr cutExpr(0);
        				for(uint j = 0; j < BCD.SSBSize; j++){
        					if(valInSSB_r_i[j] == 0) cutExpr += inSSBVar_r_i[j];
        					else cutExpr += (1 - inSSBVar_r_i[j]);
        					if(valKeySSB_r_i[j] == 0) cutExpr += keySSBVar_r_i[j];
        					else cutExpr += (1 - keySSBVar_r_i[j]);
        					if(valOutSSB_r_i[j] == 0) cutExpr += outSSBVar_r_i[j];
        					else cutExpr += (1 - outSSBVar_r_i[j]);
        				}
        				addLazy(cutExpr >= 1);
        			}
	        	}
	        }


	        if(oddSSB){
	        	//All SSB trails are odd, count the actual number of trails
	        	//Read the base model
				GRBModel m(BCD.gurobiEnv, BCD.name+".mps");
				//Add input/output constraints
				for(uint i = 0; i < BCD.blockSize; i++){
					m.addConstr(m.getVarByName(inputVar[i].get(GRB_StringAttr_VarName)) == int(round(getSolution(inputVar[i]))));
					m.addConstr(m.getVarByName(outputVar[i].get(GRB_StringAttr_VarName)) == int(round(getSolution(outputVar[i]))));

				}
				//Add key constraints
				for(uint r = 0; r < allKeyVar_size; r++){
	        		auto & allKeyVar_r = allKeyVar[r];
	        		auto & keySolution_r = keySolution[r];
	        		for(uint i = 0; i < BCD.blockSize; i++){
	        			m.addConstr(m.getVarByName(allKeyVar_r[i].get(GRB_StringAttr_VarName)) == uint(keySolution_r[i]));
	        		}
	        	}

	        	//setup params for enumerating solutions
	        	cout << endl << "counting the actual number of trails..." << endl;
	        	m.set(GRB_IntParam_PoolSearchMode, 2);
				m.set(GRB_IntParam_PoolSolutions, 2000000000); //maximum number of solution allowed by Gurobi

				//Dummy objective
				GRBLinExpr objExpr(0);
				for(uint i = 0; i < BCD.blockSize; i++){
					objExpr += m.getVarByName("x_"+to_string(nbRounds/2)+"_"+to_string(i));
				}
				m.setObjective(objExpr);
				
				//Set the counting callback and optimize
				// callbackCount cb = callbackCount(MAX_NUMBER_CBSEARCHKEY);
				// m.setCallback(&cb);
				m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_CBSEARCHKEY);
				m.set(GRB_IntParam_SolutionLimit, MAX_NUMBER_CBSEARCHKEY);
				m.set(GRB_IntParam_Threads,MAX_GRB_THREAD);
				
				m.optimize();

				if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
					uint nbTrail = m.get(GRB_IntAttr_SolCount);
					cout << endl << nbTrail << " different trails" << endl;

					//If the actual number of trails is even, we need another key
					if(nbTrail%2 == 0) oddSSB = false;
					else{
						uint weightSol = 0;
						cout << endl << "odd solution " << oddSSB << endl;
						for(auto const & ksol : keySolution){
							for(uint i = 0; i < BCD.blockSize; i++){
								cout << ksol[i];
								if(ksol[i] == 1) weightSol++;
							}
							cout << endl;
						}
						cout << "weight " << weightSol << endl;
					}
				}
				else{
					cout << endl << "Too many solutions, aborted" << endl;
					oddSSB = false;
				}

				if(!oddSSB){
			        //If we get here, the key pattern is not good, search for another one
			        //Generate a linear expression to remove he current key pattern from the solution pool
			        GRBLinExpr cutExpr(0);
			        for(uint r = 0; r < allKeyVar_size; r++){
		        		auto & allKeyVar_r = allKeyVar[r];
		        		auto & keySolution_r = keySolution[r];
		        		for(uint i = 0; i < BCD.blockSize; i++){
		        			if(keySolution_r[i] == 0) cutExpr += allKeyVar_r[i];
		        			else cutExpr += (1 - allKeyVar_r[i]);
		        		}
		        	}
		        	addLazy(cutExpr >= 1);
		        }
	        }
	        //If number of trails is odd, the solution is not removed and thus the solver finishes
        }
      } catch (GRBException e) {
        cout << "Error number: " << e.getErrorCode() << endl;
        cout << e.getMessage() << endl;
      } catch (...) {
        cout << "Error during callback" << endl;
      }
}

callbackCount::callbackCount() : 
	ctr(0),
	countLimit(0),
	wasAborted(false),
	printNbSol(false)
{}

callbackCount::callbackCount(uint64_t const xcountLimit,
							 bool xprintNbSol) :
	ctr(0),
	countLimit(xcountLimit),
	wasAborted(false),
	printNbSol(xprintNbSol)
{}

void callbackCount::callback(){
	try {
        if (where == GRB_CB_MIPSOL) { //If the solver found a solution
		    ctr++;
		    if(printNbSol)
		    	std::cout << ctr << " solutions up to now in callbackCount\r" << std::flush;
		    if(countLimit > 0 && ctr > countLimit){
		    	wasAborted = true;
		    	abort();
		    }
		}
	} catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
	} catch (...) {
	cout << "Error during callback" << endl;
	}
}


callbackDynamic::callbackDynamic(BCData const & xBCD,
								 uint const xrMiddle,
                    			 vector<uint8_t> const & xoutput,
                    			 vector<vector<uint8_t>> const & xkeyVal,
                    			 vector<GRBVar> const & xinputMiddleVar,
                    			 vector<GRBVar> const & xkeyMiddleVar,
                    			 vector<vector<vector<GRBVar>>> const & xinSSBVar,
                    			 vector<vector<vector<GRBVar>>> const & xoutSSBVar,
                    			 vector<vector<vector<GRBVar>>> const & xkeySSBVar) :
	BCD(xBCD),
	rMiddle(xrMiddle),
	output(xoutput),
	keyVal(xkeyVal),
	inputMiddleVar(xinputMiddleVar),
	keyMiddleVar(xkeyMiddleVar),
	inSSBVar(xinSSBVar),
	outSSBVar(xoutSSBVar),
	keySSBVar(xkeySSBVar)
{}

void callbackDynamic::callback(){
	try {
        if (where == GRB_CB_MIPSOL) { //If the solver found a solution
		    
        	//Get the solution for each SSB
        	vector<vector<vector<uint8_t>>> valInSSB(inSSBVar.size(),
        											vector<vector<uint8_t>>(BCD.nbSSB,
        											vector<uint8_t>(BCD.SSBSize)));
        	vector<vector<vector<uint8_t>>> valOutSSB(outSSBVar.size(),
        											vector<vector<uint8_t>>(BCD.nbSSB,
        											vector<uint8_t>(BCD.SSBSize)));
        	vector<vector<vector<uint8_t>>> valKeySSB(keySSBVar.size(),
        											vector<vector<uint8_t>>(BCD.nbSSB,
        											vector<uint8_t>(BCD.SSBSize)));
        	for(uint r = 0; r < inSSBVar.size(); r++){
        		auto & valInSSB_r = valInSSB[r];
        		auto & valOutSSB_r = valOutSSB[r];
        		auto & valKeySSB_r = valKeySSB[r];

        		auto & inSSBVar_r = inSSBVar[r];
        		auto & outSSBVar_r = outSSBVar[r];
        		auto & keySSBVar_r = keySSBVar[r];
        		for(uint i = 0; i < BCD.nbSSB; i++){
        			auto & valInSSB_r_i = valInSSB_r[i];
	        		auto & valOutSSB_r_i = valOutSSB_r[i];
	        		auto & valKeySSB_r_i = valKeySSB_r[i];

	        		auto & inSSBVar_r_i = inSSBVar_r[i];
	        		auto & outSSBVar_r_i = outSSBVar_r[i];
	        		auto & keySSBVar_r_i = keySSBVar_r[i];
        			for(uint j = 0; j < BCD.SSBSize; j++){
        				valInSSB_r_i[j] = int(round(getSolution(inSSBVar_r_i[j])));
        				valOutSSB_r_i[j] = int(round(getSolution(outSSBVar_r_i[j])));
        				valKeySSB_r_i[j] = int(round(getSolution(keySSBVar_r_i[j])));
        			}
        		}
        	}

        	//Check if this key pattern can lead to an odd number of trails by checking all SSB
        	bool oddSSB = true;
        	// cout << "Number of trails in the SSBs : ";
        	for(uint r = 0; r < inSSBVar.size(); r++){
        		auto & valInSSB_r = valInSSB[r];
        		auto & valOutSSB_r = valOutSSB[r];
        		auto & valKeySSB_r = valKeySSB[r];
        		for(uint i = 0; i < BCD.nbSSB; i++){
        			auto & valInSSB_r_i = valInSSB_r[i];
        			auto & valKeySSB_r_i = valKeySSB_r[i];
        			auto & valOutSSB_r_i = valOutSSB_r[i];
        			uint64_t nbTrail;
        			if(r == inSSBVar.size()-1 && BCD.dedicatedPRESENTlastLayer){
        				//Dedicated SSB for the last SSB if we modified the last layer for Present
        				nbTrail = BCD.countSSBTrail_dedicatedPresent(valInSSB_r_i, valKeySSB_r_i, valOutSSB_r_i);
        			}
        			else
        				nbTrail = BCD.countSSBTrail(valInSSB_r_i, valKeySSB_r_i, valOutSSB_r_i);
        			// cout << nbTrail << " ";
        			//Check if it exceeds the threshold or even
        			if(nbTrail > THRESHOLD_NBSOL_SSB_DYNAMIC || nbTrail%2 == 0){
        				oddSSB = false;
        				auto & inSSBVar_r_i = inSSBVar[r][i];
		        		auto & outSSBVar_r_i = outSSBVar[r][i];
		        		auto & keySSBVar_r_i = keySSBVar[r][i];
        				GRBLinExpr cutExpr(0);
        				for(uint j = 0; j < BCD.SSBSize; j++){
        					if(valInSSB_r_i[j] == 0) cutExpr += inSSBVar_r_i[j];
        					else cutExpr += (1 - inSSBVar_r_i[j]);
        					if(valKeySSB_r_i[j] == 0) cutExpr += keySSBVar_r_i[j];
        					else cutExpr += (1 - keySSBVar_r_i[j]);
        					if(valOutSSB_r_i[j] == 0) cutExpr += outSSBVar_r_i[j];
        					else cutExpr += (1 - outSSBVar_r_i[j]);
        				}
        				addLazy(cutExpr >= 1);
        			}
	        	}
	        }
	        // cout << endl;

	        if(oddSSB){
	        	//All SSB in the second half are odd, count the actual number of trail in the second half
	        	//Generate a model for the second half
	        	generateBaseModel(BCD.name+"_lastHalf", BCD.rMax - rMiddle, BCD.S,BCD.P,BCD.Mmodel,BCD.linAsSbox,BCD.keyAfterMC,false,true,BCD.dedicatedPRESENTlastLayer);

	        	//Read the model
	        	GRBModel m(BCD.gurobiEnv,BCD.name+"_lastHalf.mps");

	        	//Add the input/output constraints
	        	for(uint i = 0; i < BCD.blockSize; i++){
	        		m.addConstr(m.getVarByName("x_0_"+to_string(i)) == int(round(getSolution(inputMiddleVar[i]))));
	        		m.addConstr(m.getVarByName("y_"+to_string(BCD.rMax - rMiddle - 1)+"_"+to_string(i)) == output[i]);
	        	}

	        	//Add the key constraints
	        	//First key is extracted from current solution
	        	vector<uint8_t> keyMiddleVal(BCD.blockSize);
	        	for(uint i = 0; i < BCD.blockSize; i++)
	        		keyMiddleVal[i] = int(round(getSolution(keyMiddleVar[i])));
	        	for(uint i = 0; i < BCD.blockSize; i++)
	        		m.addConstr(m.getVarByName("k_0_"+to_string(i)) == keyMiddleVal[i]);
	        	//The remaining keys are already known
	        	for(uint r = rMiddle+1; r < BCD.rMax-1; r++){
	        		uint const offset = r - rMiddle;
	        		auto const & keyVal_r = keyVal[r];
	        		for(uint i = 0; i < BCD.blockSize; i++)
	        			m.addConstr(m.getVarByName("k_"+to_string(offset)+"_"+to_string(i)) == keyVal_r[i]);
	        	}

	        	//Dummy objective
				GRBLinExpr objExpr(0);
				for(uint i = 0; i < BCD.blockSize; i++){
					objExpr += m.getVarByName("x_"+to_string((BCD.rMax - rMiddle)/2)+"_"+to_string(i));
				}
				m.setObjective(objExpr);

				//setup params for enumerating solutions
	        	cout << "Counting the actual number of trails..." << endl;
	        	m.set(GRB_IntParam_PoolSearchMode, 2);
	        	m.set(GRB_IntParam_PoolSolutions, 2000000000); //maximum number of solution allowed by Gurobi
	        	uint64_t threshold = 1;
	        	if(THRESHOLD_NBSOL_DYNAMIC > 1)
	        		threshold = (BCD.rMax - rMiddle)*THRESHOLD_NBSOL_DYNAMIC;
	        	m.set(GRB_IntParam_SolutionLimit, threshold+1); //We stop if we find more solutions than the threshold
	        	m.set(GRB_IntParam_Threads,MAX_GRB_THREAD);
				m.set(GRB_DoubleParam_TimeLimit, TIME_LIMIT_CALLBACKDYNAMIC);

				m.optimize();

				if(m.get(GRB_IntAttr_Status) == GRB_OPTIMAL){
					uint nbTrail = m.get(GRB_IntAttr_SolCount);
					cout << nbTrail << " different trails" << endl;

					//If the actual number of trails is even or larger than the threshold, we need another key
					if(nbTrail > threshold || nbTrail%2 == 0) oddSSB = false;
				}
				else{
					cout << "Too many or no solutions (" << m.get(GRB_IntAttr_SolCount) << ")" << endl;
					oddSSB = false;
				}

				if(!oddSSB){
			        //If we get here, the key pattern is not good, search for another one
			        //Generate a linear expression to remove he current key pattern from the solution pool
			        GRBLinExpr cutExpr(0);
			        for(uint i = 0; i < BCD.blockSize; i++){
			        	if(keyMiddleVal[i] == 0) cutExpr += keyMiddleVar[i];
			        	else cutExpr += (1 - keyMiddleVar[i]);
			        }
		        	addLazy(cutExpr >= 1);
		        }
	        }

		}
	} catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
	} catch (...) {
	cout << "Error during callback" << endl;
	}
}


callbackLowerBound::callbackLowerBound(BCData const & xBCD,
				                       uint const xnbRounds,
				                       vector<GRBVar> const & xinputMiddleVar,
				                       vector<uint8_t> const & xoutput,
				                       vector<vector<uint8_t>> const & xkeyVal) :
	BCD(xBCD),
	nbRounds(xnbRounds),
	inputMiddleVar(xinputMiddleVar),
	output(xoutput),
	keyVal(xkeyVal)
{}


void callbackLowerBound::callback(){
	try {
        if (where == GRB_CB_MIPSOL) { //If the solver found a solution
		    //Get the value of the input
		    vector<uint8_t> input(BCD.blockSize);
		    for(uint i = 0; i < BCD.blockSize; i++)
		    	input[i] = int(round(getSolution(inputMiddleVar[i])));

		    cout << "Found input : "; BCD.printParitySetHex(input); cout << endl;

		    //Count the number of trails from input to output
		    //Generate the model
			generateBaseModel(BCD.name+"_cblb",nbRounds,BCD.S,BCD.P,BCD.Mmodel,BCD.linAsSbox,BCD.keyAfterMC,false,true,BCD.dedicatedPRESENTlastLayer);
			//Read the model
	        GRBModel m(BCD.gurobiEnv,BCD.name+"_cblb.mps");

	        //Fix the input/output
	        for(uint i = 0; i < BCD.blockSize; i++){
	        	m.addConstr(m.getVarByName("x_"+to_string(0)+"_"+to_string(i)) == input[i]);
        		m.addConstr(m.getVarByName("y_"+to_string(nbRounds-1)+"_"+to_string(i)) == output[i]);
	        }

	        //Fix the keys
	        for(uint r = 0; r < nbRounds-1; r++){
	        	auto const & keyVal_r = keyVal[r];
	        	for(uint i = 0; i < BCD.blockSize; i++)
	        		m.addConstr(m.getVarByName("k_"+to_string(r)+"_"+to_string(i)) == keyVal_r[i]);
	        }

	        m.set(GRB_IntParam_PoolSearchMode, 2);
	        m.set(GRB_IntParam_PoolSolutions, 2000000000);

	        //callback to see the current number of solution
	        callbackCount cb = callbackCount(0,false);
	        m.setCallback(&cb);

	        m.optimize();

	        //Get the number of solution
	        uint nbTrail = m.get(GRB_IntAttr_SolCount);
	        cout << nbTrail << " trails" << endl;
	        if(nbTrail%2 == 0){
	        	//Number of trail is even, remove this input
	        	GRBLinExpr cutExpr(0);
	        	for(uint i = 0; i < BCD.blockSize; i++){
	        		if(input[i] == 0) cutExpr += inputMiddleVar[i];
	        		else cutExpr += (1 - inputMiddleVar[i]);
	        	}
	        	addLazy(cutExpr >= 1);
	        }


		}
	} catch (GRBException e) {
    cout << "Error number: " << e.getErrorCode() << endl;
    cout << e.getMessage() << endl;
	} catch (...) {
	cout << "Error during callback" << endl;
	}
}
