#include"gurobi_c++.h"
#include<vector>

using namespace std;

// bit permutation of GIFT
vector<int> BP = {
	0,17,34,51,48,1,18,35,32,49,2,19,16,33,50,3,
	4,21,38,55,52,5,22,39,36,53,6,23,20,37,54,7,
	8,25,42,59,56,9,26,43,40,57,10,27,24,41,58,11,
	12,29,46,63,60,13,30,47,44,61,14,31,28,45,62,15
};

// linear inequalities for GIFT Sbox with logic condition.
// let (x0,x1,x2,x3) --> (y0,y1,y2,y3)
// the first one means  "x1 + x2 + (1-x3) + y0 + y1 + y2 >= 1"  
vector<vector<int>> table = {
	{2,0,0,1,0,0,0,2},	{2,0,1,1,0,0,2,0},	{1,0,2,2,0,1,1,0},	{2,0,1,1,0,1,0,2},
	{2,0,1,1,1,0,0,2},	{2,1,0,1,1,0,2,0},	{1,1,2,2,1,1,0,0},	{2,1,1,1,1,1,2,0},
	{0,1,2,2,0,2,0,0},	{0,2,1,2,0,0,0,2},	{1,2,2,0,0,0,2,0},	{1,2,0,0,2,2,1,0},
	{1,2,0,1,2,2,0,0},	{1,2,1,2,2,0,0,0},	{2,1,1,0,0,0,2,2},	{2,1,1,2,2,0,0,0},
	{2,2,1,0,0,2,0,1},	{1,1,1,2,0,2,2,0},	{1,1,2,1,2,0,2,0},	{1,1,2,2,2,0,0,1},
	{1,1,1,1,2,2,0,2},	{1,1,2,1,0,2,2,1},	{2,2,1,1,2,0,1,1},	{0,0,0,2,2,2,1,2},
	{0,2,0,0,2,1,2,2},	{2,0,0,0,2,2,1,2},	{0,0,2,2,1,1,2,2},	{0,2,0,2,2,1,1,2},
	{0,2,2,0,1,2,1,2},	{2,0,0,2,1,1,2,2},	{2,2,0,0,1,2,2,1},	{2,2,0,2,2,1,1,1},
	{0,0,2,2,2,2,2,1},	{0,2,2,0,2,2,2,1},	{0,2,2,2,2,2,1,1}
};

// add linear constraints for Sbox to model
void sboxConstr(GRBModel& model, vector<GRBVar> inVar, vector<GRBVar> outVar) {

	for (int i = 0; i < table.size(); i++) {
		GRBLinExpr tmp = 0;
		for (int j = 0; j < 4; j++) {
			if (table[i][j] == 0)
				tmp += inVar[j];
			else if (table[i][j] == 1)
				tmp += (1 - inVar[j]);
		}
		for (int j = 0; j < 4; j++) {
			if (table[i][4 + j] == 0)
				tmp += outVar[j];
			else if (table[i][4 + j] == 1)
				tmp += (1 - outVar[j]);
		}
		model.addConstr(tmp >= 1);
	}
}

// set each pattern
void setPattern(GRBModel& model, vector<vector<GRBVar>> var, unsigned long long int v) {
	for (int i = 0; i < 64; i++) {
		if (((v >> i) & 1) == 1) {
			model.addConstr(var[i / 4][i % 4] == 1);
		}
		else {
			model.addConstr(var[i / 4][i % 4] == 0);
		}
	}
}

// main function
int main(void) {

	// rounds that we show the lower bound
	int rounds = 9;

	// input/key/output pattern
	unsigned long long int in = 0xFFFFFFFFFFFEFFFF;
	vector<unsigned long long int> key = {
		0x0000000000000000,
		0x0000000000000000,
		0x8004000000041000,
		0x9204860410519C00,
		0x060D130120091011,
		0x6004000000000001,
		0x0000000000000000,
		0x0000000000000000
	};
	unsigned long long int out = 0x8000000000000000;

	// start gurobi 
	try {
		// create the environment
		GRBEnv env = GRBEnv();

		// enumurate solutions up to 2000000000
		env.set(GRB_IntParam_PoolSearchMode, 2);
		env.set(GRB_IntParam_PoolSolutions, 2000000000);
		env.set(GRB_DoubleParam_PoolGap, GRB_INFINITY);

		// create the model
		GRBModel model = GRBModel(env);

		// create variables
		vector<vector<vector<GRBVar>>> X(rounds, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		vector<vector<vector<GRBVar>>> Y(rounds, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		vector<vector<vector<GRBVar>>> K(rounds - 1, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		for (int r = 0; r < rounds; r++) {
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 4; j++) {
					X[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					Y[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					if (r < rounds - 1)
						K[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
				}
			}
		}

		// set pattern
		setPattern(model, X[0], in);
		setPattern(model, Y[rounds - 1], out);
		for (int r = 0; r < rounds - 1; r++)
			setPattern(model, K[r], key[r]);

		// create constraints
		for (int r = 0; r < rounds - 1; r++) {

			//sbox
			for (int i = 0; i < 16; i++)
				sboxConstr(model, X[r][i], Y[r][i]);

			// key xor and bit perm
			for (int i = 0; i < 64; i++) {
				int j = BP[i];
				model.addConstr(Y[r][i / 4][i % 4] + K[r][j / 4][j % 4] == X[r + 1][j / 4][j % 4]);
			}

		}

		// last sbox
		for (int i = 0; i < 16; i++)
			sboxConstr(model, X[rounds - 1][i], Y[rounds - 1][i]);

		// solve this model 
		model.optimize();

		int solcount = model.get(GRB_IntAttr_SolCount);
		cout << "the number of trails : " << solcount << endl;

	}
	catch (GRBException e) {
		cerr << "Error code = " << e.getErrorCode() << endl;
		cerr << e.getMessage() << endl;
	}
	catch (...) {
		cerr << "Exception during optimization" << endl;
	}

	return 0;
}