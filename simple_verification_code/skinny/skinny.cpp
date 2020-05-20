#include"gurobi_c++.h"
#include<vector>

using namespace std;

// linear inequalities for Skinny Sbox with logic condition.
vector<vector<int>> s_table = {
	{0,1,0,1,0,2,1,0},	{2,0,1,0,2,1,0,0},	{0,1,1,0,1,2,2,0},	{2,1,1,1,1,1,0,2},
	{2,2,1,2,0,0,0,0},	{2,2,2,1,0,0,0,0},	{2,0,1,1,2,0,2,0},	{1,2,1,2,2,0,0,1},
	{1,1,1,2,0,2,0,2},	{1,1,2,1,2,0,0,2},	{2,2,1,1,2,0,1,0},	{1,1,2,1,0,2,2,1},
	{1,2,1,1,2,2,1,0},	{2,0,0,0,1,2,2,2},	{1,2,2,2,0,0,2,0},	{2,1,2,2,0,0,0,2},
	{0,2,0,0,2,1,2,2},	{0,2,2,0,2,1,1,2},	{0,0,2,2,2,1,1,2},	{2,2,0,0,1,1,2,2},
	{1,2,2,1,0,2,0,2},	{1,2,1,1,2,0,2,2},	{0,2,2,2,1,2,1,1},	{2,2,2,0,2,1,1,1},
	{1,2,2,2,1,2,1,0},	{2,2,0,2,1,1,1,2},	{2,2,0,2,2,1,2,1},	{2,0,2,2,1,2,2,1},
	{2,2,2,0,1,2,1,2},	{0,2,2,2,2,1,2,1},	{2,0,2,2,1,2,1,2}
};

// linear inequalities for Skinny Lbox with logic condition.
vector<vector<int>> l_table = {
	{2,0,2,0,1,1,1,2},	{1,2,2,2,0,0,2,0},	{2,2,1,2,0,2,0,0},	{0,2,0,0,1,2,2,2},
	{1,2,1,2,2,0,0,2},	{1,2,2,1,2,0,2,0},	{2,1,1,2,0,2,2,0},	{2,2,1,1,2,2,0,0},
	{2,0,2,2,2,1,1,1},	{1,1,1,2,2,0,2,2},	{2,1,1,1,2,2,2,0},	{0,2,0,2,2,2,2,1},
	{2,0,0,2,2,2,1,2},	{0,2,2,2,2,1,2,2},	{2,1,2,2,2,2,0,2},	{2,2,2,1,0,2,2,2},
};

//
vector<int> shiftrows = { 0,13,10,7,4,1,14,11,8,5,2,15,12,9,6,3 };

// add linear constraints for Sbox to model
void sboxConstr(GRBModel& model, vector<GRBVar> inVar, vector<GRBVar> outVar, vector<vector<int>> table) {

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
			model.addConstr(var[15 - (i / 4)][i % 4] == 1);
		}
		else {
			model.addConstr(var[15 - (i / 4)][i % 4] == 0);
		}
	}
}

// main function
int main(void) {

	// rounds that we show the lower bound
	int rounds = 10;

	// input/key/output pattern
	unsigned long long int in = 0xFFFFFFFFFFFFFEFF;
	vector<unsigned long long int> key = {
		0x0000000000000000,
		0x0460000000000000,
		0x000004000000C000,
		0x80080000C4C40880,
		0x008000008000000C,
		0x08004080C00C4080,
		0x0000C040004800C0,
		0x0000400000000000,
		0x0000000000000000
	};
	unsigned long long int out = 0x1000000000000000;


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
		vector<vector<vector<GRBVar>>> Z(rounds - 1, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		vector<vector<vector<GRBVar>>> W(rounds - 1, vector<vector<GRBVar>>(16, vector<GRBVar>(4)));
		for (int r = 0; r < rounds; r++) {
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 4; j++) {
					X[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					Y[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					if (r < rounds - 1) {
						K[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
						Z[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					}
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
				sboxConstr(model, X[r][i], Y[r][i], s_table);

			//addKey
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 4; j++) {
					model.addConstr(Z[r][i][j] == Y[r][i][j] + K[r][i][j]);
				}
			}

			//shiftrows
			for (int i = 0; i < 16; i++)
				W[r][i] = Z[r][shiftrows[i]];


			//mixcolumns
			for (int col = 0; col < 4; col++) {
				for (int j = 0; j < 4; j++) {
					vector<GRBVar> inVar = { W[r][4 * col + 0][j], W[r][4 * col + 1][j], W[r][4 * col + 2][j], W[r][4 * col + 3][j] };
					vector<GRBVar> outVar = { X[r + 1][4 * col + 0][j], X[r + 1][4 * col + 1][j], X[r + 1][4 * col + 2][j], X[r + 1][4 * col + 3][j] };
					sboxConstr(model, inVar, outVar, l_table);
				}
			}

		}

		// last sbox
		for (int i = 0; i < 16; i++)
			sboxConstr(model, X[rounds - 1][i], Y[rounds - 1][i], s_table);

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