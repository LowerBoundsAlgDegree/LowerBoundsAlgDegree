#include"aes.h"

// add linear constraints for Sbox to model
void sboxConstr(GRBModel& model, vector<GRBVar> inVar, vector<GRBVar> outVar) {

	for (int i = 0; i < table.size(); i++) {
		GRBLinExpr tmp = 0;
		for (int j = 0; j < 8; j++) {
			if (table[i][j] == 0)
				tmp += inVar[j];
			else if (table[i][j] == 1)
				tmp += (1 - inVar[j]);
		}
		for (int j = 0; j < 8; j++) {
			if (table[i][8 + j] == 0)
				tmp += outVar[j];
			else if (table[i][8 + j] == 1)
				tmp += (1 - outVar[j]);
		}
		model.addConstr(tmp >= 1);
	}
}

// add linear constraints for MixColumns to model
void modelBinaryMat(GRBModel& model, vector<vector<GRBVar>> in, vector<vector<GRBVar>> out) {

	vector<GRBVar> v(32);
	vector<GRBVar> w(32);
	for (int i = 0; i < 8; i++) {
		v[i] = in[0][7 - i];
		v[8 + i] = in[1][7 - i];
		v[16 + i] = in[2][7 - i];
		v[24 + i] = in[3][7 - i];

		w[i] = out[0][7 - i];
		w[8 + i] = out[1][7 - i];
		w[16 + i] = out[2][7 - i];
		w[24 + i] = out[3][7 - i];
	}

	vector<vector<GRBVar>> element(32, vector<GRBVar>(32));
	for (int i = 0; i < 32; i++) {
		for (int j = 0; j < 32; j++) {
			if (mixcolumns[i][j] == 1) {
				element[i][j] = model.addVar(0, 1, 0, GRB_BINARY);
			}
		}
	}

	for (int j = 0; j < 32; j++) {
		if ((j % 8) == 0) {
			GRBVar tmp[11];
			int index = 0;
			for (int i = 0; i < 32; i++) {
				if (mixcolumns[i][j] == 1) {
					tmp[index] = element[i][j];
					index++;
				}
			}
			model.addGenConstrOr(v[j], tmp, 11);
		}
		else {
			GRBVar tmp[5];
			int index = 0;
			for (int i = 0; i < 32; i++) {
				if (mixcolumns[i][j] == 1) {
					tmp[index] = element[i][j];
					index++;
				}
			}
			model.addGenConstrOr(v[j], tmp, 5);
		}
	}

	for (int i = 0; i < 32; i++) {
		GRBLinExpr tmp = 0;
		for (int j = 0; j < 32; j++) {
			if (mixcolumns[i][j] == 1) {
				tmp += element[i][j];
			}
		}
		model.addConstr(w[i] == tmp);
	}


	// bypass constraints for 0xFFFFFFFF -> 0xFFFFFFFF
	GRBVar allOne = model.addVar(0, 1, 0, GRB_BINARY);
	GRBVar tmp[32] = {
		v[0], v[1], v[2], v[3], v[4], v[5], v[6], v[7], v[8], v[9], v[10], v[11], v[12], v[13], v[14], v[15],
		v[16], v[17], v[18], v[19], v[20], v[21], v[22], v[23], v[24], v[25], v[26], v[27], v[28], v[29], v[30], v[31]
	};
	model.addGenConstrAnd(allOne, tmp, 32);

	for (int i = 0; i < 32; i++) {
		model.addConstr(element[i][(i + 1) % 32] >= allOne);
	}

}

// set each pattern
void setPattern(GRBModel& model, vector<vector<GRBVar>> var, vector<int> v) {

	for (int i = 0; i < 16; i++) {
		for (int j = 0; j < 8; j++) {
			if (((v[i] >> j) & 1) == 1) {
				model.addConstr(var[i][j] == 1);
			}
			else {
				model.addConstr(var[i][j] == 0);
			}
		}
	}

}

// main function
int main(void) {

	// rounds that we show the lower bound
	int rounds = 5;

	// input/key/output pattern
	vector<int> in = { 0xFF, 0xFE, 0xFE, 0xFE, 0xFE, 0xFF, 0xEF, 0xFE, 0xFB, 0xF7, 0xFF, 0xF7, 0xFE, 0xEF, 0xF7, 0xFF };
	vector<vector<int>> key = {
		{0x00, 0x00, 0x00, 0x00, 0xBB, 0xFC, 0x3F, 0x77, 0xEF, 0xFC, 0xF2, 0x5F, 0x75, 0x7B, 0xFC, 0xFB},
		{0x77, 0xE8, 0x91, 0x9C, 0x2D, 0xE4, 0x36, 0x75, 0xC6, 0xD3, 0xE6, 0x8C, 0x4E, 0x3A, 0x4F, 0xA3},
		{0x83, 0x00, 0x00, 0x00, 0x00, 0x0D, 0x00, 0x00, 0x00, 0x00, 0x26, 0x00, 0x00, 0x00, 0x00, 0x62},
		{0xA1, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00}
	};
	vector<int> out = { 0x01, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00, 0x00 };



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
		vector<vector<vector<GRBVar>>> X(rounds, vector<vector<GRBVar>>(16, vector<GRBVar>(8)));
		vector<vector<vector<GRBVar>>> Y(rounds, vector<vector<GRBVar>>(16, vector<GRBVar>(8)));
		vector<vector<vector<GRBVar>>> Z(rounds - 1, vector<vector<GRBVar>>(16, vector<GRBVar>(8)));
		vector<vector<vector<GRBVar>>> W(rounds - 1, vector<vector<GRBVar>>(16, vector<GRBVar>(8)));
		vector<vector<vector<GRBVar>>> K(rounds - 1, vector<vector<GRBVar>>(16, vector<GRBVar>(8)));
		for (int r = 0; r < rounds; r++) {
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 8; j++) {
					X[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					Y[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
					if (r < rounds - 1) {
						W[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
						K[r][i][j] = model.addVar(0, 1, 0, GRB_BINARY);
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
				sboxConstr(model, X[r][i], Y[r][i]);

			//shiftrows
			for (int i = 0; i < 16; i++)
				Z[r][i] = Y[r][shiftrows[i]];

			//mixcolumns
			modelBinaryMat(model, { Z[r][0],Z[r][1],Z[r][2],Z[r][3] }, { W[r][0],W[r][1],W[r][2],W[r][3] });
			modelBinaryMat(model, { Z[r][4],Z[r][5],Z[r][6],Z[r][7] }, { W[r][4],W[r][5],W[r][6],W[r][7] });
			modelBinaryMat(model, { Z[r][8],Z[r][9],Z[r][10],Z[r][11] }, { W[r][8],W[r][9],W[r][10],W[r][11] });
			modelBinaryMat(model, { Z[r][12],Z[r][13],Z[r][14],Z[r][15] }, { W[r][12],W[r][13],W[r][14],W[r][15] });
			
			// addKey
			for (int i = 0; i < 16; i++) {
				for (int j = 0; j < 8; j++) {
					model.addConstr(X[r + 1][i][j] == W[r][i][j] + K[r][i][j]);
				}
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