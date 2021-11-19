//
// Created by phil on 4/15/21.
//

#include <cassert>
#include <gurobi_c++.h>


#define SILENCE(code) FILE* tmp = stdout;stdout = tmpfile(); code fclose(stdout);stdout = tmp;

using KeyPattern = std::vector<std::vector<uint64_t >>;
using KeyPatternList = std::vector<KeyPattern>;


void modelCopy(GRBModel &model, uint64_t n, GRBVar *input, GRBVar *output1, GRBVar *output2){
    for(uint64_t i = 0; i < n; i++){
        GRBVar temp[2] = {output1[i], output2[i]};
        model.addGenConstrOr(input[i], temp, 2);
    }
    model.update();
}


void modelShift(GRBModel &model, uint64_t n, GRBVar *input, GRBVar *output, uint64_t shift){
    for(uint64_t i = 0; i < n; i++){
        model.addConstr(input[i] == output[(i-shift+n) % n]);
    }
    model.update();
}


void modelAnd(GRBModel &model, uint64_t n, GRBVar *input1, GRBVar *input2, GRBVar *output){
    for(uint64_t i = 0; i < n; i++){
        model.addConstr(output[i] == input1[i]);
        model.addConstr(output[i] == input2[i]);
    }
    model.update();
}


void setOutputBit(GRBModel &model, uint64_t n, uint64_t rounds, uint64_t outputBit){
    for(uint64_t i = 0; i < 2*n; i++){
        std::string yStr;
        if(i < n){
            yStr = "Y2_" + std::to_string(rounds-1) + "_" + std::to_string(i);
        }
        else{
            yStr = "Y1_" + std::to_string(rounds-1) + "_" + std::to_string(i-n);
        }


        if(i == outputBit){
            model.addConstr(model.getVarByName(yStr) == 1);
        }
        else{
            model.addConstr(model.getVarByName(yStr) == 0);
        }
    }
}


GRBModel createSimonModel(GRBEnv &env, uint64_t n, uint64_t rounds, uint64_t shift1, uint64_t shift2,
                          uint64_t shift3, std::vector<std::vector<GRBVar>> &roundKeyVector){

    GRBModel model(env);

    // input and output variables
    GRBVar X1[rounds][n], X2[rounds][n], Y1[rounds][n], Y2[rounds][n];
    // intermediate variables
    GRBVar I[rounds][9][n];
    // inner round keys
    GRBVar RoundKey[rounds][n];

    KeyPatternList keySchedule;


    // Init all variables
    for(uint64_t r = 0; r < rounds; r++){
        for(uint64_t i = 0; i < n; i++){
            std::string str = std::to_string(r) + "_" + std::to_string(i);
            X1[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "X1_" + str);
            X2[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "X2_" + str);
            Y1[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "Y1_" + str);
            Y2[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "Y2_" + str);
            RoundKey[r][i] = model.addVar(0, 1, 0, GRB_BINARY, "K_" + str);

            if(!roundKeyVector.empty()){
                roundKeyVector[r][i] = RoundKey[r][i];
            }

            for(uint64_t j = 0; j < 9; j++){
                std::string strI = "I_" + std::to_string(r) + "_" + std::to_string(j) + "_" + std::to_string(i);
                I[r][j][i] = model.addVar(0, 1, 0, GRB_BINARY, strI);
            }
        }
    }

    // Model the rounds
    for(uint64_t r = 0; r < rounds; r++){

        // Connect rounds
        if(r > 0){
            for(uint64_t i = 0; i < n; i++){
                model.addConstr(X1[r][i] == Y2[r-1][i]);
                model.addConstr(X2[r][i] == Y1[r-1][i]);
            }
        }
        model.update();

        modelCopy(model, n, X1[r], I[r][0], I[r][1]);
        modelShift(model, n, I[r][0], I[r][2], shift1);

        modelCopy(model, n, I[r][1], I[r][3], I[r][4]);
        modelShift(model, n, I[r][3], I[r][5], shift2);

        modelAnd(model, n, I[r][2], I[r][5], I[r][6]);

        modelCopy(model, n, I[r][4], I[r][7], Y1[r]);
        modelShift(model, n, I[r][7], I[r][8], shift3);

        // Model the big XOR
        for(uint64_t i = 0; i < n; i++){
            model.addConstr(Y2[r][i] == RoundKey[r][i] + X2[r][i] + I[r][6][i] + I[r][8][i]);
        }
        model.update();
    }

    return model;
}


int64_t upperBoundAlgebraicDegree(uint64_t n, uint64_t rounds, uint64_t shift1, uint64_t shift2,
                                  uint64_t shift3, uint64_t outputBit, double timeLimit){

    int64_t upperBound = 2* int64_t(n) - 1;

    SILENCE(
    while (upperBound > 0) {
        int exitState = -1;

        std::vector<std::vector<GRBVar>> roundKeyVector;

        GRBEnv env;
        env.set(GRB_IntParam_OutputFlag, 0);
        env.set(GRB_IntParam_Threads, 1);
        env.set(GRB_DoubleParam_TimeLimit, timeLimit);

        GRBLinExpr inputSum;
        GRBLinExpr outputSum;


        auto simonModel = createSimonModel(env, n, rounds, shift1, shift2, shift3, roundKeyVector);

        for (uint64_t i = 0; i < n; i++) {
            std::string x1Str = "X1_0_" + std::to_string(i);
            std::string x2Str = "X2_0_" + std::to_string(i);

            inputSum += simonModel.getVarByName(x1Str);
            inputSum += simonModel.getVarByName(x2Str);

            std::string y1Str = "Y1_" + std::to_string(rounds - 1) + "_" + std::to_string(i);
            std::string y2Str = "Y2_" + std::to_string(rounds - 1) + "_" + std::to_string(i);

            outputSum += simonModel.getVarByName(y1Str);
            outputSum += simonModel.getVarByName(y2Str);
        }

        //simonModel.addConstr(outputSum == 1);
        simonModel.addConstr(inputSum == upperBound);
        setOutputBit(simonModel, n, rounds, outputBit);

        simonModel.optimize();

        exitState = simonModel.get(GRB_IntAttr_Status);

        if(exitState == 2){
            return upperBound;
        }
        else if(exitState == 9) {
            std::cout << "Timelimit reached" << std::endl;
            return upperBound;
        }

        assert(exitState == 3);
        upperBound--;

    })
    return -1;
}


void upperBounds(uint64_t n, uint64_t shift1, uint64_t shift2, uint64_t shift3,
                 double timeLimit){
    uint64_t rounds = 1;
    int64_t minDegree = 0;

    std::vector<std::pair<int64_t, int64_t>> boundsList;

    while(minDegree < int64_t (2*n) - 1){
        int64_t degree1 = upperBoundAlgebraicDegree(n, rounds, shift1, shift2, shift3, 0, timeLimit);
        int64_t degree2 = upperBoundAlgebraicDegree(n, rounds, shift1, shift2, shift3, n, timeLimit);

        minDegree = std::min(degree1, degree2);
        int64_t algDegree = std::max(degree1, degree2);

        std::cout << "Round: " << rounds << std::endl;
        std::cout << "MinDeg: " << minDegree;
        std::cout << "\tAlgDeg: " << algDegree << std::endl;

        boundsList.emplace_back(minDegree, algDegree);
        rounds++;
    }

    uint64_t counter = 1;
    std::cout << "MinDegree: " << std::endl << "----" << std::endl;
    for(auto const &boundPair : boundsList){
        std::cout << "(" << counter++ << ", " << boundPair.first << ")" << std::endl;
    }
    std::cout << "-----" << std::endl;

    counter = 1;
    std::cout << "Algebraic Degree: " << std::endl << "----" << std::endl;
    for(auto const &boundPair : boundsList){
        std::cout << "(" << counter++ << ", " << boundPair.second << ")" << std::endl;
    }
}

int main(int argc, char **argv){

    std::string usage = "Usage: " + std::string(argv[0]) + " n ciper timelimit\n";

    if (argc != 4){
        std::cout << usage;
        return -1;
    }

    uint64_t n = std::atoll(argv[1]);

    uint64_t shift1, shift2, shift3;

    if(std::string(argv[2]) == "simon"){
        shift1 = 1;
        shift2 = 8;
        shift3 = 2;
    }
    else if(std::string(argv[2]) == "simeck"){
        shift1 = 0;
        shift2 = 5;
        shift3 = 1;
    }
    else{
        std::cerr << "Invalid cipher" << std::endl;
        std::cout << usage;
        return -1;
    }

    double timeLimit = std::atof(argv[3]);

    std::cout << "Cipher " << std::string(argv[2]) << std::endl;
    std::cout << "n " << n << std::endl;

    upperBounds(n, shift1, shift2, shift3, timeLimit);

    return 0;
}