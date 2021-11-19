//
// Created by phil on 5/31/21.
//

#include <fstream>
#include <vector>
#include <cassert>
#include <sstream>
#include <cryptominisat5/cryptominisat.h>
#include <NTL/mat_GF2.h>

#define ASSERT(assumption, msg) if(!(assumption)){std::cerr << "ERROR" << std::endl << msg << std::endl;}

using Monomial = std::vector<uint64_t>;
using KeyPattern = std::vector<std::vector<uint64_t >>;
using KeyInputList = std::vector<std::pair<Monomial, KeyPattern>>;
using KeyInputListList = std::vector<KeyInputList>;

using TrailRow = std::vector<int64_t>;
using TrailMatrix = std::vector<TrailRow>;


// Parse input file with the key pattern
uint64_t parseFile(const std::string &fileName, const uint64_t n, KeyInputListList &keyInputListList){

    std::ifstream file(fileName, std::ios::in);
    uint64_t rounds = 0;

    if(file.is_open()){
        std::string line;

        std::getline(file, line);
        ASSERT(line == "xxxx", "read file error 1");

        while(std::getline(file, line)){
            if(line.empty()){
                file.close();
                return 0;
            }
            else{
                KeyInputList keyInputList;
                rounds++;

                do {
                    std::vector<uint64_t> input;
                    
                    for (uint64_t i = 0; i < 2 * n; i++) {
                        input.push_back(std::stoull(std::string(1, line[i])));
                    }

                    KeyPattern keyPattern;

                    for (uint64_t i = 0; i < rounds; i++) {
                        std::vector<uint64_t> roundKeyPattern;
                        std::getline(file, line);
                        for (uint64_t k = 0; k < n; k++) {
                            roundKeyPattern.push_back(std::stoull(std::string(1, line[k])));
                        }
                        keyPattern.push_back(roundKeyPattern);
                    }

                    keyInputList.emplace_back(input, keyPattern);
                    std::getline(file, line);
                    ASSERT(line == "----", "read file error 2");
                    std::getline(file, line);

                } while(line != "xxxx" && !line.empty());

                keyInputListList.push_back(keyInputList);

                if(line.empty()){
                    file.close();
                    return 0;
                }
            }
        }

    }
    else {
        std::cerr << "Error opening file!" << std::endl;
        return -1;
    }
    return 0;
}


uint64_t hammingWeight(const std::vector<uint64_t> &monomial){
    uint64_t result = 0;
    for(auto elem: monomial){
        result += elem;
    }
    return result;
}


// Model COPY for (a) -> (b_1, b_2)
inline std::string copySAT(uint64_t a, uint64_t b1, uint64_t b2){
    std::string satStr = "-" + std::to_string(a) + " " + std::to_string(b1) + " " + std::to_string(b2) + " 0\n";
    satStr += std::to_string(a) + " -" + std::to_string(b1) + " 0\n";
    satStr += std::to_string(a) + " -" + std::to_string(b2) + " 0\n";

    return satStr;
}


// Model XOR for (a_1, a_2) -> (b)
inline std::string xorSAT(uint64_t a1, uint64_t a2, uint64_t b){
    std::string satStr = std::to_string(a1) + " " + std::to_string(a2) + " -" + std::to_string(b) + " 0\n";
    satStr += std::to_string(b) + " -" + std::to_string(a1) + " 0\n";
    satStr += std::to_string(b) + " -" + std::to_string(a2) + " 0\n";
    satStr += "-" + std::to_string(a1) + " -" + std::to_string(a2) + " 0\n";

    return satStr;
}


// Model AND for (a_1, a_2) -> (b)
inline std::string andSAT(uint64_t a1, uint64_t a2, uint64_t b){
    std::string satStr = std::to_string(b) + " -" + std::to_string(a1) + " 0\n";
    satStr += std::to_string(b) + " -" + std::to_string(a2) + " 0\n";
    satStr += std::to_string(a1) + " -" + std::to_string(b) + " 0\n";
    satStr += std::to_string(a2) + " -" + std::to_string(b) + " 0\n";

    return satStr;
}


std::string generateSimonSatFormula(uint64_t n, uint64_t rounds, uint64_t shift1, uint64_t shift2, uint64_t shift3,
                                    const KeyPattern &keyPattern,
                                    const std::vector<uint64_t> &inputMonomial,
                                    uint64_t outputBit){

    uint64_t nVariables = 11*n*rounds + 2*n;

    // Clauses for set input, output and key
    uint64_t nClauses = 4*n + n*rounds;
    // Clauses for COPY
    nClauses += 3 * 3 * n * rounds;
    // Clauses for AND
    nClauses += 4 * n * rounds;
    // Clauses for XOR
    nClauses += 3 * 4 * n * rounds;

    // Round offset
    const uint64_t R_OFF = 11 * n;

    std::string cnfString = "p cnf " + std::to_string(nVariables) + " " + std::to_string(nClauses) + "\n";

    // Set input monomial
    for(uint64_t i = 0; i < 2*n; i++){
        if(inputMonomial[i] == 0){
            cnfString += "-" + std::to_string(i+1) + " 0\n";
        }
        else{
            assert(inputMonomial[i] == 1);
            cnfString += std::to_string(i+1) + " 0\n";
        }
    }


    // Propagate the rounds
    for(uint64_t r = 0; r < rounds; r++){
        // First copy
        for(uint64_t i = 1; i <= n; i++){
            uint64_t shiftedPos = r*R_OFF + 2*n + ( (i - shift1 - 1 + n) % n) + 1;
            cnfString += copySAT(r*R_OFF + i, shiftedPos, r*R_OFF + 3*n + i);
        }
        // Second copy
        for(uint64_t i = 1; i <= n; i++){
            uint64_t shiftedPos = r*R_OFF + 4*n + ( (i - shift2 - 1 + n) % n) + 1;
            cnfString += copySAT(r*R_OFF + 3*n + i, shiftedPos, r*R_OFF + 5*n + i);
        }
        // And
        for(uint64_t i = 1; i <= n; i++){
            cnfString += andSAT(r*R_OFF + 2*n + i, r*R_OFF + 4*n + i, r*R_OFF + 6*n + i);
        }
        // First Xor
        for(uint64_t i = 1; i <= n; i++){
            cnfString += xorSAT(r*R_OFF + n + i, r*R_OFF + 6*n + i, r*R_OFF + 7*n + i);
        }
        // Third copy
        for(uint64_t i = 1; i <= n; i++){
            uint64_t shiftedPos = r*R_OFF + 8*n + ( (i - shift3 - 1 + n) % n) + 1;
            cnfString += copySAT(r*R_OFF + 5*n + i, shiftedPos, (r+1)*R_OFF + n + i);
        }
        // Second Xor
        for(uint64_t i = 1; i<= n; i++){
            cnfString += xorSAT(r*R_OFF + 7*n + i, r*R_OFF + 8*n + i, r*R_OFF + 9*n + i);
        }
        // Set Key Pattern
        for(uint64_t i = 0; i < n; i++){
            if(keyPattern[r][i] == 0){
                cnfString += "-" + std::to_string(r*R_OFF + 10*n + i + 1) + " 0\n";
            }
            else{
                assert(keyPattern[r][i] == 1);
                cnfString += std::to_string(r*R_OFF + 10*n + i + 1) + " 0\n";
            }
        }
        // Third Xor
        for(uint64_t i = 1; i<= n; i++){
            cnfString += xorSAT(r*R_OFF + 9*n + i, r*R_OFF + 10*n + i, (r+1)*R_OFF + i);
        }
    }

    // Set output bit
    assert(outputBit < 2*n);
    for(uint64_t i = 0; i < 2*n; i++){
        if(i == outputBit){
            cnfString += std::to_string(rounds * R_OFF + i + 1) + " 0\n";
        }
        else{
            cnfString += "-" + std::to_string(rounds * R_OFF + i + 1) + " 0\n";
        }
    }

    return cnfString;
}


int64_t count_solutions(CMSat::SATSolver &sat_solver, bool reduced, int timeLimit){
    uint64_t solutions = 0;
    CMSat::lbool ret = sat_solver.solve();
    std::vector<CMSat::Lit> ban_solution;

    std::time_t startTime = std::time(nullptr);

    while(ret == CMSat::l_True){
        ban_solution.clear();

        if(std::time(nullptr) - startTime > timeLimit){
            return -1;
        }

        for (uint64_t var = 0; var < sat_solver.nVars(); var++) {
            if (sat_solver.get_model()[var] != CMSat::l_Undef) {
                ban_solution.emplace_back(var, (sat_solver.get_model()[var] == CMSat::l_True )) ;
            }
            else{
                std::cout << "Error " << std::endl;
            }
        }

        sat_solver.add_clause(ban_solution);
        ret = sat_solver.solve();
        solutions++;
    }
    if(reduced) {
        return solutions % 2;
    }
    else{
        return solutions;
    }
}


CMSat::SATSolver cnfToSatInstance(const std::string &cnfStr){
    CMSat::SATSolver solver;

    std::string line;
    std::stringstream sstream(cnfStr);

    //std::getline(sstream, line, '\n');
    std::getline(sstream, line, '\n');

    line = line.substr(6);
    std::size_t pos = line.find(' ');

    uint64_t noClauses = std::stoll(line.substr(pos+1));
    uint64_t noLiterals = std::stoll(line.substr(0, pos));

    uint64_t iterations = 0;

    solver.new_vars(noLiterals);

    while(std::getline(sstream, line, '\n')){

        std::istringstream istream(line);
        int64_t literal;

        std::vector<CMSat::Lit> clause;
        while( istream >> literal ) {
            if(literal != 0){

                if(literal > 0){
                    clause.emplace_back(literal-1, false);
                }
                else{
                    clause.emplace_back(-literal-1, true);
                }
            }
        }

        solver.add_clause(clause);
        iterations++;
    }

    assert(iterations == noClauses);

    return solver;
}


int64_t countTrailsSat(uint64_t n, uint64_t shift1, uint64_t shift2, uint64_t shift3, uint64_t rounds,
                       uint64_t outputBit, const std::vector<uint64_t> &inputMonomial,
                       const KeyPattern &keyPattern, bool reduced, int timeLimit){

    auto cnfStr = generateSimonSatFormula(n, rounds, shift1, shift2, shift3, keyPattern, inputMonomial, outputBit);
    auto solver = cnfToSatInstance(cnfStr);
    return count_solutions(solver, reduced, timeLimit);
}


uint64_t computeRank(const uint64_t n, const TrailMatrix &trailMatrix){
    uint64_t rank = 0, oldRank;
    std::vector<bool> unitVectors(2*n, false);
    NTL::vec_vec_GF2 matrix;

    do{
        oldRank = rank;

        for(const auto &rowRef : trailMatrix){
            auto row(rowRef);

            assert(row.size() == 2*n);
            bool usable = true;
            NTL::vec_GF2 rowVector;
            rowVector.SetLength(long(2*n));

            for(int64_t pos = 0; pos < int64_t (2*n); pos++){

                if(row[pos] == -1){
                    if(unitVectors[pos]){
                        row[pos] = 0;
                    }
                    else{
                        usable = false;
                        break;
                    }

                }
                rowVector[pos] = row[pos];
            }

            if(usable){
                matrix.append(rowVector);
            }

        }

        NTL::mat_GF2 gaussMatrix;
        NTL::conv(gaussMatrix, matrix);
        rank = NTL::gauss(gaussMatrix);

        for(int64_t row = 0; row < gaussMatrix.NumRows(); row++){
            if( NTL::weight(gaussMatrix[row]) == 1){
                matrix.append(gaussMatrix[row]);

                for(int64_t pos = 0; pos < gaussMatrix[row].length(); pos++){
                    if(gaussMatrix[row][pos] == 1){
                        unitVectors[pos] = true;
                    }
                }
            }
        }


    } while (rank != oldRank);

    return rank;
}


uint64_t rotatePosition(uint64_t n, uint64_t position, uint64_t shift){
    uint64_t half = position / n;

    if(half == 0){
        return (position+shift) % n;
    }
    else if(half == 1){
        return ((position+shift) % n) + n;
    }
    else{
        std::cerr << "rotate position error" << std::endl;
        return UINT64_MAX;
    }
}


TrailRow rotateRow(const uint64_t n, const TrailRow &trailRow, const uint64_t shift){
    TrailRow rotatedRow(2*n, -3);

    for(uint64_t i = 0; i < 2*n; i++){
        rotatedRow[i] = trailRow[rotatePosition(n, i, shift)];
    }

    return rotatedRow;
}


void computeDegree(const uint64_t n, const uint64_t shift1, const uint64_t shift2,
                   const uint64_t shift3, const KeyInputListList &keyInputListList,
                   const int64_t timeLimit){

    for(uint64_t i = 0; i < keyInputListList.size(); i++){
        uint64_t rounds = i+1;

        KeyInputList keyInputList = keyInputListList[i];

        uint64_t minDegree = UINT64_MAX;
        uint64_t algebraicDegree = 0;
        TrailMatrix trailMatrix;

        for(const auto& [input, keyPattern] : keyInputList){
            minDegree = std::min(minDegree, hammingWeight(input));
            algebraicDegree = std::max(algebraicDegree, hammingWeight(input));

            TrailRow trailRow;

            for(uint64_t outputBit = 0; outputBit < 2*n; outputBit++){
                auto trails = countTrailsSat(n, shift1, shift2, shift3, rounds, outputBit, input,
                                             keyPattern, false, int(timeLimit));
                trailRow.push_back(trails);
            }

            for(uint64_t rotation = 0; rotation < n; rotation++){
                trailMatrix.push_back(rotateRow(n, trailRow, rotation));
            }

        }

        std::cout << "Rounds: " << rounds << std::endl;
        std::cout << "MinDegree: " << minDegree << std::endl;
        std::cout << "Rank: " << computeRank(n, trailMatrix) << std::endl;
        std::cout << "Algebraic degree: " << algebraicDegree << std::endl;
        std::cout << "--------------------------------------" << std::endl;
    }
}


int main(int argc, char **argv){
    if(argc != 5){
        std::cerr << "Usage: " << argv[0] << " cipher n input timelimit" << std::endl;
        return -1;
    }

    std::string inputFile = std::string(argv[3]);
    const uint64_t n = atoll(argv[2]);
    const uint64_t timeLimit = atoll(argv[4]);

    uint64_t shift1, shift2, shift3;

    if(std::string(argv[1]) == "simon"){
        shift1 = 1;
        shift2 = 8;
        shift3 = 2;
    }
    else if(std::string(argv[1]) == "simeck"){
        shift1 = 0;
        shift2 = 5;
        shift3 = 1;
    }
    else{
        std::cerr << "Invalid cipher" << std::endl;
        return -1;
    }

    std::cout << "Verification program degree" << std::endl;
    std::cout << "Cipher: " << argv[1] << std::endl;
    std::cout << "Time limit: " << timeLimit << " seconds" << std::endl << std::endl;

    KeyInputListList keyInputListList;

    parseFile(inputFile, n, keyInputListList);
    computeDegree(n, shift1, shift2, shift3, keyInputListList, timeLimit);

    return 0;
}