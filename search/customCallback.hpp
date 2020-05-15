#ifndef H_CUSTOMCALLBACK
#define H_CUSTOMCALLBACK

#include <vector>
#include <iostream>
#include <cmath>
#include <cstdint>

#include <boost/dynamic_bitset.hpp>

#include "gurobi_c++.h"
#include "BCData.hpp"

/*
  Definition of the different callbacks used in the implementation
  Callbacks are function that can be called during the solving process of Gurobi (see its documentation for more details)
  This has a huge role in the efficiency of the algorithms in general, as it can, for example, 
  discard a solution found and look for another one without restarting the search from scratch
*/

class BCData;

//Callback class for optimization in the searchKeyPattern function
class callbackSearchKeyPattern: public GRBCallback
{
  public:

    BCData const & BCD;                                             //Some data about the block cipher
    unsigned int nbRounds;                                          //number of rounds convered by the model
    std::vector<GRBVar> inputVar;                                   //variables for the input
    std::vector<GRBVar> outputVar;                                  //variables for the output
    std::vector<std::vector<GRBVar>> allKeyVar;                     //all key variables
    std::vector<std::vector<std::vector<GRBVar>>> inSSBVar;         //input variables for the SSB
    std::vector<std::vector<std::vector<GRBVar>>> outSSBVar;        //output variables for the SSB
    std::vector<std::vector<std::vector<GRBVar>>> keySSBVar;        //key variables for the SSB
    uint64_t ctrsolutions;                                          //counter for the solutions found
    std::vector<std::vector<boost::dynamic_bitset<>>> allSolutions; //vector to keep all examined solutions

    callbackSearchKeyPattern(BCData const & xBCD,
                             unsigned int const xnbRounds,
                             std::vector<GRBVar> const & xinputVar,
                             std::vector<GRBVar> const & xoutputVar,
                             std::vector<std::vector<GRBVar>> const & xallKeyVar,
                             std::vector<std::vector<std::vector<GRBVar>>> const & xinSSBVar,
                             std::vector<std::vector<std::vector<GRBVar>>> const & xoutSSBVar,
                             std::vector<std::vector<std::vector<GRBVar>>> const & keySSBVar);

  protected:
    void callback();
};



//callback class when just counting solution, to see the progress
//Also has allows to set a limit on the number of solutions, but using the SolutionLimit parameter in the model is better
class callbackCount: public GRBCallback
{
  public:
    uint64_t ctr;
    uint64_t countLimit; 
    //Upper bound on the the number of solution to count
    //abort if it goes higher, don't limit if equal to 0
    bool wasAborted;
    bool printNbSol;

    callbackCount();

    callbackCount(uint64_t const xcountLimit,
                  bool xprintNbSol = false);

  protected:
    void callback();
};



//callback class for the improved dynamic search
class callbackDynamic: public GRBCallback
{
  public:
    BCData const & BCD;                                             //Some data about the block cipher
    uint rMiddle;
    std::vector<uint8_t> const & output;
    std::vector<std::vector<uint8_t>> const & keyVal;
    std::vector<GRBVar> const & inputMiddleVar;
    std::vector<GRBVar> const & keyMiddleVar;
    std::vector<std::vector<std::vector<GRBVar>>> const & inSSBVar;         //input variables for the SSB
    std::vector<std::vector<std::vector<GRBVar>>> const & outSSBVar;        //output variables for the SSB
    std::vector<std::vector<std::vector<GRBVar>>> const & keySSBVar;        //key variables for the SSB

    callbackDynamic(BCData const & xBCD,
                    uint const xrMiddle,
                    std::vector<uint8_t> const & xoutput,
                    std::vector<std::vector<uint8_t>> const & xkeyVal,
                    std::vector<GRBVar> const & xinputMiddleVar,
                    std::vector<GRBVar> const & xkeyMiddleVar,
                    std::vector<std::vector<std::vector<GRBVar>>> const & xinSSBVar,
                    std::vector<std::vector<std::vector<GRBVar>>> const & xoutSSBVar,
                    std::vector<std::vector<std::vector<GRBVar>>> const & xkeySSBVar);

  protected:
    void callback();
};

//callback class for the lower bound function
class callbackLowerBound: public GRBCallback
{
  public:
    BCData const & BCD;
    uint nbRounds;
    std::vector<GRBVar> const & inputMiddleVar;
    std::vector<uint8_t> const & output;
    std::vector<std::vector<uint8_t>> const & keyVal;

    callbackLowerBound(BCData const & xBCD,
                       uint const xnbRounds,
                       std::vector<GRBVar> const & xinputMiddleVar,
                       std::vector<uint8_t> const & xoutput,
                       std::vector<std::vector<uint8_t>> const & xkeyVal);

  protected:
    void callback();
};


#endif