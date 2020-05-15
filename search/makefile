CXX = g++
CXXFLAGS = -O3 -march=native -Wall -Wextra -std=c++11 -fopenmp -Wno-deprecated-copy #-g -fsanitize=address -fno-omit-frame-pointer

IFLAGS = -I $(GUROBI_HOME)/include/ -I /usr/include/x86_64-linux-gnu/
LFLAGS = -L $(GUROBI_HOME)/lib/ -lgurobi_g++5.2 -lgurobi81 -lgmpxx -lgmp

%.o: %.cpp
	$(CXX) $(IFLAGS) $(LFLAGS) $(CXXFLAGS) -c $< -o $@ 

main :  main.o BCData.o configBC.o customCallback.o SmallMatrix.o $(SHARP_OBJECTES) 
	$(CXX) $(CXXFLAGS) $(IFLAGS) -o main main.o BCData.o configBC.o customCallback.o SmallMatrix.o $(SHARP_OBJECTES) $(LFLAGS)

clean :
	rm -rf *.o