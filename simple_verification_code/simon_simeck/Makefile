CXX = g++
CXXFLAGS = -O3 -march=native


all: upperBounds verifyDegree

upperBounds: upperBounds.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -L /opt/gurobi911/linux64/lib -lntl -lcryptominisat5 -lgurobi91 -lgurobi_g++5.2

verifyDegree: verifyDegree.cpp
	$(CXX) $(CXXFLAGS) -o $@ $< -lntl -lcryptominisat5
