#include "SmallMatrix.hpp"

using namespace std;
typedef unsigned int uint;

SmallMatrix::SmallMatrix() :
	rows(),
	nrows(0),
	ncols(0)
	//Default constructor, build a 0x0 matrix
{}

SmallMatrix::SmallMatrix(uint i, uint j) :
	rows(i,0),
	nrows(i),
	ncols(j)
	//Build a i x j zero matrix
{}

SmallMatrix::SmallMatrix(vector<uint64_t> const & xrows,
						 uint xncols) :
						 rows(xrows),
						 nrows(xrows.size()),
						 ncols(xncols)						 
	//Build a matrix with rows provided by #xrows and with #xncols columns
{}

uint64_t& SmallMatrix::get(uint i){
	//Get the #i-th row
	return rows[i];
}

uint64_t const & SmallMatrix::get(uint i) const {
	return rows[i];
}

uint SmallMatrix::get(uint const i,
		 			  uint const j) const{
	//Get coefficient [i][j] of the matrix
	if(rows[i] & (uint64_t(1) << j)) return 1;
	else return 0;
}

void SmallMatrix::set(uint const i,
					  uint const j,
					  uint const x){
	//Set coefficient [#i][#j] of the matrix to #x

	if(x == 0)
		rows[i] = rows[i] & (uint64_t(-1) - (uint64_t(1) << j));
	else
		rows[i] = rows[i] | (uint64_t(1) << j);
}

void SmallMatrix::addRow(uint64_t const r){
	//Add the row #r to the matrix
	rows.emplace_back(r);
	nrows++;
}

void checkLinIndep(SmallMatrix & M, 
				   uint64_t newRow,
				   map<uint,uint> & listPivot){
	/*
		Check if #newRow is linearly independent from the rows of #M
		If so, add it to M and echelonize-reduce
		- #M is the original matrix
		- #newRow is the new row to check
		- #listPivot is a map so that listPivot[i] = j, meaning that the pivot for column i is row j of M
	*/

	if(M.nrows == 0){
		uint indexNewPivot = 0;
		while(!(newRow & (uint64_t(1) << indexNewPivot))) indexNewPivot++;
		listPivot[indexNewPivot] = M.nrows;
		M.addRow(newRow);
	}
	else{
		//For each known pivot, check if we need to reduce the new row with it
		for(auto & pivot : listPivot){
			uint indexPivot = pivot.first;
			if(newRow & (uint64_t(1) << indexPivot)){
				//We need to reduce on this pivot
				newRow ^= M.get(pivot.second);
			}
		}

		if(newRow != 0){
			//The new row was not fully cancelled by the pivots, hence it's lin indep
			//Find the column for the new pivot
			uint indexNewPivot = 0;
			while(!(newRow & (uint64_t(1) << indexNewPivot))) indexNewPivot++;
			//Add the pivot to the listof known pivots
			listPivot[indexNewPivot] = M.nrows;
			//Reduce M with this new pivot
			for(uint i = 0; i < M.nrows; i++){
				auto & Mi = M.get(i);
				if(Mi & (uint64_t(1) << indexNewPivot)) Mi ^= newRow;
			}
			//Add the new row to M
			M.addRow(newRow);
		}
	}

}

void SmallMatrix::print(){
	//Bin print
	for(auto const & r : rows){
		for(uint i = 0; i < ncols; i++){
			if(r & (uint64_t(1) << i)) cout << 1;
			else cout << 0;
		}
		cout << endl;
	}
}