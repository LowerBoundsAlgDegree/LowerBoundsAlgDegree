#ifndef H_SMALLMATRIX
#define H_SMALLMATRIX

#include <vector>
#include <cstdint>
#include <iostream>
#include <map>

/*
	Small class to do some stuff with small binary matrices (max number of columns : 64)
*/

class SmallMatrix{

	public:

		std::vector<uint64_t> rows; //The row values
		unsigned int nrows;			//Number of rows
		unsigned int ncols;			//Number of columns
		

		SmallMatrix();
		//Default constructor, build a 0x0 matrix

		SmallMatrix(uint i, uint j);
		//Build a i x j zero matrix

		SmallMatrix(std::vector<uint64_t> const & xrows,
					unsigned int const xncols);
		//Build a matrix with rows provided by #xrows and with #xncols columns

		uint64_t& get(uint i);
		uint64_t const & get(uint i) const;
		//Get the #i-th row

		unsigned int get(uint const i,
				 		 uint const j) const;
		//Get coefficient [#i][#j] of the matrix

		void set(uint const i,
				 uint const j,
				 uint const x);
		//Set coefficient [#i][#j] of the matrix to #x

		void addRow(uint64_t const r);
		//Add the row #r to the matrix

		void print();
		//Bin print

};

void checkLinIndep(SmallMatrix & M, 
						   uint64_t const newRow,
						   std::map<uint,uint> & listPivot);
/*
	Check if #newRow is linearly independent from the rows of #M
	If so, add it to M and echelonize-reduce
	- #M is the original matrix
	- #newRow is the new row to check
	- #listPivot is a map so that listPivot[i] = j, meaning that the pivot for column i is row j of M
*/

#endif