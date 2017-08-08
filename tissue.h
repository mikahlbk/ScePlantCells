//tissue.h
//=================
//Include Guards
#ifndef _TISSUE_H_INCLUDED_
#define _TISSUE_H_INCLUDED_
//=========================
//forward declarations

//=======================
//Include dependencies
#include <string>
#include <vector>
#include <fstream>

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
//=========================
// Tissue Class Declaration

class Tissue {

	private:
		// We'll need a better data structure for later
		vector<Cell*> cells;

	public:
		Tissue(string filename);
		
		void get_Cells(vector<Cell*>& cells);

		void calc_New_Forces();
		void update_Cell_Locations();
		void update_Neighbor_Cells();

		void print_Data_Output(ofstream& ofs);
		void print_VTK_File(ofstream& ofs);

		void grow_Cells(const int Ti);


};


//===========================
//End of file

#endif
