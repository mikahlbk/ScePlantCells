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
		int num_cells;
	public:
		Tissue(string filename);
		
		void get_Cells(vector<Cell*>& cells);
		void update_Life_Length();
		void update_Wall();
		void update_Cytoplasm();
		void calc_New_Forces();
		void update_Cell_Locations();
		void update_Neighbor_Cells();
		//void cell_Division(const int Ti);
		void print_Data_Output(ofstream& ofs);
		void print_VTK_File(ofstream& ofs);
		//void make_Vectors();
		void update_Adhesion();
		int get_Num_Cells() {return num_cells;}
		//void grow_Cells(const int Ti);
		//Destructor
		~Tissue();
};


//===========================
//End of file

#endif
