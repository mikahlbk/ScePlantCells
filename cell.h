//cell.h
//===================
// Inlcude Guards
#ifndef _CELL_H_INCLUDED_
#define _CELL_H_INCLUDED_
//===================
// forward declarations

//===================
// include dependencies
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "phys.h"
#include "coord.h"
#include "node.h"
//===================
// Cell Class Declaration

class Cell {

	private:
		//keep track of how many nodes the cell has created
		int num_wall_nodes;
		int num_cyt_nodes;
		// Keep info about its cyt and wall nodes
		vector<Cyt_Node*> cyt_nodes;
		vector<Wall_Node*> corners;
		// Keep info about other cells
		vector<Cell*> neigh_cells;

	public:
		// Constructors
		Cell(string filename);

		// Getters and Setters
		void get_CytNodes(vector<Cyt_Node*>& cyts);
		Wall_Node* get_WallNodes();
		void get_Neigh_Cells(vector<Cell*>& cells);

		// Keep track of neighbor cells

		// Calc Forces
		void calc_New_Forces();
		// Update Node Locations
		void update_Node_Locations();
		// Update Angles
		void update_Wall_Angles();
		// Output current frame of simulation after update locations
		void print_Data_Output(ofstream& ofs);
		void print_VTK_File(ofstream& ofs);
		// Talking to other Cells
		
		// Growth of cell
		// returns the coordinte of the wall node whose 
		// left length is the largest
		Wall_Node* find_Largest_Length(int& side);
		void add_Wall_Node(const int Ti);
		void add_Cyt_Node(const int Ti);

};


// End Cell Class
//===================

#endif



