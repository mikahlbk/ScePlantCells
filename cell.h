//cell.h
//===================
// Inlcude Guards
#ifndef _CELL_H_INCLUDED_
#define _CELL_H_INCLUDED_
//===================
// forward declarations
class Tissue;
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
		int rank;
		Tissue* my_tissue;
		int life_length;
		int num_cyt_nodes;
		Coord cell_center;
		vector<Side*> sides;
		vector<Cyt_Node*> cyt_nodes;
		vector<Cell*> neigh_cells;

	public:
		// Constructors
		Cell(int rank, Tissue* tissue);
		Cell(int rank, Coord corner, double height, double width, 
			 int Ti, Tissue* tiss);

		// Destructor
		~Cell();

		// Getters and Setters
		Coord get_Cell_Center();
		void get_Cyt_Nodes(vector<Cyt_Node*>& cyts);
		Wall_Node* get_Wall_Nodes();
		void get_Sides();
		void get_Neighbor_Cells(vector<Cell*>& cells);
		int get_Node_Count();
		
		// Keep track of neighbor cells
		void update_Neighbor_Cells();
		bool get_Reasonable_Bounds(Wall_Node* curr, Wall_Node* & A, Wall_Node* & B);

		// Forces and Positionsing
		void calc_New_Forces();
		void update_Node_Locations();
		void update_Wall_Angles();
		void udpate_Cell_Center();
		void update_Life_Length();

		//Output Functions
		void print_Data_Output(ofstream& ofs);
		void print_VTK_Points(ofstream& ofs, int& count);
		void print_VTK_Scalars(ofstream& ofs);
		void print_VTK_Vectors(ofstream& ofs);
		
		// Growth of cell
		Wall_Node* find_Largest_Length();
		void add_Wall_Node();
		void add_Cyt_Node();

};


// End Cell Class
//===================

#endif

