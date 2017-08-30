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
		Tissue* my_tissue;
		int rank;
		//keep track of when spawned
		int init_cell_time;
		//approx location
		Coord cell_center;
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
		Cell(int Ti, Tissue* tiss);
		Cell(int rank, Coord corner, double height, double width, 
			 int Ti, Tissue* tiss);
		// Destructors
		~Cell();

		// Getters and Setters
		Coord get_Cell_Center();
		void update_Cell_Center();
		void set_Wall_Cnt(int cnt);
		void set_Rank(const int id);
		void get_CytNodes(vector<Cyt_Node*>& cyts);
		void set_CytNodes(vector<Cyt_Node*>& cyts);
		void get_CornerNodes(vector<Wall_Node*>& corns);
		void set_CornerNodes(vector<Wall_Node*>& corns);
		Wall_Node* get_WallNodes();
		void get_Neighbor_Cells(vector<Cell*>& cells);
		int get_Num_Nodes();
		bool get_Reasonable_Bounds(Wall_Node* curr, Wall_Node* & A, Wall_Node* & B);

		// Calc Forces
		void calc_New_Forces();
		void update_Node_Locations();
		void update_Wall_Angles();
		void update_Neighbor_Cells();
		
		// Output current frame of simulation after update locations
		void print_Data_Output(ofstream& ofs);
		void print_VTK_Points(ofstream& ofs, int& count);
		void print_VTK_Scalars(ofstream& ofs);
		void print_VTK_Vectors(ofstream& ofs);
		
		// Growth of cell
		Wall_Node* find_Largest_Length(int& side);
		void add_Wall_Node(const int Ti);
		void add_Cyt_Node(const int Ti);

		// Division
		Cell* divide(const int Ti);
		Cell* divide_length_wise(const int Ti);
		Cell* divide_width_wise(const int Ti);

};


// End Cell Class
//===================

#endif





