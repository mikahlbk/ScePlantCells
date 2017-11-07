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
		int layer;
		int growth_rate;
		Coord cell_center;
		vector<Side*> sides;
		vector<Cyt_Node*> cyt_nodes;
		vector<Cell*> neigh_cells;
		vector<double>lengths;
	public:
		// Constructors
		Cell(int rank, Tissue* tissue);
		Cell(int rank, Coord corner, double height, double width, 
			 int Ti, Tissue* tiss, int layer);

		// Destructor
		~Cell();

		// Getters and Setters
		Coord get_Cell_Center();
		void get_Cyt_Nodes(vector<Cyt_Node*>& cyts);
		Wall_Node* get_Wall_Nodes();
		void get_Sides(vector<Side*>& sides);
		void get_Neighbor_Cells(vector<Cell*>& cells);
		int get_Node_Count();
		int get_Rank();
		void set_Rank(const int id);
		void set_Sides(vector<Side*>& sides);
		void stretching(int Ti, Coord force);
		int get_Layer(){return layer;}
		void set_Layer(int layer);
		void set_growth_rate(int growth_rate);
		int get_growth_rate(){return growth_rate;}
		void get_Lengths (vector<double>& lengths);
		// Keep track of neighbor cells
		void update_Neighbor_Cells();
		void update_adhesion_springs();
		double length();
	//	bool get_Reasonable_Bounds(Wall_Node* curr, Wall_Node* & A, Wall_Node* & B);

		// Forces and Positionsing
		void calc_New_Forces();
		void update_Node_Locations();
		void update_Wall_Angles();
		void update_Cell_Center();
		void update_Life_Length();

		//Output Functions
		void print_Data_Output(ofstream& ofs);
		int update_VTK_Indices(int& id);
		void print_VTK_Adh(ofstream& ofs);
		void print_VTK_Points(ofstream& ofs, int& count);
		void print_VTK_Scalars(ofstream& ofs);
		void print_VTK_Vectors(ofstream& ofs);
		
		// Growth of cell
		Wall_Node* find_Largest_Length();
		void add_Wall_Node();
		void add_Cyt_Node();
	
		//Division 
		Cell* divide(const int Ti);
		Cell* divide_length_wise(const int Ti);
		Cell* divide_width_wise(const int Ti);
		void add_Cyt_Node_Div();
		
};


// End Cell Class
//===================

#endif

