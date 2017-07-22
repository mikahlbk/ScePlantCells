//cell.h
//===================
// Inlcude Guards
#ifndef _CELL_H_INCLUDED_
#define _CELL_H_INCLUDED_
//===================
// forward declarations

//===================
// include dependencies
#include "coord.h"
#include "node.h"
//===================
// Cell Class Declaration

class Cell {

	private:
		//keep track of how many nodes the cell has created
		int num_wall_nodes;
		int num_cyt_nodes;
		vector<int> wall_nodes_per_frame;
		vector<int> cyt_nodes_per_frame;
		//2D vector for node locations
		vector< vector<Coord> > wall_locs;
		vector< vector<Coord> > cyt_locs;
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
		// Talking to other Cells

		// Growth of cell
		// returns the coordinte of the wall node whose 
		// left length is the largest
		Wall_Node* find_Largest_Length();
		void add_Cell_Wall_Node();
		void add_Cyt_Node();

};


// End Cell Class
//===================

#endif
