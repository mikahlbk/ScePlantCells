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
		// Keep info about its cyt and wall nodes
		vector<Cyt_Node*> cyt_nodes;
		Wall_Node* first_corner;
		// Keep info about other cells
		

	public:
		// Constructors
		Cell();
		Cell(string filename);

		// Getters and Setters
		void get_CytNodes(vector<Cyt_Node*>& cyts);
		Wall_Node* get_WallNodes();

		// Keep track of neighbor cells

		// Calc Forces
		void calc_New_Forces();
		// Update Node Locations
		void update_Node_Locations();
		// Talking to other Cells

		// Growth of cell


};


// End Cell Class
//===================

#endif
