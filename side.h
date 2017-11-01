//side.h
//======================
//Include Guards
#ifndef _SIDE_H_INCLUDED_ 
#define _SIDE_H_INCLUDED_
//=======================
//Forward Declarations
class Cell;
//=======================
//Include Declarations
#include <iostream>
#include "phys.h"
#include "coord.h"
#include "node.h"
//======================

class Side {
	protected:
		Wall_Node* end_A;
		Wall_Node* end_Z;
		int num_wall_nodes;
		Cell* my_cell;
		vector<Side*> touching_Sides;
		int side_type;
		double linear_spring;
		double bending_spring;
		vector<double>lengths;
		vector<double>forces;
	public:
		//Constructors
		Side(Coord point_a, Coord point_z, Cell* my_cell, int num_nodes);
		void connect_Ends(Side* s);
		//Getters and Setters
		//cyt nodes
		void get_Cyt_Nodes(vector<Cyt_Node*>& cyts);
		//wall nodes
		Wall_Node* get_Wall_Nodes() {return end_A;}
		Wall_Node* get_End_A() {return end_A;}
		Wall_Node* get_End_Z() {return end_Z;}
		int get_Wall_Count() {return num_wall_nodes;}
		//cell level
		void update_Touching_Sides(vector<Cell*>& neighbor_Cells);
		void get_Touching_Sides(vector<Side*>& touching_Sides);
		int get_Side_Type() {return side_type;}
		//void update_Touching_Neighbors(vector<Cell*>& touching_Neighbors);
		Cell* get_My_Cell() {return my_cell;}
		void set_My_Cell(Cell* new_cell);
		//add setter for division to set new cell
		//parameters for forces
		void set_Phys_Parameters(double kbend, double klin);
		void set_Side_Type(int type);
		double get_Linear_Spring() {return linear_spring;}
		double get_Bending_Spring() {return bending_spring;}
		void get_Lengths(vector<double>& lengths);
		void get_Forces(vector<double>& forces);
		void stretching_Test(int Ti, Coord force); 
		double length();		
		//double force_Sum();
		//Force and Positioning
		

		//Cell Growth
		void add_Wall_Node(Wall_Node* right);

		//Destructors
		~Side();
};


#endif
