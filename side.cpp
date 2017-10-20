//side.cpp
//======================
// Forward Declarations

//======================
// Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <stdlib.h>

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "side.h"
#include "cell.h"
//=====================

//Constructors
Side::Side(Coord locA, Coord locZ, Cell* cell, int num_nodes) {

	my_cell = cell;
	end_A = new Wall_Node(locA, this);
	end_Z = new Wall_Node(locZ, this);
//	cout << "made end a and z" << endl;
	Coord diff = ((locZ - locA) / (num_nodes - 1));
	Coord curr_loc = locA + diff;
	Wall_Node* curr;
	Wall_Node* prev = end_A;
 
	for (int i = 0; i < num_nodes - 2; i++) {
		curr = new Wall_Node(curr_loc, this);
		curr->set_Right_Neighbor(prev);
		prev->set_Left_Neighbor(curr);
		curr_loc += diff;
		prev = curr;
	}
//	cout << " made in between nodes" << endl;
	//link up end_Z
	end_Z->set_Right_Neighbor(prev);
	prev->set_Left_Neighbor(end_Z);

	num_wall_nodes = num_nodes;
//	cout << "linked ends " << endl;
	//set curve parameters on each end
	curr = end_A;
	for (int i = 0; i < 3; i++) {
		curr->set_Curve(true);
		curr->set_Equi_Angle(thetaCurve);
		curr = curr->get_Left_Neighbor();
	}

	curr = end_Z;
	for (int i = 0; i < 3; i++) {
		curr->set_Curve(true);
		curr->set_Equi_Angle(thetaCurve);
		curr = curr->get_Right_Neighbor();
	}
//	cout << "set curve nodes" << endl;
}

void Side::connect_Ends(Side* s) {
	//connect end_Z of curr side to end_A of side s

	Wall_Node* other_A = s->get_End_A();

	end_Z->set_Left_Neighbor(other_A);
	other_A->set_Right_Neighbor(end_Z);

	return;
}

//Destructor
Side::~Side() {

	my_cell = NULL;
	end_A = NULL;
	end_Z = NULL;

}

//=================================================================
//=======================================
//Getters and Setters
//=======================================
//=================================================================
void Side::get_Cyt_Nodes(vector<Cyt_Node*>& cyts) {
	my_cell->get_Cyt_Nodes(cyts);
	return;
}

void Side::set_Phys_Parameters(double kbend, double klin) {
	bending_spring = kbend;
	linear_spring = klin;
	return;
}

void Side::set_Side_Type(int type) {
	side_type = type;
	return;
}
void Side::set_My_Cell(Cell* new_cell) {
	this-> my_cell = new_cell;
	return;
}
void Side::get_Touching_Sides(vector<Side*>& touching_sides) {
	touching_sides = this->touching_Sides;
	return;
}

void Side::update_Touching_Sides(vector<Cell*>& neighbor_Cells) {
	Coord my_Mid_Point = (this->get_End_A()->get_Location() + this->get_End_Z()->get_Location())*0.5;
	double my_Mid_Point_X = my_Mid_Point.get_X();
	double my_Mid_Point_Y = my_Mid_Point.get_Y();
	
	vector<Cell*>touching_Neighbors;
	vector<Side*>curr_Cell_Sides;
	Cell* curr_Cell = NULL;
	Coord curr_Cell_Left;
	Coord curr_Cell_Right;
	for(int i = 0; i < neighbor_Cells.size(); i++) {
		curr_Cell = neighbor_Cells.at(i);
		curr_Cell->get_Sides(curr_Cell_Sides);
		curr_Cell_Left = curr_Cell_Sides.at(2)->get_End_Z()->get_Location();
		curr_Cell_Right = curr_Cell_Sides.at(0)->get_End_A()->get_Location();	
		//for side 0 need the y value of the top left corner
		//needs to be less than the y value of the midpoint
		if((side_type == 0) && (curr_Cell_Left.get_Y() < my_Mid_Point_Y)) {
			touching_Neighbors.push_back(curr_Cell);
		}
		//for side 1 the x value of the top left corner 
		//needs to be greater than the x value of midpoint
		else if((side_type == 1) && (curr_Cell_Left.get_X() > my_Mid_Point_X)) {
			touching_Neighbors.push_back(curr_Cell);
	 	}
		//for side 2 the y value of the bottom right corner 
		//needs to be greater than the y valu of midpoint
		else if((side_type == 2) && (curr_Cell_Right.get_Y() > my_Mid_Point_Y)) {
			touching_Neighbors.push_back(curr_Cell);
		}
		//for side 3 the x value of the bottom right corner
		//needs to be less than the x value of the midpoint
		else if((side_type == 3) && (curr_Cell_Right.get_X() < my_Mid_Point_X)) {
			touching_Neighbors.push_back(curr_Cell);
		}
	}

	vector<Side*>touching_Sides;
	Side* side = NULL;
	for(int i = 0; i < touching_Neighbors.size(); i++) {
		curr_Cell = touching_Neighbors.at(i);
		curr_Cell->get_Sides(curr_Cell_Sides);
		if(side_type == 0) {
			side = curr_Cell_Sides.at(2);
			touching_Sides.push_back(side);
		}
		else if(side_type == 1) {
			side = curr_Cell_Sides.at(3);
			touching_Sides.push_back(side);
		}
		else if(side_type == 2) {
			side = curr_Cell_Sides.at(0);
			touching_Sides.push_back(side);
		}
		else if(side_type == 3) {
			side = curr_Cell_Sides.at(1);
			touching_Sides.push_back(side);
		}
	}

	cout << "Side is: " << side_type << endl;
	cout << "With touching neighbors size: " << touching_Neighbors.size() << endl;
	this->touching_Sides = touching_Sides;
	cout << "With touching sides size: " << touching_Sides.size() << endl;
	return;
}



		
		

//=================================================================
//=======================================
//Cell
//=======================================
//=================================================================

void Side::add_Wall_Node(Wall_Node* right) {

	Wall_Node* left = right->get_Left_Neighbor();

	Coord new_loc = ((right->get_Location() + left->get_Location()) * 0.5);

	Wall_Node* new_node = new Wall_Node(new_loc, this, left, right);
	num_wall_nodes++;

	right->set_Left_Neighbor(new_node);
	left->set_Right_Neighbor(new_node);


	//fix curves
	if (right->is_Curve() && left->is_Curve()) {
		
		Wall_Node* curr = end_A;

		for (int i = 0; i < 3; i++) {
			curr->set_Curve(true);
			curr->set_Equi_Angle(thetaCurve);
			curr = curr->get_Left_Neighbor();
		}

		curr->set_Curve(false);
		curr->set_Equi_Angle(thetaFlat);


		int num_iter;
		if (right == end_Z) {
			end_Z = new_node;
		}

		curr = end_Z;
		
		for (int i = 0; i < 3; i++) {
			curr->set_Curve(true);
			curr->set_Equi_Angle(thetaCurve);
			curr = curr->get_Right_Neighbor();
		}

		curr->set_Curve(false);
		curr->set_Equi_Angle(thetaFlat);
	}
	


	return;
}







//======================
//end side.cpp
