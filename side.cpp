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

void Side::set_My_Cell(Cell* new_cell) {
	this-> my_cell = new_cell;
	return;
}
void Side::get_Touching_Neighbors(vector<Cell*>& touching_neighbors) {
	touching_neighbors = this->touching_Neighbors;
	return;
}
//void Side::update_Touching_Neighbors(vector<Cell*>* touching_Neighbors) {
//	this->touching_Neighbors = touching_Neighbors;
//	return;
//}

void Side::update_Touching_Neighbors(vector<Cell*>& neighbor_Cells) {
	Coord my_Cell_Center = my_cell->get_Cell_Center();
	double my_Cell_Center_X = my_Cell_Center.get_X();
	double my_Cell_Center_Y = my_Cell_Center.get_Y();
	Coord my_Mid_Point = (this->get_End_A()->get_Location() + this->get_End_Z()->get_Location())*0.5;
	double my_Mid_Point_X = my_Mid_Point.get_X();
	double my_Mid_Point_Y = my_Mid_Point.get_Y();
	string side;
	vector<Cell*> touching_Neighbors;
	Cell* curr_Cell = NULL;
	Coord curr_Cell_Center;
	double center_X = 0;
	double center_Y = 0;
//	cout << "Side is : " << side << endl;
	if((this->get_End_A()->get_Location().get_Y() < my_Cell_Center_Y) && (this->get_End_A()->get_Location().get_X() < my_Cell_Center_X)) {
		side = "bottom";
	}
	else if((this->get_End_A()->get_Location().get_Y() < my_Cell_Center_Y) && (this->get_End_A()->get_Location().get_X() > my_Cell_Center_X)) {
		side = "right";
	}
	else if((this->get_End_A()->get_Location().get_Y() > my_Cell_Center_Y) && (this->get_End_A()->get_Location().get_X() > my_Cell_Center_X)) {
		side = "top";
	}
	else {
		side = "left";
	}
//	cout << "Side is : " << side << endl;
//	cout << "Number of neighbor cells: " << neighbor_Cells.size() << endl;
	for(int i = 0; i < neighbor_Cells.size(); i++) {
		curr_Cell = neighbor_Cells.at(i);
		curr_Cell_Center = curr_Cell->get_Cell_Center();
		center_X = curr_Cell_Center.get_X();
		center_Y = curr_Cell_Center.get_Y();
		if((side == "bottom") && (center_Y < this->get_End_A()->get_Location().get_Y())) {
			touching_Neighbors.push_back(curr_Cell);
		}
		else if((side == "right") && (center_X > this->get_End_A()->get_Location().get_X())) {
			touching_Neighbors.push_back(curr_Cell);
		}
		else if((side == "top") && (center_Y > this->get_End_A()->get_Location().get_Y())) {
			touching_Neighbors.push_back(curr_Cell);
		}
		else if((side == "left") && (center_X < this->get_End_A()->get_Location().get_X())) {
			touching_Neighbors.push_back(curr_Cell);
		}

	}
	/*Coord curr_Side_Mid_Point;
	vector<Side*> neighbor_Sides;
	vector<Side*> curr_Cell_Sides;
	Side* curr_Side = NULL;
	double curr_len = 0;
	double smallest = 100;
	for(int j = 0; j < touching_Neighbors.size();j++) {
		curr_Cell = touching_Neighbors.at(j);
		cout << "Rank : " << curr_Cell->get_Rank() << endl;
		curr_Cell->get_Sides(curr_Cell_Sides);
		for(int i = 0; i < curr_Cell_Sides.size(); i++) {
			curr_Side = curr_Cell_Sides.at(i);
			curr_Side_Mid_Point = (curr_Side->get_End_A()->get_Location() + curr_Side->get_End_Z()->get_Location())*.5;
			curr_len = (my_Mid_Point - curr_Side_Mid_Point).length();
			if(curr_len < smallest) {
				neighbor_Sides.push_back(curr_Side);
				smallest = curr_len;
			}
		}
	}*/
	
//	cout << "Side is: " << side << endl;
//	cout << "With touching neighbors size: " << touching_Neighbors.size() << endl;
	this->touching_Neighbors = touching_Neighbors;
//	cout << "After assigment: " << this->touching_Neighbors.size() << endl;
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
