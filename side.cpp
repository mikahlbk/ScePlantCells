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
	
	Coord diff = ((locZ - locA) / (num_nodes - 1));
	Coord curr_loc = locA + diff;
	Wall_Node* curr;
	Wall_Node* prev = end_A;

	for (int i = 0; i < num_nodes - 2; i++) {
		curr = new Wall_Node(curr_loc, this);
		curr->set_Right_Neighbor(prev);
		prev->set_Left_Neighbor(curr);
		curr_loc += diff;
	}
	//link up end_Z
	end_Z->set_Right_Neighbor(prev);
	prev->set_Left_Neighbor(end_Z);

	num_wall_nodes = num_nodes;

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

//=================================================================
//=======================================
//Cell Growth
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
			curr = curr->get_Left_Neighbor();
		}

		curr->set_Curve(false);
		curr->set_Equi_Angle(thetaFlat);
	}
	

	return;
}







//======================
//end side.cpp
