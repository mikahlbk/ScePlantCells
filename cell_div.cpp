//cell_div.cpp
//========================
//Forward Declarations
//


//========================
//Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
//==========================

Cell* Cell::divide(int Ti) {
	Cell* sister = NULL;
	//calculate area
	cout << "in division function" << endl;
	double area = this->calc_Area();
	cout << "area calculated " << area << endl;
	if(area > AREA_DOUBLED) {
		if(this->layer == 0) { 
			cout << "Cell " << this->rank << "  passed area threshold for division" << endl;
			sister = this->divide_length_wise();
			cout << "divided" << endl;
		}
		if(this->layer == 3) {
			cout << "Cell " << this->rank << " passed area threshold for division" << endl;
			//sister = this->divide_width_wise();
			cout << "divided" << endl;
		}
	}
	return sister;
}


Cell* Cell::divide_length_wise() {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as the left sister
	//	this functoin will create a sister cell to the right
	//	and return it to the tissue
	cout << "Made new cell pointer" << endl;
	Cell* sister = new Cell(my_tissue);
	Wall_Node* up = NULL;
	Wall_Node* down = NULL;
	Wall_Node* left = NULL;
	Wall_Node* right = NULL;
	this->closest_node_top(up);
	if(up == NULL) {
		cout << "Top NULL" << endl;
		exit(1);
	}
	this->closest_node_bottom(down);
	if(down == NULL) {
		cout << "Bottom NULL" << endl;
		exit(1);
	}

	Wall_Node* left_start = down->get_Right_Neighbor()->get_Right_Neighbor();
	Wall_Node* right_start = down->get_Left_Neighbor()->get_Left_Neighbor();

	delete  down->get_Right_Neighbor();
	delete 	down->get_Left_Neighbor();
	delete  down;

	Wall_Node* left_end = up->get_Left_Neighbor()->get_Left_Neighbor();
	Wall_Node* right_end = up->get_Right_Neighbor()->get_Right_Neighbor();

	delete up->get_Left_Neighbor();
	delete up->get_Right_Neighbor();
	delete up;

	double left_length = (left_end->get_Location() - left_start->get_Location()).length();
	double right_length = (right_end->get_Location() - right_start->get_Location()).length();

	int total_num_left = (int) (left_length/(MembrEquLen*2))-1;
	int total_num_right = (int) (right_length/(MembrEquLen*2))-1;

	double left_space = left_length/total_num_left;
	double right_space = right_length/total_num_right;
	cout << "Total num:" << total_num_left << endl;
	Coord left_increment = Coord(0,left_space);
	Coord right_increment = Coord(0, right_space);

	cout << "make left side" << endl;
	Wall_Node* curr = NULL;
	Coord curr_coord= left_start->get_Location()+ Coord(0.015, 0);
	Wall_Node* prev = left_start;

	for(int i = 0; i< total_num_left; i++) {
		curr_coord = curr_coord + left_increment;
		curr = new Wall_Node(curr_coord, this);
		//cout << "Setting neighbors: " << i << endl;
		prev->set_Left_Neighbor(curr);
		curr->set_Right_Neighbor(prev);
		prev = curr;
	}
	curr->set_Left_Neighbor(left_end);
	left_end->set_Right_Neighbor(curr);
	left_Corner = left_start;

	cout << "make right side" << endl;;
	curr_coord = right_start->get_Location()- Coord(0.015, 0);
	prev = right_start;

	for(int i = 0; i< total_num_right; i++) {
		curr_coord = curr_coord + right_increment;
		curr = new Wall_Node(curr_coord, sister);
		//cout << "Setting neighbors: " << i << endl;
		prev->set_Right_Neighbor(curr);
		curr->set_Left_Neighbor(prev);
		prev = curr;
	}
	curr->set_Right_Neighbor(right_end);
	right_end->set_Left_Neighbor(curr);
	sister->set_Left_Corner(right_start);

	//count wall nodes
	//cout << "begin count wall nodes" << endl;
	curr = this->left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	int number_nodes_A = 0;
	//cout << "cell a counting" << endl;
	do {
		number_nodes_A++;
		curr->update_Cell(this);
		next = curr->get_Left_Neighbor();
		//cout << number_nodes_A << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while(next != orig);
	//cout << "done cell a counting" << endl;	
	this->set_Wall_Count(number_nodes_A);
	
	curr = sister->get_Left_Corner();
	next = NULL;
	orig = curr;
	int number_nodes_B = 0;
	//cout << "Cell b counting" << endl;
	do {
		number_nodes_B++;
		curr->update_Cell(sister);
		next = curr->get_Left_Neighbor();
		//cout << number_nodes_B << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while (next != orig);
	//cout << "done counting cell b" << endl;
	sister->set_Wall_Count(number_nodes_B);
	
	//cout << "updating angles" << endl;
	//update wall angles
	this->update_Wall_Angles();
	sister->update_Wall_Angles();
	this->update_Wall_Equi_Angles();
	sister->update_Wall_Equi_Angles();
	//cout << "updating center" << endl;
	//update cell center
	this->update_Cell_Center();
	sister->update_Cell_Center();

	//cout << "update layer" << endl;
	//update layer information
	sister->set_Layer(this->layer);
	//cout << "update growth rate" << endl;
	sister->set_growth_rate(this->growth_rate);
		
	//cout << "Updated Angles and cell centers"<< endl;
	
	//distribute cyt nodes between sister cells
	//cout << "deleting cyt nodes" << endl;
	int new_cyt_cnt = 15;
	//delete all old cyt nodes
	Cyt_Node* c = NULL;
	while(!cyt_nodes.empty()) {
		c= cyt_nodes.at(cyt_nodes.size() -1);
		delete c;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}	
	
	//cout << "Finished deleting old cyt nodes" << endl;
	
	cout << "get most up/down and left/right for radius" << endl;
	Wall_Node* up_sis = NULL;
	Wall_Node* down_sis = NULL;
	Wall_Node* left_sis= NULL;
	Wall_Node* right_sis = NULL;
	//Wall_Node* closest = NULL;
	//Wall_Node* closest_sis = NULL;
	this->closest_node_top(up);
	this->closest_node_bottom(down);
	sister->closest_node_top(up_sis);
	sister->closest_node_bottom(down_sis);
	this->closest_node_left(left);
	this->closest_node_right(right);
	sister->closest_node_left(left_sis);
	sister->closest_node_left(left_sis);

	double radius_x = ((right->get_Location() - left->get_Location()).length())*0.5;
	double radius_x_s = ((right_sis->get_Location() - left_sis->get_Location()).length())*0.5;
	double radius_y = ((up->get_Location() - down->get_Location()).length())*0.5;
	double radius_y_s = ((up_sis->get_Location() - down_sis->get_Location()).length())*0.5;
	//create new ones for each cell
	for(int i = 0; i < new_cyt_cnt; i++) {
		this->add_Cyt_Node_Div(radius_x,radius_y);
		sister->add_Cyt_Node_Div(radius_x_s, radius_y_s);
	}
	return sister;
}

/*Cell* Cell::divide_width_wise() {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as the left sister
	//	this functoin will create a sister cell to the right
	//	and return it to the tissue
//	cout << "Made new cell pointer" << endl;
	Cell* sister = new Cell(my_tissue);
	Wall_Node* up = NULL;
	Wall_Node* down = NULL;
	Wall_Node* left = NULL;
	Wall_Node* right = NULL;
	this->closest_node_left(left);
	if(left == NULL) {
		exit(1);
	}
	this->closest_node_right(right);
	if(right == NULL) {
		exit(1);
	}
	//	this->most_Up_Down(up, down);
	Wall_Node* top_start = left->get_Right_Neighbor()->get_Right_Neighbor();
	Wall_Node* bottom_start = left->get_Left_Neighbor()->get_Left_Neighbor();

	delete  left->get_Right_Neighbor();
	delete 	left->get_Left_Neighbor();
	delete  left;

	Wall_Node* top_end = right->get_Left_Neighbor()->get_Left_Neighbor();
	Wall_Node* bottom_end = right->get_Right_Neighbor()->get_Right_Neighbor();

	delete right->get_Left_Neighbor();
	delete right->get_Right_Neighbor();
	delete right;

	double top_length = (top_end->get_Location() - top_start->get_Location()).length();
	double bottom_length = (bottom_end->get_Location() - bottom_start->get_Location()).length();

	double top_space = top_length/(MembrEquLen*2);
	double bottom_space = bottom_length/(MembrEquLen*2);
	int total_num_top = ((int) top_space)-1;
	int total_num_bottom = ((int) bottom_space)-1;
//	cout << "Total num top:" << total_num_top << endl;
//	cout << "total num bottom: " << total_num_bottom << endl; 
	Coord top_increment = Coord(top_space,0);
	Coord bottom_increment = Coord(bottom_space,0);

//	cout << "make top side" << endl;
	Wall_Node* curr = NULL;
	Coord curr_coord = top_start->get_Location()- Coord(0,0.015);
	Wall_Node* prev = top_start;

	for(int i = 0; i< total_num_top; i++) {
		curr_coord = curr_coord + top_increment;
		curr = new Wall_Node(curr_coord, this);
		//cout << "Setting neighbors: " << i << endl;
		prev->set_Left_Neighbor(curr);
		curr->set_Right_Neighbor(prev);
		prev = curr;
	}
	curr->set_Left_Neighbor(top_end);
	top_end->set_Right_Neighbor(curr);
	left_Corner = top_start;

	//cout << "make bottom side" << endl;;
	curr_coord = bottom_start->get_Location()+ Coord(0,0.015);
	prev = bottom_start;

	for(int i = 0; i< total_num_bottom; i++) {
		curr_coord = curr_coord + bottom_increment;
		curr = new Wall_Node(curr_coord, sister);
		//cout << "Setting neighbors: " << i << endl;
		prev->set_Right_Neighbor(curr);
		curr->set_Left_Neighbor(prev);
		prev = curr;
	}
	curr->set_Right_Neighbor(bottom_end);
	bottom_end->set_Left_Neighbor(curr);
	sister->set_Left_Corner(bottom_start);

	//count wall nodes
	//cout << "begin count wall nodes" << endl;
	curr = this->left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	int number_nodes_A = 0;
	//cout << "cell a counting" << endl;
	do {
		number_nodes_A++;
		curr->update_Cell(this);
		next = curr->get_Left_Neighbor();
	//	cout << number_nodes_A << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while(next != orig);
	//cout << "done cell a counting" << endl;	
	this->set_Wall_Count(number_nodes_A);
	
	curr = sister->get_Left_Corner();
	next = NULL;
	orig = curr;
	int number_nodes_B = 0;
	//cout << "Cell b counting" << endl;
	do {
		number_nodes_B++;
		curr->update_Cell(sister);
		next = curr->get_Left_Neighbor();
	//	cout << number_nodes_B << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while (next != orig);
	//cout << "done counting cell b" << endl;
	sister->set_Wall_Count(number_nodes_B);
	
	//cout << "updating angles" << endl;
	//update wall angles
	this->update_Wall_Angles();
	sister->update_Wall_Angles();
	this->update_Wall_Equi_Angles();
	sister->update_Wall_Equi_Angles();
	//cout << "updating center" << endl;
	//update cell center
	this->update_Cell_Center();
	sister->update_Cell_Center();

	//cout << "update layer" << endl;
	//update layer information
	sister->set_Layer(this->layer);
	//cout << "update growth rate" << endl;
	sister->set_growth_rate(this->growth_rate);
		
	//cout << "Updated Angles and cell centers"<< endl;
	
	//distribute cyt nodes between sister cells
//	cout << "deleting cyt nodes" << endl;
	int new_cyt_cnt = 25;
	//delete all old cyt nodes
	Cyt_Node* c = NULL;
	while(!cyt_nodes.empty()) {
		c= cyt_nodes.at(cyt_nodes.size() -1);
		delete c;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}	
	
//	cout << "Finished deleting old cyt nodes" << endl;
	
//	cout << "get most up/down and left/right for radius" << endl;
	Wall_Node* up_sis = NULL;
	Wall_Node* down_sis = NULL;
	Wall_Node* left_sis= NULL;
	Wall_Node* right_sis = NULL;
	Wall_Node* closest = NULL;
	Wall_Node* closest_sis = NULL;
//	this->closest_node_top(up);
//	this->closest_node_bottom(down);
//	sister->closest_node_top(up_sis);
//	sister->closest_node_bottom(down_sis);
//	this->closest_node_left(left);
	this->closest_node(closest);
	sister->closest_node(closest_sis);
//	sister->closest_node_right(right_sis);
//	cout << "closest " << closest << "closest sis " << closest_sis << endl; // " up " << up << " down " << down << endl;	
	cout << "set new raddi" << endl;
	double radius = (closest->get_Location() - cell_center).length();
//	cout << "got right" << endl;
	double radius_sis = (closest_sis->get_Location() - sister->get_Cell_Center()).length();
//	cout << "got x" << endl;
	//double radius_y = ((up->get_Location() - down->get_Location()).length())*0.5;
	//double radius_y_s = ((up_sis->get_Location() - down_sis->get_Location()).length())*0.5;
	
//	cout << "create new ones for each cell" << endl;
	for(int i = 0; i < new_cyt_cnt; i++) {
		this->add_Cyt_Node_Div(radius);
		sister->add_Cyt_Node_Div(radius_sis);
	}
	return sister;
}*/


