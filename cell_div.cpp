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
//	cout << "in division function" << endl;
	double area = this->calc_Area();
//	cout << "area calculated " << area << endl;
	if(area > AREA_DOUBLED) {
		if(this->layer == 1) { 
			cout << "Cell " << this->rank << "  passed area threshold for division lengthwise" << endl;
			sister = this->divide_length_wise();
			cout << "divided" << endl;
		}
		if(this->layer == 3) {
			cout << "Cell " << this->rank << " passed area threshold for division widthwise" << endl;
			sister = this->divide_width_wise();
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
	bool islength = true;
//	cout << "Made new cell pointer" << endl;
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

	Wall_Node* left_start = down->get_Right_Neighbor()->get_Right_Neighbor()->get_Right_Neighbor();
	Wall_Node* right_start = down->get_Left_Neighbor()->get_Left_Neighbor()->get_Left_Neighbor();
	
	delete down->get_Right_Neighbor()->get_Right_Neighbor();
	delete  down->get_Right_Neighbor();
	delete down->get_Left_Neighbor()->get_Left_Neighbor();
	delete 	down->get_Left_Neighbor();
	delete  down;

	Wall_Node* left_end = up->get_Left_Neighbor()->get_Left_Neighbor()->get_Left_Neighbor();
	Wall_Node* right_end = up->get_Right_Neighbor()->get_Right_Neighbor()->get_Right_Neighbor();

	delete up->get_Left_Neighbor()->get_Left_Neighbor();
	delete up->get_Right_Neighbor()->get_Right_Neighbor();
	delete up;

	double left_slope = (left_end->get_Location().get_Y() - left_start->get_Location().get_Y())/(left_end->get_Location().get_X() - left_start->get_Location().get_X());
	double right_slope = (right_end->get_Location().get_Y() - right_start->get_Location().get_Y())/(right_end->get_Location().get_X() - right_start->get_Location().get_X());
	
	double b_left = left_end->get_Location().get_Y()-left_slope*left_end->get_Location().get_X();
	double b_right = right_end->get_Location().get_Y()-right_slope*right_end->get_Location().get_X();

	double left_length = (left_end->get_Location() - left_start->get_Location()).length();
	double right_length = (right_end->get_Location() - right_start->get_Location()).length();

	int total_num_left = (int) (left_length/(MembrEquLen*2));
	int total_num_right = (int) (right_length/(MembrEquLen*2));

	double x_left = sqrt(pow(left_end->get_Location().get_X() - left_start->get_Location().get_X(),2));
	double x_right = sqrt(pow(right_end->get_Location().get_X() - right_start->get_Location().get_X(),2));

	double delta_x_left = x_left/total_num_left;
	double delta_x_right = x_right/total_num_right;

	double start_x_left = left_start->get_Location().get_X();
	double start_x_right = right_start->get_Location().get_X();
	
//	cout << "make left side" << endl;
	Wall_Node* curr = NULL;
	Coord curr_coord;
	Wall_Node* prev = left_start;
	double curr_x = start_x_left;
	for(unsigned int i = 0; i < total_num_left; i++) {
		curr_x = curr_x + delta_x_left;
		curr_coord = Coord(curr_x + 0.04, left_slope*curr_x + b_left);
		curr = new Wall_Node(curr_coord, this);
		//cout << "Setting neighbors: " << i << endl;
		prev->set_Left_Neighbor(curr);
		curr->set_Right_Neighbor(prev);
		prev = curr;
	}
	curr->set_Left_Neighbor(left_end);
	left_end->set_Right_Neighbor(curr);
	left_Corner = left_start;

//	cout << "make right side" << endl;;
	prev = right_start;
	curr_x = start_x_right;

	for(int i = 0; i< total_num_right; i++) {
		curr_x = curr_x + delta_x_right;
		curr_coord = Coord(curr_x - .04, right_slope*curr_x + b_right);
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
	
	//delete all old cyt nodes
	int num_cyts = cyt_nodes.size();
	int new_cyt_cnt = num_cyts/2;
	Cyt_Node* c = NULL;
	while(!cyt_nodes.empty()) {
		c= cyt_nodes.at(cyt_nodes.size() -1);
		delete c;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}	
	
	//cout << "Finished deleting old cyt nodes" << endl;
	
//	cout << "get most up/down and left/right for radius" << endl;
	Wall_Node* up_sis = NULL;
	Wall_Node* down_sis = NULL;
	Wall_Node* left_sis = NULL;
	Wall_Node* right_sis = NULL;
	this->closest_node_top(up);
	this->closest_node_bottom(down);
	if(up == NULL) {

	cout << "this up null" << endl;
	}
	if(down == NULL) {
		cout << "this down null" << endl;
	}
	sister->closest_node_top(up_sis);
	sister->closest_node_bottom(down_sis);
	if(up_sis == NULL) {
		cout << "up sis null" << endl;
	}
	if(down_sis == NULL) {
		cout << "down sis null" << endl;
	}
	this->closest_node_left(left);
	this->closest_node_right(right);
	if(left == NULL) {
		cout << "left null" << endl;
	}
	if(right==NULL) {
		cout << "right null" << endl;
	}
	sister->closest_node_left(left_sis);
	sister->closest_node_right(right_sis);
	if(left_sis == NULL) {

	cout << "left sis null" << endl;
	}
	if(right_sis == NULL) {
		cout << "right sis null" << endl;
	}


	double radius_x = ((left->get_Location() - right->get_Location()).length())*0.5;
	double radius_y = ((up->get_Location() - down->get_Location()).length())*0.5;
//	cout << "first radius" << endl;
	double radius_x_s = ((left_sis->get_Location() - right_sis->get_Location()).length())*0.5;
	double radius_y_s = ((up_sis->get_Location() - down_sis->get_Location()).length())*0.5;
//	cout << "second radius" << endl;
	//create new ones for each cell
	for(int i = 0; i < new_cyt_cnt; i++) {
		this->add_Cyt_Node_Div(radius_x, radius_y,islength);
		sister->add_Cyt_Node_Div(radius_x_s, radius_y_s,islength);
	}
	return sister;
}

Cell* Cell::divide_width_wise() {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as the left sister
	//	this functoin will create a sister cell to the right
	//	and return it to the tissue
	bool islength = false;
//	cout << "Made new cell pointer" << endl;
	Cell* sister = new Cell(my_tissue);
	Wall_Node* up = NULL;
	Wall_Node* down = NULL;
	Wall_Node* left = NULL;
	Wall_Node* right = NULL;
	this->closest_node_left(left);
/*	if(up == NULL) {
		cout << "Top NULL" << endl;
		exit(1);
	}*/
	this->closest_node_right(right);
/*	if(down == NULL) {
		cout << "Bottom NULL" << endl;
		exit(1);
	}*/

	Wall_Node* top_start = left->get_Right_Neighbor()->get_Right_Neighbor()->get_Right_Neighbor();
	Wall_Node* bottom_start = left->get_Left_Neighbor()->get_Left_Neighbor()->get_Left_Neighbor();
	
	delete left->get_Right_Neighbor()->get_Right_Neighbor();
	delete left->get_Right_Neighbor();
	delete left->get_Left_Neighbor()->get_Left_Neighbor();
	delete left->get_Left_Neighbor();
	delete left;

	Wall_Node* top_end = right->get_Left_Neighbor()->get_Left_Neighbor()->get_Left_Neighbor();
	Wall_Node* bottom_end = right->get_Right_Neighbor()->get_Right_Neighbor()->get_Right_Neighbor();

	delete right->get_Left_Neighbor()->get_Left_Neighbor();
	delete right->get_Right_Neighbor()->get_Right_Neighbor();
	delete right->get_Left_Neighbor();
	delete right->get_Right_Neighbor();
	delete right;

	double top_slope = (top_end->get_Location().get_Y() - top_start->get_Location().get_Y())/(top_end->get_Location().get_X() - top_start->get_Location().get_X());
	double bottom_slope = (bottom_end->get_Location().get_Y() - bottom_start->get_Location().get_Y())/(bottom_end->get_Location().get_X() - bottom_start->get_Location().get_X());
	
	double b_top = top_end->get_Location().get_Y()-top_slope*top_end->get_Location().get_X();
	double b_bottom = bottom_end->get_Location().get_Y()-bottom_slope*bottom_end->get_Location().get_X();

	double top_length = (top_end->get_Location() - top_start->get_Location()).length();
	double bottom_length = (bottom_end->get_Location() - bottom_start->get_Location()).length();

	int total_num_top = (int) (top_length/(MembrEquLen*2));
	int total_num_bottom = (int) (bottom_length/(MembrEquLen*2));

	double x_top = sqrt(pow(top_end->get_Location().get_X() - top_start->get_Location().get_X(),2));
	double x_bottom = sqrt(pow(bottom_end->get_Location().get_X() - bottom_start->get_Location().get_X(),2));

	double delta_x_top = x_top/total_num_top;
	double delta_x_bottom = x_bottom/total_num_bottom;

	double start_x_top = top_start->get_Location().get_X();
	double start_x_bottom = bottom_start->get_Location().get_X();
	
//	cout << "make left side" << endl;
	Wall_Node* curr = NULL;
	Coord curr_coord;
	Wall_Node* prev = top_start;
	double curr_x = start_x_top;
	for(unsigned int i = 0; i < total_num_top; i++) {
		curr_x = curr_x + delta_x_top;
		curr_coord = Coord(curr_x + 0.04, top_slope*curr_x + b_top);
		curr = new Wall_Node(curr_coord, this);
		//cout << "Setting neighbors: " << i << endl;
		prev->set_Left_Neighbor(curr);
		curr->set_Right_Neighbor(prev);
		prev = curr;
	}
	curr->set_Left_Neighbor(top_end);
	top_end->set_Right_Neighbor(curr);
	left_Corner = top_start;

//	cout << "make right side" << endl;;
	prev = bottom_start;
	curr_x = start_x_bottom;

	for(int i = 0; i< total_num_bottom; i++) {
		curr_x = curr_x + delta_x_bottom;
		curr_coord = Coord(curr_x - .04, bottom_slope*curr_x + b_bottom);
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
	
	//delete all old cyt nodes
	int num_cyts = cyt_nodes.size();
	int new_cyt_cnt = num_cyts/2;
	Cyt_Node* c = NULL;
	while(!cyt_nodes.empty()) {
		c= cyt_nodes.at(cyt_nodes.size() -1);
		delete c;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}	
	
	//cout << "Finished deleting old cyt nodes" << endl;
	
//	cout << "get most up/down and left/right for radius" << endl;
	Wall_Node* up_sis = NULL;
	Wall_Node* down_sis = NULL;
	Wall_Node* left_sis = NULL;
	Wall_Node* right_sis = NULL;
	this->closest_node_top(up);
	this->closest_node_bottom(down);
	sister->closest_node_top(up_sis);
	sister->closest_node_bottom(down_sis);
	this->closest_node_left(left);
	this->closest_node_right(right);
	sister->closest_node_left(left_sis);
	sister->closest_node_right(right_sis);

	double radius_x = ((left->get_Location() - right->get_Location()).length())*0.5;
	double radius_y = ((up->get_Location() - down->get_Location()).length())*0.5;
	double radius_x_s = ((left_sis->get_Location() - right_sis->get_Location()).length())*0.5;
	double radius_y_s = ((up_sis->get_Location() - down_sis->get_Location()).length())*0.5;
	//create new ones for each cell
	for(int i = 0; i < new_cyt_cnt; i++) {
		this->add_Cyt_Node_Div(radius_x, radius_y, islength);
		sister->add_Cyt_Node_Div(radius_x_s, radius_y_s, islength);
	}
	return sister;
}

