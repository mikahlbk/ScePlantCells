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

Cell* Cell::divide() {
	Cell* sister = NULL;
	//calculate area
	double area = this->calc_Area();
	if(area > AREA_DOUBLED) {
	//	if((this->layer == 3)|| (this->layer == 4)|| (this->layer == 5)) {
			//ut << "cell will divide" << endl;
			//sister = this->divide_width_wise(Ti);		
	//	}
	//	else {
	//	cout << "Will Divide" << endl;	
		sister = this->divide_length_wise();
		//}
		//cout <<  "Divided" << endl;
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
	this->most_Up_Down();
	Wall_Node* left_start = most_down->get_Right_Neighbor()->get_Right_Neighbor();
	Wall_Node* right_start = most_down->get_Left_Neighbor()->get_Left_Neighbor();

	delete  most_down->get_Right_Neighbor();
	delete 	most_down->get_Left_Neighbor();
	delete most_down;

	Wall_Node* left_end = most_up->get_Left_Neighbor()->get_Left_Neighbor();
	Wall_Node* right_end = most_up->get_Right_Neighbor()->get_Right_Neighbor();

	delete most_up->get_Left_Neighbor();
	delete most_up->get_Right_Neighbor();
	delete most_up;

	double left_length = (left_end->get_Location() - left_start->get_Location()).length();
	double right_length = (right_end->get_Location() - right_start->get_Location()).length();

	double left_space = left_length/MembrEquLen;
	double right_space = right_length/MembrEquLen;
	int total_num = 20;
	cout << "Total num:" << total_num << endl;
	Coord left_increment = Coord(0,left_space);
	Coord right_increment = Coord(0, right_space);

	cout << "make left side" << endl;
	Wall_Node* curr = NULL;
	Coord curr_coord = left_start->get_Location();
	Wall_Node* prev = left_start;

	for(int i = 0; i< total_num; i++) {
		curr_coord = curr_coord + left_increment;
		curr = new Wall_Node(curr_coord, this);
		cout << "Setting neighbors" << endl;
		prev->set_Left_Neighbor(curr);
		curr->set_Right_Neighbor(prev);
		prev = curr;
	}
	curr->set_Left_Neighbor(left_end);
	left_end->set_Right_Neighbor(curr);
	this->left_Corner = left_start;

	cout << "make right side" << endl;
	curr_coord = right_start->get_Location();
	prev = right_start;

	for(int i = 0; i< total_num; i++) {
		curr_coord = curr_coord + right_increment;
		curr = new Wall_Node(curr_coord, sister);
		cout << "Setting neighbors" << endl;
		prev->set_Right_Neighbor(curr);
		curr->set_Left_Neighbor(prev);
		prev = curr;
	}
	curr->set_Right_Neighbor(left_end);
	left_end->set_Left_Neighbor(curr);
	sister->set_Left_Corner(right_start);

	//count wall nodes
	cout << "begin count wall nodes" << endl;
	curr = this->left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	int number_nodes_A = 0;
	cout << "cell a counting" << endl;
	do {
		number_nodes_A++;
		next = curr->get_Left_Neighbor();
		cout << number_nodes_A++ << endl;
		if(next == NULL) {
			cout << "notlinked" << endl;
			exit(1);
		}
		curr = next;
	} while(next != orig);
	cout << "done cell a counting" << endl;	
	this->set_Wall_Count(number_nodes_A);
	
	curr = sister->get_Left_Corner();
	orig = curr;
	int number_nodes_B = 0;
	cout << "Cell b counting" << endl;
	do {
		number_nodes_B++;
		next = curr->get_Left_Neighbor();
		curr = next;
	} while (next != orig);
	
	sister->set_Wall_Count(number_nodes_B);
	
	cout << "updating angles" << endl;
	//update wall angles
	this->update_Wall_Angles();
	sister->update_Wall_Angles();
	
	cout << "updating center" << endl;
	//update cell center
	this->update_Cell_Center();
	sister->update_Cell_Center();

	cout << "update layer" << endl;
	//update layer information
	sister->set_Layer(this->layer);
	
	sister->set_growth_rate(this->growth_rate);
		
	cout << "Updated Angles and cell centers"<< endl;
	
	//distribute cyt nodes between sister cells
	cout << "deleting cyt nodes" << endl;
	int new_cyt_cnt = 20;
	//delete all old cyt nodes
	Cyt_Node* c = NULL;
	while(!cyt_nodes.empty()) {
		c= cyt_nodes.at(cyt_nodes.size() -1);
		delete c;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}	
	
	cout << "Finished deleting old cyt nodes" << endl;
	
	cout << "get most up/down and left/right for radius" << endl;
	this->most_Up_Down();
	sister->most_Up_Down();
	this->most_Left_Right();
	sister->most_Left_Right();

	double radius_x = ((most_right->get_Location() - most_left->get_Location()).length())*0.5;
	double radius_x_s = ((sister->get_most_right()->get_Location() - sister->get_most_left()->get_Location()).length())*0.5;
	double radius_y = ((most_up->get_Location() - most_down->get_Location()).length())*0.5;
	double radius_y_s = ((sister->get_most_up()->get_Location() - sister->get_most_down()->get_Location()).length())*0.5;
	//create new ones for each cell
	for(int i = 0; i < new_cyt_cnt; i++) {
		this->add_Cyt_Node_Div(radius_x, radius_y);
		sister->add_Cyt_Node_Div(radius_x_s, radius_y_s);
	}
	cout << "Cell A" << endl;
	/*vector<Cyt_Node*>CellA;
	this->get_Cyt_Nodes(CellA);
	int counter = 0;
	for(int i = 0;i<CellA.size();i++) {
		if(CellA.at(i)->get_My_Cell() == this) {
			counter++;
		}
	}*/
	//cout << "Number cyt nodes assigned to Cell A is: " << counter << endl;
	
	//cout << "Cell B" << endl;
/*	vector<Cyt_Node*>CellB;
	sister->get_Cyt_Nodes(CellB);
	counter = 0;
	for(int i = 0;i<CellB.size();i++) {
		if(CellB.at(i)->get_My_Cell() == sister) {
			counter++;
		}
	}*/
	//cout << "Number cyt nodes assigned to Cell B is: " << counter << endl;
		
	return sister;
}

/*
Cell* Cell::divide_width_wise(const int Ti) {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as the left sister
	//	this functoin will create a sister cell to the top
	//	and return it to the tissue
	//  cout << "Number Cells : " << my_tissue->get_Num_Cells() << endl;
//	int new_Rank = my_tissue->get_Num_Cells();
    //cout << new_Rank << " is new rank " << endl;
	Cell* sister = new Cell(Ti, my_tissue);
//	sister->set_Rank(new_Rank);
	vector<Side*> sister_sides;
	//Find division point
	//divide at midpoint
	Coord mid_left = ((sides.at(3)->get_End_A()->get_Location()+ sides.at(3)->get_End_Z()->get_Location())*0.5); 
	Coord mid_right = ((sides.at(1)->get_End_A()->get_Location()+ sides.at(1)->get_End_Z()->get_Location())*0.5); 
	
	Coord divider = ((mid_left + mid_right)/2);
	double divide_Y = divider.get_Y();

	//cout << "Begin Iterating Through Wall Nodes" << endl;

	Wall_Node* curr = sides.at(1)->get_End_A();
	Wall_Node* next = NULL;
	Wall_Node* prev = NULL;
	bool passed_div_line = false;
	int num_nodes = 0;
	Wall_Node* A1_End_Z = NULL;
	Wall_Node* B1_End_A = NULL;
		
	//make sides A1 and B1	
	do {
		next = curr->get_Left_Neighbor();
		prev = curr->get_Right_Neighbor();
		if(curr->get_Location().get_Y() > divide_Y) {
			//Cell B
			A1_End_Z = prev;
			B1_End_A = curr;
			passed_div_line = true;
		}
		curr = next;
		num_nodes++;
	} while( !passed_div_line );
	Side* A1 = new Side(sides.at(1)->get_End_A()->get_Location(),A1_End_Z->get_Location()-Coord(0,0.2),this, num_nodes-1);
	A1->set_Phys_Parameters(kBendHigh,kLinearLow);
	A1->set_Side_Type(1);
	Side* B1 = new Side(B1_End_A->get_Location()+Coord(0,0.2), sides.at(1)->get_End_Z()->get_Location(), sister, sides.at(1)->get_Wall_Count()-num_nodes-1);
	B1->set_Phys_Parameters(kBendHigh,kLinearLow);
	B1->set_Side_Type(1);
	//delete previous A1
	//cout << "deleting side 1" << endl;
	delete sides.at(1);
	//set new A1 and B1
	sides.at(1) = A1;

	//reassign A2 to B2
	Side* B2 = this->sides.at(2);
	B2->set_My_Cell(sister);
	
	//make sides A3 and B3
	curr = sides.at(3)->get_End_A();
	passed_div_line = false;
	Wall_Node* A3_End_A = NULL;
	Wall_Node* B3_End_Z = NULL;
			
	num_nodes = 0;
	//make sides A3 and B3
	do {
		next = curr->get_Left_Neighbor();
		prev = curr->get_Right_Neighbor();
		if(curr->get_Location().get_Y() < divide_Y) {
			//Cell A
			A3_End_A = curr;
			B3_End_Z = prev;
			passed_div_line = true;
		}
		curr = next;
		num_nodes++;
	} while(!passed_div_line);
	
	Side* A3 = new Side(A3_End_A->get_Location()-Coord(0,0.2),sides.at(3)->get_End_Z()->get_Location(),this, num_nodes-1);
	A3->set_Phys_Parameters(kBendHigh,kLinearLow);
	A3->set_Side_Type(3);
	Side* B3 = new Side(sides.at(3)->get_End_A()->get_Location(),B3_End_Z->get_Location()+Coord(0,0.2), sister, sides.at(3)->get_Wall_Count()-num_nodes-1);
	B3->set_Phys_Parameters(kBendHigh,kLinearLow);
	B3->set_Side_Type(3);
	//delete current A3
	delete sides.at(3);
	//set new A3
	sides.at(3) = A3;

	//make new cell wall which is A2 and B0
	Side* A2 = new Side(sides.at(1)->get_End_Z()->get_Location()+Coord(-0.04,0.04),sides.at(3)->get_End_A()->get_Location() + Coord(0.04,0.04), this, sides.at(0)->get_Wall_Count());
	Side* B0 = new Side(B3->get_End_Z()->get_Location()+ Coord(.04,-.04), B1->get_End_A()->get_Location()-Coord(0.04,0.04), sister, sides.at(0)->get_Wall_Count());
	A2->set_Phys_Parameters(kBendLow,kLinearHigh);
	A2->set_Side_Type(2);
	B0->set_Phys_Parameters(kBendLow,kLinearHigh);
	B0->set_Side_Type(0);
	//Assign A0 and all B
	sides.at(2) = A2;
	sister_sides.push_back(B0);
	sister_sides.push_back(B1);
	sister_sides.push_back(B2);
	sister_sides.push_back(B3);
	//cout<< "reassign sides in a" << endl;
	//connect the sides cell A
	this->sides.at(0)->connect_Ends(sides.at(1));
	this->sides.at(1)->connect_Ends(sides.at(2));
	this->sides.at(2)->connect_Ends(sides.at(3));
	this->sides.at(3)->connect_Ends(sides.at(0));  

	//cout << "reassign sides in b" << endl;
	//connect the sides cell B	
	sister_sides.at(0)->connect_Ends(sister_sides.at(1));
	sister_sides.at(1)->connect_Ends(sister_sides.at(2));
	sister_sides.at(2)->connect_Ends(sister_sides.at(3));
	sister_sides.at(3)->connect_Ends(sister_sides.at(0));
	
	//changes sides vector in cell B aka sister
	sister->set_Sides(sister_sides);

	int layer = this->get_Layer();
	sister->set_Layer(layer);
	//cout << "updating angles" << endl;
	//update angles and cell centers for both cells
	//cout << "cell A update wall angles" << endl;
	this->update_Wall_Angles();
	//cout << "cell B update wall angles" << endl;
	sister->update_Wall_Angles();
	//cout << "update centers" << endl;
	this->update_Cell_Center();
	sister->update_Cell_Center();
	
	double rate = (-.25*sister->get_Cell_Center().length() + 11.7)*2000;
	sister->set_growth_rate(rate);
	//cout << "Updated Angles and cell centers"<< endl;
	
	//distribute cyt nodes between sister cells
	int new_cyt_cnt = 20;
	//delete all old cyt nodes
	Cyt_Node* c = NULL;
	while(!cyt_nodes.empty()) {
		c= cyt_nodes.at(cyt_nodes.size() -1);
		delete c;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}	

	//cout << "Finished deleting old cyt nodes" << endl;
	//create new ones for each cell
	for(int i = 0; i < new_cyt_cnt; i++) {
		this->add_Cyt_Node_Div();
		sister->add_Cyt_Node_Div();
	}
	//cout << "Cell A" << endl;
	vector<Cyt_Node*>CellA;
	this->get_Cyt_Nodes(CellA);
	int counter = 0;
	for(int i = 0;i<CellA.size();i++) {
		if(CellA.at(i)->get_My_Cell() == this) {
			counter++;
		}
	}
	//cout << "Number cyt nodes assigned to Cell A is: " << counter << endl;
	
	//cout << "Cell B" << endl;
	vector<Cyt_Node*>CellB;
	sister->get_Cyt_Nodes(CellB);
	counter = 0;
	for(int i = 0;i<CellB.size();i++) {
		if(CellB.at(i)->get_My_Cell() == sister) {
			counter++;
		}
	}
	//cout << "Number cyt nodes assigned to Cell B is: " << counter << endl;
		
	return sister;
}


*/
