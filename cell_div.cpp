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
#include "side.h"
#include "tissue.h"
//==========================


Cell* Cell::divide(const int Ti) {
	Cell* sister = NULL;
	//calculate area
	double area = 0;
	Coord x_vec;
	Coord y_vec;
	double x_length;
	double y_length;
	x_vec = sides.at(0)->get_End_Z()->get_Location() - sides.at(0)->get_End_A()->get_Location();
	y_vec = sides.at(3)->get_End_A()->get_Location() - sides.at(3)->get_End_Z()->get_Location();
	x_length = x_vec.length();
	y_length = y_vec.length();
	area = y_length*x_length;
	if(area > 4) {
		if((this->layer == 3)|| (this->layer == 4)|| (this->layer == 5)) {
			//cout << "cell will divide" << endl;
			sister = this->divide_width_wise(Ti);		
		}
		else {
			sister = this->divide_length_wise(Ti);
		}
	}
	return sister;

}

Cell* Cell::divide_length_wise(const int Ti) {
	//current cell will split into two daughter cells
	//	-"this" will keep its entity as the left sister
	//	this functoin will create a sister cell to the right
	//	and return it to the tissue
//  cout << "Number Cells : " << my_tissue->get_Num_Cells() << endl;
//	int new_Rank = my_tissue->get_Num_Cells();
//	cout << new_Rank << " is new rank " << endl;
	Cell* sister = new Cell(Ti, my_tissue);
//	sister->set_Rank(new_Rank);
	vector<Side*> sister_sides;
	//Find division point
	//divide at midpoint
	Coord mid_top = ((sides.at(2)->get_End_A()->get_Location()+ sides.at(2)->get_End_Z()->get_Location())*0.5); 
	Coord mid_bot = ((sides.at(0)->get_End_A()->get_Location()+ sides.at(0)->get_End_Z()->get_Location())*0.5); 
	
	Coord divider = ((mid_top + mid_bot)/2);
	double divide_X = divider.get_X();

	cout << "Begin Iterating Through Wall Nodes" << endl;

	Wall_Node* curr = sides.at(0)->get_End_A();
	Wall_Node* next = NULL;
	Wall_Node* prev = NULL;
	bool passed_div_line = false;
	int num_nodes = 0;
	Wall_Node* A0_End_Z = NULL;
	Wall_Node* B0_End_A = curr;
		
	//make sides A0 and B0	
	do {
		next = curr->get_Left_Neighbor();
		prev = curr->get_Right_Neighbor();
		if(curr->get_Location().get_X() > divide_X) {
			//Cell B
			A0_End_Z = prev;
			B0_End_A = curr;
			passed_div_line = true;
		}
		curr = next;
		num_nodes++;
	} while( !passed_div_line );
	Side* A0 = new Side(sides.at(0)->get_End_A()->get_Location(),A0_End_Z->get_Location()-Coord(.2,0),this, num_nodes-1);
	A0->set_Phys_Parameters(kBendLow,kLinearHigh);
	A0->set_Side_Type(0);
	Side* B0 = new Side(B0_End_A->get_Location()+Coord(.2,0), sides.at(0)->get_End_Z()->get_Location(), sister, sides.at(0)->get_Wall_Count()-num_nodes-1);
	B0->set_Phys_Parameters(kBendLow,kLinearHigh);
	B0->set_Side_Type(0);
	//delete previous A0
	//cout << "deleting side 0" << endl;
	delete sides.at(0);
	//set new A0 and B0
	sides.at(0) = A0;
	sister_sides.push_back(B0);

	curr = sides.at(2)->get_End_A();
	passed_div_line = false;
	Wall_Node* A2_End_A = NULL;
	Wall_Node* B2_End_Z = NULL;
			
	num_nodes = 0;
	//make sides A1 and B3
	do {
		next = curr->get_Left_Neighbor();
		prev = curr->get_Right_Neighbor();
		if(curr->get_Location().get_X() < divide_X) {
			//Cell A
			A2_End_A = curr;
			B2_End_Z = prev;
			passed_div_line = true;
		}
		curr = next;
		num_nodes++;
	} while(!passed_div_line);
	
	Side* A2 = new Side(A2_End_A->get_Location()-Coord(.2,0),sides.at(2)->get_End_Z()->get_Location(),this, num_nodes-1);
	A2->set_Phys_Parameters(kBendLow,kLinearHigh);
	A2->set_Side_Type(2);

	Side* B2 = new Side(sides.at(2)->get_End_A()->get_Location(),B2_End_Z->get_Location()+Coord(.2,0), sister, sides.at(2)->get_Wall_Count()-num_nodes-1);
	B2->set_Phys_Parameters(kBendLow,kLinearHigh);
	B2->set_Side_Type(2);
	//delete current A2
	delete sides.at(2);
	//set new A2
	sides.at(2) = A2;
	//reassign A1 to B1
	Side* B1 = this->sides.at(1);
	B1->set_My_Cell(sister);
	sister_sides.push_back(B1);
	//Assign B2
	sister_sides.push_back(B2);
	//make new cell wall which is A1 and B3
	Side* A1 = new Side(sides.at(0)->get_End_Z()->get_Location()+Coord(0.04,0.04),sides.at(2)->get_End_A()->get_Location() + Coord(0.04,-0.04), this, sides.at(3)->get_Wall_Count());
	Side* B3 = new Side(sister_sides.at(2)->get_End_Z()->get_Location() - Coord(.04,.04), sister_sides.at(0)->get_End_A()->get_Location() + Coord(-0.04,0.04), sister, sister_sides.at(1)->get_Wall_Count());
	A1->set_Phys_Parameters(kBendHigh,kLinearLow);
	A1->set_Side_Type(1);
	B3->set_Phys_Parameters(kBendHigh,kLinearLow);
	B3->set_Side_Type(3);
	//Assign A1 and B3
	sides.at(1) = A1;
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



