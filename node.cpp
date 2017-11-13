//node.cpp
//=========================
#include <iostream>
#include <vector>
#include <cmath>
//=========================
#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
//=========================

//========================================
/** class Node Functions **/
Node::Node(Coord loc) {
	my_loc = loc;
	new_force = Coord();
}

Coord Node::get_Location() {
    return my_loc;
}

Coord Node::get_Force() {
	return new_force;
}

void Node::update_Location() {
	my_loc += new_force * dt;
    return;
}

Node::~Node() {}
//========================================
/**class Cyt Node Functions**/
Cyt_Node::Cyt_Node(Coord loc, Cell* my_cell) : Node(loc) {
	this->my_cell = my_cell;

}

void Cyt_Node::calc_Forces() {
	//for cytoplasm, just need morse potential for int-int and int-membr
	Coord Fii = calc_Morse_II();
	Coord Fmi = calc_Morse_MI(my_cell->get_Wall_Nodes());
    new_force = Fmi + Fii;
	
	return;
}

// Needs to have access:
//		-all the other cyt nodes of cell
//		-all the membr nodes of cell

Coord Cyt_Node::calc_Morse_II() {
	//calc force for II
	Coord Fii; //initialized to zero

	vector<Cyt_Node*>cyts;
	if (my_cell==NULL) {
		cout << "Error: Trying to access NULL Pointer. Aborting!" << endl;
		exit(1);
	}
	my_cell->get_Cyt_Nodes(cyts);

	for (unsigned int j = 0; j < cyts.size(); j++) {
		//don't calculate yourself
		if (cyts.at(j) != this) {
			//calc morse between this node and node j
			Fii += morse_Equation(cyts.at(j));
		}
	}

	return Fii;
}

Coord Cyt_Node::calc_Morse_MI(Wall_Node* orig) {
	//calc force for IM
	Coord Fmi;

	Wall_Node* curr = orig;
	
	do {
		Fmi += morse_Equation(curr);
		//update curr_wall
		curr = curr->get_Left_Neighbor();

	} while (curr != orig); 

	return Fmi;
}

Coord Cyt_Node::morse_Equation(Cyt_Node* cyt) {
	
	if (cyt == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
		exit(1);
	}
	
	//use Int-Int variables
    Coord Fii;
    Coord diff_vect = cyt->get_Location() - my_loc;
    double diff_len = diff_vect.length();
    double attract = (U_II/xsi_II)*exp(diff_len*(-1)/xsi_II);
    double repel = (W_II/gamma_II)*exp(diff_len*(-1)/gamma_II);
    
	Fii = diff_vect*((-attract + repel)/diff_len);
	return Fii;
}

Coord Cyt_Node::morse_Equation(Wall_Node* wall) {

	if (wall == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
		exit(1);
	} 

	//use Mem-Int variables
	Coord Fmi;
	Coord diff_vect = wall->get_Location() - my_loc; 
	double diff_len = diff_vect.length();
	double attract = (U_MI/xsi_MI)*exp(diff_len*(-1)/xsi_MI);
    double repel = (W_MI/gamma_MI)*exp(diff_len*(-1)/gamma_MI);
    
    Fmi = diff_vect*((-attract + repel)/diff_len);
	return Fmi;
}

Cyt_Node::~Cyt_Node() {
	my_cell = NULL;
}


//======================================================
/** class Wall Node Functions **/

// Constructors-----------------
Wall_Node::Wall_Node(Coord loc, Cell* my_cell) : Node(loc) {
	this->my_cell = my_cell;
}

Wall_Node::Wall_Node(Coord loc, Cell* my_cell, Wall_Node* left, Wall_Node* right) : Node(loc)   {
	this->left = left;
    this->right = right;
	this-> my_cell = my_cell;
	this-> pull = false;
	this->F_ext = Coord(0,0);
	update_Angle();
}

Wall_Node::~Wall_Node() {
	my_cell = NULL;	
	left = NULL;
	right = NULL;
	//closest = NULL;
}

//  Getters and Setters--------------------
void Wall_Node::set_Equi_Angle(double angle) {
	equi_angle = angle;
	return;
}

void Wall_Node::set_Left_Neighbor(Wall_Node* new_Left) {
	this->left = new_Left;
	return;
}

void Wall_Node::set_Right_Neighbor(Wall_Node* new_Right) {
	this->right = new_Right;
	return;
}

void Wall_Node::update_Angle() {
	Coord left_vect = get_Left_Neighbor()->get_Location() - get_Location();
	Coord right_vect = get_Right_Neighbor()->get_Location() - get_Location();

	
	double left_len = left_vect.length();
	double right_len = right_vect.length();

	double costheta = left_vect.dot(right_vect) / (left_len * right_len);
	double theta = acos( min( max(costheta,-1.0), 1.0) );

	double crossProd = left_vect.cross(right_vect);

	if (crossProd < 0.0) {
		theta = 2 * pi - theta;
	}
	
	//update protected member variables
	my_angle = theta;
	cross_Prod = crossProd;

	return;
}

void Wall_Node::update_Equi_Angle(double new_theta) {
	equi_angle = new_theta;

	return;
}

void Wall_Node::set_Closest(Wall_Node*  closest, double closest_len) {
	this->closest = closest;
	this->closest_len = closest_len;
	return;
}

// Calc Force Functions -----------------------
void Wall_Node::calc_Forces() {
	// Initialize force sum to zero by default constructor
	Coord sum;
	
	sum += calc_Morse_SC();
//	cout << "SC success" << endl;	
	cyt_force = sum;

	sum += calc_Morse_DC();
	//cout << "DC" << calc_Morse_DC() << endl;
	sum += calc_Linear();
//	cout << "linear" << endl;
	sum += calc_Bending();
//	cout << "bending" << endl;
	
	if(pull == true) {
		//cout << "External Force is" << calc_External() << endl;
		sum+= calc_External();
		pull = false;
	}
	// Update new_force variable for location updating
	new_force = sum;

	return;
}

//morse potential between wall node i and every cyt node in cell
Coord Wall_Node::calc_Morse_SC() {
	vector<Cyt_Node*> cyt_nodes;
	my_cell->get_Cyt_Nodes(cyt_nodes);

	Coord Fmi;
	
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Fmi += this->morse_Equation(cyt_nodes.at(i));
	}
	//cout << "	morse_sc: " << Fmi << endl;	
	return Fmi;
}


//morse potential between wall node i and every cyt node in cell
Coord Wall_Node::calc_Morse_DC() {
	Coord Fdc;
	vector<Cell*> cells;
	my_cell->get_Neighbor_Cells(cells);	
	Wall_Node* curr = NULL;
	Wall_Node* orig = NULL;
	for (unsigned int i = 0; i < cells.size(); i++) {
		curr = cells.at(i)->get_Wall_Nodes();
		orig = curr;
		do {
			Fdc += morse_Equation(curr);
			curr = curr->get_Left_Neighbor();
		} while (curr != orig);
	}

	if(this->closest != NULL){
		Fdc += linear_Equation_ADH(this->closest);
	}
	
	return Fdc;
}


//bending force of node
Coord Wall_Node::calc_Bending() {
	Coord F_bend;

	F_bend += bending_Equation_Center();
	F_bend += bending_Equation_Left();
	F_bend += bending_Equation_Right();
	
	if (cross_Prod < 0.0) {
		F_bend = F_bend*(-1);
	}	
	//cout << "	bending: " << F_bend << endl;
	return F_bend;
}

//spring force of neighboring springs
Coord Wall_Node::calc_Linear() {
	Coord F_lin;

//	cout << "calc left" << endl;
	F_lin += linear_Equation(left);

//	cout << "calc right" << endl;
	F_lin += linear_Equation(right);
	
//	cout << "	linear: " << F_lin << endl;
	return F_lin;
}

Coord Wall_Node::calc_External() {

	F_ext += Coord(EXTERNAL_FORCE*dt,0);
	
	if(this->get_Location().get_X() > my_cell->get_Cell_Center().get_X()) {
		return F_ext;
	}
	else if(this->get_Location().get_X() < my_cell->get_Cell_Center().get_X()) {
		return F_ext*-1;
	}
}

void Wall_Node::pull_node() {
	this->pull = true;
	return;
}


//===========================================================
// Mathematical force calculations


Coord Wall_Node::morse_Equation(Cyt_Node* cyt) {
	if (cyt == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
	}
	
	//use Membr - int variables
	Coord Fmi;
	Coord diff_vect = cyt->get_Location() - my_loc;
	double diff_len = diff_vect.length();
	double attract = (U_MI/xsi_MI)*exp(diff_len*(-1)/xsi_MI);
	double repel = (W_MI/gamma_MI)*exp(diff_len*(-1)/gamma_MI);

	Fmi = diff_vect * ((-attract + repel) / diff_len);

	return Fmi;
}

Coord Wall_Node::morse_Equation(Wall_Node* wall) {
	if (wall == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
	}

	//use Mem-Mem variables
    Coord Fmmd;
    Coord diff_vect = wall->get_Location() - my_loc;
    double diff_len = diff_vect.length();
    double attract = (U_MM/xsi_MM)*exp(diff_len*(-1)/xsi_MM);
    double repel = (W_MM/gamma_MM)*exp(diff_len*(-1)/gamma_MM);
    
    Fmmd = diff_vect * ((-attract + repel) / diff_len);
	return Fmmd;
}



Coord Wall_Node::bending_Equation_Center() {
	Coord F_center;
	double self_Constant; 
	
	double eps = 0.0001;

	if (abs(my_angle - pi) < eps) {
		return F_center;
	}
	else {
		self_Constant = K_BEND*(my_angle - equi_angle)/(sqrt(1-pow(cos(my_angle),2)));
	}

	Coord left_vect = left->get_Location() - my_loc;
	Coord right_vect = right->get_Location() - my_loc;
	double left_len = left_vect.length();
	double right_len = right_vect.length();
	Coord term_l1 = (left_vect*(-1))/(left_len*right_len);
	Coord term_l2 = left_vect*cos(my_angle)/pow(left_len,2);
	Coord term_r1 = (right_vect*(-1))/(left_len*right_len);
	Coord term_r2 = right_vect*cos(my_angle)/pow(right_len,2);

	F_center = (term_l1 + term_l2 + term_r1 + term_r2) * self_Constant;
	
	//cout << "Bending center: " << F_center << endl;	
	return F_center;
}

Coord Wall_Node::bending_Equation_Left() {
	Coord F_left;
	//double left_k_bend = left->get_Bending_Spring();
	double left_equi_angle = left->get_Equi_Angle();
	double left_angle = left->get_Angle();
	double left_Constant;
	
	double eps = 0.0001;

	if (abs(left_angle - pi) < eps) {
		return F_left;
	}
	else {
		left_Constant = K_BEND*(left_angle - left_equi_angle)/(sqrt(1-pow(cos(left_angle),2)));
	}

	Coord left_vect = left->get_Location() - my_loc;
	double left_len = left_vect.length();
	Coord left_left_vect = left->get_Left_Neighbor()->get_Location()-left->get_Location();
	double left_left_len = left_left_vect.length();
	Coord left_left_term1 = left_left_vect/(left_left_len*left_len);
	Coord left_term2 = left_vect*cos(left_angle)/pow(left_len,2);

	F_left = (left_left_term1 + left_term2) * left_Constant;
	
	//cout << "Bending left: " << F_left << endl;
	return F_left;
}

Coord Wall_Node::bending_Equation_Right() {
	Coord F_right;
	//double right_k_bend = right->get_Bending_Spring();
	double right_equ_angle = right->get_Equi_Angle();
	double right_angle = right->get_Angle();
	double right_Constant;
	
	double eps = 0.0001;

	if (abs(right_angle - pi) < eps) {
		return F_right;
	}
	else {
		right_Constant = K_BEND*(right_angle-right_equ_angle)/(sqrt(1-pow(cos(right_angle),2)));
	}

	Coord right_vect = right->get_Location() - my_loc;
	double right_len = right_vect.length();
	Coord right_right_vect = right->get_Right_Neighbor()->get_Location()-right->get_Location();
	double right_right_len = right_right_vect.length();
	Coord right_right_term1 = right_right_vect/(right_right_len*right_len);
	Coord right_term2 = right_vect*cos(right_angle)/pow(right_len,2);

	F_right = (right_right_term1 + right_term2)*right_Constant;
	
	//cout << "Bending right: " << F_right << endl;
	return F_right;
}

Coord Wall_Node::linear_Equation(Wall_Node* wall) {
	if (wall == NULL) {
		cout << "ERROR: Trying to access NULL pointer. Aborting!" << endl;
	}
	
	//use spring constant variables
	Coord F_lin;
	Coord diff_vect = wall->get_Location() - my_loc;
	double diff_len = diff_vect.length();
	Coord scaled_k = K_LINEAR*(diff_len - MembrEquLen);
	F_lin = (diff_vect/diff_len).distribute(scaled_k);

	return F_lin;	
}

Coord Wall_Node::linear_Equation_ADH(Wall_Node* wall) {
	if (wall == NULL) {
		cout << "Problems for days" << endl;
	}
	Coord F_lin;
	Coord diff_vect = wall->get_Location() - my_loc;
	double diff_len = diff_vect.length();
	F_lin = (diff_vect/diff_len)*(K_ADH*(diff_len - MembrEquLen_ADH));

	return F_lin;
}

//==========================================================
//Adhesion functions
Wall_Node* Wall_Node::find_Closest_Node(vector<Cell*>& neighbors) {
	Wall_Node* curr = NULL;
	Wall_Node* orig = NULL;
	Wall_Node* next = NULL;
	Cell* curr_cell = NULL;
	Wall_Node* closest = NULL;
	double curr_dist = 0;
	double smallest = 100;
	for(int i = 0; i < neighbors.size(); i++) {
		curr_cell = neighbors.at(i);
		//find the closest node on curr_Side
		curr = curr_cell->get_Left_Corner();
		orig = curr;
		do{
			next = curr->get_Left_Neighbor();
			curr_dist = (this->my_loc - curr->get_Location()).length();
			if(curr_dist < ADHThresh) {
				if(curr_dist < smallest) {
					closest = curr;
					smallest = curr_dist;
				}
			}
			curr = next;
		} while (next != orig);
	}
	return closest;
}

Wall_Node* Wall_Node::find_Closest_Node_Beg(vector<Cell*>& neighbors) {
	Wall_Node* curr = NULL;
	Wall_Node* orig = NULL;
	Wall_Node* next = NULL;
	Cell* curr_cell = NULL;
	Wall_Node* closest = NULL;
	double curr_dist = 0;
	double smallest = 100;
	for(int i = 0; i < neighbors.size(); i++) {
		curr_cell = neighbors.at(i);
		//find the closest node on curr_Side
		curr = curr_cell->get_Left_Corner();
		orig = curr;
		do{
			next = curr->get_Left_Neighbor();
			curr_dist = (this->my_loc - curr->get_Location()).length();
			if(curr_dist < ADHThreshBeg) {
				if(curr_dist < smallest) {
					closest = curr;
					smallest = curr_dist;
				}
			}
			curr = next;
		} while (next != orig);
	}
	return closest;
}

void Wall_Node::make_Connection(Wall_Node* curr_Closest) {
	double curr_dist = 0;
	if(curr_Closest != NULL) {
		//cout << "Making connection" << endl;
		if(curr_Closest->get_Closest() != NULL) {
			curr_dist = (this->get_Location() - curr_Closest->get_Location()).length();
			if(curr_dist < this->closest_len) {
				this->closest_len = curr_dist;
				this->closest = curr_Closest;
				curr_Closest->set_Closest(this, curr_dist);
			}
			else if (curr_dist > this->closest_len) {
				//do nothing
			}
		}

		else {
			this->closest_len = curr_dist;
			this->closest = curr_Closest;
			curr_Closest->set_Closest(this, curr_dist);
		}
	}
}
//==========================================================
// End of node.cpp

