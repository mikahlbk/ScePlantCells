//node.cpp
//=========================
#include <iostream>
#include <vector>
#include <math.h>
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

void Node::update_Location() {
    my_loc += new_force * dt;
    return;
}

//========================================
/** class Cyt Node Functions **/
Cyt_Node::Cyt_Node(Coord loc) : Node(loc) {};

void Cyt_Node::calc_Forces(Cell* my_cell) {
	//for cyt, just need morse potential for int-int and int-membr
	vector<Cyt_Node*> cyts;
	my_cell->get_CytNodes(cyts);

	Coord Fii = calc_Morse_II(cyts);

	Coord Fmi = calc_Morse_MI(my_cell->get_WallNodes());
	
    new_force = Fmi + Fii;

	return;
}

// Needs to have access:
//		-all the other cyt nodes of cell
//		-all the membr nodes of cell
Coord Cyt_Node::calc_Morse_II(vector<Cyt_Node*>& cyt_nodes) {
	//calc force for II
	Coord Fii; //need to initialize to zero

	for (unsigned int j = 0; j < cyt_nodes.size(); j++) {
		//don't calculate yourself
		if (cyt_nodes.at(j) != this) {
			//calc morse between this node and node j
			Fii += morse_Equation(cyt_nodes.at(j));
		}
	}

	return Fii;
}

Coord Cyt_Node::calc_Morse_MI(Wall_Node* orig) {
	//calc force for IM
	Coord Fmi;
	Wall_Node* curr = orig;
	
	do {
		Fmi = morse_Equation(curr);
		//update curr_wall
		curr = curr->get_Left_Neighbor();

	} while (curr != orig); 

	return Fmi;
}

Coord Cyt_Node::morse_Equation(Cyt_Node* cyt) {
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
	//use Mem-Int variables
	Coord Fmi;
	Coord diff_vect = wall->get_Location() - my_loc; 
	double diff_len = diff_vect.length();
	double attract = (U_MI/xsi_MI)*exp(diff_len*(-1)/xsi_MI);
    double repel = (W_MI/gamma_MI)*exp(diff_len*(-1)/gamma_MI);
    
    Fmi = diff_vect*((-attract + repel)/diff_len);
	return Fmi;
}


//======================================================
/** class Wall Node Functions **/

// Constructors-----------------
Wall_Node::Wall_Node(Coord loc) : Node(loc) {};

Wall_Node::Wall_Node(Coord loc, Wall_Node* left, Wall_Node* right) : Node(loc) {
    this->left = left;
    this->right = right;

	update_Angle();
}

//  Getters and Setters--------------------
double Wall_Node::get_Angle() {
	return my_angle;
}

Wall_Node* Wall_Node::get_Left_Neighbor() {
	return left;
}

Wall_Node* Wall_Node::get_Right_Neighbor() {
	return right;
}

void Wall_Node::set_Left_Neighbor(Wall_Node* new_Left){
	//change the pointer
	this->left = new_Left; 
}

void Wall_Node::set_Right_Neighbor(Wall_Node* new_Right) {
	//change the pointer
	this->right = new_Right;
}

// Calc Force Functions -----------------------
void Wall_Node::calc_Forces(Cell* my_cell) {
	// Initialize force sum to zero by default constructor
	Coord sum;
	// gather cyt nodes of your cell for morse calc
	vector<Cyt_Node*> cyts;
	my_cell->get_CytNodes(cyts);
//	sum += calc_Morse_SC(cyts);

	//will be implimented later
	//sum += calc_Morse_DC(my_cell->get_Neigh_Cells());
	
//	sum += calc_Linear();
	sum += calc_Bending();

	// Update new_force variable for location updating
	new_force = sum;;
}

void Wall_Node::update_Angle() {
	Coord left_vect = get_Left_Neighbor()->get_Location() - get_Location();
	Coord right_vect = get_Right_Neighbor()->get_Location() - get_Location();
	
	double left_len = left_vect.length();
	double right_len = right_vect.length();

	double costheta = left_vect.dot(right_vect) / (left_len * right_len);
	double theta = acos( min( max(costheta,-1.0), 1.0) );

	double crossProd = left_vect.cross(right_vect);


	if (crossProd < 0) {
		theta = 2 * pi - theta;
	}
	
	//update protected member variables
	my_angle = theta;
	cross_Prod = crossProd;

	return;
}

//morse potential between wall node i and every cyt node in cell
Coord Wall_Node::calc_Morse_SC(vector<Cyt_Node*>& cyt_nodes) {
	Coord Fmi;
	
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Fmi += this->morse_Equation(cyt_nodes.at(i));
	}
	
	return Fmi;
}

//probably need vector of relatively close cells
Coord Wall_Node::calc_Morse_DC(vector<Cell*>& cells) {
	Coord Fdc;
	/*
	//iterate through each cell
	for (int i = 0; i < cells.size(); i++) {
		//find which nodes from cell.at(i) you will need

		//iterate through membrane nodes of each cell
		for () {
			
			Fdc += this->morse_Function();

		}

	}
	*/
	return Fdc;
}

Coord Wall_Node::calc_Bending() {
	Coord F_bend;

	F_bend += bending_Equation_Center();
	F_bend += bending_Equation_Left();
	F_bend += bending_Equation_Right();
	
	if (cross_Prod < 0) {
		F_bend = F_bend*(-1);
	}	
    cout << F_bend << endl;
	
	return F_bend;
}

Coord Wall_Node::morse_Equation(Cyt_Node* cyt) {
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
	//use Int-Int variables
    Coord Fmmd;
    Coord diff_vect = wall->get_Location() -  my_loc;
    double diff_len = diff_vect.length();
    double attract = (U_MM/xsi_MM)*exp(diff_len*(-1)/xsi_MM);
    double repel = (W_MM/gamma_MM)*exp(diff_len*(-1)/gamma_MM);
    
    Fmmd = diff_vect * ((-attract + repel) / diff_len);

	return Fmmd;
}

Coord Wall_Node::linear_Equation(Wall_Node* wall, double k_linear) {
	//use spring constant variables
	Coord Flin;
	Coord diff_vect = wall->get_Location() - my_loc;
	double diff_len = diff_vect.length();
	
	Flin = (diff_vect/diff_len)*(k_linear * (diff_len - MembrEquLen));

	return Flin;	
}

Coord Wall_Node::bending_Equation_Center() {
	Coord F_center;
	double k_bend = get_bendingSpring();
	double equ_angle = get_Equi_Angle();
//	cout << my_angle << endl;
//	cout << equ_angle << endl;
	double self_Constant; 
	
//	double eps = 0.0001;

//	if (my_angle - pi < eps) {
//		return F_center;
//	}
//	else {
		self_Constant = k_bend*(my_angle - equ_angle)/(sqrt(1-pow(cos(my_angle),2)) + .001);
//	}

	Coord left_vect = left->get_Location() - my_loc;
	Coord right_vect = right->get_Location() - my_loc;
	double left_len = left_vect.length();
	double right_len = right_vect.length();
	Coord term_l1 = left_vect*(-1)/(left_len*right_len);
	Coord term_l2 = left_vect*(cos(my_angle))/pow(left_len,2);
	Coord term_r1 = right_vect*(-1)/(left_len*right_len);
	Coord term_r2 = right_vect*(cos(my_angle))/pow(right_len,2);

	F_center = (term_l1 + term_l2 + term_r1 + term_r2) * self_Constant;
	
	return F_center;
}

Coord Wall_Node::bending_Equation_Left() {
	Coord F_left;
	double left_k_bend = left->get_bendingSpring();
	double left_equ_angle = left->get_Equi_Angle();
//	cout << left_equ_angle << endl;
	double left_angle = left->get_Angle();
//	cout << left_angle << endl;
	double left_Constant;
	
//	double eps = 0.001;

//	if (left_angle - pi < eps) {
//		return F_left;
//	}
//	else {
		left_Constant = left_k_bend*(left_angle - left_equ_angle)/(sqrt(1-pow(cos(left_angle),2)) + .001);
//	}

	Coord left_vect = left->get_Location() - my_loc;
	double left_len = left_vect.length();
	Coord left_left_vect = left->get_Left_Neighbor()->get_Location()-left->get_Location();
	double left_left_len = left_left_vect.length();
	Coord left_left_term1 = left_left_vect/(left_left_len*left_len);
	Coord left_term2 = left_vect*cos(left_angle)/pow(left_len,2);

	F_left = (left_left_term1 + left_term2) * left_Constant;
	
	return F_left;
}

Coord Wall_Node::bending_Equation_Right() {
	Coord F_right;
	double right_k_bend = right->get_bendingSpring();
	double right_equ_angle = right->get_Equi_Angle();
//	cout << right_equ_angle << endl;
	double right_angle = right->get_Angle();
//	cout << right_angle << endl;
	double right_Constant;
	
//	double eps = 0.001;

//	if (right_angle - pi < eps) {
//		return F_right;
//	}
//	else {
		right_Constant = right_k_bend*(right_angle-right_equ_angle)/(sqrt(1-pow(cos(right_angle),2)) + .001);
//	}

	Coord right_vect = right->get_Location() - my_loc;
	double right_len = right_vect.length();
	Coord right_right_vect = right->get_Right_Neighbor()->get_Location()-right->get_Location();
	double right_right_len = right_right_vect.length();
	Coord right_right_term1 = right_right_vect/(right_right_len*right_len);
	Coord right_term2 = right_vect*cos(right_angle)/pow(right_len,2);

	F_right = (right_right_term1 + right_term2)*right_Constant;
	return F_right;
}



//========================================================
/** class Corner Node Functions **/
Corner_Node::Corner_Node(Coord loc) : Wall_Node(loc) {}

Corner_Node::Corner_Node(Coord loc, Wall_Node* left, Wall_Node* right) 
    : Wall_Node(loc, left, right) {}

double Corner_Node::get_Equi_Angle() {
	return thetaCorner;
}

//should never call this function
double Corner_Node::get_linearSpring() {
	return 100000;
}

double Corner_Node::get_bendingSpring() {
	return kBendEnd;
}

Coord Corner_Node::calc_Linear() {
	//as a corner, have to find out which spring is end and 
	//	which is flank

	//calc left spring force
	Coord F_left = this->linear_Equation(left, left->get_linearSpring());

	//calc right spring force
	Coord F_rt = this->linear_Equation(right, right->get_linearSpring());

	return F_left + F_rt;
}

bool Corner_Node::is_Corner() {
	return true;
}


//===========================================================
/** class Flank Node function **/
Flank_Node::Flank_Node(Coord loc) : Wall_Node(loc) {};

Flank_Node::Flank_Node(Coord loc, Wall_Node* left, Wall_Node* right) 
	: Wall_Node(loc, left, right) {}

double Flank_Node::get_Equi_Angle() {
	return thetaFlank;
}

Coord Flank_Node::calc_Linear() {
	//as a flank node, both springs on either side have flank constants

	//calc left spring force
	Coord F_left = this->linear_Equation(left, kLinearFlank);

	//calc right spring force
	Coord F_rt = this->linear_Equation(right, kLinearFlank);

	return F_left + F_rt;
}

double Flank_Node::get_linearSpring() {
	return kLinearFlank;
}

double Flank_Node::get_bendingSpring() {
	return kBendEnd;
}

bool Flank_Node::is_Corner() {
	return false;
}


//==============================================================
/** class End Node function **/
End_Node::End_Node(Coord loc) : Wall_Node(loc) {};

End_Node::End_Node(Coord loc, Wall_Node* left, Wall_Node* right)
    : Wall_Node(loc, left, right) {}

Coord End_Node::calc_Linear() {
	//as end node, both springs on either sied have end constants

	//calc left spring force
	Coord F_left = this->linear_Equation(left, kLinearEnd);

	//calc right spring force
	Coord F_rt = this->linear_Equation(right, kLinearEnd);

	return F_left + F_rt;
}

double End_Node::get_Equi_Angle() {
	return thetaEnd;
}

double End_Node::get_linearSpring() {
	return kLinearEnd;
}

double End_Node::get_bendingSpring() {
	return kBendEnd;
}

bool End_Node::is_Corner() {
	return false;
}

//==========================================================
// End of node.cpp

	
