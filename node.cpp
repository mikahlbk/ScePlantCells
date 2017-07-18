//node.cpp
#include "coord.h"
#include "node.h"

/** class Node Functions **/
Node::Node(Coord loc) {
	my_loc = loc;
	new_Force = Coord();
}

Coord Node::get_location() {
    return my_loc;
}

void Node::update_Location() {
    my_loc += new_force*dt;
    return;
}


/** class Cyt Node Functions **/
Cyt_Node::Cyt_Node(Coord loc) : Node(loc) {};

void Cyt_Node::calc_Forces(Cell* my_cell) {
	//for cyt, just need morse potential for int-int and int-membr
	Coord Fii = calc_Morse_II(my_cell->get_CytNodes());

	Coord Fmi = calc_Morse_MI(my_cell->get_WallNodes());
	
    new_Force = Fmi + Fii;
}

// Needs to have access:
//		-all the other cyt nodes of cell
//		-all the membr nodes of cell
Coord Cyt_Node::calc_Morse_II(const vector<Cyt_Node*>& cyt_nodes) {
	//calc force for II
	Coord Fii; //need to initialize to zero

	for (int j = 0; j < cyt_nodes.size(); j++) {
		//don't calculate yourself
		if (cyt_nodes.at(j) != this) {
			//calc morse between this node and node j
			Fii += morse_Equation(cyt_nodes.at(i));
		}
	}

	return Fii;
}

Coord Cyt_Node::calc_Morse_MI(Wall_Node* orig)
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
    Coord diff_vect = my_loc - cyt->get_Location();
    double diff_len = diff_vect.length();
    double attract = (U_II/xsi_II)*exp(diff_len*(-1)/xsi_II);
    double repel = (W_II/gamma_II)*exp(diff_len*(-1)/gamma_II);
    
	Fii = diff_vect*(-attract + repel)/diff_len;
	return Fii;

}

Coord Cyt_Node::morse_Equation(Wall_Node* wall) {
	//use Mem-Int variables
	Coord Fmi;
	Coord diff_vect = my_loc - cyt->get_Location();
	double diff_len = diff_vect.length();
	double attract = (U_MI/xsi_MI)*exp(diff_len*(-1)/xsi_MI);
    double repel = (W_MI/gamma_MI)*exp(diff_len*(-1)/gamma_MI);
    
    Fmi = diff_vect*(-attract + repel)/diff_len;
	return Fmi;
}



/** class Wall Node Functions **/
Wall_Node::Wall_Node(Coord loc) : Node(loc) {};

Wall_Node::Wall_Node(Coord loc, Wall_Node* left, Wall_Node* right, double angle) : Node(loc) {
    this->left = left;
    this->right = right;
	my_angle = angle;
}

double Wall_Node::get_Angle() {
	return my_angle;
}

Wall_Node* Wall_Node::get_Left_Neighbor() {
	return left;
}

Wall_Node* Wall_Node::get_Right_Neighbor() {
	return right;
}

void Wall_Node::calc_Forces(Cell* my_cell) {
	
	Coord sum;

	sum += calc_Morse_SC(my_cell->get_CytNodes);
	sum += calc_Morse_DC(my_cell->get_Neigh_Cells());
	sum += calc_Linear();
	sum += calc_Bending();

	new_Force = sum;;
}
//morse potential between wall node i and every cyt node in cell
Coord Wall_Node::calc_Morse_SC(vector<Cyt_Node*> cyt_nodes) 	Coord Fmi;
	
	for (int i = 0; i < cyt_nodes.size(); i++) {
		Fmi += this->morse_Equation(cyt_nodes.at(j));
	}
	
	return Fmi;
}
//probably need vector of relatively close cells
Coord Wall_Node::calc_Morse_DC(vector<Cell*>& cells) {
	Coord Fdc;
	//iterate through each cell
	for (int i = 0; i < cells.size(); i++) {
		//find which nodes from cell.at(i) you will need

		//iterate through membrane nodes of each cell
		for () {
			
			Fdc += this->morse_Function();

		}

	}

	return Fdc;
}

Coord Wall_Node::morse_Equation(Cyt_Node* cyt) {
	//use Membr - int variables
	Coord Fmi
	Coord diff_vect = my_loc - cyt->get_Location();
	double diff_len = diff_vect.length();
	double attract = (U_MI/xsi_MI)*exp(diff_len*(-1)/xsi_MI);
	double repel = (W_MI/gamma_MI)*exp(diff_len*(-1)/gamma_MI);
	
}

Coord Wall_Node::morse_Equation(Wall_Node* wall) {
	//use Int-Int variables
    Coord Fmmd;
    Coord diff_vect = my_loc - wall->get_Location();
    double diff_len = diff_vect.length();
    double attract = (U_MMD/xsi_MMD)*exp(diff_len*(-1)/xsi_MMD);
    double repel = (W_MMD/gamma_MMD)*exp(diff_len*(-1)/gamma_MMD);
    
    Fmmd = diff_vect*(-attract + repel)/diff_len;
}

Coord Wall_Node::linear_Equation(Wall_Node* wall, k_linear) {
	//use spring constant  variables
	Coord Flin;
	Coord diff_vect = wall->get_location - my_loc;
	double diff_len = diff_vect.length();
	
	Fmi = diff_vect*(k_linear*(diff_len-MembrEquLen)/diff_len);
	return Fmi;	
	
}

Coord Wall_Node::bending_Equation_Center() {
	Coord F_center;
	double k_bend = get_bendingSpring();
	double equ_angle = get_Equi_Angle();
	double self_Constant = k_bend*(my_angle - equ_angle)/sqrt(1-pow(cos(my_angle),2));
	Coord left_vect = left->get_Location()-my_loc;
	Coord right_vect = right->get_Location()-my_loc;
	double left_len = left_vect.length();
	double right_len = right_vect.length();
	Coord left_term1 = left_vect*(-1/(left_len*right_len));
	Coord left_term2 = left_vect*(cos(my_angle)/pow(left_len,2));
	Coord left_term3 = right_vect*(-1/(left_len*right_len));
	Coord left_term4 = right_vect*(cos(my_angle)/pow(right_len,2));

	F_center = (left_term1 + left_term2 + left_term3 + left_term4)*self_Constant;
	return F_center;
}

Coord Wall_Node::bending_Equation_Left() {
	Coord F_left;
	double left_k_bend = left->get_bendingSpring();
	double left_equ_angle = left->get_Equi_Angle();
	double left_angle = left->get_Angle();
	double left_Constant = k_bend*(left_angle - left_equ_angle)/sqrt(1-pow(cos(left_angle),2));
	Coord left_vect = left->get_Location() - my_loc;
	double left_len = left_vect.length();
	Coord left_left_vect = left->get_Left_Neighbor()->get_Location()-left->get_Location();
	double left_left_len = left_left_vect.length();
	Coord left_left_term1 = left_left_vect/(left_left_len*left_len);
	Coord left_term2 = left_vect*((-1)*cos(left_angle)/pow(left_len,2));

	F_left = (left_left_term1 + left_term2)*left_Constant;
	return F_left;
}

Coord Wall_Node::bending_Equation_right() {
	Coord F_right;
	double right_k_bend = right->get_bendingSpring();
	double right_equ_angle = rightt->get_Equi_Angle();
	double right_angle = right->get_Angle();
	double right_Constant = k_bend*(right_angle -right_equ_angle)/sqrt(1-pow(cos(right_angle),2));
	Coord right_vect = right->get_Location() - my_loc;
	double right_len = right_vect.length();
	Coord right_right_vect = right->get_Right_Neighbor()->get_Location()-right->get_Location();
	double right_right_len = right_right_vect.length();
	Coord right_right_term1 = right_right_vect/(right_right_len*right_len);
	Coord right_term2 = right_vect*((-1)*cos(right_angle)/pow(right_len,2));

	F_right = (right_right_term1 + right_term2)*right_Constant;
	return F_right;
}

Coord Wall_Node::calc_Bending() {
	Coord F_bend;
	F_center = bending_Equation_Center();
	F_left = bending_Equation_Left();
	F_right = bending_Equation_Right();

	F_bend = F_center + F_left + F_right;

	return F_bend;
}


/** class Corner Node Functions **/
Corner_Node::Corner_Node(Coord loc) : Wall_Node(loc) {};

Corner_Node::Corner_Node(Coord loc, Wall_Node* left, Wall_Node* right, double angle) 
    : Wall_Node(loc, left, right, angle) {}

double Corner_Node::get_Equi_Angle() {
	return thetaCorner;
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



/** class Flank Node function **/
Flank_Node::Flank_Node(Coord loc) : Wall_Node(loc) {};

Flank_Node::Flank_Node(Coord loc, Wall_Node* left, Wall_Node* right, double angle) 
	: Wall_Node(loc, left, right, angle) {}

Coord Flank_Node::calc_Linear() {
	//as a flank node, both springs on either side have flank constants

	//calc left spring force
	Coord F_left = this->linear_Equation(left, );

	//calc right spring force
	Coord F_rt = this->linear_Equation(right, );

	return F_left + F_rt;
}

double Flank_Node::get_Equi_Angle() {
	return thetaFlank;
}

double Flank_Node::get_linearSpring() {
	return kLinearFlank;
}

/** class Edge Node function **/
End_Node::End_Node(Coord loc) : Wall_Node(loc) {};

End_Node::End_Node(Coord loc, Node* left, Node* right, double angle)
    : Wall_Node(loc, left, right, angle) {}

Coord End_Node::calc_Linear() {
	//as end node, both springs on either sied have end constants

	//calc left spring force
	Coord F_left = this->linear_Equation(left, );

	//calc right spring force
	Coord F_rt = this->linear_Equation(right, );

	return F_left + F_rt;
}

double End_Node::get_Equi_Angle() {
	return thetaEnd;
}

double End_Node::get_linearSpring() {
	return kLinearEnd;
}





	
