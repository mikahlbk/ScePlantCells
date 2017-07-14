//node.cpp
#include "coord.h"
#include "node.h"

/** class Node Functions **/
Node::Node(Location loc) {
	my_loc = loc;
}

Location Node::get_location() {
    return my_loc;

}

Force Node::morse_Equation() {

}

/** class Cyt Node Functions **/
Cyt_Node::Cyt_Node(Location loc) : Node(loc) {};

Force Cyt_Node::calc_Forces(Cell* my_cell) {
	//for cyt, just need morse potential for int-int and int-membr
	Force Fii = calc_Morse_II(my_cell->get_CytNodes());

	Force Fmi = calc_Morse_MI(my_cell->get_WallNodes());
	
    return Fmi + Fii;
}

// Needs to have access:
//		-all the other cyt nodes of cell
//		-all the membr nodes of cell
Force Cyt_Node::calc_Morse_II(vector<Cyt_Node*>& cyt_nodes) {
	//calc force for II
	Force FII; //need to initialize to zero

	for (int j = 0; j < cyt_nodes.size(); j++) {
		//don't calculate yourself
		if (cyt_nodes.at(j) != this) {
			//calc morse between this node and node j
			FII += this->morse_Function(cyt_nodes.at(j), Uii, Wii, Zii, Gii);
		}
	}

	return FII;
}

Force Cyt_Node::calc_Morse_MI(Wall_Node* curr)
	//calc force for IM
	Force FMI;
	Wall_Node* orig = curr;
	
	do {
		FMI = this->morse_Function(curr, Umi, Wmi, Zmi, Gii);
		//update curr_wall
		curr = curr->get_Left_Neighbor();

	} while (curr != orig); 

	return FMI;
}



/** class Wall Node Functions **/
Wall_Node::Wall_Node(Location loc) : Node(loc) {};

Wall_Node::Wall_Node(Location loc, Wall_Node* left, Wall_Node* right) : Node(loc) {
    this->left = left;
    this->right = right;
}

double Wall_Node::get_Angle() {
	return my_angle;
}

Force Wall_Node::calc_Forces(Cell* my_cell) {
	
	Force sum;

	sum += calc_Morse_SC(my_cell->get_CytNodes);
	sum += calc_Morse_DC();
	sum += calc_Linear();
	sum += calc_Bending();

	return sum;
}
//morse potential between wall node i and every cyt node in cell
Force Wall_Node::calc_Morse_SC(vector<Cyt_Node*> cyt_nodes) {
	Force Fmi;
	
	for (int i = 0; i < cyt_nodes.size(); i++) {
		Fmi += this->morse_Function(cyt_nodes.at(j), Umi, Wmi, Zmi, Gii);
	}
	
	return Fmi;
}
//probably need vector of relatively close cells
Force Wall_Node::calc_Morse_DC(vector<Cell*>& cells) {
	Force Fdc;
	//iterate through each cell
	for (int i = 0; i < cells.size(); i++) {
		//iterate through membrane nodes of each cell
		for () {
			
			Fdc += this->morse_Function(   , Ummd, Wmmd, Zmmd, Gmmd);

		}

	}

	return Fdc;
}

Force Wall_Node::calc_Bending() {
	Force F_bend;

	return F_bend;
}
/** class Corner Node Functions **/
Corner_Node::Corner_Node(Location loc) : Wall_Node(loc) {};

Corner_Node::Corner_Node(Location loc, Wall_Node* left, Wall_Node* right, double angle) 
    : Wall_Node(loc, left, right, angle) {}

Force Corner_Node::calc_Forces() {
	Force sum;

	//calc Morse_SC
	sum += calc_Morse_SC();

	//calc Morse_DC
	//calc linear

	return sum;
}

Force Corner_Node::calc_Linear() {
	//as a corner, have to find out which spring is end and 
	//	which is flank

	//calc left spring force
	Force F_left = this->linear_Equation(left, left->get_linearSpring());

	//calc right spring force
	Force F_rt = this->linear_Equation(right, right->get_linearSpring());

	return F_left + F_rt;
}

/** class Flank Node function **/
Flank_Node::Flank_Node(Location loc) : Wall_Node(loc) {};

Flank_Node::Flank_Node(Location loc, Wall_Node* left, Wall_Node* right, double angle) 
	: Wall_Node(loc, left, right, angle) {}

Force Flank_Node::calc_Linear() {
	//as a flank node, both springs on either side have flank constants

	//calc left spring force
	Force F_left = this->linear_Equation(left, );

	//calc right spring force
	Force F_rt = this->linear_Equation(right, );

	return F_left + F_rt;
}

/** class Edge Node function **/
End_Node::End_Node(Location loc) : Wall_Node(loc) {};

End_Node::End_Node(Location loc, Node* left, Node* right, double angle)
    : Wall_Node(loc, left, right, angle) {}

Force End_Node::calc_Linear() {
	//as end node, both springs on either sied have end constants

	//calc left spring force
	Force F_left = this->linear_Equation(left, );

	//calc right spring force
	Force F_rt = this->linear_Equation(right, );

	return F_left + F_rt;
}



