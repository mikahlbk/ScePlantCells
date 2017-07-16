//node.cpp
#include "coord.h"
#include "node.h"

/** class Node Functions **/
Node::Node(Coord loc) {
	my_loc = loc;
}

Coord Node::get_location() {
    return my_loc;

}

Coord Node::morse_Equation() {

}

/** class Cyt Node Functions **/
Cyt_Node::Cyt_Node(Coord loc) : Node(loc) {};

Force Cyt_Node::calc_Forces(Cell* my_cell) {
	//for cyt, just need morse potential for int-int and int-membr
	Coord Fii = calc_Morse_II(my_cell->get_CytNodes());

	Coord Fmi = calc_Morse_MI(my_cell->get_WallNodes());
	
    return Fmi + Fii;
}

// Needs to have access:
//		-all the other cyt nodes of cell
//		-all the membr nodes of cell
Coord Cyt_Node::calc_Morse_II(vector<Cyt_Node*>& cyt_nodes) {
	//calc force for II
	Coord FII; //need to initialize to zero

	for (int j = 0; j < cyt_nodes.size(); j++) {
		//don't calculate yourself
		if (cyt_nodes.at(j) != this) {
			//calc morse between this node and node j
			FII += this->morse_Function(cyt_nodes.at(j), Uii, Wii, Zii, Gii);
		}
	}

	return FII;
}

Coord Cyt_Node::calc_Morse_MI(Wall_Node* curr)
	//calc force for IM
	Coord FMI;
	Wall_Node* orig = curr;
	
	do {
		FMI = this->morse_Function(curr, Umi, Wmi, Zmi, Gii);
		//update curr_wall
		curr = curr->get_Left_Neighbor();

	} while (curr != orig); 

	return FMI;
}



/** class Wall Node Functions **/
Wall_Node::Wall_Node(Coord loc) : Node(loc) {};

Wall_Node::Wall_Node(Coord loc, Wall_Node* left, Wall_Node* right) : Node(loc) {
    this->left = left;
    this->right = right;
}

double Wall_Node::get_Angle() {
	return my_angle;
}

Coord Wall_Node::calc_Forces(Cell* my_cell) {
	
	Coord sum;

	sum += calc_Morse_SC(my_cell->get_CytNodes);
	sum += calc_Morse_DC();
	sum += calc_Linear();
	sum += calc_Bending();

	return sum;
}
//morse potential between wall node i and every cyt node in cell
Coord Wall_Node::calc_Morse_SC(vector<Cyt_Node*> cyt_nodes) {
	Coord Fmi;
	
	for (int i = 0; i < cyt_nodes.size(); i++) {
		Fmi += this->morse_Function(cyt_nodes.at(j), Umi, Wmi, Zmi, Gii);
	}
	
	return Fmi;
}
//probably need vector of relatively close cells
Coord Wall_Node::calc_Morse_DC(vector<Cell*>& cells) {
	Coord Fdc;
	//iterate through each cell
	for (int i = 0; i < cells.size(); i++) {
		//iterate through membrane nodes of each cell
		for () {
			
			Fdc += this->morse_Function(   , Ummd, Wmmd, Zmmd, Gmmd);

		}

	}

	return Fdc;
}

Coord Wall_Node::calc_Bending() {
	Coord F_bend;

	return F_bend;
}
/** class Corner Node Functions **/
Corner_Node::Corner_Node(Coord loc) : Wall_Node(loc) {};

Corner_Node::Corner_Node(Coord loc, Wall_Node* left, Wall_Node* right, double angle) 
    : Wall_Node(loc, left, right, angle) {}

Coord Corner_Node::calc_Forces() {
	Coord sum;

	//calc Morse_SC
	sum += calc_Morse_SC();

	//calc Morse_DC
	//calc linear

	return sum;
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



