//cell.cpp
//===================
// Forward Declarations

//===================
// Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "coord.h"
#include "node.h"
#include "cell.h"
//===================

// Cell Class Member functions

// Constructors

Cell::Cell(string filename) {
	
	ifstream ifs(filename.c_str());

	if(!ifs) {
		cout << filename << " is not available" << endl;
		return 1;
	}

	stringstream ss;
	string line;
	string temp;
	char comma;
	double x, y;
	Wall_Node* prev_wall = NULL;
	Wall_Node* curr_wall = NULL;
	Cyt_Node* cyt = NULL;

	while (getline(ifs,line)) {
		ss.str(line);

		getline(ss,temp,':');
		
		if (temp == "Wall_Nodes") {
			cout << "Started taking in Wall Node Locations" << endl;
		}
		else if (temp == "Cyt_Nodes") {
			//wrap back around to first wall node
			first_corner->set_Right_Neighbor(prev_wall);
			prev_wall->set_Left_Neighbor(first_corner);

			//stuff for cyt_nodes

		}
		else if (temp == "Corner") {
			ss >> x;
			ss >> comma;
			ss >> y;
			Coord temp(x,y);
			double angle = pi / 2;
			
			curr_wall = new Corner_Node(temp,angle);

			if (prev_wall == NULL) {
				first_corner = curr_wall;
			}
			else {
				prev_wall->set_Left_Neighbor(curr_wall);
				curr_wall->set_Right_Neighbor(prev_wall);
			}

			prev_wall = curr_wall;
		}
		else if (temp == "End") {
			ss >> x;
			ss >> comma;
			ss >> y;
			Coord temp(x,y);
			double angle = pi;
			
			curr_wall = new End_Node(temp,angle);

			prev_wall->set_Left_Neighbor(curr_wall);
			curr_wall->set_Right_Neighbor(prev_wall);

			prev_wall = curr_wall;
		}
		else if (temp == "Flank") {
			ss >> x;
			ss >> comma;
			ss >> y;
			Coord temp(x,y);
			double angle = pi;

			curr_wall = new Flank_Node(temp,angle);

			prev_wall->set_Left_Neighbor(curr_wall);
			curr_wall->set_Right_Neighbor(prev_wall);

			prev_wall = curr_wall;
		}
		else if (temp == "Cyt") {
			ss >> x;
			ss >> comma;
			ss >> y;
			Coord temp(x,y);
			
			cyt = new Cyt_Node(temp);
			cyt_nodes.push_back(cyt);
		}
		
		ss.clr();
	}

}

// Getters and Setters

void Cell::get_CytNodes(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
}

Wall_Node* Cell::get_WallNodes() {
	return first_corner;
}

// Calc Force

void Cell::calc_New_Forces() {
	//calc forces on cyt nodes
	for (int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->calc_Forces();
	}

	//calc forces on wall nodes
	Wall_Node* curr = first_corner;
	
	do {
		curr->calc_Forces();
		curr = curr->get_Left_Neighbor();
	
	} while(curr != first_corner);

	return;
}

// Update Node Locations

void Cell::update_Node_Positions() {
	
	//update cyt nodes
	for (int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->update_Location();
	}

	//update wall nodes
	Wall_Node* curr = first_corner;
	
	do {
		curr->update_Location();
		curr = curr->get_Left_Neighbor();
	
	} while(curr != first corner);
	
	return;
}

Wall_Node* Cell:: find_Largest_Length() {
	Wall_Node* curr = first_corner;
	Wall_Node* biggest = first_corner;
	Coord left_Neighb_loc;
	Coord curr_Loc;
	Coord diff_vect;
	double curr_len = 0;
	double len;
	//loop through all possible Cell Wall 'links' to find biggest
	for (i = first_corner;i = first_corner->get_Right_Neighbor();i = curr->get_Left_Neighbor()) {
		//finding current lengths and comparing
		curr = i;
		left_Neighb_loc = curr->get_Left_Neighbor()->get_Location();
		curr_Loc = curr->get_Location();
		diff_vect = left_Neighb_loc - curr_Loc;
		len = diff_vect.length();
		if(len > curr_len) {
			curr_len = len;
			biggest = curr;
		}
	}
	return biggest;
}

void Cell::add_Cell_Wall_Node() {
	//we will add the node based on what the function find_Largest_Length() returns
	//find_Largest_Length() returns a pointer to the node where the largest
	//Cell Wall link is to the left of that node
	//first we apply find_Largest_Length() to the cell to get a pointer to the right
	//of where the new node is added
	Wall_Node*  right_Node = find_Largest_Length();
	//to the left of this node will be the node to the left of the new node
	Wall_Node*  left_Node = right_Node->get_Left_neighbor();
	//now we find the coords of each of these nodes to use to find the new coords
	Coord right_Coords = right_Node->get_Location();  
	Coord left_Coords = left_Node->get_Location();
	//add it halfway between these two coords
	Coord new_Coords = (right_Coords + left_Coords)*(1/2);
	//make the new wall node
	Node new_Node(new_Coords);
	new_Node.set_Left_Neighbor(left_Node);
	new_Node.set_Right_Neighbor(right_Node);
	right_Node->set_Left_Neighbor(new_Node*);
	left_Node->set_Right_Neighbor(new_Node*);
}

void Cell::add_Cyt_Node() {
	double len;
	double width;
	Coord len_vect;
	Coord width_vect;
	Coord new_Coords;
	Wall_Node* second_corner;
	Wall_Node* fourth_corner;
	Wall_Node* curr = first_corner;
	double i = 1;
	do {
		curr = curr->get_Left_Neighbor();
		//if(curr is a corner node) {
			i++;
		}
	} while(i<2);
	second_corner = curr;
	curr = first_corner;
	double i = 1;
	do {
		curr = curr->get_Left_Neighbor();
		//if(curr is a corner node) {
			i++;
		}
	} while(i<3);
	fourth_corner = curr;
	len_vect = first_corner->get_Location()-fourth_corner->get_Location();
	width_vect = first_corner->get_Location() - second_corner->get_Location();
	//new_Coords_x = random number between zero and one multiplied by the length	
	//new_Coords_y = same	
	//new_Coords = (x,y)
	Cyt_Node new_Cyt_Node(new_Coords);
}
	






