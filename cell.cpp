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

	//Initialize node counters to 0
	num_wall_nodes = 0;
	num_cyt_nodes = 0;
	vector<Coord> init_walls;
	vector<Coord> init_cyts;
	
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
			
			curr_wall = new Corner_Node(temp);

			if (prev_wall != NULL) {
				prev_wall->set_Left_Neighbor(curr_wall);
				curr_wall->set_Right_Neighbor(prev_wall);
			}
			corners.push_back(curr_wall);
			prev_wall = curr_wall;
			
			init_walls.push_back(temp);
			num_wall_nodes++;
		}
		else if (temp == "End") {
			ss >> x;
			ss >> comma;
			ss >> y;
			Coord temp(x,y);
			
			curr_wall = new End_Node(temp);

			prev_wall->set_Left_Neighbor(curr_wall);
			curr_wall->set_Right_Neighbor(prev_wall);
			prev_wall = curr_wall;

			init_walls.push_back(temp);
			num_wall_nodes++;
		}
		else if (temp == "Flank") {
			ss >> x;
			ss >> comma;
			ss >> y;
			Coord temp(x,y);

			curr_wall = new Flank_Node(temp);

			prev_wall->set_Left_Neighbor(curr_wall);
			curr_wall->set_Right_Neighbor(prev_wall);
			prev_wall = curr_wall;

			init_walls.push_back(temp);
			num_wall_nodes++;
		}
		else if (temp == "Cyt") {
			ss >> x;
			ss >> comma;
			ss >> y;
			Coord temp(x,y);
			
			cyt = new Cyt_Node(temp);
			cyt_nodes.push_back(cyt);

			init_cyts.push_back(temp);
			num_cyt_nodes++;
		}
		
		ss.clr();
	}

	ifs.close();

	// keep track of initial locations
	wall_node_locs.push_back(init_walls);
	cyt_node_locs.push_back(init_cyts);
	wall_nodes_per_frame.push_back(num_wall_nodes);
	cyt_nodes_per_frame.push_back(num_cyt_nodes);
	// update angles of wall nodes
	update_Wall_Angles();
}

// Getters and Setters

void Cell::get_CytNodes(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
	return;
}

Wall_Node* Cell::get_WallNodes() {
	return first_corner;
}

void get_Neigh_Cells(vector<Cell*>& cells) {
	cells = neigh_cells;
	return;
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

	//
	
	return;
}

void Cell::update_Wall_Angles() {

	Wall_Node* curr = corners.at(0);
	
	do {
		curr->update_Angle();
		curr = curr->get_Left_Neighbor();
	} while (curr != corners.at(0));	

	return;
}

Wall_Node* Cell::find_Largest_Length(int& side) {
	side = 0; //know that we start on bottom flank
	/* We know we start at first entry of corner vector. 
		Each time we pass another corner, increment side by 1  */
	// side = 1 or 3 => end
	// side = 2 or 4 => flank
	Wall_Node* curr = first_corner;
	Wall_Node* biggest = NULL;
	Coord left_Neighb_loc;
	Coord curr_Loc;
	Coord diff_vect;
	double max_len = 0;
	double len;
	//loop through all possible Cell Wall 'links' to find biggest
	do {
		//if encounter a corner, increment side to know if end or flank
		if (curr->is_Corner()) {
			side++;
		}

		//finding current lengths and comparing
		left_Neighb_loc = curr->get_Left_Neighbor()->get_Location();
		curr_Loc = curr->get_Location();
		diff_vect = left_Neighb_loc - curr_Loc;
		len = diff_vect.length();
		if(len > max_len) {
			max_len = len;
			biggest = curr;
		}
		curr = curr->get_Left_Neighbor();

	} while (curr != corners.at(0));

	return biggest;
}

void Cell::add_Cell_Wall_Node() {
	//we will add the node based on what the function find_Largest_Length() returns
	//find_Largest_Length() returns a pointer to the node where the largest
	//Cell Wall link is to the left of that node
	//first we apply find_Largest_Length() to the cell to get a pointer to the right
	//of where the new node is added
	Wall_Node* right_Node = find_Largest_Length();
	//to the left of this node will be the node to the left of the new node
	Wall_Node* left_Node = right_Node->get_Left_neighbor();
	//now we find the coords of each of these nodes to use to find the new coords
	Coord right_Coords = right_Node->get_Location();  
	Coord left_Coords = left_Node->get_Location();
	//add it halfway between these two coords

	Coord new_Coords = (right_Coords + left_Coords)*(1/2);
	
	if ( ) { //end
		Wall_Node* new_Node = new End_Node(new_Coords, left_Node, right_Node);
	}
	else { //flank
		Wall_Node* new_Node = new Flank_Node(new_Coords, left_Node, right_Node);
	}
	
	right_Node->set_Left_Neighbor(new_Node);
	left_Node->set_Right_Neighbor(new_Node);
}

void Cell::add_Cyt_Node() {

	Coord len_vect = corners.at(0)->get_Location() - corners.at(3)->get_Location();
	Coord width_vect = corners.at(0)->get_Location() - corners.at(1)->get_Location();
	double new_x = (rand() / RAND_MAX) * width_vect.length();  
	double new_y = (rand() / RAND_MAX) * len_vect.length();

	Coord new_Coord(new_x, new_y);
	new_Coord += corners.at(0);

	Cyt_Node* cyt = new Cyt_Node(new_Coord);
	cyt_nodes.push_back(cyt);

	return;
}
	






