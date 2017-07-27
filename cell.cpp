//cell.cpp
//===================
// Forward Declarations

//===================
// Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>

#include "phys.h"
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
			corners.at(0)->set_Right_Neighbor(prev_wall);
			prev_wall->set_Left_Neighbor(corners.at(0));

			//stuff for cyt_nodes

		}
		else if (temp == "Corner") {
			ss >> x;
			cout << "X: " << x << endl;
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

			num_wall_nodes++;
		}
		else if (temp == "Cyt") {
			ss >> x;
			ss >> comma;
			ss >> y;
			Coord temp(x,y);
			
			cyt = new Cyt_Node(temp);
			cyt_nodes.push_back(cyt);

			num_cyt_nodes++;
		}
		
		ss.clear();
	}

	ifs.close();

	// update angles of wall nodes
	update_Wall_Angles();
}

// Getters and Setters

void Cell::get_CytNodes(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
	return;
}

Wall_Node* Cell::get_WallNodes() {
	return corners.at(0);
}

void Cell::get_Neigh_Cells(vector<Cell*>& cells) {
	cells = neigh_cells;
	return;
}

// Calc Force

void Cell::calc_New_Forces() {
	//calc forces on cyt nodes
	
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->calc_Forces(this);
	}
	
	//calc forces on wall nodes
	Wall_Node* curr = corners.at(0);
	int i = 0;
	
	do {
		cout << "CellRank: " << i << endl;
		curr->calc_Forces(this);
		curr = curr->get_Left_Neighbor();
		i++;
	
	} while(curr != corners.at(0));

	return;
}

// Update Node Locations

void Cell::update_Node_Locations() {
	
	//update cyt nodes

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->update_Location();
	}

	//update wall nodes
	Wall_Node* curr = corners.at(0);
	
	do {
		curr->update_Location();
		curr = curr->get_Left_Neighbor();
	
	} while(curr != corners.at(0));

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

void Cell::print_Data_Output(ofstream& ofs) {
	
	Wall_Node* curr = corners.at(0);
	int i = 0;

	do {
		Coord loc = curr->get_Location();
		Coord f = curr->get_New_Forces();

		ofs << "Node " << i << ':' << endl;
		ofs << "	Loc: " << loc << endl;
		ofs << "	Force: " << f << endl;

		curr = curr->get_Left_Neighbor();
		i++;

	} while(curr != corners.at(0));

	return;
}

void Cell::print_VTK_File(ofstream& ofs) {

	ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Point representing Sub_cellular elem model" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;
	ofs << "POINTS " << num_wall_nodes + num_cyt_nodes << " float" << endl;

	Wall_Node* curr_wall = corners.at(0);
	do {
		Coord loc = curr_wall->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
	} while(curr_wall != corners.at(0));

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord loc = cyt_nodes.at(i)->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
	}

	ofs << endl;

	ofs << "CELLS " << 1 << ' ' << num_wall_nodes + num_cyt_nodes + 1  << endl;

	ofs << num_wall_nodes + num_cyt_nodes;

	for (int i = 0; i < num_wall_nodes + num_cyt_nodes; i++) {
		ofs << ' ' << i;
	}

	ofs << endl << endl;

	ofs << "CELL_TYPES " << 1 << endl;
	ofs << 2 << endl;
	
	return;
}



Wall_Node* Cell::find_Largest_Length(int& side) {
	side = 0; //know that we start on bottom flank
	/* We know we start at first entry of corner vector. 
		Each time we pass another corner, increment side by 1  */
	// side = 1 or 3 => end
	// side = 2 or 4 => flank
	Wall_Node* curr = corners.at(0);
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
	int side = 0;
	Wall_Node* right_Node = find_Largest_Length(side);
	//to the left of this node will be the node to the left of the new node
	Wall_Node* left_Node = right_Node->get_Left_Neighbor();
	//now we find the coords of each of these nodes to use to find the new coords
	Coord right_Coords = right_Node->get_Location();  
	Coord left_Coords = left_Node->get_Location();
	//add it halfway between these two coords

	Coord new_Coords = (right_Coords + left_Coords)*(1/2);
	Wall_Node* new_Node;
	
	if (side % 2 == 1 ) { //end
		new_Node = new End_Node(new_Coords, left_Node, right_Node);
	}
	else { //flank
		new_Node = new Flank_Node(new_Coords, left_Node, right_Node);
	}
	
	right_Node->set_Left_Neighbor(new_Node);
	left_Node->set_Right_Neighbor(new_Node);

	num_wall_nodes++;
	return;
}

void Cell::add_Cyt_Node() {

	Coord len_vect = corners.at(0)->get_Location() - corners.at(3)->get_Location();
	Coord width_vect = corners.at(0)->get_Location() - corners.at(1)->get_Location();
	double new_x = (rand() / RAND_MAX) * width_vect.length();  
	double new_y = (rand() / RAND_MAX) * len_vect.length();

	Coord new_Coord(new_x, new_y);
	new_Coord += corners.at(0)->get_Location();

	Cyt_Node* cyt = new Cyt_Node(new_Coord);
	cyt_nodes.push_back(cyt);

	num_cyt_nodes++;
	return;
}
	






