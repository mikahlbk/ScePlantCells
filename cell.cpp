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


