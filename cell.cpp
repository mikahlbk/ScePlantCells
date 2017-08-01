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
#include "tissue.h"
//===================

// Cell Class Member functions

// Constructors

Cell::Cell(int rank, Coord corner, double height, 
			double width, int Ti, Tissue* tiss)    {
	num_wall_nodes = 0;
	num_cyt_nodes = 0;	
	this->my_tissue = tiss;
	this->rank = rank;

	init_cell_time = Ti;

	int num_Init_Wall_Nodes = 100;
	double perim = (height * 2) + (width * 2);
	//space between wall nodes
	double space = perim / num_Init_Wall_Nodes;
	int num_end_nodes = (width / space) - 1;
	int num_flank_nodes = (height / space) - 1;

	/*
	int num_end_nodes = 19;
	int num_flank_nodes = 34;

	double f_space = height / num_flank_nodes;
	double e_space = width / num_end_nodes;
	*/
	double curr_X;
	double curr_Y;
	Coord location;

	Wall_Node* currW;
	Wall_Node* prevW;

	// Create first corner
	prevW = new Corner_Node(corner, this);
	corners.push_back(prevW);
	num_wall_nodes++;

	//create lower end
	curr_X = corner.get_X() + space;
	curr_Y = corner.get_Y();

	for (int i = 0; i < num_end_nodes; i++) {
		location = Coord(curr_X, curr_Y);
		currW = new End_Node(location, this);

		// Set neighbor relationships
		currW->set_Right_Neighbor(prevW);
		prevW->set_Left_Neighbor(currW);
		
		//update for next iteration
		prevW = currW;
		curr_X += space;
		num_wall_nodes++;
	}

	//create second corner
	location = Coord(curr_X, curr_Y);
	currW = new Corner_Node(location, this);
	currW->set_Right_Neighbor(prevW);
	prevW->set_Left_Neighbor(currW);
	corners.push_back(currW);
	prevW = currW;
	num_wall_nodes++;

	//create right flank
	//   curr_X should be good
	curr_Y += space;
	
	for (int i = 0; i < num_flank_nodes; i++) {
		location = Coord(curr_X, curr_Y);
		currW = new Flank_Node(location, this);

		// Set neighbor relationships
		currW->set_Right_Neighbor(prevW);
		prevW->set_Left_Neighbor(currW);
		
		//update for next iteration
		prevW = currW;
		curr_Y += space;
		num_wall_nodes++;
	}

	//create third corner
	location = Coord(curr_X, curr_Y);
	currW = new Corner_Node(location, this);
	currW->set_Right_Neighbor(prevW);
	prevW->set_Left_Neighbor(currW);
	corners.push_back(currW);
	prevW = currW;
	num_wall_nodes++;

	//create upper end
	curr_X -= space;
		//curr_Y should be good
	
	for (int i = 0; i < num_end_nodes; i++) {
		location = Coord(curr_X, curr_Y);
		currW = new End_Node(location, this);

		// Set neighbor relationships
		currW->set_Right_Neighbor(prevW);
		prevW->set_Left_Neighbor(currW);
		
		//update for next iteration
		prevW = currW;
		curr_X -= space;
		num_wall_nodes++;
	}

	//create fourth corner
	location = Coord(curr_X, curr_Y);
	currW = new Corner_Node(location, this);
	currW->set_Right_Neighbor(prevW);
	prevW->set_Left_Neighbor(currW);
	corners.push_back(currW);
	prevW = currW;
	num_wall_nodes++;

	//create left flank
	//	remember that we already have first corner. 
	//	don't make it again.
	
	//   curr_X should be good
	curr_Y -= space;
	
	for (int i = 0; i < num_flank_nodes; i++) {
		location = Coord(curr_X, curr_Y);
		currW = new Flank_Node(location, this);

		// Set neighbor relationships
		currW->set_Right_Neighbor(prevW);
		prevW->set_Left_Neighbor(currW);
		
		//update for next iteration
		prevW = currW;
		curr_Y -= space;
		num_wall_nodes++;
	}

	//connect first and last nodes
	prevW->set_Left_Neighbor(corners.at(0));
	corners.at(0)->set_Right_Neighbor(prevW);

	//initialize angles
	update_Wall_Angles();

	//make cytoplasm node
	// For rectangular distribution
	double corner_offset = 0.1;
	double scal_x_offset = 0.9;
	double scal_y_offset = 0.9;
	Coord corn_off(corner_offset, corner_offset);

	int init_Cyt_Count = 20;
	/* For circular distribution
	double init_Cell_Radius = .9375;
	double random_angle;
	double random_radius;
	*/
	double x;
	double y;

	for(int i = 0; i < init_Cyt_Count; i++) {
		
		// USING POSITIONS OF CORNER NODES FOR CYT NODE ALLOCATION
		// ---distributes more evenly throughout start cell
		x = (static_cast<double>(rand()) / RAND_MAX) * scal_x_offset * (width - corner_offset);
		y = (static_cast<double>(rand()) / RAND_MAX) * scal_y_offset * (height - corner_offset); 

		Coord rand_off(x,y);
		
		location = (corner + corn_off + rand_off);

		/* USING ANGLE AND RADIUS FOR CYT NODE ALLOCATION

		random_angle = 2*pi*(static_cast<double>(rand()) / RAND_MAX);
		random_radius = init_Cell_Radius*(static_cast<double>(rand()) / RAND_MAX);
		x = (random_radius * cos(random_angle)) + corner.get_X() + (width / 2);
		//also need to add the x and y value of corner node
		//of this cell to put these values in the right place
		y = (random_radius * sin(random_angle)) + corner.get_Y() + (height / 2);
		*/

		//create cyt node
		Cyt_Node* cyt = new Cyt_Node(location,this);
		cyt_nodes.push_back(cyt);
		num_cyt_nodes++;
	}
}

// Getters and Setters

void Cell::get_CytNodes(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
	return;
}

Wall_Node* Cell::get_WallNodes() {
	return corners.at(0);
}

void Cell::get_Neighbor_Cells(vector<Cell*>& cells) {
	my_tissue->get_Cells(cells);
	return;
}

int Cell::get_Num_Nodes() {
	return num_wall_nodes + num_cyt_nodes;
}

// Calc Force

void Cell::calc_New_Forces() {
	//calc forces on cyt nodes
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->calc_Forces();
	}

	//calc forces on wall nodes
	Wall_Node* curr = corners.at(0);
	
	do {
		curr->calc_Forces();
		curr = curr->get_Left_Neighbor();
	
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
	
	ofs << "This is where data output goes" << endl;

	return;
}

void Cell::print_VTK_Points(ofstream& ofs, int& count) {

	Wall_Node* curr_wall = corners.at(0);
	do {
		Coord loc = curr_wall->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		count++;
	} while(curr_wall != corners.at(0));

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord loc = cyt_nodes.at(i)->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
		count++;
	}

	return;
}



Wall_Node* Cell::find_Largest_Length(int& side) {
	side = 0; //know that we start on bottom flank
	// We know we start at first entry of corner vector. 
	//Each time we pass another corner, increment side by 1 
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

void Cell::add_Wall_Node(const int Ti) {

	if ((Ti - init_cell_time) % ADD_WALL_TIMER != 0) {
		return;
	}

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
//	cout << "Right side coordinates" << endl;
//	cout << right_Coords << endl;
	Coord left_Coords = left_Node->get_Location();
//	cout << "Left side coordinates" << endl;
//	cout << left_Coords << endl;
	//add it halfway between these two coords

	Coord new_Coords = (right_Coords + left_Coords)*(.5);
//	cout << "New Coords: " << new_Coords << endl;
	Wall_Node* new_Node;
	
	if (side % 2 == 1 ) { //end
		new_Node = new End_Node(new_Coords, this, left_Node, right_Node);
	}
	else { //flank
		new_Node = new Flank_Node(new_Coords, this, left_Node, right_Node);
	}
	
	right_Node->set_Left_Neighbor(new_Node);
	left_Node->set_Right_Neighbor(new_Node);

	num_wall_nodes++;
	return;
}

void Cell::add_Cyt_Node(const int Ti) {

	if ((Ti - init_cell_time) % ADD_CYT_TIMER != 0) {
		return;
	}

	// Puts cyt node anywhere in cell with a min dist away from walls
	double min_x = max(corners.at(0)->get_Location().get_X(), corners.at(3)->get_Location().get_X());
	double max_x = min(corners.at(1)->get_Location().get_X(), corners.at(2)->get_Location().get_X()); 
	double min_y = max(corners.at(0)->get_Location().get_Y(), corners.at(1)->get_Location().get_Y());
	double max_y = min(corners.at(2)->get_Location().get_Y(), corners.at(3)->get_Location().get_Y());

	double corn_offset = 0.08;
	Coord init_off(corn_offset, corn_offset);
	Coord min_corner(min_x, min_y);

	double scal_x_offset = 0.75;
	double scal_y_offset = 0.9;
	double scal_x = (static_cast<double>(rand()) / RAND_MAX) * scal_x_offset * (max_x - (min_x + corn_offset));
	double scal_y = (static_cast<double>(rand()) / RAND_MAX) * scal_y_offset * (max_y - (min_y + corn_offset));
	Coord scal_C(scal_x, scal_y);

	/* Puts new cyt node in dead center of cell
	Coord len_mid = (corners.at(0)->get_Location() + corners.at(3)->get_Location()) * 0.5;
	Coord width_mid = (corners.at(0)->get_Location() + corners.at(1)->get_Location()) * 0.5;
	double new_x = width_mid.get_X();  
	double new_y = len_mid.get_Y();
	*/

	Coord new_Coord = (min_corner + init_off + scal_C);

	Cyt_Node* cyt = new Cyt_Node(new_Coord, this);
	cyt_nodes.push_back(cyt);

	num_cyt_nodes++;
	return;
}





