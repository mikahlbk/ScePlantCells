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

	Coord corn_off(0.05 * width, 0.05 * height);
	double scal_x_offset = 0.9;
	double scal_y_offset = 0.9;

	int init_Cyt_Count = 20;

	double x;
	double y;

	for(int i = 0; i < init_Cyt_Count; i++) {
		
		// USING POSITIONS OF CORNER NODES FOR CYT NODE ALLOCATION
		// ---distributes more evenly throughout start cell
		x = (static_cast<double>(rand()) / RAND_MAX) * scal_x_offset * width;
		y = (static_cast<double>(rand()) / RAND_MAX) * scal_y_offset * height; 

		Coord rand_off(x,y);
		
		location = (corner + corn_off + rand_off);

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


void Cell::find_Big_Gaps(vector<Wall_Node*>& walls, vector<int>& sides) {
	
	int side = 0;
	// We know we start at first entry of corner vector. 
	//Each time we pass another corner, increment side by 1 
	// side = 1 or 3 => end
	// side = 2 or 4 => flank
	Wall_Node* curr = corners.at(0);

	Coord left_loc;
	Coord curr_loc;
	Coord diff_vect;
	double cut_off = 0.15;
	double len;

	do {
		//if encounter a corner, increment side to know if end or flank
		if (curr->is_Corner()) {
			side++;
		}

		//finding current lengths and comparing
		left_loc = curr->get_Left_Neighbor()->get_Location();
		curr_loc = curr->get_Location();
		diff_vect = left_loc - curr_loc;
		len = diff_vect.length();
		if(len > cut_off) {
			//Push wall node and it's side ID onto vectors
			walls.push_back(curr);
			sides.push_back(side);
		}

		curr = curr->get_Left_Neighbor();

	} while (curr != corners.at(0));

	return;
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
	vector<Wall_Node*> walls;
	vector<int> sides;

	find_Big_Gaps(walls, sides);

	cout << "Find_Big_Gaps size: " << walls.size() << endl;

	Wall_Node* right = NULL;
	Wall_Node* left = NULL;
	Coord newCoord;
	Wall_Node* newNode = NULL;

	for (unsigned int i = 0; i < walls.size(); i++) {
		//Get nodes to right and left of new node
		right = walls.at(i);
		left = right->get_Left_Neighbor();
		// New node's location
		newCoord = ((right->get_Location() + left->get_Location()) * 0.5);
		
		if (sides.at(i) % 2 == 1) {
			newNode = new End_Node(newCoord, this, left, right);
		}
		else { //flank
			newNode = new Flank_Node(newCoord, this, left, right);
		}
		//increment wall node counter
		num_wall_nodes++;
		//update neighbors of left and right
		right->set_Left_Neighbor(newNode);
		left->set_Right_Neighbor(newNode);
	}

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

	Coord init_off((max_x - min_x) * 0.05, (max_y - min_y) * 0.05);
	Coord min_corner(min_x, min_y);

	double scal_x_offset = 0.9;
	double scal_y_offset = 0.9;
	double scal_x = (static_cast<double>(rand()) / RAND_MAX) * scal_x_offset * (max_x - min_x);
	double scal_y = (static_cast<double>(rand()) / RAND_MAX) * scal_y_offset * (max_y - min_y);
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

// Gets called by update_Node_Locations()
Cell* Cell::divide_length_wise() {
	// Current cell splits into two daughter cells
	//    -"this" will keep its entity as the left sister, 
	//		create a sister cell to the right of it, and return it 
	//      to the tissue.

	//ask tissue for new id num
//	int new_id = my_tissue->get();
//	Cell* sister = new Cell(new_id, my_tissue);

	// Find midpoint
	Coord mid_top = ((corners.at(3)->get_Location() + corners.at(2)->get_Location()) * 0.5);
	Coord mid_bot = ((corners.at(0)->get_Location() + corners.at(1)->get_Location()) * 0.5);
	Coord divider = (mid_top - mid_bot);
	double divide_X = divider.get_X();

	vector<Wall_Node*> corners_A;
	vector<Wall_Node*> corners_B;

	Wall_Node* curr = corners.at(0);
	Wall_Node* next = NULL;
	Wall_Node* prev = NULL;
	
	int side = 0;
	corners_A.push_back(curr);
	// distribute wall nodes between sister cells
	do { 

		next  = curr->get_Left_Neighbor();
			

	} while(next->get_Location().get_X() < divide_x);

	// curr is new corner1
	corners_A.push_back(curr);
	corners_B.push_back(next);
	corners_B.push_back(corners.at(1));
	corners_B.push_back(corners.at(2));
	
	curr = corners.at(3);
	do { 

		next  = curr->get_Right_Neighbor();
			

	} while(next->get_Location().get_X() < divide_x);
	corners_A.push_back(curr);
	corners_A.push_back(corners.at(3));
	corners_B.push_back(next);
	
	//reset corner vector of cell A and add side to nodes
	int num_nodes = 0;
	curr = corners.at(1);
	do {
		next = curr->get_Left_Neighbor();
		num_nodes++;
	} while(next->isCorner() ==  False)
	
	Coord total_dist_coord = corners_A.at(2)-corners_A.at(1);
	double total_dist  = total_dist_coord.length();
	double space = total_dist/(num_nodes-1);
	prevw = corners_A.at(1);
	curr get x x coord
	curr get y y coord
	
	for (int i = 0; i < num_nodes-1; i++) {
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

	prevw = corners_B.at(0);
	curr get x x coord
	curr get y y coord
	
	for (int i = 0; i < num_nodes-1; i++) {
		location = Coord(curr_X, curr_Y);
		currW = new Flank_Node(location, this);

		// Set neighbor relationships
		currW->set_Left_Neighbor(prevW);
		prevW->set_Right_Neighbor(currW);
		
		//update for next iteration
		prevW = currW;
		curr_Y += space;
		num_wall_nodes++;	
	}


	//add side to nodes for cell b
	//same as above
	//split up cytoplasm
	//reassign cyt vector to A
	//make sister cell by feeding in starter corner and cyt node
	//create new cell with inputs as the corner vec and cyt vec
	//

	

	return sister;
}


// Gets called by update_Node_Locations()
Cell* Cell::divide_width_wise() {
	// Current cell splits into two daughter cells
	//    -"this" will keep its entity, create a sister cell, and return it 
	//      to the tissue.
	Cell* sister = NULL;

	// Find midpoint


	return sister;
}








/*
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
*/

