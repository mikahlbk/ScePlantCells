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

Cell::Cell(int rank, Coord corner, double height, double width, int Ti, Tissue* tiss)    {
	
	num_wall_nodes = 0;
	num_cyt_nodes = 0;	
	this->my_tissue = tiss;
	this->rank = rank;

	cell_center = Coord(corner.get_X() + (width / 2), corner.get_Y() + (height / 2));

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

	//update cell_Center
	double x_val = ((corners.at(0)->get_Location().get_X()) + (corners.at(1)->get_Location().get_X())) / 2;
	double y_val = ((corners.at(0)->get_Location().get_Y()) + (corners.at(3)->get_Location().get_Y())) / 2;

	cell_center = Coord(x_val, y_val);

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
Coord Cell::get_Cell_Center() {
	return cell_center;
}

void Cell::get_CytNodes(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
	return;
}

Wall_Node* Cell::get_WallNodes() {
	return corners.at(0);
}

void Cell::get_CornerNodes(vector<Wall_Node*>& corns) {
	corns = this->corners;
	return;
}

void Cell::update_Neighbor_Cells() {
	//clear prev vector of neigh cells
	neigh_cells.clear();
	//grab all cells from tissue
	vector<Cell*> all_Cells;
	my_tissue->get_Cells(all_Cells);

	
	// Empty variables for holding info about other cells
	Cell* curr = NULL;
	Coord curr_Cent;
	double curr_cX, curr_cY;
	double curr_minX, curr_maxX, curr_minY, curr_maxY;
	vector<Wall_Node*> curr_corners;
	
	// All necessary info about my cell location
	double my_cX = cell_center.get_X();
	double my_cY = cell_center.get_Y();
	double my_minX = corners.at(0)->get_Location().get_X();
	double my_maxX = corners.at(1)->get_Location().get_X();
	double my_minY = corners.at(0)->get_Location().get_Y();
	double my_maxY = corners.at(3)->get_Location().get_Y();

	double prelim_threshold = 5.0;
	double sec_threshold = 0.5;

	bool checkA = false;
	bool checkB = false;

	// iterate through all cells
	for (unsigned int i = 0; i < all_Cells.size(); i++) {
		curr = all_Cells.at(i);
		//reset boolean variables
		checkA = false;
		checkB = false;
		//check to make sure not pointing at yourself
		if (curr != this) {
			curr_Cent = curr->get_Cell_Center();
			// Check if cell centers are close enough together
			if ( (curr_Cent - cell_center).length() < prelim_threshold ) {
				
				curr_cX = curr_Cent.get_X();
				curr_cY = curr_Cent.get_Y();

				curr->get_CornerNodes(curr_corners);


				if (curr_cX < my_cX) {
					if (curr_cY < my_cY) {
						//check if curr_maxX and curr_maxY are in range
						// of my_minX and my_minY
						curr_maxX = curr_corners.at(2)->get_Location().get_X();
						curr_maxY = curr_corners.at(2)->get_Location().get_Y();
						
						if ( (curr_maxX > my_minX) || ( (my_minX - curr_maxX) < sec_threshold) ) {
							checkA = true;
						}

						if ( (curr_maxY > my_minY) || ( (my_minY - curr_maxY) < sec_threshold) ) {
							checkB = true;
						}

					}
					else {
						//check if curr_maxX and curr_minY are in range
						// of my_minX and my_maxY
						curr_maxX = curr_corners.at(2)->get_Location().get_X();
						curr_minY = curr_corners.at(0)->get_Location().get_Y();
						
						if ( (curr_maxX > my_minX) || ( (my_minX - curr_maxX) < sec_threshold) ) {
							checkA = true;
						}

						if ( (curr_minY < my_maxY) || ( (curr_minY - my_maxY) < sec_threshold) ) {
							checkB = true;
						}

					}
				}
				else { // curr_cX > my_cX
					if (curr_cY < my_cY) {
						//check if curr_minX and curr_maxY are in range
						// of my_maxX and my_minY
						curr_minX = curr_corners.at(0)->get_Location().get_X();
						curr_maxY = curr_corners.at(2)->get_Location().get_Y();
						
						if ( (curr_minX < my_maxX) || ( (curr_minX - my_maxX) < sec_threshold) ) {
							checkA = true;
						}

						if ( (curr_maxY > my_minY) || ( (my_minY - curr_maxY) < sec_threshold) ) {
							checkB = true;
						}

					}
					else {
						//check if curr_minX and curr_minY are in range
						// of my_maxX and my_maxY
						curr_minX = curr_corners.at(0)->get_Location().get_X();
						curr_minY = curr_corners.at(0)->get_Location().get_Y();
						
						if ( (curr_minX < my_maxX) || ( (curr_minX - my_maxX) < sec_threshold) ) {
							checkA = true;
						}

						if ( (curr_minY < my_maxY) || ( (curr_minY - my_maxY) < sec_threshold) ) {
							checkB = true;
						}

					}
				}

			}
			//else already too far away
				
			// if both checks come out true, then add curr to neigh cells
			if ( checkA && checkB) {
				neigh_cells.push_back(curr);
			}
			
		}
		//else you're pointing at yourself and shouldnt do anything

	}
	
	cout << "Cell: " << rank << " -- neighbors: " << neigh_cells.size() << endl;

	return;
}

void Cell::get_Neighbor_Cells(vector<Cell*>& cells) {
	//copy over neighboring cells
	cells = neigh_cells;

	return;
}

int Cell::get_Num_Nodes() {
	return num_wall_nodes + num_cyt_nodes;
}

bool Cell::get_Reasonable_Bounds(Wall_Node* curr, Wall_Node* & A, Wall_Node* & B) {

	bool close_enough = true;

	//Is the curr wall_node below, above, left or right of the cell

	int grid = 0;
	double threshold = 0.5;

	double min_x = corners.at(0)->get_Location().get_X();
	double min_y = corners.at(0)->get_Location().get_Y();
	double max_x = corners.at(2)->get_Location().get_X();
	double max_y = corners.at(2)->get_Location().get_Y();

	double x_val = curr->get_Location().get_X();
	double y_val = curr->get_Location().get_Y();

	if (y_val < min_y) {
		if (x_val < min_x) { 
			grid = 1;
			A = corners.at(0);
			B = A;
		}
		else if (x_val > max_x) { 
			grid = 3;
			A = corners.at(1);
			B = A;
		}
		else { 
			grid = 2;
			A = corners.at(0);
			B = corners.at(1);
		}
	}
	else if (y_val > max_y) {
		if (x_val < min_x) { 
			grid = 7; 
			A = corners.at(3);
			B = A;
		}
		else if (x_val > max_x) { 
			grid = 5; 
			A = corners.at(2);
			B = A;
		}
		else { 
			grid = 6; 
			A = corners.at(2);
			B = corners.at(3);
		}
	}
	else {
		if (x_val < min_x) { 
			grid = 8; 
			A = corners.at(3);
			B = corners.at(0);
		}
		else { 
			grid = 4; 
			A = corners.at(1);
			B = corners.at(2);
		}
	}

	
	// Check if you're even close enough

	if (grid % 2 == 0) {
		
		if ( (curr->get_Location() - A->get_Location()).length() > threshold ) {
			close_enough = false;
		}
		else {
			close_enough = true;
		}
	}
	else {
		Coord midpoint = (A->get_Location() + B->get_Location()) / 2;
		if ( (curr->get_Location() - midpoint).length() < threshold ) {
			close_enough = true;
		}
		else if ( (curr->get_Location() - A->get_Location()).length() < threshold) {
			close_enough = true;
		}
		else if ( (curr->get_Location() - B->get_Location()).length() < threshold) {
			close_enough = true;
		}
		else {
			close_enough = false;
		}
	}


	return close_enough;
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
		//cout << " calc force for this node" << endl;
		curr->calc_Forces();
		curr = curr->get_Left_Neighbor();
	
	} while (curr != corners.at(0));

	//cout << "Finished Cell::calc_NewForces" << endl;

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

	//update cell_Center
	double x_val = ((corners.at(0)->get_Location().get_X()) + (corners.at(1)->get_Location().get_X())) / 2;
	double y_val = ((corners.at(0)->get_Location().get_Y()) + (corners.at(3)->get_Location().get_Y())) / 2;

	cell_center = Coord(x_val, y_val);
	
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

void Cell::print_VTK_Scalars(ofstream& ofs) {

	Wall_Node* curr_wall = corners.at(0);
	do {
		Coord force = curr_wall->get_Force();
		ofs << force.length() << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while(curr_wall != corners.at(0));

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord force = cyt_nodes.at(i)->get_Force();
		ofs << force.length() << endl;
	}

	return;
}

void Cell::print_VTK_Vectors(ofstream& ofs) {

	Wall_Node* curr_wall = corners.at(0);
	do {
		Coord force = curr_wall->get_CytForce();
		ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while(curr_wall != corners.at(0));

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord force = cyt_nodes.at(i)->get_Force();
		ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;
	}

	return;
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





