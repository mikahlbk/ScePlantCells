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
#include "side.h"
#include "tissue.h"
//===================

// Cell Class Member functions

// Constructors
Cell::Cell(int rank, Tissue* tissue) {
	this->rank = rank;
	my_tissue = tissue;
	num_cyt_nodes = 0;
	life_length = 0;
}

Cell::Cell(int rank, Coord corner, double height, double width, int Ti, Tissue* tiss)    {
	
	num_cyt_nodes = 0;	
	this->my_tissue = tiss;
	this->rank = rank;

	Coord cell_center = Coord(corner.get_X() + (width / 2), corner.get_Y() + (height / 2));

	life_length = 0;

	//rough estimates for cell sizing
	int num_Init_Wall_Nodes = 100;
	double perim = (height * 2) + (width * 2);
	//space between wall nodes
	double space = perim / num_Init_Wall_Nodes;
	int num_end_nodes = (width / space);
	int num_flank_nodes = (height / space);
	
	// Side A
	Coord locA = ( corner + Coord(0.04,0.0) );
	Coord locZ = ( corner + Coord((width - 0.08),0.0) ); 
	Side* s = new Side(locA, locZ, this, num_end_nodes);
//	cout << "made side" << endl;
	s->set_Phys_Parameters(kBendLow, kLinearHigh);
	if ((this->get_Rank() == 3) || (this->get_Rank() == 4)) {
		s->set_Phys_Parameters(kBendHigh, kLinearLow);
    }
//	cout << "set params" << endl;
	sides.push_back(s);
//	cout << "I made side A" << endl;
	//Side B
	locA = locZ + Coord(0.04,0.04);
	locZ = locA + Coord(0.0, (height - 0.08));
	s = new Side(locA, locZ, this, num_flank_nodes);
	s->set_Phys_Parameters(kBendHigh, kLinearLow);
	if ((this->get_Rank() == 3) || (this->get_Rank() == 4)) {
		s->set_Phys_Parameters(kBendLow, kLinearHigh);
    }
	sides.push_back(s);
//	cout << "I made side B" << endl;
	//side C
	locA = locZ + Coord(-0.04, 0.04);
	locZ = locA + Coord(-(width - 0.08), 0.0);
	s = new Side(locA, locZ, this, num_end_nodes);
	s->set_Phys_Parameters(kBendLow, kLinearHigh);
	if ((this->get_Rank() == 3) || (this->get_Rank()==4)) {
		s->set_Phys_Parameters(kBendHigh, kLinearLow);
    }
	sides.push_back(s);
//	cout << "I made side C" << endl;
	//side D
	locA = locZ + Coord(-0.04, -0.04);
	locZ = corner + Coord(0.0, 0.04);
	s = new Side(locA, locZ, this, num_flank_nodes);
	s->set_Phys_Parameters(kBendHigh, kLinearLow);
	if ((this->get_Rank() == 3) || (this->get_Rank() ==4)) {
		s->set_Phys_Parameters(kBendLow, kLinearHigh);
    }
	sides.push_back(s);
//	cout << "I made side D" << endl;
	//Connect the four disjoint sides
	sides.at(0)->connect_Ends(sides.at(1));
	sides.at(1)->connect_Ends(sides.at(2));
	sides.at(2)->connect_Ends(sides.at(3));
	sides.at(3)->connect_Ends(sides.at(0));
	
	//update wall angles
	update_Wall_Angles();
	
//	cout << "I connected the sides" << endl;
	//Insert cytoplasm nodes
	num_cyt_nodes = 20;
	
	Coord corn_off(0.1 * width, 0.1 * height);
	double scal_x_offset = 0.8;
	double scal_y_offset = 0.8;
	Coord location;
	Cyt_Node* cyt;
	double x;
	double y;

	for (int i = 0; i < num_cyt_nodes; i++) {
		// USING POSITIONS OF CORNER NODES FOR CYT NODE ALLOCATION
		// ---distributes more evenly throughout start cell
		x = (static_cast<double>(rand()) / RAND_MAX) * scal_x_offset * width;
		y = (static_cast<double>(rand()) / RAND_MAX) * scal_y_offset * height; 
		Coord rand_off(x,y);
		location = (corner + corn_off + rand_off);

		cyt = new Cyt_Node(location,this);
		cyt_nodes.push_back(cyt);
	}
	
}

// Destructor
Cell::~Cell() {

	// Delete Cyt Nodes
	Cyt_Node* cyt = NULL;
	while ( !cyt_nodes.empty()) {
		cyt = cyt_nodes.at(cyt_nodes.size() - 1);
		delete cyt;
		cyt_nodes.pop_back();
	}
	// Delete Wall Nodes
	Wall_Node* curr = sides.at(0)->get_Wall_Nodes();
	Wall_Node* next = NULL;
	Wall_Node* last = curr->get_Right_Neighbor();
	last->set_Left_Neighbor(NULL);

	while (curr != NULL) {
		next = curr->get_Left_Neighbor();
		delete curr;
		curr = next;
	}

	my_tissue = NULL;
}

//=============================================================
//========================================
// Getters and Setters
//========================================
//=============================================================
Coord Cell::get_Cell_Center() {
	return cell_center;
}
int Cell::get_Rank(){
	return rank;
}

void Cell::set_Rank(const int id) {
	rank = id;
	return;
}

void Cell::get_Cyt_Nodes(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
	return;
}

Wall_Node* Cell::get_Wall_Nodes() {
	return sides.at(0)->get_End_A();
}

void Cell::get_Sides(vector<Side*>& sides) {
	sides = this->sides;
	return;
}

void Cell::set_Sides(vector<Side*>& sides) {
	this->sides = sides;
	return;
}

void Cell::get_Neighbor_Cells(vector<Cell*>& cells) {
	cells = neigh_cells;
	return;
}

int Cell::get_Node_Count() {
	int wall_cnt = 0;
	for (unsigned int i = 0; i < sides.size(); i++) {
		wall_cnt += sides.at(i)->get_Wall_Count();
	}

	return (wall_cnt + num_cyt_nodes);
}

//=============================================================
//=========================================
// Keep Track of neighbor cells
//=========================================
//=============================================================

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
	vector<Side*> curr_sides;
	
	// All necessary info about my cell location
	double my_cX = cell_center.get_X();
	double my_cY = cell_center.get_Y();
	double my_minX = sides.at(0)->get_End_A()->get_Location().get_X();
	double my_maxX = sides.at(0)->get_End_Z()->get_Location().get_X();
	double my_minY = sides.at(3)->get_End_Z()->get_Location().get_Y();
	double my_maxY = sides.at(3)->get_End_A()->get_Location().get_Y();

	double prelim_threshold = 5.0;
	double sec_threshold = 2;

	bool checkA = false;
	bool checkB = false;
	
//	cout<<"number cells in system " << all_Cells.size() << endl;
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

				curr->get_Sides(curr_sides);
				curr_maxX = curr_sides.at(0)->get_End_Z()->get_Location().get_X();
				curr_maxY = curr_sides.at(3)->get_End_A()->get_Location().get_Y();
				curr_minX = curr_sides.at(0)->get_End_A()->get_Location().get_X();
				curr_minY = curr_sides.at(3)->get_End_Z()->get_Location().get_Y();
			

				if (curr_cX < my_cX) {
					if (curr_cY < my_cY) {
						//check if curr_maxX and curr_maxY are in range
						// of my_minX and my_minY
					
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
					
						if ( (curr_maxX > my_minX) || ( (my_minX - curr_maxX) < sec_threshold) ) {
							checkA = true;
						}

						if ( (curr_minY < my_maxY) || ( (curr_minY - my_maxY) < sec_threshold) ) {
							checkB = true;
						}

					}
					//cout << "right A: " << checkA << "right B: " << checkB << endl;
				}
				else { // curr_cX > my_cX
					if (curr_cY < my_cY) {
						//check if curr_minX and curr_maxY are in range
						// of my_maxX and my_minY
					
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
					
						if ( (curr_minX < my_maxX) || ( (curr_minX - my_maxX) < sec_threshold) ) {
							checkA = true;
						}

						if ( (curr_minY < my_maxY) || ( (curr_minY - my_maxY) < sec_threshold) ) {
							checkB = true;
						}

					}
				}
				//cout << "left A: " << checkA << "left B: " << checkB << endl;

			}
			//else already too far away
				
			// if both checks come out true, then add curr to neigh cells
			if ( checkA && checkB) {
				neigh_cells.push_back(curr);
				//cout << rank << "has neighbor" << curr->get_Rank() << endl;
			}
			
		}
		//else you're pointing at yourself and shouldnt do anything

	}
	Side* side = NULL;
	for(int i = 0; i<sides.size();i++) {
		side = sides.at(i);
		side->update_Neighbor_Sides(neigh_cells);
	}
		
	cout << "Cell: " << rank << " -- neighbors: " << neigh_cells.size() << endl;

	return;
}

void Cell::update_adhesion_springs() {
	vector<Cell*> neighbor_Cells;
	this->get_Neighbor_Cells(neighbor_Cells);
	vector<Side*> neighbor_Sides;
	Side* curr_side = NULL;
	Wall_Node* curr_Node = NULL;
	Wall_Node* next_Node = NULL;
	Wall_Node* curr_Closest = NULL;
	double curr_len = 0;
	for(int i = 0; i<sides.size();i++) {
		curr_side = sides.at(i);
		curr_side->get_Neighbor_Sides(neighbor_Sides);
		curr_Node = curr_side->get_End_A();
		do {
			next_Node = curr_Node->get_Left_Neighbor();
			curr_Closest = curr_Node->find_Closest_Node(neighbor_Sides);
			curr_Node->make_Connection(curr_Closest);
			curr_Node = next_Node;
		} while(next_Node->get_My_Side() != curr_side);
	}
}
	
/*
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
*/

//===============================================================
//============================
//  Forces and Positioning
//============================
//===============================================================
void Cell::calc_New_Forces() {
	//calc forces on cyt nodes
//	cout << "entered calc forces" << endl;
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
//		cout << "entered loop" << endl;
		cyt_nodes.at(i)->calc_Forces();
	}
//	cout << "cyt nodes forces calced" << endl;
	//calc forces on wall nodes
	Wall_Node* curr = sides.at(0)->get_Wall_Nodes();
	Wall_Node* orig = curr;
	
	do {
		curr->calc_Forces();
		curr = curr->get_Left_Neighbor();
	
	} while (curr != orig);

//	cout << "wall nodes forces calced" << endl;
	return;
}


void Cell::update_Node_Locations() {
	
	//update cyt nodes
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->update_Location();
	}

	//update wall nodes
	Wall_Node* curr = sides.at(0)->get_Wall_Nodes();
	Wall_Node* orig = curr;

	do {
		curr->update_Location();
		curr = curr->get_Left_Neighbor();
	
	} while(curr != orig);

	//update cell_Center
	update_Cell_Center();
	//update wall_angles
	update_Wall_Angles();
	
	return;
}

void Cell::update_Wall_Angles() {
	
	Wall_Node* curr = sides.at(0)->get_Wall_Nodes();
	Wall_Node* orig = curr;
	do {
		curr->update_Angle();
		curr = curr->get_Left_Neighbor();
	} while (curr != orig);	

	return;
}

void Cell::update_Cell_Center() {
/*
	double min_x = (sides.at(0)->get_Wall_Nodes()->get_Location()).get_X();
	double max_x = (sides.at(2)->get_Wall_Nodes()->get_Location()).get_X();
	double min_y = (sides.at(1)->get_Wall_Nodes()->get_Location()).get_Y();
	double max_y = (sides.at(3)->get_Wall_Nodes()->get_Location()).get_Y();

	double x = (min_x + max_x) / 2;
	double y = (min_y + max_y) / 2;

	Coord center(x,y);

	cell_center = center;
*/

	Coord low = sides.at(0)->get_Wall_Nodes()->get_Location();
	Coord high = sides.at(2)->get_Wall_Nodes()->get_Location();

	Coord mid = ((low + high) * 0.5);
	cell_center = mid;

	return;
}

void Cell::update_Life_Length() {
	life_length++;

	//check if cell can add a cyt or wall node

	if (life_length % ADD_WALL_TIMER == ADD_WALL_TIMER - 1) {
		add_Wall_Node();
	}

	if (life_length % ADD_CYT_TIMER == ADD_CYT_TIMER - 1) {
		add_Cyt_Node();
	}

	return;
}

//===========================================================
//==================================
// Output Functions
//==================================
//===========================================================

void Cell::print_Data_Output(ofstream& ofs) {
	
	ofs << "This is where data output goes" << endl;

	return;
}

void Cell::print_VTK_Points(ofstream& ofs, int& count) {

	Wall_Node* curr_wall = sides.at(0)->get_Wall_Nodes();
	do {
		Coord loc = curr_wall->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		count++;
	} while (curr_wall != sides.at(0)->get_Wall_Nodes());

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord loc = cyt_nodes.at(i)->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
		count++;
	}

	return;
}

void Cell::print_VTK_Scalars(ofstream& ofs) {

	Wall_Node* curr_wall = sides.at(0)->get_Wall_Nodes();
	do {
		Coord force = curr_wall->get_CytForce();
		ofs << force.length() << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != sides.at(0)->get_Wall_Nodes());

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord force = cyt_nodes.at(i)->get_Force();
		ofs << force.length() << endl;
	}

	return;
}

void Cell::print_VTK_Vectors(ofstream& ofs) {

	Wall_Node* curr_wall = sides.at(0)->get_Wall_Nodes();
	do {
		Coord force = curr_wall->get_CytForce();
		ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while(curr_wall != sides.at(0)->get_Wall_Nodes());

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord force = cyt_nodes.at(i)->get_Force();
		ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;
	}

	return;
}

//=====================================================================
//==========================================
// Growth of Cell
//==========================================
//=====================================================================


Wall_Node* Cell::find_Largest_Length() {
	
	Wall_Node* curr = sides.at(0)->get_Wall_Nodes();
	Wall_Node* biggest = NULL;
	Wall_Node* orig = curr;
	Coord left_Neighb_loc;
	Coord curr_Loc;
	Coord diff_vect;
	double max_len = 0;
	double len;
	int big_gaps = 0;
	//loop through all possible Cell Wall 'links' to find biggest
	do {
		//finding current lengths and comparing
		left_Neighb_loc = curr->get_Left_Neighbor()->get_Location();
		curr_Loc = curr->get_Location();
		diff_vect = left_Neighb_loc - curr_Loc;
		len = diff_vect.length();
		if (len > MEMBR_THRESH_LENGTH) {
			big_gaps++;
			if(len > max_len) {
				max_len = len;
				biggest = curr;
			}
		}
		curr = curr->get_Left_Neighbor();

	} while (curr != orig);

	cout << "Cell " << rank << " -- big gaps: " << big_gaps << endl;

	return biggest;
}


void Cell::add_Wall_Node() {
	cout << "add wall node check" << endl;
	//find node to the right of largest spring
	Wall_Node* right = find_Largest_Length();
  cout << "LL found" << endl;
	//find which side this node belongs to
	if(right != NULL) {
		Side* s = right->get_My_Side();
		//tell that side to add a new node to the left of curr pointer
		s->add_Wall_Node(right);
		cout << "wall node added" << endl;
	}
	return;
}

void Cell::add_Cyt_Node() {

	Cyt_Node* cyt = new Cyt_Node(cell_center, this);
	cyt_nodes.push_back(cyt);

	num_cyt_nodes++;
	return;
}

void Cell::add_Cyt_Node_Div() {
	Side* s0 = sides.at(0);
	Side* s1 = sides.at(1);
	double x_length = s0->get_End_Z()->get_Location().get_X()-s0->get_End_A()->get_Location().get_X();
	double y_length = s1->get_End_Z()->get_Location().get_Y()-s1->get_End_A()->get_Location().get_Y();
	Coord init_off(x_length*0.05, y_length*0.05);
	Coord min_corner(s0->get_End_A()->get_Location().get_X(),s0->get_End_A()->get_Location().get_Y());

	double x_coord_offset = 0.9;
	double y_coord_offset = 0.9;
	double x_coord_scaled = (static_cast<double>(rand())/RAND_MAX)*x_coord_offset*x_length;
	double y_coord_scaled = (static_cast<double>(rand())/RAND_MAX)*y_coord_offset*y_length;
	Coord scal(x_coord_scaled, y_coord_scaled); 
	
	Coord rand_Coord = (min_corner+init_off + scal);
		
	Cyt_Node* cyt = new Cyt_Node(rand_Coord,this);
	cyt_nodes.push_back(cyt);

	num_cyt_nodes++;
	return;
}





