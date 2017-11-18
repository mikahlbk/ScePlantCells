//cell.cpp
//===================
// Forward Declarations

//===================
// Include Dependencies
#include <iostream>
#include <sstream>
#include <fstream>
#include <vector>
#include <cstdlib>
#include <ctime>

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
//===================

// Cell Class Member functions

// Constructors
Cell::Cell(Tissue* tissue) {
	my_tissue = tissue;
	//rank assigned in tissue class
	//num_cyt_nodes assigned
//	most_up = NULL;
//	most_down = NULL;
//	most_left = NULL;
//	most_right = NULL;
	//layer inherited
	//growth rate inherited
	//cell center computed
	//num_wall nodes computer
	//left_Corner assigned
	life_length = 0;
}


Cell::Cell(int rank, Coord center, double radius, Tissue* tiss, int layer)    {

	this->rank = rank;
	this->my_tissue = tiss;
	num_cyt_nodes = 0;
	this->layer = layer;
//	int init_radius = radius;
	this->cell_center = center;
	double rate = (-0.25*cell_center.length() + 11.7)*2000;
	this->set_growth_rate(rate);
	num_wall_nodes = 0;
	
	life_length = 0;
//	most_up = NULL;
//	most_down = NULL;
//	most_left = NULL;
//	most_right = NULL;
	//rough estimates for cell sizing
	int num_Init_Wall_Nodes = Init_Wall_Nodes;
	double angle_increment = (2*pi)/num_Init_Wall_Nodes;
	
	//make all wall nodes
	double curr_X;
	double curr_Y;
	Coord location = this->cell_center;;
	Wall_Node* currW;
	Wall_Node* prevW;
	Wall_Node* orig;
	double curr_theta = 0;
	curr_X = cell_center.get_X() + radius*cos(curr_theta);
	curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
	location = Coord(curr_X,curr_Y);
	//make the first note
	prevW = new Wall_Node(location,this);
	num_wall_nodes++;
	orig = prevW;
	//this will be the "starter" node
	this->left_Corner = orig;

	//make successive nodes
	for(int i = 0; i<num_Init_Wall_Nodes-1; i++) {
		curr_theta = curr_theta + angle_increment;
		curr_X = cell_center.get_X() + radius*cos(curr_theta);
		curr_Y = cell_center.get_Y() + radius*sin(curr_theta);
		location = Coord(curr_X,curr_Y);
		currW = new Wall_Node(location, this);
		num_wall_nodes++;
		prevW->set_Left_Neighbor(currW);
		currW->set_Right_Neighbor(prevW);
		prevW = currW;
	}
	
	//connect last node to starter node
	currW->set_Left_Neighbor(orig);
	orig->set_Right_Neighbor(currW);

	//insert cytoplasm nodes
	int	num_init_cyt_nodes = Init_Num_Cyt_Nodes;
	double scal_x_offset = 0.8;
	//Coord location;
	Cyt_Node* cyt;
	double x;
	double y;
	for (int i = 0; i < num_init_cyt_nodes; i++) {
		// USING POSITIONS OF CELL CENTER FOR CYT NODE ALLOCATION
		// ---distributes more evenly throughout start cell
		double rand_radius = (static_cast<double>(rand()) / RAND_MAX)*scal_x_offset*radius;
		double rand_angle = (static_cast<double>(rand()) / RAND_MAX)*2*pi;
		x = cell_center.get_X()+ rand_radius*cos(rand_angle);
		y = cell_center.get_Y()+ rand_radius*sin(rand_angle);
		location = Coord(x,y);

		cyt = new Cyt_Node(location,this);
		cyt_nodes.push_back(cyt);
		num_cyt_nodes++;
	}
	
	//update equilibrium angle
	update_Wall_Equi_Angles();
	//update wall angles
	update_Wall_Angles();
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
	Wall_Node* curr = left_Corner;
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
void Cell::get_Cyt_Nodes(vector<Cyt_Node*>& cyts) {
	cyts = cyt_nodes;
	return;
}

void Cell::set_Rank(const int id) {
	this->rank = id;
	return;
}

void Cell::set_Layer(int layer) {
	this->layer = layer;
	return;
}
void Cell::set_growth_rate(double growth_rate) {
	this->growth_rate = growth_rate;
	return;
}
void Cell::get_Neighbor_Cells(vector<Cell*>& cells) {
	cells = neigh_cells;
	return;
}

void Cell::get_Strain(vector<double>& strain) {
	strain = this->strain_vec;
	return;
}

void Cell::get_Stress(vector<double>& stress) {
	stress = this->stress_vec;
	return;
}

int Cell::get_Node_Count() {
	int node_count = 0;

	node_count = num_wall_nodes + num_cyt_nodes;

	return node_count;
}

void Cell::set_Left_Corner(Wall_Node*& new_left_corner) {
	this->left_Corner = new_left_corner;
	return;
}

void Cell::set_Wall_Count(int& number_nodes) {
	this->num_wall_nodes = number_nodes;
	return;
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
	Coord distance;
	
	double prelim_threshold = 10;
	//double sec_threshold = 1;

	// iterate through all cells
	for (unsigned int i = 0; i < all_Cells.size(); i++) {
		curr = all_Cells.at(i);
		if (curr != this) {
			curr_Cent = curr->get_Cell_Center();
			// Check if cell centers are close enough together
			distance = this->cell_center - curr_Cent;
			//cout << "Distance = " << distance << endl;
			if ( distance.length() < prelim_threshold ) {
				neigh_cells.push_back(curr);
				//cout << rank << "has neighbor" << curr->get_Rank() << endl;
			}
			
		}
		//else you're pointing at yourself and shouldnt do anything
	}	
	
//	cout << "Cell: " << rank << " -- neighbors: " << neigh_cells.size() << endl;

	return;
}

void Cell::update_adhesion_springs(int Ti) {
	int time = Ti;
	vector<Cell*>neighbors;
	this->get_Neighbor_Cells(neighbors);
	Wall_Node* curr_Node = NULL;
	Wall_Node* next_Node = NULL;
	Wall_Node* curr_Closest = NULL;
	double curr_len = 0;
	for(int i = 0; i<num_wall_nodes;i++) {
		curr_Node = left_Corner;
		do {
			next_Node = curr_Node->get_Left_Neighbor();
//			if(time < Relaxation_Time) {
//				curr_Closest = curr_Node->find_Closest_Node_Beg(neighbors);
//			}
//			else {
			curr_Closest = curr_Node->find_Closest_Node(neighbors);
//			}
			//if(curr_Closest == NULL) {
				//cout << "Did not find a curr closest" << endl;
			//}
			//else {
				//cout << curr_Closest << endl;
			//}
			curr_Node->make_Connection(curr_Closest);
			curr_Node = next_Node;
		} while(next_Node != left_Corner);
	}
}

//===============================================================
//============================
//  Forces and Positioning
//============================
//===============================================================
void Cell::calc_New_Forces() {
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->calc_Forces();
	}

	//calc forces on wall nodes
	Wall_Node* curr = left_Corner; 
	Wall_Node* orig = curr;
	int counter = 0;
	do {
		counter++;
//		cout << "Wall node number: " << counter << endl;
		curr->calc_Forces();
		curr = curr->get_Left_Neighbor();
	
	} while (curr != orig);

	return;
}

void Cell::update_Node_Locations() {
	//update cyt nodes
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		cyt_nodes.at(i)->update_Location();
	}

	//update wall nodes
	Wall_Node* curr = left_Corner;
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
	
	Wall_Node* curr = left_Corner;
	Wall_Node* orig = curr;
	do {
		curr->update_Angle();
		curr = curr->get_Left_Neighbor();
	} while (curr != orig);	

	return;
}

void Cell::update_Wall_Equi_Angles() {
	double new_equi_angle = (num_wall_nodes-2)*pi/num_wall_nodes;
	Wall_Node* curr = left_Corner;
	Wall_Node* orig = left_Corner;
	
	do {
		curr->update_Equi_Angle(new_equi_angle);
		curr = curr->get_Left_Neighbor();
	} while (curr != orig);	
	
	return;
}

void Cell::update_Cell_Center() {
	Wall_Node* curr = left_Corner;
	Wall_Node* orig = curr;
	Coord total_location = Coord();

	do {
		total_location += curr->get_Location();
		curr = curr->get_Left_Neighbor();
	} while(curr != orig);

	cell_center = total_location*(1.0/num_wall_nodes);
	
	return;
}

void Cell::update_Life_Length() {
	life_length++;
	return;
}
void Cell::wall_Node_Check() {
	
	//check if cell can add a cyt or wall node

	if (life_length % ADD_WALL_TIMER == ADD_WALL_TIMER-1) {
	//	cout << "adding a wall node" << endl;
		add_Wall_Node();
	}
	return;
}
void Cell::cytoplasm_Check() {

	if (life_length % ADD_CYT_TIMER == ADD_CYT_TIMER-1) {
		//cout << "adding cyt node" << endl;
		add_Cyt_Node();
	}

	return;
}

void Cell::stretch() {
	//stretch the cell
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	double y_coord = cell_center.get_Y();
	double curr_coord;
	do {
		curr_coord = curr->get_Location().get_Y();
		next = curr->get_Left_Neighbor();
		curr->pull_node();
	//	cout << "pulling" << endl;
		curr = next;
			
	} while (next != orig);
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

int Cell::update_VTK_Indices(int& id) {
//	cout << "ID before: " << id << endl;
	int rel_cnt = 0;

	Wall_Node* curr_wall = left_Corner;
	do { 
		curr_wall->update_VTK_Id(id);
		id++;
		if(curr_wall->get_Closest()!= NULL) {
			rel_cnt++;
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while (curr_wall != left_Corner);
	
	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
			cyt_nodes.at(i)->update_VTK_Id(id);
			id++;
	}
//	cout << "ID after: " << id << endl;
	return rel_cnt;
}
void Cell::print_VTK_Adh(ofstream& ofs) {

	int my_id, nei_id;
	Wall_Node* neighbor = NULL;
	Wall_Node* curr_wall = left_Corner;

	do {
		neighbor = curr_wall->get_Closest();
		if(neighbor != NULL) {
			my_id = curr_wall->get_VTK_Id();
			nei_id = neighbor->get_VTK_Id();
			ofs << 2 << ' ' << my_id << ' ' << nei_id << endl;
		}
		curr_wall = curr_wall->get_Left_Neighbor();
	} while(curr_wall != left_Corner);
	return;
}
void Cell::print_VTK_Points(ofstream& ofs, int& count) {

	Wall_Node* curr_wall = left_Corner;
	Wall_Node* orig = curr_wall;
	//	cout << "knows left corner" << endl;
	do {
		Coord loc = curr_wall->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
		//cout<< "maybe cant do left neighbor" << endl;
		curr_wall = curr_wall->get_Left_Neighbor();
		count++;
		//cout << "did it  " << count << endl;
	} while (curr_wall != orig);
	
	//cout << "walls worked" << endl;
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord loc = cyt_nodes.at(i)->get_Location();
		ofs << loc.get_X() << ' ' << loc.get_Y() << ' ' << 0 << endl;
		count++;
	}
	//cout << "points worked" << endl;
	return;
}

void Cell::print_VTK_Scalars(ofstream& ofs) {

	Wall_Node* curr_wall = left_Corner;
	do {
		Coord force = curr_wall->get_CytForce();
		ofs << force.length() << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);

	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		Coord force = cyt_nodes.at(i)->get_Force();
		ofs << force.length() << endl;
	}

	return;
}

void Cell::print_VTK_Vectors(ofstream& ofs) {

	Wall_Node* curr_wall = left_Corner;
	do {
		Coord force = curr_wall->get_CytForce();
		ofs << force.get_X() << ' ' << force.get_Y() << ' ' << 0 << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while(curr_wall != left_Corner);

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
	
	Wall_Node* curr = left_Corner;
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

//	cout << "Cell " << rank << " -- big gaps: " << big_gaps << endl;

	return biggest;
}


void Cell::add_Wall_Node() {
	//find node to the right of largest spring
	Wall_Node* right = find_Largest_Length();
	Wall_Node* left = NULL;
	Coord location;
	Wall_Node* added_node = NULL;
	if(right != NULL) {
	//	cout << "wasnt null" << endl;
		left = right->get_Left_Neighbor();
		location  = (right->get_Location() + left->get_Location())*0.5;
		added_node = new Wall_Node(location, this, left, right);
		//cout << "made new node" << endl;
		right->set_Left_Neighbor(added_node);
		left->set_Right_Neighbor(added_node);
		num_wall_nodes++;
		update_Wall_Equi_Angles();
		update_Wall_Angles();
	}
	else {
		//cout << "null" << endl;
	}
	
	return;
}

void Cell::add_Cyt_Node() {

	Cyt_Node* cyt = new Cyt_Node(cell_center, this);
	cyt_nodes.push_back(cyt);

	num_cyt_nodes++;
	return;
}

void Cell::add_Cyt_Node_Div(double& radius_x, double& radius_y) {
	//USING POSITIONS OF CELL CENTER FOR CYT NODE ALLOCATION
	// ---distributes more evenly throughout start cell
	double offset = 0.3;
	Coord location;
	Cyt_Node* cyt;
	double x;
	double y;
	
	double rand_radius_x = (static_cast<double>(rand()) / RAND_MAX)*offset*radius_x;
	double rand_radius_y = (static_cast<double>(rand()) / RAND_MAX)*offset*radius_y;
	double rand_angle = (static_cast<double>(rand()) / RAND_MAX)*2*pi;
	x = cell_center.get_X()+ rand_radius_x*cos(rand_angle);
	y = cell_center.get_Y()+ rand_radius_y*sin(rand_angle);
	location = Coord(x,y);
	cyt = new Cyt_Node(location,this);
	cyt_nodes.push_back(cyt);
	num_cyt_nodes++;
	return;
}

void Cell::most_Left_Right(Wall_Node*& left, Wall_Node*& right) {
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	
	double x_Coord = cell_center.get_X();
	double curr_Coord;
	double x_max = x_Coord;
	double x_min = x_Coord;
	
	do {
		curr_Coord = curr->get_Location().get_X();
		//cout << curr_Coord << "is curr coord" << endl;
		next = curr->get_Left_Neighbor();
		if(curr_Coord > x_max) {
			x_max = curr_Coord;
			right = curr;
			//cout << "right is: " << right << endl;

		}
		else if(curr_Coord < x_min) {
			x_min = curr_Coord;
			left = curr;

		}
		curr = next;
	} while (next != orig);
	return;
}

	
double Cell::extensional_Length() {
 
	double length = (most_right->get_Location() - most_left->get_Location()).length();

	return length;
}

void Cell::closest_node_top(Wall_Node*& up) {
	double y_coord = cell_center.get_Y();
	double x_coord = cell_center.get_X();

	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;

	double curr_dist;
	double curr_y;
	double curr_x;
	double window = .15;
	double smallest_dist = 100;
	Coord curr_coord;
	Wall_Node* closest = NULL;
	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		if((curr_x < x_coord + window) && (curr_x > x_coord - window)) {
			if(curr_y > y_coord) {
				curr_dist = (cell_center - curr_coord).length();
				if(curr_dist < smallest_dist) {
					smallest_dist = curr_dist;
					closest = curr;
				}
			}
		}
		curr = next;
	} while (next != orig);

	up = closest;
	return;
}

void Cell::closest_node_bottom(Wall_Node*& down) {
	double y_coord = cell_center.get_Y();
	double x_coord = cell_center.get_X();

	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;

	double curr_dist;
	double curr_y;
	double curr_x;
	double window = .15;
	double smallest_dist = 100;
	Coord curr_coord;
	Wall_Node* closest = NULL;
	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		if((curr_x < x_coord + window) && (curr_x > x_coord - window)) {
			if(curr_y < y_coord) {
				curr_dist = (cell_center - curr_coord).length();
				if(curr_dist < smallest_dist) {
					smallest_dist = curr_dist;
					closest = curr;
				}
			}
		}
		curr = next;
	} while (next != orig);

	down = closest;
	return;
}

void Cell::closest_node_left(Wall_Node*& left) {
	double y_coord = cell_center.get_Y();
	double x_coord = cell_center.get_X();

	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;

	double curr_dist;
	double curr_y;
	double curr_x;
	double window = .15;
	double smallest_dist = 100;
	Coord curr_coord;
	Wall_Node* closest = NULL;
	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		if((curr_y < y_coord + window) && (curr_y > y_coord - window)) {
			if(curr_x < x_coord) {
				curr_dist = (cell_center - curr_coord).length();
				if(curr_dist < smallest_dist) {
					smallest_dist = curr_dist;
					closest = curr;
				}
			}
		}
		curr = next;
	} while (next != orig);

	left = closest;
	return;
}

void Cell::closest_node_right(Wall_Node*& right) {
	double y_coord = cell_center.get_Y();
	double x_coord = cell_center.get_X();

	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;

	double curr_dist;
	double curr_y;
	double curr_x;
	double window = .15;
	double smallest_dist = 100;
	Coord curr_coord;
	Wall_Node* closest = NULL;
	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		if((curr_y < y_coord + window) && (curr_y > y_coord - window)) {
			if(curr_x > x_coord) {
				curr_dist = (cell_center - curr_coord).length();
				if(curr_dist < smallest_dist) {
					smallest_dist = curr_dist;
					closest = curr;
				}
			}
		}
		curr = next;
	} while (next != orig);

	right = closest;
	return;
}
void Cell::most_Up_Down(Wall_Node*& up, Wall_Node*& down) {
	update_Cell_Center();
	
	double y_Coord = cell_center.get_Y();
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	double y_max = y_Coord;
	double y_min = y_Coord;
	double curr_Coord;
	
	do {
		curr_Coord = curr->get_Location().get_Y();
		next = curr->get_Left_Neighbor();
		if(curr_Coord > y_max) {
			y_max = curr_Coord;
			up = curr;
		}
		else if(curr_Coord < y_min) {
			y_min = curr_Coord;
			down = curr;
		}
		curr = next;
	} while (next != orig);
	
	return;
}

double Cell::tensile_Length() {
	double length = 0;
	Wall_Node* curr = most_up;
	Wall_Node* next = NULL;
	do {

		next = curr->get_Left_Neighbor();
		length += (next->get_Location()-curr->get_Location()).length();
		curr = next;
	} while (next != most_down);
	
	return length;
}

double Cell::total_Force() {
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	Coord force_sum;
	double force;

	do {
		force_sum += curr->get_f_EXT();
		next = curr->get_Left_Neighbor();
		curr = next;
	} while (next != curr);

	force = force_sum.length();

	return force;
}

void Cell::add_strain(double& new_length) {
	strain_vec.push_back(new_length);
	return;
}

void Cell::add_stress(double& new_length, double& new_force) {
	double curr_stress = new_force/new_length;
	stress_vec.push_back(curr_stress);
	return;
}

double Cell::calc_Area() {
	Wall_Node* left = NULL;
	Wall_Node* right = NULL;
	Wall_Node* up = NULL;
	Wall_Node* down = NULL;
	this->closest_node_top(up);
	this->closest_node_bottom(down);
	//cout << "got up down" << endl;
	this->closest_node_left(left);
	this->closest_node_right(right);
	//cout << "got left right" << endl;
	//cout << "right " << right << endl;
	//cout << "left " << left << endl;
	//cout << "up" << up << endl;
	//cout << "down " << down << endl;
	//if(right == NULL) {
		//exit(1);
	//}
	double width = (right->get_Location() - left->get_Location()).length();
	//cout << "width computed" << endl;
	double length = (up->get_Location() - down->get_Location()).length();
	//cout << "length computed" << endl;
	double area = pi*width*length;

	return area;
}



	




