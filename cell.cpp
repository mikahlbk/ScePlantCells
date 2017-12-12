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
#include <random>

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
	//just divided so reset life length
	life_length = 0;
	//rank assigned at time of division
	//num_cyt_nodes assigned at time of division
	//starts at 0
	num_cyt_nodes = 0;
	//layer inherited from parent cell
	//growth rate determined based on new cell center
	//at time of division
	//cell center computed at time of division
	//num_wall nodes computed at time of division
	num_wall_nodes = 0;
	//left_Corner assigned at time of division
	//WUS computed based on new cell center
	wuschel = 0;
	//cytokinin computed based on new cell center
	cytokinin = 0;
	Cell_Progress = 0;
	Cell_Progress_add_node = 0;
}


Cell::Cell(int rank, Coord center, double radius, Tissue* tiss, int layer)    {

	this->rank = rank;
//	cout << "rank " << rank << endl; 
	this->my_tissue = tiss;
	num_cyt_nodes = 0;
	this->layer = layer;
//	int init_radius = radius;
	this->cell_center = center;
//	cout << "length: " << cell_center.get_Y() << endl;
	num_wall_nodes = 0;
	life_length = 0;
	Cell_Progress = 0;
	Cell_Progress_add_node = 0;
	this->calc_WUS();
	this->calc_CYT();
	double K_LINEAR_Y = .1540*pow(wuschel,3) + -4.8350*pow(wuschel,2) + 54.2901*wuschel + -50.7651;
	double K_LINEAR_X = -13.2177*wuschel + 473.7440;
	if(K_LINEAR_Y > 1000) {
		K_LINEAR_Y = 1000;
	}
	if(K_LINEAR_X > 1000) {
		K_LINEAR_X = 1000;
	}
	this->K_LINEAR = Coord(K_LINEAR_X, K_LINEAR_Y);
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
	//set damping for cells that act as anchor points
	if(layer == 6) {
		this->damping = .01;
	}
	else {
		this->damping = 1;
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

void Cell::set_Damping(double& new_damping) {
	this->damping = new_damping;
	return;
}
void Cell::set_Layer(int layer) {
	this->layer = layer;
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
	return num_wall_nodes + num_cyt_nodes;
}

void Cell::set_Left_Corner(Wall_Node*& new_left_corner) {
	this->left_Corner = new_left_corner;
	return;
}

void Cell::set_Wall_Count(int& number_nodes) {
	this->num_wall_nodes = number_nodes;
	return;
}
/*void Cell::set_Area(double& new_area){
	area = new_area;
	return;
}*/
void Cell::reset_Cell_Progress(){
	this->Cell_Progress = 0;

	return;
}
void Cell::set_Cell_Progress_add_node(double& new_area) {
	this->Cell_Progress_add_node = new_area;
	return;
}
void Cell::calc_WUS() {
	this->wuschel = 0.0246*pow(cell_center.length(),2) + -1.8659*cell_center.length() + 28.5785;
	//cout << "wuschel: " << wuschel << endl;
	return;
}
void Cell::calc_CYT() {
	this->cytokinin = -.0713*pow(cell_center.length()*.1,2) + .10761*cell_center.length()*.1 + 10.6624; 
	//	cout << "cytokinin: " << cytokinin << endl;
	return;
}

//=============================================================
//=========================================
// Keep Track of neighbor cells and Adhesion springs
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

void Cell::update_adhesion_springs() {
//	int time = Ti;
	vector<Cell*>neighbors;
	this->get_Neighbor_Cells(neighbors);
	Wall_Node* curr_Node = NULL;
	Wall_Node* next_Node = NULL;
	Wall_Node* curr_Closest = NULL;
	double curr_len = 0;
//	for(int i = 0; i<num_wall_nodes;i++) {
	curr_Node = left_Corner;
	do {
		next_Node = curr_Node->get_Left_Neighbor();
		curr_Closest = curr_Node->find_Closest_Node(neighbors);
		curr_Node->make_Connection(curr_Closest);
		curr_Node = next_Node;
	} while(next_Node != left_Corner);
	//}
	return;
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
	double new_damping = 0;
	for (unsigned int i = 0; i < cyt_nodes.size(); i++) {
		new_damping = cyt_nodes.at(i)->get_My_Cell()->get_Damping();
		cyt_nodes.at(i)->update_Location(new_damping);

	}

	//update wall nodes
	Wall_Node* curr = left_Corner;
	Wall_Node* orig = curr;

	do {
		new_damping = curr->get_My_Cell()->get_Damping();
		curr->update_Location(new_damping);
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
//	cout << life_length << endl;
	return;
}
void Cell::wall_Node_Check() {
	//cout << "adding a wall node" << endl;
	add_Wall_Node();
	return;
}
//======================================================
//==================================
//Functions for calibration
//==================================
//======================================================
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

double Cell::extensional_Length() {
 	Wall_Node* right = NULL;
	Wall_Node* left = NULL;
	closest_node_right(right);
	closest_node_left(left);
	double length = (right->get_Location() - left->get_Location()).length();

	return length;
}

double Cell::tensile_Length() {
	double length = 0;
	Wall_Node* up = NULL;
	Wall_Node* bottom = NULL;
	closest_node_bottom(bottom);
	closest_node_top(up);
	Wall_Node* curr = up;
	Wall_Node* next = NULL;
	do {

		next = curr->get_Left_Neighbor();
		length += (next->get_Location()-curr->get_Location()).length();
		curr = next;
	} while (next != bottom);
	
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

void Cell::print_VTK_Scalars_Force(ofstream& ofs) {

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

void Cell::print_VTK_Scalars_WUS(ofstream& ofs) {

	//double concentration = 0;
	Wall_Node* curr_wall = left_Corner;
	do {
		double concentration = curr_wall->get_My_Cell()->get_WUS_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);


	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		double concentration = cyt_nodes.at(i)->get_My_Cell()->get_WUS_concentration();
		ofs << concentration << endl;
	}
	return;
}

void Cell::print_VTK_Scalars_CYT(ofstream& ofs) {

//	double concentration = 0;

	Wall_Node* curr_wall = left_Corner;
	do {
		double concentration = curr_wall->get_My_Cell()->get_CYT_concentration();
		ofs << concentration << endl;

		curr_wall = curr_wall->get_Left_Neighbor();
		
	} while (curr_wall != left_Corner);


	for(unsigned int i = 0; i < cyt_nodes.size(); i++) {
		double concentration = cyt_nodes.at(i)->get_My_Cell()->get_CYT_concentration();
		ofs << concentration << endl;
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
void Cell::find_Largest_Length(Wall_Node* first, Wall_Node* second) {
	
	Wall_Node* curr = left_Corner;
	Wall_Node* biggest = NULL;
	Wall_Node* orig = curr;
	Coord left_Neighb_loc;
	Coord curr_Loc;
	Coord diff_vect;
	double max_len = 0;
	double second_max = 0;
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
			if(len > second_max) {
				second_max = len;
				second = curr;
				if(len > max_len) {
					max_len = len;
					first = curr;
				}
			}
		}
		curr = curr->get_Left_Neighbor();

	} while (curr != orig);

//	cout << "Cell " << rank << " -- big gaps: " << big_gaps << endl;

	return;
}


void Cell::add_Wall_Node() {
	//find node to the right of largest spring
	Wall_Node* right = NULL;
	Wall_Node* second = NULL;
	find_Largest_Length(right, second);
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

void Cell::add_Cyt_Node_Div(double& radius_x, double& radius_y, bool islength) {
	//USING POSITIONS OF CELL CENTER FOR CYT NODE ALLOCATION
	// ---distributes more evenly throughout start cell
	double offsetx = 0;
	double offsety = 0;
	if(islength) {
		offsetx = 0.4;
		offsety = 0.8;
	}
	else {
		offsetx = 0.8;
		offsety = 0.4;
	}
	
	Coord location;
	Cyt_Node* cyt;
	double x;
	double y;
	double rand_radius_x = (static_cast<double>(rand()) / RAND_MAX)*offsetx*radius_x;
	double rand_radius_y = (static_cast<double>(rand())/ RAND_MAX)*offsety*radius_y; 
	double rand_angle = (static_cast<double>(rand()) / RAND_MAX)*2*pi;
	x = cell_center.get_X()+ rand_radius_x*cos(rand_angle);
	y = cell_center.get_Y()+ rand_radius_y*sin(rand_angle);
	location = Coord(x,y);
	cyt = new Cyt_Node(location,this);
	cyt_nodes.push_back(cyt);
	num_cyt_nodes++;
	return;
}

void Cell::closest_node_top(Wall_Node*& up) {
	//finds node with highest y value that is
	//in a certain window around x value of center
	double y_coord = cell_center.get_Y();
	double x_coord = cell_center.get_X();

	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;

	double curr_dist;
	double curr_y;
	double curr_x;
	double smallest_dist = 500;
//	cout << "smallest dist top" << endl;
	Coord curr_coord;
	Wall_Node* closest = NULL;

	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		curr_dist = sqrt(pow(x_coord -curr_x,2));
		if(curr_dist < smallest_dist) {
		//	cout << "passed smallest check" << endl;
			if(curr_y > y_coord) {
		//		cout << "passed y coord check" << endl;	
				//cout << "passed smallest check" << endl;
					smallest_dist = curr_dist;
					//cout << "curr " << curr << endl;
					closest = curr;
				}
		}
		curr = next;
	} while (next != orig);
	//cout << "closest" << closest << endl;
	up = closest;
	//cout << "up in function" << up << endl;
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
	double smallest_dist = 500;
	Coord curr_coord;
	Wall_Node* closest = NULL;

	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		curr_dist = sqrt(pow(x_coord - curr_x, 2));
		if(curr_dist < smallest_dist) {
			if(curr_y < y_coord) {
					smallest_dist = curr_dist;
					closest = curr;
			}
		}
		curr = next;
	} while (next != orig);
//	cout << "bottom" << endl;
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
	double smallest_dist = 500;
	Coord curr_coord;
	Wall_Node* closest = NULL;

	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		curr_dist = sqrt(pow(y_coord - curr_y,2));
		if(curr_dist < smallest_dist) {
			//cout << "passed window check" << endl;
			if(curr_x < x_coord){ 
					smallest_dist = curr_dist;
					closest = curr;
			}
		}
		curr = next;
	} while (next != orig);
//	cout << "left" << endl;
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
	double smallest_dist = 500;
	Coord curr_coord;
	Wall_Node* closest = NULL;
	do {
		curr_coord = curr->get_Location();
		next = curr->get_Left_Neighbor();
		curr_y = curr_coord.get_Y();
		curr_x = curr_coord.get_X();
		curr_dist = sqrt(pow(y_coord-curr_y,2));
		if(curr_dist < smallest_dist) {
			if(curr_x > x_coord) {
					smallest_dist = curr_dist;
					closest = curr;
			}
		}
		curr = next;
	} while (next != orig);
//	cout << "right" << endl;
	right = closest;
	return;
}

void Cell::update_Cell_Progress(int Ti) {
	Cell* new_Cell= NULL;
	vector<Cell*> cells;
	this->update_Life_Length();
	this->my_tissue->get_Cells(cells);
	int number_cells = cells.size();
	double sigma = (((double) rand()/(double) RAND_MAX))*.004 - .002;
//	cout << "Sigma: " << sigma << endl;
	double rate = (.002+sigma)*exp(.08*dt*life_length);
	Cell_Progress = Cell_Progress + rate*dt;
	//cout << "growth rate: " << rate << endl;
	//cout<< "new: " << Cell_Progress << endl;
	//cout << "Old: " << Cell_Progress_add_node << endl;
	if(this->calc_Area() >= 8) {
		new_Cell = this->divide();
		new_Cell->set_Rank(number_cells);
		my_tissue->update_Num_Cells(new_Cell);
		this->update_Neighbor_Cells();
		this->update_adhesion_springs();
		new_Cell->update_Neighbor_Cells();
		new_Cell->update_adhesion_springs();
	}
	else {
		if(Cell_Progress - Cell_Progress_add_node > 0.05) {
			this->add_Cyt_Node();
			//cout << "added cyt node" << endl;
			Cell_Progress_add_node = Cell_Progress;
		}	
		//this->Cell_Progress = new_Cell_Progress;
		this->update_Life_Length();
	}
	//cout << "Cell Prog: " << Cell_Progress_add_node << endl;
	return;
}

double Cell::calc_Area() {
	Wall_Node* curr = left_Corner;
	Wall_Node* next = NULL;
	Wall_Node* orig = curr;
	Coord a_i;
	Coord a_j;
	double area = 0;
	double curr_area = 0;
	do {
		next = curr->get_Left_Neighbor();
		a_i = curr->get_Location() - cell_center;
		a_j = next->get_Location() - cell_center;
		curr_area = 0.5*sqrt(pow(a_i.cross(a_j),2));
		area += curr_area;
		curr = next;
	} while(next != orig);

	return area;
}
///=================================
///==========================
//Measurements for Division
///==========================
///==================================

double Cell::compute_pressure(Wall_Node* start, Wall_Node* end) {
	Wall_Node* curr_wall = start;
	Wall_Node* next = NULL;
	Wall_Node* orig = end;
	Coord force;
	double curr_length = 0;
	double total_length = 0;
	double pressure = 0;
	do {
		next = curr_wall->get_Left_Neighbor();
		force += curr_wall->get_CytForce();
		curr_length = (curr_wall->get_Location() - next->get_Location()).length();
		total_length += curr_length;
		curr_wall = next;
	} while (next != orig);
	
	pressure = force.length()/total_length;

	return pressure;
}

double Cell::compute_sigma_long() {
	Wall_Node* top = NULL;
	this->closest_node_top(top);
	double length = ((cell_center - top->get_Location()).length())*(0.5);
	double y_cut_off = cell_center.get_Y() + length;
	Wall_Node* left = NULL;
	Wall_Node* right = NULL;
	this->closest_node_left(left);
	this->closest_node_right(right);
	Wall_Node* curr = left;
	Wall_Node* start = NULL;
	Wall_Node* end = NULL;
	bool possible = true;
	do {
		if(curr->get_Location().get_Y() > y_cut_off) {
			start = curr;
			possible = false;
		}  
		curr = curr->get_Right_Neighbor();
	} while (possible);

	curr = right;
	possible = true;
	do {
		if(curr->get_Location().get_Y() > y_cut_off) {
			end = curr;
			possible = false;
		}
		curr = curr->get_Left_Neighbor();
	} while (possible);
	double pressure_top = this->compute_pressure(start, end);
	
	Wall_Node* bottom = NULL;
	this->closest_node_bottom(bottom);
	length = ((cell_center - bottom->get_Location()).length())*(0.5);
	y_cut_off = cell_center.get_Y() - length;
//	Wall_Node* left = NULL;
//	Wall_Node* right = NULL;
//	this->most_left_node(left);
//	this->most_right_node(right);
	curr = left;
	start = NULL;
	end = NULL;
	possible = true;
	do {
		if(curr->get_Location().get_Y() < y_cut_off) {
			start = curr;
			possible = false;
		}  
		curr = curr->get_Left_Neighbor();
	} while (possible);

	curr = right;
	possible = true;
	do {
		if(curr->get_Location().get_Y() < y_cut_off) {
			end = curr;
			possible = false;
		}
		curr = curr->get_Right_Neighbor();
	} while (possible);
	double pressure_bottom = this->compute_pressure(start, end);

	double sigma_long = pressure_top + pressure_bottom;
	return sigma_long;
}

double Cell::compute_sigma_trans() {

	Wall_Node* left = NULL;
	this->closest_node_left(left);
	double length = ((cell_center - left->get_Location()).length())*(0.5);
	double x_cut_off = cell_center.get_X() + length;
	Wall_Node* top = NULL;
	Wall_Node* bottom = NULL;
	this->closest_node_top(top);
	this->closest_node_bottom(bottom);
	Wall_Node* curr = top;
	Wall_Node* start = NULL;
	Wall_Node* end = NULL;
	bool possible = true;
	do {
		if(curr->get_Location().get_X() < x_cut_off) {
			start = curr;
			possible = false;
		}  
		curr = curr->get_Left_Neighbor();
	} while (possible);

	curr = bottom;
	possible = true;
	do {
		if(curr->get_Location().get_X() < x_cut_off) {
			end = curr;
			possible = false;
		}
		curr = curr->get_Right_Neighbor();
	} while (possible);
	double pressure_left = this-> compute_pressure(start, end);

	Wall_Node* right = NULL;
	this->closest_node_right(right);
	length = ((cell_center - right->get_Location()).length())*(0.5);
	x_cut_off = cell_center.get_X() + length;
	//Wall_Node* top = NULL;
	//Wall_Node* bottom = NULL;
	//this->most_top_node(top);
	//this->most_bottom_node(bottom);
	curr = bottom;
	start = NULL;
	end = NULL;
	possible = true;
	do {
		if(curr->get_Location().get_X() < x_cut_off) {
			start = curr;
			possible = false;
		}  
		curr = curr->get_Left_Neighbor();
	} while (possible);

	curr = top;
	possible = true;
	do {
		if(curr->get_Location().get_X() < x_cut_off) {
			end = curr;
			possible = false;
		}
		curr = curr->get_Right_Neighbor();
	} while (possible);
	double pressure_right = this->compute_pressure(start, end);


	double sigma_trans = pressure_left + pressure_right;
	return sigma_trans;
}



