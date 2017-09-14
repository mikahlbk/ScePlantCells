//cell_div.cpp
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

Cell* Cell::divide(const int Ti) {

	Cell* sister = NULL;
    //calculate area
	double area = 0;
	Coord x_vec;
	Coord y_vec;
	double x_length = 0;
	double y_length = 0;
	y_vec = corners.at(3)->get_Location() - corners.at(0)->get_Location();
	x_vec = corners.at(1)->get_Location() - corners.at(0)->get_Location();
	y_length = y_vec.length();
	x_length = x_vec.length();
	area = y_length*x_length;
	if ((area >  6)) {
			if ((this->rank == 3) || (this->rank == 4) || (this->rank == 5) || (this->rank  == 6)) {
				sister = this -> divide_length_wise(Ti);
			}
			else {
				if (y_length > x_length) {
					sister = this->divide_width_wise(Ti);
			    }		
				else {
					sister = this->divide_length_wise(Ti);
				}		
			}			
	}


	return sister;
}

Cell* Cell::divide_length_wise(const int Ti) {
	// Current cell splits into two daughter cells
	//    -"this" will keep its entity as the left sister, 
	//		create a sister cell to the right of it, and return it 
	//      to the tissue.

	//ask tissue for new id num
	Cell* sister = new Cell(Ti, my_tissue);

	// Find midpoint
	Coord mid_top = ((corners.at(3)->get_Location() + corners.at(2)->get_Location()) * 0.5);
	Coord mid_bot = ((corners.at(0)->get_Location() + corners.at(1)->get_Location()) * 0.5);
	Coord divider = ((mid_top + mid_bot) / 2);
	double divide_X = divider.get_X();
	// vectors and counters for each cell
	vector<Wall_Node*> corners_A;
	vector<Wall_Node*> corners_B;
	int wall_cnt_A = 0;
	int wall_cnt_B = 0;

	// For iterating through all wall_nodes
	Wall_Node* curr = corners.at(0);
	Wall_Node* next = NULL;
	Wall_Node* prev = NULL;
	int side = 0;

	cout << "BEGIN ITERATING THROUGH WALL NODES" << endl;

	// distribute wall nodes between sister cells
	do { 

		next = curr->get_Left_Neighbor();
		
		if (curr->get_Location().get_X() < divide_X) {
			//Cell A
			
			if (curr->is_Corner()) {
				corners_A.push_back(curr);
				side++;
			}

			if (side == 1) {
				//on bottom end
				
				if (next->get_Location().get_X() < divide_X) {
					prev = curr;
					wall_cnt_A++;
				}
				else {
					// Create new corner for cell A
					Wall_Node* prev_prev = prev->get_Right_Neighbor();
					Wall_Node* temp = new Corner_Node(prev->get_Location(), this);
					temp->set_Right_Neighbor(prev_prev);
					prev_prev->set_Left_Neighbor(temp);
					corners_A.push_back(temp);
					// delete curr and prev
					delete curr;
					delete prev;
					// deleted curr and prev but curr was never counted
					// and added temp so dont incrememnt
				}

			}
			else if (side == 3) {
				//on top end
				if (corners_A.size() < 3) {
					//need to delete curr and next, then insert a corner
					//WARNING: prev is set to null, do not try to access
					Wall_Node* next_next = next->get_Left_Neighbor();
					Wall_Node* temp = new Corner_Node(next->get_Location(), this);
					temp->set_Left_Neighbor(next_next);
					next_next->set_Right_Neighbor(temp);
					corners_A.push_back(temp);
					//delete curr and prev
					delete curr;
					delete next;
					// reset curr and next
					prev = temp;
					next = next_next;
					// increment because added temp. never counted curr or next
					wall_cnt_A++;
				}
				else {
					prev = curr;
					wall_cnt_A++;
				}
			}
			else if (side == 2) {
				cout << "ERROR: DIVISION. CAN'T BE ON THIS FLANK" << endl;
			}
			else { //side == 4
				prev = curr;
				wall_cnt_A++;
			}
		}
		else { //on right side of dividing line
			//Cell B
			if (curr->is_Corner()) {
				corners_B.push_back(curr);
				side++;
			}

			if (side == 1) {
				if (corners_B.size() < 1) {
					//need to delete curr and next, then insert a corner
					//WARNING: prev is undefined, do not try to access
					Wall_Node* next_next = next->get_Left_Neighbor();
					Wall_Node* temp = new Corner_Node(next->get_Location(), sister);
					temp->set_Left_Neighbor(next_next);
					next_next->set_Right_Neighbor(temp);
					corners_B.push_back(temp);
					//delete
					delete curr;
					delete next;
					// reset next
					prev = temp;
					next = next_next;
					// Never added curr or next, only temp so add 1
					wall_cnt_B++;

				}
				else {
					curr->set_My_Cell(sister);
					prev = curr;
					wall_cnt_B++;
				}
			}
			else if (side == 2) {
				curr->set_My_Cell(sister);
				prev = curr;
				wall_cnt_B++;
			}
			else if (side == 3) {
				//on top end		
				if (next->get_Location().get_X() > divide_X) {
					curr->set_My_Cell(sister);
					prev = curr;
					wall_cnt_B++;
				}
				else {
					// delete curr and prev, then add corner
					Wall_Node* prev_prev = prev->get_Right_Neighbor();
					Wall_Node* temp = new Corner_Node(prev->get_Location(), sister);
					temp->set_Right_Neighbor(prev_prev);
					prev_prev->set_Left_Neighbor(temp);
					corners_B.push_back(temp);
					// delete
					delete curr;
					delete prev;
					// incr because of temp, but decremetn because deleted prev
					//dont change count
				}	
			}
			else { //side == 4
				cout << "ERROR: DIVISION. CAN'T BE ON THIS FLANK" << endl;
			}

		}

		curr = next;

	} while(curr != corners.at(0));

	cout << "FINISHED ITERATING THROUGH WALL NODES" << endl;

	// create new flank walls for both cells
	double spacing = 0.08;
	
	//Cell A
	Wall_Node* x = corners_A.at(1);
	Wall_Node* y = corners_A.at(2);
	Coord diff = (y->get_Location() - x->get_Location());
	double total_length = diff.length();
	int num_needed_nodes = (total_length / spacing) - 1;
	Coord slope_vect = ((diff / total_length) * spacing);
	Coord curr_pos = (x->get_Location() + slope_vect);

	for (int i = 0; i < num_needed_nodes; i++) {
		curr = new Flank_Node(curr_pos, this);
		wall_cnt_A++;
		// set neighbors
		curr->set_Right_Neighbor(x);
		x->set_Left_Neighbor(curr);
		// update x and new pos
		curr_pos += slope_vect;
		x = curr;
	}
	// finish the wall
	y->set_Right_Neighbor(x);
	x->set_Left_Neighbor(y);

	
	// Cell B
	x = corners_B.at(3);
	y = corners_B.at(0);
	diff = (y->get_Location() - x->get_Location());
	total_length = diff.length();
	num_needed_nodes = (total_length / spacing) - 1;
	slope_vect = ((diff / total_length) * spacing);
	curr_pos = (x->get_Location() + slope_vect);

	for (int i = 0; i < num_needed_nodes; i++) {
		curr = new Flank_Node(curr_pos, sister);
		wall_cnt_B++;
		// set neighbors
		curr->set_Right_Neighbor(x);
		x->set_Left_Neighbor(curr);
		// update x and new pos
		curr_pos += slope_vect;
		x = curr;
	}
	// finish the wall
	y->set_Right_Neighbor(x);
	x->set_Left_Neighbor(y);

	cout << "FINISHED TWO EXTRA FLANK WALLS" << endl;

	//update each cell's private member fields
	//   This
	num_wall_nodes = wall_cnt_A;
	corners = corners_A;

	//   Sister
	sister->set_Wall_Cnt(wall_cnt_B);
	sister->set_CornerNodes(corners_B);

	//update angles and cell centers for both cells
	this->update_Wall_Angles();
	cout << "update THIS wall angles" << endl;
	sister->update_Wall_Angles();
	cout << "updated wall angles" << endl;
	this->update_Cell_Center();
	sister->update_Cell_Center();

	cout << "updated angles and cell centers" << endl;

	// distribute cyt nodes between sister cells
	int new_cyt_cnt = 18;
	// Delete all old cyt nodes
	Cyt_Node* c = NULL;
	while ( !cyt_nodes.empty() ) {
		c = cyt_nodes.at(cyt_nodes.size() - 1);
		delete c;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}
	cout << "FINISHED DELETING OLD CYT NODES" << endl;
	//create new ones for each cell
	for (int i = 0; i < new_cyt_cnt; i++) {
		//add_Cyt_Node take in int value to check if
		// it's proper time to add a node
		// pass in times to work
		this->add_Cyt_Node(init_cell_time);
		sister->add_Cyt_Node(Ti);
	}

	
	return sister;
}


// Gets called by update_Node_Locations()
Cell* Cell::divide_width_wise(const int Ti) {
	// Current cell splits into two daughter cells
	//    -"this" will keep its entity, create a sister cell, and return it 
	//      to the tissue.
	Cell* sister = new Cell(Ti, my_tissue);

	// Find midpoint
	Coord mid_left = ((corners.at(0)->get_Location() + corners.at(3)->get_Location()) * 0.5);
	Coord mid_right = ((corners.at(1)->get_Location() + corners.at(2)->get_Location()) * 0.5);
	Coord divider = ((mid_left + mid_right) / 2);
	double divide_Y = divider.get_Y();
	// vectors and counters for each cell
	vector<Wall_Node*> corners_A;
	vector<Wall_Node*> corners_B = {NULL, NULL, NULL, NULL};
	int wall_cnt_A = 0;
	int wall_cnt_B = 0;

	// For iterating through all wall_nodes
	Wall_Node* curr = corners.at(0);
	Wall_Node* next = NULL;
	Wall_Node* prev = NULL;
	int side = 0;

	cout << "BEGIN ITERATING THROUGH WALL NODES" << endl;

	// distribute wall nodes between sister cells
	do { 

		next = curr->get_Left_Neighbor();
		
		if (curr->get_Location().get_Y() < divide_Y) {
			//Cell A
			
			if (curr->is_Corner()) {
				corners_A.push_back(curr);
				side++;
			}

			if (side == 1) {
				//on bottom end
				prev = curr;
				wall_cnt_A++;	
			}
			else if (side == 3) {
				//on top end
				cout << "ERROR: DIVISION. CAN'T BE ON THIS FLANK" << endl;
			}
			else if (side == 2) {

				if (next->get_Location().get_Y() < divide_Y) {
					prev = curr;
					wall_cnt_A++;
				}
				else {
					// Create new corner for cell A
					Wall_Node* prev_prev = prev->get_Right_Neighbor();
					Wall_Node* temp = new Corner_Node(prev->get_Location(), this);
					temp->set_Right_Neighbor(prev_prev);
					prev_prev->set_Left_Neighbor(temp);
					corners_A.push_back(temp);
					// delete curr and prev
					delete curr;
					delete prev;
					// deleted curr and prev but curr was never counted
					// and added temp so dont incrememnt
				}


			}
			else { //side == 4
				if (corners_A.size() < 4) {
					//need to delete curr and next, then insert a corner
					//WARNING: prev is set to null, do not try to access
					Wall_Node* next_next = next->get_Left_Neighbor();
					Wall_Node* temp = new Corner_Node(next->get_Location(), this);
					temp->set_Left_Neighbor(next_next);
					next_next->set_Right_Neighbor(temp);
					corners_A.push_back(temp);
					//delete curr and prev
					delete curr;
					delete next;
					// reset curr and next
					prev = temp;
					next = next_next;
					// increment because added temp. never counted curr or next
					wall_cnt_A++;
				}
				else {
					prev = curr;
					wall_cnt_A++;
				}
			}
		}
		else { //above the dividing line
			//Cell B
			if (side == 1) {
				// bottom end
				cout << "ERROR: DIVISION. CAN'T BE ON THIS FLANK" << endl;
			}
			else if (side == 2) {
				//right flank
				if (corners_B.at(1) == NULL) {
					//need to delete curr and next, then insert a corner
					//WARNING: prev is undefined, do not try to access
					Wall_Node* next_next = next->get_Left_Neighbor();
					Wall_Node* temp = new Corner_Node(next->get_Location(), sister);
					temp->set_Left_Neighbor(next_next);
					next_next->set_Right_Neighbor(temp);
					corners_B.at(1) = temp;
					//delete
					delete curr;
					delete next;
					// reset next
					prev = temp;
					next = next_next;
					// Never added curr or next, only temp so add 1
					wall_cnt_B++;
				}
				else {
					if (curr->is_Corner()) {
						corners_B.at(2) = curr;
						side++;
					}
					curr->set_My_Cell(sister);
					prev = curr;
					wall_cnt_B++;
				}
			}
			else if (side == 3) {
				// top end
				if (curr->is_Corner()) {
					corners_B.at(3) = curr;
					side++;
				}
				curr->set_My_Cell(sister);
				prev = curr;
				wall_cnt_B++;
			}
			else { //side == 4
				// left flank
				if (next->get_Location().get_Y() > divide_Y) {
					curr->set_My_Cell(sister);
					prev = curr;
					wall_cnt_B++;
				}
				else {
					// delete curr and prev, then add corner
					Wall_Node* prev_prev = prev->get_Right_Neighbor();
					Wall_Node* temp = new Corner_Node(prev->get_Location(), sister);
					temp->set_Right_Neighbor(prev_prev);
					prev_prev->set_Left_Neighbor(temp);
					corners_B.at(0) = temp;
					// delete
					delete curr;
					delete prev;
					// incr because of temp, but decremetn because deleted prev
					//dont change count
				}	
			}
		}

		curr = next;

	} while(curr != corners.at(0));

	cout << "FINISHED ITERATING THROUGH WALL NODES" << endl;

	cout << "CornersA.size() = " << corners_A.size() << endl;
	for (unsigned int i = 0; i < corners_A.size(); i++) {
		if (corners_A.at(i) == NULL) {
			cout << "   A(" << i << ") is NULL" << endl;
		}
		else {
			cout << "   A(" << i << ") is good" << endl;
		}
	}
	cout << "CornersB.size() = " << corners_B.size() << endl;
	for (unsigned int i = 0; i < corners_B.size(); i++) {
		if (corners_B.at(i) == NULL) {
			cout << "   B(" << i << ") is NULL" << endl;
		}
		else {
			cout << "   B(" << i << ") is good" << endl;
		}
	}


	// create new flank walls for both cells
	curr = corners.at(0)->get_Left_Neighbor();
	int num_needed_nodes = 0;
	do {
		num_needed_nodes++;
		curr = curr->get_Left_Neighbor();
	} while ( !(curr->is_Corner()) );

	cout << "starting end walls" << endl;
	//Cell A
	Wall_Node* x = corners_A.at(2);
	Wall_Node* y = corners_A.at(3);
	Coord diff = (y->get_Location() - x->get_Location());
	double total_length = diff.length();
	double spacing = total_length / (num_needed_nodes + 1);
	Coord slope_vect = ((diff / total_length) * spacing);
	Coord curr_pos = (x->get_Location() + slope_vect);

	for (int i = 0; i < num_needed_nodes; i++) {
		curr = new End_Node(curr_pos, this);
		wall_cnt_A++;
		// set neighbors
		curr->set_Right_Neighbor(x);
		x->set_Left_Neighbor(curr);
		// update x and new pos
		curr_pos += slope_vect;
		x = curr;
	}
	// finish the wall
	y->set_Right_Neighbor(x);
	x->set_Left_Neighbor(y);

	cout << "finished A" << endl;

	// Cell B
	x = corners_B.at(0);
	y = corners_B.at(1);
	diff = (y->get_Location() - x->get_Location());
	total_length = diff.length();
	spacing = total_length / (num_needed_nodes + 1);
	slope_vect = ((diff / total_length) * spacing);
	curr_pos = (x->get_Location() + slope_vect);

	for (int i = 0; i < num_needed_nodes; i++) {
		curr = new End_Node(curr_pos, sister);
		wall_cnt_B++;
		// set neighbors
		curr->set_Right_Neighbor(x);
		x->set_Left_Neighbor(curr);
		// update x and new pos
		curr_pos += slope_vect;
		x = curr;
	}
	// finish the wall
	y->set_Right_Neighbor(x);
	x->set_Left_Neighbor(y);

	cout << "FINISHED TWO EXTRA END WALLS" << endl;

	//update each cell's private member fields
	//   This
	num_wall_nodes = wall_cnt_A;
	corners = corners_A;
	cout << "Check A" << endl;

	//   Sister
	sister->set_Wall_Cnt(wall_cnt_B);
	cout << "check here" << endl;
	sister->set_CornerNodes(corners_B);
	cout << "Check B" << endl;

	//update angles and cell centers for both cells
	this->update_Wall_Angles();
	cout << "update THIS wall angles" << endl;
	sister->update_Wall_Angles();
	cout << "updated wall angles" << endl;
	this->update_Cell_Center();
	sister->update_Cell_Center();

	cout << "updated angles and cell centers" << endl;

	// distribute cyt nodes between sister cells
	int new_cyt_cnt = 18;
	// Delete all old cyt nodes
	Cyt_Node* c = NULL;
	while ( !cyt_nodes.empty() ) {
		c = cyt_nodes.at(cyt_nodes.size() - 1);
		delete c;
		cyt_nodes.pop_back();
		num_cyt_nodes--;
	}
	cout << "FINISHED DELETING OLD CYT NODES" << endl;
	//create new ones for each cell
	for (int i = 0; i < new_cyt_cnt; i++) {
		//add_Cyt_Node take in int value to check if
		// it's proper time to add a node
		// pass in times to work
		this->add_Cyt_Node(init_cell_time);
		sister->add_Cyt_Node(Ti);
	}


	return sister;
}



//end file
