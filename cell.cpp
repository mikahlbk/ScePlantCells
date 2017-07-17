//cell.cpp
//===================
// Forward Declarations

//===================
// Include Dependencies

//===================

// Cell Class Member functions

// Constructors

Cell::Cell() {}

Cell::Cell(string filename) {


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


