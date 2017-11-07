//tissue.cpp
//=========================

//=========================
//Include Dependencies
#include <vector>
#include <fstream>
#include <sstream>
#include <string>

#include "phys.h"
#include "coord.h"
#include "cell.h"
#include "tissue.h"
#include "side.h"
//=========================
// Public Member Functions for Tissue.cpp

Tissue::Tissue(string filename) {
	num_cells = 0;
	ifstream ifs(filename.c_str());

	if(!ifs) {
		cout << filename << " is not available" << endl;
		return;
	}

	stringstream ss;
	string line;
	string temp;
	char trash;
	int rank;
	int layer;
	double height, width;
	Coord corner;
	double x, y;
	Cell* curr;

	while (getline(ifs,line)) {
		ss.str(line);

		getline(ss,temp,':');

		if (temp == "CellRank") {
			ss >> rank;
		}
		else if (temp == "Corner") {
			ss >> x >> trash >> y;
			Coord loc(x,y);
			corner = loc;
		}
		else if (temp == "Height") {
			ss >> height;
		}
		else if (temp == "Width") {
			ss >> width;
		}
		else if (temp == "Layer") {
			ss >> layer;
		}
		else if (temp == "End_Cell") {
			//create new cell with collected data and push onto vector 
			curr = new Cell(rank, corner, height, width, 0, this, layer);
			num_cells++;
			cells.push_back(curr);
		}

		ss.clear();
	}

	ifs.close();
	this->num_cells = num_cells;
}

Tissue::~Tissue() {
	
	Cell* curr = NULL;
	while ( !cells.empty() ) {
		curr = cells.at(cells.size() - 1);
		delete curr;
		cells.pop_back();
	}
}

void Tissue::get_Cells(vector<Cell*>& cells) {
	cells = this->cells;
	return;
}

void Tissue::update_Life_Length() {
//	cout << "tissue life length check" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Life_Length();
	}
	return;
}

void Tissue::calc_New_Forces() {

	for (unsigned int i = 0; i < cells.size(); i++) {
//		cout << "calc new forces" << endl;
		cells.at(i)->calc_New_Forces();
	}

	return;
//	cout<<"calculated forces"<<endl;

}

void Tissue::stretch(int Ti) {
	Coord force = Coord(0,0);
	for (unsigned int i = 0; i< cells.size();i++) {
		cout << "tissue stretch" << endl;
		cells.at(i)->stretching(Ti, force);
	}
	return;
}

void Tissue::update_Cell_Locations() {
	
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Node_Locations();
	}

	return;
//	cout<<"Updated locations"<<endl;
}
	

void Tissue::update_Neighbor_Cells() {
	//update vectors of neighboring cells
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Neighbor_Cells();
	}
	
	return;
}

void Tissue::update_Adhesion() {
	Wall_Node* curr = NULL;
	Wall_Node* orig = NULL;
	Wall_Node* next = NULL;
	for (unsigned int i = 0;i<cells.size();i++) {
		//cout<< "clearing current closest info" << endl;
		curr = cells.at(i)->get_Wall_Nodes();
		orig = curr;
		do {
			next = curr->get_Left_Neighbor();
			curr->set_Closest(NULL, 100);
			curr = next;
		} while(next != orig);
	}
	
	for(unsigned int i=0;i<cells.size();i++) {
		//cout << "Updating adhesion for cell" << endl;
		cells.at(i)->update_adhesion_springs();
	}
}

void Tissue::cell_Division(const int Ti) {
	bool divided = false;
	Cell* new_cell = NULL;
//	cout << "number cells: " << cells.size()<< endl;
	for(unsigned int i = 0; i < cells.size();i++) {
		//cout << "current divide cell: " << i << endl;
		new_cell = cells.at(i)->divide(Ti);
		if (new_cell !=NULL) {
			divided = true;
			new_cell->set_Rank(num_cells);
			num_cells++;
			cells.push_back(new_cell);
		}

	}
	if(divided) {
		this->update_Neighbor_Cells();
		this->update_Adhesion();
	}
	return;
}

void Tissue::make_Vectors() {
	Cell* curr= NULL;
	vector<double>curr_lengths;
	vector<Side*>sides;
	vector<double>curr_forces;
	Side* curr_side = NULL;
	ofstream myfile1("cell_vec.txt");
	ofstream myfile2("side_vecs.txt");
	if(myfile1.is_open()) {
		for(int i = 0; i < cells.size();i++) {
			curr = cells.at(i);
			curr->get_Lengths(curr_lengths);
			for(int j = 0; j< curr_lengths.size();j++) {
				myfile1 << curr_lengths.at(j) << endl;
			}

		}
		myfile1.close();
	}
	else {
		cout << "unable to open file" << endl;
	}
	
	if(myfile2.is_open()) {
		for(int i = 0; i < cells.size();i++) {
			curr = cells.at(i);
			curr->get_Sides(sides);
			curr_side = sides.at(3);
			curr_side->get_Lengths(curr_lengths); 
			curr_side->get_Forces(curr_forces);
		//	myfile2<< "First the Lengths" << endl;
			for(int j = 0; j < curr_lengths.size();j++) {
				myfile2 << curr_lengths.at(j) << endl;
			}
		//	myfile2 << "Now the Forces" << endl;
			for(int j = 0; j< curr_forces.size();j++) {
				myfile2 << curr_forces.at(j) << endl;
			}
		}
		myfile2.close();
	}
	else {
		cout << "unable to open file" << endl;
	}
	return;
}
	
void Tissue::print_Data_Output(ofstream & ofs) {
	return;
}

int Tissue::update_VTK_Indices() {
	
	int id = 0;
	int rel_cnt = 0;

	//iterate through cells to reassign vtk id's - starting at 0
	for (unsigned int i = 0; i < cells.size(); i++) {
		//iterate through 
		rel_cnt += cells.at(i)->update_VTK_Indices(id);
	}

	return rel_cnt;
}

void Tissue::print_VTK_File(ofstream& ofs) {

	int rel_cnt = update_VTK_Indices();
		
	ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Point representing Sub_cellular elem model" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;
	// Good up to here

	//Need total number of points for all cells
	int num_Points = 0;
	for (unsigned int i = 0; i < cells.size(); i++) {
		num_Points += cells.at(i)->get_Node_Count();
	}

	ofs << "POINTS " << num_Points << " float" << endl;
	
	vector<int> start_points;
	vector<int> end_points;
	int count = 0;
	for (unsigned int i = 0; i < cells.size(); i++) {
		start_points.push_back(count);
		cells.at(i)->print_VTK_Points(ofs, count);
		end_points.push_back(count - 1);
	}
	ofs << endl;

	ofs << "CELLS " << cells.size() << ' ' << num_Points + start_points.size() << endl;
	// output list of indices for each our our plant cells
	for (unsigned int i = 0; i < cells.size(); i++) {
		ofs << cells.at(i)->get_Node_Count();

		for (int k = start_points.at(i); k <= end_points.at(i); k++) {
			ofs << ' ' << k;
		}
		ofs << endl;
	}
	// output pairs of node indices to draw adh line
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Adh(ofs);
	}
	ofs << endl;

	ofs << "CELL_TYPES " << start_points.size() << endl;
	for (unsigned int i = 0; i < start_points.size(); i++) {
		// type for entire cell relationship
		ofs << 2 << endl;
	}
	for (int i = 0; i < rel_cnt; i++) {
		// type for adh spring relationship
		ofs << 3 << endl;
	}

	ofs << endl;

	
	ofs << "POINT_DATA " << num_Points << endl;
	ofs << "SCALARS magnitude float " << 1 << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars(ofs);
	}

	ofs << endl;

	ofs << "VECTORS force float" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Vectors(ofs);
	}

	

	return;
}
/*
void Tissue::grow_Cells(const int Ti) {
	
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->add_Wall_Node(Ti);
		cells.at(i)->add_Cyt_Node(Ti);
	}

	return;
}
*/

//=========================
//End of tissue.cpp
