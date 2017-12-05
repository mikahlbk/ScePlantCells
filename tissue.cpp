//tissue.cpp
//=========================

//=========================
//Include Dependencies
#include <vector>
#include <fstream>
#include <sstream>
#include <string>
#include <iomanip>

#include "phys.h"
#include "coord.h"
#include "cell.h"
#include "tissue.h"
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
	double radius;
	Coord center;
	double x, y;
	Cell* curr;

	while (getline(ifs,line)) {
		ss.str(line);

		getline(ss,temp,':');

		if (temp == "CellRank") {
			ss >> rank;
		}
		else if (temp == "Center") {
			ss >> x >> trash >> y;
			Coord loc(x,y);
			center = loc;
		}
		else if (temp == "Radius") {
			ss >> radius;
		}
		else if (temp == "Layer") {
			ss >> layer;
		}
		else if (temp == "End_Cell") {
			//create new cell with collected data and push onto vector 
			curr = new Cell(rank, center, radius, this, layer);
			num_cells++;
			cells.push_back(curr);
		}

		ss.clear();
	}

	ifs.close();
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

void Tissue::update_Cell_Cycle(int Ti) {
	for (unsigned int i = 0; i < cells.size(); i++) {
		//cout << "updating cell i" << endl;
		cells.at(i)->update_Cell_Progress(Ti);
	}
	return;
}

/*void Tissue::update_Cytoplasm() {
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->cytoplasm_Check();
	}
	return;
}*/

void Tissue::update_Wall() {
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->wall_Node_Check();
		//cout<< "Wall Count Cell " << i << ": " << cells.at(i)->get_Wall_Count() << endl;
	}



	return;
}
void Tissue::update_Num_Cells(Cell*& new_Cell) {
	num_cells++;
	cells.push_back(new_Cell);
	return;
}

void Tissue::calc_New_Forces() {
	Cell* curr = NULL;
	for (unsigned int i = 0; i < cells.size(); i++) {
		curr = cells.at(i);
	//	cout << "calc forces for cell: " << i << endl;
		curr->calc_New_Forces();
	}
	return;


}

void Tissue::update_Cell_Locations() {
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Node_Locations();
		//updates cell centers
	}
	return;
}

void Tissue::update_Neighbor_Cells() {
	//update vectors of neighboring cells
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Neighbor_Cells();
	}
	return;
}

void Tissue::stretching_Test() {
	//stretch the cell
	for(unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->stretch();
	}
	return;
}

void Tissue::cell_strain() {
	double curr_length;
	Cell* curr_cell;
	for(unsigned int i = 0; i < cells.size(); i++) {
		curr_cell = cells.at(i);
		curr_length = curr_cell->extensional_Length();
		curr_cell->add_strain(curr_length);
	}
	return;
}

void Tissue::cell_stress() {
	double curr_length;
	double curr_force;
	Cell* curr_cell;
	for(unsigned int i = 0; i < cells.size(); i ++) {
		curr_cell = cells.at(i);
		curr_length = curr_cell->tensile_Length();
		curr_force = curr_cell->total_Force();
		curr_cell->add_stress(curr_length, curr_force);
	}
	return;
}

void Tissue::update_Adhesion(int Ti) {
	int time = Ti;
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

void Tissue::make_Vectors() {
	Cell* curr = NULL;
	vector<double> strain;
	vector<double> stress;
	ofstream myfile1("strain_vec.txt");
	ofstream myfile2("stress_vec.txt");
	if(myfile1.is_open()){
		for(int i = 0; i < cells.size();i++) {
			curr = cells.at(i);
			curr->get_Strain(strain);
			for(int j= 0; j < strain.size(); j++) {
					myfile1 << strain.at(j) << endl;
			}
		}
		myfile1.close();
	}
	else {
		cout << "unable to open file" << endl;
	}
	if(myfile2.is_open()){
		for(int i = 0; i < cells.size(); i++) {
			curr = cells.at(i);
			curr->get_Stress(stress);
			for(int j = 0;j < stress.size(); j++) {
				myfile2 << stress.at(j) << endl;
			}
		}
		myfile2.close();
	}
	else {
		cout << "unable to open file" << endl;
	}
	return;
}

void Tissue::print_Data_Output(ofstream& ofs) {
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

//	cout << "final ID: " << id << endl;
//	cout << "rel_cnt: " << rel_cnt << endl;

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
	ofs << "CELLS " << cells.size()+rel_cnt << ' ' << (num_Points + start_points.size()) + (rel_cnt*3) << endl;

	for (unsigned int i = 0; i < cells.size(); i++) {
		ofs << cells.at(i)->get_Node_Count();

		for (int k = start_points.at(i); k <= end_points.at(i); k++) {
			ofs << ' ' << k;
		}
		ofs << endl;
	}
	
//	//output pairs of node indices to draw adh line
	for(unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Adh(ofs);
	}

	ofs << endl;

	ofs << "CELL_TYPES " << start_points.size() + rel_cnt << endl;
	for (unsigned int i = 0; i < start_points.size(); i++) {
		ofs << 2 << endl;
	}

	for(unsigned int i = 0; i < rel_cnt; i++) {
//		//type for adh relationship
		ofs << 3 << endl;
	}

	ofs << endl;


	ofs << "POINT_DATA " << num_Points << endl;
	ofs << "SCALARS magnitude double " << 1 << endl;
	ofs << "LOOKUP_TABLE default" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Scalars_WUS(ofs);
	}

	ofs << endl;

	ofs << "VECTORS force float" << endl;
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->print_VTK_Vectors(ofs);
	}

	

	return;
}


//=========================
//End of tissue.cpp
