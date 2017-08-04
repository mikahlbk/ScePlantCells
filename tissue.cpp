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
//=========================
// Public Member Functions for Tissue.cpp

Tissue::Tissue(string filename) {

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
		else if (temp == "End_Cell") {
			//create new cell with collected data and push onto vector 
			curr = new Cell(rank, corner, height, width, 0, this);
			cells.push_back(curr);
		}

		ss.clear();
	}

	ifs.close();
}

void Tissue::get_Cells(vector<Cell*>& cells) {
	cells = this->cells;
	return;
}

void Tissue::calc_New_Forces(const int Ti) {

	string init = "DataOutput/wall_forces_";
	string end = ".txt";
	string num = to_string(Ti);

	string filename = init + num + end;

	ofstream ofs(filename.c_str());

	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->calc_New_Forces(ofs);
	}

	ofs.close();

	return;
}

void Tissue::update_Cell_Locations() {
	
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Node_Locations();
	}

	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->update_Wall_Angles();
	}

	return;
}

void Tissue::print_Data_Output(ofstream & ofs) {
	return;
}

void Tissue::print_VTK_File(ofstream& ofs) {
		
	ofs << "# vtk DataFile Version 3.0" << endl;
	ofs << "Point representing Sub_cellular elem model" << endl;
	ofs << "ASCII" << endl << endl;
	ofs << "DATASET UNSTRUCTURED_GRID" << endl;
	// Good up to here

	//Need total number of points for all cells
	int num_Points = 0;
	for (unsigned int i = 0; i < cells.size(); i++) {
		num_Points += cells.at(i)->get_Num_Nodes();
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

	ofs << "CELLS " << cells.size() << ' ' << num_Points + start_points.size() << endl;

	for (unsigned int i = 0; i < cells.size(); i++) {
		ofs << cells.at(i)->get_Num_Nodes();

		for (int k = start_points.at(i); k <= end_points.at(i); k++) {
			ofs << ' ' << k;
		}
		ofs << endl;
	}

	ofs << endl;

	ofs << "CELL_TYPES " << start_points.size() << endl;
	for (unsigned int i = 0; i < start_points.size(); i++) {
		ofs << 2 << endl;
	}

	return;
}

void Tissue::grow_Cells(const int Ti) {
	
	for (unsigned int i = 0; i < cells.size(); i++) {
		cells.at(i)->add_Wall_Node(Ti);
		cells.at(i)->add_Cyt_Node(Ti);
	}

	return;
}

//=========================
//End of tissue.cpp
