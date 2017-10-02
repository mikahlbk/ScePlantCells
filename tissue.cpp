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

void Tissue::cell_Division(const int Ti) {
	bool divided = false;
	Cell* new_cell = NULL;
//	cout << "number cells: " << cells.size()<< endl;
	for(unsigned int i = 0; i < 1;i++) {
		//cout << "current divide cell: " << i << endl;
		new_cell = cells.at(i)->divide(Ti);
		if(new_cell !=NULL) {
			divided = true;
			new_cell->set_Rank(num_cells);
			num_cells++;
			cells.push_back(new_cell);
		}

	}
	if(divided) {
		this->update_Neighbor_Cells();
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

	for (unsigned int i = 0; i < cells.size(); i++) {
		ofs << cells.at(i)->get_Node_Count();

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

	ofs << endl;

	/*
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

	*/

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
