// Main.cpp
//============================

//===========================
// Include Dependencies
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>

#include "phys.h"
#include "coord.h"
//#include "node.h"
#include "cell.h"
//==============================

using namespace std;

//============================

int main() {

	string init_cell = "cell_start.txt";
	//make a new cell object

	Cell* growing_Cell = new Cell(init_cell);
	cout << "Finished creating Cell" << endl;
	//parameters for time step
    double numSteps = 10;

	// Variable for dataoutput
	int digits;
	string format = ".vtk";
	string Number;
	string initial = "Animation/Plant_Cell_";
	string Filename;
	ofstream ofs;
	int out = 0; //counter for creating output/vtk files

	//loop for time steps
	for(int Ti = 0; Ti * dt < numSteps; Ti++) {
		//loop through all cells
		//for now only one cell

		//Print to dataOutput and VTK files

		if (Ti % 100 == 0) {
	
			digits = ceil(log10(out + 1));
			if (digits == 1 || digits == 0) {
				Number = "0000" + to_string(out);
			}
			else if (digits == 2) {
				Number = "000" + to_string(out);
			}
			else if (digits == 3) {
				Number = "00" + to_string(out);
			}
			else if (digits == 4) {
				Number = "0" + to_string(out);
			}

			Filename = initial + Number + format;

			ofs.open(Filename.c_str());
			growing_Cell->print_VTK_File(ofs);
			ofs.close();
		
			out++;
		}

		if (Ti % 1000 == 0) {
			cout << "Simulation still running. Ti: " << Ti << endl;
		}

		//growth

//		if (Ti % 1000 == 0) {
//			cout << "Adding cell wall node" << endl;
//			cout << "Ti : " << Ti << endl;
//			growing_Cell->add_Cell_Wall_Node();
//			cout << "Completed adding cell wall node" << endl;
//		}
		if (Ti % 5000 == 0) {
//			cout << "Adding cyt node" << endl;
			growing_Cell->add_Cyt_Node();
//			cout << "Completed adding cyt node" << endl;
		}

		growing_Cell->calc_New_Forces();
		growing_Cell->update_Node_Locations();
		growing_Cell->update_Wall_Angles();
		
		

	}
		
	return 0;
}











