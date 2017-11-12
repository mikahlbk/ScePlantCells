// Main.cpp
//============================

//===========================
// Include Dependencies
#include <iostream>
#include <math.h>
#include <string>
#include <vector>
#include <fstream>
#include <ctime>

#include "phys.h"
#include "coord.h"
#include "node.h"
#include "cell.h"
#include "tissue.h"
//==============================

using namespace std;

//============================

int main() {

	int start = clock();
	
	string init_tissue = "cell_start.txt";
	
	//make new cell objects in tissue
	Tissue growing_Tissue(init_tissue);

	cout << "Finished creating Cells" << endl;
	
	//parameters for time step
    double numSteps = 500;

	// Variable for dataoutput
	int digits;
	string format = ".vtk";
	string Number;
	string initial = "Animation/Plant_Cell_";
	string Filename;
	ofstream ofs;
	int out = 0; //counter for creating output/vtk files

	//loop for time steps
	for(int Ti = 0; Ti*dt< numSteps; Ti++) {
		//loop through all cells
		//for now only one cell
		//cout << "Ti = " << Ti << endl;
		//Print to dataOutput and VTK files
		if (Ti % 100  == 0) {
			
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
			growing_Tissue.print_VTK_File(ofs);
			ofs.close();	
			out++;
		}
	
		if (Ti % 1000 == 0) {
			cout << "Simulation still running. Ti: " << Ti << endl;
		}	
		
		
		// Update Each cell's neighboring cells
		if (Ti % 10  == 0) {
			//cout << "Find Neighbors" << endl;
			growing_Tissue.update_Neighbor_Cells();
			//cout << "Make Adhesion" << endl;
			growing_Tissue.update_Adhesion();
    	}
		

		// Tissue Growth
		growing_Tissue.update_Life_Length();

		// Add cyt node/ wall node 
		if(Ti > 35000) {
			growing_Tissue.update_Cytoplasm();
		}

		growing_Tissue.update_Wall();
//		cout << "updated life length" << endl;
		//Calculate new forces on cells and nodes
		growing_Tissue.calc_New_Forces();
//		cout << "calculated forces" << endl;
		//Update node positions
		growing_Tissue.update_Cell_Locations();
//		cout << "updated node positions" << endl;
		//Division if necessary
//		growing_Tissue.cell_Division(Ti);
//		cout << "Division" << endl;	
	}

	int stop = clock();

	cout << "Time: " << (stop - start) / double(CLOCKS_PER_SEC) * 1000 << endl;
		
	return 0;
}











