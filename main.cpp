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

int main(int argc, char* argv[]) {

	string anim_folder = argv[1];

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
	string initial = "/Plant_Cell_";
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

			Filename = anim_folder +  initial + Number + format;

			ofs.open(Filename.c_str());
			growing_Tissue.print_VTK_File(ofs);
			ofs.close();	
			out++;
		}
		
		if (Ti % 100 == 0) {
			cout << "Simulation still running. Ti: " << Ti << endl;
		}
		
	
		// Tissue Growth
		//cout << "update life length" << endl;

		growing_Tissue.update_Life_Length();
		//cout << "life length updated" << endl;	
		//Add cyt node/ wall node 
		if(Ti > 500) {
		//	cout << "update cytoplasm" << endl;
			growing_Tissue.update_Cytoplasm();
		}

		if(Ti > 500) {
		//	cout << "Updated Wall" << endl;
			growing_Tissue.update_Wall();
		}
		
		if((Ti > 200)) {
			if (Ti% 100  == 0 ) {
			//	cout << "Find Neighbors" << endl;
				growing_Tissue.update_Neighbor_Cells();
		//	}
		//	if(Ti%10 == 0) {
			//	cout << "Make Adhesion" << endl;
				growing_Tissue.update_Adhesion(Ti);
			}
		}
		//Division if necessary
		if(Ti > 1500) {
		//	if(Ti% 10 == 0) {
			cout << "Check if cells need to divide" << endl;
				growing_Tissue.cell_Division(Ti);
			//}
		}
		
		//cout << "division success" << endl;
		/*if((Ti > 268583)) {
			growing_Tissue.stretching_Test();
			growing_Tissue.cell_stress();
		}
		if((Ti > 268582)) {
			growing_Tissue.cell_strain();
		}*/
		// Update Each cell's neighboring cells
		
		//cout << "do forces" << endl;
		//Calculate new forces on cells and nodes
		growing_Tissue.calc_New_Forces();
		//cout << "calculated forces" << endl;
		//Update node positions
		growing_Tissue.update_Cell_Locations();
		//updates cell centers
		//cout << "updated node positions" << endl;
		
	}

	//for calibration of elastic modulus
	//growing_Tissue.make_Vectors();
	
	int stop = clock();

	cout << "Time: " << (stop - start) / double(CLOCKS_PER_SEC) * 1000 << endl;
		
	return 0;
}











