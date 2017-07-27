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
	int numSteps = 4;

	// Variable for dataoutput
	int digits;
	string format1 = ".vtk";
	string format2 = ".txt";
	string Number;
	string initial_1 = "Animation/Plant_Cell_";
	string initial_2 = "DataOutput/Plant_Cell_";
	string Filename;
	ofstream ofs;

	//loop for time steps
	for(int Ti = 0; Ti < numSteps; Ti++) {
		//loop through all cells
		//for now only one cell
		cout << "Ti: " << Ti << endl;
		//Print to dataOutput and VTK files

		digits = ceil(log10(Ti + 1));
		if (digits == 1 || digits == 0) {
			Number = "0000" + to_string(Ti);
		}
		else if (digits == 2) {
			Number = "000" + to_string(Ti);
		}
		else if (digits == 3) {
			Number = "00" + to_string(Ti);
		}
		else if (digits == 4) {
			Number = "0" + to_string(Ti);
		}

		Filename = initial_1 + Number + format1;
		ofs.open(Filename.c_str());
		growing_Cell->print_VTK_File(ofs);
		ofs.close();

		Filename = initial_2 + Number + format2;
		ofs.open(Filename.c_str());
		growing_Cell->print_Data_Output(ofs);
		ofs.close();
		
		growing_Cell->calc_New_Forces();
		growing_Cell->update_Node_Locations();
		growing_Cell->update_Wall_Angles();

	}
	//=================
	/*
	// Five locations for Wall Nodes
	Coord locA(0.1, 0.1); //first_corner
	Coord locB(0.2, 0.1); //end node to left of A
	Coord locC(0.3, 0.1); //end node to left of B
	Coord locD(0.1, 0.2); //flank node to right of A
	Coord locE(0.1, 0.3); //flank node to right of D

	Wall_Node* A = new Corner_Node(locA);
	Wall_Node* B = new End_Node(locB);
	Wall_Node* C = new End_Node(locC);
	Wall_Node* D = new Flank_Node(locD);
	Wall_Node* E = new Flank_Node(locE);

	A->set_Left_Neighbor(B);
	B->set_Right_Neighbor(A);
	B->set_Left_Neighbor(C);
	C->set_Right_Neighbor(B);
	
	A->set_Right_Neighbor(D);
	D->set_Left_Neighbor(A);
	D->set_Right_Neighbor(E);
	E->set_Left_Neighbor(D);

	A->update_Angle();
	B->update_Angle();
	D->update_Angle();

	// Two locations for Cyt Nodes
	Coord locF(0.15, 0.2);
	Coord locG(0.2, 0.15);

	Cyt_Node* F = new Cyt_Node(locF);
	Cyt_Node* G = new Cyt_Node(locG);
	
	vector<Cyt_Node*> cyts;
	cyts.push_back(F);
	cyts.push_back(G);


	// Perform Calculations

	//Bending
	cout << "A's Bending Force" << endl;
	Coord bend = A->calc_Bending();
	cout << "	Total Bending Force: " << bend << endl << endl;

	cout << "A's Morse Force" << endl;
	Coord morse = A->calc_Morse_SC(cyts);
	cout << "	Total Morse Force: " << morse << endl << endl;
	
	cout << "A's Linear Force" << endl;
	Coord lin = A->calc_Linear();
	cout << "	Total Spring Force: " << lin << endl << endl;

	Coord total = bend + morse + lin;
	cout << "A's Total Force: " << total << endl;


	*/
	
	return 0;
}



