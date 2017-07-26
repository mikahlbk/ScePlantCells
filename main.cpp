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

//	cout << "A's Morse Force" << endl;
//	Coord morse = A->calc_Morse_SC(cyts);
//	cout << "	Total Morse Force: " << morse << endl << endl;
	
	cout << "A's Linear Force" << endl;
	Coord lin = A->calc_Linear();
	cout << "	Total Spring Force: " << lin << endl << endl;

	cout << "B's Linear Force" << endl;
	Coord linB = B->calc_Linear();
	cout << " Total Spring Force: " << linB << endl << endl;

	Coord total = bend +lin; //morse 
	cout << "A's Total Force: " << total << endl;



	
	return 0;
}



