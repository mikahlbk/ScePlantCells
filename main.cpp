#include <iostream>

using namespace std;

int main() {
	//make a new cell object
	Cell* growing_Cell = new Cell();

	//parameters for time step
	double dt = .0001;
	double numSteps = 100;

	//loop for time steps
	for(i=0;i=numSteps;i+=dt) {
		//loop through all cells
		//for now only one cell
		//growth
		if (i % 200 == 0) {
			growing_Cell->add_Cell_Wall_Node();
		}
		if (i % 500 ==0) {
			growing_Cell->add_Cyt_Node();
		}
		
	growing_Cell->calc_New_Forces();
	growing_Cell->update_Node_Locations();
	growing_Cell->update_Wall_Angels();
		
	return 0;
}
