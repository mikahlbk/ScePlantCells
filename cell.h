//cell.h
//===================
// Inlcude Guards
#ifndef _CELL_H_INCLUDED_
#define _CELL_H_INCLUDED_
//===================
// forward declarations
class Tissue;
//===================
// include dependencies
#include <string>
#include <vector>
#include <fstream>
#include <iostream>

#include "phys.h"
#include "coord.h"
#include "node.h"
//===================
// Cell Class Declaration

class Cell {

	private:
		int rank;
		Tissue* my_tissue;
		int life_length;
		int num_cyt_nodes;
		int layer;
		int growth_rate;
		Coord cell_center;
		int num_wall_nodes;
		vector<double> strain_vec;
		vector<double> stress_vec;
		vector<Cyt_Node*> cyt_nodes;
		vector<Cell*> neigh_cells;
		Wall_Node* most_up;
		Wall_Node* most_down;
		Wall_Node* most_right;
		Wall_Node* most_left;
		Wall_Node* left_Corner;	
	public:
		
		Cell(Tissue* tissue);
		Cell(int rank, Coord center, double radius, Tissue* tiss, int layer);

		// Destructor
		~Cell();

		// Getters and Setters
		int get_Rank() {return rank;}
		int get_Life_Length() {return life_length;}
		int get_Cytoplasm_Count() {return num_cyt_nodes;}
		int get_Layer() {return layer;}
		int get_Growth_Rate() {return growth_rate;}
		Coord get_Cell_Center() {return cell_center;}
		int get_Wall_Count() {return num_wall_nodes;}
		int get_Node_Count();
		void get_Cyt_Nodes(vector<Cyt_Node*>& cyts);
		Tissue* get_Tissue() {return my_tissue;}
		Wall_Node* get_Wall_Nodes() {return left_Corner;}
		void set_Rank(const int id);
		void set_Layer(int layer);
		void set_growth_rate(double rate);
		void get_Neighbor_Cells(vector<Cell*>& cells);
		void get_Strain(vector<double>& strain);
		void get_Stress(vector<double>& stress);
		Wall_Node* get_Left_Corner() {return left_Corner;}
		void set_Left_Corner(Wall_Node*& new_left_corner);
		void set_Wall_Count(int& number_nodes);
		//Wall_Node* get_most_right() {return most_right;}
		//Wall_Node* get_most_left() {return most_left;}
		//Wall_Node* get_most_up() {return most_up;}
		//Wall_Node* get_most_down() {return most_down;}
		// Keep track of neighbor cells
		void update_Neighbor_Cells();
		void update_adhesion_springs(int Ti);
	
		// Forces and Positionsing
		void calc_New_Forces();
		void update_Node_Locations();
		void update_Wall_Angles();
		void update_Wall_Equi_Angles();
		void update_Cell_Center();
		void update_Life_Length();
		void wall_Node_Check();
		void cytoplasm_Check();
		void stretch();
		//Output Functions
		void print_Data_Output(ofstream& ofs);
		int update_VTK_Indices(int& id);
		void print_VTK_Adh(ofstream& ofs);
		void print_VTK_Points(ofstream& ofs, int& count);
		void print_VTK_Scalars(ofstream& ofs);
		void print_VTK_Vectors(ofstream& ofs);
		
		// Growth of cell
		Wall_Node* find_Largest_Length();
		void add_Wall_Node();
		void add_Cyt_Node();
		double total_Force();
		void most_Left_Right(Wall_Node*& left, Wall_Node*& right);
		double tensile_Length();
		double extensional_Length();
		void most_Up_Down(Wall_Node*& up, Wall_Node*& down);
		void add_stress(double& new_length, double& new_force);
		void add_strain(double& new_length);
		double calc_Area();
		//Division 
		Cell* divide(int Ti);
		Cell* divide_length_wise();
		//Cell* divide_width_wise(const int Ti);
		void add_Cyt_Node_Div(double& radius_x, double& radius_y);
		
};


// End Cell Class
//===================

#endif

