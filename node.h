//node.h
//=====================
// Include Guards
#ifndef _NODE_H_INCLUDED_  //if node.h hasn't been included yet
#define _NODE_H_INCLUDED_  //    define it so the compiler knows 
//=====================
// Forward Declarations 
class Wall_Node;
class Side;
class Cell;
//=====================
// Include Declarations
#include <iostream>
#include <vector>
#include "phys.h"
#include "coord.h"
//=====================

class Node {
    protected:
    //variables that will be shared by all nodes
		Coord my_loc;
		Coord new_force;
    public:
    //functions that you will want performed on all nodes
        //Constructor
        Node(Coord loc);
        //some functions you can define in base class because 
        //    all nodes will use the exact same function
        virtual Coord get_Location();
		virtual Coord get_Force();
		virtual void update_Location();
        //other functions might be executed differently based on
        //    which node you are. Thus define as "pure virtual" and 
        //    properly define them in a derived class
        virtual void calc_Forces() = 0;

		virtual ~Node();
};

class Cyt_Node: public Node {
    private: 
    //if don't need to keep any more information then just leave blank
		Cell* my_cell;
    public:
        Cyt_Node(Coord loc, Cell* my_cell);
        virtual void calc_Forces();
		Coord calc_Morse_II();
		Coord calc_Morse_MI(Wall_Node* orig);
		Coord morse_Equation(Cyt_Node* cyt);
		Coord morse_Equation(Wall_Node* wall);
		~Cyt_Node();
};

class Wall_Node: public Node {
    protected:
    //variables that will be shared by all wall nodes
        Wall_Node* left;
        Wall_Node* right;
        double my_angle;
		double equi_angle;
		double cross_Prod;
		Coord cyt_force;
		bool on_curve;
		Side* my_side;
    public:
    //function that you want performed on all wall nodes
		// Constructors
        Wall_Node(Coord loc, Side* my_side);
        Wall_Node(Coord loc, Side* my_side, Wall_Node* left, Wall_Node* right);
        //maybe could define them here if corner and edge both perform
        //    these functions identically

		// Getters and Setters
		double get_Angle() {return my_angle;}
		double get_Equi_Angle() {return equi_angle;}
		void set_Equi_Angle(double angle);
		bool is_Curve() {return on_curve;}
		void set_Curve(bool on_curve);
		Coord get_CytForce() {return cyt_force;}
        Wall_Node* get_Left_Neighbor() {return left;}
		Wall_Node* get_Right_Neighbor() {return right;}
		double get_Linear_Spring();
		double get_Bending_Spring();
		void set_Left_Neighbor(Wall_Node* new_Left);
		void set_Right_Neighbor(Wall_Node* new_Right);
		Side* get_My_Side() {return my_side;}

		// Force Calculations
		virtual void calc_Forces();
		Coord calc_Morse_SC();
		Coord calc_Morse_DC();
		Coord calc_Bending();
		Coord calc_Linear();
		void update_Angle();

		// Mathematical Force Equations
		Coord morse_Equation(Cyt_Node* cyt);
		Coord morse_Equation(Wall_Node* wall);
		Coord bending_Equation_Center();
		Coord bending_Equation_Left();
		Coord bending_Equation_Right();
		Coord linear_Equation(Wall_Node* wall);

		~Wall_Node();
};


//===========================
#endif  
