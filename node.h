//node.h
//=====================
// Include Guards
#ifndef NODE_H  //if node.h hasn't been included yet
#define NODE_H  //    define it so the compiler knows 
//=====================
// Forward Declarations 
class Wall_Node;
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
		Cell* my_cell;
    public:
    //functions that you will want performed on all nodes
        //Constructor
        Node(Coord loc, Cell* my_cell);
        //some functions you can define in base class because 
        //    all nodes will use the exact same function
        virtual Coord get_Location();
		virtual Coord get_Force();
		virtual void update_Location();
        //other functions might be executed differently based on
        //    which node you are. Thus define as "pure virtual" and 
        //    properly define them in a derived class
        virtual void calc_Forces() = 0;
};

class Cyt_Node: public Node {
    private: 
    //if don't need to keep any more information then just leave blank

    public:
        Cyt_Node(Coord loc, Cell* my_cell);
        virtual void calc_Forces();
		Coord calc_Morse_II(vector<Cyt_Node*>& cyt_nodes);
		Coord calc_Morse_MI(Wall_Node* curr);
		Coord morse_Equation(Cyt_Node* cyt);
		Coord morse_Equation(Wall_Node* wall);
};

class Wall_Node: public Node {
    protected:
    //variables that will be shared by all wall nodes
        Wall_Node* left;
        Wall_Node* right;
        double my_angle;
		double cross_Prod;
		Coord cyt_force;

    public:
    //function that you want performed on all wall nodes
		// Constructors
        Wall_Node(Coord loc, Cell* my_cell);
        Wall_Node(Coord loc, Cell* my_cell, Wall_Node* left, Wall_Node* right);
        //maybe could define them here if corner and edge both perform
        //    these functions identically

		// Getters and Setters
		virtual double get_Angle();
		virtual Coord get_CytForce();
        virtual Wall_Node* get_Left_Neighbor();
		virtual Wall_Node* get_Right_Neighbor();
		virtual void set_Left_Neighbor(Wall_Node* new_Left);
		virtual void set_Right_Neighbor(Wall_Node* new_Right);
		// Force Calculations
		virtual void calc_Forces();
		virtual void update_Angle();
		virtual Coord calc_Morse_SC();
		virtual Coord calc_Morse_DC();
		virtual Coord calc_Bending();
		virtual Coord morse_Equation(Cyt_Node* cyt);
		virtual Coord morse_Equation(Wall_Node* wall);
		virtual Coord linear_Equation(Wall_Node* wall, double k_Linear);
		virtual Coord bending_Equation_Center();
		virtual Coord bending_Equation_Left();
		virtual Coord bending_Equation_Right();

        //otherwise set as pure virtual
		virtual Coord calc_Linear() = 0;
		virtual double get_Equi_Angle() = 0;
		virtual double get_linearSpring() = 0;
		virtual double get_bendingSpring() = 0;
		virtual bool is_Corner() = 0;
};

class Corner_Node: public Wall_Node {
    public:
        Corner_Node(Coord loc, Cell* my_cell);
        Corner_Node(Coord loc, Cell* my_cell, Wall_Node* left, Wall_Node* right);
		virtual double get_Equi_Angle();
		virtual double get_linearSpring();
		virtual double get_bendingSpring();
		virtual Coord calc_Linear();
		virtual bool is_Corner();
};

class Flank_Node: public Wall_Node {
    public:
        Flank_Node(Coord loc, Cell* my_cell);
        Flank_Node(Coord loc, Cell* my_cell, Wall_Node* left, Wall_Node* right);
		virtual Coord calc_Linear();
		virtual double get_Equi_Angle();
		virtual double get_linearSpring();
		virtual double get_bendingSpring();
		virtual bool is_Corner();
};

class End_Node: public Wall_Node {
	public:
		End_Node(Coord loc, Cell* my_cell);
		End_Node(Coord loc, Cell* my_cell, Wall_Node* left, Wall_Node* right);
		virtual Coord calc_Linear();
		virtual double get_Equi_Angle();
		virtual double get_linearSpring();
		virtual double get_bendingSpring();
		virtual bool is_Corner();
};

//===========================
#endif  
