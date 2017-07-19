//node.h
//=====================
// Include Guards
#ifndef NODE_H  //if node.h hasn't been included yet
#define NODE_H  //    define it so the compiler knows 
//=====================
// Forward Declarations 

//=====================
// Include Declarations
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
		virtual void update_Location();
        //other functions might be executed differently based on
        //    which node you are. Thus define as "pure virtual" and 
        //    properly define them in a derived class
        virtual void calc_Forces(Cell* my_cell) = 0;
};

class Cyt_Node: public Node {
    private: 
    //if don't need to keep any more information then just leave blank

    public:
        Cyt_Node(Coord loc);
        virtual void calc_Forces(Cell* my_cell);
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

    public:
    //function that you want performed on all wall nodes
        Wall_Node(Coord loc, double angle);
        Wall_Node(Coord loc, Wall_Node* left, Wall_Node* right, double angle);
        //maybe could define them here if corner and edge both perform
        //    these functions identically
		virtual double get_Angle();
        virtual Wall_Node* get_Left_Neighbor();
		virtual Wall_Node* get_Right_Neighbor();
		virtual void set_Left_Neighbor();
		virtual void set_Right_Neighbor();
		virtual Coord calc_Forces(Cell* my_cell);
		virtual Coord calc_Morse_SC(vector<Cyt_Node*>& cyt_nodes);
		virtual Coord calc_Morse_DC(vector<Cell*>& cells);
		virtual Coord morse_Equation(Cyt_Node* cyt);
		virtual Coord morse_Equation(Wall_Node* wall);
		virtual Coord linear_Equation(Wall_Node* wall, k_Linear);
		virtual Coord bending_Equation_Center();
		virtual Coord bending_Equation_Left();
		virtual Coord bending_Equation_Right();
		virtual Coord calc_Bending();

        //otherwise set as pure virtual
		virtual Coord calc_Linear() = 0;
		virtual double get_Equi_Angle() = 0;
		virtual double get_linearSpring() = 0;
		virtual double get_bendingSpring() = 0;
};

class Corner_Node: public Wall_Node {
    public:
        Corner_Node(Coord loc);
        Corner_Node(Coord loc, Node* left, Node* right, double angle);
		virtual double get_Equi_Angle();
		virtual Coord calc_Linear();
};

class Flank_Node: public Wall_Node {
    public:
        Flank_Node(Coord loc);
        Flank_Node(Coord loc, Node* left, Node* right, double angle);
		virtual Coord calc_Linear();
		virtual double get_Equi_Angle();
		virtual double get_linearSpring();
};

class End_Node: public Wall_Node {
	public:
		End_Node(Coord loc);
		End_Node(Coord loc, Node* left, Node* right, double angle);
		virtual Coord calc_Linear();
		virtual Coord calc_Bending();
		virtual double get_Equi_Angle();
		virtual double get_linearSpring()
};

//===========================
#endif  
