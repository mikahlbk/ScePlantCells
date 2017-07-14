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
		Location my_loc;
		Force new_Force;
    public:
    //functions that you will want performed on all nodes
        //Constructor
        Node(Location loc);
        //some functions you can define in base class because 
        //    all nodes will use the exact same function
        virtual Location get_location();
		virtual morse_Equation(Node* node, double U, double W,
							double Z, double G);
        //other functions might be executed differently based on
        //    which node you are. Thus define as "pure virtual" and 
        //    properly define them in a derived class
        virtual Force calc_Forces() = 0;
};

class Cyt_Node: public Node {
    private: 
    //if don't need to keep any more information then just leave blank

    public:
        Cyt_Node(Location loc);
        virtual Force calc_Forces(Cell* my_cell);
		virtual Force calc_Morse_II(vector<Cyt_Node*>& cyt_nodes);
		virtual Force calc_Morse_MI(Wall_Node* curr);
};

class Wall_Node: public Node {
    protected:
    //variables that will be shared by all wall nodes
        Wall_Node* left;
        Wall_Node* right;
        double my_angle;

    public:
    //function that you want performed on all wall nodes
        Wall_Node(Location loc);
        Wall_Node(Location loc, Wall_Node* left, Wall_Node* right, double angle);
        //maybe could define them here if corner and edge both perform
        //    these functions identically
		virtual double get_Angle();
        virtual Force calc_Forces();
		virtual Force calc_Morse_SC();
		virtual Force calc_Morse_DC();
		virtual Force linear_Equation();
		virtual Force bending_Equation();

        //otherwise set as pure virtual
		virtual Force calc_Linear() = 0;
		virtual Force calc_Bending() = 0;
		virtual double get_Equi_Angle() = 0;
		virtual double get_linearSpring() = 0;
};

class Corner_Node: public Wall_Node {
    public:
        Corner_Node(Location loc);
        Corner_Node(Location loc, Node* left, Node* right, double angle);
		virtual double get_Equi_Angle();
		virtual Force calc_Linear();
		virtual Force calc_Bending();
};

class Flank_Node: public Wall_Node {
    public:
        Flank_Node(Location loc);
        Flank_Node(Location loc, Node* left, Node* right, double angle);
		virtual Force calc_Linear();
		virtual Force calc_Bending();
		virtual double get_Equi_Angle();
		virtual double get_linearSpring();
};

class End_Node: public Wall_Node {
	public:
		End_Node(Location loc);
		End_Node(Location loc, Node* left, Node* right, double angle);
		virtual Force calc_Linear();
		virtual Force calc_Bending();
		virtual double get_Equi_Angle();
		virtual double get_linearSpring()
};

//===========================
#endif  
