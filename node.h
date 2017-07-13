//node.h
//=====================
// Include Guards
#ifndef NODE_H  //if node.h hasn't been included yet
#define NODE_H  //    define it so the compiler knows 
//=====================
// Forward Declarations 

//=====================
// Include Declarations
#include "loc.h"
//=====================

class Node {
    private:
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
        //other functions might be executed differently based on
        //    which node you are. Thus define as "pure virtual" and 
        //    properly define them in a derived class
        virtual Force calc_Forces() = 0;
		virtual Force calc_Morse() = 0;
		virtual Force calc_Linear() = 0;
		virtual Force calc_Bending() = 0;
        virtual void set_New_Location() = 0;
};

class Wall_Node: public Node {
    private:
    //variables that will be shared by all wall nodes
        Node* left_neighbor;
        Node* right_neighbor;
        double my_angle;

    public:
    //function that you want performed on all wall nodes
        Wall_Node(Location loc);
        Wall_Node(Location loc, Node* left, Node* right, double angle);
        //maybe could define them here if corner and edge both perform
        //    these functions identically
		virtual double get_Angle();
        virtual Force calc_Forces();
		virtual Force calc_Morse();
		virtual Force calc_Linear();
		virtual Force calc_Bending();
        virtual void set_Location();
        //otherwise set as pure virtual

};

class Cyt_Node: public Node {
    private: 
    //if don't need to keep any more information then just leave blank

    public:
        Cyt_Node(Location loc);

        virtual Force calc_Forces();
		virtual Force calc_Morse();
        virtual void set_Location();

};

class Corner_Node: public Wall_Node {
    private:

    public:
        Corner_Node(Location loc);
        Corner_Node(Location loc, Node* left, Node* right, double angle);
		virtual Force calc_Forces();
		virtual Force 
};

class Flank_Node: public Wall_Node {
    private:

    public:
        Flank_Node(Location loc);
        Flank_Node(Location loc, Node* left, Node* right, double angle);

};

class End_Node: public Wall_Node {
	private:

	public:
		End_Node(Location loc);
		End_Node(Location loc, Node* left, Node* right, double angle);
		
};

//===========================
#endif  
