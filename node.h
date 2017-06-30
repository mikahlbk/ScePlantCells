//node.h

#ifndef NODE_H  //if node.h hasn't been included yet
#define NODE_H  //    define it so the compiler knows 

class Node {
    private:
    //variables that will be shared by all nodes
        double x_coord;
        double y_coord;
    public:
    //functions that you will want performed on all nodes
        //Constructor
        Node(double x, double y);
        //some functions you can define in base class because 
        //    all nodes will use the exact same function
        virtual void get_location(double& x, double& y);
        //other functions might be executed differently based on
        //    which node you are. Thus define as "pure virtual" and 
        //    properly define them in a derived class
        virtual void set_Location() = 0;
        virtual double calc_Forces() = 0;
};

class Wall_Node: public Node {
    private:
    //variables that will be shared by all wall nodes
        Node* left_neighbor;
        Node* right_neighbor;
        double my_angle;

    public:
    //function that you want performed on all wall nodes
        Wall_Node(double x, double y);
        Wall_Node(double x, double y, Node* left, Node* right);
        //maybe could define them here if corner and edge both perform
        //    these functions identically
        virtual void set_Location();
        virtual double calc_Forces();
        //otherwise set as pure virtual
        virtual bool is_Corner() = 0;

};

class Cyt_Node: public Node {
    private: 
    //if don't need to keep any more information then just leave blank

    public:
        Cyt_Node(double x, double y);
        virtual void set_Location();
        virtual double calc_Forces();

};

class Corner_Node: public Wall_Node {
    private:

    public:
        Corner_Node(double x, double y);
        Corner_Node(double x, double y, Node* left, Node* right);
        virtual bool is_Corner();

};

class Edge_Node: public Wall_Node {
    private:

    public:
        Edge_Node(double x, double y);
        Edge_Node(double x, double y, Node* left, Node* right);
        virtual bool is_Corner();

};


#endif  
