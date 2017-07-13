//node.cpp

#include "node.h"

/** class Node Functions **/
Node::Node(double x, double y) {
    this->x_coord = x;
    this->y_coord = y;
}
void Node::get_location(double& x, double& y) {
    x = this->x_coord;
    y = this->y_coord;
    return;
}
/** class Wall Node Functions **/
Wall_Node::Wall_Node(double x, double y) : Node(x,y) {};

Wall_Node::Wall_Node(double x, double y, Node* left, Node* right) : Node(x,y) {
    this->left_neighbor = left;
    this->right_neighbor = right;
}



Force Wall_Node::get_angle() {

}

Force Wall_Node::calc_Forces() {
    double number;

    return number;
}

Force Wall_Node::calc_Morse() {

}

Force Wall_Node::calc_Linear() {

}

Force Wall_Node::calc_Bending() {

}

/** class Cyt Node Functions **/
Cyt_Node::Cyt_Node(Location loc) : Node(loc) {};

Force Cyt_Node::calc_Forces() {
    double number;

    return number;
}

Force Cyt_node::calc_Morse() {


}

/** class Corner Node Functions **/
Corner_Node::Corner_Node(double x, double y) : Wall_Node(x,y) {};

Corner_Node::Corner_Node(double x, double y, Node* left, Node* right) 
    : Wall_Node(x,y,left,right) {}

bool Corner_Node::is_Corner() {
    return true;
}

/** class Edge Node function **/
Edge_Node::Edge_Node(double x, double y) : Wall_Node(x,y) {};

Edge_Node::Edge_Node(double x, double y, Node* left, Node* right)
    : Wall_Node(x,y,left,right) {}

bool Edge_Node::is_Corner() {
    return false;
}




