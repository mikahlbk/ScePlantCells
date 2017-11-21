// phys.h

//===================
// Include Guard
#ifndef _PHYS_H_INCLUDED_
#define _PHYS_H_INCLUDED_
//===================
// Forward declarations

//====================
// Include dependencies
#include <math.h>
#include "coord.h"
//=====================

// Simulation Constants
const double dt = .001;
const int Relaxation_Time = 3000;
const int ADD_WALL_TIMER = 100;
const int ADD_CYT_TIMER = 500;

const int MAX_NUM_WALL = 200;
const int MAX_NUM_CYT = 60;

const int Init_Num_Cyt_Nodes = 20;
const int Init_Wall_Nodes = 100;
const double AREA_DOUBLED = 4;
// Global Physics Constants
const double pi = acos(-1.0);
const double window = .2;
///// Cell wall mechanical parameters

const double K_BEND = .001;
const Coord K_LINEAR_LONG = Coord(450, 100);
const Coord K_LINEAR_WIDE = Coord(100, 450);
const double K_ADH = 20;
const double MembrEquLen_ADH = 0.6;

const double EXTERNAL_FORCE = 2;
//linear spring equilibrium length
const double MembrEquLen = .0625; 
const double MEMBR_THRESH_LENGTH = 0.1; 
const double ADHThreshBeg = .78125;
const double ADHThresh = .78125;
const double damping = 1;
///// Subcellular element parameters for membrane - membrane interactions

const double U_MM =  3.90625;
const double W_MM =  0;
const double xsi_MM = 0.125;
const double gamma_MM = 1.5625;	

///// Subcellular element parameters for membrane  - internal interactions

const double U_MI = 0.78125;  
const double W_MI = 0.00;
const double xsi_MI = 0.125;
const double gamma_MI = 0.625;

///// Subcellular element parameters for internal - internal interactions

const double U_II = 0.488;
const double W_II = 0.146484;
const double xsi_II = 0.3125;
const double gamma_II = 1.25;


//=====================
#endif
