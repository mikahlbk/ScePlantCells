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
//=====================

// Simulation Constants
const double dt = 0.00003;
const int ADD_WALL_TIMER = 1500;
const int ADD_CYT_TIMER = 1000;

const int MAX_NUM_WALL = 200;
const int MAX_NUM_CYT = 60;


// Global Physics Constants

const double pi = acos(-1.0);

///// Cell wall mechanical parameters

const double thetaFlat = pi;
const double thetaCurve = pi * (11.0 / 12.0);

const double kBendHigh = 2;
const double kBendLow = 2;

const double kLinearHigh = 10;
const double kLinearLow = 8;

//linear spring equilibrium length
const double MembrEquLen = .0625; 
const double MEMBR_THRESH_LENGTH = 0.15; 

///// Subcellular element parameters for membrane - membrane interactions

const double U_MM = .2;
const double W_MM = .1;
const double xsi_MM = 0.06;
const double gamma_MM = 0.07;	

///// Subcellular element parameters for membrane  - internal interactions

const double U_MI = 0.78125;   //0.78125  
const double W_MI = 0.00;    //0.00;
const double xsi_MI = 0.125; //0.125;
const double gamma_MI = 0.425;  //0.625;

///// Subcellular element parameters for internal - internal interactions

const double U_II = 0.488;
const double W_II = 0.146484;
const double xsi_II = 0.3125;
const double gamma_II = 1.25;


//=====================
#endif
