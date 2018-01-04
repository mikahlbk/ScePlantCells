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
const double dt = .00075;
const int ADD_WALL_TIMER = 2;
//const int ADD_CYT_TIMER = needs to be written to follow exponential growth;

const int Init_Num_Cyt_Nodes = 20;
const int Init_Wall_Nodes = 250;
const double AREA_THRESH = 48;
const double growth_rate = .0001;
// Global Physics Constants
const double pi = acos(-1.0);

///// Cell wall mechanical parameters
const double K_BEND = .0001;

//Adhesion spring mechanical params
const double K_ADH = 45;
const double MembrEquLen_ADH = 0.6;
const double ADHThresh = .78525;

//Calibration parameters
const double EXTERNAL_FORCE = 2;
const double INIT_AREA = .018;
const double GROWTH = 0.018;

//linear spring equilibrium length
const double MembrEquLen = .0625; 
const double MEMBR_THRESH_LENGTH = 0.1; 


///// Subcellular element parameters for membrane - membrane interactions
const double U_MM =  3.90625;
const double W_MM =  0;
const double xsi_MM = 0.125;
const double gamma_MM = 1.5625;	

///// Subcellular element parameters for membrane  - internal interactions
const double U_MI = 8;
const double W_MI = 0;
const double xsi_MI = 0.27;
const double gamma_MI = .625;

///// Subcellular element parameters for internal - internal interactions
const double U_II = 21;
const double W_II = 6;
const double xsi_II = 0.58;
const double gamma_II = 1.34;


//=====================
#endif
