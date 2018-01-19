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
#include "rand.h"
//=====================

// Simulation Constants
const double dt = .00075;

const int Init_Num_Cyt_Nodes = 20;
const int Init_Wall_Nodes = 300;
const double AREA_THRESH = 64;
const double GROWTH_RATE = .000001;
// Global Physics Constants
const double pi = acos(-1.0);

///// Cell wall mechanical parameters
const double K_BEND = 3;

//Adhesion spring mechanical params
const double K_ADH = 20;
const double MembrEquLen_ADH = 0.6;
const double ADHThresh = .65;

//Calibration parameters
const double EXTERNAL_FORCE = 1;
const double INIT_AREA = .018;
const double GROWTH = 0.018;

//linear spring equilibrium length
const double MembrEquLen = .0625; 
const double MEMBR_THRESH_LENGTH = 0.095; 


///// Subcellular element parameters for membrane - membrane interactions
const double U_MM =  3.90625;
const double W_MM =  0;
const double xsi_MM = 0.125;
const double gamma_MM = 1.5625;	

///// Subcellular element parameters for membrane  - internal interactions
const double U_MI = 2.78288;
const double W_MI = 2.46484;
const double xsi_MI = .125;
const double gamma_MI = 1.55;

///// Subcellular element parameters for internal - internal interactions
const double U_II = 5.288;
const double W_II = 8.46484;
const double xsi_II = .8125;
const double gamma_II = 3;
const double ALPHA = 1.75;
const double DELTA = .35;
const double MORSE_EQ = 2;
const double ALPHA_MI = 1.2;
const double DELTA_MI = .15;
const double MORSE_EQ_MI = 1.5;
//=====================
#endif
