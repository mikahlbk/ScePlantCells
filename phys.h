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
const double dt = 0.0001;


// Global Physics Constants

const double pi = acos(-1.0);

///// Cell wall mechanical parameters

//equilibrium angle for end node
const double thetaEnd = pi;

//equilibrium angle for corner node
const double thetaCorner = pi / 2;

//equilibrium angle for flank node
const double thetaFlank = pi;

//rotational spring constant for end node
const double kBendEnd = 6;

//rotational spring constant for flank node
const double kBendFlank = 0.6;

//linear spring constant for end node
const double kLinearEnd = 30;

//linear spring constant for flank node
const double kLinearFlank = 6;

//linear spring equilibrium length
const double MembrEquLen = 0.0625;

///// Subcellular element parameters for membrane - membrane interactions

const double U_MM = 3.90625;
const double W_MM = 3.90625;
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
