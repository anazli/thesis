#ifndef GUARD_PARAMETERS_H
#define GUARD_PARAMETERS_H

#include<cmath>
#include "vec3.h"

const unsigned int chunkSZ = 8; //omp for schedule static chunksize
const double PI = acos(-1.);
const size_t DIM = 3; //coords dim
const double L = 10.E-7;//m //box length
//energy cutoff distance
const double RCW = 2.E-8;//for van der waals
const double RCM = 6.E-8; //for magnetic dipolar
const double Rlj = 3.E-8; //for lenard jones
const double NM = 1.E-9; // one nano meter
const double MODAL = 7.; //modal diameter
const double SIGMA = 0.24; //dispersion
const double LEFT = 1.; //Left boundary of diameter's distribution
const double RIGHT = 10.; //Right boundary
//#####################################################################
const double Eo = 8.854E-12; //m^-3*kg^-1*s^4*A^2 permittivity of free space
const double Er = 80.1; //relative permittivity of water
const double charge = 1.6022E-19; //C elementary charge
const double Mo = 4. * PI * 1.E-7; //kg*s^-2*A^-2 permeability of free space
const double Mr = 1.; //relative permeability of water
//#####################################################################
const double Kb = 1.38064852E-23; //boltzmann constant
const double ALPHA = 5.E-20; //hamaker constant
const double grafting = 1.E+18; //nm;
const double TEMP = 1.0; //temperature
const double BETA = 1./(Kb * TEMP);
const double KANIS = 9000.;  // Anisotropy Constant
const double MSAT = 446000.; // Saturation Magnetization
const Vec3 field(0., 0., -0.01); //magnetic field vector 

#endif
