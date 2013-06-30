#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "Physicalconstants.h"

Physicalconstants::Physicalconstants()
{
    hbarC = 0.197327053;  //GeV*fm
    e_sq = 1./137.*4.*M_PI;
    q_sq = pow(2./3., 2) + pow(-1./3., 2) + pow(-1./3., 2); //u, d, s
    N_c = 3.;   // number of color
    C_F = (N_c*N_c - 1)/(2.*N_c);
    alpha_EM = 1./137.;
    g_s = 2.0;

}

Physicalconstants::~Physicalconstants()
{

}
