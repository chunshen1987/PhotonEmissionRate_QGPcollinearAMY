#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include "Physicalconstants.h"

Physicalconstants::Physicalconstants(ParameterReader* paraRdr_in)
{
    paraRdr = paraRdr_in;

    hbarC = 0.197327053;  //GeV*fm
    e_sq = 1./137.*4.*M_PI;
    N_c = 3;   // number of color
    N_F = paraRdr_in->getVal("N_flavor");
    C_F = (N_c*N_c - 1)/(2.*N_c);
    alpha_EM = 1./137.;
    g_s = 2.0;

}

Physicalconstants::~Physicalconstants()
{

}

double Physicalconstants::get_q_sq()
{
   if(N_F == 1)
      q_sq = pow(2./3., 2);
   else if(N_F == 2)
      q_sq = pow(2./3., 2) + pow(-1./3., 2);
   else if(N_F == 3)
      q_sq = pow(2./3., 2) + 2.*pow(-1./3., 2);
   else if(N_F == 4)
      q_sq = 2.*pow(2./3., 2) + 2.*pow(-1./3., 2);
   else if(N_F == 5)
      q_sq = 2.*pow(2./3., 2) + 3.*pow(-1./3., 2);
   else if(N_F == 6)
      q_sq = 3.*pow(2./3., 2) + 3.*pow(-1./3., 2);
   else 
   {
      cout << "N_F is not physical!" << endl;
      exit(1);
   }
   return(q_sq);
}
