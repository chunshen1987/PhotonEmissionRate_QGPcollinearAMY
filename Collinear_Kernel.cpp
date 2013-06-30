#include<iostream>
#include<sstream>
#include<fstream>
#include<cmath>
#include<iomanip>
#include<string>
#include<fstream>

#include <gsl/gsl_errno.h>
#include <gsl/gsl_matrix.h>
#include <gsl/gsl_odeiv2.h>
#include <gsl/gsl_sf_bessel.h>

#include "Collinear_Kernel.h"
#include "ParameterReader.h"
#include "Arsenal.h"
#include "Physicalconstants.h"

using namespace std;

int Diffeq_functions(double t, const double y[], double f[], void *params)
{
   double gammaEuler = 0.577216;
   double kappa = ((double *)params)[0];
   double coeff = ((double *)params)[1];
   double K0 = gsl_sf_bessel_K0(-t);
   double Cprime = 4.*kappa*1./(2.*M_PI)*(gammaEuler + log(-t/2.) + K0);
   f[0] = y[2];
   f[1] = y[3];
   f[2] = - 3.*y[2]/t + kappa*y[0] + coeff*Cprime*y[1];
   f[3] = - 3.*y[3]/t + kappa*y[1] - coeff*Cprime*y[0];
   return GSL_SUCCESS;
}

int Diffeq_Jacobian(double t, const double y[], double *dfdy, double dfdt[], void *params)
{
   double gammaEuler = 0.577216;
   double kappa = ((double *)params)[0];
   double coeff = ((double *)params)[1];
   double K0 = gsl_sf_bessel_K0(-t);
   double K1 = gsl_sf_bessel_K1(-t);
   double Cprime = 4.*kappa*1./(2.*M_PI)*(gammaEuler + log(-t/2.) + K0);
   double dCprimedt = 4.*kappa*1./(2.*M_PI)*(1./t + K1);
   gsl_matrix_view dfdy_mat 
     = gsl_matrix_view_array (dfdy, 4, 4);
   gsl_matrix * m = &dfdy_mat.matrix; 
   gsl_matrix_set (m, 0, 0, 0.0);
   gsl_matrix_set (m, 0, 1, 0.0);
   gsl_matrix_set (m, 0, 2, 1.0);
   gsl_matrix_set (m, 0, 3, 0.0);
   gsl_matrix_set (m, 1, 0, 0.0);
   gsl_matrix_set (m, 1, 1, 0.0);
   gsl_matrix_set (m, 1, 2, 0.0);
   gsl_matrix_set (m, 1, 3, 1.0);
   gsl_matrix_set (m, 2, 0, kappa);
   gsl_matrix_set (m, 2, 1, coeff*Cprime);
   gsl_matrix_set (m, 2, 2, -3./t);
   gsl_matrix_set (m, 2, 3, 0.0);
   gsl_matrix_set (m, 3, 0, - coeff*Cprime);
   gsl_matrix_set (m, 3, 1, kappa);
   gsl_matrix_set (m, 3, 2, 0.0);
   gsl_matrix_set (m, 3, 3, -3./t);
   dfdt[0] = 0.0;
   dfdt[1] = 0.0;
   dfdt[2] = 3.*y[2]/(t*t) + coeff*dCprimedt*y[1];
   dfdt[3] = 3.*y[3]/(t*t) - coeff*dCprimedt*y[0];
   return GSL_SUCCESS;
}

Collinear_Kernel::Collinear_Kernel(ParameterReader* paraRdr_in)
{
    paraRdr = paraRdr_in;

    npt_k = paraRdr->getVal("npt_k");
    double kT_min = paraRdr->getVal("k_min");
    double kT_max = paraRdr->getVal("k_max");
    double dkT = (kT_max - kT_min)/(npt_k-1);
    kT = new double [npt_k];
    for(int i = 0; i < npt_k ; i++)
       kT[i] = kT_min + i*dkT;

    npt_T = paraRdr->getVal("npt_T");
    double T_min = paraRdr->getVal("T_min");
    double T_max = paraRdr->getVal("T_max");
    double dT = (T_max - T_min)/(npt_T - 1);
    temperature = new double [npt_T];
    for(int i = 0; i < npt_T ; i++)
       temperature[i] = T_min + i*dT;
        
    rateTable = new double* [npt_k];
    for(int i = 0; i < npt_k ; i++)
       rateTable[i] = new double [npt_T];

    npt_ktilde = paraRdr->getVal("npt_ktilde"); 
    double ktilde_min = paraRdr->getVal("ktilde_min");
    double ktilde_max = paraRdr->getVal("ktilde_max");
    double dktilde = (ktilde_max - ktilde_min)/(npt_ktilde - 1);
    fTilde = new double* [npt_ktilde];
    ktilde = new double [npt_ktilde];
    rawRatetable = new double [npt_ktilde];
    
    npt_p_plus = paraRdr->getVal("npt_p_plus");
    for(int i = 0; i < npt_ktilde ; i++)
    {
       fTilde[i] = new double [npt_p_plus];
       ktilde[i] = ktilde_min + i*dktilde;
    }
}

Collinear_Kernel::~Collinear_Kernel()
{
    for(int i = 0; i < npt_ktilde; i++)
       delete [] fTilde[i];
    delete [] ktilde;
    delete [] fTilde;
    delete [] rawRatetable;
    
    delete [] kT;
    delete [] temperature;
    for(int i = 0; i < npt_k; i++)
       delete [] rateTable[i];
    delete [] rateTable;
}

void Collinear_Kernel::generateEmissionrateTable()
{
    Physicalconstants physconst;
    double alpha_EM = physconst.get_alphaEM();
    double C_F = physconst.get_C_F();
    double q_sq = physconst.get_q_sq();
    double g_s = physconst.get_g_s_const();

    for(int i = 0; i < npt_T; i++)
    {
       double m_inf_sq = C_F*g_s/4.*temperature[i]*temperature[i];
       double prefactor = 2.*alpha_EM*(3*q_sq)*m_inf_sq/pow(2*M_PI, 3);
       for(int j = 0; j < npt_k; j++)
       {
           double var = kT[j]/temperature[i];
           double result;
           interpolation1D_linear(ktilde, rawRatetable, &var, &result, npt_ktilde, 1);
           rateTable[j][i] = result*prefactor;
       }
    }
}

void Collinear_Kernel::outputEmissionrateTable(string filename)
{
    ofstream of;
    of.open(filename.c_str());
    for(int i = 0; i < npt_T; i++)
    {
       for(int j = 0; j < npt_k; j++)
          of << scientific << setw(20) << setprecision(8)
             << rateTable[j][i] << "   ";
       of << endl;
    }
    of.close();
}


double Collinear_Kernel::SovleDiffeq(double ktilde, double p_plustilde)
{
   Physicalconstants Phycons;
   double C_F = Phycons.get_C_F();
   double N_C = Phycons.get_N_c();
   double N_F = 3.0;
   double kappa = 3.*C_F/(4.*(N_C + N_F/2.));
   double *params = new double [3];
   params[0] = kappa;
   params[1] = 2.*p_plustilde*(p_plustilde + ktilde)/ktilde;
   gsl_odeiv2_system sys = {Diffeq_functions, Diffeq_Jacobian, 4, params};

   double relerr = 1e-8;
   double abserr = 0.0;
   double hstart = 1e-4;
   gsl_odeiv2_driver * d = gsl_odeiv2_driver_alloc_y_new (&sys, gsl_odeiv2_step_rk8pd, hstart, abserr, relerr);
   double t = -5.0;
   double y[4] = {0.01, 0.01, -0.1, -0.1};

   double t_end[2] = {-5e-3, -5e-4};
   double Re_h[2], Im_h[2];
   
   for(int i = 0; i < 2; i++)
   {
      int status = gsl_odeiv2_driver_apply (d, &t, t_end[i], y);
      if(status != GSL_SUCCESS)
      {
         cout << "ODE solver error!" << endl;
         exit(1);
      }
      //printf ("%.5e %.5e %.5e %.5e %.5e\n", -t, y[0], y[1], y[2], y[3]);
      Re_h[i] = y[0];
      Im_h[i] = y[1];
   }
   double temp = -sqrt(kappa)*t_end[0];
   double A_11 = - gsl_sf_bessel_K1(temp)/temp;
   double A_12 = - gsl_sf_bessel_I1(temp)/temp;
   temp = -sqrt(kappa)*t_end[1];
   double A_21 = - gsl_sf_bessel_K1(temp)/temp;
   double A_22 = - gsl_sf_bessel_I1(temp)/temp;

   double denominator = A_11*A_22 - A_21*A_12;
   double Re_C1 = (Re_h[0]*A_22 - Re_h[1]*A_12)/denominator;
   double Im_C1 = (Im_h[0]*A_22 - Im_h[1]*A_12)/denominator;
   double Re_C2 = (Re_h[1]*A_11 - Re_h[0]*A_21)/denominator;
   double Im_C2 = (Im_h[1]*A_11 - Im_h[0]*A_21)/denominator;

   double result = 2./M_PI*params[1]*kappa*(Re_C1*Im_C2 - Im_C1*Re_C2)/(Re_C1*Re_C1 + Im_C1*Im_C1);
   delete [] params;
   return(result);
}

void Collinear_Kernel::calculateGluonselfEnergy(double q_0, double q, gluonSelfenergy* gluon_ptr)
{
    double temp_log;
    temp_log = log((q + q_0)/(q - q_0));
    gluon_ptr->Re_PI_L = -1. + q_0/(2.*q)*temp_log;
    gluon_ptr->Im_PI_L = - M_PI * q_0/(2.*q);
    gluon_ptr->Re_PI_T = q_0*q_0/(2.*q*q) + (q*q - q_0*q_0)*q_0/(4.*q*q*q)*temp_log;
    gluon_ptr->Im_PI_T = - (q*q - q_0*q_0)*q_0/(4.*q*q*q)*M_PI;
}

double Collinear_Kernel::fermiDistribution(double ptilde)
{
    double result;
    result = 1./(exp(ptilde) + 1.);
    return(result);
}
