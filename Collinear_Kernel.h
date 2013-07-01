#ifndef Collinear_Kernel_h
#define Collinear_Kernel_h

#include "ParameterReader.h"
#include "Physicalconstants.h"

struct gluonSelfenergy
{
   double Re_PI_L;
   double Im_PI_L;
   double Re_PI_T;
   double Im_PI_T;
};

class Collinear_Kernel
{
   private:
      ParameterReader *paraRdr;

      Physicalconstants *Phycons;

      int npt_p_plus_1, npt_p_plus_2;
      double *p_plus_pt_1, *wp_plus_1, *p_plus_pt_standard_1, *wp_plus_standard_1;
      double *p_plus_pt_2, *wp_plus_2, *p_plus_pt_standard_2, *wp_plus_standard_2;

      int npt_ktilde;
      double *ktilde;
      double *rawRatetable;

      int npt_k, npt_T;
      double *kT, *temperature;
      double **rateTable;

   public:
      Collinear_Kernel(ParameterReader* paraRdr_in);
      ~Collinear_Kernel();

      void generateEmissionrateTable();
      void outputEmissionrateTable(string filename);

      void calRawEmissionTable();
      double Integral_p_plus(double ktilde);
      double SolveDiffeq(double ktilde, double p_plustilde);
      void calculateGluonselfEnergy(double q_0, double q, gluonSelfenergy* gluon_ptr);

      double fermiDistribution(double ptilde);

};

#endif
