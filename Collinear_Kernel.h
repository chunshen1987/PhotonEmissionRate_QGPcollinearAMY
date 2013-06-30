#ifndef Collinear_Kernel_h
#define Collinear_Kernel_h

#include "ParameterReader.h"

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

      int npt_p_plus;

      int npt_ktilde;
      double **fTilde;
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

      double SovleDiffeq(double ktilde, double p_plustilde);
      void calculateGluonselfEnergy(double q_0, double q, gluonSelfenergy* gluon_ptr);

      double fermiDistribution(double ptilde);

};

#endif
