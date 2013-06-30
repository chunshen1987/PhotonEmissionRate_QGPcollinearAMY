#ifndef Physicalconstants_H
#define Physicalconstants_H

class Physicalconstants
{
   private:
      double hbarC;
      double alpha_EM;
      double e_sq;
      double q_sq;
      double N_c;
      double C_F;

      double g_s;

   public:
      Physicalconstants();
      ~Physicalconstants();

      double get_hbarC() {return(hbarC);};
      double get_alphaEM() {return(alpha_EM);};
      double get_e_sq() {return(e_sq);};
      double get_q_sq() {return(q_sq);};
      double get_C_F() {return(C_F);};
      double get_N_c() {return(N_c);};
      double get_g_s_const() {return(g_s);};
};

#endif
