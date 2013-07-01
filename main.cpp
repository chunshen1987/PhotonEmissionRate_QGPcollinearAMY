#include <iostream>
#include <fstream>
#include <sstream>
#include <cmath>
#include <iomanip>
#include <stdio.h>
#include <stdlib.h>
#include <time.h>

#include "ParameterReader.h"
#include "Collinear_Kernel.h"
#include "Stopwatch.h"
#include "Arsenal.h"

using namespace std;

int main(int argc, char** argv)
{
   Stopwatch sw; 
   sw.tic();

   ParameterReader* paraRdr = new ParameterReader();
   paraRdr->readFromFile("parameters.dat");
   paraRdr->readFromArguments(argc, argv);
   
   Collinear_Kernel test(paraRdr);
   //Collinear emission
   string filename = "CollinearLPM";
   test.calculatePhotonEmissionRates(filename);
   test.calRawEmissionTable();
   //double result = test.SolveDiffeq(50.0, -25.0);
   //cout << result << endl;
   //test.calRawEmissionTable();

   sw.toc();
   cout << "totally takes : " << sw.takeTime() << "sec." << endl;
   return 0;
}
