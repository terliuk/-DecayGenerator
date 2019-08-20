#include <iostream>
#include "DecayGenerator.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>


int main (int argc, char *argv[])
{
  //std::cout<<argv[1]<<std::endl;
  std::string m = (argc > 1) ? std::string(argv[1]) : std::string("MM");
  DecayGenerator *gen = new DecayGenerator(std::string(m));
  gen->printModel();
    
  int nev = 50; 
  double T1[nev], T2[nev], costhetas[nev];
  gen->GenerateEvents(nev, T1, T2, costhetas);

  for(int i=0 ; i<nev; i++){
    std::cout<< i<<std::endl;
    std::cout << "T1 = " << T1[i] << std::endl;
    std::cout << "T2 = " << T2[i] << std::endl;
    std::cout << "CT = " << costhetas[i] << std::endl;
    std::cout << T1[i] + T2[i] << std::endl;
  }
  return 0;
}

