#include <iostream>
#include "DecayGenerator.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>


int main (int argc, char *argv[])
{
  
  DecayGenerator *gen = new DecayGenerator();
  gen->printModel();
  std::cout<<gen->rho_MM(100.0, 0.75)<<std::endl;
  std::cout<<gen->rho_MM(500.0, 0.75)<<std::endl; 
  std::cout<<gen->rho_MM(1220.0,0.75)<<std::endl; 
  return 0;
}

