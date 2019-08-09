// A simple program that computes the square root of a number
#include <iostream>
#include "DecayGenerator.h"
#include <gsl/gsl_sf_bessel.h>
#include <gsl/gsl_sf_result.h>

int main (int argc, char *argv[])
{
  
  DecayGenerator *gen = new DecayGenerator();
  gen->printModel();
  std::cout<<"F(0.001) = " << gen->F(0.001)<<std::endl;

  double x = 5.0;
  double y = gsl_sf_bessel_J0 (x);
  gsl_sf_result res;
  std::cout<<"gsl_sf_bessel_J0(5)"<<y<<std::endl;
  return 0;
}

