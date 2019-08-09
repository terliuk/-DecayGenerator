#include "DecayGenerator.h"

DecayGenerator::DecayGenerator(   std::string model_, // model name, e.g. MM 
                                  int zdaughter_,  // Z of the daughter, 56 for Xe->Ba decay
                                  double q_  // Q value of the decay
                                ){
    m_e = 510.998;  // keV
    alpha = 1.0/137.036; 
    setModel(model_);
    setZdaughter(zdaughter_);
    setQ(q_);
}

void DecayGenerator::setModel(std::string m){model = m; }
std::string DecayGenerator::getModel(){return model; }

void DecayGenerator::printModel(){
    std::cout<<"=== Generating decay mode === "<<std::endl;
    std::cout<<"Current model : "<<model<<std::endl;
    std::cout<<"Daughter Z    : "<<Z_d<<std::endl; 
    std::cout<<"Q value       : "<<Q<<" keV"<<std::endl;
    std::cout<<"============================= "<<std::endl;
}

void DecayGenerator::setZdaughter(int z_){ Z_d = z_;}
int DecayGenerator::getZdaughter(){return Z_d;}

void DecayGenerator::setQ(double q_){ Q = q_;}
double DecayGenerator::getQ(){return Q;}

// momentum of electron in units of electron mass
double DecayGenerator::p(double t_){ return sqrt(t_*(t_+2) ); } 
// speed of the electron in units of c
double DecayGenerator::beta(double t_){ return (p(t_)/(t_ +1)); } 

double DecayGenerator::F(double t_){return F(t_, Z_d);}

double DecayGenerator::F(double t_, int z_){ 
    long double s = sqrt(1.0 - pow((alpha*Z_d),2));
    long double u = alpha*Z_d*(t_ +1)/p(t_); 
    
    //std::cout<<"s = "<<s<<std::endl;
    //std::cout<<"u = "<<u<<std::endl;
    gsl_sf_result *lnf = new gsl_sf_result ; 
    long double result = (pow(p(t_), 2*s -u )*
                     exp(M_PI*u)  );
    return double(result);
}
