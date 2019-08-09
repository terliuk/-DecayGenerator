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
// This is a Fermi function that 
double DecayGenerator::F(double t_, int z_){ 
    double s = sqrt(1.0 - pow((alpha*Z_d),2));
    double u = alpha*Z_d*(t_ +1)/p(t_); 
    gsl_sf_result lnf, arg; 
    int k = gsl_sf_lngamma_complex_e(s, u, &lnf, &arg); 
    double result = (pow(p(t_), 2*s -2 )*
                     exp(M_PI*u + 2.0*lnf.val)  );
    return result;
}
// 
double DecayGenerator::rho_MM(double T1, double costheta){
    double t1 = T1 / m_e;
    double t2 = (Q - T1) / m_e;
    return ( (t1 + 1.0)*p(t1) * 
             (t2 + 1.0)*p(t2) * 
             F(t1) * F(t2) * 
             (1 - beta(t1)*beta(t2)*costheta) 
           );
    }
// 
double DecayGenerator::rho_RHC(double T1, double costheta){
    double t1 = T1 / m_e;
    double t2 = (Q - T1) / m_e;
    return ( (t1 + 1.0)*p(t1) * 
             (t2 + 1.0)*p(t2) * 
             F(t1) * F(t2) * 
             pow(t2 - t1, 2) *
             (1 + beta(t1)*beta(t2)*costheta) 
           );
    }
// 
double DecayGenerator::rho_2vbb(double T1, double T2, double costheta){
    if ( (Q - T1 - T2) < 0.0  ){return 0.0;} 
    double t1 = T1 / m_e ; 
    double t2 = T2 / m_e ; 
    
    return ( (t1 + 1.0)*p(t1) * 
             (t2 + 1.0)*p(t2) * 
             F(t1) * F(t2) * 
             pow( (Q / m_e - t1 - t2), 5)*
             (1 - beta(t1)*beta(t2)*costheta) 
           );
    }

