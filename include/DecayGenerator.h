#include <string>
#include <iostream>
#define _USE_MATH_DEFINES
#include <math.h> 
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>

class DecayGenerator{
    
    public: 
        DecayGenerator(std::string s= "MM", 
                      int z= 56, double q = 2457.8);
        // 
        void setModel(std::string);
        std::string getModel();
        void printModel(); 
        //
        void setZdaughter(int);
        int getZdaughter();
        // 
        void setQ(double);
        double getQ();
        // 
        double p(double);
        double beta(double);      
        double F(double);      
        double F(double, int);    
        //
        double rho_MM(double, double);
        double rho_RHC(double, double); 
        double rho_2vbb(double, double, double); 
            
    private:
        std::string model; // Model: MM, RHC, 2vbb etc. 
        int Z_d; // int of daughter Z
        double Q;
        double alpha; 
        double m_e;
        
};
