#include <string>
#include <iostream>
#include <tuple>

#define _USE_MATH_DEFINES
#include <math.h> 
#include <gsl/gsl_sf_result.h>
#include <gsl/gsl_sf_gamma.h>
#include <random>
#include <assert.h>

class DecayGenerator{
    
    public: 
        
        
        DecayGenerator(std::string s= "MM", 
                      int z= 56, double q = 2457.8, uint64_t seed = 5489u);
        //DecayGenerator(std::string);
        // 
        void initialize(std::string, int, double, uint64_t);
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
        void setMax(double);
        double getMax(); 
        //
        double p(double);
        double beta(double);      
        double F(double);      
        double F(double, int);    
        //
        double rho_MM(double, double);
        double rho_RHC(double, double); 
        double rho_2vbb(double, double, double); 
        void GenerateEvents(int, double*, double*, double*);
        std::tuple<double,double,double> GenerateOneEvent();   
     
    private:
        // variables for overloading 
        
        ///
        std::string model; // Model: MM, RHC, 2vbb etc. 
        int Z_d; // int of daughter Z
        double Q;
        double alpha; 
        double m_e;
        double maxval;
        std::mt19937_64 rng;
        std::uniform_real_distribution<double> unif; 
};
