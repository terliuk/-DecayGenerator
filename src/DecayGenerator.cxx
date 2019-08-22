#include "DecayGenerator.h"

DecayGenerator::DecayGenerator(   std::string model_, // model name, e.g. MM 
                                  int zdaughter_,  // Z of the daughter, 56 for Xe->Ba decay
                                  double q_, 
                                  uint64_t seed  // Q value of the decay
                                ){initialize(model_,zdaughter_,q_,seed);}

//DecayGenerator::DecayGenerator(std::string model_){initialize(model_,56,2457.8,5489u);}

/// Initialization function 
void DecayGenerator::initialize(std::string model_, 
                                  int zdaughter_,
                                  double q_, 
                                  uint64_t seed ){
    m_e = 510.998;  // keV
    alpha = 1.0/137.036; 
    setModel(model_);
    setZdaughter(zdaughter_);
    setQ(q_);
    rng = std::mt19937_64(seed);
    unif = std::uniform_real_distribution<double>(0.0, 1.0);
}

void DecayGenerator::setModel(std::string m){
    model = m; 
    
    if (model==std::string("MM") ){         maxval = 1515.0; } 
    else if (model==std::string("RHC") ){ maxval = 9950.0;} 
    else if (model==std::string("2vbb") ) { maxval = 48600; }
    else{maxval = -1.0;}
    assert (("Cannot find decay model" ,maxval > 0.0) ) ;
}
std::string DecayGenerator::getModel(){return model; }

void DecayGenerator::printModel(){
    std::cout<<"=== Generating decay mode === "<<std::endl;
    std::cout<<"Current model : "<<model<<std::endl;
    std::cout<<"Daughter Z    : "<<Z_d<<std::endl; 
    std::cout<<"Q value       : "<<Q<<" keV"<<std::endl;
    std::cout<<"Max rho val   : "<<maxval<<std::endl;
    std::cout<<"============================= "<<std::endl;
}

void DecayGenerator::setZdaughter(int z_){ Z_d = z_;}
int DecayGenerator::getZdaughter(){return Z_d;}

void DecayGenerator::setQ(double q_){ Q = q_;}
double DecayGenerator::getQ(){return Q;}

void DecayGenerator::setMax(double mv_){ maxval = mv_;}
double DecayGenerator::getMax(){return maxval;}
// momentum of electron in units of electron mass
double DecayGenerator::p(double t_){ return sqrt(t_*(t_+2) ); } 
// speed of the electron in units of c
double DecayGenerator::beta(double t_){ return (p(t_)/(t_ +1)); } 

double DecayGenerator::F(double t_){return F(t_, Z_d);}

// This is a Fermi function that 
double DecayGenerator::F(double t_, int z_){ 
    double s = sqrt(1.0 - pow((alpha*z_),2));
    double u = alpha*z_*(t_ +1)/p(t_); 
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
//

void DecayGenerator::GenerateEvents(int nevents, 
                                    double* T1s, 
                                    double* T2s, 
                                    double* costhetas){
    for (int i =0; i < nevents; i++){
        double t1_ , t2_, cth_; 
        std::tie(t1_,t2_, cth_) = GenerateOneEvent();
        T1s[i] = t1_;
        T2s[i] = t2_; 
        costhetas[i] = cth_;
    }    
    return;
}
//  
std::tuple<double,double,double> DecayGenerator::GenerateOneEvent(){
    while (true){  // XXX: what about adding a safeguard against infinite loops? 
        double t1 = unif(rng)*Q;
        double t2 = (model == std::string("2vbb") ) ? unif(rng)*Q : t2 = Q - t1;
        double cth = (-1 + 2*unif(rng));
        double check = unif(rng)*maxval;
        double fval = -1;
        if (model==std::string("MM") )        {  fval = rho_MM(t1, cth);      } 
        else if (model==std::string("RHC") )  {  fval = rho_RHC(t1, cth);     } 
        else if (model==std::string("2vbb") ) {  fval = rho_2vbb(t1,t2, cth); }
        assert(( "rho value has wrong value, are you using a defined decay model? "&& fval >= 0.0)) ; 
        assert( ( "Rho value is larger than maximally allowed. You can use setMaximum() function, but this is strange" && fval < maxval ));
        if (check < fval ){ return std::make_tuple(t1,t2,cth);} 
    }
}
//
boost::python::tuple DecayGenerator::GenerateOneEventPy(){
    double t1_ , t2_, cth_; 
    std::tie(t1_,t2_, cth_) = GenerateOneEvent();
    return boost::python::make_tuple(t1_,t2_, cth_);
}
//

boost::python::numpy::ndarray DecayGenerator::GenerateEventsPy(int nev){
    namespace p = boost::python;
    namespace np = boost::python::numpy;
    double events_c[nev][3];

    for(int i=0; i < nev ; i++){
        std::tie(events_c[i][0], events_c[i][1], events_c[i][2]) = GenerateOneEvent();
        std::cout<<i<<"\t"<<events_c[i][0]<<"\t"<<events_c[i][1]<<"\t"<<events_c[i][2]<<std::endl;
    }
    np::ndarray events = np::from_data(events_c,
                                    np::dtype::get_builtin<double>(), 
                                    p::make_tuple(nev,3), // shape is N x 3
                                    p::make_tuple(3*sizeof(double),1*sizeof(double)), //stride, each of new row is 3 x 1
                                    p::object()); 
    return events.copy();
}
