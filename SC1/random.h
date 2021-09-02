#ifndef i_random
#define i_random
#include <random>
#include <vector>
using namespace std;
extern std::mt19937 generator;

struct random
{

    static int randomInteger (int); 
    static int randomInteger(int a,int b);
    static double randomDouble (double);
    static double randomDouble(double fMin, double fMax);
    static int randomNormalint (int a,int b);
    static double randomWeibull(double a, double b);
    static double randomNormaldouble (double a, double b);


    static double cell_cycle_time (double,int);

};
double inverse_erf(double x);

class gauss_randomize {
public:
    gauss_randomize();
    gauss_randomize(short dataset);
    ~gauss_randomize();
    double get_distribution_value();
private:
    void gauss_initialize(double, double, int, double, double);
    double gauss(double&,double&,double&);
    void cyster_initialize();
    double cyster07angle_wt(int angle);
    vector <double> field;
    int arraydim;
    
    ////§§§    21-03-2017
    static const double Ny;
};


#endif
