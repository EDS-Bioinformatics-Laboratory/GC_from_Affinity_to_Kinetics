#include "bcr.h"
#include "random.h"
#include <cmath>
#include <sstream>
using namespace std;
#define PI 3.14159265
//vector <pair<double ,double >> initial_affinity_seeds;
vector <vector <int>> initial_affinity_seeds;
vector <long double> delta_distribution_x;
// Constructor of BCR for the shape-space method
BCR::BCR(parameters& p)
{   /// Initialize mutation counters
    nMutFromGermline.beneficial = 0.;
    nMutFromGermline.harmful = 0.;
    nMutFromGermline.silent = 0.;
    nMutFromGermline.key_mutations_beneficial = 0.;
    nMutFromGermline.lethal = 0.;
    P_on=0.0;
    P_off=0.0;
    ///Initialize distance
    distance=5.0;
    checked_initial_dissociation=false;
    BCReceptor.resize(p.par[BCR_Length],0);
    double tmp_rnd = random::randomDouble(1);
    if (tmp_rnd<0.33)
    {BCReceptor={4,5,4,2};}
    else if (tmp_rnd<0.66)
    {BCReceptor={4,5,2,4};}
    else BCReceptor={2,5,1,3};
    m_theta=0;
    delta_Pon=0;
    delta_Poff=0;
    delta_mtheta=0;
    delta_Affinity=0;
}

void BCR::calculate_on_off_P(parameters &p)
{
     double dist = 0.;
    for( int i = 0; i < BCReceptor.size(); i++){
        dist += fabs(double (3. - BCReceptor[i])); //   now the target antigen is 3333, modify later
    }
    distance=dist;
    
    double d_a = distance * sin(m_theta*PI/180.0);
    P_on  = exp((-1.0*powf(d_a, 2))/powf(2.8, 2));
    
    double d_d = distance * cos(m_theta*PI/180.0);
    P_off = exp((-1.0*powf(d_d, 2))/powf(2.8, 2));

    if (P_off>1)
    {P_off=1.0;
        cout<<"Error: Poff bigger than one"<<endl;
    }
    if (P_off<0)
    {P_off=0.0;
        cout<<"Error: Poff smaller than zero"<<endl;
    }
    if (P_on>1.0)
    {P_on=1.0;
        cout<<"Error: Pon bigger than one"<<endl;
    }
    if (P_on<0.0)
    {P_on=0.0;
        cout<<"Error: Pon smaller than zero"<<endl;
    }
}
void BCR::calculate_koff(parameters &p)
{
//    double dist = 0.;
//   for( int i = 0; i < BCReceptor.size(); i++){
//       dist += fabs(double (3. - BCReceptor[i])); //   now the target antigen is 3333, modify later
//   }
//    distance=dist;
//       k_off = exp(distance-7.0);
}
void initialize_Seeds (parameters &p)
{
    initial_affinity_seeds.resize(p.par[BCR_pool],vector<int>(p.par[BCR_Length]));
       for (int j=0;j<p.par[BCR_pool];j++){
           for (int i =0; i< p.par[BCR_Length];i++){
               initial_affinity_seeds[j][i]=random::randomInteger(p.par[Nmax]+1);
           }}
    //Initialize delta_kon/koff distribution
    delta_distribution_x.reserve(10000);
    for (int i=0;i<10000;i++){
        long double number = random::randomWeibull(15,10);
        delta_distribution_x.push_back((number-11)/2.); }
}
void BCR::mutateBCR(parameters& p)
{
    delta_mtheta=0;
    delta_Poff=0;
    delta_Pon=0;
    delta_Affinity=0;
    calculate_on_off_P(p);
    double tmp_on_1 =P_on;
    double tmp_off_1 =P_off;
    double tmp_mtheta_1 =m_theta;
    double tmp_aff1 = getMyAffinity4Ag(p);
    
    if(random::randomDouble(1) < pMut)
    {

    double step=0;
    step= 1.0;
    short int randomL; //random location
    short int randomC; //random change
    randomL = random::randomInteger(4);
    randomC = random::randomInteger(2);
    while ((BCReceptor[randomL]==9 && randomC==1)|| (BCReceptor[randomL]== 0 && randomC==0))
    {
    randomL = random::randomInteger(4);
    randomC = random::randomInteger(2);
    }
    if (randomC==1)
    {BCReceptor[randomL] += step;}
    if (randomC==0)
    {BCReceptor[randomL] -= step;}
    }
    calculate_on_off_P(p);
    delta_Affinity = getMyAffinity4Ag(p) - tmp_aff1;
    delta_mtheta=m_theta - tmp_mtheta_1;
    delta_Pon= P_on - tmp_on_1;
    delta_Poff = P_off - tmp_off_1;
    
    if (delta_Affinity>0){
    nMutFromGermline.beneficial++;}
    else if (delta_Affinity<0){
    nMutFromGermline.harmful++; //Increase numer of mutations from the beginning founder
    }
    //Silent
    else{nMutFromGermline.silent = nMutFromGermline.silent + 1;}
    return; // Returns
}
//This function calculates affinity using shape-space concept by computing the distance between Ag and bcr
double BCR::getMyAffinity4Ag (parameters &p)
{
//    calculate_koff(p);
//    if (k_off==0)
//    {cout<<"Error, Koff=0, infinit affinity."<<endl;
//        exit(1);
//    }
//    if (checked_initial_dissociation)
//    {return 0.0;}
    double dist = 0.;
       for( int i = 0; i < BCReceptor.size(); i++){
           dist += fabs(double (3. - BCReceptor[i])); //   now the target antigen is 3333, modify later
       }
       distance=dist;
    return (long double)(exp((-1.0*powf(distance, 2))/powf(2.8, 2))); 
}

string BCR::print_BCR(){
    stringstream res;
    int L = (int) BCReceptor.size();
    for(int i = 0; i < L; ++i)
    {
//        res<<BCReceptor[i];
    }
    return res.str();
}




