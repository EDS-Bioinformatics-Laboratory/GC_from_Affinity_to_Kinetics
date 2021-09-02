#ifndef OUTPUT_H
#define OUTPUT_H
//#include "observer.h"
#include "vector3d.h"
#include <iostream>
#include <fstream>
#include <time.h>
#include <stdlib.h>
#include <cmath>
#include <sstream>
#include <vector>
#include "cell.h"

#ifdef _WIN32
#include <windows.h>
#endif

#ifdef __linux__
#include <sys/stat.h>
#endif

#ifdef __APPLE__
#include <sys/param.h>
#include <sys/stat.h>
#endif

class simulation;
class parameters;
class cell;
enum {event_born, event_divide, event_die, event_CBdiff2CC, event_bind_FDC,event_pick_Ag,event_FDC_selected, event_start_TC_contact, event_stop_TC_contact, event_start_TC_signaling, event_stop_TC_signaling, event_tc_selected,event_MC_differentiation,event_PC_differentiation,event_mutation,evetn_recycling, events_counts};

class output
{
public:
    string parfname;
    output(string _parfname = string());
    ~output();
    void record_output_time_step(double currentTime, simulation& currentSim, parameters &p);
    void recordEvent(double currentTime  , double value);
    short Output_ID;
    void initialize_fileds();
    void clear_fileds();
    void write2file_time(double currentTime, simulation& currentSim, parameters&p);
    void write_event( cell* Cellx, stringstream &sim_output);
    void close_event(B_cell* Cellx, stringstream &sim_output,double time,int SIM_ID,int concept_id, int Sim_seed);

    void write_event_2file(stringstream &sim_output);

    void Plasma_output(double currentTime, simulation& currentSim, parameters&p);
    void createFolder(string folderName,parameters &p);
    string output_path;
    vector<string> storage;

};


struct ministat
{
    
    ministat()
    {
        clear_ministat();
    }
    void clear_ministat()
    {
        N=0;
        sum=0;
        sumsq=0;
        values.clear();
    }
    void add(double v)
    {
        if(not(std::isnan(v))){
            sum += v;
            sumsq += v*v;
            N++;}
        else {
            cout<<"NAN value in output!"<<endl;
        }
        values.push_back(v);
    }
    int N;
    vector<double> values;
    double sum;
    double sumsq;
    double average()
    {
        // if(N==NAN){return 0;}
        if (N<=0){return 0;}
        if(std::isnan(sum)){
            cout<<"Error in output calculation, NAN value, ERP01"<<endl;
            return 0;}
        else{return double(sum / (double) N);}
    }
    
    double stddev()
    {

        if (N<=1){return 0;}
        if(sum == 0){return 0;} 
        if (std::isnan(sum)) {
            cout<<"Error in output calculation, NAN value, ERP02"<<endl;
            return 0;}
        else
        {   if (N>1)
            
        {
            double std = (double) sqrt( (double) (sumsq / (double) N) - ((sum / (double) N) * (sum / (double)(N))));
            if (std::isnan(std))
            {return 0.0;}
            return std;
        }
        else {
            return 0;
        }
            
        }
    }
    string print()
    {
        stringstream res;
        res << "[" << N << "]";
        for(int i = 0; i < N; ++i)
        {
            res << "\t" << values[i];
        }
        return res.str();
    }
};

//ministats
extern vector<ministat> Bcell_counts;   //number of Bcells
extern ministat Plasma_counts;

//Mutations
extern ministat Bcell_mutation;
extern ministat Plasma_mutations;

//Affinities
extern ministat Bcell_affinity;
extern ministat Plasma_affinity;

//Antigens
extern vector<ministat> Bcell_Antigen;

//Antibodies
extern ministat Antibody;

//CentroBlasts cycle and divisions
extern ministat Bcell_cycle;
extern ministat Bcell_divisions;
extern ministat Plasma_antigen;
extern ministat Plasma_divisions;


#endif // OUTPUT_H
