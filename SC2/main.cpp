#include "mafalda.h"
#include "random.h"
//#include "GC3D.h"
#include <iostream>
#include <vector>
#include <string>
#include <getopt.h>
//#include <sys/stat.h>
using namespace std;
string outputFolder=string();
#include <time.h>
#ifdef _WIN32
#include <windows.h>
#endif
#ifdef __linux__
#include <sys/stat.h>
#endif
#ifdef __APPLE__
#include <sys/param.h>
//#include <boost/filesystem.hpp>
#endif
std::mt19937 generator;

unsigned int rdtsc(){
    unsigned int lo,hi;
    __asm__ __volatile__ ("rdtsc" : "=a" (lo), "=d" (hi));
    return ((unsigned int)hi << 16) | lo;
}

int main(int argc, char** argv){
 

// Parsing arguments given from command line
// See https://www.gnu.org/software/libc/manual/html_node/Getopt-Long-Option-Example.html#Getopt-Long-Option-Example
   string parfname = string();
// Using hyohasma parameter file
    int takeHyphasmaFile = false;
// Using specific seed for random number generator
    int requested_seed = 1;
//Simulation id
    short simulation_id = 1;
    short Concept_id = 109;
    string additional_parameter_path;
// will be the remaining number of arguments during parsing
    int c = 0;
    while (c != -1){
        // Definition of the list of possible arguments: either "a_word" or "-X argument".
        // Other arguments (not from the predefined list), like parameter files, will also be retrieved at the end.
        static struct option long_options[] = {
            {"hypster",     no_argument,       0, 'h'},   // These options don’t set a flag.
            {"seed",    required_argument, 0, 's'},
            {"outputFolder",    required_argument, 0, 'o'},
            {"simulation_id", required_argument,0,'i'},
            {"additional_parameters", required_argument,0,'p'},
            {0, 0, 0, 0}
        };
    
// Parses the arguments one by one
// The index of the identified argument will be put inside.
        int option_index = 0;
        c = getopt_long (argc, argv, "hs:o:i:p:", long_options, &option_index);  // "as:" means, expects -a without argument or -s with argument
        switch (c)
        {
        case 0:
            {if (long_options[option_index].flag != 0)  // i.e. if there was a flag associated. Nothing to do
                cout << long_options[option_index].name << endl;
            if (optarg)
                cout << " with arg " << optarg << endl;
                break;}
        case 'h':
            {
//                cout << "Detected   -h  -> will read the parameter file as hyphasma file" << endl;
            takeHyphasmaFile = true;
                break;}
        case 's':
            { requested_seed = atoi(optarg); //turn string into int
//            cout << "Detecetd   -s  -> using seed: " << requested_seed << endl;
//printf ("option -s with value `%s'\n", optarg);
                break;}
        case 'o':
            { outputFolder = optarg;
                
//            cout << "Detected  -o -> Using output folder: " << outputFolder << endl;
//printf ("option -s with value `%s'\n", optarg);
                break;}
                
            case 'i':
            {    requested_seed = atoi(optarg); //
                
//                cout << "Detected  -i  -> simulation ID: " << simulation_id << endl;
                //printf ("option -s with value `%s'\n", optarg);
                break;}
            
            case 'p':
            {     additional_parameter_path = optarg; //
                
                cout << "Detected  -p  -> Additional parameters: " << additional_parameter_path.c_str()<< endl;
                //printf ("option -s with value `%s'\n", optarg);
                break;}
                
        case '?':
            { // getopt_long should have printed an error message.
                break;}
        default:
            {// no arguments within the list. Other arguments might be given (see next loop)
            //cout << "no argument given" << endl;
                break;}
        }
    }
    // Additional arguments, that are not in the 'official' list options (the parameter file for instance).
    int nArguments = argc - optind;
   
    if((parfname.size() > 0) && (!parfname.substr(parfname.size()-4,4).compare(string(".par")))){
        parfname = parfname.substr(0, parfname.size()-4);
        cout << " The parameter file contained .par at the end => cutted it into " << parfname << endl;
    }

// If the output path is not decleared or set as automatic it sets the output folder name as the parameter file name
    if((outputFolder.size() == 0) || (!outputFolder.compare(string("auto")))){
        outputFolder = parfname;
    }
    
    
    
    unsigned  int seed_zero= requested_seed;
    srand(seed_zero);
    generator.seed(seed_zero);
    
    cout<<"Simulation seed (id) : "<<seed_zero<<endl;
// Print all parameters that is used to analyze file
// Create simulation instance
    cout<< "Reading internal parameters";
    parameters currentParameterSet;//(takeHyphasmaFile,parfname);
    cout<<", done."<<endl;
//      currentParameterSet.writeparameters("paramparam.txt");
        currentParameterSet.convert_parameters();
        currentParameterSet.writeparameters(outputFolder,simulation_id);
        simulation Sim(currentParameterSet);
        Sim.Simulation_ID=simulation_id;
        Sim.Sim_Seed=requested_seed;
        Sim.Concept_ID=Concept_id;
        Sim.currentOutput->Output_ID=requested_seed;
//      initGC3D(argc, argv);
    
        if (additional_parameter_path.size()>0)
        {
            currentParameterSet.Set_additional_parameters(additional_parameter_path);
        }
    
        Sim.simulate(*Sim.currentLattice, currentParameterSet);
    

    return 0;

}
