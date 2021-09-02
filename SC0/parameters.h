#ifndef PARAMETERS_H
#define PARAMETERS_H
#include <string>
#include <iostream>
#include <vector>
#include "setparam.h"

using namespace std;
enum { P_silent,P_affinity_change,P_harmful,P_lethal,P_beneficial,P_key_beneficial,initial_koff,beneficial_change,harmful_change,key_beneficial_change};

//Name of parameters
enum {
//A
    agSaturation ,
    AgAmountperFDC ,
    Ag_threshold,
//B
    Bcell_speed ,
    Bcell_tp ,
    Bcell_stddev ,
    Bcell_tp_stddev ,
    BCR_pool,
    BCR_Length ,

//C
    chemo_dx ,
    chemmax ,
    chemosteep ,
    chemohalf ,
    CXCL12crit ,
    CXCL13crit ,
    CXCL12recrit ,
    CXCL13recrit ,
    c_G1 ,
    c_S ,
    c_G2 ,
    c_M ,
    c_G1_stddev ,
    c_S_stddev ,
    c_G2_stddev ,
    c_M_stddev ,
    collectionFDCperiod ,
    Ccdif_delay_stddev ,
    CB_radius,
//D
    difDelay ,
    dimension,
    dx ,
    dt ,
    DendriteLength ,
    DeleteAgInFreshCC ,
//E
    expMin ,
    expMax ,
    eta ,
//G
    Gamma ,
//I
    InitialNumberSC ,
    InitialNumberFDC ,
    InitialNumberTC ,
    InitialNumberCB ,
//K
    kon ,
    koff,

//M
//    Memorycell_tp ,
//    Memorycell_speed ,
//    Memorycell_tp_stddev ,
    macrophage ,
//N
    Nmax ,
    nDiv ,
    nDiv_stddev ,
    nDivinflow ,
    NoMutFounderCells ,
//P
    Plasmacell_tp ,
    Plasmacell_speed ,
    Plasmacell_tp_stddev ,
    pSel ,
    pApoCC ,
    pApoCB ,
    p_dif ,
    pMHCdepHill ,
    pMHCdepMin ,
    pMHCdepMax ,
    pMHCdepK ,
    pmutB4StartMut ,
    pmutAfterStartMut ,
    pmutAfterSelection ,
    pmutAffinityExponent ,
    pDivideAgAssymetric ,
    polarityIndex ,
//R
    radius ,
    rateCBinflow ,
//S
    StartMutation ,
    smoothnessStopCBinflow ,
//T
    timeStopCBinflow ,
    Tcell_speed ,
    Tcell_tp ,
    Tcell_stddev ,
    Tcell_tp_stddev ,
    tcTime ,
    tcRescueTime ,
    tmax ,
    tolight ,
    testDelay ,
//W
    widthPI ,
//Z
    zoneRatioGC,
//------------------------End---------
    N_par // Parameter counter, should be at the end always
};

class parameters
{
public:
    parameters(int ,string );
    parameters();
    ~parameters();
    vector <double> par;
    vector <string> names;
    void init_param();
    void writeparameters(string fname, short simulation_id_tmp);
    bool readparameters(string fname);
    void Set_additional_parameters(string additional_parameters_path);
    void convert_parameters ();
    string parameter_file_name;
    void matchFromHyphasma(hyphasmaParameter &hypar);
    vector<double> additional_parameters;
    double Avogadro_constant;
    string print();
};


#endif // PARAMETERS_H
