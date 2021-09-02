#ifndef CELL_H
#define CELL_H
#include <vector>
#include <string>
#include "vector3d.h"
#include "parameters.h"
#include "bcr.h"
#include <utility>
#include <sstream>
//The number Pi
 #define PI 3.141592654

using namespace std;

//Declearations
class lattice;
class T_cell;
class output;
class simulation;

//Cell types
enum celltype {empty,Stromalcell,FDCell, TFHC, Centroblast, Centrocyte, Plasmacell,Memorycell, border, cell_type_counter };

//Cell cycles
enum CellCycleState {cycle_G1, cycle_S, cycle_G2, cycle_M, cycle_Divide, cycle_G0, cycle_Ncellstates};

//Cell states
enum Cellstate {founder,unselected,contact_FDC,FDC_selected,contact_TC,TC_selected,recycled,apoptosis,TC_free,TC_connected,Plasma_Out,Plasma_in_GC,cell_state_counter};

class cell
{
public:
    cell(); //cell with new ID
    cell(cell* copied_cell);    //Constructs a new cell and copies fields from another cell to it
    virtual ~cell(){}   //Deconstructor
    int ID;     // ID of current cell
    int MID;    // ID of the mother cell
    double Born_time;
    Cellstate cell_state; //Status of the cell
    celltype cell_type; //Cell type
    vector3D position; // Current position
    vector3D polarity; // Current direction of movement
    double persistence_time; // Time left for next turn
    double speed;   //Speed of cell
    stringstream event; //records info of cell
    bool can_move; //A switch to turn moving on/off
    CellCycleState cyclestate;
    void getNewPersistentTime(parameters& p); //Update time left to calculate next polarity
    void getRandomPolarity(parameters& p, lattice& l); //Get a random polarity
    virtual void getNewPolarity(parameters& p, lattice& l); //Get a new polarity (not random)
    string print_neighbours(lattice &l); //Prints neighbours of the cell on current lattice l
    void move(parameters &p,  lattice& l , vector<vector3D> &redo_list); //move
    string printcell(); //Prints cell info
    double set_speed(parameters &p); // To change cell speed
    int get_new_clonal_id();
};

class B_cell: public cell {
public:

    B_cell(parameters& p); // Get new B_cell with random BCR.
    B_cell(parameters& p, B_cell* Mom_cell); // Copy everything from mother cellto daughter.
    ~B_cell(){}
    BCR myBCR;
    int total_number_of_divisions;
    int nFDCcontacts; // number of FDC contacts, not binding.
    double MyAffinity ;
    double pMHC_dependent_number_of_divisions;
    double cycle_state_time; //Time that has passed in the current cycle state
    double time_of_cycle_state_switch;   //Total time of current cycle state that cell has to pass to go to the next cycle//s
    double Bc_Tc_interaction_clock;
    double Recycling_delay; // Time it takes for B cell to recycle after getting selected by T cells
    double BC_FDC_interaction_clock; // Time since a B_cell becomes CC_free (in sec).
    double TC_selected_clock; // Time since a B_cell becomes selected by a Tcell
    double clock;
    //    double FDCinteractiontime; //Time since CC became in contact Ag. Swich on when using network
    double retained_Ag; // Ag internalized from interaction with FDC.
    bool IamHighAg; //Recycled cell will become output
    bool   Selected_by_FDC; //Indicates if B-cell rescued by FDC
    bool   Selected_by_TC;
    T_cell* interactingTC; // ID of T cell with who Bcell is interacting.
    CellCycleState cyclestate;
    int nDivisions2do; //Number of divisions left for cell to do
    double TCsignalDuration; // Acumulated signal from currently interacting TC (in sec).(As imput to ODE)
    double fdc_interaction_time_history;
    int fdc_binding_nums;
    int fdc_binding_nums_no_ag;

    pair<double ,double> Tc_interaction_history;
    vector3D fdc_pos={0,0,0};
    bool isResponsive2CXCL12;
    bool isResponsive2CXCL13;
    
    void setMyAffinity(parameters &p);
    void transmit_CCdelay2cycle(parameters &p);
    void ContinueCellCycle(parameters &p);
    void clockreset();//reset clocks if CCrecycled but not differentiated to output
    void timeleft2recycle(parameters & p); //When finished dividing CBs calculate a remining time to differentiate to CC_free.
    void set_Retained_Ag(parameters& p); 
    void mutate(parameters& p); //Mutation happens in BCR
    void proliferate(parameters &p,lattice &l,double time,vector<B_cell*> &ListB_cell, output &currentoutput, simulation &curent_sim);
    void Resensitize2Chemokines(parameters& p, lattice& l);
    void getNewPolarity(parameters& p, lattice& l);
    bool IsCycling() ;
    
    int clonal_id;
    int m_clonal_id;
    string printBcell(); 
};

void redo_move(vector <vector3D> &redo_list,lattice &l,double time);

class T_cell: public cell {
public:
    T_cell(parameters& p);
    ~T_cell(){}
    Cellstate cell_state;
    vector<B_cell*> interactingCC;
    int nIncontactCCs;
    void liberateCC_TC(B_cell *bc);
    string printTcell();
    void getNewPolarity(parameters& p, lattice& l);

};

class FDC: public cell {
public:
     FDC() : cell() {volume = 0; AgperDendrite = 0;}
      ~FDC(){}
    vector<vector3D> occupiedPositions;
    int volume;
    double AgperDendrite;
    bool can_move;
};

class Stromal_cell: public cell {
public:
    Stromal_cell() : cell() {}
    bool can_move;
     ~Stromal_cell(){}
};

class Memory_cell: public cell {
public:
    Memory_cell(parameters & p);
    Memory_cell(parameters & p , B_cell* Bcell);
    Memory_cell(B_cell* Bcell); // Copy fields from Bcell into M_cell.
     ~Memory_cell(){}
    BCR myBCR;
    double MyAffinity ;
    double BLIMP1;
    bool can_move;
    double IRF4;
};

class Plasma_cell: public cell {
public:
    Plasma_cell(parameters & p);
    Plasma_cell(parameters & p , B_cell* Bcell);
     ~Plasma_cell(){}
    BCR myBCR;
    bool   Selected_by_FDC; //Indicates if B-cell rescued by FDC
    bool   Selected_by_TC;
    double MyAffinity ;
    bool can_move;
    int total_number_of_divisions;
    double  birth_time;
    double fdc_interaction_time_history;
    int fdc_binding_nums;
    int fdc_binding_nums_no_ag;
    int nFDCcontacts; // number of FDC contacts.
    pair<double ,double> Tc_interaction_history;
    double retained_Ag;
    bool isResponsive2CXCL12;
    bool isResponsive2CXCL13;
    void getNewPolarity(parameters& p, lattice& l);
     
    int clonal_id; 
    int m_clonal_id;

};


#endif // CELL_H


