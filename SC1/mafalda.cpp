#include "mafalda.h"
#include "random.h"
#include "lattice.h"
#include "cell.h"
//#include "GC3D.h"
#include <iostream>
#include <vector>
#include <list>
#include <cmath>
#include "bcr.h"
#include "output.h"
#include "performance.h"
#include <algorithm>
bool pause1 = false;
using namespace std;


/// Simulation initializer
/// @param p is list of parameters
simulation::simulation(parameters& p) {
  
  currentOutput = new output();
  currentLattice =
      new lattice(p, "cxcl12_3d_5micron.sig", "cxcl13_3d_5micron.sig");
  ListB_cell.reserve(50000);
  ListT_cell.reserve(50000);
  ListP_cell.reserve(50000);
  ListM_cell.reserve(50000);
}
/// Destructor of Simulation class
simulation::~simulation() {
  delete currentLattice;
  delete currentOutput;
  int NSC = (int)ListSC.size();
  for (int i = 0; i < NSC; ++i) {
    delete ListSC[i];
  }
  int NFDC = (int)ListFDC.size();
  for (int i = 0; i < NFDC; ++i) {
    delete ListFDC[i];
  }
  int NBC = (int)ListB_cell.size();
  for (int i = 0; i < NBC; ++i) {
    delete ListB_cell[i];
  }
  int NTC = (int)ListT_cell.size();
  for (int i = 0; i < NTC; ++i) {
    delete ListT_cell[i];
  }
  int NPC = (int)ListP_cell.size();
  for (int i = 0; i < NPC; ++i) {
    delete ListP_cell[i];
  }
  ListB_cell.clear();
  ListP_cell.clear();
  ListSC.clear();
  ListFDC.clear();
  ListT_cell.clear();
  ListM_cell.clear();
}

/// Initialize cells
/// @param l is lattice of simulation
/// @param p is parameter list
void simulation::InitialCells(lattice& l, parameters& p) {
    /*Initialize Stromal cells,randomley placed in dark zone,immobile, not transparent */
  for (unsigned int j = 0; j < p.par[InitialNumberSC]; j = j + 1) {
    Stromal_cell* newSC = new Stromal_cell();
    newSC->position = l.getFreePosition(p.par[zoneRatioGC], 1);
    ListSC.push_back(newSC);
    newSC->cell_type = Stromalcell;
    newSC->can_move = false;
    l.putcellat(newSC);
  }
  cerr << p.par[InitialNumberSC] << " SC generated " << endl;

  /*Initialize FDC network */
  /* If the dendrit goes beyond the simulation space the ag amount will be distributed over the other dendrities. To work with FDCs you need to work with the amount of Ag in each lattice node because FDCs are transparent*/
  l.set_initial_fdc_position();

  for (unsigned int j = 0; j < p.par[InitialNumberFDC]; j++) {
    FDC* newFDC = new FDC();
//    newFDC->position = l.getFreePosition(0, p.par[zoneRatioGC]);
    newFDC->position = l.get_fdc_position(j);
      if (not(l.insideBorders(newFDC->position)))
      {cout<<"Error, FDC position out of borders.";
          exit(1);
      }
      
    newFDC->can_move = false;
    newFDC->occupiedPositions.reserve(6 * p.par[DendriteLength]);
    newFDC->volume = 1;
    newFDC->cell_type = FDCell;
    newFDC->occupiedPositions.push_back(newFDC->position);
    int ttmp[3], tmp[3];
    ttmp[0] = newFDC->position.X;
    ttmp[1] = newFDC->position.Y;
    ttmp[2] = newFDC->position.Z;
    tmp[0] = ttmp[0];
    tmp[1] = ttmp[1];
    tmp[2] = ttmp[2];
    for (int i = 1; i <= (p.par[DendriteLength]); i++) {
      for (int j = 0; j < 3; j++) {
        // positive direction
        tmp[j] = ttmp[j] + i;
        if (l.insideBorders(vector3D(tmp[0], tmp[1], tmp[2]))) {
          newFDC->volume += 1;
          newFDC->occupiedPositions.push_back(vector3D(tmp[0], tmp[1], tmp[2]));
        }
        // negative direction
        tmp[j] = ttmp[j] - i;
        if (l.insideBorders(vector3D(tmp[0], tmp[1], tmp[2]))) {
          newFDC->volume += 1;
          newFDC->occupiedPositions.push_back(vector3D(tmp[0], tmp[1], tmp[2]));
        }
        tmp[j] = ttmp[j];
      }
    }
            
    // Amount of Ag per dendrite
    newFDC->AgperDendrite = p.par[AgAmountperFDC] / newFDC->volume;
    for (int i = 0; i < newFDC->volume; i++) {
      l.putAgFDCat(newFDC->occupiedPositions.at(i), newFDC,
                   newFDC->AgperDendrite);
      l.AddTotalAmountAginLattice(newFDC->AgperDendrite);
    }
    ListFDC.push_back(newFDC);
  }
  cerr << p.par[InitialNumberFDC] << " FDC generated" << endl;

  /// Initialize (100) affinity seeds for incoming CBs
  initialize_Seeds(p);

  /// Initialize the Seeder B cells, placed randomley in dark zone. B cells start the cycle in G1 phase with an affinity from initial seeds pool
    //Theta=90 -> low association/dissociation theta=0 high association/dissociation
    double  m_theta_values_founder[3]={90.0,45.0,0.0};
  for (unsigned int j = 0; j < p.par[InitialNumberCB]; j = j + 1) {
    B_cell* newB_cell = new B_cell(p);
    newB_cell->myBCR.m_theta=m_theta_values_founder[j];
    newB_cell->myBCR.calculate_on_off_P(p);
    newB_cell->clonal_id=j;
    newB_cell->m_clonal_id=j;
    newB_cell->myBCR.calculate_on_off_P(p);
    newB_cell->cell_state = founder;
    newB_cell->cell_type = Centroblast;
    newB_cell->persistence_time = p.par[Bcell_tp];  // Time left for next turn
    newB_cell->speed = p.par[Bcell_speed];
    newB_cell->Born_time=0.0;
    newB_cell->can_move = true;  // A switch to turn moving on/off
    newB_cell->setMyAffinity(p);
    newB_cell->time_of_cycle_state_switch = random::cell_cycle_time(p.par[c_G1], cycle_G1);
    newB_cell->cyclestate = cycle_G1;
    newB_cell->position = l.getFreePosition(p.par[zoneRatioGC], 1);
    newB_cell->polarity = l.get_random_direction();
    newB_cell->myBCR.checked_initial_dissociation=false;
    if (not(l.insideBorders(newB_cell->position))) {
      cerr << "Influx of B-cell at border position: " << newB_cell->printcell()
           << endl;
      exit(1);
    }
    newB_cell->nDivisions2do = p.par[nDiv];
    newB_cell->getNewPersistentTime(p);
    newB_cell->isResponsive2CXCL12 = true;
    newB_cell->myBCR.pMut = p.par[pmutAfterStartMut];
    newB_cell->isResponsive2CXCL13 = false;
    l.putcellat(newB_cell);
    newB_cell->event<<newB_cell->ID<<","<<0.0<<","<<newB_cell->MID<<","<<newB_cell->m_clonal_id<<","<<newB_cell->clonal_id<<","; // 0.0 is time of creation
    ListB_cell.push_back(newB_cell);
  }
  cerr << p.par[InitialNumberCB] << " CB generated" << endl;

  /// Initialize T follicular helper cells, mobile when they are not interacting with B cells, don't divide*/
  for (unsigned int j = 0; j < p.par[InitialNumberTC]; j = j + 1) {
    T_cell* newT_cell = new T_cell(p);
    newT_cell->cell_state = TC_free;
    newT_cell->cell_type = TFHC;
    newT_cell->can_move = true;
    newT_cell->position = l.getFreePosition(0, p.par[zoneRatioGC]);
    newT_cell->polarity = l.get_random_direction();
    newT_cell->getNewPersistentTime(p);
    newT_cell->persistence_time = p.par[Tcell_tp];
    newT_cell->speed = p.par[Tcell_speed];
    l.putcellat(newT_cell);
    ListT_cell.push_back(newT_cell);
  }
  cerr << p.par[InitialNumberTC] << " TC generated" << endl;
}

/// Simulation function
/// @param l is lattice
/// @param p is parameter list
void simulation::simulate(lattice& l, parameters& p) {
  double now_wall = get_wall_time();
  double now_cpu = get_cpu_time();
  // Create simulation output folder
  currentOutput->createFolder(outputFolder, p);
  // Cell initiation
  InitialCells(l, p);
    
//Record FDC_CXCL_AG distribution
static bool fdc_cxcl_data=false;
    if (fdc_cxcl_data) {
         fdc_cxcl_data = false;
        // do the initialization part
         FILE *FDC_CXCL;
         string folder1 = outputFolder + "FDC_CXCL_seed_" + to_string(Sim_Seed)+".csv";
         FDC_CXCL=fopen(folder1.c_str(), "w");
         fprintf(FDC_CXCL,"%s\n","Index,x1,x2,x3,type,cxcl12,cxcl13,Cag");

            for (int i=0; i<l.grid.size();i++)
            {
                for (int j=0; j<l.grid.size();j++)
                {
                    for (int k=0; k<l.grid.size();k++)
                    {
                        if (l.insideBorders(vector3D(i,j,k))){
                            fprintf(FDC_CXCL,"%d,%d,%d,%d,%d,%.16G,%.16G,%f\n",0,i,j,k,l.celltypeat(i, j, k),l.chemoat(CXCL12, i, j, k),l.chemoat(CXCL13, i, j, k),l.getAgat(i, j, k));

                        }
                        else {
                            fprintf(FDC_CXCL,"%d,%d,%d,%d,%d,%.16G,%.16G,%f\n",0,i,j,k,-1,l.chemoat(CXCL12, i, j, k),l.chemoat(CXCL13, i, j, k),l.getAgat(i, j, k));
                        }
                    }
                }
            }
        fclose(FDC_CXCL);
}
    ////
    
  int Total_time_steps = p.par[tmax] / p.par[dt];
  cout << "\nSimulation time= " << p.par[tmax] << endl;
  // Time loop
  for (int counter = 0; counter <= Total_time_steps; counter++)
  {
    // Time:
    double t = double(counter) * p.par[dt];  // in hours
    double recording_time_period = 2.;  // in hour
    int recording_time_steps = double(recording_time_period / p.par[dt]);
    if (not(pause1)) {
      // Record Output every fix ouhr
      if (fmod(counter, recording_time_steps) < p.par[dt]) {
        currentOutput->record_output_time_step(t, *this, p);
      }
        if (fmod(t, 24.0) ==0.0)
            cout << "time= " << int(t) << " Conecpt_ID="<< Concept_ID<<" SIM_id= " << Simulation_ID <<" seed="<< Sim_Seed<< endl;
        
        /// Redo the movment for thoes cells which can not move due to cell trafficking
      vector<vector3D> redo_list;
      redo_list.reserve(12000);
    /// Vector that stores the dead cells to remove later
      going_to_delet.reserve(10000);
      /// Shuffle the cell lists
      if (ListB_cell.size() > 0)
          std::shuffle(ListB_cell.begin(), ListB_cell.end(),generator);
      if (ListT_cell.size() > 0)
        std::shuffle(ListT_cell.begin(), ListT_cell.end(),generator);
        if (ListP_cell.size()>0)
        std::shuffle (ListP_cell.begin(),ListP_cell.end(),generator );
//        if (ListM_cell.size()>0)
      //  std::random_shuffle ( ListM_cell.begin(), ListM_cell.end() );
      // Calculations for Output cells at each time step
      Calc_Out(t, p, l, redo_list);
      // Calculations for T cells at each time step
      Calc_TC(p, l, redo_list);
      // Calculations for B cells at each time step
      Calc_BC(t, p, l, redo_list, going_to_delet);
      // Transfer newly differentiated Plasma cells from B cell list to avoid
      // interference in output files.
      transfer_plasma_from_Bcell_list(t, p, l, redo_list);
      // Influx of B cells to GC as an option
      // BCinflux(t,p,l);
      // Display simulation
      // Visualise(t,p);
      /// Redo move
      bool allow_exchange = true;
      if (allow_exchange) {
        if (redo_list.size() > 0) {
          redo_move(redo_list, l,t);
        }
      }
      // Remove dead cells
      clean_dead_cells(l);
    } else {
      counter--;
      t = double(counter) * p.par[dt];  // in hour
        //Visualise(t,p);
    }
  }
    
  cerr << "Simulation " << Sim_Seed <<" finished" << endl;
    ///    write Events
  for (unsigned int i = 0; i < ListB_cell.size(); i++) {
    B_cell* Bcell = ListB_cell.at(i);
    currentOutput->close_event(Bcell, sim_output,p.par[tmax] + 2,Simulation_ID,Concept_ID,Sim_Seed);
    currentOutput->write_event(Bcell, sim_output);
  }
  cerr << "Writing Events" << endl;
  currentOutput->write_event_2file(sim_output);
  cerr << "Writing Output files" << endl;
  currentOutput->Plasma_output(p.par[tmax] + 1, *this, p);
  double then_wall = get_wall_time();
  double then_cpu = get_cpu_time();

  cout << "Wall-time: " << then_wall - now_wall
       << " CPU-time: " << then_cpu - now_cpu << " seed: " << Sim_Seed <<endl<< endl;
}


/// Influx of B cells into the GC B cells influx to GC by a probability that can change with time, they find a random position in the whole GC to enter
/// @param time time
/// @param p parameters
/// @param l lattice

void simulation::BCinflux(double time, parameters& p, lattice& l) {
    
  double pBCinflux = double((p.par[rateCBinflow] * p.par[dt])) /
                     double((1.0 + exp((time - p.par[timeStopCBinflow]) /
                                       p.par[smoothnessStopCBinflow])));
  if (random::randomDouble(1.) < pBCinflux) {
    B_cell* Bcell = new B_cell(p);
    // Initialize
    Bcell->cell_state = founder;
    Bcell->cell_type = Centroblast;
    Bcell->persistence_time = p.par[Bcell_tp];
    Bcell->speed = p.par[Bcell_speed];
    Bcell->can_move = true;
    Bcell->setMyAffinity(p);
    Bcell->cyclestate = cycle_G1;
    Bcell->time_of_cycle_state_switch =
        random::cell_cycle_time(p.par[c_G1], cycle_G1);
    Bcell->position = l.getFreePosition(0, 1);  // Take free position
    if (not(l.insideBorders(Bcell->position))) {
      cerr << "Influx of B-cell at border position: " << Bcell->printcell()
           << endl;
      exit(1);
    }
    Bcell->polarity = l.get_random_direction();
    Bcell->nDivisions2do = p.par[nDivinflow];
    Bcell->getNewPersistentTime(p);
    Bcell->myBCR.pMut = p.par[pmutAfterStartMut];
    Bcell->isResponsive2CXCL12 = true;
    Bcell->isResponsive2CXCL13 = false;
    Bcell->event<<Bcell->ID<<","<<time<<","<<Bcell->MID<<","<<Bcell->m_clonal_id<<","<<Bcell->clonal_id <<",";
    l.putcellat(Bcell);
    ListB_cell.push_back(Bcell);
  }
}

/// Calculation of T cells
void simulation::Calc_TC(parameters& p, lattice& l,
                         vector<vector3D>& redo_list) {
  int N_T_cell = int(ListT_cell.size());
  for (int i = 0; i < N_T_cell; i++) {
    T_cell* Tcell = (T_cell*) ListT_cell.at(i);
    switch (Tcell->cell_state) {
      case TC_free: {
        Tcell->can_move = true;
        Tcell->move(p, l, redo_list);
        break;
      }
      case TC_connected: {
        Tcell->can_move = false;
        // Tcells in contact to CCs do not move!
        sort(Tcell->interactingCC.begin(), Tcell->interactingCC.end(),
             [](const B_cell* x, const B_cell* y) {
               return (x->retained_Ag > y->retained_Ag);
             });  
        Tcell->polarity.X =
            double(Tcell->interactingCC[0]->position.X - Tcell->position.X);
        Tcell->polarity.Y =
            double(Tcell->interactingCC[0]->position.Y - Tcell->position.Y);
        Tcell->polarity.Z =
            double(Tcell->interactingCC[0]->position.Z - Tcell->position.Z);
        break;
      }
    }
  }
}

// Visualization function
// void simulation::Visualise(double t, parameters &p)
//{
//
//      if (t<=0)
//        {
////            glm::vec3 cameraPosition(10.0f, 20.0f, 10.0f+ cameraDistance);
////            gluLookAt(100.0,100.0,0.0, 0.0, 0.0, 0.0, 0.0, 0.0, 0.0);
//            nextToDisplay(&ListB_cell, &ListT_cell, &ListFDC, &ListP_cell,
//            NULL, &ListSC, currentLattice, t);
//            display();
//            glutMainLoopEvent();
//
//        }
//    else if (t>0 && t<=720)
//        {
////            gluPerspective(60.0, 1.0, 1, 500); // Note : deph test works
///only if the first plane is > 0
//
////            gluLookAt(10, 10, 10,  0, 0, 0, 0, 1, 0);
//
//            if(fmod(t+1e-9,1)< p.par[dt])
//                {
//                    nextToDisplay(&ListB_cell, &ListT_cell, &ListFDC,
//                    &ListP_cell, NULL, &ListSC, currentLattice, t);
//                    display();
//                     glutMainLoopEvent();
//
//                }
//        }
//    else if (t>720 && t< 7200)
//        {
//            nextToDisplay(&ListB_cell, &ListT_cell, &ListFDC, &ListP_cell,
//            NULL, &ListSC, currentLattice, t);
//            display();
//             glutMainLoopEvent();
//        }
//    else if (t>=7200)
//        {
//            if(fmod(t+1e-9,1)< p.par[dt])
//                {
//                    nextToDisplay(&ListB_cell, &ListT_cell, &ListFDC,
//                    &ListP_cell, NULL, &ListSC, currentLattice, t);
//                    display();
//                     glutMainLoopEvent();
//                }
//        }
//}

/// Calculation of B cells
void simulation::Calc_BC(double t, parameters& p, lattice& l,
                         vector<vector3D>& redo_list,
                         vector<int>& going_to_delet) {
    
 int N_B_cell= int (ListB_cell.size());
 for ( int i = 0; i < N_B_cell; i++) {
    B_cell* Bcell = ListB_cell.at(i);
      
///Centrocytes_________________________________
    if (Bcell->cell_type == Centrocyte) {
      switch (Bcell->cell_state) {
        case unselected: {
          Bcell->clock += 1;
          Bcell->Resensitize2Chemokines(p, l);
          Bcell->BC_FDC_interaction_clock += p.par[dt];
            /// CCs finished FDC interaction windows
          if (Bcell->BC_FDC_interaction_clock > p.par[collectionFDCperiod]) {
            if (Bcell->retained_Ag > 20) {
              Bcell->retained_Ag = 20;
            }
            if (Bcell->retained_Ag <= 0) {
              Bcell->cell_state = apoptosis;
              Bcell->isResponsive2CXCL12 = false;
              Bcell->isResponsive2CXCL13 = true;
            } else {
              /// CCs selected by FDCs here
              Bcell->cell_state = FDC_selected;
              Bcell->Selected_by_FDC = true;
              Bcell->can_move = true;
            }
              static int tmp_count=0;
              if (tmp_count==0)
              {
              FILE *Ag_collection_events;
              string folder2 = outputFolder + "Ag_collection_events_seed_" + to_string(Sim_Seed)+".csv"  ;
                 
                 Ag_collection_events = fopen(folder2.c_str(), "w");
                 fprintf(Ag_collection_events, "%s","Born_time,ID,MID,M_clone_id,clone_id,States,Distances,Affinity,delta_Affinity,P_on,delta_Pon,P_off,delta_Poff,m_theta,delta_mtheta,N_Ags,N_divisions,N_Mutations_beneficial,N_Mutations_harmful,FDC_interaction_nums,FDC_binding_nums,FDC_binding_no_Ag,FDC_interaction_time_total,FDC_Selected,SIM_ID,Concept_ID,Sim_seed\n");
                  fclose(Ag_collection_events);
                  tmp_count++;
              }
              FILE *Ag_collection_events;
                           string folder2 = outputFolder + "Ag_collection_events_seed_" + to_string(Sim_Seed)+".csv"  ;
                              Ag_collection_events = fopen(folder2.c_str(), "a");
              fprintf(Ag_collection_events, "%f,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d,%d,%d,%d,%f,%d,%d,%d,%d\n",Bcell->Born_time,Bcell->ID,Bcell->MID,Bcell->m_clonal_id,Bcell->clonal_id,Bcell->cell_state,Bcell->myBCR.distance,Bcell->myBCR.getMyAffinity4Ag(p),Bcell->myBCR.delta_Affinity,Bcell->myBCR.P_on,Bcell->myBCR.delta_Pon,Bcell->myBCR.P_off,Bcell->myBCR.delta_Poff,Bcell->myBCR.m_theta,Bcell->myBCR.delta_mtheta,Bcell->retained_Ag,Bcell->nDivisions2do,Bcell->myBCR.nMutFromGermline.beneficial,Bcell->myBCR.nMutFromGermline.harmful,Bcell->nFDCcontacts,Bcell->fdc_binding_nums,Bcell->fdc_binding_nums_no_ag,Bcell->fdc_interaction_time_history,Bcell->Selected_by_FDC,Simulation_ID,Concept_ID,Sim_Seed);
                               fclose(Ag_collection_events);
              
          }
            /// CCs are in  FDC interaction window
          else if (Bcell->clock > p.par[testDelay]) {
            vector3D fdc_position(-1, -1, -1);
            vector<vector3D> neighbours = l.getNeighbour_nn(Bcell->position);
            for (unsigned int j = 0; j < neighbours.size(); j++) {
              if (l.insideBorders(neighbours[j])) {
                if (l.getAgat(neighbours[j]) > 0.) {
                  fdc_position = neighbours[j];
                }
              }
            }
              ///Binding CC to FDC
            if (fdc_position.X != -1) {
              bool suppress_next_interaction = false;
              Bcell->myBCR.calculate_on_off_P(p);
              Bcell->setMyAffinity(p);
              double binding_probability = Bcell->myBCR.P_on;
              Bcell->nFDCcontacts += 1;
              if (random::randomDouble(1) <= binding_probability)
                  
              {
                Bcell->fdc_pos = fdc_position;
                short success;
                if ((l.getAgat(fdc_position) >= 1.) && (random::randomDouble(1) < (double)(l.getAgat(fdc_position) /p.par[agSaturation])))  
                {
                  success = 1;
                } else {
                  success = 0;
                }
                if (success == 1) {
                  Bcell->cell_state = contact_FDC;
                  l.removeAgAt(Bcell->fdc_pos);
                  Bcell->can_move = false;
                  Bcell->fdc_binding_nums += 1;
                } else {
                  suppress_next_interaction = true;
                    Bcell->fdc_binding_nums_no_ag += 1;
                }
              } else {
                suppress_next_interaction = true;
              }
              if (suppress_next_interaction) {
                Bcell->clock = 0;
              }
            }
          }
          if (Bcell->cell_state == unselected) {
            Bcell->move(p, l, redo_list);  
          }
          break;
        }
              
        case contact_FDC: {
          Bcell->BC_FDC_interaction_clock += p.par[dt];
          Bcell->fdc_interaction_time_history += p.par[dt];
            Bcell->myBCR.calculate_on_off_P(p);
            long double Ag_collcetion_prob = Bcell->myBCR.P_off ;

///Antigen collection, SC-1: With interruptions.
           if (not(Bcell->myBCR.checked_initial_dissociation)){
            	if (random::randomDouble(1.0)<Ag_collcetion_prob)
                    {
                     	Bcell->myBCR.checked_initial_dissociation=true;
                    }
            	else {
                	l.returnAg(Bcell->fdc_pos);
               		Bcell->cell_state = unselected;
                	Bcell->can_move = true;
                	Bcell->clock = 0;
                      }
            }
           
	 if (Bcell->myBCR.checked_initial_dissociation)
            {
             	if (random::randomDouble(1.0)<Ag_collcetion_prob)
                {
                    if (random::randomDouble(1.0) <= 0.04)
                		{
                   		Bcell->retained_Ag += 1.0;
                    		Bcell->cell_state = unselected;
                    		Bcell->can_move = true;
                    		Bcell->clock = 0;
                    		Bcell->myBCR.checked_initial_dissociation=false;
                		}
                }
                else {
               		l.returnAg(Bcell->fdc_pos);
                	Bcell->cell_state = unselected;
               		Bcell->can_move = true;
                	Bcell->clock = 0;
                    Bcell->myBCR.checked_initial_dissociation=false;
                }
            }

          Bcell->isResponsive2CXCL13 = false;
          break;
        }

        case FDC_selected: {
          Bcell->Resensitize2Chemokines(p, l);
          vector3D tc_position(-1, -1, -1);
          vector<vector3D> tmp_neighbours = l.getNeighbour_nn(Bcell->position);
          vector<vector3D> TC_neighbours;
          for (unsigned int k = 0; k < tmp_neighbours.size(); k++) {
            if (l.insideBorders(tmp_neighbours[k])) {
              if (l.celltypeat(tmp_neighbours[k]) == TFHC) {
                TC_neighbours.push_back(tmp_neighbours[k]);
              }
            }
          }
          if (TC_neighbours.size() > 0) {
            short x = random::randomInteger(int(TC_neighbours.size()));
            tc_position = TC_neighbours[x];
          }
          if (tc_position.X == -1) {
            Bcell->move(p, l, redo_list);
          } else {
            /// bind CC to TC
            Bcell->cell_state = contact_TC;
            Bcell->can_move = false;
            Bcell->Bc_Tc_interaction_clock = 0;
            T_cell* TC = (T_cell*)l.cellat(tc_position);
            Bcell->interactingTC = TC;
            TC->nIncontactCCs += 1;
            TC->cell_state = TC_connected;
            TC->interactingCC.push_back(Bcell);
          }
          break;
        }
        case contact_TC: {
          Bcell->isResponsive2CXCL13 = false;
          // duration of contact
          Bcell->Tc_interaction_history.first += p.par[dt];
          Bcell->Bc_Tc_interaction_clock += p.par[dt];
          T_cell* TC = (T_cell*)Bcell->interactingTC;
          if (TC == NULL) {
            cout << "Accessing null Tcell from CC" << endl;
          } else if (TC->ID != Bcell->interactingTC->ID) {
            cout << "Accessing wrong TC from CC" << endl;
          } else {
            vector3D CC_neighbour = l.get_nn_directed2(TC);
            if (CC_neighbour.X != -1) {
              cell* cellthere = l.grid.at(CC_neighbour.X)
                                    .at(CC_neighbour.Y)
                                    .at(CC_neighbour.Z);
              if (cellthere != NULL) {
                B_cell* neighbour_CC = (B_cell*)l.cellat(CC_neighbour);
                if (neighbour_CC->ID == Bcell->ID) {
                  // TC and CC are face2face
                  Bcell->TCsignalDuration += p.par[dt];
                  // duration of signal
                  Bcell->Tc_interaction_history.second += p.par[dt];
                }
              } else {
                cout << "time= " << t
                     << " There is no CC in the directed position, TC pos="
                     << TC->position.print()
                     << " CC pos=" << Bcell->position.print()
                     << " Directedpos=" << CC_neighbour.print()
                     << " CCID=" << Bcell->ID << " TCID=" << TC->ID
                     << " cell_state=" << Bcell->cell_state << endl;
              }
            } else {
              //  cout<<"The CC neighbour is negative -1"<<endl;
            }
          }
          if (Bcell->TCsignalDuration > p.par[tcRescueTime]) {
            // CC_TC Selection
            if (Bcell->retained_Ag > 20) {
              Bcell->retained_Ag = 20;
            }
            Bcell->timeleft2recycle(p);
            double pMHC = Bcell->retained_Ag;
            double ag_factor = pow(pMHC, p.par[pMHCdepHill]);
            // Record selected CC mutation frequencies
            // Number of pmhc dpendent divisions
            Bcell->pMHC_dependent_number_of_divisions =
                p.par[pMHCdepMin] +
                (p.par[pMHCdepMax] - p.par[pMHCdepMin]) * ag_factor /
                    (ag_factor + pow(p.par[pMHCdepK], p.par[pMHCdepHill]));
            double ndivtmp =2.0;
            //                            if(Bcell->pMHC_dependent_number_of_divisions>=0)
            //                            {
            //                                ndivtmp=Bcell->pMHC_dependent_number_of_divisions;
            //                            }
            Bcell->nDivisions2do = int(ndivtmp);
            //                            ndivtmp -=
            //                            double(Bcell->nDivisions2do);
            //                            if (random::randomDouble(1) < ndivtmp)
            //                            {
            //                                ++Bcell->nDivisions2do;
            //                            }

            Bcell->cell_state = TC_selected;

            Bcell->Selected_by_TC = true;

            Bcell->can_move = true;
            Bcell->TC_selected_clock = 0.0;
          }
          if (not(Bcell->cell_state == TC_selected) &&
              Bcell->Bc_Tc_interaction_clock > p.par[tcTime]) {
            // apoptotic cells mutation frequency
            Bcell->cell_state = apoptosis;
            Bcell->can_move = true;
            Bcell->isResponsive2CXCL13 = true;
            Bcell->isResponsive2CXCL12 = false;
          }

          if (Bcell->cell_state == TC_selected ||
              Bcell->cell_state == apoptosis) {
            if (not(TC == NULL)) {
              TC->liberateCC_TC(Bcell);
            }
            Bcell->can_move = true;
          }
          break;
        }
        case TC_selected: {
          Bcell->isResponsive2CXCL13 = false;
          Bcell->TC_selected_clock += p.par[dt];
          if (Bcell->TC_selected_clock > Bcell->Recycling_delay) {
            if (random::randomDouble(1) < p.par[p_dif])
            {
              // Here cells recycle
              Bcell->isResponsive2CXCL12 = true;
              Bcell->isResponsive2CXCL13 = false;
              Bcell->cell_type = Centroblast;
              
              if (Bcell->retained_Ag > 0.) {
                Bcell->IamHighAg = true;
              }

             
              Bcell->cell_state = recycled;
              Bcell->setMyAffinity(p);
              Bcell->cyclestate = cycle_G1;
              Bcell->transmit_CCdelay2cycle(p);
              if (Bcell->nDivisions2do <= 0) {
                Bcell->cyclestate = cycle_G0;
              }
              Bcell->Recycling_delay = 0;
            } else {
              Bcell->move(p, l, redo_list);
            }
          } else {
            Bcell->move(p, l, redo_list);
          }
          break;
        }

        case apoptosis: {
          if (random::randomDouble(1) < p.par[macrophage]) {
            //#record_event
            currentOutput->close_event(Bcell, sim_output, t,Simulation_ID,Concept_ID,Sim_Seed);
            currentOutput->write_event(Bcell, sim_output);
            going_to_delet.push_back(Bcell->ID);
          } else {
            Bcell->isResponsive2CXCL12 = false;
            Bcell->isResponsive2CXCL13 = true;
            Bcell->Resensitize2Chemokines(p, l);
            Bcell->move(p, l, redo_list);
          }
          break;
        }
        default:
          break;
      }
    }
    //___________________________________________________________________________________________________________

    // Centroblasts_______________________________________________________________________________________________
    if (Bcell->cell_type == Centroblast) {
      // Increase cell cycle time
      Bcell->cycle_state_time += p.par[dt];
      // Resensitize
      Bcell->Resensitize2Chemokines(p, l);
      // Switch cycle state
      if ((Bcell->cycle_state_time >= Bcell->time_of_cycle_state_switch) &&
          (Bcell->cyclestate < cycle_Divide)) {
        Bcell->ContinueCellCycle(p);
        Bcell->cycle_state_time = 0.0;
      }
      if (Bcell->cyclestate == cycle_Divide) {
        Bcell->proliferate(p, l, t, ListB_cell, *currentOutput, *this);
      } else if (Bcell->cyclestate == cycle_G0) {
        if (random::randomDouble(1) <p.par[p_dif])  
        {
          if (Bcell->retained_Ag > 0. && Bcell->IamHighAg) {
            // Output cells differentiation
            Bcell->cell_type = Plasmacell;
            // record frequency of  output cells and their mutations
          }
          // centrocytes
          else {
            // CCs created here
            Bcell->isResponsive2CXCL12 = false;
            Bcell->isResponsive2CXCL13 = true;
            Bcell->cell_type = Centrocyte;
            Bcell->cell_state = unselected;
            
            Bcell->nFDCcontacts = 0;
            Bcell->retained_Ag = 0;
            Bcell->Selected_by_FDC = false;
            Bcell->Selected_by_TC = false;
            Bcell->clock = 0;  // not sure
            Bcell->interactingTC = NULL;
            Bcell->Bc_Tc_interaction_clock = 0.;
            Bcell->Recycling_delay = 0.;  // Time spent inside CBgoingLZ (time
                                          // for moving to Light Zone)
            Bcell->BC_FDC_interaction_clock =
                0.;  // Time since a B_cell became CC_free (in sec).
            Bcell->TC_selected_clock = 0.0;
            Bcell->Selected_by_FDC = false;
            Bcell->Selected_by_TC = false;
            Bcell->TCsignalDuration =
                0.;  // Acumulated signal from currently interacting TC (in
                     // sec).(As imput to ODE)
            Bcell->Tc_interaction_history.first = 0.0;
            Bcell->Tc_interaction_history.second = 0.0;
            Bcell->fdc_interaction_time_history = 0.0;
            Bcell->fdc_binding_nums_no_ag = 0.0;
            Bcell->fdc_binding_nums = 0.0;

          }
        }
      }
      if ((Bcell->cell_type == Centroblast) && (Bcell->cyclestate != cycle_M)) {
        Bcell->move(p, l, redo_list);
      }
    }
    //______________________________________________________________________________________
  }
}

void simulation::clean_dead_cells(lattice& l) {

  long tmp_size_1 = ListB_cell.size();
  for (unsigned int j = 0; j < going_to_delet.size(); j++) {
    if (ListB_cell.size() == 1) {
      if (going_to_delet[j] == ListB_cell.at(0)->ID) {
        l.removecellat(ListB_cell.at(0)->position);
        delete ListB_cell[0];
        ListB_cell.pop_back();
      }
      if (ListB_cell.size() > 0) {
        cout << "Error in deleting dead B cells (1)." << endl;
        exit(1);
      }
    } else if (ListB_cell.size() > 1) {
      for (unsigned int ks = 0; ks < ListB_cell.size(); ks++) {
        long tmp_size_2 = ListB_cell.size();
        if (ListB_cell.at(ks) != NULL) {
          if (ListB_cell.at(ks)->ID == going_to_delet[j]) {
            l.removecellat(ListB_cell.at(ks)->position);  // remove from lattice
            delete ListB_cell.at(ks);
            ListB_cell.at(ks) = NULL;
            tmp_size_2--;
          }
        }
      }
    }
  }

  ListB_cell.erase(remove_if(ListB_cell.begin(), ListB_cell.end(),
                             [](const B_cell* x) { return (x == NULL); }),
                   ListB_cell.end());
  if (tmp_size_1 != (going_to_delet.size() + ListB_cell.size())) {
    cout << "Error in deleting dead B cells (2)." << endl;
    exit(1);
  }
  going_to_delet.clear();
}

void simulation::Calc_Out(double t, parameters& p, lattice& l,
                          vector<vector3D>& redo_list) {
  for (unsigned int i = 0; i < ListP_cell.size(); i++) {
    Plasma_cell* Plasma = ListP_cell.at(i);
    if (not(Plasma->cell_state == Plasma_Out)) {
      if (l.is_at_border(Plasma->position)) {
        l.removecellat(Plasma->position);
        Plasma->cell_state = Plasma_Out;
      } else {
        Plasma->move(p, l, redo_list);
      }
    }
  }

  //
  //_______________________________________________________________________________________________________________
}
void simulation::transfer_plasma_from_Bcell_list(double t, parameters& p,
                                                 lattice& l,
                                                 vector<vector3D>& redo_list) {
  
  // Output
  // cells___________________________________________________________________________________________________

  for (unsigned int i = 0; i < ListB_cell.size(); i++) {
    B_cell* Bcell = ListB_cell.at(i);
    if (Bcell->cell_type == Plasmacell) {
      Plasma_cell* new_Plasma = new Plasma_cell(p, Bcell);
      new_Plasma->birth_time = t;
      new_Plasma->ID = Bcell->ID;
      new_Plasma->MID = Bcell->MID;
      l.removecellat(Bcell->position);
      l.putcellat(new_Plasma);
      ListP_cell.push_back(new_Plasma);
      delete ListB_cell.at(i);
      ListB_cell.at(i) = NULL;
    }
  }

  ListB_cell.erase(remove_if(ListB_cell.begin(), ListB_cell.end(),
                             [](const B_cell* x) { return (x == NULL); }),
                   ListB_cell.end());

  for (unsigned int i = 0; i < ListB_cell.size(); i++) {
    B_cell* Bcell = ListB_cell.at(i);

    if (Bcell == NULL) {
      cout << "Erorr, NULL BC" << endl;
    }
  }

  for (unsigned int i = 0; i < ListP_cell.size(); i++) {
    Plasma_cell* Plasma = ListP_cell.at(i);

    if (Plasma == NULL) {
      cout << "Erorr, NULL PC" << endl;
    }
  }
}
