#include "output.h"
#include <set>
#include <numeric>
#include "mafalda.h"
#include "cell.h"
#include <sys/types.h>
#include <sys/stat.h>
#ifdef __linux__
#include <sys/types.h>
#include <sys/stat.h>
#endif
using namespace std;

// ministats
    vector<ministat> Bcell_counts;  // number of Bcells
    ministat Bcell_affinity;
    ministat Bcell_distances;
    ministat Bcell_P_on;
    ministat Bcell_P_off;
    ministat Bcell_Ag;
    ministat B_cell_m_theta;

    vector<ministat>Clon_Pon;
    vector<ministat>Clon_Poff;
    vector<ministat>Clon_Affinity;
    vector<ministat>Clon_Distance;
    vector<ministat>Clon_counts;
    vector<ministat>Clon_Ag;
    vector<ministat>Clone_m_theta;

    vector<ministat>Outcell_Pon;
    vector<ministat>Outcell_Poff;
    vector<ministat>Outcell_Affinity;
    vector<ministat>Outcell_Distance;
    vector<ministat>Outcell_counts;
    vector<ministat>Outcell_m_theta;

    ministat total_outcell_affinity;
    ministat total_outcell_distances;
    ministat total_outcell_P_on;
    ministat total_outcell_P_off;
    ministat total_outcell_count;
    ministat total_outcell_m_theta;

//Initialize ministats
    output::output(string _parfname) : parfname(_parfname) { initialize_fileds(); }
    output::~output() {}

// Create folder
void output::createFolder(string folderName, parameters &p)  {

        output_path = folderName;
        struct stat info;
        if( stat( output_path.c_str(), &info ) != 0 )
        printf( "Output folder: cannot access %s\n", output_path.c_str() );
        else if( info.st_mode & S_IFDIR )
        printf( "Output folder: already exists %s\n", output_path.c_str() );
        else
        {
        const char *c = output_path.c_str();
        #ifdef __linux__
        mkdir(c, S_IRWXU | S_IRWXG | S_IROTH | S_IXOTH);
        #endif
        #ifdef __APPLE__
        char cstr[output_path.size() + 1];
        strcpy(cstr, output_path.c_str());
        mkdir(cstr, 0777);
        #endif
        }
   
// Define files and Headers

    //B-cell-time.csv
    FILE *Bcell_time;
    string folder1 = output_path + "Bcell_time_seed_" + to_string(Output_ID)+".csv" ;
    Bcell_time = fopen(folder1.c_str(), "w");
    fprintf(Bcell_time, "%s","time,BCs,CBs,CCs,founder,unselected,contact_FDC,FDC_Selected,contact_TC,TC_Selected,recycled,apoptosis,distances_avg,distances_sd,affinity_avg,affinity_sd,Pon_avg,Pon_sd,Poff_avg,Poff_sd,m_theta_avg,m_theta_sd,N_Ag_avg,N_Ag_sd,total_consumed_Ag,SIM_ID,Concept_ID,Sim_seed\n");
    fclose(Bcell_time);
    
    //Outcell-time.csv
    FILE *Outcell_time;
    string folder0 = output_path + "Outcell_time_seed_" + to_string(Output_ID)+".csv" ;
    Outcell_time = fopen(folder0.c_str(), "w");
    fprintf(Outcell_time, "%s","time,N_clone1,N_clone2,N_clone3,clone1_Distance_avg,clone1_Distance_sd,clone2_Distance_avg,clone2_Distance_sd,clone3_Distance_avg,clone3_Distance_sd,clone1_Affinity_avg,clone1_Affinity_sd,clone2_Affinity_avg,clone2_Affinity_sd,clone3_Affinity_avg,clone3_Affinity_sd,clone1_Pon_avg,clone1_Pon_sd,clone2_Pon_avg,clone2_Pon_sd,clone3_Pon_avg,clone3_Pon_sd,clone1_Poff_avg,clone1_Poff_sd,clone2_Poff_avg,clone2_Poff_sd,clone3_Poff_avg,clone3_Poff_sd,clone1_m_theta_avg,clone1_m_theta_sd,clone2_m_theta_avg,clone2_m_theta_sd,clone3_m_theta_avg,clone3_m_theta_sd,total_output,total_Distance_avg,total_Distance_sd,total_Affinity_avg,total_Affinity_sd,total_Pon_avg,total_Pon_sd,total_Poff_avg,total_Poff_sd,total_m_theta_avg,total_m_theta_sd,SIM_ID,Concept_ID,Sim_seed\n");
    fclose(Outcell_time);
    
    for (int i=0;i<3;i++){
    //Clone_X_time.csv
    FILE *Clone_X;
    string folder3 = output_path + "Clone_"+to_string(i+1)+"_seed_" + to_string(Output_ID)+".csv" ;
    Clone_X = fopen(folder3.c_str(), "w");
    fprintf(Clone_X, "%s","time,Total,founder,unselected,contact_FDC,FDC_Selected,contact_TC,TC_Selected,recycled,apoptosis,Distance_avg,Distance_sd,Affinity_avg,Affinity_sd,Pon_avg,Pon_sd,Poff_avg,Poff_sd,m_theta_avg,m_theta_sd,N_Ag_avg,N_Ag_sd,clone_consumed_Ag,SIM_ID,Concept_ID,Sim_seed\n");
    fclose(Clone_X);
    }
    
    //Event-time.csv
    FILE *Event_time;
    string folder2 = output_path + "Event_data_seed_" + to_string(Output_ID)+".csv" ;
    Event_time = fopen(folder2.c_str(), "w");
    fprintf(Event_time, "%s","ID,Born_time,MID,M_clone_id,D_clone_id,States,Distances,Affinity,delta_Affinity,P_on,delta_Pon,P_off,delta_Poff,m_theta,delta_mtheta,N_Ag,N_divisions,N_Mutations_beneficial,N_Mutations_harmful,FDC_interaction_nums,FDC_binding_nums,FDC_binding_nums_no_Ag,FDC_interaction_time_total,TC_interaction_time,TC_signaling_time,Selected_by_FDC,Selected_by_TC,Leave_time,SIM_ID,Concept_ID,Sim_seed\n");
    fclose(Event_time);
     }

void output::initialize_fileds() {

// 0-9 -> CC states (7) + CB (2)

for (int i = 0; i < 8; i++) {
    Bcell_counts.push_back(ministat());  // number of Bcells
  }
    
    for (int i = 0; i < 3; i++) {
         Clon_Pon.push_back(ministat());
         Clon_Poff.push_back(ministat());
         Clon_Affinity.push_back(ministat());
         Clon_Distance.push_back(ministat());
         Clon_counts.push_back(ministat());
         Clon_Ag.push_back(ministat());
         Clone_m_theta.push_back(ministat());
         Outcell_Pon.push_back(ministat());
         Outcell_Poff.push_back(ministat());
         Outcell_Affinity.push_back(ministat());
         Outcell_Distance.push_back(ministat());
         Outcell_counts.push_back(ministat());  // number of Bcells
         Outcell_m_theta.push_back(ministat());
    }
}

void output::clear_fileds() {

    for (int i = 0; i < 8; i++) {
      Bcell_counts.at(i).clear_ministat();   //.at(i).clear_ministat();   //number of Bcells
        }
      Bcell_affinity.clear_ministat();
      Bcell_distances.clear_ministat();
      Bcell_P_on.clear_ministat();
      Bcell_P_off.clear_ministat();
      Bcell_affinity.clear_ministat();
      Bcell_distances.clear_ministat();
      Bcell_Ag.clear_ministat();
    
      for (int i = 0; i < 3; i++) {
      Clon_Pon.at(i).clear_ministat();
      Clon_Poff.at(i).clear_ministat();
      Clon_Affinity.at(i).clear_ministat();
      Clon_Distance.at(i).clear_ministat();
      Clon_counts.at(i).clear_ministat();
      Clon_Ag.at(i).clear_ministat();
      Clone_m_theta.at(i).clear_ministat();
      Outcell_Pon.at(i).clear_ministat();
      Outcell_Poff.at(i).clear_ministat();
      Outcell_Affinity.at(i).clear_ministat();
      Outcell_Distance.at(i).clear_ministat();
      Outcell_counts.at(i).clear_ministat();  // number of Bcells
      Outcell_m_theta.at(i).clear_ministat();
        }

     total_outcell_affinity.clear_ministat();
     total_outcell_distances.clear_ministat();
     total_outcell_P_on.clear_ministat();
     total_outcell_P_off.clear_ministat();
     total_outcell_count.clear_ministat();
     total_outcell_m_theta.clear_ministat();
}

//Recording temporal data from simulation
// Take fields from simulation into master observer variable to create file.
void output::record_output_time_step(double currentTime, simulation &currentSim,
                                     parameters &p) {
    
// Bcell data
    for (int i = 0; i < currentSim.ListB_cell.size(); i++) {
    B_cell *Bcell = currentSim.ListB_cell.at(i);
    if (Bcell->cell_state > 7){
      cout << "Error, wrong cell in BC list," << Bcell->cell_state << endl;
      exit(1);
    }
    Bcell_counts[Bcell->cell_state].add(1);
    Bcell->setMyAffinity(p);
    Bcell->myBCR.calculate_on_off_P(p);
    Bcell_affinity.add(Bcell->MyAffinity);
    Bcell_distances.add(Bcell->myBCR.distance);
    Bcell_P_on.add(Bcell->myBCR.P_on);
    Bcell_P_off.add(Bcell->myBCR.P_off);
    Bcell_Ag.add(Bcell->retained_Ag);
    B_cell_m_theta.add(Bcell->myBCR.m_theta);
   }
    //All BCs CBs and CCs
    double BC_counts=0;
    for (int i=0; i<8;i++)
         BC_counts += Bcell_counts[i].sum;
    //0 is founder state and 6 is recycled
    double CBs_counts = Bcell_counts[0].sum+Bcell_counts[6].sum;
    // 1,2,3,4,5,7 are states of CCs
    //Without apoptosis
    double CCs_counts = Bcell_counts[1].sum+Bcell_counts[2].sum+Bcell_counts[3].sum+Bcell_counts[4].sum+Bcell_counts[5].sum;
    //Write in file
    FILE *Bcell_time_data;
    string folder1 = output_path + "Bcell_time_seed_" + to_string(Output_ID)+".csv" ;
    Bcell_time_data = fopen(folder1.c_str(), "a");
    fprintf(Bcell_time_data, "%f,%f,%f,%f,", currentTime,BC_counts,CBs_counts,CCs_counts);
    for (int i = 0; i < 8; i++) {
    fprintf(Bcell_time_data, "%f,", Bcell_counts[i].sum);
    }
    fprintf(Bcell_time_data, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d\n",Bcell_distances.average(),Bcell_distances.stddev(),Bcell_affinity.average(),Bcell_affinity.stddev(),Bcell_P_on.average(),Bcell_P_on.stddev(),Bcell_P_off.average(),Bcell_P_off.stddev(),B_cell_m_theta.average(),B_cell_m_theta.stddev(),Bcell_Ag.average(),Bcell_Ag.stddev(),currentSim.currentLattice->TotalAmountAginLattice-Bcell_Ag.sum,currentSim.Simulation_ID,currentSim.Concept_ID,currentSim.Sim_Seed);
    fclose(Bcell_time_data);
    
        
//Output cells
    
    //Calculations:
    for (int i = 0; i < currentSim.ListP_cell.size(); i++) {
    Plasma_cell *Pcell = currentSim.ListP_cell.at(i);
    Pcell->myBCR.calculate_on_off_P(p);
    Outcell_Pon.at(Pcell->m_clonal_id).add(Pcell->myBCR.P_on);
    Outcell_Poff.at(Pcell->m_clonal_id).add(Pcell->myBCR.P_off);
    Outcell_Affinity.at(Pcell->m_clonal_id).add(Pcell->myBCR.getMyAffinity4Ag(p));
    Outcell_Distance.at(Pcell->m_clonal_id).add(Pcell->myBCR.distance);
    Outcell_counts.at(Pcell->m_clonal_id).add(1);
    Outcell_m_theta.at(Pcell->m_clonal_id).add(Pcell->myBCR.m_theta);
    total_outcell_affinity.add(Pcell->myBCR.getMyAffinity4Ag(p));
    total_outcell_distances.add(Pcell->myBCR.distance);
    total_outcell_P_on.add(Pcell->myBCR.P_on);
    total_outcell_P_off.add(Pcell->myBCR.P_off);
    total_outcell_m_theta.add(Pcell->myBCR.m_theta);
    total_outcell_count.add(1);
    }
    //Write to file
    FILE *Outcell_time_data;
    string folder0 = output_path + "Outcell_time_seed_" + to_string(Output_ID)+".csv" ;
    Outcell_time_data = fopen(folder0.c_str(), "a");
    fprintf(Outcell_time_data, "%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d\n",currentTime,Outcell_counts.at(0).sum,Outcell_counts.at(1).sum,Outcell_counts.at(2).sum,
            Outcell_Distance.at(0).average(),Outcell_Distance.at(0).stddev(),Outcell_Distance.at(1).average(),Outcell_Distance.at(1).stddev(),Outcell_Distance.at(2).average(),Outcell_Distance.at(2).stddev(),
            Outcell_Affinity.at(0).average(),Outcell_Affinity.at(0).stddev(),Outcell_Affinity.at(1).average(),Outcell_Affinity.at(1).stddev(),Outcell_Affinity.at(2).average(),Outcell_Affinity.at(2).stddev(),
            Outcell_Pon.at(0).average(),Outcell_Pon.at(0).stddev(),Outcell_Pon.at(1).average(),Outcell_Pon.at(1).stddev(),Outcell_Pon.at(2).average(),Outcell_Pon.at(2).stddev(),
            Outcell_Poff.at(0).average(),Outcell_Poff.at(0).stddev(),Outcell_Poff.at(1).average(),Outcell_Poff.at(1).stddev(),Outcell_Poff.at(2).average(),Outcell_Poff.at(2).stddev(),Outcell_m_theta.at(0).average(),Outcell_m_theta.at(0).stddev(),Outcell_m_theta.at(1).average(),Outcell_m_theta.at(1).stddev(),Outcell_m_theta.at(2).average(),Outcell_m_theta.at(2).stddev(),
            total_outcell_count.sum,total_outcell_distances.average(),total_outcell_distances.stddev(),total_outcell_affinity.average(),total_outcell_affinity.stddev(),total_outcell_P_on.average(),total_outcell_P_on.stddev(),total_outcell_P_off.average(),total_outcell_P_off.stddev(),total_outcell_m_theta.average(),total_outcell_m_theta.stddev(),currentSim.Simulation_ID,currentSim.Concept_ID,currentSim.Sim_Seed);
       fclose(Outcell_time_data);
       
//Clones
        //Calculations:
        //initialize vector of CC states for clones 
        vector<vector<int> > Clone_State_counts( 3 , vector<int> (8, 0));
        for (int i = 0; i < currentSim.ListB_cell.size(); i++) {
        B_cell *Clone_cell = currentSim.ListB_cell.at(i);
        Clone_cell->myBCR.calculate_on_off_P(p);
        Clone_State_counts.at(Clone_cell->m_clonal_id).at(Clone_cell->cell_state) +=1;
          if (Clone_cell->cell_state!=apoptosis){
	Clon_Pon.at(Clone_cell->m_clonal_id).add(Clone_cell->myBCR.P_on);
        Clon_Poff.at(Clone_cell->m_clonal_id).add(Clone_cell->myBCR.P_off);
        Clon_Affinity.at(Clone_cell->m_clonal_id).add(Clone_cell->myBCR.getMyAffinity4Ag(p));
        Clon_Distance.at(Clone_cell->m_clonal_id).add(Clone_cell->myBCR.distance);
        Clon_counts.at(Clone_cell->m_clonal_id).add(1);
        Clon_Ag.at(Clone_cell->m_clonal_id).add(Clone_cell->retained_Ag);
        Clone_m_theta.at(Clone_cell->m_clonal_id).add(Clone_cell->myBCR.m_theta);
             }
        }
    
        //Write to file
        for (int i=0; i<3;i++)
        {
          FILE *Clone_X;
          string folderX = output_path + "Clone_"+to_string(i+1)+"_seed_" + to_string(Output_ID)+".csv" ; 
          Clone_X = fopen(folderX.c_str(), "a");
            fprintf(Clone_X, "%f,%f,%d,%d,%d,%d,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d\n",currentTime,Clon_counts.at(i).sum,Clone_State_counts.at(i).at(0),Clone_State_counts.at(i).at(1),Clone_State_counts.at(i).at(2),Clone_State_counts.at(i).at(3),Clone_State_counts.at(i).at(4),Clone_State_counts.at(i).at(5),Clone_State_counts.at(i).at(6),Clone_State_counts.at(i).at(7),Clon_Distance.at(i).average(),Clon_Distance.at(i).stddev(),Clon_Affinity.at(i).average(),Clon_Affinity.at(i).stddev(),Clon_Pon.at(i).average(),Clon_Pon.at(i).stddev(),Clon_Poff.at(i).average(),Clon_Poff.at(i).stddev(),Clone_m_theta.at(i).average(),Clone_m_theta.at(i).stddev(),Clon_Ag.at(i).average(),Clon_Ag.at(i).stddev(),Clon_Ag.at(i).sum,currentSim.Simulation_ID,currentSim.Concept_ID,currentSim.Sim_Seed);
            fclose(Clone_X);
        }
    
//clear stats
    Clone_State_counts.clear();
    Clone_State_counts.shrink_to_fit();
    clear_fileds();
}

void output::write_event(cell *Cellx, stringstream &sim_output) {
  sim_output << Cellx->event.str() << endl;
}

void output::write_event_2file(stringstream &sim_output) {
  FILE *event_data;
  string folder1 = output_path + "Event_data_seed_" + to_string(Output_ID)+".csv" ;
  event_data = fopen(folder1.c_str(), "a");
  fprintf(event_data, "%s", sim_output.str().c_str());
  fclose(event_data);
}
// for B cells
//cell states: 0-founder,
//1-unselected,
//2-contact_FDC,
//3-FDC_selected,
//4-contact_TC,
//5-TC_selected,
//6-recycled,
//7-apoptosis,
//8-TC_free,
//9-TC_connected,
//10-Plasma_Out,
//11-Plasma_in_GC,
//12-cell_state_counter
void output::close_event(B_cell *Cellx, stringstream &sim_output, double time,int SIM_ID, int concept_id,int Sim_seed) {
    Cellx->event << Cellx->cell_state << "," << Cellx->myBCR.distance<<","<<Cellx->MyAffinity<<","<< Cellx->myBCR.delta_Affinity<<","<<Cellx->myBCR.P_on<<","<<Cellx->myBCR.delta_Pon<<","<<Cellx->myBCR.P_off<<","<<Cellx->myBCR.delta_Poff<<","<<Cellx->myBCR.m_theta<<","<<Cellx->myBCR.delta_mtheta<<","<< Cellx->retained_Ag << "," << Cellx->total_number_of_divisions << "," << Cellx->myBCR.nMutFromGermline.beneficial << "," <<Cellx->myBCR.nMutFromGermline.harmful << "," << Cellx->nFDCcontacts << ","<< Cellx->fdc_binding_nums<<","<<Cellx->fdc_binding_nums_no_ag<<","<< Cellx->fdc_interaction_time_history << ","<< Cellx->Tc_interaction_history.first << "," << Cellx->Tc_interaction_history.second << ","<< Cellx->Selected_by_FDC << "," << Cellx->Selected_by_TC << "," << time<<","<<SIM_ID<<","<<concept_id<<","<<Sim_seed;
}

void output::Plasma_output(double endtime, simulation &currentSim,
                           parameters &p) {
  FILE *Plasma_cells_data;
  string folder1 = output_path + "Event_data_seed_" + to_string(Output_ID)+".csv";
  Plasma_cells_data = fopen(folder1.c_str(), "a");
  // Plasma data
  for (int j = 0; j < currentSim.ListP_cell.size(); j++) {
    Plasma_cell *Plasma = currentSim.ListP_cell.at(j);
      fprintf(Plasma_cells_data, "%d,%f,%d,%d,%d,%d,%f,%f,%f,%f,%f,%f,%f,%f,%f,%f,%d,%d,%d,%d,%d,%d,%f,%f,%f,%d,%d,%f,%d,%d,%d\n",
              Plasma->ID, Plasma->birth_time, Plasma->MID,Plasma->m_clonal_id,Plasma->clonal_id,Plasma->cell_state,Plasma->myBCR.distance,Plasma->MyAffinity,Plasma->myBCR.delta_Affinity,Plasma->myBCR.P_on,Plasma->myBCR.delta_Pon,Plasma->myBCR.P_off,Plasma->myBCR.delta_Poff,Plasma->myBCR.m_theta,Plasma->myBCR.delta_mtheta,Plasma->retained_Ag, Plasma->total_number_of_divisions, Plasma->myBCR.nMutFromGermline.beneficial, Plasma->myBCR.nMutFromGermline.harmful,Plasma->nFDCcontacts, Plasma->fdc_binding_nums,Plasma->fdc_binding_nums_no_ag,Plasma->fdc_interaction_time_history, Plasma->Tc_interaction_history.first, Plasma->Tc_interaction_history.second, Plasma->Selected_by_FDC,Plasma->Selected_by_TC,endtime,currentSim.Simulation_ID,currentSim.Concept_ID,currentSim.Sim_Seed);
  }
  fclose(Plasma_cells_data); 
}

