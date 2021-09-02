#ifndef BCR_H
#define BCR_H
#include <vector>
#include "parameters.h"
#include <string>
void initialize_Seeds (parameters &p);

struct mutation_array
{
    int silent =0;
    int beneficial=0;
    int harmful=0;
    int key_mutations_beneficial=0;
    int lethal=0;

};
class BCR
{
    //This class includes the implementation of B-cell receptor
public:
    //Shape-space method
    BCR(parameters &p);
    vector<int> BCReceptor;
    double pMut; // Probability of mutation
    mutation_array nMutFromGermline; // Number of mutations from the beginning founder cell
    void mutateBCR(parameters &p);
    double getMyAffinity4Ag(parameters& p);
    double distance;
    string print_BCR();
    bool checked_initial_dissociation;
    void calculate_koff(parameters &p);
    void calculate_on_off_P(parameters &p);
    double P_on;
    double P_off;
    double m_theta; // This angle  is used to specify the contribution of Pon/Poff according to each mutational distance.
    double delta_Pon;
    double delta_Poff;
    double delta_mtheta;
    double delta_Affinity;
};


#endif // BCR_H
