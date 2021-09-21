# Publication
This software accompanies the publication

**From affinity selection to kinetic selection in germinal center modeling**

Danial Lashgari1, Elena Merino Tjiero1, Michael Meyer-Hermann2,3, Marit J. van Gils5, Huub Hoefsloot5,#, Antoine H.C. van Kampen1,5,#,*

1 Bioinformatics Laboratory, Epidemiology and Data Science, Amsterdam Public Health research institute, Amsterdam Institute for Infection and Immunity, Amsterdam, the Netherlands.
2 Department for Systems Immunology and Braunschweig Integrated Centre of Systems Biology, Helmholtz Centre for Infection Research, Braunschweig, Germany.
3 Institute for Biochemistry, Biotechnology and Bioinformatics, Technische Universität Braunschweig, Braunschweig, Germany.
4 Laboratory of Experimental Virology, Department of Medical Microbiology, Center for Infection and Immunity Amsterdam, Academic Medical Center, University of Amsterdam, the Netherlands
5 Biosystems Data Analysis, Swammerdam Institute for Life Sciences, University of Amsterdam, Amsterdam, the Netherlands.

# Contributed equally 
*Corresponding author: 
Antoine H.C. van Kampen
Bioinformatics Laboratory
Epidemiology & Data Science
Amsterdam University Medical Centers
Meibergdreef 9, 1105 AZ Amsterdam, the Netherlands
a.h.vankampen@amsterdamumc.nl (AvK)
tel. +31-20-5667096
![image](https://user-images.githubusercontent.com/68376494/134162160-1b5b8a45-d2a3-4d36-9a5e-ec6d1f4caaea.png)


# Project
## GC_from_Affinity_to_Kinetics
This repository includes code ONLY. The Agent-based model is based on [Hyphasma](https://www.helmholtz-hzi.de/en/research/research-topics/immune-response/systems-immunology/our-research/) (e.g., Michael-Meyer Hermann, 2012) and ´.  
The different scenarios can be modelled by codes provided in different folders:
   SC0: The Reference Scenario in the manuscript
   SC1: Scenario-1 in the manuscript
   SC2: SCenario-2 in the manuscript

## Software
All software is written in C++. All the parameters used are borrowed from the original model [Meyer-Hermann et al., 2012] except parameter Theta that is explained in the manuscript. To change the parameters, use parameters.cpp

### Run simulations 
#Make project on Unix-based systems using KONKOFF.pro and qmake in following steps
# qmake KONKOFF.pro
# make
# run ./KONKOFF 
#### Options 
  # -s "random seed" 
  # -o "Path of output folder"

## References
* Meyer-Hermann, M., et al. (2012). "A theory of germinal center B cell selection, division, and exit." Cell Rep 2(1): 162-174.
* Robert, P. A., et al. (2017). "How to Simulate a Germinal Center." Methods Mol Biol 1623: 303-334.
	
