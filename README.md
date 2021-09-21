

# Publication

This software accompanies the publication
 **From affinity selection to kinetic selection in germinal center modeling**

<sup>
Danial Lashgari1, Elena Merino Tjiero<sup>1</sup>, Michael Meyer-Hermann<sup>2,3</sup>, Marit J. van Gils <sup>4</sup>, Huub Hoefsloot<sup>5#</sup>, Antoine H.C. van Kampen<sup>1,5,#*</sup></sup>

<sup> 1 Bioinformatics Laboratory, Epidemiology and Data Science, Amsterdam Public Health research institute, Amsterdam Institute for Infection and Immunity, Amsterdam, the Netherlands.
2 Department for Systems Immunology and Braunschweig Integrated Centre of Systems Biology, Helmholtz Centre for Infection Research, Braunschweig, Germany.
3 Institute for Biochemistry, Biotechnology and Bioinformatics, Technische Universität Braunschweig, Braunschweig, Germany.
4 Laboratory of Experimental Virology, Department of Medical Microbiology, Center for Infection and Immunity Amsterdam, Academic Medical Center, University of Amsterdam, the Netherlands
5 Biosystems Data Analysis, Swammerdam Institute for Life Sciences, University of Amsterdam, Amsterdam, the Netherlands.
</sup>
<sup>
	
#Contributed equally 
*Corresponding author: 
Antoine H.C. van Kampen
Bioinformatics Laboratory
Epidemiology & Data Science
Amsterdam University Medical Centers
Meibergdreef 9, 1105 AZ Amsterdam, the Netherlands
a.h.vankampen@amsterdamumc.nl (AvK)
tel. +31-20-5667096 
</sup>


# Software
This repository includes codes and nesseccary files for runing simulatios. The Agent-based model is based on [Meyer-Hermann et al., 2012] and [Robert et al., 2017].

Different scenarios can be modelled by codes provided in different folders:
   SC0: The Reference Scenario in the manuscript
   SC1: Scenario-1 in the manuscript
   SC2: SCenario-2 in the manuscript

All software is written in C++. All the parameters used are borrowed from the original model [Meyer-Hermann et al., 2012] except parameter $\theta$ that is explained in the manuscript. To change the parameters, use parameters.cpp

### Run simulations   
Run simulations Unix-based systems using following steps
qmake `KONKOFF.pro` 
make
run ./KONKOFF </br>
 
#### Options 
 -s "random seed" 
 -o "Path of output folder"

## References
* Meyer-Hermann, M., et al. (2012). "A theory of germinal center B cell selection, division, and exit." Cell Rep 2(1): 162-174.
* Robert, P. A., et al. (2017). "How to Simulate a Germinal Center." Methods Mol Biol 1623: 303-334.
	
