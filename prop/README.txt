 
Files:
  - prop.c, rand.c xsum.c, propjob.c, prop.config Makefile
  
Configurations:
  1. External file configuration: prop.config
    - Change parameters values in the prop.config
    - Save and it is ready to portray change in program. No need to recompile program
    
  2. Internal configuration: prop.c 
    - Parameters:
      - LAST: number of time steps from end of simulation to be considered as summary for values
      - TestPlot: 0 or 1; If 0, does not produce test group 'TestG' behavior plot
      - TestG: Index of group in a system to be checked for its individual group behavior
      
      - UsVsNature_Them: 1 or 2; 1 means Us vs. Nature games and 2 means Us vs. Them games
      - PIB_MFA: 1 or 2; If 1, program use update strategy of predicting individual behavior for leaders for option 3 update strategy by leaders; 
			 If 2, program use update strategy of mean field approximation for leaders for option 3 update strategy by leaders
      - ALLDATAFILE : If 1, generate all data files for individual runs and summary too; if 0, generate data files for only individual runs and averaged run but not for summary
      - AVERAGE_GRAPH_ONLY: 0 or 1; If 0, produces all individual run graphs also; if 1, produces only average graph for all runs
      - GRAPHS: not used
      - INIT_COM_EFFORT: inital commoners' effort (0 or 1)
      - INIT_LEAD_PUN_EFFORT: initial leaders' punishment effort between 0 and 1
      - INIT_LEAD_NORM_EFFORT: initial leaders' norm promoting effort
      - UPDATE_COM: 0 or 1; if 0, does not update commoners' efforts, if 1, updates commoners' effort
      - UPDATE_LEAD_PUN_EFFORT: 0 or 1; if 0, does not update leaders' punishment efforts, if 1, updates leaders' punishment effort
      - UPDATE_LEAD_NORM_EFFORT: 0 or 1; if 0, does not update leaders' norm promoting efforts, if 1, updates leaders' norm promoting effort
      - CLUSTER: 0 or 1; if 0, program assumes it is running in pc mode and displays graphs in gnuplot; if 1, program assumes it is running in cluster and does not produces graph
      - SKIP: Number of time steps to skip before taking data points for statistics
      
### PC simulation (running in personal machine):
  Usage:
    - make (compiles code)
    - ./prop prop.config (executes program)
	
  Output:
    Dynamics graphs: 
    - Outputs 4 graphs: efforts (x, y, z), payoff (pi_0, pi_1), utility function value (uc, ul) and group propduction P versus time steps (with SKIP as multiplier) in x-axis
    - x is efforts by commoners
    - y is punishment effort by leader
    - z is norm promoting effort by leader
    - color code: red - commoner related traits, green and blue: leader related traits 
	
### Cluster simulation (running multiple jobs in cluster)  
  Usage:
    - Edit parameter list in propjob.c
    - Edit parameters in prop.c
    - make (compiles code)
    - ./propjob 	
	
  Graph plotting:
    - change parameters in xsum.c
    - Vsim parameter: options - (1,2,3); indicates which strategy update option program will plot graphs for
    - make
    - ./xsum
  
  Graphs plotted by xsum.c:
    - Graphs plot are summarized as traits value versus benefit b are plotted 
    - In major y-axis, graphs for X, y, z and P are plotted against benefit b. 
      X: group effort by commoners
      Y: punishment effort by leader
      Z: norm promoting effort by leader
      P: group production
      b: benefit per individual
    - Graphs for different values of 'a' where Theta = a/(n+a)  
    - Graphs plotted for each y axis parameter are in same y-scale for comparison
  
  
Note: 
In this version, leader assumes commoners update with probability 'Upc' where Upc is sum of probability of updating commoners from configuration file.
To make leader assume all commoners update, uncomment old implementation (just one line statement) in method F() enclosed in block PIB_MFA = 1 
and comment out new implementation block. Comments in F() method should help identify new and old implementation blocks.

      
      
      

    