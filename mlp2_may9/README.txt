## Read Me file for mlp.c
 
 Compiling and Running:
  - make
  - ./mlp mlp.config
 
 Files required:
  - mlp.c
  - mlp.config
  - rand.c
  - Makefile

 *********** Configuration file ***************
 
 Parameters on which simulation are run can be configured in configuration file mlp.config

 Seed:   seed to initialize psuedorandom number generator
 Runs:   number of runs for each simulation
 T:      time period of simulation for each run
 n:      number of commoners in a group
 G:      number of groups in a polity
 H:      number of polities in a system;
 b:      benefit per commoner
 B:      benefit per chief
 cx:     cost coefficient for commoners effort
 cy:     cost coeffficient for leaders effort
 cz:     cost coefficient for chiefs effort 
 L:      number of leaders with least efforts to be punished by chief in a group
 k:      punishment to commoner by leader
 K:      punishment to leader by chief
 delta:  punishment cost to leader for punishing commoner
 DELTA:  punishment cost to chief for punishing chief
 Theta_u:tax factor imposed by leader on commoner
 Theta_d:reward factor given by leader to commoner
 Eta_u:  tax factor imposed by chief on leader
 Eta_d:  reward factor given by chief to group
 Vop:    options - {1, 2, 3}; indicates set of probabilities of method to update commoner, leader and chief's strategies; for example, 1 indicates first set of options {...} in array Vc[][3], Vl[][3] and Vcf[][3] in mlp.c 
 Rho:    probability that polity plays to us vs them game
 m:      probability of migration (probability of copying strategy from other group or polity)
 x0:     half effort parameter of commoners
 y0:     half effort parameter of leaders
 z0:     half effort parameter of chiefs
 Y0:     half strength equivalent of polity
 e:      efficiency of leader
 E:      efficiency of chief
 Sigma:  standard deviation of effort distribution (normal)
 Sigma_t: standard deviation of tax distribution (normal)

 
 *********** program settings defined in mlp.c ***************
 
  CLUSTER:     			simulation on cluster or not (1 or 0)
  SKIP:        			time interval between calculation with snapshot of states
  STU:         			time period from last which is considered as summary of simulation
  ALLDATAFILE: 			if 1, generate all data files for individual runs and summary too 
  AVERAGE_GRAPH_ONLY : 		if 1, it suppresses graphs of individual runs while running in pc (CLUSTER = 0) mode and produces only averaged graphs only
  GRAPHS:      			if 1, saves graphs as png files 
  UPDATE_COM:  			1 or 0; if 1, update strategy of commoner; if 0, does not update strategy of commoner
  UPDATE_LEAD_EFFORT: 		1 or 0; if 1, update effort made by leader; if 0, does not update effort 
  UPDATE_CHIEF_EFFORT: 		1 or 0; if 1, update strategy of chief; if 0, does not update effort
  UPDATE_LEAD_PUN_EFFORT: 	1 or 0; if 1, update punishment effort 'p' of leader; if 0, do not update
  UPDATE_CHIEF_PUN_EFFORT: 	1 or 0; if 1, update punishment effort 'q' of chief; if 0, do not update
  INIT_COM_EFFORT: 		initial effort of commoners
  INIT_LEAD_EFFORT: 		initial effort of leaders
  INIT_CHIEF_EFFORT: 		initial effort of chiefs 
  INIT_LEAD_PUN_EFFORT: 	initial punishment effort of leaders
  INIT_CHIEF_PUN_EFFORT: 	initial punishment effort of chiefs    
  Vc[][3]: 			set of options (indicated by Vop in mlp.config file) to be used by commoner during strategy update; each option set denotes probablities of using random mutation, selective copy and myopic optimization method in format {random_mutation, selective_copy, optimization}
  Vl[][3]: 			set of options (indicated by Vop in mlp.config file) to be used by leader during strategy update; each option set denotes probablities of using random mutation, selective copy and myopic optimization method in format {random_mutation, selective_copy, optimization}
  Vcf[][3]: 			set of options (indicated by Vop in mlp.config file) to be used by chief during strategy update; each option set denotes probablities of using random mutation, selective copy and myopic optimization method in format {random_mutation, selective_copy, optimization}
  
  
 ************************************** Output ******************************************
  Dynamics graphs: 
  - Outputs 4 graphs: efforts (x, y, z), payoff (pi_0, pi_1, pi_2), punishment efforts (p, q), and group production and polity strength (P & Q) versus time steps (with SKIP as multiplier) in x-axis
    Effort graph:
      - red line (x) is average effort per commoner 
      - green line (y) is effort by leader
      - blue line (z) is effort by chief
    Payoff graph:
      - red line (pi_0) is payoff of commoner 
      - green line (pi_1) is payoff of leader
      - blue line (pi_2) is payoff of chief
    Punishment graph:
      - green line (p) is punishment effort by leader
      - blue line (q) is punishment effort by chief
    Production and polity strength graph:
      - red line (P) is group production 
      - green line (Q) is polity strength
  
  
  
  
  
  
 
  
