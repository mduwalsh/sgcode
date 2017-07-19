// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<time.h>


/* sets of parameters to be simulated on */
int Vop[]      = {1, 2};
int n[]        = {4, 8, 16};
int Lambda[]   = {4, 8, 16};
double b[]     = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
double B[]     = {1};
double delta[] = {0.05, 0.1, 0.2};
double DELTA[] = {0};
double k[]     = {0.125, 0.25, 0.5};
double K[]     = {0};

double Theta_ua[] = {0, 1, 2};
double Eta_ua[]   = {0};
double Theta_da[] = {0};
double Eta_da[]   = {0};

double cx[]    = {1};
double cy[]    = {1};
double cz[]    = {1};
double cp[]    = {0.125};
double cq[]    = {1};

double x_0[]   = {1};    // x0
double y_0[]   = {1};          // y0
double z_0[]   = {1};          // z0
double X0[]    = {2, 4, 8};
double Y0[]    = {1};
double e[]     = {1};
double E[]     = {1};
double Rho[]   = {0.0};

double Dummy[] = {1.0};

/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 20;  

int      G         = 8;   
int      H         = 100;
int      T         = 1500;   
int      L         = 2;
double   m         = 0.1;

double   Sigma     = 0.1; 
double   Sigma_t   = 0.1; 


int NoOfJobsInBatch = 50;             // numbers of jobs to run before submitting other number of jobs after some sleep time
int SleepTime       = 30;             // sleep time in seconds; time before next batch of jobs scheduled

/* end of other basic parameters */

int _Vop, _n, _Lambda;
double _b, _B, _delta, _DELTA, _K, _k, _Theta_ua, _Theta_da, _Eta_ua, _Eta_da, _cx, _cy, _cz, _cp, _cq, _x_0, _y_0, _z_0, _X0, _Y0, _e, _E, _Rho;

void createConfig(int j)
{
  FILE *fp;
  char cfile[20];
  sprintf(cfile, "mlp%d.config", j);
  if(!(fp = fopen(cfile, "w"))){
    printf("Error!! %s couldn't be created.\n", cfile);
    return;
  } 
  fprintf(fp, "unsigned Seed      = %u; \n", Seed);
  fprintf(fp, "int      Runs      = %d; \n", Runs);  
  fprintf(fp, "int      T         = %d; \n", T);
  fprintf(fp, "int      Lambda    = %d; \n", _Lambda);
  fprintf(fp, "\n");
  fprintf(fp, "int      n         = %d; \n", _n);
  fprintf(fp, "int      G         = %d; \n", G);
  fprintf(fp, "int      H         = %d; \n", H);
  fprintf(fp, "\n");
  fprintf(fp, "double   b         = %lf; \n", _b);
  fprintf(fp, "double   B         = %lf; \n", _b);
  fprintf(fp, "double   cx        = %lf; \n", _cx);
  fprintf(fp, "double   cy        = %lf; \n", _cy);
  fprintf(fp, "double   cz        = %lf; \n", _cz);  
  fprintf(fp, "double   cp        = %lf; \n", _cp);
  fprintf(fp, "double   cq        = %lf; \n", _cq);  
  fprintf(fp, "\n");  
  fprintf(fp, "unsigned L         = %d; \n", L);
  fprintf(fp, "double   k         = %lf; \n", _k);
  fprintf(fp, "double   K         = %lf; \n", _n*_k);  
  fprintf(fp, "double   delta     = %lf; \n", _delta);
  fprintf(fp, "double   DELTA     = %lf; \n", _n*_delta);
  fprintf(fp, "\n");  
  fprintf(fp, "double   Theta_ua  = %lf; \n", _Theta_ua);
  fprintf(fp, "double   Theta_da  = %lf; \n", _Theta_da);
  fprintf(fp, "double   Eta_ua    = %lf; \n", _Eta_ua);
  fprintf(fp, "double   Eta_da    = %lf; \n", _Eta_da);
  fprintf(fp, "\n");
  fprintf(fp, "unsigned Vop       = %d; \n", _Vop);
  fprintf(fp, "double   m         = %lf; \n", m); 
  fprintf(fp, "double   Rho       = %lf; \n", _Rho);
  fprintf(fp, "\n");  
  fprintf(fp, "double   x0        = %lf; \n", _x_0);
  fprintf(fp, "double   y0        = %lf; \n", _y_0);
  fprintf(fp, "double   z0        = %lf; \n", _z_0);
  fprintf(fp, "double   X0        = %lf; \n", _X0);
  fprintf(fp, "double   Y0        = %lf; \n", _Y0);
  fprintf(fp, "\n");
  fprintf(fp, "double   e        = %lf; \n", _e);
  fprintf(fp, "double   E        = %lf; \n", _E);  
  fprintf(fp, "\n");  
  fprintf(fp, "double   Sigma     = %lf; \n", Sigma);
  fprintf(fp, "double   Sigma_t   = %lf; \n", Sigma_t);
  fclose(fp);
}

void createScript(int j, int _n)
{
  char sfile[20];
  //char path[256];
  FILE *fp;
  sprintf(sfile, "mlpscript.%d.sh", j);
  fp = fopen(sfile, "w");
  fprintf(fp, "#PBS -l nodes=1\n");
  fprintf(fp, "#PBS -N ml.%d.%d\n", _n, j);
  fprintf(fp, "#PBS -l walltime=24:00:00\n" );
  fprintf(fp, "cd $PBS_O_WORKDIR\n");
  //getcwd(path, 256);
  fprintf(fp, "./mlp mlp%d.config\n", j);
  fclose(fp);  
}
int main()
{
  int i[26], is[26], j, sc;
  char scall[50];
  char *logfile = "log_job_count.txt";
  FILE *fp = fopen(logfile, "a");  
  time_t ltime;
  ltime = time(NULL);
  struct tm *tm;
  tm = localtime(&ltime);  
  fprintf(fp, "\n%s\n", asctime(tm));        // time stamp to log file  
  fclose(fp);
  
  is[0]  = sizeof(n)/sizeof(int);
  is[1]  = sizeof(b)/sizeof(double);
  is[2]  = sizeof(B)/sizeof(double);
  is[3]  = sizeof(delta)/sizeof(double);
  is[4]  = sizeof(DELTA)/sizeof(double);
  is[5]  = sizeof(K)/sizeof(double);
  is[6]  = sizeof(k)/sizeof(double);
  is[7]  = sizeof(Theta_ua)/sizeof(double);
  is[8]  = sizeof(Theta_da)/sizeof(double);
  is[9]  = sizeof(Eta_ua)/sizeof(double);
  is[10]  = sizeof(Eta_da)/sizeof(double);
  is[11]  = sizeof(cx)/sizeof(double);
  is[12]  = sizeof(cy)/sizeof(double);
  is[13]  = sizeof(cz)/sizeof(double);
  is[14]  = sizeof(Vop)/sizeof(int);
  is[15]  = sizeof(Rho)/sizeof(double);
  is[16]  = sizeof(Lambda)/sizeof(int);
  is[17]  = sizeof(x_0)/sizeof(double);
  is[18]  = sizeof(y_0)/sizeof(double);
  is[19]  = sizeof(z_0)/sizeof(double);
  is[20]  = sizeof(X0)/sizeof(double);
  is[21]  = sizeof(Y0)/sizeof(double);
  is[22]  = sizeof(e)/sizeof(double);
  is[23]  = sizeof(E)/sizeof(double);
  is[24]  = sizeof(cp)/sizeof(double);
  is[25]  = sizeof(cq)/sizeof(double);
  
  for(j = 1, i[0] = 0; i[0] < is[0]; i[0]++){  // n
    _n = n[i[0]];
    for(i[1] = 0; i[1] < is[1]; i[1]++){  //b
      _b = b[i[1]];
      for(i[2] = 0; i[2] < is[2]; i[2]++){  // B
	_B = B[i[2]];
	for(i[3] = 0; i[3] < is[3]; i[3]++){  // delta
	  _delta = delta[i[3]];
	  for(i[4] = 0; i[4] < is[4]; i[4]++){  // DELTA
	    _DELTA = DELTA[i[4]];
	    for(i[5] = 0; i[5] < is[5]; i[5]++){  // K
	      _K = K[i[5]];
	      for(i[6] = 0; i[6] < is[6]; i[6]++){  // k
		_k = k[i[6]];
		for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_ua
		  _Theta_ua = Theta_ua[i[7]];
		  for(i[8] = 0; i[8] < is[8]; i[8]++){  // Theta_da
		    _Theta_da = Theta_da[i[8]];
		    for(i[9] = 0; i[9] < is[9]; i[9]++){  // Eta_ua
		      _Eta_ua = Eta_ua[i[9]];
		      for(i[10] = 0; i[10] < is[10]; i[10]++){  // Eta_da
			_Eta_da = Eta_da[i[10]];
			for(i[11] = 0; i[11] < is[11]; i[11]++){  // cx
			  _cx = cx[i[11]];
			  for(i[12] = 0; i[12] < is[12]; i[12]++){  // cy
			    _cy = cy[i[12]];
			    for(i[13] = 0; i[13] < is[13]; i[13]++){  // cz
			      _cz = cz[i[13]];
			      for(i[14] = 0; i[14] < is[14]; i[14]++){  // Vop
				_Vop = Vop[i[14]];
				for(i[15] = 0; i[15] < is[15]; i[15]++){  // Rho
				  _Rho = Rho[i[15]];
				  for(i[16] = 0; i[16] < is[16]; i[16]++){  // Lambda
				    _Lambda = Lambda[i[16]];
				    for(i[17] = 0; i[17] < is[17]; i[17]++){  // x_0
				      _x_0 = x_0[i[17]];
				      for(i[18] = 0; i[18] < is[18]; i[18]++){  // y_0
					_y_0 = y_0[i[18]];
					for(i[19] = 0; i[19] < is[19]; i[19]++){  // z_0
					  _z_0 = z_0[i[19]];
					  for(i[20] = 0; i[20] < is[20]; i[20]++){  // X0
					    _X0 = X0[i[20]];
					    for(i[21] = 0; i[21] < is[21]; i[21]++){  // Y0
					      _Y0 = Y0[i[21]];
					      for(i[22] = 0; i[22] < is[22]; i[22]++){  // e
						_e = e[i[22]];
						for(i[23] = 0; i[23] < is[23]; i[23]++){  // E
						  _E = E[i[23]];
						  for(i[24] = 0; i[24] < is[24]; i[24]++){  // cp
						    _cp = cp[i[24]];
						    for(i[25] = 0; i[25] < is[25]; i[25]++, j++){ // cq 
						      _cq = cq[i[25]];
						      
						      createConfig(j);
						      createScript(j, _n);   
						      sprintf(scall, "qsub mlpscript.%d.sh", j);
						      sc = system(scall);
						      if(sc) printf("Error submitting jobs!!\n");
						      
						      if(j % NoOfJobsInBatch  == 0){         // jobs in batch check		
							printf("scheduled jobs count: %d ", j);
							sleep(SleepTime);                    // sleep before submitting next batch of jobs
							// log into file every batch completed
							fp = fopen(logfile, "a");
							fprintf(fp, "%d \t", j);
							fclose(fp);
							sleep(1);
						      }						      
						    }
						  }
						}
					      }
					    }
					  }
					}
				      }
				    }
				  }
				}
			      }
			    }
			  }
			}
		      }
		    }
		  }
		}
	      }
	    }
	  }
	}
      }
    }
  }
  printf("\n%d\n", j);
  return 0;
}




/** Usage:
    compile : gcc -o mljob mljob.c
    run     : ./mljob
**/

