// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

/* sets of parameters to be simulated on */
int n[]        = {1};
int K[]     = {4};

double b[]     = {1};
double delta[] = {0.05, 0.2};
double k[]     = {0.125, 0.5};
double Theta[] = {1, 2, 4};

double cx[]    = {1};
double cy[]    = {0.125};
double cz[]    = {1};

double Vc1    = 0.01;
double Vc2    = 0.24;
double Vc3    = 0.; 
double Vl1    = 0.01;
double Vl2    = 0.24;
double Vl3    = 0.; 

double X0[]    = {4, 8, 16};    // x0
double s0[]   = {0.0,  0.5};          // s0
double s1[]     = {0.0,  0.5};
double e[]    = {1};
double z_0[]  = {1};


/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 20;  

int      G         = 500;   
int      H         = 1;
int      T         = 1500;   
double   Lambda    = 4;
double   m         = 0.5;

double   Sigma     = 0.1; 
double   Sigma_t   = 0.1; 


/* end of other basic parameters */

void createConfig(int j, int _n, double _b, int _K, double _delta, double _k, double _Theta, double _cx, double _cy, double _cz, double _X0, double _s0, double _s1, double _e, double _z_0)
{
  FILE *fp;
  char cfile[200];
  sprintf(cfile, "prop%d.config", j);
  if(!(fp = fopen(cfile, "w"))){
    printf("Error!! %s couldn't be created.\n", cfile);
    return;
  } 
 fprintf(fp, "unsigned Seed      = %u;\n", Seed);
 fprintf(fp, "int      Runs      = %d;\n", Runs);
 fprintf(fp, "int      T         = %d;\n", T);
  
 fprintf(fp, "int      n         = %d;\n", _n);
 fprintf(fp, "int      G         = %d;\n", G); 
 fprintf(fp, "int      K         = %d;\n", _K);
  
 fprintf(fp, "double   b         = %lf;\n", _b);  
 fprintf(fp, "double   cx        = %lf;\n", _cx);
 fprintf(fp, "double   cy        = %lf;\n", _cy);
 fprintf(fp, "double   cz        = %lf;\n", _cz);     
  
 fprintf(fp, "double   k         = %lf;\n", _k);
 fprintf(fp, "double   delta     = %lf;\n", _delta);  
 fprintf(fp, "double   Theta     = %lf;\n", _Theta);
 
 fprintf(fp, "double   Vc1       = %lf;\n", Vc1);
 fprintf(fp, "double   Vc2       = %lf;\n", Vc2);
 fprintf(fp, "double   Vc3       = %lf;\n", Vc3); 
 fprintf(fp, "double   Vl1       = %lf;\n", Vl1);
 fprintf(fp, "double   Vl2       = %lf;\n", Vl2);
 fprintf(fp, "double   Vl3       = %lf;\n", Vl3); 
 fprintf(fp, "double   m         = %lf;\n", m);
  
 fprintf(fp, "double   X0        = %lf;\n", _X0);  
 fprintf(fp, "double   z0        = %lf;\n", _z_0);
 fprintf(fp, "double   e         = %lf;\n", _e);
  
 fprintf(fp, "double   s0        = %lf;\n", _s0);  
 fprintf(fp, "double   s1        = %lf;\n", _s1); 
 fprintf(fp, "double   Lambda    = %lf;\n", Lambda);    
  
 fprintf(fp, "double   Sigma     = %lf;\n", Sigma);

  fclose(fp);
}

void createScript(int j, int _n, double _b, int _K, double _delta, double _k, double _Theta, double _cx, double _cy, double _cz, double _X0, double _s0, double _s1, double _e, double _z_0)
{
  char sfile[200];
  //char path[256];
  FILE *fp;
  sprintf(sfile, "propscript.%d.sh", j);
  if(!(fp = fopen(sfile, "w"))){
    printf("Error!! %s couldn't be created.\n", sfile);
    return;
  }

  fprintf(fp, "#PBS -l nodes=1\n");
  fprintf(fp, "#PBS -N p.%d.%.2lf.%.2lf.%.2lf.%.2lf%.2lf.%.2lf.%.2lf.%.2lf%.2lf\n", _K, _k, _b, _Theta, _cx, _cy, _cz, _X0, _s0, _s1);
  fprintf(fp, "#PBS -l walltime=24:00:00\n" );
  fprintf(fp, "cd $PBS_O_WORKDIR\n");
  //getcwd(path, 256);
  fprintf(fp, "./propsim prop%d.config\n", j);
  fclose(fp);  
}
int main()
{
  int i[14], is[14], j, sc;
  char scall[50];
  is[0]  = sizeof(n)/sizeof(int);
  is[1]  = sizeof(b)/sizeof(double);
  is[2]  = sizeof(K)/sizeof(int);
  is[3]  = sizeof(delta)/sizeof(double);
  is[4]  = sizeof(k)/sizeof(double);
  is[5]  = sizeof(Theta)/sizeof(double);
  is[6]  = sizeof(cx)/sizeof(double);
  is[7]  = sizeof(cy)/sizeof(double);
  is[8]  = sizeof(cz)/sizeof(double);
  is[9]  = sizeof(X0)/sizeof(double);
  is[10]  = sizeof(s0)/sizeof(double);
  is[11]  = sizeof(s1)/sizeof(double);
  is[12]  = sizeof(e)/sizeof(double);
  is[13]  = sizeof(z_0)/sizeof(double);
    
  for(j = 1, i[0] = 0; i[0] < is[0]; i[0]++)
    for(i[1] = 0; i[1] < is[1]; i[1]++)
      for(i[2] = 0; i[2] < is[2]; i[2]++)
	for(i[3] = 0; i[3] < is[3]; i[3]++)
	  for(i[4] = 0; i[4] < is[4]; i[4]++)
	    for(i[5] = 0; i[5] < is[5]; i[5]++)
	      for(i[6] = 0; i[6] < is[6]; i[6]++)
		for(i[7] = 0; i[7] < is[7]; i[7]++)
		  for(i[8] = 0; i[8] < is[8]; i[8]++)
		    for(i[9] = 0; i[9] < is[9]; i[9]++)
		      for(i[10] = 0; i[10] < is[10]; i[10]++)
			for(i[11] = 0; i[11] < is[11]; i[11]++)
			  for(i[12] = 0; i[12] < is[12]; i[12]++)
			    for(i[13] = 0; i[13] < is[13]; i[13]++, j++){
			      if(j>0){			     
				createConfig(j, n[i[0]], b[i[1]], K[i[2]], delta[i[3]], k[i[4]], Theta[i[5]], cx[i[6]], cy[i[7]], cz[i[8]], X0[i[9]], s0[i[10]], s1[i[11]], e[i[12]], z_0[i[13]]);
				createScript(j, n[i[0]], b[i[1]], K[i[2]], delta[i[3]], k[i[4]], Theta[i[5]], cx[i[6]], cy[i[7]], cz[i[8]], X0[i[9]], s0[i[10]], s1[i[11]], e[i[12]], z_0[i[13]]);   
				sprintf(scall, "qsub propscript.%d.sh", j);
				sc = system(scall); 
				//printf("a scl \n");
				if(sc) printf("Error submitting jobs!!\n");
    // 			    if(j%150 == 0) {
    // 			      printf("submitted # %d\n", j); 
    // 			      sleep(00);
    // 			    }
			      }

			    }
  printf("\n%d\n", j);
  return 0;
}




/** Usage:
    compile : gcc -o propjob propjob.c
    run     : ./propjob
**/

 
