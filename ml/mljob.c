// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

/* sets of parameters to be simulated on */
int N[]        = {4, 8, 16};

double K[]     = {1, 1.5, 2};
double B2[]    = {1, 2, 4 };
double Theta[] = {0, 0.1, 0.2, 0.3};
double Eta[]   = {0, 0.1, 0.2, 0.3};
/*
double B0[]    = {4, 8};
double B1[]    = {1.0};
*/
double C0[]    = {1};
double C1[]    = {1};
double C2[]    = {1};
double C3[]    = {1};
double C4[]    = {1};
double V1[]    = {0.1};
double V2[]    = {0.2};
double V3[]    = {0.};
double x0[]    = {1};    // x0
double yy0[]   = {1};          // y0
double z0[]    = {1};          // z0
double X0[]    = {1};
double e[]     = {2, 4, 6};
double E[]     = {1};

/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 10;  

int      G         = 50;   
int      H         = 50;
int      T         = 1000;   

double   Sigma     = 0.1; 
double   Sigma_t   = 0.1; 


/* end of other basic parameters */


void createConfig(int i, int n, double k, double B2, double theta, double eta, double c0, double c1, double c2, double c3, double c4, double v1, double v2, double v3, double x, double ee, double EE, double x_0, double y_0, double z_0);
void createScript(int i, int n, double k, double B2, double theta, double eta, double c0, double c1, double c2, double c3, double c4, double v1, double v2, double v3, double x, double ee, double EE, double x_0, double y_0, double z_0);

int main()
{
  int n, k, b2, th, et, c0, c1, c2, c3, c4, v1, v2, v3, sc, x, ee, EE, x_0, y_0, z_0;
  int ns, ks, b2s, ths, ets, c0s, c1s, c2s, c3s, c4s, v1s, v2s, v3s, xs, es, Es, x0s, y0s, z0s;
  char scall[50];
  ns  = sizeof(N)/sizeof(int);
  ks  = sizeof(K)/sizeof(double);
  b2s  = sizeof(B2)/sizeof(double);
  ths  = sizeof(Theta)/sizeof(double);
  ets  = sizeof(Eta)/sizeof(double);
  c0s  = sizeof(C0)/sizeof(double);
  c1s  = sizeof(C1)/sizeof(double);
  c2s  = sizeof(C2)/sizeof(double);
  c3s  = sizeof(C3)/sizeof(double);
  c4s  = sizeof(C4)/sizeof(double);
  v1s  = sizeof(V1)/sizeof(double);
  v2s  = sizeof(V2)/sizeof(double);
  v3s  = sizeof(V3)/sizeof(double);   
  xs  = sizeof(X0)/sizeof(double);
  es  = sizeof(e)/sizeof(double);
  Es  = sizeof(E)/sizeof(double);
  x0s = sizeof(x0)/sizeof(double);
  y0s = sizeof(yy0)/sizeof(double);
  z0s = sizeof(z0)/sizeof(double);
  
  int i = 1;  
  for(n = 0; n < ns; n++)
    for(x = 0; x < xs; x++)
      for(k = 0; k < ks; k++)
	for(b2 = 0; b2 < b2s; b2++)
	  for(th = 0; th < ths; th++)
	    for(et = 0; et < ets; et++)
	      for(c0 = 0; c0 < c0s; c0++)
		for(c1 = 0; c1 < c1s; c1++)
		  for(c2 = 0; c2 < c2s; c2++)
		    for(c3 = 0; c3 < c3s; c3++)
		      for(c4 = 0; c4 < c4s; c4++)
			for(v1 = 0; v1 < v1s; v1++)
			  for(v2 = 0; v2 < v2s; v2++)
			    for(v3 = 0; v3 < v3s; v3++)
			      for(ee = 0; ee < es; ee++)
				for(EE = 0; EE < Es; EE++)
				  for(x_0 = 0; x_0 < x0s; x_0++)
				    for(y_0 = 0; y_0 < y0s; y_0++)
				      for(z_0 = 0; z_0 < z0s; z_0++, i++){				      
					createConfig(i, N[n], K[k], B2[b2], Theta[th], Eta[et], C0[c0], C1[c1], C2[c2], C3[c3], C4[c4], V1[v1], V2[v2], V3[v3], X0[x], e[ee], E[EE], x0[x_0], yy0[y_0], z0[z_0]);
					createScript(i, N[n], K[k], B2[b2], Theta[th], Eta[et], C0[c0], C1[c1], C2[c2], C3[c3], C4[c4], V1[v1], V2[v2], V3[v3], X0[x], e[ee], E[EE], x0[x_0], yy0[y_0], z0[z_0]);   
					sprintf(scall, "qsub mlscript.%d.sh", i);
					sc = system(scall);
					if(sc) printf("Error submitting jobs!!\n");
				      }
  printf("\n%d\n", i);
  return 0;
}

void createConfig(int i, int n, double k, double b2, double theta, double eta, double c0, double c1, double c2, double c3, double c4, double v1, double v2, double v3, double x, double ee, double EE, double x_0, double y_0, double z_0 )
{
  FILE *fp;
  char cfile[20];
  sprintf(cfile, "ml%d.config", i);
  if(!(fp = fopen(cfile, "w"))){
    printf("Error!! %s couldn't be created.\n", cfile);
    return;
  } 
  fprintf(fp, "unsigned Seed      = %u; \n", Seed);
  fprintf(fp, "int      Runs      = %d; \n", Runs);  
  fprintf(fp, "int      T         = %d; \n", T);
  fprintf(fp, "\n");
  fprintf(fp, "int      N         = %d; \n", n);
  fprintf(fp, "int      G         = %d; \n", G);
  fprintf(fp, "int      H         = %d; \n", H);
  fprintf(fp, "\n");
  fprintf(fp, "double   K         = %lf; \n", k);
  fprintf(fp, "double   B2        = %lf; \n", b2);
  fprintf(fp, "double   Theta     = %lf; \n", theta);
  fprintf(fp, "double   Eta       = %lf; \n", eta);  
  fprintf(fp, "\n");  
  fprintf(fp, "double   C0        = %lf; \n", c0);
  fprintf(fp, "double   C1        = %lf; \n", c1);
  fprintf(fp, "double   C2        = %lf; \n", c2);  
  fprintf(fp, "double   C3        = %lf; \n", c3);
  fprintf(fp, "double   C4        = %lf; \n", c4);
  fprintf(fp, "\n");  
  fprintf(fp, "double   V1        = %lf; \n", v1);
  fprintf(fp, "double   V2        = %lf; \n", v2);
  fprintf(fp, "double   V3        = %lf; \n", v3);  
  fprintf(fp, "\n");  
  fprintf(fp, "double   x0        = %lf; \n", x_0);
  fprintf(fp, "double   y0        = %lf; \n", y_0);
  fprintf(fp, "double   z0        = %lf; \n", z_0);
  fprintf(fp, "double   X0        = %lf; \n", x);
  fprintf(fp, "\n");
  fprintf(fp, "double   e        = %lf; \n", ee);
  fprintf(fp, "double   E        = %lf; \n", EE);  
  fprintf(fp, "\n");  
  fprintf(fp, "double   Sigma     = %lf; \n", Sigma);
  fprintf(fp, "double   Sigma_t   = %lf; \n", Sigma_t);
  fclose(fp);
}

void createScript(int i, int n, double k, double b2, double theta, double eta, double c0, double c1, double c2, double c3, double c4, double v1, double v2, double v3, double x, double ee, double EE, double x_0, double y_0, double z_0)
{
  char sfile[20];
  //char path[256];
  FILE *fp;
  sprintf(sfile, "mlscript.%d.sh", i);
  fp = fopen(sfile, "w");
  fprintf(fp, "#PBS -l nodes=1\n");
  fprintf(fp, "#PBS -N ml.%d.%.2lf%.2lf.%.2lf.%.2lf.%.2lf%.2lf.%.2lf.%.2lf.%.2lf.%.2lf.%.2lf.%.2lf.%.2lf.%.2lf.%.2lf.%.2lf.%.2lf.%.2lf\n", n, k, b2, theta, eta, c0, c1, c2, c3, c4, x, ee, EE, v1, v2, v3, x_0, y_0, z_0);
  fprintf(fp, "#PBS -l walltime=24:00:00\n" );
  fprintf(fp, "cd $PBS_O_WORKDIR\n");
  //getcwd(path, 256);
  fprintf(fp, "./ml ml%d.config\n", i);
  fclose(fp);  
}


/** Usage:
    compile : gcc -o mljob mljob.c
    run     : ./mljob
**/

