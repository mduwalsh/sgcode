#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

/* sets of parameters to be simulated on */
int G[]        = {10, 20, 40};
int N[]        = {4, 8, 16};
double K[]     = {1, 1.5, 2};
double B2[]    = {1, 2, 4};
double Theta[] = {0., 0.1, 0.2, 0.3};
double Eta[]   = {0, 0.1, 0.2, 0.3};
double C0[]    = {1};
double C1[]    = {1};
double C2[]    = {1};
double C3[]    = {1};
double C4[]    = {1};
double V1[]    = {0.};
double V2[]    = {0.};
double V3[]    = {0.3};

double x_0[]   = {1};
double y_0[]   = {0.25};
double z_0[]   = {0.1};
double X0[]    = {1};
double e[]     = {1, 2, 4};
double E[]     = {1};

/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 4;  

int      H         = 100;
int      T         = 1000;   

double   Sigma     = 0.1; 
double   Sigma_t   = 0.1; 
/* end of other basic parameters */

int main()
{
  int g, n, k, b2, th, et, c0, c1, c2, c3, c4, v1, v2, v3, sc, x, ee, EE, x00, y00, z00;
  int gs, ns, ks, b2s, ths, ets, c0s, c1s, c2s, c3s, c4s, v1s, v2s, v3s, xs, es, Es, x00s, y00s, z00s;
  char scall[200];
  gs  = sizeof(G)/sizeof(int);
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
  x00s = sizeof(x_0)/sizeof(double);
  y00s = sizeof(y_0)/sizeof(double);
  z00s = sizeof(z_0)/sizeof(double);
  
  int i = 1;  
  for(g = 0; g < gs; g++)
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
				    for(z00 = 0; z00 < z00s; z00++)
				      for(x00 = 0; x00 < x00s; x00++)
					for(y00 = 0; y00 < y00s; y00++)
					  for(k = 0; k < ks; k++, i++){
					    //y_0[y00] = x_0[x00];
					    //z_0[y00] = x_0[x00];
					    X0[y00] = N[n]*x_0[x00];
					    E[EE] = e[ee];
					    sprintf(scall, "./fplot %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", G[g], N[n], K[k], B2[b2], Theta[th], Eta[et], C0[c0], C1[c1], C2[c2], C3[c3], C4[c4], V1[v1], V2[v2], V3[v3], X0[x], e[ee], E[EE], x_0[x00], y_0[y00], z_0[z00], Runs);
					    sc = system(scall);
					    if(sc) printf("Error submitting jobs!!\n");
					  }      
  printf("\n%d\n", i);
  return 0;
}


/** Usage:
    compile : gcc -o allplot allplot.c
    run     : ./allplot
**/



