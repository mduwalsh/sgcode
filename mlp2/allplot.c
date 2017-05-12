#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>

/* sets of parameters to be simulated on */
int G[]        = {10, 20, 40};
int n[]        = {4, 8, 16};
double b[]     = {1, 1.5, 2};
double B[]    = {1, 2, 4};
double Theta_u[] = {0., 0.1, 0.2, 0.3};
double Theta_d[] = {0., 0.1, 0.2, 0.3};
double Eta_u[]   = {0, 0.1, 0.2, 0.3};
double Eta_d[]   = {0, 0.1, 0.2, 0.3};
double cx[]    = {1};
double cy[]    = {1};
double cz[]    = {1};

double V1[]    = {0.};
double V2[]    = {0.};
double V3[]    = {0.3};

double x_0[]   = {1};
double y_0[]   = {0.25};
double z_0[]   = {0.1};
double Y0[]    = {1};
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
  int g, nn, bb, bB, tu, td, eu, ed, ccx, ccy, ccz, v1, v2, v3, sc, yY0, ee, EE, x00, y00, z00;
  int gs, nns, bbs, bBs, tus, tds, eus, eds, cxs, cys, czs, v1s, v2s, v3s, yY0s, es, Es, x00s, y00s, z00s;
  char scall[200];
  gs  = sizeof(G)/sizeof(int);
  nns  = sizeof(n)/sizeof(int);
  bbs  = sizeof(b)/sizeof(double);
  bBs  = sizeof(B)/sizeof(double);
  tus  = sizeof(Theta_u)/sizeof(double);
  tds  = sizeof(Theta_d)/sizeof(double);
  eus  = sizeof(Eta_u)/sizeof(double);
  eds  = sizeof(Eta_d)/sizeof(double);
  cxs  = sizeof(cx)/sizeof(double);
  cys  = sizeof(cy)/sizeof(double);
  czs  = sizeof(cz)/sizeof(double);
  v1s  = sizeof(V1)/sizeof(double);
  v2s  = sizeof(V2)/sizeof(double);
  v3s  = sizeof(V3)/sizeof(double);   
  yY0s  = sizeof(Y0)/sizeof(double);
  es  = sizeof(e)/sizeof(double);
  Es  = sizeof(E)/sizeof(double);
  x00s = sizeof(x_0)/sizeof(double);
  y00s = sizeof(y_0)/sizeof(double);
  z00s = sizeof(z_0)/sizeof(double);
  
  int i = 1;  
  for(g = 0; g < gs; g++)
    for(nn = 0; nn < nns; n++)
      for(yY0 = 0; yY0 < yY0s; yY0++)
	for(bb = 0; bb < bbs; bb++)
	  for(bB = 0; bB < bBs; bB++)
	    for(tu = 0; tu < tus; tu++)
	      for(td = 0; td < tds; td++)
		for(eu = 0; eu < eus; eu++)
		  for(ed = 0; ed < eds; ed++)
		    for(ccx = 0; ccx < cxs; ccx++)
		      for(ccy = 0; ccy < cys; ccy++)
			for(ccz = 0; ccz < czs; ccz++)
			  for(v1 = 0; v1 < v1s; v1++)
			    for(v2 = 0; v2 < v2s; v2++)
			      for(v3 = 0; v3 < v3s; v3++)
				for(ee = 0; ee < es; ee++)
				  for(EE = 0; EE < Es; EE++)
				    for(z00 = 0; z00 < z00s; z00++)
				      for(x00 = 0; x00 < x00s; x00++)
					for(y00 = 0; y00 < y00s; y00++, i++){    
					    sprintf(scall, "./fplot %d %d %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %lf %d", G[g], n[nn], b[bb], B[bB], Theta_u[tu], Theta_d[td], Eta_u[eu], Eta_d[td], cx[ccx], cy[ccy], cz[ccz], V1[v1], V2[v2], V3[v3], e[ee], E[EE], Y0[yY0], x_0[x00], y_0[y00], z_0[z00], Runs);
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



