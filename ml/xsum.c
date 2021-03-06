#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h> 
#include<stdbool.h>

#define DEBUG 1
#if DEBUG
#define __USE_GNU
#include <fenv.h>       /* enable floating exceptions */
unsigned long Where;    // debugging counter
#endif

/* Parameters for plotting */
double G[]     = {10};
double N[]     = {4, 8, 16};
double K[]     = { 1, 1.5, 2};
double B2[]    = { 1, 2, 4};
double Theta[] = {0.0, 0.1, 0.2, 0.3};
double Eta[]   = {0.0, 0.1, 0.2, 0.3};
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
double e[]     = {1, 2, 4};



// Parameters to specify type of plots
#define WINPOPUP 0           // 0 or 1; if 0, png file is created; if 1, window pop up is generated
#define STACKED 0            //  1 or 0; if 1 VBAR should be defined

#define TRAIT "effort"       // options are: effort, payoff, tax

#define yAXIS "N"            // variation in plots in inner rows in a graph
#define YAXIS "Theta"            // variation in plots in outer rows in a graph
#define xAXIS1 "e"           // horizontal 1st inner x axis
#define xAXIS2 "K"            // horizontal 1st inner x axis in a plot in a graph
#define XAXIS "Eta"            // variation in plots in columns in a graph
#define VBAR  "K"        // variation in vertical stacked bars in a cluster of a histogram

#define NP 18

double Pn[NP], *P[NP];        // Pn stores value of parameter currently under multiple loops of parameters for each indexed parameter; P array points to reference of parameter used in plotting
char Pa[NP][10];               // stores alphbet strings array used to denote specific parameter in creating data file 
int Ps[NP];                   // Ps stores no. of values under one parameter indexed by its index.

// methods declaration
void set_param();
void prep_infile(char *fname, char *appnd);  
void prep_outfile(char *file, char *appnd);
void plot_graphs();
void prepare_data_files();


/* Parameters for plotting */
#define EPS 0                     // 1: eps image; 0: png image

#define BOXWIDTH 0.75


#define PREC "7.4" /* printing precision */
 
// write data in d array to file fp in n columns
/*
 * fp : file pointer
 * d  : pointer to data to be written
 * ba : alphabet representation of parameter for which data is to be written
 * bi : value of parameter for wihic data is to be written
 * n  : no. of columns of data to be written
 * pad: 1 or 0; 1 if no. of columms of data for each row is different. it pads columns with zero to max possible columns; in this case, data is collected for different individuals in a group. group size may vary
 */
void write_data(FILE **fp, double *d, char *ba, double bi, int n, int pad)
/*
 * d: data array 
 */
{
  int i, j, k, nx;
  /*for(j=3; j--;){
    fprintf(fp[j], "%s=%.2f ", ba, bi);
  }*/
  for(i = 0; i < 3; i++){
    fprintf(fp[i], "%"PREC"lf \t", d[i]);
  }
  /*
  if(pad){
    for(nx = 0, i = 0; i < sizeof(N)/sizeof(double); i++){
      if(nx < N[i]) nx = N[i];
    }
    if(n < nx){
      for(k = n; k < nx; k++)
	fprintf(fp, "%.4f ", 0.0);
    }
  }
  fprintf(fp, "\n");
  */
}

// calc averages of data in columns 
/*
 * fp: file pointer
 * d : pointer array in which data to be stored
 * n : no. of columns of data to be read/written
 * toNormal : normalize data or not; 1 or 0
 */
int calc_data(FILE *fp, double *d, int n, bool toNormal)
{
  int pos, nc = 0; 
  double val, sum = 0.0;
  char line[1024], *str;
  while( fgets(line, 1024, fp) !=NULL ) {
      // Just search for the latest line, do nothing in the loop
  } 
  // just read first line and do nothing
  str = line;
  while(sscanf(str, "%lf%n", &val, &pos)==1){
    d[nc] = val;
    if(toNormal){
      sum += d[nc];
    }
    str += pos;
    nc++;   
  }
  if(nc == n){
    if(toNormal){          // normalization of values
      while(nc--){
	d[nc] /= sum;
      }
    }
    return 1;
  }
  else 
    return 0;
}


int main(int argc, char **argv)
{      
  if(STACKED){
    if(strcmp(VBAR, "") == 0 || !VBAR){
      printf("Please define VBAR parameter! \n");
      exit(1);
    }
  }  
  set_param();                                            // sets P array with pointers to array of parameters like B, N, Beta,...
  prepare_data_files();                                   // readies data files for plotting graphs  
  plot_graphs();  
  return 0;
}

void set_param(){
  int i, h = -1, p = -1, z = -1, y = -1, x=-1, v = -1;
  double *Param[NP];
  char ParamStr[NP][10];
  char ParamA[NP][10];
  int ParamS[NP];
  int gs, ns, ks, b2s, ths, ets, c0s, c1s, c2s, c3s, c4s, v1s, v2s, v3s, xs, es, Es, x0s, y0s, z0s;
  
  // calculating size of arrays
  gs  = sizeof(G)/sizeof(double);
  ns  = sizeof(N)/sizeof(double);
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
  //xs  = sizeof(X0)/sizeof(double);
  es  = sizeof(e)/sizeof(double);
  x0s  = sizeof(x_0)/sizeof(double);
  y0s  = sizeof(y_0)/sizeof(double);
  z0s  = sizeof(z_0)/sizeof(double);
  
  
  // copying array address
  Param[0] = G;
  Param[1] = N; 
  Param[2] = K; 
  Param[3] = B2;
  Param[4] = Theta;
  Param[5] = Eta;
  Param[6] = C0;
  Param[7] = C1;
  Param[8] = C2;
  Param[9] = C3;
  Param[10] = C4;
  Param[11] = V1;
  Param[12] = V2;
  Param[13] = V3;
  Param[14] = e;
  Param[15] = x_0;
  Param[16] = y_0;
  Param[17] = z_0;
  
  // array size
  ParamS[0] = gs;
  ParamS[1] = ns; 
  ParamS[2] = ks; 
  ParamS[3] = b2s;
  ParamS[4] = ths;
  ParamS[5] = ets;
  ParamS[6] = c0s;
  ParamS[7] = c1s;
  ParamS[8] = c2s;
  ParamS[9] = c3s;
  ParamS[10] = c4s;
  ParamS[11] = v1s;
  ParamS[12] = v2s;
  ParamS[13] = v3s;
  ParamS[14] = es;
  ParamS[15] = x0s;
  ParamS[16] = y0s;
  ParamS[17] = z0s;
  // setting array name as strings
  sprintf(ParamStr[0], "G");
  sprintf(ParamStr[1], "N");
  sprintf(ParamStr[2], "K");
  sprintf(ParamStr[3], "B2");
  sprintf(ParamStr[4], "Theta");
  sprintf(ParamStr[5], "Eta");
  sprintf(ParamStr[6], "C0");
  sprintf(ParamStr[7], "C1");
  sprintf(ParamStr[8], "C2");
  sprintf(ParamStr[9], "C3");
  sprintf(ParamStr[10], "C4");
  sprintf(ParamStr[11], "V1");
  sprintf(ParamStr[12], "V2");
  sprintf(ParamStr[13], "V3");
  sprintf(ParamStr[14], "e");
  sprintf(ParamStr[15], "x_0");
  sprintf(ParamStr[16], "y_0");
  sprintf(ParamStr[17], "z_0");
  // setting alphabets representing parameters
  sprintf(ParamA[0], "g");
  sprintf(ParamA[1], "n");
  sprintf(ParamA[2], "k");
  sprintf(ParamA[3], "b2");
  sprintf(ParamA[4], "theta");
  sprintf(ParamA[5], "eta");
  sprintf(ParamA[6], "c0");
  sprintf(ParamA[7], "c1");
  sprintf(ParamA[8], "c2");
  sprintf(ParamA[9], "c3");
  sprintf(ParamA[10], "c4");
  sprintf(ParamA[11], "v1");
  sprintf(ParamA[12], "v2");
  sprintf(ParamA[13], "v3");
  sprintf(ParamA[14], "e");
  sprintf(ParamA[15], "x_0");
  sprintf(ParamA[16], "y_0");
  sprintf(ParamA[17], "z_0");
  
  // assigning array reference according to alphabets for easy location of main axes's parameter
  for(i = 0; i < NP; i++){
    if( strcmp(ParamStr[i], xAXIS1) == 0){               // minor x-axis 1
      h = i;
      P[NP-1] = Param[i];
      sprintf(Pa[NP-1],"%s", ParamA[i]);
      Ps[NP-1] = ParamS[i];
    }
    if( strcmp(ParamStr[i], xAXIS2) == 0){               // minor x-axis 2
      v = i;
      P[NP-2] = Param[i];
      sprintf(Pa[NP-2],"%s", ParamA[i]);
      Ps[NP-2] = ParamS[i];
    }
    else if( strcmp(ParamStr[i], XAXIS) == 0 ){         // major x-axis
      p = i;
      P[NP-3] = Param[i];
      sprintf(Pa[NP-3], "%s", ParamA[i]);
      Ps[NP-3] = ParamS[i];         
    }
    else if( strcmp(ParamStr[i], yAXIS) == 0 ){         // minor y-axis
      x = i;
      P[NP-4] = Param[i];
      sprintf(Pa[NP-4], "%s", ParamA[i]);
      Ps[NP-4] = ParamS[i];
    }
    else if( strcmp(ParamStr[i], YAXIS) == 0 ){         // major y-axis
      z = i;
      P[NP-5] = Param[i];
      sprintf(Pa[NP-5], "%s", ParamA[i]);
      Ps[NP-5] = ParamS[i];
    }
/*#if STACKED
    else if( strcmp(ParamStr[i], VBAR) == 0){
      y = i;
      P[NP-4] = Param[i];
      sprintf(Pa[NP-4], "%s", ParamA[i]);
      Ps[NP-4] = ParamS[i];
    }    
#endif*/
  }
  if( strcmp(XAXIS, "") == 0 ){
    p = -1;
  }
  int k;
  for(k = 0, i = 0; i < NP; i++){        
    if( i != h && i != p && i != z && i != y && i != x && i != v){            
      P[k] =  Param[i];
      sprintf(Pa[k],"%s", ParamA[i]);
      Ps[k] = ParamS[i];
      k++;      
    }
  }  
} 

// prepares name of file from which data to be read
void prep_infile(char *fname, char *appnd)
{
  char str[20];
  strcpy(str, "");
  strcpy(fname, "");
  int i, n = 0;
  
  // for g
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "g") == 0){
      sprintf(str,"g%02.0f", Pn[i]);
      break;
    }
  }  
  strcat(fname, str);
  
  // for n
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "n") == 0){
      sprintf(str,"n%02.0f", Pn[i]);
      n = (int)Pn[i];
      break;
    }
  }  
  strcat(fname, str);
  
  // for k
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "k") == 0){
      sprintf(str,"k%0.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);  
  
  // for b2
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "b2") == 0){
      sprintf(str,"b2%0.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str); 
  
  // for theta
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "theta") == 0){
      sprintf(str,"theta%.2f", Pn[i]);
      break;
    }    
  }     
  strcat(fname, str);
  
  // for eta
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "eta") == 0){
      sprintf(str,"eta%.2f", Pn[i]);
      break;
    }    
  }   
  strcat(fname, str);
  
  // for c0
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "c0") == 0){
      sprintf(str,"c0%.2f", Pn[i]);
      break;
    }    
  }     
  strcat(fname, str);
  
  // for c1
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "c1") == 0){
      sprintf(str,"c1%.2f", Pn[i]);
      break;
    }    
  }   
  strcat(fname, str);
  
  // for c2
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "c2") == 0){
      sprintf(str,"c2%.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);
  
  // for c3
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "c3") == 0){
      sprintf(str,"c3%.2f", Pn[i]);
      break;
    }    
  }   
  strcat(fname, str);
  
  // for c4
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "c4") == 0){
      sprintf(str,"c4%.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);
  
  // for v1
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "v1") == 0){
      sprintf(str,"v1%.2f", Pn[i]);
      break;
    }    
  }    
  strcat(fname, str);
  
  // for v2
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "v2") == 0){
      sprintf(str,"v2%.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);
  
  strcpy(str, "");
  // for v3
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "v3") == 0){
      sprintf(str,"v3%.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);

  // for X0
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "x_0") == 0){
      sprintf(str,"X0%.2f", n*Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);  
  
  // for e
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "e") == 0){
      sprintf(str,"e%.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);  
  
  // for E
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "e") == 0){
      sprintf(str,"E%.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);    
  
  // for x0
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "x_0") == 0){
      sprintf(str,"x0%.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);    
  
  // for y0
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "y_0") == 0){
      sprintf(str,"y0%.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);    
  
  // for z0
  for(i = 0; i < NP; i++){
    if(strcmp(Pa[i], "z_0") == 0){
      sprintf(str,"z0%.2f", Pn[i]);
      break;
    }    
  }  
  strcat(fname, str);    
  
  strcat(fname, appnd);  
} 

// prepares file to write data for graphs
void prep_outfile(char *file, char *appnd)
{
#if !STACKED
  sprintf(file, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s", Pa[0], Pn[0], Pa[1], Pn[1], Pa[2], Pn[2], Pa[3], Pn[3], Pa[4], Pn[4], Pa[5], Pn[5], Pa[6], Pn[6], Pa[7], Pn[7], Pa[8], Pn[8], Pa[9], Pn[9], Pa[10], Pn[10], Pa[11], Pn[11], Pa[12], Pn[12], Pa[NP-5], Pn[NP-5], Pa[NP-4], Pn[NP-4], appnd);
/*#else
  sprintf(file, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s", Pa[0], Pn[0], Pa[1], Pn[1], Pa[2], Pn[2], Pa[3], Pn[3], Pa[4], Pn[4], Pa[5], Pn[5], Pa[6], Pn[6], Pa[7], Pn[7], Pa[8], Pn[8], Pa[9], Pn[9], Pa[10], Pn[10], Pa[11], Pn[11], Pa[12], Pn[12], Pa[13], Pn[13], Pa[14], Pn[14], Pa[NP-3], Pn[NP-3], Pa[NP-2], Pn[Np-2], appnd);*/
#endif
}

void prepare_data_files()
{
  char infile[300], edata[3][300], str[300];
  double *d = NULL, v = 0;
  char a[20];
  FILE *fp1, *fp2[3];
  int j, pi, qi, ai, bi, ci, di, ei, fi, gi, hi, ri, si, ti, ui, vi, wi, xi, yi, zi, ni = 0, pad = 0, ft = 0;   
  d = malloc(3*sizeof(double));
  ni = 3;                                                               // no. of data columns
  if(strcmp(TRAIT, "tax") == 0) ni = 2;                                 // only two columns of data for tax
  
  // looping through all parameter arrays
  for(ri = 0; ri < Ps[0]; ri++){    
    Pn[0] = P[0][ri];
      for(si = 0; si < Ps[1]; si++){
	Pn[1] = P[1][si];
	for(ti = 0; ti < Ps[2]; ti++){
	  Pn[2] = P[2][ti];
	  for(ui = 0; ui < Ps[3]; ui++){    	    
	    Pn[3] = P[3][ui];
	    for(vi = 0; vi < Ps[4]; vi++){	      
	      Pn[4] = P[4][vi];
	      for(wi = 0; wi < Ps[5]; wi++){
		Pn[5] = P[5][wi];		
		for(xi = 0; xi < Ps[6]; xi++){
		  Pn[6] = P[6][xi];		
		  for(ai = 0; ai < Ps[7]; ai++){
		    Pn[7] = P[7][ai];
		    for(bi = 0; bi < Ps[8]; bi++){
		      Pn[8] = P[8][bi];
		      for(ci = 0; ci < Ps[9]; ci++){
			Pn[9] = P[9][ci];
			for(di = 0; di < Ps[10]; di++){
			  Pn[10] = P[10][di];
			  for(ei = 0; ei < Ps[11]; ei++){
			    Pn[11] = P[11][ei];
			    for(fi = 0; fi < Ps[12]; fi++){
			      Pn[12] = P[12][fi];
			      //for(gi = 0; gi < Ps[13]; gi++){
				//Pn[13] = P[13][gi];
				for(hi = 0; hi < Ps[NP-5]; hi++){
				  Pn[NP-5] = P[NP-5][hi];
#if !STACKED
				  for(yi = 0; yi < Ps[NP-4]; yi++){
				    Pn[NP-4] = P[NP-4][yi];        
#endif
				    for(j = 0; j < 3; j++){
				      if(strcmp(TRAIT, "effort") == 0) sprintf(str, "ehist_%d.dat", j);
				      else if(strcmp(TRAIT, "tax") == 0) sprintf(str, "thist_%d.dat", j);
				      else if(strcmp(TRAIT, "payoff") == 0) sprintf(str, "phist_%d.dat", j);			 
				      
				      prep_outfile(edata[j], str);                                                               // data for histogram of property of interest
				      fp2[j] = fopen(edata[j], "w");
				      // write headers to histogram data file
				      fprintf(fp2[j], "%s's ", Pa[NP-3]);     
				      for(qi = 0; qi < Ps[NP-2]; qi++){
					for(pi = 0; pi < Ps[NP-1]; pi++){ 
					  fprintf(fp2[j], "\t%.2lf/%.2lf  ", P[NP-1][pi],P[NP-2][qi]);   
					}
				      }
				      fprintf(fp2[j], "\n");
				    }
				    for(zi = 0; zi < Ps[NP-3]; zi++){
				      Pn[NP-3] = P[NP-3][zi];
				      #if !STACKED  
					sprintf(a,"%s", Pa[NP-3]); v = P[NP-3][zi];					
				      #endif
				      if(strcmp(a, "theta") == 0) sprintf(a, "{/Symbol q}");
				      if(strcmp(a, "eta") == 0) sprintf(a, "{/Symbol h}");
				      fprintf(fp2[0], "\"%s=%.2f\" ", a, v);
				      fprintf(fp2[1], "\"%s=%.2f\" ", a, v);
				      fprintf(fp2[2], "\"%s=%.2f\" ", a, v);
				      for(qi = 0; qi < Ps[NP-2]; qi++){
					Pn[NP-2] = P[NP-2][qi];					
					
					ft = 0;
					for(pi = 0; pi < Ps[NP-1]; pi++){ 					  
					  Pn[NP-1] = P[NP-1][pi];					  
					  // files to be read for any property of interest for set of paramters
					  if(strcmp(TRAIT, "effort") == 0) prep_infile(infile, "xsum.dat");
					  else if(strcmp(TRAIT, "tax") == 0) prep_infile(infile, "tsum.dat");
					  else if(strcmp(TRAIT, "payoff") == 0) prep_infile(infile, "psum.dat");
					  
					  if(!(fp1 = fopen(infile, "r"))){
					    printf("Could not open file '%s'\n",infile);
					    ft = 1;
					    break; //exit(1);
					  }
					  if(calc_data(fp1, d, ni, 0)){   // calculate data and not normalize			    
					    fclose(fp1);
					    write_data(fp2, d, a, v, ni, pad);   // write data to histogram data file for efforts // j is passed as delta for linear spacing    
					  }
					  else{
					    fclose(fp1);
					    printf("error in data read '%s': column count does not match!!\n", infile);
					    exit(1);
					  }    
					}// pi loop ends					
				      }// qi loop ends 
				      for(j=3; j--;){
					fprintf(fp2[j],"\n");
				      }				      			      
				    } // zi loop ends
				    for(j=3; j--;){
				      fclose(fp2[j]);
				      if(ft) remove(edata[j]);     // if input file not read, destroy output file opened
				    }	
				  }  // yi loop ends
				}  // hi loop ends
			      //}   // gi loop ends
			    }   // fi loop ends
			  }  // ei loop ends
			}  // di loop ends
		      }  // ci loop ends
		    } // bi loop ends
		  }  // ai loop ends
	      }  // xi loop ends
	    }  // wi loop ends
	  }  // vi loop ends
	}  // ui loop ends
      }  // ti loop ends
    }  // si loop ends       
  }  // ri loop ends
   // reading data and prepare histogram data file for multiple values of b and delta but same other parameters completed
}

void hist_sum(char *edata, int n, char *title, char *outputfile, int wid)
{
  int i, j, k, c, l;  
  char albt[3][20];
  FILE *fp1;
  /*if(strcmp(TRAIT, "effort") == 0) j = 1;
  else if(strcmp(TRAIT, "payoff") == 0) j = 3;
  else if(strcmp(TRAIT, "tax") == 0) j = 2;  
  */
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode  
#if STACKED 
  int pi;
#endif
  char data[300]; char str[20];
  //int bs  = sizeof(B)/sizeof(double); // no. of Bs in array corresponds to no. of clusters
#if WINPOPUP 
  fprintf(gp, "set term x11 %d \n", wid);
#else
  if(EPS){
    fprintf(gp, "set term postscript eps enhanced color solid font \"Helvetica,25\" \n");
    fprintf(gp, "set output '%s.eps' \n", outputfile);
  }
  else{
    fprintf(gp, "set term pngcairo size 1600,900 enhanced color solid font \"Helvetica,8\" \n");
    fprintf(gp, "set output '%s.png' \n", outputfile);  
  }
#endif
  
  fprintf(gp, "set key autotitle columnheader\n");
  fprintf(gp, "set key outside \n");
  fprintf(gp, "set style data histogram \n");  
#if STACKED
  fprintf(gp, "set style histogram rowstacked gap 1 title offset 0.0, -1 \n");
#else
  fprintf(gp, "set style histogram cluster gap 1 title offset 0.0, 0\n");
  fprintf(gp, "set offset -0.5, -0.5, 0, 0\n");
#endif
  fprintf(gp, "set style fill solid border -1\n");
  fprintf(gp, "set boxwidth %lf relative\n", BOXWIDTH);  
  
  fprintf(gp, "set lmargin 10 \n");  
  
  //fprintf(gp, "set multiplot layout %d, %d title '%s' offset 0.02, 0 \n", Ps[NP-4]*Ps[NP-3], Ps[NP-2], "");                        // set subplots layout    
  fprintf(gp, "set multiplot layout %d, %d title '%s' offset 0.02, 0 \n", Ps[NP-4]*Ps[NP-3], 3, "");                        // set subplots layout   
  
#if STACKED
    fprintf(gp, "unset xtics \n");
#else
    fprintf(gp, "set auto x \n");    
#endif  
  
  fprintf(gp, "set style arrow 1 nohead lw 2 \n");
  fprintf(gp, "set ylabel ''\n");
  fprintf(gp, "set label 10 '%s' at screen %.4f,%.4f center font \",12\"\n", TRAIT, 0.02, 1.0-0.008);  
      
  if(strcmp(XAXIS, "") != 0){                                          // when xaxis is specified
    
    fprintf(gp, "ymax = 0.0\n");
    
    for(c = 0, j = 0; j < Ps[NP-5]; j++){                                        // major y axis
      
      for(k = 0; k < Ps[NP-4]; k++){                                      // minor y axis      
	
	for(i = 0; i < 3; i++){                                    // major xaxis	
	  if(strcmp(TRAIT, "effort") == 0) sprintf(str, "ehist_%d.dat", i);
	  else if(strcmp(TRAIT, "payoff") == 0) sprintf(str, "phist_%d.dat", i);    
	  else if(strcmp(TRAIT, "tax") == 0) sprintf(str, "thist_%d.dat", i);    
	  sprintf(data, "%s%s%.2f%s%.2f%s", edata, Pa[NP-5], P[NP-5][j], Pa[NP-4], P[NP-4][k], str);
	  if(!(fp1 = fopen(data, "r"))){
	    printf("Could not open file '%s'\n",data);
	    break; //exit(1);
	  }
	  fclose(fp1);
	  /*fprintf(gp, "stats '%s' using 2 prefix 'A%d' nooutput\n", data, c);
	  fprintf(gp, "stats '%s' using 3 prefix 'B%d' nooutput\n", data, c);	  
	  fprintf(gp, "if(A%d_max > ymax) {ymax = A%d_max}\n", c, c);
	  fprintf(gp, "if(B%d_max > ymax) {ymax = B%d_max}\n", c, c);
	  if(strcmp(TRAIT, "tax") !=0){
	    fprintf(gp, "stats '%s' using 4 prefix 'C%d' nooutput\n", data, c);
	    fprintf(gp, "if(C%d_max > ymax) {ymax = C%d_max}\n", c, c);
	  }	  
	  */
	}
      }
    }
    //fprintf(gp, "ymax = ymax+0.01\n");
    //fprintf(gp, "set format y \"%%.2f\"\n");
    //fprintf(gp, "set ytics ymax/3\n");
    //fprintf(gp, "set yrange [0:ymax]\n");
    //fprintf(gp, "set autoscale y \n");
      
    for(c = 0, j = 0; j < Ps[NP-5]; j++){                  // major y axis
      if(strcmp(Pa[NP-5], "xx") == 0)
	sprintf(albt[2], "X0");
      else 
	sprintf(albt[2], "%s", Pa[NP-5]);
      if(strcmp(albt[2], "theta") == 0) sprintf(albt[2], "{/Symbol q}");
      if(strcmp(albt[2], "eta") == 0) sprintf(albt[2], "{/Symbol h}");
      fprintf(gp, "set label 1 '%s = %g' at screen %.4f,%.4f rotate by 90 center font \",12\"\n", albt[2], P[NP-5][j], 0.01, 1.0 + 0.5/(double)(Ps[NP-5])  - ((1.0-0.00)/(double)(Ps[NP-5]))*(j+1));
      
      for(k = 0; k < Ps[NP-4]; k++){                                      // minor y axis     
	if(strcmp(Pa[NP-4], "xx") == 0)
	  sprintf(albt[1], "X0");
	else 
	  sprintf(albt[1], "%s", Pa[NP-4]);
	if(strcmp(albt[1], "theta") == 0) sprintf(albt[1], "{/Symbol q}");
	if(strcmp(albt[1], "eta") == 0) sprintf(albt[1], "{/Symbol h}");
	fprintf(gp, "set label 2 '%s = %g' at screen %.4f,%.4f rotate by 90 center font \",10\"\n", albt[1], P[NP-4][k], 0.02, 1.0 -0.005 + 0.5/(double)(Ps[NP-4]*Ps[NP-5]) - ( 1.0/((Ps[NP-4]*Ps[NP-5])) )*((j*Ps[NP-4])+k+1));
	fprintf(gp, "set lmargin 10.0 \n");
	
	for(i = 0; i < 3; i++, c++){                                    // major xaxis	
	  if(strcmp(TRAIT, "effort") == 0) sprintf(str, "ehist_%d.dat", i);
	  else if(strcmp(TRAIT, "payoff") == 0) sprintf(str, "phist_%d.dat", i);    
	  else if(strcmp(TRAIT, "tax") == 0) sprintf(str, "thist_%d.dat", i);    
	  sprintf(data, "%s%s%.2f%s%.2f%s", edata, Pa[NP-5], P[NP-5][j], Pa[NP-4], P[NP-4][k], str);
	  if(!(fp1 = fopen(data, "r"))){
	    printf("Could not open file '%s'\n",data);
	    break; //exit(1);
	  }
	  fclose(fp1);
	  if(i == 0){ sprintf(albt[0], "Commoner");}
	  else if(i == 1){ sprintf(albt[0], "Leader");}
	  else if(i == 2){ sprintf(albt[0], "Chief");}
	  //fprintf(gp, "set label 3 '%s' at screen %.4f,%.4f center font \",10\"\n", albt[0], (1.0/(Ps[NP-2]+1))*(i+1), 1.0-0.01);	  
	  fprintf(gp, "set label 3 '%s' at screen %.4f,%.4f center font \",10\"\n", albt[0], (1.0/(3+1))*(i+1), 1.0-0.01);	  
  #if !STACKED
	  sprintf(data, "%s%s%.2f%s%.2f%s", edata, Pa[NP-5], P[NP-5][j], Pa[NP-4], P[NP-4][k], str);  
	  if(!(fp1 = fopen(data, "r"))){
	    printf("Could not open file '%s'\n",data);
	    break; //exit(1);
	  }
	  fclose(fp1);
	  // set yrange
	  fprintf(gp, "ymax = 0.0\n");
	  fprintf(gp, "set autoscale y\n");
	  for(l = 0; l < (Ps[NP-1]*Ps[NP-2]); l++){
	    fprintf(gp, "stats '%s' using %d prefix 'A%d%d' nooutput\n", data, l+2, l, c);
	    fprintf(gp, "if(A%d%d_max > ymax) {ymax = A%d%d_max}\n", l,c, l, c);
	  }
	  //fprintf(gp, "stats '%s' using 2 prefix 'A%d' nooutput\n", data, c);
	  //fprintf(gp, "stats '%s' using 3 prefix 'B%d' nooutput\n", data, c);	  
	  //fprintf(gp, "if(A%d_max > ymax) {ymax = A%d_max}\n", c, c);
	  //fprintf(gp, "if(B%d_max > ymax) {ymax = B%d_max}\n", c, c);
	  /*if(strcmp(TRAIT, "tax") !=0){
	    fprintf(gp, "stats '%s' using 4 prefix 'C%d' nooutput\n", data, c);
	    fprintf(gp, "if(C%d_max > ymax) {ymax = C%d_max}\n", c, c);
	  }*/
	  fprintf(gp, "unset autoscale y\n");
	  fprintf(gp, "ymax = ymax+0.01\n");
	  fprintf(gp, "set format y \"%%.2f\"\n");
	  fprintf(gp, "set ytics ymax/3\n");
	  fprintf(gp, "set yrange [0:ymax]\n"); 
    
	  fprintf(gp, "plot \\\n");  
	  if(strcmp(TRAIT, "tax") == 0){
	    fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) lc col notitle \\\n", Ps[NP-1]*Ps[NP-2]+1, data);
	  }
	  else{
	    fprintf(gp, "for[col=2:%d] '%s' u col:xtic(1) lc (col-2)/%d+1 notitle \\\n", Ps[NP-1]*Ps[NP-2]+1, data, Ps[NP-2]);  
	  }
	  fprintf(gp, "\n"); 
	  fprintf(gp, "set lmargin 0.0 \n");
  #endif  
	}
	fprintf(gp, "unset label 3\n");
      }
      fprintf(gp, "unset label 2\n");
    }
    fprintf(gp, "unset label 1\n");
  }  
  fprintf(gp, "unset label 10\n");
  fprintf(gp, "unset multiplot \n");
  
  fflush(gp); 
  pclose(gp);
}

void plot_graphs()
{
  char str[300], title[300], edata[300];   
  int ai, bi, ci, di, ei, fi, gi, hi, ri, si, ti, ui, vi, wi, xi, ni = 0;    
  
  // for simulation of ml.c only 
  int cnt=0;
  if(strcmp(TRAIT, "effort") == 0) cnt = 1;
  else if(strcmp(TRAIT, "payoff") == 0) cnt = 3;
  else if(strcmp(TRAIT, "tax") == 0) cnt = 2; 
  
  int wid = 0;
  
  // plotting summary histograms all at a time // the one with cluster base is removed from for loop (like for B in original case)
  
  ni = 3;
  if(strcmp(TRAIT, "tax")==0) ni = 2;
  for(ri = 0; ri < Ps[0]; ri++){          
    for(si = 0; si < Ps[1]; si++){   
     
      for(ti = 0; ti < Ps[2]; ti++){   
	for(ui = 0; ui < Ps[3]; ui++){    
	  for(vi = 0; vi < Ps[4]; vi++){  
	    for(wi = 0; wi < Ps[5]; wi++){      
	      for(xi = 0; xi < Ps[6]; xi++){    
		for(ai = 0; ai < Ps[7]; ai++){  
		  for(bi = 0; bi < Ps[8]; bi++){  
		    for(ci = 0; ci < Ps[9]; ci++){          
		      for(di = 0; di < Ps[10]; di++){  
			for(ei = 0; ei < Ps[11]; ei++){  
			  for(fi = 0; fi < Ps[12]; fi++){
			    //for(gi = 0; gi < Ps[13]; gi++){
			      //for(hi = 0; hi < Ps[14]; hi++){
#if !STACKED		
				if(strcmp(XAXIS, "") != 0){                      // when xaxis is specified
				  sprintf(edata, "%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f%s%.2f", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][ai], Pa[8], P[8][bi], Pa[9], P[9][ci], Pa[10], P[10][di], Pa[11], P[11][ei], Pa[12], P[12][fi]);   // data for histogram of efforts
				  // plot histograms      
				  sprintf(title, "%s= %g, %s= %g, %s= %g, %s= %.g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g, %s= %g", Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][ai], Pa[8], P[8][bi], Pa[9], P[9][ci], Pa[10], P[10][di], Pa[11], P[11][ei], Pa[12], P[12][fi]);
				  //sprintf(str, "sum_%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3], P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[7], P[7][ai], Pa[8], P[8][bi], Pa[9], P[9][ci], Pa[10], P[10][di], Pa[11], P[11][ei], Pa[12], P[12][fi], Pa[13], P[13][gi], Pa[14], P[14][hi]);
				  sprintf(str, "sum_%s%04.1f%s%04.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f%s%.1f_%d_v%d%d%d",Pa[0], P[0][ri], Pa[1], P[1][si], Pa[2], P[2][ti], Pa[3],  P[3][ui], Pa[4], P[4][vi], Pa[5], P[5][wi], Pa[6], P[6][xi], Pa[10], P[10][di], Pa[11], P[11][ei], Pa[12], P[12][fi], cnt, (int)(P[7][ai]*10) , (int)(P[8][bi]*10), (int)(P[9][ci]*10));
				  
				  hist_sum(edata, ni, title, str, ++wid);   
				} 
	      
#endif
			      //}  // hi
			    //}  // gi
			  }  // fi
			}  // ei
		      }  // di
		    }  // ci
		  }  // bi
		}  // ai
	      }  // xi
	    }  // wi
	  }  // vi
	}  // ui
      }  // ti
    }  // si
  }  // ri   
  
}
  
/** Usage:
    compile : gcc -o xsum xsum.c -lm
    run     : ./xsum
**/
 
 
