// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

#define LAST 500
#define SKIP 1

#define PIB_MFA 1

/* sets of parameters to be simulated on */
int Vsim[] = {1};
double *Vc, *Vl;
double Vc1[]     = {0.01, 0.24, 0.0};
double Vc2[]     = {0.0, 0.01, 0.24};
double Vc3[]     = {0.0, 0.01, 0.24};

double Vl1[]     = {0.01, 0.24, 0.0};
double Vl2[]     = {0.01, 0.24, 0.0};
double Vl3[]     = {0.0, 0.01, 0.48};

int n[]        = {4, 8, 16};
double b[]     = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
int K[]     = {4};

double delta[] = {0.05, 0.1, 0.2};
double k[]     = {0.125, 0.25, 0.5};
double Theta[] = {0, 0.1, 0.2};

double cx[]    = {1};
double cy[]    = {0.125};
double cz[]    = {1};

double X0[]    = {2, 4, 8};    // X0
double z_0[]   = {1.0/16.0, 1.0/8.0, 1.0/4.0};          // z0
double e[]     = {1, 2, 4};

/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 10;  

int      G         = 500;   
int      H         = 1;
int      T         = 1500;   
double   Lambda    = 25; 
double   m         = 0.5;

double   Sigma     = 0.1; 
double   Sigma_t   = 0.1; 

/* end of other basic parameters */
int i[24], is[24];
int _n, _K; 
double _b, _delta, _k, _Theta, _cx, _cy, _cz, _z_0, _X0, _e;


void prep_file_in(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "g%02dn%02dK%02db%0.2fk%0.2fdelta%.2fTheta%0.2fcx%.2fcy%.2fcz%.2fe%.2fX0%.2fz0%.2fvc1%.2fvc2%.2fvc3%.2fvl1%.2fvl2%.2fvl3%.2f%s", G, _n, _K, _b, _k, _delta, _Theta, _cx, _cy, _cz, _e, _X0, _z_0, Vc[0], Vc[1], Vc[2], Vl[0], Vl[1], Vl[2], apndStr);
}


void write_param_data(FILE *fp, int vsim)
{
  fprintf(fp, "%g \t %1d \t %02d \t %0.2f \t %0.2f \t %0.2f \t %0.2f \t %0.2f \t %0.2f \t %0.2f \t %0.2f \t %0.2f \t %0.2f \t", Lambda, vsim, _n, _b, _k, _delta, _Theta, _cx, _cy, _cz, _e, _X0, _z_0);
}



void prep_data_file(FILE *fpox, int vsim)
{
  double x, y, z, P, wc, wl;
  char str[1000], fin[500];
  FILE *fpix;    
  // // create files (B, x and X0), (B, y and X0), and (B, P and X0) and payofss Ws  
  for(i[0] = 0; i[0] < is[0]; i[0]++){  // n    
    _n = n[i[0]];	
    for(i[1] = 0; i[1] < is[1]; i[1]++){  //b
      _b = b[i[1]];
      for(i[2] = 0; i[2] < is[2]; i[2]++){  // K
	_K = K[i[2]];
	for(i[3] = 0; i[3] < is[3]; i[3]++){  // delta
	  _delta = delta[i[3]];
	  for(i[4] = 0; i[4] < is[4]; i[4]++){  // k
	    _k = k[i[4]];
	    for(i[5] = 0; i[5] < is[5]; i[5]++){  //  Theta
	      _Theta = Theta[i[5]];
	      for(i[6] = 0; i[6] < is[6]; i[6]++){  // cx
		_cx = cx[i[6]];
		for(i[7] = 0; i[7] < is[7]; i[7]++){  // cy
		  _cy = cy[i[7]];
		  for(i[8] = 0; i[8] < is[8]; i[8]++){  // cz
		    _cz = cz[i[8]];		
		      for(i[10] = 0; i[10] < is[10]; i[10]++){  // z_0
			_z_0 = z_0[i[10]];		    
			for(i[11] = 0; i[11] < is[11]; i[11]++){  // e
			  _e = e[i[11]];     
			    for(i[9] = 0; i[9] < is[9]; i[9]++)  // X0						
			    {
			      _X0 = X0[i[9]];			
			      write_param_data(fpox, vsim);
			      // read from input files
			      prep_file_in(fin, "std.dat");
			      if((fpix = fopen(fin, "r"))){
				if(fgets(str, 1000, fpix) == NULL) printf("fgets read error!\n"); // read 1st line header 				
				// read data x_v, y_v and z_v for avg. std dev of 10 overall runs
				if(fscanf(fpix, "%lf %lf %lf", &x, &y, &z) == 3){
				  // write data
				  fprintf(fpox, "%.4f\t%.4f\t%.4f\t", x, y, z);				  
				}
				// read data x_v, y_v and z_v for std dev of 10 average values
				if(fscanf(fpix, "%lf %lf %lf", &x, &y, &z) == 3){
				  // write data
				  fprintf(fpox, "%.4f\t%.4f\t%.4f\n", x, y, z);				  
				}
				fclose(fpix);
			      }
			      else{
				printf("Can't open %s!\n", fin);    
				exit(1);
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

int main()
{
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
  is[10]  = sizeof(z_0)/sizeof(double);
  is[11]  = sizeof(e)/sizeof(double);
  
  int j;
  FILE *fpox;
  char fout[500];
  for(j = 0; j < sizeof(Vsim)/sizeof(int); j++){
    if(Vsim[j] == 1){ Vc = Vc1; Vl = Vl1;}
    else if(Vsim[j] == 2){ Vc = Vc2; Vl = Vl2;}
    else if(Vsim[j] == 3){ Vc = Vc3; Vl = Vl3;}  
    if(Vsim[j] == 3){
      if(PIB_MFA == 1){
	sprintf(fout, "stdout_%d_pib.dat", Vsim[j]);
      }
      else if(PIB_MFA == 2){
	sprintf(fout, "stdout_%d_mfa.dat", Vsim[j]);
      }
    }
    else{
      sprintf(fout, "stdout_%d.dat", Vsim[j]);
    }
    fpox = fopen(fout, "w"); 
    // write headers
    fprintf(fpox, "Lambda \toption\t n \t b \t k \t delta \t Theta \t cx \t cy \t cz \t e \t X0 \t z_0 \t xv1 \t\t yv1 \t\t zv1 \t\t xv2 \t\t yv2 \t\t zv2\n");
    prep_data_file(fpox, Vsim[j]);     
    fclose(fpox);
  }
  return 0;
}




/** Usage:
    compile : gcc -o dsum dsum.c
    run     : ./dsum
**/


 

