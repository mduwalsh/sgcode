// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

/* sets of parameters to be simulated on */
int Vsim[] = {3};         // option 1, 2, 3; strategy update  
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

double delta[] = {0.1};
double k[]     = {0.25};
double Theta[] = {1, 2, 4};

double cx[]    = {1};
double cy[]    = {0.125};
double cz[]    = {1};

double X0[]    = {4, 8, 16};    // X0
double s0[]   = {0.25};          // s0
double s1[]     = {0.25};
double e[]     = {1};
double z_0[]   = {1};

/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 10;  

int      G         = 500;   
int      H         = 1;
int      T         = 1500;   
double   Lambda    = 4; 
double   m         = 0.5;

double   Sigma     = 0.1; 
double   Sigma_t   = 0.1; 

/* end of other basic parameters */
int i[24], is[24];
int _n, _K; 
double _b, _delta, _k, _Theta, _cx, _cy, _cz, _X0, _s0, _s1, _e, _z_0;


void prep_file_in(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "g%02dn%02dK%02db%0.2fk%0.2fdelta%.2fTheta%0.2fcx%.2fcy%.2fcz%.2fX0%.2fe%.2fz0%.2fs0%.2fs1%.2fvc1%.2fvc2%.2fvc3%.2fvl1%.2fvl2%.2fvl3%.2f%s", G, _n, _K, _b, _k, _delta, _Theta, _cx, _cy, _cz, _X0, _e, _z_0, _s0, _s1, Vc[0], Vc[1], Vc[2], Vl[0], Vl[1], Vl[2], apndStr);
}

void prep_file_out(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "n%02dK%02dth%.2fk%.2fdl%.2fcx%.2fcy%.2fcz%.2fe%.2fz0%.2fs0%.2fs1%.2fvc1%.2fvc2%.2fvc3%.2fvl1%.2fvl2%.2fvl3%.2f%s",  _n, _K, _Theta, _k, _delta, _cx, _cy, _cz, _e, _z_0, _s0, _s1, Vc[0], Vc[1], Vc[2], Vl[0], Vl[1], Vl[2], apndStr);
}

void prep_file_plot(char *fname, char *prpndStr, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "%sK%02dcx%.2fcy%.2fcz%.2fvc1%.2fvc2%.2fvc3%.2fvl1%.2fvl2%.2fvl3%.2fk%.2fdl%.2fe%.2fz0%.2fs0%.2fs1%.2f%s", prpndStr, _K, _cx, _cy, _cz, Vc[0], Vc[1], Vc[2], Vl[0], Vl[1], Vl[2], _k, _delta, _e, _z_0, _s0, _s1, apndStr);
}

void plot(FILE *gp, char *title, char *type, int data_cols)
{      
    // one plot of graph
    char str[30], fdata[300];
    if(strcmp(type, "X") == 0)
      sprintf(str, "_xdata.dat");
    else if(strcmp(type, "Y") == 0)
      sprintf(str, "_ydata.dat");
    else if(strcmp(type, "Z") == 0)
      sprintf(str, "_zdata.dat");
    else if(strcmp(type, "P") == 0)
      sprintf(str, "_Pdata.dat");    
    else if(strcmp(type, "W_c") == 0)
      sprintf(str, "_wcdata.dat");
    else if(strcmp(type, "W_l") == 0)
      sprintf(str, "_wldata.dat");    
    fprintf(gp, "set title ''\n");
    fprintf(gp, "set format y \"%%.2f\"\n");
    fprintf(gp, "set ylabel '%s' font \",12\" \n", type);
    fprintf(gp, "set xlabel '%s'\n", "");
    fprintf(gp, "plot \\\n"); 
    for(i[0] = 0; i[0] < is[0]; i[0]++){  // n    
      _n = n[i[0]];
      prep_file_out(fdata, str);
      if(i[0] == 0){
	fprintf(gp, "for [col=2:%d] '%s' using 1:col with lines lw %d lc col-1 title columnheader \\\n", data_cols, fdata, i[0]+1);
      }
      else{
	fprintf(gp, ", for [col=2:%d] '%s' using 1:col with lines lw %d lc col-1 notitle \\\n", data_cols, fdata, i[0]+1);
      }
      					      
    }
    fprintf(gp, "\n");    
}

void plot_data_points()
{
  char fout[300], fdata[300];
  int count = 0;
  FILE *gp;  
  gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  fprintf(gp, "xmax = 0.1\n");
  fprintf(gp, "ymax = 0.1\n");
  fprintf(gp, "zmax = 0.1\n");
  fprintf(gp, "wcmax = 0.1 \n");
  fprintf(gp, "wlmax = 0.1 \n");
  
  for(i[2] = 0; i[2] < is[2]; i[2]++){  // K
    _K = K[i[2]];
    for(i[3] = 0; i[3] < is[3]; i[3]++){  // delta
      _delta = delta[i[3]];
      for(i[4] = 0; i[4] < is[4]; i[4]++){  // k
	_k = k[i[4]];	
	for(i[6] = 0; i[6] < is[6]; i[6]++){  // cx
	  _cx = cx[i[6]];
	  for(i[7] = 0; i[7] < is[7]; i[7]++){  // cy
	    _cy = cy[i[7]];
	    for(i[8] = 0; i[8] < is[8]; i[8]++){  // cz
	      _cz = cz[i[8]];		
	      for(i[10] = 0; i[10] < is[10]; i[10]++){  // s1
		_s1 = s1[i[10]];		    
		for(i[11] = 0; i[11] < is[11]; i[11]++){  // s0
		  _s0 = s0[i[11]];  		  		    
		  for(i[12] = 0; i[12] < is[12]; i[12]++){  // e
		    _e = e[i[12]];		    
		    for(i[13] = 0; i[13] < is[13]; i[13]++){  // z0
		      _z_0 = z_0[i[13]];  		
		      for(i[5] = 0; i[5] < is[5]; i[5]++){  //  Theta
			_Theta = Theta[i[5]];
			for(i[0] = 0; i[0] < is[0]; i[0]++){  // n    
			  _n = n[i[0]];		      
			  for(i[9] = 0; i[9] < is[9]; i[9]++){  // X0						  		
			    //for x
			    prep_file_out(fdata, "_xdata.dat");
			    fprintf(gp, "stats '%s' using %d prefix 'A%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
			    fprintf(gp, "if(A%d%d_max > xmax) {xmax = A%d%d_max}\n", i[9], count, i[9], count);
			    //for y
			    prep_file_out(fdata, "_ydata.dat");
			    fprintf(gp, "stats '%s' using %d prefix 'B%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
			    fprintf(gp, "if(B%d%d_max > ymax) {ymax = B%d%d_max}\n", i[9], count, i[9], count);
			    //for z
			    prep_file_out(fdata, "_zdata.dat");
			    fprintf(gp, "stats '%s' using %d prefix 'C%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
			    fprintf(gp, "if(C%d%d_max > zmax) {zmax = C%d%d_max}\n", i[9], count, i[9], count);
			    //for wc
			    prep_file_out(fdata, "_wcdata.dat");
			    fprintf(gp, "stats '%s' using %d prefix 'D%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
			    fprintf(gp, "if(D%d%d_max > wcmax) {wcmax = D%d%d_max}\n", i[9], count, i[9], count);
			    //for wl
			    prep_file_out(fdata, "_wldata.dat");
			    fprintf(gp, "stats '%s' using %d prefix 'E%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
			    fprintf(gp, "if(E%d%d_max > wlmax) {wlmax = E%d%d_max}\n", i[9], count, i[9], count);
			  }
			  count++;
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
  for(i[2] = 0; i[2] < is[2]; i[2]++){  // K
    _K = K[i[2]];
    for(i[3] = 0; i[3] < is[3]; i[3]++){  // delta
      _delta = delta[i[3]];
      for(i[4] = 0; i[4] < is[4]; i[4]++){  // k
	_k = k[i[4]];	
	for(i[6] = 0; i[6] < is[6]; i[6]++){  // cx
	  _cx = cx[i[6]];
	  for(i[7] = 0; i[7] < is[7]; i[7]++){  // cy
	    _cy = cy[i[7]];
	    for(i[8] = 0; i[8] < is[8]; i[8]++){  // cz
	      _cz = cz[i[8]];		
	      for(i[10] = 0; i[10] < is[10]; i[10]++){  // s1
		_s1 = s1[i[10]];		    
		for(i[11] = 0; i[11] < is[11]; i[11]++){  // s0
		  _s0 = s0[i[11]];  		 
		  for(i[12] = 0; i[12] < is[12]; i[12]++){  // e
		    _e = e[i[12]];		    
		    for(i[13] = 0; i[13] < is[13]; i[13]++){  // z0
		      _z_0 = z_0[i[13]];  		
		      fprintf(gp, "set term pngcairo size 1600,900 enhanced color solid font \"Helvetica,8\" \n");
		      //multiplot begins for efforts					  
		      prep_file_plot(fout, "ef", ".png");
		      fprintf(gp, "set output '%s' \n", fout);
		      fprintf(gp, "set key autotitle columnheader\n");
		      fprintf(gp, "set key inside left top \n");
		      fprintf(gp, "set bmargin 3 \n");
		      fprintf(gp, "set tmargin 1 \n");
		      fprintf(gp, "set lmargin 6 \n");
		      fprintf(gp, "set rmargin 12 \n");
		      fprintf(gp, "set multiplot layout %d, %d title '%s' offset 0.02, 0 \n", 4, is[5], "");		  
		      //for x		  
		      // major x axis
		      fprintf(gp, "set yrange [0:xmax] \n");
		      for(i[5] = 0; i[5] < is[5]; i[5]++){  //  Theta
			_Theta = Theta[i[5]];
			// set x axis label at top
			fprintf(gp, "set label 1 '%s=%g' at screen %.4f,%.4f center font \",10\"\n", "a", _Theta, (1.0/(is[5]))*(i[5]+1) - 0.5/is[5], 1.0-0.02);      
			plot(gp, "", "X", is[5]+1);
			fprintf(gp, "unset label 1\n");
		      }
		      //for y
		      fprintf(gp, "set yrange [0:ymax] \n");
		      // major x axis
		      for(i[5] = 0; i[5] < is[5]; i[5]++){  //  Theta
			_Theta = Theta[i[5]];
			plot(gp, "", "Y", is[5]+1);
		      }
		      //for z
		      fprintf(gp, "set yrange [0:zmax] \n");
		      // major x axis
		      for(i[5] = 0; i[5] < is[5]; i[5]++){  //  Theta
			_Theta = Theta[i[5]];
			plot(gp, "", "Z", is[5]+1);
		      }
		      //for P
		      fprintf(gp, "set yrange [0:1.1] \n");
		      // major x axis
		      for(i[5] = 0; i[5] < is[5]; i[5]++){  //  Theta
			_Theta = Theta[i[5]];
			fprintf(gp, "set label 2 '%s' at screen %.4f,%.4f center font \",12\"\n", "b",  (1.0/(is[5]))*(i[5]+1) - 0.5/is[5], 0.02);
			plot(gp, "", "P", is[5]+1);
			fprintf(gp, "unset label 2\n");
		      }
		      
		      fprintf(gp, "unset multiplot\n");
		      // multiplot ends
		      
		      
		      prep_file_plot(fout, "py", ".png");
		      fprintf(gp, "set output '%s' \n", fout);
		      fprintf(gp, "set key autotitle columnheader\n");
		      fprintf(gp, "set key inside left top \n");
		      fprintf(gp, "set bmargin 3 \n");
		      fprintf(gp, "set tmargin 1 \n");
		      fprintf(gp, "set lmargin 6 \n");
		      fprintf(gp, "set rmargin 12 \n");
		      fprintf(gp, "set multiplot layout %d, %d title '%s' offset 0.02, 0 \n", 2, is[5], "");
		      fprintf(gp, "set autoscale y \n");
		      //for wc
		      fprintf(gp, "set yrange [:wcmax] \n");
		      // major x axis
		      for(i[5] = 0; i[5] < is[5]; i[5]++){  //  Theta
			_Theta = Theta[i[5]];
			// set x axis label at top
			fprintf(gp, "set label 1 '%s=%g' at screen %.4f,%.4f center font \",10\"\n", "a", _Theta, (1.0/(is[5]))*(i[5]+1) - 0.5/is[5], 1.0-0.02);      
			plot(gp, "", "W_c", is[5]+1);
			fprintf(gp, "unset label 1\n");
		      }					  
		      //for wl
		      fprintf(gp, "set yrange [:wlmax] \n");
		      // major x axis
		      for(i[5] = 0; i[5] < is[5]; i[5]++){  //  Theta
			_Theta = Theta[i[5]];
			// set x axis label at bottom
			fprintf(gp, "set label 2 '%s' at screen %.4f,%.4f center font \",12\"\n", "b",  (1.0/(is[5]))*(i[5]+1) - 0.5/is[5], 0.02);
			plot(gp, "", "W_l", is[5]+1);					   					    
			fprintf(gp, "unset label 2\n");
		      }			
		      
		      fprintf(gp, "unset multiplot\n");
		      // multiplot ends		  
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
  fflush(gp);
  pclose(gp);
}

void prep_data_file()
{
  double x, y, z, P, wc, wl;
  char str[1000], fout[500], fin[500];
  FILE *fpix, *fpiP, *fpiw, *fpox, *fpoy, *fpoz, *fpoP, *fpowc, *fpowl;    
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
		    for(i[10] = 0; i[10] < is[10]; i[10]++){  // s1
		      _s1 = s1[i[10]];		    
		      for(i[11] = 0; i[11] < is[11]; i[11]++){  // s0
			_s0 = s0[i[11]];  
			for(i[12] = 0; i[12] < is[12]; i[12]++){  // e
			  _e = e[i[12]];		    
			  for(i[13] = 0; i[13] < is[13]; i[13]++){  // z0
			    _z_0 = z_0[i[13]];  
			    for(i[0] = 0; i[0] < is[0]; i[0]++){  // n    
			      _n = n[i[0]];
			      // // create files (B, x and X0), (B, y and X0), and (B, P and X0) and payofss Ws
			      prep_file_out(fout, "_xdata.dat");
			      fpox = fopen(fout, "w");
			      prep_file_out(fout, "_ydata.dat");
			      fpoy = fopen(fout, "w");
			      prep_file_out(fout, "_zdata.dat");
			      fpoz = fopen(fout, "w");
			      prep_file_out(fout, "_Pdata.dat");					      
			      fpoP = fopen(fout, "w");			
			      prep_file_out(fout, "_wcdata.dat");
			      fpowc = fopen(fout, "w");
			      prep_file_out(fout, "_wldata.dat");					      
			      fpowl = fopen(fout, "w");			
			      // write headers
			      fprintf(fpox, "b "); fprintf(fpoy, "b "); fprintf(fpoz, "b "); 
			      fprintf(fpoP, "b "); 
			      fprintf(fpowc, "b "); fprintf(fpowl, "b ");			
			      for(i[9] = 0; i[9] < is[9]; i[9]++)  // X0						
			      {
				fprintf(fpox, "X0=%.1f  ", X0[i[9]]);
				fprintf(fpoy, "X0=%.1f  ", X0[i[9]]);
				fprintf(fpoz, "X0=%.1f  ", X0[i[9]]);
				fprintf(fpoP, "X0=%.1f  ", X0[i[9]]);			  
				fprintf(fpowc, "X0=%.1f  ", X0[i[9]]);
				fprintf(fpowl, "X0=%.1f  ", X0[i[9]]);			  
			      }
			      fprintf(fpox, "\n"); fprintf(fpoy, "\n"); fprintf(fpoz, "\n"); 
			      fprintf(fpoP, "\n"); 
			      fprintf(fpowc, "\n"); fprintf(fpowl, "\n"); 
			      for(i[1] = 0; i[1] < is[1]; i[1]++){  //b
				_b = b[i[1]];			  
				fprintf(fpox, "%.1lf ", _b);
				fprintf(fpoy, "%.1lf ", _b);
				fprintf(fpoz, "%.1lf ", _b);
				fprintf(fpoP, "%.1lf ", _b);			  
				fprintf(fpowc, "%.1lf ", _b);
				fprintf(fpowl, "%.1lf ", _b);			  
				for(i[9] = 0; i[9] < is[9]; i[9]++)  // X0						
				{
				  _X0 = X0[i[9]];
				  // read from input files
				  prep_file_in(fin, "xsum.dat");
				  if((fpix = fopen(fin, "r"))){
				    if(fgets(str, 1000, fpix) == NULL) printf("fgets read error!\n"); // read 1st line header 
				    // read data x and y
				    if(fscanf(fpix, "%lf %lf %lf", &x, &y, &z) == 3){
				      // write data
				      fprintf(fpox, "%.4f  ", x*_n);
				      fprintf(fpoy, "%.4f  ", y);						
				      fprintf(fpoz, "%.4f  ", z);
				    }
				    fclose(fpix);
				  }
				  else{
				    printf("Can't open %s!\n", fin);    
				  }
				  
				  prep_file_in(fin, "Psum.dat");
				  if((fpiP = fopen(fin, "r"))){
				    if(fgets(str, 1000, fpiP) == NULL) printf("fgets read error!\n");  // read 1st line header 
				    // read data P
				    if(fscanf(fpiP, "%lf", &P) == 1){
				      // write data
				      fprintf(fpoP, "%.4f  ", P);						      						
				    }
				  }
				  else{
				    printf("Can't open %s!\n", fin);
				  }
				  fclose(fpiP);    
				  
				  // payoffs w
				  prep_file_in(fin, "psum.dat");
				  if((fpiw = fopen(fin, "r"))){
				    if(fgets(str, 1000, fpiw) == NULL) printf("fgets read error!\n"); // read 1st line header 
				    // read data wc and wl
				    if(fscanf(fpiw, "%lf %lf", &wc, &wl) == 2){
				      // write data
				      fprintf(fpowc, "%.4f  ", wc);
				      fprintf(fpowl, "%.4f  ", wl);					
				    }
				    fclose(fpiw);
				  }
				  else{
				    printf("Can't open %s!\n", fin);    
				  }
				  
				}
				fprintf(fpox, "\n"); fprintf(fpoy, "\n"); fprintf(fpoz, "\n"); 
				fprintf(fpoP, "\n"); 
				fprintf(fpowc, "\n"); fprintf(fpowl, "\n");
			      }
			      fclose(fpox); fclose(fpoy); fclose(fpoz); fclose(fpoP); fclose(fpowc); fclose(fpowl);
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
  is[10]  = sizeof(s1)/sizeof(double);
  is[11]  = sizeof(s0)/sizeof(double);
  is[12]  = sizeof(e)/sizeof(double);
  is[13]  = sizeof(z_0)/sizeof(double);
  
  int j;
  for(j = 0; j < sizeof(Vsim)/sizeof(int); j++){
    if(Vsim[j] == 1){ Vc = Vc1; Vl = Vl1;}
    else if(Vsim[j] == 2){ Vc = Vc2; Vl = Vl2;}
    else if(Vsim[j] == 3){ Vc = Vc3; Vl = Vl3;}  
    prep_data_file();    
    plot_data_points();
  }
  return 0;
}




/** Usage:
    compile : gcc -o mljob mljob.c
    run     : ./mljob
**/


 
