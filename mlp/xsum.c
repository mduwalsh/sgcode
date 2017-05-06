// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

#define STU 500
#define SKIP 10

/* sets of parameters to be simulated on */
int n[]        = {4, 8, 16};

double b[]     = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
//double b[]     = {3, 4};
//double b[]     = {5, 6};
//double b[]     = {7, 8};
//double b[]     = {9, 10};
//double b[]     = {11, 12};
//double b[]     = {13, 14};
//double b[]     = {15, 16};
//double b[]     = {17, 18};
//double b[]     = {19, 20 };
double B[]     = {1};
double delta[] = {0};
double DELTA[] = {0};
double K[]     = {0};
double k[]     = {0};

double Theta_u[] = {0, 0.1, 0.2};
double Eta_u[]   = {0.0, 0.1, 0.2};
double Theta_d[] = {1};
double Eta_d[]   = {1};

double cx[]    = {1};
double cy[]    = {1};
double cz[]    = {1};


double V1[]    = {0.0};
double V2[]    = {0.};
double V3[]    = {0.25};
double x_0[]   = {1};    // x0
double y_0[]   = {0.0625, 0.125, 0.25};          // y0
double z_0[]   = {0.0625, 0.125, 0.25};          // z0
double X0[]    = {2, 4, 8};
double Y0[]    = {3.2};
double e[]     = {1, 2, 4};
double E[]     = {1, 2, 4};

/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 10;  

int      G         = 8;   
int      H         = 100;
int      T         = 1500;   
int      L         = 5;
double   m         = 0.1;

double   Sigma     = 0.1; 
double   Sigma_t   = 0.1; 

/* end of other basic parameters */
int i[24], is[24];
int _n; 
double _b, _B, _delta, _DELTA, _K, _k, _Theta_u, _Theta_d, _Eta_u, _Eta_d, _cx, _cy, _cz, v1, v2, v3, _x_0, _y_0, _z_0, _X0, _Y0, _e, _E;

void prep_file_in(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "g%02dn%02db%0.2fB%.2ftu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fv1%.2fv2%.2fv3%.2fY0%.2fX0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s", G, _n, _b, _B, _Theta_u, _Theta_d, _Eta_u, _Eta_d, _cx, _cy, _cz, v1, v2, v3, _Y0, _X0, _e, _E, _x_0, _y_0, _z_0, apndStr);
}

void prep_file_out(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "n%02dtu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fv1%.2fv2%.2fv3%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s",  _n, _Theta_u, _Theta_d, _Eta_u, _Eta_d, _cx, _cy, _cz, v1, v2, v3, _e, _E, _x_0, _y_0, _z_0, apndStr);
}

void prep_file_plot(char *fname, char *prpndStr, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "%seu%.2fcx%.2fcy%.2fcz%.2fv1%.2fv2%.2fv3%.2fE%.2fx0%.2fz0%.2fe%.2fy0%.2f%s", prpndStr, _Eta_u,_cx, _cy, _cz, v1, v2, v3,_E, _x_0, _z_0, _e, _y_0, apndStr);
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
    else if(strcmp(type, "Q") == 0)
      sprintf(str, "_Qdata.dat");
    else if(strcmp(type, "W_c") == 0)
      sprintf(str, "_wcdata.dat");
    else if(strcmp(type, "W_l") == 0)
      sprintf(str, "_wldata.dat");
    else if(strcmp(type, "W_cf") == 0)
      sprintf(str, "_wcfdata.dat");
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
  char fout[300];
  FILE *gp;  
  
  for(i[2] = 0; i[2] < is[2]; i[2]++){  // B
    //_B = B[i[2]];
    for(i[3] = 0; i[3] < is[3]; i[3]++){  // delta
      _delta = delta[i[3]];
      for(i[4] = 0; i[4] < is[4]; i[4]++){  // DELTA
	_DELTA = DELTA[i[4]];
	for(i[5] = 0; i[5] < is[5]; i[5]++){  //  K
	  _K = K[i[5]];
	  for(i[6] = 0; i[6] < is[6]; i[6]++){  // k
	    _k = k[i[6]];
	    for(i[8] = 0; i[8] < is[8]; i[8]++){  // Theta_d
	      //_Theta_d = Theta_d[i[8]];
	      for(i[9] = 0; i[9] < is[9]; i[9]++){  // Eta_u
		_Eta_u = Eta_u[i[9]];
		for(i[10] = 0; i[10] < is[10]; i[10]++){  // Eta_d
		  //_Eta_d = Eta_d[i[10]];
		  _Eta_d = _Eta_u;
		  for(i[11] = 0; i[11] < is[11]; i[11]++){  // cx
		    _cx = cx[i[11]];
		    for(i[12] = 0; i[12] < is[12]; i[12]++){  // cy
		      _cy = cy[i[12]];
		      for(i[13] = 0; i[13] < is[13]; i[13]++){  // cz
			_cz = cz[i[13]];
			for(i[14] = 0; i[14] < is[14]; i[14]++){  // V1
			  v1 = V1[i[14]];
			  for(i[15] = 0; i[15] < is[15]; i[15]++){  // V2
			    v2 = V2[i[15]];
			    for(i[16] = 0; i[16] < is[16]; i[16]++){  // V3
			      v3 = V3[i[16]];
			      for(i[17] = 0; i[17] < is[17]; i[17]++){  // x_0
				_x_0 = x_0[i[17]];
				for(i[18] = 0; i[18] < is[18]; i[18]++){  // y_0
				  _y_0 = y_0[i[18]];
				  for(i[19] = 0; i[19] < is[19]; i[19]++){  //z_0 	
				    _z_0 = z_0[i[19]];
				    for(i[21] = 0; i[21] < is[21]; i[21]++){  // Y0
				      //_Y0 = Y0[i[21]];
				      for(i[22] = 0; i[22] < is[22]; i[22]++){  // e
					_e = e[i[22]];
					for(i[23] = 0; i[23] < is[23]; i[23]++){  // E 
					  _E = E[i[23]]; 
					  gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
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
					  fprintf(gp, "set multiplot layout %d, %d title '%s' offset 0.02, 0 \n", 5, is[7], "");
					  //for x
					  // major x axis
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    _Theta_d = _Theta_u;
					    // set x axis label at top
					    fprintf(gp, "set label 1 '%s=%g' at screen %.4f,%.4f center font \",10\"\n", "{/Symbol q}", _Theta_u, (1.0/(is[7]))*(i[7]+1) - 0.5/is[7], 1.0-0.02);      
					    plot(gp, "", "X", is[7]+1);
					    fprintf(gp, "unset label 1\n");
					  }
					  //for y
					  // major x axis
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    _Theta_d = _Theta_u;
					    plot(gp, "", "Y", is[7]+1);
					  }
					  //for z
					  // major x axis
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    _Theta_d = _Theta_u;
					    plot(gp, "", "Z", is[7]+1);
					  }
					  //for P
					  // major x axis
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    _Theta_d = _Theta_u;					    				   					    
					    plot(gp, "", "P", is[7]+1);
					  }
					  //for Q
					  // major x axis
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    _Theta_d = _Theta_u;
					    // set x axis label at bottom
					    fprintf(gp, "set label 2 '%s' at screen %.4f,%.4f center font \",12\"\n", "b",  (1.0/(is[7]))*(i[7]+1) - 0.5/is[7], 0.02);
					    plot(gp, "", "Q", is[7]+1);					   					    
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
					  fprintf(gp, "set multiplot layout %d, %d title '%s' offset 0.02, 0 \n", 3, is[7], "");
					  //for x
					  // major x axis
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    _Theta_d = _Theta_u;
					    // set x axis label at top
					    fprintf(gp, "set label 1 '%s=%g' at screen %.4f,%.4f center font \",10\"\n", "{/Symbol q}", _Theta_u, (1.0/(is[7]))*(i[7]+1) - 0.5/is[7], 1.0-0.02);      
					    plot(gp, "", "W_c", is[7]+1);
					    fprintf(gp, "unset label 1\n");
					  }					  
					  //for wl
					  // major x axis
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    _Theta_d = _Theta_u;
					    // set x axis label at bottom
					    plot(gp, "", "W_l", is[7]+1);					   					    
					    fprintf(gp, "unset label 2\n");
					  }
					  //for wcf
					  // major x axis
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    _Theta_d = _Theta_u;
					    // set x axis label at bottom
					    fprintf(gp, "set label 2 '%s' at screen %.4f,%.4f center font \",12\"\n", "b",  (1.0/(is[7]))*(i[7]+1) - 0.5/is[7], 0.02);
					    plot(gp, "", "W_cf", is[7]+1);					   					    
					    fprintf(gp, "unset label 2\n");
					  }
					  
					  fprintf(gp, "unset multiplot\n");
					  // multiplot ends
					  
					  
					  fflush(gp);
					  pclose(gp);
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

void prep_data_file()
{
  double x, y, z, P, Q, wc, wl, wcf;
  char str[1000], fout[500], fin[500];
  FILE *fpix, *fpiP, *fpiQ, *fpiw, *fpox, *fpoy, *fpoz, *fpoP, *fpoQ, *fpowc, *fpowl, *fpowcf;  
  
  for(i[2] = 0; i[2] < is[2]; i[2]++){  // B
    //_B = B[i[2]];
    for(i[3] = 0; i[3] < is[3]; i[3]++){  // delta
      _delta = delta[i[3]];
      for(i[4] = 0; i[4] < is[4]; i[4]++){  // DELTA
	_DELTA = DELTA[i[4]];
	for(i[5] = 0; i[5] < is[5]; i[5]++){  //  K
	  _K = K[i[5]];
	  for(i[6] = 0; i[6] < is[6]; i[6]++){  // k
	    _k = k[i[6]];
	    for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
	      _Theta_u = Theta_u[i[7]];
	      for(i[8] = 0; i[8] < is[8]; i[8]++){  // Theta_d
		//_Theta_d = Theta_d[i[8]];
		_Theta_d = _Theta_u;
		for(i[9] = 0; i[9] < is[9]; i[9]++){  // Eta_u
		  _Eta_u = Eta_u[i[9]];
		  for(i[10] = 0; i[10] < is[10]; i[10]++){  // Eta_d
		    _Eta_d = Eta_d[i[10]];
		    _Eta_d = _Eta_u;
		    for(i[11] = 0; i[11] < is[11]; i[11]++){  // cx
		      _cx = cx[i[11]];
		      for(i[12] = 0; i[12] < is[12]; i[12]++){  // cy
			_cy = cy[i[12]];
			for(i[13] = 0; i[13] < is[13]; i[13]++){  // cz
			  _cz = cz[i[13]];
			  for(i[14] = 0; i[14] < is[14]; i[14]++){  // V1
			    v1 = V1[i[14]];
			    for(i[15] = 0; i[15] < is[15]; i[15]++){  // V2
			      v2 = V2[i[15]];
			      for(i[16] = 0; i[16] < is[16]; i[16]++){  // V3
				v3 = V3[i[16]];
				for(i[17] = 0; i[17] < is[17]; i[17]++){  // x_0
				  _x_0 = x_0[i[17]];
				  for(i[18] = 0; i[18] < is[18]; i[18]++){  // y_0
				    _y_0 = y_0[i[18]];
				    for(i[19] = 0; i[19] < is[19]; i[19]++){  //z_0 	
				      _z_0 = z_0[i[19]];
				      for(i[21] = 0; i[21] < is[21]; i[21]++){  // Y0
					_Y0 = Y0[i[21]];					
					for(i[22] = 0; i[22] < is[22]; i[22]++){  // e
					  _e = e[i[22]];
					  for(i[23] = 0; i[23] < is[23]; i[23]++){  // E 
					    _E = E[i[23]]; 
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
					      prep_file_out(fout, "_Qdata.dat");					      
					      fpoQ = fopen(fout, "w");
					      prep_file_out(fout, "_wcdata.dat");
					      fpowc = fopen(fout, "w");
					      prep_file_out(fout, "_wldata.dat");					      
					      fpowl = fopen(fout, "w");
					      prep_file_out(fout, "_wcfdata.dat");
					      fpowcf = fopen(fout, "w");
					      // write headers
					      fprintf(fpox, "b "); fprintf(fpoy, "b "); 
					      fprintf(fpoP, "b "); fprintf(fpoQ, "b "); 
					      fprintf(fpowc, "b "); fprintf(fpowl, "b ");
					      fprintf(fpoz, "b "); fprintf(fpowcf, "b ");
					      for(i[20] = 0; i[20] < is[20]; i[20]++)  // X0						
					      {
						fprintf(fpox, "X0=%.1f  ", X0[i[20]]);
						fprintf(fpoy, "X0=%.1f  ", X0[i[20]]);
						fprintf(fpoz, "X0=%.1f  ", X0[i[20]]);
						fprintf(fpoP, "X0=%.1f  ", X0[i[20]]);
						fprintf(fpoQ, "X0=%.1f  ", X0[i[20]]);
						fprintf(fpowc, "X0=%.1f  ", X0[i[20]]);
						fprintf(fpowl, "X0=%.1f  ", X0[i[20]]);
						fprintf(fpowcf, "X0=%.1f  ", X0[i[20]]);
					      }
					      fprintf(fpox, "\n"); fprintf(fpoy, "\n"); fprintf(fpoz, "\n"); 
					      fprintf(fpoP, "\n"); fprintf(fpoQ, "\n"); 
					      fprintf(fpowc, "\n"); fprintf(fpowl, "\n"); fprintf(fpowcf, "\n");
					      for(i[1] = 0; i[1] < is[1]; i[1]++){  //b
						_b = b[i[1]];
						_B = _b;
						fprintf(fpox, "%.1lf ", _b);
						fprintf(fpoy, "%.1lf ", _b);
						fprintf(fpoz, "%.1lf ", _b);
						fprintf(fpoP, "%.1lf ", _b);
						fprintf(fpoQ, "%.1lf ", _b);
						fprintf(fpowc, "%.1lf ", _b);
						fprintf(fpowl, "%.1lf ", _b);
						fprintf(fpowcf, "%.1lf ", _b);
						for(i[20] = 0; i[20] < is[20]; i[20]++)  // X0						
						{
						  _X0 = X0[i[20]];
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
						      fprintf(fpoP, "%.4f  ", P/(STU/SKIP +1));						      						
						    }
						  }
						  else{
						    printf("Can't open %s!\n", fin);
						  }
						  fclose(fpiP);
						  
						  prep_file_in(fin, "Qsum.dat");
						  if((fpiQ = fopen(fin, "r"))){
						   if(fgets(str, 1000, fpiQ) == NULL) printf("fgets read error!\n");  // read 1st line header 
						    // read data P
						    if(fscanf(fpiQ, "%lf", &Q) == 1){
						      // write data
						      fprintf(fpoQ, "%.4f  ", Q);						      						
						    }
						  }
						  else{
						    printf("Can't open %s!\n", fin);
						  }
						  fclose(fpiQ);
						  
						  // payoffs w
						  prep_file_in(fin, "psum.dat");
						  if((fpiw = fopen(fin, "r"))){
						    if(fgets(str, 1000, fpiw) == NULL) printf("fgets read error!\n"); // read 1st line header 
						    // read data wc and wl
						    if(fscanf(fpiw, "%lf %lf %lf", &wc, &wl, &wcf) == 3){
						      // write data
						      fprintf(fpowc, "%.4f  ", wc);
						      fprintf(fpowl, "%.4f  ", wl);	
						      fprintf(fpowcf, "%.4f  ", wcf);
						    }
						    fclose(fpiw);
						  }
						  else{
						    printf("Can't open %s!\n", fin);    
						  }
						  
						}
						fprintf(fpox, "\n"); fprintf(fpoy, "\n"); fprintf(fpoz, "\n"); 
						fprintf(fpoP, "\n"); fprintf(fpoQ, "\n"); 
						fprintf(fpowc, "\n"); fprintf(fpowl, "\n"); fprintf(fpowcf, "\n");
					      }
					      fclose(fpox); fclose(fpoy); fclose(fpoz); fclose(fpoP); fclose(fpoQ); fclose(fpowc); fclose(fpowl); fclose(fpowcf);
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

int main()
{
  is[0]  = sizeof(n)/sizeof(int);
  is[1]  = sizeof(b)/sizeof(double);
  is[2]  = sizeof(B)/sizeof(double);
  is[3]  = sizeof(delta)/sizeof(double);
  is[4]  = sizeof(DELTA)/sizeof(double);
  is[5]  = sizeof(K)/sizeof(double);
  is[6]  = sizeof(k)/sizeof(double);
  is[7]  = sizeof(Theta_u)/sizeof(double);
  is[8]  = sizeof(Theta_d)/sizeof(double);
  is[9]  = sizeof(Eta_u)/sizeof(double);
  is[10]  = sizeof(Eta_d)/sizeof(double);
  is[11]  = sizeof(cx)/sizeof(double);
  is[12]  = sizeof(cy)/sizeof(double);
  is[13]  = sizeof(cz)/sizeof(double);
  is[14]  = sizeof(V1)/sizeof(double);
  is[15]  = sizeof(V2)/sizeof(double);
  is[16]  = sizeof(V3)/sizeof(double);
  is[17]  = sizeof(x_0)/sizeof(double);
  is[18]  = sizeof(y_0)/sizeof(double);
  is[19]  = sizeof(z_0)/sizeof(double);
  is[20]  = sizeof(X0)/sizeof(double);
  is[21]  = sizeof(Y0)/sizeof(double);
  is[22]  = sizeof(e)/sizeof(double);
  is[23]  = sizeof(E)/sizeof(double);
  
  prep_data_file();
  
  plot_data_points();
  
  return 0;
}




/** Usage:
    compile : gcc -o mljob mljob.c
    run     : ./mljob
**/

/*
 * createConfig(j, n[i[0]], b[i[1]], b[i[1]], delta[i[3]], DELTA[i[4]], K[i[5]], k[i[6]], Theta_u[i[7]], Theta_d[i[8]], Eta_u[i[9]], Eta_d[i[10]], cx[i[11]], cy[i[12]], cz[i[13]], V1[i[14]], V2[i[15]], V3[i[16]], x_0[i[17]], y_0[i[18]], z_0[i[19]], X0[i[20]], X0[i[20]], e[i[22]], E[i[23]]);
						  createScript(j, n[i[0]], b[i[1]], b[i[1]], delta[i[3]], DELTA[i[4]], K[i[5]], k[i[6]], Theta_u[i[7]], Theta_d[i[8]], Eta_u[i[9]], Eta_d[i[10]], cx[i[11]], cy[i[12]], cz[i[13]], V1[i[14]], V2[i[15]], V3[i[16]], x_0[i[17]], y_0[i[18]], z_0[i[19]], X0[i[20]], X0[i[20]], e[i[22]], E[i[23]]);   
						  sprintf(scall, "qsub mlpscript.%d.sh", j);
						  sc = system(scall);
						  if(sc) printf("Error submitting jobs!!\n");
*/
