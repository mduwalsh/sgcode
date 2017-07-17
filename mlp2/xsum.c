// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>

#define STU 500
#define SKIP 10

/* sets of parameters to be simulated on */
int Vop[]      = {1};
int n[]        = {4, 8, 16};
double b[]     = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
double B[]     = {1};
double delta[] = {0.05, 0.1, 0.2};
double DELTA[] = {0};
double K[]     = {0};
double k[]     = {0.125, 0.25, 0.5};

double Theta_u[] = {0, 0.1, 0.2};
double Eta_u[]   = {0, 0.1, 0.2};
double Theta_d[] = {1};
double Eta_d[]   = {1};

double cx[]    = {1};
double cy[]    = {1};
double cz[]    = {1};

double x_0[]   = {1};    // x0
double y_0[]   = {1};          // y0
double z_0[]   = {1};          // z0
double X0[]    = {2, 4, 8};
double Y0[]    = {1};
double e[]     = {1};
double E[]     = {1};
double Rho[]   = {0.0, 0.25, 0.5, 0.75, 1};


/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 10;  

int      G         = 8;   
int      H         = 100;
int      T         = 1500;   
int      L         = 2;
double   m         = 0.1;

double   Sigma     = 0.1; 
double   Sigma_t   = 0.1; 

double Vc[][3] = { {0.01, 0.24, 0.0}, {0.00, 0.01, 0.24}, {0.00, 0.01, 0.24} }; // commoner strategy update method probability sets for different options (1, 2, 3)
double Vl[][3] = { {0.01, 0.24, 0.0}, {0.01, 0.24, 0.0}, {0.00, 0.01, 0.48} };  // leader strategy update method probability sets for different options (1, 2, 3)
double Vcf[][3] = { {0.01, 0.24, 0.0}, {0.01, 0.24, 0.0}, {0.00, 0.01, 0.48} };  // chief strategy update method probability sets for different options (1, 2, 3)

/* end of other basic parameters */
int i[23], is[23];
int _n, _Vop; 
double _b, _B, _delta, _DELTA, _K, _k, _Theta_u, _Theta_d, _Eta_u, _Eta_d, _cx, _cy, _cz, _x_0, _y_0, _z_0, _Rho, _X0, _Y0, _e, _E;



// data file to read data from
void prep_file_in(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "g%02dn%02db%0.2fB%.2fk%.2fdl%.2ftu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fVop%dRho%.2fY0%.2fX0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s", G, _n, _b, _B, _k, _delta, _Theta_u, _Theta_d, _Eta_u, _Eta_d, _cx, _cy, _cz, _Vop, _Rho, _X0, _X0, _e, _E, _x_0, _y_0, _z_0, apndStr);
}

// data file prepared to plot graph from
void prep_file_out(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "n%02dk%.2fdl%.2ftu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fVop%dRho%.2fY0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s",  _n, _k, _delta, _Theta_u, _Theta_d, _Eta_u, _Eta_d, _cx, _cy, _cz, _Vop, _Rho, _Y0, _e, _E, _x_0, _y_0, _z_0, apndStr);
}

void prep_file_plot(char *fname, char *prpndStr, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "%sk%.2fdl%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fVop%dRho%.2fY0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s", prpndStr, _k, _delta, _Theta_d, _Eta_u, _Eta_d, _cx, _cy, _cz, _Vop, _Rho, _Y0, _e, _E, _x_0, _y_0, _z_0, apndStr);
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
    else if(strcmp(type, "p") == 0)
      sprintf(str, "_pdata.dat");
    else if(strcmp(type, "q") == 0)
      sprintf(str, "_qdata.dat");
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

void set_max_yrange(FILE *gp)
{
  int count = 0;
  char fdata[300];
  fprintf(gp, "xmax = 0.1\n");
  fprintf(gp, "ymax = 0.1\n");
  fprintf(gp, "zmax = 0.1\n");
  fprintf(gp, "pmax = 0.1\n");
  fprintf(gp, "qmax = 0.1\n");
  fprintf(gp, "Pmax = 0.1\n");
  //fprintf(gp, "Qmax = 0.1\n");
  fprintf(gp, "wcmax = 0.1 \n");
  fprintf(gp, "wlmax = 0.1 \n");
  fprintf(gp, "wcfmax = 0.1 \n");
  for(i[0] = 0; i[0] < is[0]; i[0]++){  // n
    _n = n[i[0]];
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
		_Theta_d = Theta_d[i[8]];
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
			  for(i[14] = 0; i[14] < is[14]; i[14]++){  // Vop
			    _Vop = Vop[i[14]];
			    for(i[15] = 0; i[15] < is[15]; i[15]++){  // Rho
			      _Rho = Rho[i[15]];
			      for(i[16] = 0; i[16] < is[16]; i[16]++){  // x_0
				_x_0 = x_0[i[16]];
				for(i[17] = 0; i[17] < is[17]; i[17]++){  // y_0 
				  _y_0 = y_0[i[17]];
				  for(i[18] = 0; i[18] < is[18]; i[18]++){ // z_0 
				    _z_0 = z_0[i[18]];
				    for(i[20] = 0; i[20] < is[20]; i[20]++){  // Y0
				      _Y0 = Y0[i[20]];
				      for(i[21] = 0; i[21] < is[21]; i[21]++){  // e
					_e = e[i[21]];
					for(i[22] = 0; i[22] < is[22]; i[22]++){  // E
					  _E = E[i[22]];  
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					  
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
					    //for y
					    prep_file_out(fdata, "_pdata.dat");
					    fprintf(gp, "stats '%s' using %d prefix 'D%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
					    fprintf(gp, "if(D%d%d_max > pmax) {pmax = D%d%d_max}\n", i[9], count, i[9], count);
					    //for z
					    prep_file_out(fdata, "_qdata.dat");
					    fprintf(gp, "stats '%s' using %d prefix 'E%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
					    fprintf(gp, "if(E%d%d_max > qmax) {qmax = E%d%d_max}\n", i[9], count, i[9], count);
					    //for wc
					    prep_file_out(fdata, "_Pdata.dat");
					    fprintf(gp, "stats '%s' using %d prefix 'F%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
					    fprintf(gp, "if(F%d%d_max > Pmax) {Pmax = F%d%d_max}\n", i[9], count, i[9], count);
					    //for wl
					    prep_file_out(fdata, "_wcdata.dat");
					    fprintf(gp, "stats '%s' using %d prefix 'G%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
					    fprintf(gp, "if(G%d%d_max > wcmax) {wcmax = G%d%d_max}\n", i[9], count, i[9], count);
					    //for wc
					    prep_file_out(fdata, "_wcfdata.dat");
					    fprintf(gp, "stats '%s' using %d prefix 'H%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
					    fprintf(gp, "if(H%d%d_max > wcfmax) {wcfmax = H%d%d_max}\n", i[9], count, i[9], count);
					    //for wl
					    prep_file_out(fdata, "_wldata.dat");
					    fprintf(gp, "stats '%s' using %d prefix 'I%d%d' nooutput\n", fdata, i[9]+2, i[9], count);
					    fprintf(gp, "if(I%d%d_max > wlmax) {wlmax = I%d%d_max}\n", i[9], count, i[9], count);
					  
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
		}
	      }
	    }
	  }
	}
      }
    }
  }
}
  

void plot_data_points()
{
  char fout[300];
  FILE *gp;  
  gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  // gets and sets maximum values in yrange for each traits
  set_max_yrange(gp);
  
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
	      _Theta_d = Theta_d[i[8]];
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
			for(i[14] = 0; i[14] < is[14]; i[14]++){  // Vop
			  _Vop = Vop[i[14]];
			  for(i[15] = 0; i[15] < is[15]; i[15]++){  // Rho
			    _Rho = Rho[i[15]];
			    for(i[16] = 0; i[16] < is[16]; i[16]++){  // x_0
			      _x_0 = x_0[i[16]];
			      for(i[17] = 0; i[17] < is[17]; i[17]++){  // y_0 
				_y_0 = y_0[i[17]];
				for(i[18] = 0; i[18] < is[18]; i[18]++){ // z_0 
				  _z_0 = z_0[i[18]];
				  for(i[20] = 0; i[20] < is[20]; i[20]++){  // Y0
				    _Y0 = Y0[i[20]];
				    for(i[21] = 0; i[21] < is[21]; i[21]++){  // e
				      _e = e[i[21]];
				      for(i[22] = 0; i[22] < is[22]; i[22]++){  // E
					_E = E[i[22]]; 					
					  
					  
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
					  fprintf(gp, "set multiplot layout %d, %d title '%s' offset 0.02, 0 \n", 4, is[7], "");
					  //for x
					  // major x axis
					  fprintf(gp, "set yrange [0:xmax] \n");
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    //_Theta_d = _Theta_u;
					    // set x axis label at top
					    fprintf(gp, "set label 1 '%s=%g' at screen %.4f,%.4f center font \",10\"\n", "{/Symbol q}", _Theta_u, (1.0/(is[7]))*(i[7]+1) - 0.5/is[7], 1.0-0.02);      
					    plot(gp, "", "X", is[7]+1);
					    fprintf(gp, "unset label 1\n");
					  }
					  //for y
					  // major x axis
					  fprintf(gp, "set yrange [0:ymax] \n");
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    //_Theta_d = _Theta_u;
					    plot(gp, "", "Y", is[7]+1);
					  }
					  //for z
					  // major x axis
					  /*fprintf(gp, "set yrange [0:zmax] \n");
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    //_Theta_d = _Theta_u;
					    plot(gp, "", "Z", is[7]+1);
					  }*/
					  //for p
					  // major x axis
					  fprintf(gp, "set yrange [0:pmax] \n");
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    //_Theta_d = _Theta_u;
					    plot(gp, "", "p", is[7]+1);
					  }
					  //for q
					  // major x axis
					  /*fprintf(gp, "set yrange [0:qmax] \n");
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    //_Theta_d = _Theta_u;
					    plot(gp, "", "q", is[7]+1);
					  }*/
					  //for P
					  // major x axis
					  fprintf(gp, "set yrange [0:Pmax] \n");
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    //_Theta_d = _Theta_u;					    				   					    
					    plot(gp, "", "P", is[7]+1);
					  }

					  //for Q
					  // major x axis
					  fprintf(gp, "set yrange [0:Qmax] \n");
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
					  fprintf(gp, "set multiplot layout %d, %d title '%s' offset 0.02, 0 \n", 2, is[7], "");
					  //for x
					  // major x axis
					  fprintf(gp, "set yrange [0:wcmax] \n");
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    //_Theta_d = _Theta_u;
					    // set x axis label at top
					    fprintf(gp, "set label 1 '%s=%g' at screen %.4f,%.4f center font \",10\"\n", "{/Symbol q}", _Theta_u, (1.0/(is[7]))*(i[7]+1) - 0.5/is[7], 1.0-0.02);      
					    plot(gp, "", "W_c", is[7]+1);
					    fprintf(gp, "unset label 1\n");
					  }					  
					  //for wl
					  // major x axis
					  fprintf(gp, "set yrange [0:wlmax] \n");
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    //_Theta_d = _Theta_u;
					    // set x axis label at bottom
					    plot(gp, "", "W_l", is[7]+1);					   					    
					    fprintf(gp, "unset label 2\n");
					  }
					  //for wcf
					  // major x axis
					  /*fprintf(gp, "set yrange [0:wcfmax] \n");
					  for(i[7] = 0; i[7] < is[7]; i[7]++){  // Theta_u
					    _Theta_u = Theta_u[i[7]];
					    //_Theta_d = _Theta_u;
					    // set x axis label at bottom
					    fprintf(gp, "set label 2 '%s' at screen %.4f,%.4f center font \",12\"\n", "b",  (1.0/(is[7]))*(i[7]+1) - 0.5/is[7], 0.02);
					    plot(gp, "", "W_cf", is[7]+1);					   					    
					    fprintf(gp, "unset label 2\n");
					  }*/
					  
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
  double x, y, z, p, q, P, Q, wc, wl, wcf;
  char str[1000], fout[500], fin[500];
  FILE *fpix, *fpip, *fpiP, *fpiQ, *fpiw, *fpox, *fpoy, *fpoz, *fpop, *fpoq, *fpoP, *fpoQ, *fpowc, *fpowl, *fpowcf;  
  
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
		_Theta_d = Theta_d[i[8]];
		//_Theta_d = _Theta_u;
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
			  for(i[14] = 0; i[14] < is[14]; i[14]++){  // Vop
			    _Vop = Vop[i[14]];
			    for(i[15] = 0; i[15] < is[15]; i[15]++){  // Rho
			      _Rho = Rho[i[15]];
			      for(i[16] = 0; i[16] < is[16]; i[16]++){  // x_0
				 _x_0 = x_0[i[16]];
				for(i[17] = 0; i[17] < is[17]; i[17]++){  // y_0
				  _y_0 = y_0[i[17]];
				  for(i[18] = 0; i[18] < is[18]; i[18]++){  //z_0 
				    _z_0 = z_0[i[18]];
				    for(i[20] = 0; i[20] < is[20]; i[20]++){  // Y0	
				      _Y0 = Y0[i[20]];	
				      for(i[21] = 0; i[21] < is[21]; i[21]++){  // e
					_e = e[i[21]];
					for(i[22] = 0; i[22] < is[22]; i[22]++){  // E 
					  _E = E[i[22]];     
					    for(i[0] = 0; i[0] < is[0]; i[0]++){  // n    
					      _n = n[i[0]];
					      // // create files (B, x and X0), (B, y and X0), and (B, P and X0) and payofss Ws
					      prep_file_out(fout, "_xdata.dat");
					      fpox = fopen(fout, "w");
					      prep_file_out(fout, "_ydata.dat");
					      fpoy = fopen(fout, "w");
					      prep_file_out(fout, "_zdata.dat");
					      fpoz = fopen(fout, "w");
					      prep_file_out(fout, "_pdata.dat");
					      fpop = fopen(fout, "w");
					      prep_file_out(fout, "_qdata.dat");
					      fpoq = fopen(fout, "w");
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
					      fprintf(fpop, "b "); fprintf(fpoq, "b ");
					      for(i[19] = 0; i[19] < is[19]; i[19]++)  // X0						
					      {
						fprintf(fpox, "X0=%.1f  ", X0[i[19]]);
						fprintf(fpoy, "X0=%.1f  ", X0[i[19]]);
						fprintf(fpoz, "X0=%.1f  ", X0[i[19]]);
						fprintf(fpop, "X0=%.1f  ", X0[i[19]]);
						fprintf(fpoq, "X0=%.1f  ", X0[i[19]]);
						fprintf(fpoP, "X0=%.1f  ", X0[i[19]]);
						fprintf(fpoQ, "X0=%.1f  ", X0[i[19]]);
						fprintf(fpowc, "X0=%.1f  ", X0[i[19]]);
						fprintf(fpowl, "X0=%.1f  ", X0[i[19]]);
						fprintf(fpowcf, "X0=%.1f  ", X0[i[19]]);
					      }
					      fprintf(fpox, "\n"); fprintf(fpoy, "\n"); fprintf(fpoz, "\n"); 
					      fprintf(fpop, "\n"); fprintf(fpoq, "\n");
					      fprintf(fpoP, "\n"); fprintf(fpoQ, "\n"); 
					      fprintf(fpowc, "\n"); fprintf(fpowl, "\n"); fprintf(fpowcf, "\n");
					      for(i[1] = 0; i[1] < is[1]; i[1]++){  //b
						_b = b[i[1]];
						_B = _b;
						fprintf(fpox, "%.1lf ", _b);
						fprintf(fpoy, "%.1lf ", _b);
						fprintf(fpoz, "%.1lf ", _b);
						fprintf(fpop, "%.1lf ", _b);
						fprintf(fpoq, "%.1lf ", _b);
						fprintf(fpoP, "%.1lf ", _b);
						fprintf(fpoQ, "%.1lf ", _b);
						fprintf(fpowc, "%.1lf ", _b);
						fprintf(fpowl, "%.1lf ", _b);
						fprintf(fpowcf, "%.1lf ", _b);
						for(i[19] = 0; i[19] < is[19]; i[19]++)  // X0						
						{
						    _X0 = X0[i[19]];
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
						    
						    // p & q
						    prep_file_in(fin, "punsum.dat");
						    if((fpip = fopen(fin, "r"))){
						    if(fgets(str, 1000, fpip) == NULL) printf("fgets read error!\n");  // read 1st line header 
						      // read data P
						      if(fscanf(fpip, "%lf %lf", &p, &q) == 2){
							// write data
							fprintf(fpop, "%.4f  ", p);
							fprintf(fpoq, "%.4f  ", q);	
						      }
						    }
						    else{
						      printf("Can't open %s!\n", fin);
						    }
						    fclose(fpip);
						    
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
						    
						    /*prep_file_in(fin, "Qsum.dat");
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
						    fclose(fpiQ); */
						    
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
						fprintf(fpop, "\n"); fprintf(fpoq, "\n"); 
						fprintf(fpoP, "\n"); fprintf(fpoQ, "\n"); 
						fprintf(fpowc, "\n"); fprintf(fpowl, "\n"); fprintf(fpowcf, "\n");
					      }
					      fclose(fpox); fclose(fpoy); fclose(fpoz); fclose(fpop); fclose(fpoq); fclose(fpoP); fclose(fpoQ); fclose(fpowc); fclose(fpowl); fclose(fpowcf);
					   
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
  is[14]  = sizeof(Vop)/sizeof(int);
  is[15]  = sizeof(Rho)/sizeof(double);  
  is[16]  = sizeof(x_0)/sizeof(double);
  is[17]  = sizeof(y_0)/sizeof(double);
  is[18]  = sizeof(z_0)/sizeof(double);
  is[19]  = sizeof(X0)/sizeof(double);
  is[20]  = sizeof(Y0)/sizeof(double);
  is[21]  = sizeof(e)/sizeof(double);
  is[22]  = sizeof(E)/sizeof(double);
  
  prep_data_file();
  
  plot_data_points();
  
  return 0;
}




/** Usage:
    compile : gcc -o xsum xsum.c -lm
    run     : ./xsum
**/


