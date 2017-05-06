#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h> 
#include<stdbool.h>

/* Parameters for plotting */
int      Runs;
int      G;
int      N;       
double   K, B2, Theta, Eta;       
double   C0, C1, C2, C3, C4;        
double   V1, V2, V3;
double   X0;
double x_0, y_0, z_0, X0;            // x_0: half effort of commoner; y_0: half effort of leader; z_0: half effort of chief; X0: half strength of polity
double e, E;                      // e: leader efficiency; E: chief efficiency

/* Parameters for plotting */

#define SKIP 10
#define PREC "7.4" /* printing precision */

void prep_file(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "g%02dn%02dk%.2fb2%.2ftheta%.2feta%.2fc0%.2fc1%.2fc2%.2fc3%.2fc4%.2fv1%.2fv2%.2fv3%.2fX0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s", G, N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, V1, V2, V3, X0, e, E, x_0, y_0, z_0, apndStr);
}

// plot data using lines for all columns as y axis data
void plotef(int wid, char *xdata, char *pdata, char *tdata, char *gpdata, char *sdata, int datacolumn, char *title, char *outputfile, int r)
{    
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  
  fprintf(gp, "set key outside vertical\n");  
  fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,8\" \n");
  fprintf(gp, "set output '%s' \n", outputfile);  
  fprintf(gp, "set lmargin at screen 0.1 \n");
  fprintf(gp, "set rmargin at screen 0.8 \n");
  
  fprintf(gp, "stats '%s' using 2:3 prefix 'A' nooutput\n", xdata);
  fprintf(gp, "stats '%s' using 4 prefix 'B' nooutput\n", xdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'C' nooutput\n", pdata);
  fprintf(gp, "stats '%s' using 4 prefix 'D' nooutput\n", pdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'E' nooutput\n", gpdata);  
  fprintf(gp, "stats '%s' using 2:3 prefix 'F' nooutput\n", sdata);

  fprintf(gp, "set xlabel 'Time (X%d)' \n", SKIP);
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 4,1 title '%s' \n", title);    // set subplots layout
  if(r == -1){ // average
    fprintf(gp, "set label 1 'Average' at screen 0.12,0.98 center font \",13\"\n");  // label to indicate average
  }
  else{
    fprintf(gp, "set label 1 'run: %d' at screen 0.13,0.98 center font \",13\"\n", r);  // label to indicate individual run index
  }
  
  fprintf(gp, "set ylabel 'efforts (com & lead)' \n");
  fprintf(gp, "set y2label 'efforts (chief)' \n");
  
  fprintf(gp, "unset autoscale y\n");
  fprintf(gp, "ymax = A_max_x\n");
  fprintf(gp, "if(A_max_y > A_max_x) {ymax = A_max_y}\n");  
  fprintf(gp, "if(ymax <= 0.0) {ymax = 1.0}\n");
  fprintf(gp, "y2max = B_max\n");    
  fprintf(gp, "if(y2max == 0.0) {y2max = 1.0}\n");  
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set format y2 \"%%.2f\"\n");  
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "set y2tics y2max/3 nomirror \n");  
  fprintf(gp, "set yrange [0:ymax+ymax/10]\n");  
  fprintf(gp, "set y2range [0:y2max+y2max/10]\n");  
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader axes x1y1, '%s' using 1:4 with lines lw 2 title columnheader axes x1y2 \n", datacolumn-1, xdata, xdata);
  fprintf(gp, "unset y2tics\n set ytics mirror\n");
  
  fprintf(gp, "unset y2label \n");
  fprintf(gp, "set ylabel 'tax' \n");
  fprintf(gp, "set ytics 0.5 nomirror\n");  
  fprintf(gp, "set yrange [0.00:1+0.1]\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 lt col title columnheader \n", datacolumn-1, tdata);
  
  fprintf(gp, "set ylabel 'payoff' \n");
  fprintf(gp, "ymax = C_max_x\n");
  fprintf(gp, "if(C_max_y > C_max_x) {ymax = C_max_y}\n");
  fprintf(gp, "if(D_max > ymax) ymax = D_max\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "set autoscale y\n");  
  fprintf(gp, "set yrange [:ymax+ymax/10]\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, pdata);
  fprintf(gp, "unset autoscale y\n");
  
  fprintf(gp, "set ylabel 'P' \n");
  fprintf(gp, "set y2label 'S' \n");  
  fprintf(gp, "ymax = E_max_x\n");
  fprintf(gp, "y2max = F_max_x\n");
  fprintf(gp, "if(E_max_y > E_max_x) {ymax = E_max_y}\n");
  fprintf(gp, "if(F_max_y > F_max_x) {y2max = F_max_y}\n");
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  fprintf(gp, "if(y2max < 0.05) {y2max = 1.0}\n");
  fprintf(gp, "set yrange [0:ymax+0.01]\n");
  fprintf(gp, "set y2range [0:y2max+0.01]\n");  
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set format y2 \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "set y2tics y2max/3 nomirror \n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader axes x1y1,  for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader axes x1y2 \n", 2, gpdata, 2, sdata);
  fprintf(gp, "unset y2tics\n set ytics mirror\n");
  fprintf(gp, "unset y2label \n");
  fprintf(gp, "unset label 1\n");
    
  fprintf(gp, "unset multiplot \n");
  fflush(gp); 
  pclose(gp);
}


int main(int argc, char **argv)
{
  if(argc ^ 22){ printf("Usage: ./fplot g n k B2 theta eta c0 c1 c2 c3 c4 v1 v2 v3 X0 e E x0 y0 z0 runs\n"); return 0;}
  
  G = atoi(argv[1]);
  N = atoi(argv[2]);
  K = atof(argv[3]);
  B2 = atof(argv[4]);
  Theta = atof(argv[5]);
  Eta = atof(argv[6]);
  C0 = atof(argv[7]);
  C1 = atof(argv[8]);
  C2 = atof(argv[9]);
  C3 = atof(argv[10]);
  C4 = atof(argv[11]);
  V1 = atof(argv[12]);
  V2 = atof(argv[13]);
  V3 = atof(argv[14]);
  X0 = atof(argv[15]);
  e  = atof(argv[16]);
  E  = atof(argv[17]);
  x_0 = atof(argv[18]);
  y_0 = atof(argv[19]);
  z_0 = atof(argv[20]);
  Runs = atoi(argv[21]);
	 
  //open file for writing average efforts per rank per generation
  int i;
  char xdata[200], pdata[200], tdata[200], tstr[200], gpdata[200], sdata[200];
  prep_file(xdata, "x.dat"); prep_file(pdata,"p.dat"); prep_file(tdata, "t.dat"); prep_file(sdata, "s.dat");
  prep_file(gpdata, "gp.dat");
  
  char title[200], xpng[200];
  
  // plot average run
  sprintf(title, "g= %d, n= %d, k=%.2f, b2=%.2f, theta= %.2f, eta= %.2f, c=%.1f%.1f%.1f%.1f%.1f, v=%d%d%d, X0=%.1f, e=%.1f, E=%.1f, x0=%.1f, y0=%.1f, z0=%.1f", G, N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, (int)(V1*10), (int)(V2*10), (int)(V3*10), X0, e, E, x_0, y_0, z_0);
  sprintf(xpng, "g%02dn%02dk%0.2fb2%.2fth%.2fet%.2fc%.1f%.1f%.1f%.1f%.1fe%.1fE%.1fX0%.1fx0%.1fy0%.1fz0%.1fv%d%d%d.png", G, N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, e, E, X0, x_0, y_0, z_0, (int)(V1*10), (int)(V2*10), (int)(V3*10)); 
  plotef(0, xdata, pdata, tdata, gpdata, sdata, 3+1, title, xpng, -1);    
  // plot for individual runs
  for(i = 0; i < Runs; i++){
    sprintf(tstr, "x%d.dat", i);
    prep_file(xdata, tstr);
    sprintf(tstr, "p%d.dat", i);
    prep_file(pdata, tstr);
    sprintf(tstr, "t%d.dat", i);
    prep_file(tdata, tstr);
    sprintf(tstr, "gp%d.dat", i);
    prep_file(gpdata, tstr);  
    sprintf(tstr, "s%d.dat", i);
    prep_file(sdata, tstr);
    sprintf(title, "g= %d, n= %d, k= %.2f, b2= %.2f, theta= %.2f, eta= %.2f, c=%.1f%.1f%.1f%.1f%.1f, v=%d%d%d, x0=%.1f, e=%.1f, E=%.1f, x0=%.1f, y0=%.1f, z0=%.1f", N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, (int)(V1*10), (int)(V2*10), (int)(V3*10), X0, e, E, x_0, y_0, z_0);
    sprintf(xpng, "g%02dn%02dk%0.2fb2%.2fth%.2fet%.2fc%.1f%.1f%.1f%.1f%.1fe%.1fE%.1fX0%.1fx0%.1fy0%.1fz0%.1fv%d%d%d_%d.png", G, N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, e, E, X0, x_0, y_0, z_0, (int)(V1*10), (int)(V2*10), (int)(V3*10), i+1); 
    plotef(0, xdata, pdata, tdata, gpdata, sdata, 3+1, title, xpng, i);        
  }  

  return 0;
}



/** Usage:
    compile : gcc -o fplot fplot.c -lm
    run     : ./fplot n k b0 b1 b2 c0 c1 c2 v0 v1 v2 runs
**/

