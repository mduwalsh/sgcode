#include<stdlib.h>
#include<stdio.h>
#include<math.h>
#include<string.h> 
#include<stdbool.h>

/* Parameters for plotting */
int      Runs;
int      n, G;       
double   b, B, Theta_u, Theta_d, Eta_u, Eta_d;       
double   cx, cy, cz;        
double   V1, V2, V3;
double   x_0, y_0, z_0, Y0;            // x_0: half effort of commoner; y_0: half effort of leader; z_0: half effort of chief; X0: half strength of polity
double   e, E;                      // e: leader efficiency; E: chief efficiency

/* Parameters for plotting */

#define SKIP 10
#define PREC "7.4" /* printing precision */


void prep_file(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "g%02dn%02db%0.2fB%.2ftu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fv1%.2fv2%.2fv3%.2fY0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s", G, n, b, B, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, V1, V2, V3, Y0, e, E, x_0, y_0, z_0, apndStr);
}

// plot data using lines for all columns as y axis data
void plotef(int wid, char *xdata, char *pdata, char *tdata, char *gpdata, char *pundata, int datacolumn, char *title, char *outputfile, int r)
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
  fprintf(gp, "stats '%s' using 2:3 prefix 'F' nooutput\n", pundata);

  fprintf(gp, "set xlabel 'Time (X%d)' \n", SKIP);
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 4,1 title '%s' \n", title);    // set subplots layout
  if(r == -1){ // average
    fprintf(gp, "set label 1 'Average' at screen 0.12,0.99 center font \",13\"\n");  // label to indicate average
  }
  else{
    fprintf(gp, "set label 1 'run: %d' at screen 0.13,0.99 center font \",13\"\n", r);  // label to indicate individual run index
  }
  
  fprintf(gp, "unset autoscale y\n");
  fprintf(gp, "ymax = A_max_x\n");
  fprintf(gp, "if(A_max_y > A_max_x) {ymax = A_max_y}\n");  
  fprintf(gp, "if(ymax < B_max) {ymax = B_max}\n");
  fprintf(gp, "if(ymax <= 0.0) {ymax = 1.0}\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "set yrange [0:ymax+ymax/10]\n");    
  fprintf(gp, "set key outside vertical height -1\n");  
  fprintf(gp, "set ylabel 'efforts' \n");  
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, xdata);
  fprintf(gp, "unset y2tics\n set ytics mirror\n");

  fprintf(gp, "set key outside vertical height -1\n");  
  fprintf(gp, "set ylabel 'payoff' \n");
  fprintf(gp, "ymax = C_max_x\n");
  fprintf(gp, "ymin = C_min_x\n");
  fprintf(gp, "if(C_max_y > C_max_x) {ymax = C_max_y}\n");
  fprintf(gp, "if(D_max > ymax) ymax = D_max\n");
  fprintf(gp, "if(C_min_y < C_min_x) {ymin = C_min_y}\n");
  fprintf(gp, "if(D_min < ymin) ymin = D_min\n");
  fprintf(gp, "yr = ymax-ymin\n");
  fprintf(gp, "set ytics yr/3 nomirror\n");
  fprintf(gp, "set autoscale y\n");  
  fprintf(gp, "set yrange [ymin:ymax+ymax/2]\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, pdata);
  fprintf(gp, "unset autoscale y\n");
  
  fprintf(gp, "set key outside vertical height -1\n");  
  fprintf(gp, "set ylabel 'p & q' \n");
  fprintf(gp, "ymax = F_max_x\n");
  fprintf(gp, "if(F_max_y > F_max_x) {ymax = F_max_y}\n");
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  fprintf(gp, "set yrange [0:ymax+0.1]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 lc col title columnheader \n", 3, pundata);  
  
  fprintf(gp, "set ylabel 'P & Q' \n");
  fprintf(gp, "ymax = E_max_x\n");
  fprintf(gp, "if(E_max_y > E_max_x) {ymax = E_max_y}\n");
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  fprintf(gp, "set yrange [0:ymax+0.01]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", 3, gpdata);
  
  fprintf(gp, "unset label 1\n");
  fprintf(gp, "unset multiplot \n"); 
  
  fflush(gp); 
  pclose(gp);
}


int main(int argc, char **argv)
{
  if(argc ^ 22){ printf("Usage: ./fplot G n b B theta_u theta_d eta_u eta_d cx cy cz v1 v2 v3 e E Y0 x0 y0 z0 runs\n"); return 0;}
  
  G = atoi(argv[1]);
  n = atoi(argv[2]);
  b = atof(argv[3]);
  B = atof(argv[4]);
  Theta_u = atof(argv[5]);
  Theta_d = atof(argv[6]);
  Eta_u = atof(argv[7]);
  Eta_d = atof(argv[8]);
  cx = atof(argv[9]);
  cy = atof(argv[10]);
  cz = atof(argv[11]);
  V1 = atof(argv[12]);
  V2 = atof(argv[13]);
  V3 = atof(argv[14]);
  e  = atof(argv[15]);
  E  = atof(argv[16]);
  Y0 = atof(argv[17]);
  x_0 = atof(argv[18]);
  y_0 = atof(argv[19]);
  z_0 = atof(argv[20]);
  Runs = atoi(argv[21]);
  
  //open file for writing average efforts per rank per generation
  int i;
  char xdata[200], pdata[200], tdata[200], tstr[200], gpdata[200], pundata[200];
  prep_file(xdata, "x.dat"); prep_file(pdata,"p.dat"); prep_file(tdata, "t.dat"); prep_file(pundata, "pun.dat");
  prep_file(gpdata, "gp.dat");
  
  char title[200], xpng[200];
  
  // plot average run
  sprintf(title, "b:%.1f, B:%.1f t_u:%.1f, t_d:%.1f, e_u:%.1f, e_d%.1f, c:%.1f,%.1f,%.1f, v:%d%d%d, e:%.1f, E:%.1f, x0:%.2f, y0:%.2f, z0:%.2f, Y0:%.2f", b, B, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, (int)(V1*10), (int)(V2*10), (int)(V3*10), e, E, x_0, y_0, z_0, Y0);
  sprintf(xpng, "g%02dn%02db%0.2fB%.2ftu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fY0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2fv%d%d%d.png", G, n, b, B, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, Y0, e, E, x_0, y_0, z_0, (int)(V1*10), (int)(V2*10), (int)(V3*10)); 
  plotef(0, xdata, pdata, tdata, gpdata, pundata, 3+1, title, xpng, -1);    
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
    sprintf(tstr, "pun%d.dat", i);
    prep_file(pundata, tstr);
    sprintf(title, "b:%.1f, B:%.1f t_u:%.1f, t_d:%.1f, e_u:%.1f, e_d%.1f, c:%.1f,%.1f,%.1f, v:%d%d%d, e:%.1f, E:%.1f, x0:%.2f, y0:%.2f, z0:%.2f, Y0:%.2f", b, B, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, (int)(V1*10), (int)(V2*10), (int)(V3*10), e, E, x_0, y_0, z_0, Y0);
    sprintf(xpng, "g%02dn%02db%0.2fB%.2ftu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fY0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2fv%d%d%d_%d.png", G, n, b, B, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, Y0, e, E, x_0, y_0, z_0, (int)(V1*10), (int)(V2*10), (int)(V3*10), i+1); 
    plotef(0, xdata, pdata, tdata, gpdata, pundata, 3+1, title, xpng, i);        
  }  

  return 0;
}



/** Usage:
    compile : gcc -o fplot fplot.c -lm
    run     : ./fplot g n k B2 theta eta c0 c1 c2 c3 c4 v1 v2 v3 X0 e E x0 y0 z0 runs
**/

