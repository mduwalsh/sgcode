#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <float.h>
#include <time.h>
#include <unistd.h>
 
#define DEBUG 1
#if DEBUG
#define __USE_GNU
#include <fenv.h>       /* enable floating exceptions */
unsigned long Where;    // debugging counter
#endif

#define NB_MULTI_SIM 1

#include "rand.c"

#define LAST 500
#define TestPlot 0            // plot test group statistics or not
#define TestG 2               // test group index

#define UsVsNature_Them 1    // 1: Us vs Nature game; 2: Us vs Them game
#define PIB_MFA 1            // 1: Predicting individual behavior; 2: Mean-field approximation // forecast method

#define AVERAGE_GRAPH_ONLY 0  // if 1, generate only average graphs and supress individual run graphs
#define ALLDATAFILE 1        // if 1, generate all data files for individual runs and summary too 
#define GRAPHS      0        // if 1, saves graphs as png files 

#define INIT_COM_EFFORT rnd(2)              // 0 or 1
#define INIT_LEAD_PUN_EFFORT U01()
#define INIT_LEAD_NORM_EFFORT U01()

// change values to 0 if update strategy for any one of commoner or lead is to be turned off
#define UPDATE_COM 1
#define UPDATE_LEAD_PUN_EFFORT 1
#define UPDATE_LEAD_NORM_EFFORT 1

#define CLUSTER 1            // is this simulation on cluster yes or no (1 or 0)
#define SKIP 1           // time interval between snapshot of states

int Nsim[] = {4, 8, 16};
double Bsim[] = {1, 2, 3, 4, 5, 6, 7, 8, 9, 10};
int Vsim[] = {1, 2, 3};

double LOG_FLT_MAX; 

void *Malloc(size_t size)
{
  void *p = malloc(size);

  if (p) return p;
  printf("malloc failed\n");
  _exit(2);
}

#define malloc(x) Malloc(x) 

void *Calloc(size_t nmemb, size_t size)
{
  void *p = calloc(nmemb, size);

  if (p) return p;
  printf("calloc failed\n");
  _exit(3);
}

#define calloc(x,y) Calloc(x,y) 

#define MAX(x,y)         (((x)>(y))? (x): (y))
#define MIN(x,y)         (((x)<(y))? (x): (y))

// system structure
typedef struct
{
  unsigned int  x;            // 1 or 0 
  double        pi;           // payoff   
  double        uc;           // utility function
} Commoner;

typedef struct
{
  double        y;            // punishment effort
  double        z;            // norm promoting effort
  double        pi;           // payoff 
  double        ul;           // utility function
} Leader;


typedef struct
{
  double        P;            // group production P(x) value
  double        X;            // sum of efforts from commoners in group
  double        theta;        // theta: tax on commoners in a group 
  
  unsigned int  nc;           // number of commoners in a group
  Leader        *lead;        // leader
  Commoner      *com;         // group members / commoners; 
} group;

typedef struct
{
  unsigned int ng;            // no. of groups in a polity  
  double SX;                  // sum of group efforts
  group        *g;            // polity members 
} polity;

polity *Polity;               // "whole system"

// system structure ends
// Global variables
/** configuration parameters **/
unsigned Seed;                    // Seed: Seed of psuedorandom number from config file
unsigned long Seed_i;             // Seed_i: seed of each run
unsigned int Runs;                         // no. of runs
unsigned int G;                            // no. of groups
unsigned int n;                            // no. of commoners
double b;                                  // benefit factor
unsigned int T;                            // time period of simulation
double Theta, Theta_a;                     // taxes / rewards
double cx, cy, cz;                // cost parameter
double Sigma;                     // standard deviation for distribution of mutation for contribution
double Vc1, Vc2, Vc3;                // probability with which participants choose method of changing efforts 
double Vl1, Vl2, Vl3;
double X0, z_0;            // X0: half effort of commoner; z_0: half coordination effort by leaders
double k, delta;        // punishments parameter
double m;                   // probability of migration
double e;                   // efficiency parameter
double s0, s1;
double Lambda;
unsigned int K;                   // candidate strategies
double Upc;

// statistical variables
double *xmean;          // average effort by individuals
double *ymean;            // average effort by leaders
double *zmean;            // average effort by chiefs
double *pi0mean;                  // average payoff of commoners
double *pi1mean;                  // average payoff of leaders
double *Pmean;                    // group production
double *ucmean;                   // commoners utility function
double *ulmean;                   // leader utility function
double xvmean, yvmean, zvmean;    // standard deviation
double xmlv, ymlv, zmlv;          // standard deviation of average values for 10 runs

double xm_l, xv_l, ym_l, yv_l, zm_l, zv_l, p0m_l, p1m_l, ucm_l, ulm_l, Pm_l;

dist *VCdist;                      // 'dist' is defined in rand.c; Probability distribution of choosing method of changing efforts
dist *VLdist;

#define EXPECT(a,b,c) if ((a) != fscanf(f, " " b "%*[^\n]",c)){ fclose(f); printf("Error: %s\n",b); return 1; }

int read_config(char *file_name)
{
  FILE *f;

  if (!(f = fopen(file_name,"r"))) return 1;
  EXPECT(1, "unsigned Seed      = %u;", &Seed);
  EXPECT(1, "int      Runs      = %d;", &Runs);
  EXPECT(1, "int      T         = %d;", &T);
  
  EXPECT(1, "int      n         = %d;", &n);
  EXPECT(1, "int      G         = %d;", &G); 
  EXPECT(1, "int      K         = %d;", &K);
  
  EXPECT(1, "double   b         = %lf;", &b);  
  EXPECT(1, "double   cx        = %lf;", &cx);
  EXPECT(1, "double   cy        = %lf;", &cy);
  EXPECT(1, "double   cz        = %lf;", &cz);     
  
  EXPECT(1, "double   k         = %lf;", &k);
  EXPECT(1, "double   delta     = %lf;", &delta);  
  EXPECT(1, "double   Theta     = %lf;", &Theta_a);
 
  EXPECT(1, "double   Vc1       = %lf;", &Vc1);
  EXPECT(1, "double   Vc2       = %lf;", &Vc2);
  EXPECT(1, "double   Vc3       = %lf;", &Vc3); 
  EXPECT(1, "double   Vl1       = %lf;", &Vl1);
  EXPECT(1, "double   Vl2       = %lf;", &Vl2);
  EXPECT(1, "double   Vl3       = %lf;", &Vl3); 
  EXPECT(1, "double   m         = %lf;", &m);
  
  EXPECT(1, "double   X0        = %lf;", &X0);  
  EXPECT(1, "double   z0        = %lf;", &z_0);
  EXPECT(1, "double   e         = %lf;", &e);
  
  EXPECT(1, "double   s0        = %lf;", &s0);  
  EXPECT(1, "double   s1        = %lf;", &s1); 
  EXPECT(1, "double   Lambda    = %lf;", &Lambda);    
  
  EXPECT(1, "double   Sigma     = %lf;", &Sigma);

  fclose(f);   
  
  if (Runs < 1) exit(1);
  return 0;
}
#undef EXPECT

void prep_file(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "g%02dn%02dK%02db%0.2fk%0.2fdelta%.2fTheta%0.2fcx%.2fcy%.2fcz%.2fX0%.2fe%.2fz0%.2fs0%.2fs1%.2fvc1%.2fvc2%.2fvc3%.2fvl1%.2fvl2%.2fvl3%.2f%s", G, n, K, b, k, delta, Theta_a, cx, cy, cz, X0, e, z_0, s0, s1, Vc1, Vc2, Vc3, Vl1, Vl2, Vl3, apndStr);
}

// generates random number following exponential distribution
double randexp(double lambda)
/*
 * lambda: rate parameter of exponential distribution
 */
{
  return -log(1.0-U01()) / lambda;  
}

// merge sort
void merge (double *a, int n, int m) {
    int i, j, k;
    double *x = malloc(n * sizeof(double));
    for (i = 0, j = m, k = 0; k < n; k++) {
        x[k] = j == n      ? a[i++]
             : i == m      ? a[j++]
             : a[j] < a[i] ? a[j++]
             :               a[i++];
    }
    for (i = n; i--;) {
        a[i] = x[i];
    }
    free(x);
}
 
void merge_sort (double *a, int n) {
    if (n < 2)
        return;
    int m = n>>1;   // divide by 2
    merge_sort(a, m);
    merge_sort(a + m, n - m);
    merge(a, n, m);
}

void merge_key_value (int *a, double *v,  int n, int m) {
    int i, j, k;
    int *x = malloc(n * sizeof(int));
    for (i = 0, j = m, k = 0; k < n; k++) {
        x[k] = j == n      ? a[i++]
             : i == m      ? a[j++]
             : v[a[j]] < v[a[i]] ? a[j++]
             :               a[i++];
    }
    for (i = n; i--;) {
        a[i] = x[i];
    }
    free(x);
}

void merge_sort_key_value (int *a, double *v, int n) {
    if (n < 2)
        return;
    int m = n>>1;   // divide by 2
    merge_sort_key_value(a, v, m);
    merge_sort_key_value(a + m, v, n - m);
    merge_key_value(a, v, n, m);
} 

// allocate memory for system
void setup()
{
  int j;
  group *g;
  // allocate memory for Polity and groups and commoners
  Polity = malloc(sizeof(polity));                                // one polity 
  Polity->g = malloc(G*sizeof(group));                            // allocate memory for groups in each polity    
  for(j = 0; j < G; j++){                                         // through all groups in a polity
    g = Polity->g+j;                                              // get reference to group pointer
    g->lead = malloc(sizeof(Leader));                             // allocate memory for leader in a group
    g->com = malloc(n*sizeof(Commoner));                          // allocate memory for commoners in a group
  }
  VCdist = allocdist(4);                                           // allocate memory for strategy update method probability distribution 
  VLdist = allocdist(4); 
}

void cleanup()                                                       // cleans Polity system memory
{
  int j;
  group *g;  
  for(j = 0; j < G; j++){                                          // through all groups in a polity
    g = Polity->g+j;                                               // get reference to group pointer
    free(g->lead); 
    free(g->com);
  }
  free(Polity->g);  
  free(Polity);
  freedist(VCdist); freedist(VLdist);
}


void allocStatVar()
{
  // allocate memory to store statistics points of snapshot along time period of simulaltion
  xmean = calloc((int)(T/SKIP+1), sizeof(double));         // for effort by commoners
  ymean = calloc((int)(T/SKIP+1), sizeof(double));         // for pun effort by leaders
  zmean = calloc((int)(T/SKIP+1), sizeof(double));         // for norm effort by leaders
  pi0mean = calloc((int)(T/SKIP+1), sizeof(double));                   // for payoff of commoners
  pi1mean = calloc((int)(T/SKIP+1), sizeof(double));                   // for payoff of leaders
  ucmean = calloc((int)(T/SKIP+1), sizeof(double));                     // for utility funciton of commoners 
  ulmean = calloc((int)(T/SKIP+1), sizeof(double));                     // for utility function by leaders  
  Pmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for group production   
}

void clearStatVar()
{
  free(xmean); free(ymean); free(zmean); free(pi0mean); free(pi1mean); 
  free(ucmean); free(ulmean);  
  free(Pmean); 
}

// calculates all the stats
void calcStat(int d, int r) 
{
  double xm = 0, ym = 0, zm = 0, p0m = 0, p1m = 0, ucm = 0, ulm = 0, Pm = 0;
  double xv = 0, yv = 0, zv = 0;	//SG standard deviations
  int j, i, ncm, nl, l;    
  group *g;
  Commoner *com;
  Leader *ld; 

  static FILE **fp = NULL;  // file pointers for individual run  
  if(d < 0){
    //for (l = 4; l--; fclose(fp[l]));
    free(fp); 
    return;
  }     
  
  if(!d){   
    xm_l = 0, ym_l = 0, xv_l = 0, yv_l = 0, zm_l = 0, zv_l = 0, p0m_l = 0, p1m_l = 0, ucm_l = 0, ulm_l = 0, Pm_l = 0;         // reset last summary variables
    char ixdata[200], ipdata[200], iudata[200], igpdata[200], tstr[100];
    fp = malloc(4*sizeof(FILE *));    
    sprintf(tstr, "x%d.dat", r); prep_file(ixdata, tstr);    
    sprintf(tstr, "p%d.dat", r); prep_file(ipdata, tstr);    
    sprintf(tstr, "u%d.dat", r); prep_file(iudata, tstr);    
    sprintf(tstr, "gp%d.dat", r); prep_file(igpdata, tstr); 
    /*fp[0] = fopen(ixdata, "w");
    fp[1] = fopen(ipdata, "w");
    fp[2] = fopen(iudata, "w");  
    fp[3] = fopen(igpdata, "w");    
    
    // write headers       
    fprintf(fp[0], "0 \t x\ty\tz\n");
    fprintf(fp[1], "0 \t com\tlead\n");
    fprintf(fp[2], "0 \t uc\t ul\n");
    fprintf(fp[3], "0 \tP\n");   */     
  }
 
  // calculate stat of traits   
  for(j = 0; j < G; j++){                                                         // through all groups
    g = Polity->g+j;
    // commoners
    for(i = 0; i < n; i++){                                                       // through all commoners
      com = g->com+i;
      xm += com->x;                                                                // sum of commoner's effort
      p0m += com->pi;                                                                // sum of commoner's payoff
      ucm += com->uc;
      xv += (com->x)*(com->x);
    }      
    // leader stats
    ld = g->lead;                                                                // ind == leader
    ym += ld->y;                                                                // sum of leader's efforts
    zm += ld->z;
    
    yv += (ld->y)*(ld->y);		// SG
    zv += (ld->z)*(ld->z);		// SG
    
    p1m += ld->pi;                                                                 // sum of leader's payoff
    ulm += ld->ul;                                                                   // sum of leader's punishment effort 

    Pm += g->P;                                                                     // sum P over all groups        
  } 
  // compute averages of traits
  ncm = G*n;         // total no. of commoners in system
  nl  = G;           // total no. of leaders in system  
  
  xm /= ncm; 
  p0m /= ncm;
  ucm  /= ncm;
  
  xv = xv/ncm - xm*xm;
  if(xv>=0){ xv = sqrt(xv);} 
  
  ym /= nl;
  zm /= nl;
  
  yv = yv/nl - ym*ym;		//SG
  if (yv>=0){
    yv=sqrt(yv);
  }
  zv = zv/nl - zm*zm;		//SG
  if (zv>=0){
    zv=sqrt(zv);
  }
  
  p1m /= nl;
  ulm /= nl;  
  Pm /= nl;      
  
  // store in averages in variable as sum to average it after multiple runs
  xmean[d] += xm;
  ymean[d] += ym;
  zmean[d] += zm;
  pi0mean[d] += p0m;
  pi1mean[d] += p1m;  
  ucmean[d] += ucm;
  ulmean[d] += ulm;  
  Pmean[d] += Pm;  
       
  // write data for individual runs
  /*fprintf(fp[0], "%d  %.4lf  %.4lf  %.4lf\n", d, xm, ym, zm);
  fprintf(fp[1], "%d  %.4lf  %.4lf\n", d, p0m, p1m);
  fprintf(fp[2], "%d  %.4lf  %.4lf\n", d, ucm, ulm);  
  fprintf(fp[3], "%d  %.4lf\n", d, Pm); */  

  // print final values for each run
  if(d > (T/SKIP - LAST) ){
    xm_l += xm; 
    xv_l += xv;
    ym_l += ym; 
    yv_l += yv;
    zm_l += zm;
    zv_l += zv;
    p0m_l += p0m;
    p1m_l += p1m;
    ucm_l += ucm;
    ulm_l += ulm;
    Pm_l += Pm;
  }
  if(d == T/SKIP){
    xm_l /= LAST; 
    xv_l /= LAST;
    ym_l /= LAST; 
    yv_l /= LAST;
    zm_l /= LAST;
    zv_l /= LAST;
    p0m_l /= LAST;
    p1m_l /= LAST;
    ucm_l /= LAST;
    ulm_l /= LAST;
    Pm_l /= LAST;
    // mean of std deviations
    xvmean += xv_l;
    yvmean += yv_l;
    zvmean += zv_l;
    // std deviation of last mean values
    xmlv += xm_l*xm_l;
    ymlv += ym_l*ym_l;
    zmlv += zm_l*zm_l;
#if !CLUSTER
    //printf("run#%d \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%lu\n", r, xm, ym, zm, p0m, p1m, ucm, ulm, Pm, Seed_i); 
    printf("run#%d \t%.3lf \t%.3lf (%.3lf)\t%.3lf (%.3lf)\t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%lu\n", r, xm_l, ym_l, yv_l, zm_l, zv_l, p0m_l, p1m_l, ucm_l, ulm_l, Pm_l, Seed_i); 
#endif
  }  

}

void writeDataToFile()
{
  int j;
  char xdata[200], pdata[200], udata[200], gpdata[200];
  FILE **fp = malloc(8*sizeof(FILE *));    
  prep_file(xdata, "x.dat");    
  prep_file(pdata, "p.dat");    
  prep_file(udata, "u.dat");   
  prep_file(gpdata, "gp.dat");  
//   fp[0]= fopen(xdata, "w");
//   fp[1]= fopen(pdata, "w");
//   fp[2]= fopen(udata, "w");
//   fp[3]= fopen(gpdata, "w");  

  int stu = T/SKIP + 1 - LAST;
  double xsum = 0, ysum = 0, zsum = 0, p0sum = 0, p1sum = 0, ucsum = 0, ulsum = 0, Psum = 0;  

  // write headers
//   fprintf(fp[0], "0 \t x\ty\tz\n");
//   fprintf(fp[1], "0 \t com\tlead\n");
//   fprintf(fp[2], "0 \t uc \t ul\n");
//   fprintf(fp[3], "0 \t P\n");
  
  for(j = 0; j < (int)(T/SKIP)+1; j++){
    // average accumulated mean values for multiple runs
    xmean[j] /= (double)Runs;
    ymean[j] /= (double)Runs;
    zmean[j] /= (double)Runs;
    pi0mean[j] /= (double)Runs;
    pi1mean[j] /= (double)Runs;   
    ucmean[j] /= (double)Runs;
    ulmean[j] /= (double)Runs;    
    Pmean[j] /= (double)Runs;   
    // write data to file
//     fprintf(fp[0], "%d  %.4lf  %.4lf  %.4lf\n", j, xmean[j], ymean[j], zmean[j]);
//     fprintf(fp[1], "%d  %.4lf  %.4lf\n", j, pi0mean[j], pi1mean[j]);
//     fprintf(fp[2], "%d  %.4lf  %.4lf\n", j, ucmean[j], ulmean[j]); 
//     fprintf(fp[3], "%d  %.4lf\n", j, Pmean[j]);    
    
// LAST points summary values
    if(j < stu) continue;
    xsum += xmean[j];
    ysum += ymean[j];
    zsum += zmean[j];
    p0sum += pi0mean[j];
    p1sum += pi1mean[j];
    
    ucsum += ucmean[j];
    ulsum += ulmean[j]; 
    Psum += Pmean[j];

  }
  

  double st = LAST;
  xsum /= st; ysum /= st; zsum /= st; p0sum /= st; p1sum /= st;  
  ucsum /= st; ulsum /= st; Psum /= st;                                     // computing average over summary period time
  
  // mean of std deviations  
  xvmean /= (double)Runs;
  yvmean /= (double)Runs;
  zvmean /= (double)Runs;
  
  // std deviations of means over individual runs
  xmlv = xmlv/Runs - xsum*xsum;
  ymlv = ymlv/Runs - ysum*ysum;
  zmlv = zmlv/Runs - zsum*zsum;  
  if(xmlv > 0) xmlv = sqrt(xmlv);
  if(ymlv > 0) ymlv = sqrt(ymlv);
  if(zmlv > 0) zmlv = sqrt(zmlv);
   
  
 #if ALLDATAFILE
  // write data to file
  char  xsumdata[200], psumdata[200], usumdata[200], Psumdata[200], stddata[200];
  prep_file(stddata, "std.dat");
  FILE *fpstd = fopen(stddata, "w");
  prep_file(xsumdata, "xsum.dat");    
  prep_file(psumdata, "psum.dat");    
  prep_file(usumdata, "usum.dat");    
  prep_file(Psumdata, "Psum.dat");
  fp[4]= fopen(xsumdata, "w");
  fp[5]= fopen(psumdata, "w");
  fp[6]= fopen(usumdata, "w");
  fp[7]= fopen(Psumdata, "w");
  
  // write headers  
  fprintf(fp[4], "x\ty\tz\n");
  fprintf(fp[5], "com\tlead\n");
  fprintf(fp[6], "uc\t ul\n");
  fprintf(fp[7], "P\n");
     
  fprintf(fp[4], "%.4lf  %.4lf  %.4lf\n", xsum, ysum, zsum);
  fprintf(fp[5], "%.4lf  %.4lf\n", p0sum, p1sum);
  fprintf(fp[6], "%.4lf  %.4lf\n", ucsum, ulsum);    
  fprintf(fp[7], "%.4lf \n", Psum);
  
  // write std deviations
  fprintf(fpstd, "x_v  y_v  z_v\n");
  fprintf(fpstd, "%.4lf  %.4lf  %.4lf\n", xvmean, yvmean, zvmean);
  fprintf(fpstd, "%.4lf  %.4lf  %.4lf\n", xmlv, ymlv, zmlv);
  fclose(fpstd);
#endif  
  
// free file pointers
//   for(j = 0; j < 4; j++){
//     fclose(fp[j]);
//   } 
#if ALLDATAFILE
  for(j = 4; j < 8; j++){
    fclose(fp[j]);
  }
#endif 
  free(fp);
  
#if !CLUSTER
  // print averaged final values  
  //printf("\nAvg: \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \n", xmean[T/SKIP], ymean[T/SKIP], zmean[T/SKIP],  pi0mean[T/SKIP], pi1mean[T/SKIP], ucmean[T/SKIP], ulmean[T/SKIP], Pmean[T/SKIP]);  
  printf("\nAvg: \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \t%.3lf \n", xsum, ysum, zsum,  p0sum, p1sum, ucsum, ulsum, Psum);
#endif
}

void plotall(int m)
/*
 * m: 0 or 1; to save graphs as image file or not
 */
{ 
  // write data to file
  char xdata[200], pdata[200], udata[200], gpdata[200];
  char title[200], xpng[200], str[100];
  int datacolumn = 2+1;
  sprintf(str, "x.dat"); prep_file(xdata, str);
  sprintf(str, "p.dat"); prep_file(pdata, str);
  sprintf(str, "u.dat"); prep_file(udata, str); 
  sprintf(str, "gp.dat"); prep_file(gpdata, str);  
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "n:%d, b:%.1f, K:%d, k:%.1f, delta:%.1f, Theta:%.1f, cx:%.1f, cy:%.1f, cz:%.1f, X0:%.2f, vc:%d%d%d, vl:%d%d%d", n, b, K, k, delta, Theta, cx, cy, cz, X0, (int)(Vc1*10), (int)(Vc2*10), (int)(Vc3*10), (int)(Vl1*10), (int)(Vl2*10), (int)(Vl3*10));
  fprintf(gp, "set key outside vertical  spacing 1 width 1 height -2\n");        
  if(m){ // save graphs as file
    sprintf(xpng, "g%02dn%02dK%02db%0.2fk%.2fdelta%.2ftheta%.2fcx%.2fcy%.2fcz%.2fX0%.2fvc%d%d%dvl:%d%d%d.png", G, n, K, b, k, delta, Theta, cx, cy, cz, X0, (int)(Vc1*10), (int)(Vc2*10), (int)(Vc3*10), (int)(Vl1*10), (int)(Vl2*10), (int)(Vl3*10)); 
    fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,8\" \n");
    fprintf(gp, "set output '%s' \n", xpng);
  }
  else{
    fprintf(gp, "set term x11 dashed %d \n", 0);
  }
  fprintf(gp, "set lmargin at screen 0.1 \n");
  fprintf(gp, "set rmargin at screen 0.8 \n");
  
  fprintf(gp, "stats '%s' using 2:3 prefix 'A' nooutput\n", xdata);
  fprintf(gp, "stats '%s' using 4 prefix 'B' nooutput\n", xdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'C' nooutput\n", pdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'D' nooutput\n", udata);
  fprintf(gp, "stats '%s' using 2 prefix 'E' nooutput\n", gpdata);  
  

  fprintf(gp, "set xlabel 'Time(X%lu)' \n", (unsigned long)SKIP);
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 4,1 title '%s' \n", title);    // set subplots layout
  fprintf(gp, "set label 1 'Average' at screen 0.12,0.99 center font \",13\"\n");  // label to indicate average
  
  
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
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn+1, xdata);
  fprintf(gp, "set ytics mirror\n");

  fprintf(gp, "set key outside vertical height -1\n");  
  fprintf(gp, "set ylabel 'payoff' \n");
  fprintf(gp, "ymax = C_max_x\n");
  fprintf(gp, "ymin = C_min_x\n");
  fprintf(gp, "if(C_max_y > C_max_x) {ymax = C_max_y}\n");  
  fprintf(gp, "if(C_min_y < C_min_x) {ymin = C_min_y}\n");  
  fprintf(gp, "yr = ymax-ymin\n");
  fprintf(gp, "set ytics yr/3 nomirror\n");
  fprintf(gp, "set autoscale y\n");  
  fprintf(gp, "set yrange [ymin:ymax+ymax/2]\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, pdata);
  fprintf(gp, "unset autoscale y\n");
  
  fprintf(gp, "set key outside vertical height -1\n");  
  fprintf(gp, "set ylabel 'uc & ul' \n");
  fprintf(gp, "ymax = D_max_x\n");
  fprintf(gp, "if(D_max_y > D_max_x) {ymax = D_max_y}\n");
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  fprintf(gp, "ymin = D_min_x\n");
  fprintf(gp, "if(D_min_y < D_min_x) {ymin = D_min_y}\n");
  fprintf(gp, "set yrange [ymin:ymax+0.1]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "yr = ymax-ymin\n");
  fprintf(gp, "set ytics yr/3 nomirror\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, udata);  
  
  fprintf(gp, "set ylabel 'P' \n");
  fprintf(gp, "ymax = E_max\n");  
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  fprintf(gp, "set yrange [0:ymax+0.01]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn-1, gpdata);
  fprintf(gp, "unset label 1\n");
  fprintf(gp, "unset multiplot \n");
  
  fflush(gp); 
  pclose(gp);  
  
#if !CLUSTER
#if !ALLDATAFILE
  remove(xdata);
  remove(pdata);
  remove(udata);
  remove(gpdata);  
#endif
#endif
}

void plotallIndividualRun(int r, int m)
/*
 * r: run
 * m: 0 or 1; to save graphs as image file or not
 */
{
  // write data to file
  char xdata[200], pdata[200], udata[200], gpdata[200];
  char title[200], xpng[200], str[100];
  int datacolumn = 2+1;
  sprintf(str, "x%d.dat", r); prep_file(xdata, str);
  sprintf(str, "p%d.dat", r); prep_file(pdata, str);
  sprintf(str, "u%d.dat", r); prep_file(udata, str); 
  sprintf(str, "gp%d.dat", r); prep_file(gpdata, str);  
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "n:%d, b:%.1f, K:%d, k:%.1f, delta:%.1f, Theta:%.1f, cx:%.1f, cy:%.1f, cz:%.1f, X0:%.2f, vc:%d%d%d, vl:%d%d%d", n, b, K, k, delta, Theta, cx, cy, cz, X0, (int)(Vc1*10), (int)(Vc2*10), (int)(Vc3*10), (int)(Vl1*10), (int)(Vl2*10), (int)(Vl3*10));
  fprintf(gp, "set key outside vertical  spacing 1 width 1 height -2\n");        
  if(m){ // save graphs as file
    sprintf(xpng, "g%02dn%02dK%02db%0.2fk%.2fdelta%.2ftheta%.2fcx%.2fcy%.2fcz%.2fX0%.2fvc%d%d%dvl:%d%d%d_%d.png", G, n, K, b, k, delta, Theta, cx, cy, cz, X0, (int)(Vc1*10), (int)(Vc2*10), (int)(Vc3*10), (int)(Vl1*10), (int)(Vl2*10), (int)(Vl3*10), r); 
    fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,8\" \n");
    fprintf(gp, "set output '%s' \n", xpng);
  }
  else{
    fprintf(gp, "set term x11 dashed %d \n", 0);
  }
  fprintf(gp, "set lmargin at screen 0.1 \n");
  fprintf(gp, "set rmargin at screen 0.8 \n");
  
  fprintf(gp, "stats '%s' using 2:3 prefix 'A' nooutput\n", xdata);
  fprintf(gp, "stats '%s' using 4 prefix 'B' nooutput\n", xdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'C' nooutput\n", pdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'D' nooutput\n", udata);
  fprintf(gp, "stats '%s' using 2 prefix 'E' nooutput\n", gpdata);  
  

  fprintf(gp, "set xlabel 'Time(X%lu)' \n", (unsigned long)SKIP);
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 4,1 title '%s' \n", title);    // set subplots layout
  fprintf(gp, "set label 1 'run: %d' at screen 0.13,0.98 center font \",13\"\n", r);  // label to indicate individual run index
  
  
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
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn+1, xdata);
  fprintf(gp, "unset y2tics\n set ytics mirror\n");

  fprintf(gp, "set key outside vertical height -1\n");  
  fprintf(gp, "set ylabel 'payoff' \n");
  fprintf(gp, "ymax = C_max_x\n");
  fprintf(gp, "ymin = C_min_x\n");
  fprintf(gp, "if(C_max_y > C_max_x) {ymax = C_max_y}\n");  
  fprintf(gp, "if(C_min_y < C_min_x) {ymin = C_min_y}\n");  
  fprintf(gp, "yr = ymax-ymin\n");
  fprintf(gp, "set ytics yr/3 nomirror\n");
  fprintf(gp, "set autoscale y\n");  
  fprintf(gp, "set yrange [ymin:ymax+ymax/2]\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, pdata);
  fprintf(gp, "unset autoscale y\n");
  
  fprintf(gp, "set key outside vertical height -1\n");  
  fprintf(gp, "set ylabel 'uc & ul' \n");
  fprintf(gp, "ymax = D_max_x\n");
  fprintf(gp, "if(D_max_y > D_max_x) {ymax = D_max_y}\n");
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  fprintf(gp, "set yrange [:ymax+0.1]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "yr = ymax-ymin\n");
  fprintf(gp, "set ytics yr/3 nomirror\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn, udata);  
  
  fprintf(gp, "set ylabel 'P' \n");
  fprintf(gp, "ymax = E_max\n");  
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  fprintf(gp, "set yrange [0:ymax+0.01]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", datacolumn-1, gpdata);
  fprintf(gp, "unset label 1\n");
  fprintf(gp, "unset multiplot \n");
  
  fflush(gp); 
  pclose(gp);  
}


void testGPlot(char *data, int r)
{   
  FILE *gp = popen("gnuplot -persistent", "w");
  //fprintf(gp, "set key outside vertical spacing 1 width 1 height 4\n");       
  fprintf(gp, "set key outside vertical\n");   
  fprintf(gp, "set term x11 dashed %d \n", 0);
  fprintf(gp, "set lmargin at screen 0.1 \n");
  fprintf(gp, "set rmargin at screen 0.86 \n");
  //set tmargin at screen 0.97
  //set bmargin at screen 0.6  
  fprintf(gp, "stats '%s' using 2 prefix 'A' nooutput\n", data);
  fprintf(gp, "stats '%s' using 3 prefix 'B' nooutput\n", data);
  fprintf(gp, "stats '%s' using 4 prefix 'C' nooutput\n", data);
  fprintf(gp, "stats '%s' using 5 prefix 'D' nooutput\n", data); 
  fprintf(gp, "stats '%s' using 6 prefix 'E' nooutput\n", data);
  fprintf(gp, "stats '%s' using 7 prefix 'F' nooutput\n", data);
  
  fprintf(gp, "set xlabel 'Time(X%lu)' \n", (unsigned long)SKIP);
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 3,1 title '%s, G:%d' \n", "Group dynamics", TestG);    // set subplots layout
  fprintf(gp, "set label 1 'run: %d' at screen 0.13,0.98 center font \",13\"\n", r);  // label to indicate individual run index
  
  fprintf(gp, "unset autoscale y\n");
  fprintf(gp, "ymax = A_max\n");
  fprintf(gp, "if(B_max > A_max) {ymax = B_max}\n");  
  fprintf(gp, "if(ymax < C_max) {ymax = C_max}\n");
  fprintf(gp, "if(ymax <= 0.0) {ymax = 1.0}\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 mirror\n");
  fprintf(gp, "set yrange [0:ymax+ymax/10]\n");    
  fprintf(gp, "set key outside vertical height -1\n");  
  fprintf(gp, "set ylabel 'efforts' \n");  
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", 4, data);
  
  fprintf(gp, "set ylabel 'P' \n");
  fprintf(gp, "ymax = D_max\n");  
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  fprintf(gp, "set yrange [0:ymax+0.01]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "plot '%s' using 1:%d with lines lw 2 title columnheader \n",data, 5);
  
  fprintf(gp, "set ylabel 'payoffs' \n");
  fprintf(gp, "ymax = E_max\n");  
  fprintf(gp, "if(ymax < F_max) {ymax = F_max}\n");
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  fprintf(gp, "set yrange [0:ymax+0.01]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  fprintf(gp, "plot for [col=6:7] '%s' using 1:col with lines lw 2 title columnheader \n",data);
  
  fprintf(gp, "unset label 1\n");
  fprintf(gp, "unset multiplot\n");
  
  fflush(gp);
  pclose(gp);  
}

double P(double z, double X, double rSX)
/*
 * X : group effort
 * SX: sum of group efforts in whole group
 */
{    
//   double r = 1.0 + (e*z)/(n*z_0 + z);
#if UsVsNature_Them == 1  // us vs nature 
//   return r*X/( r*X + X0 );
  return X/( X + X0 );
#endif
#if UsVsNature_Them == 2  // us vs them
  if(rSX <= 0.0) return 1.0;
//   return (r*X*G) / rSX;
  return (X*G) / rSX;
#endif
}

// returns frequency of contributors in group
double q( double X)
{
  return X/n;
}

// return strength of norm internalization
double s(double X)
{
  return s0 + s1*q(X);  
}

// utility function
double u(unsigned int x1, unsigned int x0, double X, double y, double z, double rSX)
{  
  double _s = s(X);
//   double r = 1.0 + (e*z)/(n*z_0 + z);
//   double pi_c = (1.0-Theta)*b*P(z, X-x0+x1,rSX+r*(x1-x0)) - cx*x1 - k*y*(1.0-x1);            // expected material payoff for commoners
  double pi_c = (1.0-Theta)*b*P(z, X-x0+x1,rSX+(x1-x0)) - cx*x1 - k*y*(1.0-x1);            // expected material payoff for commoners
  return (1.0-_s)*pi_c + _s*z*x1;
}

// returns x' using myopic optimization for commoners
double f(int x, double X, double y, double z, double rSX)
{
  double u0, u1;
  u0 = u(0, x, X, y, z, rSX);                               // utility function when x' = 0
  u1 = u(1, x, X, y, z, rSX);                               // utility function when x' = 1

  //return ( u0 < u1 )? 1: 0;  				//SG
 
  double lpd = Lambda*(u0-u1);
  if(lpd < -LOG_FLT_MAX){
    return 1;
  }
  if(lpd > LOG_FLT_MAX){
    return 0;
  }
  else{
    return (U01() < 1.0/( 1.0+exp(lpd) ) )? 1 : 0;  
  }
}

#if PIB_MFA == 1            // Predicting individual behavior
// returns new forecasted sum of efforts from commoners given y and z
double F(double X, double y, double z, unsigned int *x0, unsigned int *x1, double rSX, unsigned int *ux)
/*
 * X: sum of efforts from commoners in group
 * y: punishment effort
 * z: norm effort
 * x0: efforts array of commoners
 * x1: new forecasted efforts array of commoners
 * ux: array indicating commoners updating or not as assumed by leader
 */
{
  int i;
  double s;
  //printf("\n");
  for(s = 0, i = 0; i < n; i++){
   // commented out because now leaders assume only some commoners update efforts with some probability
    //x1[i] = f(x0[i], X, y, z, rSX);                // new forecasted efforts of commoners
    
    // new implementation of leader assumption of commoners updating efforts with some probability
    if(ux[i]){//if(U01() < Upc){
       x1[i] = f(x0[i], X, y, z, rSX);                // new forecasted efforts of commoners
    }
    else{
      x1[i] = x0[i];
    }
    // end of new implementation of leader assumption of commoners updating efforts with some probability
    s += x1[i];       
  }
  return s;
}
#endif
#if PIB_MFA == 2         // Mean-field approximation
// returns new forecasted sum of efforts from commoners given y and z
double F(double X, double y, double z)
{
  double C;
  double _s = s(X);
  if(_s == 1.0){
    if(z > 0.0){ 
      return n;         // all commoners will contribute
    }
    else{
      C = cx - k*y;
    }   
  }
  else{
    C = cx - k*y - (_s*z)/(1.0-_s);
  }
//   double r = 1.0 + (e*z)/(n*z_0 + z);
//   C = C/r;
  if(C <= 0.0) return n;                      // to avoid division by zero and negative C
  double B = (1.0 - Theta)*b;
#if UsVsNature_Them == 1  
  double R = B/(C*X0);
  double X1;
  if(R > 1.0){
    X1 = X0*(sqrt(R) - 1.0);
    if(X1 > n) return n;                    // X' cannot be greater than n
    return X1;
  }
  else
    return 0;
#endif
#if UsVsNature_Them == 2
  double X1;
  X1 = B/C;
  if(X1 > n) return n;
  return X1;
#endif
}

#endif

// returns utility function leader trying to maximize
double u_l(double y, double z, double X1, double X2, double rSX2)
{
  return -cy*n*y - cz*z - delta*(n-X1)*y + Theta*n*b*P(z, X2, rSX2); 
}

// calculates group effort X
void calcX(int j)
{
  int i;  
  double sx;
  group *g = Polity->g+j;
  for(sx = 0, i = 0; i < n; i++){                       // through all commoners
    sx += g->com[i].x;                                   // sum production efforts from commoners
  }
  g->X = sx;                                                // sum of efforts from commoners  
}

// calculates group production P
void calcP(int j, double rSX)
{
  group *g = Polity->g+j;
  g->P = P(g->lead->z, g->X, rSX);
}

// calculates payoff of a group j
void calcPayoff(int j)
{
  int i;   
  group *g;
  Leader *ld;
  Commoner *com;
  double gpc, cp;    
  
  g = Polity->g+j;
  ld = g->lead;
  gpc = (1.0 - g->theta)*b*g->P;                       // group production share to each commoner 
  // update payoff of commoners 
  for(cp = 0, i = 0; i < n; i++){
    com = g->com+i;
    com->pi = gpc - cx*com->x;                               // payoff due to group production - cost of x 
    if(com->x == 0){
      if(U01() < ld->y){
	com->pi -= k;                                        // payoff after punishment from leader with probability y
	cp += delta;
      }
    }    
  }
  // update payoff of leader
  ld->pi = g->theta*n*b*g->P -cy*n*ld->y -cz*ld->z -cp; // payoff due to group production - cost of y - cost of punishing     
}

void calcUtilityFunction(int j, double rSX)
{
  int i;
  Commoner *com;
  Leader *ld;  
  group *g;
  g = Polity->g+j;
  ld = g->lead;
  for(i = 0; i < n; i++){
    com = g->com+i;
    com->uc = u(com->x, com->x, g->X, ld->y, ld->z, rSX);
  }
  ld->ul = u_l(ld->y, ld->z, g->X, g->X, rSX);
}

// returns if it copys from other individual or group 
int sc(double p0, double p1)
/*
 * p0: payoff of updating member
 * p1: payoff of target member to copy
 */
{
  // return (p0 < p1)? 1 : 0;
  
  double lpd = Lambda*(p0-p1);
  if(lpd < -LOG_FLT_MAX){
    return 1;
  }
  if(lpd > LOG_FLT_MAX){
    return 0;
  }
  else{
    return (U01() < 1.0/( 1.0+exp(lpd) ) )? 1 : 0;  
  }
}

// updates strategy of commoner
void updateStrategyCommoner(int v, int j, int i, double y, double z, double rSX)
/*
 * v: type of strategy update
 * j: group index
 * i: commoner index
 * y: leader's punishing effort
 * z: leader's norm promoting effort
 * SX: sum of group efforts across all groups
 */
{
  group *g;
  Commoner *com;
  g = Polity->g+j;
  com = g->com+i;
  // random mutation
  if(v == 1){
    com->x = 1-com->x;  // 0 or 1    
  }
  // selective copying
  else if(v == 2){
    int a;
    if(U01() < (1.0-m)){                                       // copy from same group      
      do{ a = rnd(n); }while( a == i);                        // select another individual in group
      if( sc(com->pi, g->com[a].pi) ){ 
	com->x = g->com[a].x;
      }   
    }
    else{                                                     // choose commoner from another group c in same polity
      int c = j;
      if(G > 1){                                              // if more than one group in polity, choose different group
	do{ c = rnd(G); }while(c == j);	                  // select group c (another group to copy from)
      }      
      a = rnd(n);                                             // select another individual in group b
      if( sc(com->pi, Polity->g[c].com[a].pi) ){ 
	com->x = Polity->g[c].com[a].x;
      }
    }    
  }
  // myopic optimization
  else if(v == 3){    
    com->x = f(com->x, g->X, y, z, rSX);  
  }
}

int updateY0rZ()
{
  int rp;
  double yzp;
  yzp = U01();
  if(yzp < 0.25){
    rp = 0;
  }
  else if(yzp < 0.5){
    rp = 1;
  }
  else if(yzp < 0.75){
    rp = 2;
  }
  else{
    rp = 3;
  }
  return rp;
}

// updates strategy of leader // random mutation
void updateStrategyLeader_V1(int j)
{
  group *g;
  Leader *ld;
  g = Polity->g+j;
  ld = g->lead;    
  int rp = updateY0rZ();
  if(UPDATE_LEAD_PUN_EFFORT == 1 && (rp == 0 || rp == 2)){
    ld->y = normal(ld->y, Sigma);
    ld->y = MAX(ld->y, 0.0);
    ld->y = MIN(ld->y, 1.00);
  }
  if(UPDATE_LEAD_NORM_EFFORT == 1 && (rp == 1 || rp == 2)){
    ld->z = normal(ld->z, Sigma);
    ld->z = MAX(ld->z, 0.0);
    //ld->z = MIN(ld->z, 1.00);
  }
}

// updates strategy of leader // selective copying
void updateStrategyLeader_V2(int j)
{
  group *g;
  Leader *ld;
  g = Polity->g+j;
  ld = g->lead;
  int rp = updateY0rZ();
  int a = j;
  if(G > 1){
    do{ a = rnd(G); }while( a == j);                                   // select another group in polity
    if( sc(ld->pi, Polity->g[a].lead[0].pi) ){                              // if leader in selected group has higher payoff copy strategy
      if( UPDATE_LEAD_PUN_EFFORT == 1 && (rp == 0 || rp == 2)){
	ld->y = (Polity->g+a)->lead->y;                          // copy punishment effort
      }
      if(UPDATE_LEAD_NORM_EFFORT == 1 && (rp == 1 || rp == 2)){
	ld->z = (Polity->g+a)->lead->z;                          // copy norm effort
      }
    }
  }
}


// updates strategy of leader // myopic optimization
// Note: this method changes value of x0 array 
void updateStrategyLeader_V3(int j, unsigned int *x0, double X, double rSX)
/*
 * j: index of group of leader
 * x0: array of unchanged this round efforts from commoners in group j
 * X: group effort of this round
 * SX: sum of group efforts across all groups
 */
{
  group *g;
  Leader *ld;
  g = Polity->g+j;
  ld = g->lead;  
  
  int i, rp;
  double X1, rSX1, s;
  // choose to update y or z or both or neither
  rp = updateY0rZ();
  if(rp == 3) return;     // no update
    
#if PIB_MFA == 1
  unsigned int *x1 = malloc(n*sizeof(unsigned int));             // new forecasted efforts by commoners which sum up to X1 group effort
  unsigned int *ux = malloc(n*sizeof(unsigned int));
  for(i = 0; i < n; i++){
    if(U01() < Upc)
      ux[i] = 1;
    else 
      ux[i] = 0;
  }
#endif
  double *X2 = malloc((K+1)*sizeof(double));                    // forecasted Xs for (K+1) candidate y and z 
  double *y = calloc((K+1),sizeof(double));
  double *z = calloc((K+1),sizeof(double));
  double *ul = malloc((K+1)*sizeof(double));                    // utility function array    
  dist *ul_dist = allocdist(K+1);                               // utility function distribution array  
  double ul_max;
  
  // candidate strategies
  y[0] = ld->y;
  z[0] = ld->z;    
  for(i = 1; i < K+1; i++){    
    if(UPDATE_LEAD_PUN_EFFORT == 1 && (rp == 0 || rp == 2)){
      y[i] = normal(ld->y, Sigma);
      y[i] = MAX(y[i], 0.0);
      y[i] = MIN(y[i], 1.00);
    }
    else{
      y[i] = ld->y;
    }    
    if(UPDATE_LEAD_NORM_EFFORT == 1 && (rp == 1 || rp == 2)){
      z[i] = normal(ld->z, Sigma);
      z[i] = MAX(z[i], 0.0);
      //z[i] = MIN(z[i], 1.0);
    }
    else{
      z[i] = ld->z;
    }     
  }
#if PIB_MFA == 1      // Predicting individual behavior
   X1 = F(X, ld->y, ld->z, x0, x1, rSX, ux);                      // forecasted X for present y and z
#else                 // Mean-field approximation
  X1 = F(X, ld->y, ld->z);
#endif
//   double r = 1.0 + (e*z[0])/(n*z_0 + z[0]);
//   rSX1 = rSX + r*(X1 - X);                                        // new sum of group efforts across all groups
  rSX1 = rSX + (X1 - X);                                        // new sum of group efforts across all groups
  ul[0] = u_l(y[0], z[0], X1, X1, rSX1);                    // utility function
  ul_max = ul[0];
  for(i = 1; i < K+1; i++){      
#if PIB_MFA == 1
    X2[i] = F(X1, y[i], z[i], x1, x0, rSX1, ux);                 // next round forecasted group effort by commoners      
#else
    X2[i] = F(X1, y[i], z[i]);
#endif
//     r = 1.0 + (e*z[i])/(n*z_0 + z[i]);
//     ul[i] = u_l(y[i], z[i], X1, X2[i], rSX1+r*(-X1+X2[i]));                    // utility function
    ul[i] = u_l(y[i], z[i], X1, X2[i], rSX1+(-X1+X2[i]));                    // utility function
    if(ul[i] > ul_max) ul_max = ul[i]; 
  }    
 for(s = 0, i = 0; i < K+1; i++){
    ul_dist->p[i] = exp( (ul[i]-ul_max)*Lambda );
    s += ul_dist->p[i];
  }
  if(s > 0.0){
    initdist(ul_dist, s);                                      // initialize distribution
    // update y and z
    i = drand(ul_dist);
    ld->y = y[i];
    ld->z = z[i];      
  }
  
#if PIB_MFA == 1
  free(x1); free(ux);
#endif
  free(X2);
  free(y); free(z);
  free(ul); freedist(ul_dist);  
}


void updateStrategy(int j, double rSX)
{
  group *g = Polity->g+j;
  Leader *ld = g->lead;
  int i, v;
  double y, z;
  y = ld->y;
  z = ld->z;
#if UPDATE_LEAD_PUN_EFFORT || UPDATE_LEAD_NORM_EFFORT
  v = drand(VLdist);
  if(v == 1){
    updateStrategyLeader_V1(j);
  }
  else if(v == 2){
    updateStrategyLeader_V2(j);
  }
  else if(v == 3){
    unsigned int *x0 = malloc(n*sizeof(unsigned int));
    for(i = 0; i < n; i++){
      x0[i] = (g->com+i)->x;
    }
    updateStrategyLeader_V3(j, x0, g->X, rSX);        // note: elements' value of x0 changes
    free(x0);
  }
#endif
  // update commoners
#if UPDATE_COM
  for(i = 0; i < n; i++){
    v = drand(VCdist);
    if(v){
      updateStrategyCommoner(v, j, i, y, z, rSX);
    }
  }
#endif
}

// play game
void playGame()
{
  int j;
  double rSX;
  group *g;
  // update X for every group in this round
  for(rSX = 0, j = 0; j < G; j++){
    g = Polity->g+j;
    calcX(j);                       // calculates and assigns X 
//     r = 1.0 + (e*g->lead->z)/(n*z_0 + g->lead->z);
//     rSX += r*g->X;
    rSX += g->X;
  }
  Polity->SX = rSX;                  // sum of group efforts for this round  
  // update group production for every group in this round
  for(j = 0; j < G; j++){    
    calcP(j, rSX);                   // updates group production value for each group
    
    calcPayoff(j);                  // updates payoff of groups
    
    calcUtilityFunction(j, rSX);
    updateStrategy(j, rSX);             // updates strategy for group j    
  }
  
}

// initialize variables
void init(int r, int Runs)	//SG
{
  group *g;
  Leader *ld;
  Commoner *com;
  int i, j;
  for(j = 0; j < G; j++){
    g = Polity->g+j;
    ld = g->lead;
    for(i = 0; i < n; i++){
      com = g->com+i;
      com->x = INIT_COM_EFFORT;
      com->pi = 0;
      com->uc = 0;
    }
    if (UPDATE_LEAD_PUN_EFFORT) {
      ld->y = (r+U01())/Runs;
    }
    else{
      ld->y = INIT_LEAD_PUN_EFFORT;  
    }
    // SG 
    // ld->y = INIT_LEAD_PUN_EFFORT;
    ld->z = INIT_LEAD_NORM_EFFORT;
    ld->pi = 0;
    ld->ul = 0;
    g->nc = n;
    g->theta = Theta;
    g->X = 0;
    g->P = 0;
  }
  Polity->SX = 0;
  Polity->ng = G;
  // initialize strategy update method probability distribution
  VCdist->p[0] = 1.0-Vc1-Vc2-Vc3;
  VCdist->p[1] = Vc1; 
  VCdist->p[2] = Vc2;
  VCdist->p[3] = Vc3;
  initdist(VCdist, 1.0);  
  
  VLdist->p[0] = 1.0-Vl1-Vl2-Vl3;
  VLdist->p[1] = Vl1; 
  VLdist->p[2] = Vl2;
  VLdist->p[3] = Vl3;
  initdist(VLdist, 1.0); 
  Upc = Vc1 + Vc2 + Vc3;
  
}

int main(int argc, char **argv)
{
#if DEBUG
  feenableexcept(FE_DIVBYZERO| FE_INVALID|FE_OVERFLOW); // enable exceptions
#endif
  if(argc ^ 2){
    printf("Usage: ./prop prop.config\n");
    exit(1);
  }
  if(read_config(argv[1])){           // read config
    printf("READDATA: Can't process %s \n", argv[1]);
    return 1;
  }    
  if(TestPlot && TestG >= G){
    printf("\n Program terminated! Error: TestG (test group index) must be less than G. Note: [index starts from 0] \n");
    exit(1);
  }
  if(s0 > 1.0){
    printf("s0 value should be less than 1.0\n");
    exit(1);
  }
  if((s0+s1) > 1.0){
    printf("s0+s1 value should be less than 1.0\n");
    exit(1);
  }
  LOG_FLT_MAX = log(FLT_MAX);
  
#if NB_MULTI_SIM == 1
  int nsim, bsim, vsim;
  for(nsim = 0; nsim < (sizeof(Nsim)/sizeof(int)); nsim++){
    n = Nsim[nsim];
    
    for(bsim = 0; bsim < (sizeof(Bsim)/sizeof(double)); bsim++){
      b = Bsim[bsim];
      for(vsim = 0; vsim < (sizeof(Vsim)/sizeof(int)); vsim++){
	if(Vsim[vsim] == 1) {
	  Vc1 = 0.01; Vc2 = 0.24; Vc3 = 0.0;
	  Vl1 = 0.01; Vl2 = 0.24; Vl3 = 0.0;
	}
	else if(Vsim[vsim] == 2) {
	  Vc1 = 0.00; Vc2 = 0.01; Vc3 = 0.24;
	  Vl1 = 0.01; Vl2 = 0.24; Vl3 = 0.00;
	}
	else if(Vsim[vsim] == 3) {
	  Vc1 = 0.00; Vc2 = 0.01; Vc3 = 0.24;
	  Vl1 = 0.00; Vl2 = 0.01; Vl3 = 0.48;
	}
	
#endif
  xvmean = 0; yvmean = 0; zvmean = 0;    // standard deviation
  xmlv = 0; ymlv = 0; zmlv  = 0;          // standard deviation of average values for 10 runs
  initrand(Seed);
  allocStatVar();                                             // allocate memory for global statistic variables
  
#if TestPlot
  FILE *fgr;
  group *tG;  
  double s_pi;
#endif  
  
  int r, i, j;  
#if !ALLDATAFILE
  char xdata[200], str[30];
#endif
  time_t now;     
  // print headers for values to be displayed on std output
#if !CLUSTER
  //printf("\nValues:\t x\t y\t z\t pi_0\t pi_1\t uc\t ul\t P\t seed\n");
  printf("\nValues:\t x\t y    (st.d)\t z    (st.d)\t pi_0\t pi_1\t uc\t ul\t P\t seed\n");		//SG
#endif
  
  Theta = Theta_a/(n + Theta_a);                              // modified implementation of tax Theta
  
  for( r = 0; r < Runs; r++){                                 // through all sets of runs    
    // initialize pseudorandom generator with seed
    if(Seed == 0){      
      now = time(0);
      Seed_i = ((unsigned long)now) + r;  
    }
    else{
      Seed_i = Seed;
    }
    initrand(Seed_i);      
    // uncomment line below to get seed in case you got core dump and want to get seed to replicate bug
    //printf("\n run# %d seed: %lu\n",r+1, Seed_i);
    
    // setup and initialize variables
    setup();                                                  // allocates memory for all polity system variables
    init(r, Runs);                                                   // intialize polity system state values
    calcStat(0, r);
#if TestPlot
    tG = Polity->g + TestG;
    fgr = fopen("testgroup.dat", "w");
    fprintf(fgr, "T X/n y z P Pi_c Pi_l\n");    
    for(s_pi = 0, j = 0; j < n; j++) s_pi += tG->com[j].pi;
    fprintf(fgr, "%d %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n", 0, (tG->X)/n, (tG->lead)->y, (tG->lead)->z, tG->P, s_pi/n, tG->lead->pi);
#endif
    for( i = 0; i < T; i++){                                  // through all time points of simulation
      playGame();                                             // play us vs nature and us vs them game and update strategies
      if( (i+1) % SKIP == 0){                                 // every SKIP time, take snapshot of states of traits	
	calcStat((i+1)/SKIP, r);                              // calculate statistics and write individual runs data to file
#if TestPlot	
	for(s_pi = 0.0, j = 0; j < n; j++) s_pi += tG->com[j].pi;
	fprintf(fgr, "%d %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf\n", (i+1)/SKIP, (tG->X)/n, (tG->lead)->y, (tG->lead)->z, tG->P, s_pi/n, tG->lead->pi);
#endif
      }   
    }
    calcStat(-1, -1);                                        // free file pointers for individual run data files
#if TestPlot	
    fclose(fgr);	
#endif

#if !CLUSTER
    
  #if !AVERAGE_GRAPH_ONLY
      if(Runs > 1)
	//plotallIndividualRun(r, 0); 
  #endif
      
  #if GRAPHS
      if(Runs > 1)
	//plotallIndividualRun(r, 1);
  #endif
      // remove individual datafiles if ALLDATAFILE set to 0
  #if !ALLDATAFILE        
	sprintf(str, "x%d.dat", r); prep_file(xdata, str); remove(xdata);
	sprintf(str, "p%d.dat", r); prep_file(xdata, str); remove(xdata);
	sprintf(str, "u%d.dat", r); prep_file(xdata, str); remove(xdata);
	sprintf(str, "gp%d.dat", r); prep_file(xdata, str); remove(xdata);           
  #endif
  #if TestPlot	
	testGPlot("testgroup.dat", r);              // plots test group dynamics	
  #endif
      
#endif
   
    cleanup();                                                // free all memory allocated for polity system
  }
  writeDataToFile();  
  clearStatVar();                                             // free other memory allocated for statistics variables

#if !CLUSTER
  //plotall(0);
#if GRAPHS
  //plotall(1);
#endif
#endif    
  #if TestPlot	
	remove("testgroup.dat");              // plots test group dynamics	
  #endif
#if NB_MULTI_SIM == 1
      }
    }
  }
#endif
  return 1;
}


/* 
gcc -Wall -O2 -march=native -pipe -o csi csi.c -lm
./csi csi.config
valgrind -v --track-origins=yes --leak-check=full --show-leak-kinds=all ./csi csi.config
gcc -g -o csi csi.c -lm
gdb csi
run csi.config

//profiling code
gcc -Wall -O2 -march=native -pipe -pg csi.c -o csi -lm
./csi csi.config
 gprof csi gmon.out > analysis.txt  
*/





