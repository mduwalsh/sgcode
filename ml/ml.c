#include <stdlib.h>
#include <stdio.h>
#include <math.h>
#include <time.h>
#include <unistd.h>

#define DEBUG 1
#if DEBUG
#define __USE_GNU
#include <fenv.h>       /* enable floating exceptions */
unsigned long Where;    // debugging counter
#endif

#include "rand.c"

#define TestLeader 0              // turn on or off printing strategy and payoff for leader
#define TestCommoner 0            // turn on or off printing strategy and payoff for Commoner
#define TestH 0                  // test polity interested in
#define TestG 0                  // test group interested in
#define TestI 2                   // test individual interested in

#define CLUSTER 1            // is this simulation on cluster yes or no (1 or 0)
#define SKIP 10           // time interval between snapshot of states
#define STU 200             // summary time period for simulation

#define AVERAGE_GRAPH_ONLY 0  // if 1, generate only average graphs and supress individual run graphs
#define TURNOFF_TEST_PLOT 1 // if 1, does not produce frequency distribution and traits graphs for test group

#define ALLDATAFILE 1        // if 1, generate all data files for individual runs and summary too 
#define GRAPHS      0        // if 1, saves graphs as png files 

#define TRAITS 1         // no. of strategic traits
#define EVENTS 3         // no. of events type

#define INIT_COM_EFFORT U01()*1.0
#define INIT_LEAD_EFFORT U01()*1.0
#define INIT_CHIEF_EFFORT U01()*1.0
// change values to 0 if update strategy for any one of commoner or lead or chief is to be turned off
#define UPDATE_COM 1
#define UPDATE_LEAD_EFFORT 0
#define UPDATE_CHIEF_EFFORT 1

//#define INIT_THETA 0.5
//#define INIT_ETA 0.5
#define UPDATE_LEAD_TAX 0
#define UPDATE_CHIEF_TAX 0


#define GP 1 // 1, 2, 3; uses different calculation for group production
#define PS 1 // 1, 2, 3; uses different calculation for polity strength

FILE *fgr;

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

typedef struct
{
  double        pi;           // payoff commoner
  double        t;            // tax
  double        x[TRAITS];    // trait vector   
} individual;

typedef struct
{
  double        p;            // group production
  double        X;            // sum of efforts from commoners in group
  double        theta;        // tax on commoners in a group
  int           nc;           // number of commoners in a group
  int           nl;           // number of leaders in a group; //for now it is 1
  individual    *lead;        // leaders
  individual    *com;         // group members / commoners; 
} group;

typedef struct
{
  int          ng;            // no. of groups in a polity   
  double       eta;           // tax on leaders
  double       S;             // polity strength
  group        *g;            // polity members  
  individual   *chief;        // chief;
  
} polity;

polity *Polity;             // "system structure"

// Global variables
/** configuration parameters **/
unsigned Seed;                    // Seed: Seed of psuedorandom number passed
unsigned long Seed_i;             // Seed_i: seed of each run
int Runs;                         // no. of runs
int G;                            // no. of groups
int N;                            // no. of individuals
int H;                            // no. of politites
double K;                         // benefit factor
int T;                            // max time of simulation
double Theta, Eta;                // taxes
double C0, C1, C2, C3, C4;        // cost parameter
double B0, B1, B2;                // expected benefit for a game
double Sigma;                     // standard deviation for distribution of mutation for contribution
double Sigma_t;                   // standard deviation for distribution of mutation for tax levels
double V1, V2, V3;                // probability with which participants choose method of changing efforts and tax levels
double x_0, y_0, z_0, X0;            // x_0: half effort of commoner; y_0: half effort of leader; z_0: half effort of chief; X0: half strength of polity
double e, E;                      // e: leader efficiency; E: chief efficiency


// statistical variables
double (*Xmean)[TRAITS];          // average effort by individuals
double (*Ymean)[TRAITS];            // average effort by leaders
double (*Zmean)[TRAITS];            // average effort by chiefs
double *Pi0mean;                  // average payoff of commoners
double *Pi1mean;                  // average payoff of leaders
double *Pi2mean;                  // average payoff of chiefs
double *Tmean;                    // average tax on commoners
double *Emean;                    // average tax on leaders
double *GPmean;                   // group production
double *Smean;                    // average polity strength
double *Qmean;                    // production to leader per each commoner
double *Umean;                    // production to chief per each commoner

dist *Mdist;                      // 'dist' is defined in rand.c; Probability distribution of choosing method of changing efforts and tax levels 

//FILE *ftp;

void prep_file(char *fname, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "g%02dn%02dk%0.2fb2%.2ftheta%.2feta%.2fc0%.2fc1%.2fc2%.2fc3%.2fc4%.2fv1%.2fv2%.2fv3%.2fX0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s", G, N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, V1, V2, V3, X0, e, E, x_0, y_0, z_0, apndStr);
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
inline void merge (double *a, int n, int m) {
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
 
inline void merge_sort (double *a, int n) {
    if (n < 2)
        return;
    int m = n>>1;   // divide by 2
    merge_sort(a, m);
    merge_sort(a + m, n - m);
    merge(a, n, m);
}



#define EXPECT(a,b,c) if ((a) != fscanf(f, b"%*[ ^\n]\n",c)){ fclose(f); printf("Error: %s\n",b); return 1; }

int read_config(char *file_name)
{
  FILE *f;

  if (!(f = fopen(file_name,"r"))) return 1;
  EXPECT(1, "unsigned Seed      = %u;", &Seed);
  EXPECT(1, "int      Runs      = %d;", &Runs);
  EXPECT(1, "int      T         = %d;", &T);
  
  EXPECT(1, "int      N         = %d;", &N);
  EXPECT(1, "int      G         = %d;", &G);  
  EXPECT(1, "int      H         = %d;", &H);
  
  EXPECT(1, "double   K         = %lf;", &K);  
  EXPECT(1, "double   B2        = %lf;", &B2); 
  EXPECT(1, "double   Theta     = %lf;", &Theta);
  EXPECT(1, "double   Eta       = %lf;", &Eta);
  
  EXPECT(1, "double   C0        = %lf;", &C0);
  EXPECT(1, "double   C1        = %lf;", &C1);
  EXPECT(1, "double   C2        = %lf;", &C2);   
  EXPECT(1, "double   C3        = %lf;", &C3);
  EXPECT(1, "double   C4        = %lf;", &C4);   

  EXPECT(1, "double   V1        = %lf;", &V1);
  EXPECT(1, "double   V2        = %lf;", &V2);
  EXPECT(1, "double   V3        = %lf;", &V3); 
  
  EXPECT(1, "double   x0        = %lf;", &x_0);
  EXPECT(1, "double   y0        = %lf;", &y_0);
  EXPECT(1, "double   z0        = %lf;", &z_0);
  EXPECT(1, "double   X0        = %lf;", &X0);  
  
  EXPECT(1, "double   e         = %lf;", &e);
  EXPECT(1, "double   E         = %lf;", &E);
  
  EXPECT(1, "double   Sigma     = %lf;", &Sigma);
  EXPECT(1, "double   Sigma_t   = %lf;", &Sigma_t);    
  
  /*EXPECT(1, "double   B0        = %lf;", &B0);
  EXPECT(1, "double   B1        = %lf;", &B1);
  EXPECT(1, "double   B2        = %lf;", &B2);  
  */
  

  fclose(f);   
  
  if (Runs < 1) exit(1);
  return 0;
}
#undef EXPECT

// allocate memory for Polity system
void setup()
{
  int h, j;
  group *gr;
  // allocate memory for Polity and groups and commoners
  Polity = malloc(H*sizeof(polity));
  for(h = 0; h < H; h++){                                            // through all polities
    (Polity+h)->chief = malloc(sizeof(individual));                 // allocate memory for chief
    (Polity+h)->g = malloc(G*sizeof(group));                          // allocate memory for groups in each polity    
    for(j = 0; j < G; j++){                                          // through all groups in a polity
      gr = (Polity+h)->g+j;                                             // get reference to group pointer
      gr->lead = malloc(1*sizeof(individual));                       // allocate memory for leader in a group
      gr->com = malloc(N*sizeof(individual));                          // allocate memory for commoners in a group
    }
  }
}

void allocStatVar()
{
  // allocate memory to store statistics points of snapshot along time period of simulaltion
  Xmean = calloc((int)(T/SKIP+1), sizeof(double [TRAITS]));         // for effort by commoners
  Ymean = calloc((int)(T/SKIP+1), sizeof(double [TRAITS]));         // for effort by leaders
  Zmean = calloc((int)(T/SKIP+1), sizeof(double [TRAITS]));         // for effort by chiefs
  Pi0mean = calloc((int)(T/SKIP+1), sizeof(double));                   // for payoff of commoners
  Pi1mean = calloc((int)(T/SKIP+1), sizeof(double));                   // for payoff of leaders
  Pi2mean = calloc((int)(T/SKIP+1), sizeof(double));                   // for payoff of chiefs
  Tmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for tax paid by commoners to leaders
  Emean = calloc((int)(T/SKIP+1), sizeof(double));                    // for tax paid by leaders to chiefs
  GPmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for group production 
  Smean = calloc((int)(T/SKIP+1), sizeof(double));                    // for polity strength
  Qmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for group production 
  Umean = calloc((int)(T/SKIP+1), sizeof(double));                    // for polity strength
  
  Mdist = allocdist(4);                                              // allocate memory for distribution for update strategy method
}

void cleanup()                                                       // cleans Polity system memory
{
  int h, j;
  group *gr;
  for(h = 0; h < H; h++){                                            // through all polities
    for(j = 0; j < G; j++){                                          // through all groups in a polity
      gr = (Polity+h)->g+j;                                          // get reference to group pointer
      free(gr->lead); 
      free(gr->com);
    }
    free((Polity+h)->g);
    free((Polity+h)->chief);    
  }
  free(Polity);
}

void clearStatVar()
{
  free(Xmean); free(Ymean); free(Zmean); free(Pi0mean); free(Pi1mean); free(Pi2mean); 
  free(Tmean); free(Emean); free(GPmean); free(Smean); free(Qmean); free(Umean);
  freedist(Mdist);  
}

/***************************************Statistics calcualtions and writing to file*******************************************************/
/*****************************************************************************************************************************************/

// calculates all the stats
void calcStat(int k, int r) 
{
  double xm = 0, ym = 0, zm = 0, p0m = 0, p1m = 0, p2m = 0, tm = 0, em = 0, gpm = 0, sm = 0, um = 0, qm =  0, tbp = 0;
  int h, j, i, ncm = 0, nl = 0, ncf = 0;  
  polity *p;
  group *g;
  individual *ind;  
#if !CLUSTER
  FILE *fx, *fy, *fz, *fp0, *fp1, *fp2, *ft, *fe;
  char fname[100];
  
  if(k== (T/SKIP)){                                // data files for freq. distribution
    sprintf(fname, "x_%d.dat", r);
    fx = fopen(fname, "w");   
    sprintf(fname, "y_%d.dat", r);
    fy = fopen(fname, "w");
    sprintf(fname, "z_%d.dat", r);
    fz = fopen(fname, "w");
    sprintf(fname, "t_%d.dat", r);
    ft = fopen(fname, "w");
    sprintf(fname, "e_%d.dat", r);
    fe = fopen(fname, "w");
    sprintf(fname, "p0_%d.dat", r);
    fp0 = fopen(fname, "w");
    sprintf(fname, "p1_%d.dat", r);
    fp1 = fopen(fname, "w");
    sprintf(fname, "p2_%d.dat", r);
    fp2 = fopen(fname, "w");
    
  }
#endif

  static FILE **fp = NULL;  // file pointers for individual run  
  if(k < 0){
    for (h = 5; h--; fclose(fp[h]));
    free(fp); 
    return;
  }     
  
  if(!k){   
    char ixdata[200], ipdata[200], itdata[200], igpdata[200], isdata[200], tstr[100];
    fp = malloc(5*sizeof(FILE *));    
    sprintf(tstr, "x%d.dat", r); prep_file(ixdata, tstr);    
    sprintf(tstr, "p%d.dat", r); prep_file(ipdata, tstr);    
    sprintf(tstr, "t%d.dat", r); prep_file(itdata, tstr);    
    sprintf(tstr, "gp%d.dat", r); prep_file(igpdata, tstr); 
    sprintf(tstr, "s%d.dat", r); prep_file(isdata, tstr);
    fp[0] = fopen(ixdata, "w");
    fp[1] = fopen(ipdata, "w");
    fp[2] = fopen(itdata, "w");  
    fp[3] = fopen(igpdata, "w");
    fp[4] = fopen(isdata, "w");
    
    // write headers   
    for( j = 0; j < 3; j++){
      fprintf(fp[j], "%d\t", 0); 
      if(j < 2){
	fprintf(fp[j], "com\tlead\tchief\n");
      }
      else{
	fprintf(fp[j], "theta\teta\n");
      }
    }     
    fprintf(fp[3], "0\tP\tQ\n");
    fprintf(fp[4], "0\tS\tU\n");
  }
 
  // calculate stat of traits
  for(h = 0; h < H; h++){                                                               // through all polities
    p = Polity+h;
    for(j = 0; j < p->ng; j++){                                                         // through all groups
      g = p->g+j;
      for(i = 0; i < g->nc; i++){                                                       // through all commoners
	ind = g->com+i;
	xm += ind->x[0];                                                                // sum of commoner's effort
	p0m += ind->pi;                                                                 // sum of commoner's payoff
	#if !CLUSTER
	if(k== (T/SKIP)){
	  fprintf(fx, "%.6lf\n", ind->x[0]);
	  fprintf(fp0, "%.6lf\n", ind->pi);
	}
#endif
      }
      ncm += g->nc;                                                                     // sum of commoner's number
      for(i = 0; i < g->nl; i++){                                                       // through all leaders
	ind = g->lead+i;
	ym += ind->x[0];                                                                // sum of leader's efforts
	p1m += ind->pi;                                                                 // sum of leader's payoff
#if !CLUSTER
	if(k== (T/SKIP)){
	  fprintf(fy, "%.6lf\n", ind->x[0]);
	  fprintf(fp1, "%.6lf\n", ind->pi);
	}
#endif
      }
      tm += g->theta;                                                                   // sum of tax factor theta on commoners
#if !CLUSTER
      if(k== (T/SKIP)){
	fprintf(ft, "%.6lf\n", g->theta);
      }
#endif
      nl += g->nl;                                                                     // sum of no. of leaders
      gpm += g->p;                                                                     // sum P over all groups
      tbp = g->theta*g->p;                                             	       // theta*b*P over of a group
      qm += (1.0-p->eta)*tbp;                                             	       // sum (1-eta)*theta*b*P over all groups
      um += p->eta*tbp;                                                                // sum eta*theta*b*P over all groups      
    }
    zm += (p->chief)->x[0];                                                               // sum of chief's efforts
    p2m += (p->chief)->pi;
    em += p->eta;                                                                       // sum of tax on leaders by chiefs in polity
    sm += p->S;                                                                         // sum of polity strength in each polity
    ncf++;                                                                              // sum of number of chiefs in system
    
#if !CLUSTER
    if(k== (T/SKIP)){
      fprintf(fz, "%.6lf\n", (p->chief)->x[0]);
      fprintf(fp2, "%.6lf\n", (p->chief)->pi);
      fprintf(fe, "%.6lf\n", p->eta);
    }
#endif
  }
  // compute averages of traits
  xm /= ncm; 
  p0m /= ncm;
  ym /= nl;
  p1m /= nl;
  tm /= nl;
  zm /= ncf;
  p2m /= ncf;
  em /= ncf;  
  gpm /= (nl);        
  sm /= (ncf);      
  qm /= nl;
  um /= nl;   // nl = G*H
  // store in averages in variable as sum to average it after multiple runs
  Xmean[k][0] += xm;
  Ymean[k][0] += ym;
  Zmean[k][0] += zm;
  Pi0mean[k] += p0m;
  Pi1mean[k] += p1m;
  Pi2mean[k] += p2m;
  Tmean[k] += tm;
  Emean[k] += em;
  GPmean[k] += gpm;
  Smean[k] += sm;  
  Qmean[k] += qm;
  Umean[k] += um;
       
  // write data for individual runs
  fprintf(fp[0], "%d  %.4lf  %.4lf  %.4lf\n", k, xm, ym, zm);
  fprintf(fp[1], "%d  %.4lf  %.4lf  %.4lf\n", k, p0m, p1m, p2m);
  fprintf(fp[2], "%d  %.4lf  %.4lf\n", k, tm, em);  
  fprintf(fp[3], "%d  %.4lf  %.4lf\n", k, gpm, qm); 
  fprintf(fp[4], "%d  %.4lf  %.4lf\n", k, sm, um); 
#if !CLUSTER
  // print final values for each run
  if(k == T/SKIP){
    printf("run#%d \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %lu\n", r, xm, ym, zm, tm, em, p0m, p1m, p2m, gpm, sm, qm, um, Seed_i);  
    fclose(fx); fclose(fy); fclose(fz); fclose(ft); fclose(fe); fclose(fp0); fclose(fp1); fclose(fp2);
  }
  
#endif
}

void writeDataToFile()
{
  int j;
  char xdata[200], pdata[200], tdata[200], gpdata[200], sdata[200];
  FILE **fp = malloc(8*sizeof(FILE *));    
  prep_file(xdata, "x.dat");    
  prep_file(pdata, "p.dat");    
  prep_file(tdata, "t.dat");   
  prep_file(gpdata, "gp.dat");
  prep_file(sdata, "s.dat");
  fp[0]= fopen(xdata, "w");
  fp[1]= fopen(pdata, "w");
  fp[2]= fopen(tdata, "w");
  fp[6]= fopen(gpdata, "w");
  fp[7]= fopen(sdata, "w");
#if ALLDATAFILE
  int stu = (T-STU)/SKIP;
  double xsum = 0, ysum = 0, zsum = 0, p0sum = 0, p1sum = 0, p2sum =0, tsum = 0, esum = 0;
  char  xsumdata[200], psumdata[200], tsumdata[200];
  prep_file(xsumdata, "xsum.dat");    
  prep_file(psumdata, "psum.dat");    
  prep_file(tsumdata, "tsum.dat");  
#endif  
  // write headers
  for( j = 0; j < 3; j++){
    fprintf(fp[j], "%d\t", 0); 
    if(j < 2){
      fprintf(fp[j], "com\tlead\tchief\n");
    }
    else{
      fprintf(fp[j], "theta\teta\n");
    }
  }  
  fprintf(fp[6], "0  P  Q\n");
  fprintf(fp[7], "0  S  U\n");
  for(j = 0; j < (int)(T/SKIP)+1; j++){
    // average accumulated mean values for multiple runs
    Xmean[j][0] /= (double)Runs;
    Ymean[j][0] /= (double)Runs;
    Zmean[j][0] /= (double)Runs;
    Pi0mean[j] /= (double)Runs;
    Pi1mean[j] /= (double)Runs;
    Pi2mean[j] /= (double)Runs;
    Tmean[j] /= (double)Runs;
    Emean[j] /= (double)Runs;
    GPmean[j] /= (double)Runs;
    Smean[j] /= (double)Runs;
    Qmean[j] /= (double)Runs;
    Umean[j] /= (double)Runs;
    // write data to file
    fprintf(fp[0], "%d  %.4lf  %.4lf  %.4lf\n", j, Xmean[j][0], Ymean[j][0], Zmean[j][0]);
    fprintf(fp[1], "%d  %.4lf  %.4lf  %.4lf\n", j, Pi0mean[j], Pi1mean[j], Pi2mean[j]);
    fprintf(fp[2], "%d  %.4lf  %.4lf\n", j, Tmean[j], Emean[j]); 
    fprintf(fp[6], "%d  %.4lf  %.4lf\n", j, GPmean[j], Qmean[j]);
    fprintf(fp[7], "%d  %.4lf  %.4lf\n", j, Smean[j], Umean[j]);    
    
#if ALLDATAFILE
    if(j < stu) continue;
    xsum += Xmean[j][0];
    ysum += Ymean[j][0];
    zsum += Zmean[j][0];
    p0sum += Pi0mean[j];
    p1sum += Pi1mean[j];
    p2sum += Pi2mean[j];
    tsum += Tmean[j];
    esum += Emean[j];
#endif
  }
  for(j = 0; j < 3; j++){
    fclose(fp[j]);
  }
  fclose(fp[6]);
  fclose(fp[7]);
#if ALLDATAFILE
  double st = STU/SKIP + 1;
  xsum /= st; ysum /= st; zsum /= st; p0sum /= st; p1sum /= st; p2sum /= st; tsum /= st; esum /= st;                // computing average over summary period time
  // write data to file
  fp[3]= fopen(xsumdata, "w");
  fp[4]= fopen(psumdata, "w");
  fp[5]= fopen(tsumdata, "w");
   // write headers
  for( j = 3; j < 6; j++){
    if(j < 5){
      fprintf(fp[j], "com\tlead\tchief\n");
    }
    else{
      fprintf(fp[j], "lead\tchief\n");
    }
  }   
  fprintf(fp[3], "%.4lf  %.4lf  %.4lf\n", xsum, ysum, zsum);
  fprintf(fp[4], "%.4lf  %.4lf  %.4lf\n", p0sum, p1sum, p2sum);
  fprintf(fp[5], "%.4lf  %.4lf\n", tsum, esum); 
  for(j = 3; j < 6; j++){
    fclose(fp[j]);
  }
#endif
  free(fp);
#if !CLUSTER
  // print averaged final values
  printf("\nAvg: \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf\n", Xmean[T/SKIP][0], Ymean[T/SKIP][0], Zmean[T/SKIP][0], Tmean[T/SKIP], Emean[T/SKIP], Pi0mean[T/SKIP], Pi1mean[T/SKIP], Pi2mean[T/SKIP], GPmean[T/SKIP], Smean[T/SKIP], Qmean[T/SKIP], Umean[T/SKIP]);  
#endif
}

void plotall(int m)
/*
 * m: 0 or 1; to save graphs as image file or not
 */
{
    // write data to file
  char xdata[200], pdata[200], tdata[200], gpdata[200], sdata[200];
  char title[200], xpng[200];
  int datacolumn = 3+1;
  prep_file(xdata, "x.dat"); prep_file(pdata,"p.dat"); prep_file(tdata, "t.dat"); prep_file(gpdata, "gp.dat"); prep_file(sdata, "s.dat");
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "n= %d, k=%.2f, b2=%.2f, theta=%.2f, eta=%.2f, c=%.1f%.1f%.1f%.1f%.1f, v=%d%d%d", N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, (int)(V1*10), (int)(V2*10), (int)(V3*10));
  
  //fprintf(gp, "set key outside vertical spacing 1 width 1 height 4\n");       
  fprintf(gp, "set key outside vertical height -1\n");        
  if(m){
    sprintf(xpng, "ept_n%02dk%.2fb2%.2ftheta%.2feta%.2fc%.1f%.1f%.1f%.1f%.1fv%d%d%d.png", N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, (int)(V1*10), (int)(V2*10), (int)(V3*10)); 
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
  fprintf(gp, "stats '%s' using 4 prefix 'D' nooutput\n", pdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'E' nooutput\n", gpdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'F' nooutput\n", sdata);
  
  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 4,1 title '%s' \n", title);    // set subplots layout
  fprintf(gp, "set label 1 'Average' at screen 0.12,0.98 center font \",13\"\n");  // label to indicate average
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
  /*fprintf(gp, "set ylabel 'strength' \n");
  fprintf(gp, "ymax = F_max\n");  
  fprintf(gp, "set yrange [-0.1:ymax+ymax/10+0.1]\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader \n", 2, sdata);
  */

  fprintf(gp, "unset multiplot \n");
  
  fflush(gp); 
  pclose(gp);
#if !CLUSTER
#if !ALLDATAFILE
  remove(xdata);
  remove(pdata);
  remove(tdata);
  remove(gpdata);
  remove(sdata);
#endif
#endif
}

void plotTraitsDist(int r)
{  
  char fx[100], fy[100], fz[100], ft[100], fe[100], fp0[100], fp1[100], fp2[100], title[200];  
  sprintf(fx, "x_%d.dat", r);
  
  sprintf(fy, "y_%d.dat", r);
  
  sprintf(fz, "z_%d.dat", r);
  
  sprintf(ft, "t_%d.dat", r);
  
  sprintf(fe, "e_%d.dat", r);
  
  sprintf(fp0, "p0_%d.dat", r);
  
  sprintf(fp1, "p1_%d.dat", r);
  
  sprintf(fp2, "p2_%d.dat", r);
  
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "traits distribution");  
  fprintf(gp, "set key outside vertical\n");        
  fprintf(gp, "set term x11 dashed %d \n", 0);
 
  
  fprintf(gp, "stats '%s' using 1 prefix 'A' nooutput\n", fx);
  fprintf(gp, "stats '%s' using 1 prefix 'B' nooutput\n", fy);
  fprintf(gp, "stats '%s' using 1 prefix 'C' nooutput\n", fz);
  fprintf(gp, "stats '%s' using 1 prefix 'D' nooutput\n", fp0);
  fprintf(gp, "stats '%s' using 1 prefix 'E' nooutput\n", fp1);  
  fprintf(gp, "stats '%s' using 1 prefix 'F' nooutput\n", fp2);

 
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 3,3 title '%s' \n", title);    // set subplots layout
  fprintf(gp, "set ylabel 'Freq' \n"); 
  fprintf(gp, "set tics out nomirror \n");
  fprintf(gp, "set style fill solid 0.5\n"); //fillstyle
  
  fprintf(gp, "set autoscale x \n");
  fprintf(gp, "set xlabel 'x' \n");
  fprintf(gp, "n = 100\n"); //number of intervals
  fprintf(gp, "max = A_max\n");
  fprintf(gp, "min = A_min\n");
  fprintf(gp, "width=(max-min)/n \n");
  fprintf(gp, "if(max-min==0) set xrange[0:max+0.5]\n");
  fprintf(gp, "if(width==0) width = 0.01\n");
  fprintf(gp, "hist(x,width)=width*floor(x/width)+width/2.0\n");
  fprintf(gp, "set boxwidth width*0.9\n");
  fprintf(gp, "plot '%s' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'red' notitle\n", fx);   
  
  fprintf(gp, "set autoscale x \n");
  fprintf(gp, "set xlabel 'y' \n");
  fprintf(gp, "n = 100\n"); //number of intervals
  fprintf(gp, "max = B_max\n");
  fprintf(gp, "min = B_min\n");
  fprintf(gp, "width=(max-min)/n \n");
  fprintf(gp, "if(max-min==0) set xrange[0:max+0.5]\n");
  fprintf(gp, "if(width==0) width = 0.01\n");  
  fprintf(gp, "hist(x,width)=width*floor(x/width)+width/2.0\n");
  fprintf(gp, "set boxwidth width*0.9\n");
  
  fprintf(gp, "plot '%s' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'green' notitle\n", fy); 
  
  fprintf(gp, "set autoscale x \n");
  fprintf(gp, "set xlabel 'z' \n");
  fprintf(gp, "n = 100\n"); //number of intervals
  fprintf(gp, "max = C_max\n");
  fprintf(gp, "min = C_min\n");
  fprintf(gp, "width=(max-min)/n \n");
  fprintf(gp, "if(max-min==0) set xrange[0:max+0.5]\n");
  fprintf(gp, "if(width==0) width = 0.01\n");
  fprintf(gp, "hist(x,width)=width*floor(x/width)+width/2.0\n");
  fprintf(gp, "set boxwidth width*0.9\n");
  fprintf(gp, "plot '%s' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'blue' notitle\n", fz); 
  
  fprintf(gp, "set autoscale x \n");
  fprintf(gp, "set xlabel 'pi_0' \n");
  fprintf(gp, "n = 100\n"); //number of intervals
  fprintf(gp, "max = D_max\n");
  fprintf(gp, "min = D_min\n");
  fprintf(gp, "width=(max-min)/n \n");
  fprintf(gp, "if(max-min==0) set xrange[0:max+0.5]\n");
  fprintf(gp, "if(width==0) width = 0.01\n");
  fprintf(gp, "hist(x,width)=width*floor(x/width)+width/2.0\n");
  fprintf(gp, "set boxwidth width*0.9\n");
  fprintf(gp, "plot '%s' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'red' notitle\n", fp0);   
  
  fprintf(gp, "set autoscale x \n");
  fprintf(gp, "set xlabel 'pi_1' \n");
  fprintf(gp, "n = 100\n"); //number of intervals
  fprintf(gp, "max = E_max\n");
  fprintf(gp, "min = E_min\n");
  fprintf(gp, "width=(max-min)/n \n");
  fprintf(gp, "if(max-min==0) set xrange[0:max+0.5]\n");
  fprintf(gp, "if(width==0) width = 0.01\n");
  fprintf(gp, "hist(x,width)=width*floor(x/width)+width/2.0\n");
  fprintf(gp, "set boxwidth width*0.9\n");
  fprintf(gp, "plot '%s' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'green' notitle\n", fp1); 
  
  fprintf(gp, "set autoscale x \n");
  fprintf(gp, "set xlabel 'pi_2' \n");
  fprintf(gp, "n = 100\n"); //number of intervals
  fprintf(gp, "max = F_max\n");
  fprintf(gp, "min = F_min\n");
  fprintf(gp, "width=(max-min)/n \n");
  fprintf(gp, "if(max-min==0) set xrange[0:max+0.5]\n");
  fprintf(gp, "if(width==0) width = 0.01\n");
  fprintf(gp, "hist(x,width)=width*floor(x/width)+width/2.0\n");
  fprintf(gp, "set boxwidth width*0.9\n");
  fprintf(gp, "plot '%s' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'blue' notitle\n", fp2); 
  
  
  
  fprintf(gp, "set autoscale x \n");
  fprintf(gp, "set xlabel 'theta' \n");
  fprintf(gp, "n = 100\n"); //number of intervals
  fprintf(gp, "max = 1.0\n");
  fprintf(gp, "min = 0.0\n");
  fprintf(gp, "width=(max-min)/n \n");
  fprintf(gp, "hist(x,width)=width*floor(x/width)+width/2.0\n");
  fprintf(gp, "set boxwidth width*0.9\n");
  fprintf(gp, "set xrange [0:1]\n");
  fprintf(gp, "plot '%s' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'green' notitle\n", ft); 
  
  fprintf(gp, "set autoscale x \n");
  fprintf(gp, "set xlabel 'eta' \n");
  fprintf(gp, "n = 100\n"); //number of intervals
  fprintf(gp, "max = 1.0\n");
  fprintf(gp, "min = 0.0\n");
  fprintf(gp, "width=(max-min)/n \n");
  fprintf(gp, "hist(x,width)=width*floor(x/width)+width/2.0\n");
  fprintf(gp, "set boxwidth width*0.9\n");
  fprintf(gp, "set xrange [0:1]\n");
  fprintf(gp, "plot '%s' u (hist($1,width)):(1.0) smooth freq w boxes lc rgb'blue' notitle\n", fe); 
  
  fprintf(gp, "unset multiplot \n");
  
  fflush(gp); 
  pclose(gp);
}

void plotallIndividualRun(int r, int m)
/*
 * r: run
 * m: 0 or 1; to save graphs as image file or not
 */
{
    // write data to file
  char xdata[200], pdata[200], tdata[200], gpdata[200], sdata[200];
  char title[200], xpng[200], str[100];
  int datacolumn = 3+1;
  sprintf(str, "x%d.dat", r); prep_file(xdata, str);
  sprintf(str, "p%d.dat", r); prep_file(pdata, str);
  sprintf(str, "t%d.dat", r); prep_file(tdata, str); 
  sprintf(str, "gp%d.dat", r); prep_file(gpdata, str); 
  sprintf(str, "s%d.dat", r); prep_file(sdata, str); 
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "n= %d, k=%.2f, B2=%.2f theta=%.2f, eta=%.2f, c=%.1f%.1f%.1f%.1f%.1f, v=%d%d%d", N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, (int)(V1*10), (int)(V2*10), (int)(V3*10));
  
  fprintf(gp, "set key outside vertical height -1\n");        
  if(m){ // save graphs as file
    sprintf(xpng, "ept_n%02dk%.2fb2%.2ftheta%.2feta%.2fc%.1f%.1f%.1f%.1f%.1fv%d%d%d_%d.png", N, K, B2, Theta, Eta, C0, C1, C2, C3, C4, (int)(V1*10), (int)(V2*10), (int)(V3*10), r); 
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
  fprintf(gp, "stats '%s' using 4 prefix 'D' nooutput\n", pdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'E' nooutput\n", gpdata);  
  fprintf(gp, "stats '%s' using 2:3 prefix 'F' nooutput\n", sdata);

  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 4,1 title '%s' \n", title);    // set subplots layout
  fprintf(gp, "set label 1 'run: %d' at screen 0.13,0.98 center font \",13\"\n", r);  // label to indicate individual run index
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
  /*
  gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "theta plot");  
  fprintf(gp, "set key outside vertical\n");   
  fprintf(gp, "set term x11 dashed %d \n", 2);
  */

}

void plotLines(FILE *gp, char *data, int col, char *title)
{   
  //fprintf(gp, "set key outside vertical spacing 1 width 1 height 4\n");       
  fprintf(gp, "set key outside vertical\n");   
  fprintf(gp, "set term x11 dashed %d \n", 0);
    
  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title '%s'\n", "");
  /*fprintf(gp, "set multiplot layout 3,1 title '%s' \n", title);    // set subplots layout
  fprintf(gp, "set ylabel 'efforts' \n");  
  fprintf(gp, "set autoscale y\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader\n", 4, data);
  
  fprintf(gp, "set ylabel 'tax' \n");  
  fprintf(gp, "set autoscale y\n");
  fprintf(gp, "plot for [col=5:%d] '%s' using 1:col with lines lw 2 lt col-3 title columnheader\n", 6, data);
  
  fprintf(gp, "set ylabel 'payoff' \n");  
  fprintf(gp, "set autoscale y\n");
  fprintf(gp, "plot for [col=7:%d] '%s' using 1:col with lines lw 2 lt col-5 title columnheader\n", 8, data);
  
  fprintf(gp, "unset autoscale y\n");
  fprintf(gp, "unset multiplot\n");
  */
  fprintf(gp, "set title '%s'\n", "traits");
  fprintf(gp, "set ylabel 'X & theta & eta' \n");  
  fprintf(gp, "set autoscale y\n");
  fprintf(gp, "plot '%s' using 1:2 with lines lw 2 title columnheader, ", data);
  fprintf(gp, "'%s' using 1:5 with lines lw 2 title columnheader\n", data);
  //fprintf(gp, "'%s' using 1:6 with lines lw 2 title columnheader\n", data);
  fprintf(gp, "unset autoscale y\n");
}

void plotTest(char *data, int col, char *title)
{   
  FILE *gp = popen("gnuplot -persistent", "w");
  //fprintf(gp, "set key outside vertical spacing 1 width 1 height 4\n");       
  fprintf(gp, "set key outside vertical\n");   
  fprintf(gp, "set term x11 dashed %d \n", 0);
  fprintf(gp, "set lmargin at screen 0.1 \n");
  fprintf(gp, "set rmargin at screen 0.86 \n");
  //set tmargin at screen 0.97
  //set bmargin at screen 0.6  
  fprintf(gp, "stats '%s' using 2:3 prefix 'A' nooutput\n", data);
  fprintf(gp, "stats '%s' using 4 prefix 'B' nooutput\n", data);
  fprintf(gp, "stats '%s' using 5:6 prefix 'C' nooutput\n", data);
  fprintf(gp, "stats '%s' using 7:8 prefix 'D' nooutput\n", data);
  fprintf(gp, "stats '%s' using 9 prefix 'E' nooutput\n", data);
  fprintf(gp, "stats '%s' using 10:11 prefix 'F' nooutput\n", data);
  fprintf(gp, "stats '%s' using 12 prefix 'G' nooutput\n", data);
  
  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title '%s'\n", "");
  fprintf(gp, "set multiplot layout 4,1 title 'H:%d, G:%d, %s' \n", TestH, TestG, title);    // set subplots layout
  fprintf(gp, "set ylabel 'efforts' \n");  
  fprintf(gp, "ymax = A_max_x\n");
  fprintf(gp, "if(A_max_y > ymax) ymax = A_max_y \n");
  fprintf(gp, "if(B_max > ymax) ymax = B_max \n");
  fprintf(gp, "set yrange[0:ymax]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3\n");
  fprintf(gp, "plot for [col=2:%d] '%s' using 1:col with lines lw 2 title columnheader\n", 4, data);
  
  fprintf(gp, "set ylabel 'tax' \n");  
  fprintf(gp, "ymax = C_max_x\n");
  fprintf(gp, "if(C_max_y > ymax) ymax = C_max_y \n");
  fprintf(gp, "set yrange[0:ymax]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3\n");
  fprintf(gp, "plot for [col=5:%d] '%s' using 1:col with lines lw 2 lt col-3 title columnheader\n", 6, data);
  
  fprintf(gp, "set ylabel 'payoff' \n");  
  fprintf(gp, "ymax = D_max_x\n");
  fprintf(gp, "if(D_max_y > ymax) ymax = D_max_y \n");
  fprintf(gp, "if(E_max > ymax) ymax = E_max \n");
  fprintf(gp, "set autoscale y\n");
  fprintf(gp, "set yrange[:ymax]\n");
  fprintf(gp, "set format y \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3\n");
  fprintf(gp, "plot for [col=7:%d] '%s' using 1:col with lines lw 2 lt col-6 title columnheader\n", 9, data);
  
  fprintf(gp, "set ylabel 'P & Q & U' \n");  
  //fprintf(gp, "set y2label 'U'\n");
  fprintf(gp, "ymax = F_max_x\n");
  //fprintf(gp, "y2max = G_max\n");
  fprintf(gp, "if(F_max_y > ymax) {ymax = F_max_y}\n");
  fprintf(gp, "if(G_max > ymax) {ymax = G_max}\n");
  //fprintf(gp, "if(G_max > y2max) {y2max = G_max}\n");
  fprintf(gp, "if(ymax < 0.05) {ymax = 1.0}\n");
  //fprintf(gp, "if(y2max < 0.05) {y2max = 1.0}\n");
  fprintf(gp, "set yrange [0:ymax+0.01]\n");
  //fprintf(gp, "set y2range [0:y2max+0.01]\n");  
  fprintf(gp, "set format y \"%%.2f\"\n");
  //fprintf(gp, "set format y2 \"%%.2f\"\n");
  fprintf(gp, "set ytics ymax/3 nomirror\n");
  //fprintf(gp, "set y2tics y2max/3 nomirror \n");
  //fprintf(gp, "plot for [col=10:%d] '%s' using 1:col with lines lw 2 lt col-9 title columnheader axes x1y1,  for [col=12:%d] '%s' using 1:col with lines lw 2 lt col-9 title columnheader axes x1y2 \n", 11, data, 12, data);  
  fprintf(gp, "plot for [col=10:%d] '%s' using 1:col with lines lw 2 lt col-9 title columnheader \n", 12, data); 
  //fprintf(gp, "unset y2tics\n set ytics mirror\n");
  //fprintf(gp, "unset y2label\n");
  fprintf(gp, "unset multiplot\n");
  
  fflush(gp);
  pclose(gp);  
}

/*****************************************************************************************************************************************/
/*****************************************************************************************************************************************/


/*****************************************************************************************************************************************/
/****************************************************methods implementing model***********************************************************/

double groupProduction(double sx, double y, int n)
/*
 * sx: sum of efforts by commoner
 * y: effort by leader
 * n: no. of commoners in group / group size
 */
{   
#if GP==1  
  //return y*B0*sx;
  double r = ( 1.0 + (e*y)/(n*y_0+y) );
  return sx/( sx + (n*x_0)/r );
#endif
#if GP==2  
  return pow(sx,y);
#endif
#if GP==3  
  return pow(sx,y);               // sx here is sum of (x)^(1/y)
#endif
}

double polityStrength(double sx, double z, int g)
/*
 * sx: sum of group production reaching to chief
 * z: effort by chief
 * g: polity size /  no. of groups in a polity
 */
{  
#if PS==1  
  return sx/( sx + (g*X0)/( 1.0 + (E*z)/( z+(g*z_0) ) ) );
#endif
#if PS==2  
  return pow(sx,z);
#endif
#if PS==3  
  return pow(sx,z);               // sx here is sum of (X)^(1/z)
#endif
}

// updates group production p
void updateGroupProduction(int h, int j)
{
  int i;
  double s;
  group *gr;
  gr = Polity[h].g+j;
#if GP==3
  double i_y = 1.0/gr->lead[0].x[0];
#endif
  for(s = 0, i = 0; i < gr->nc; i++){                       // through all commoners
#if GP==3
    s += pow(gr->com[i].x[0],i_y);     
#else
    s += gr->com[i].x[0];                                   // sum production efforts from commoners
#endif
  }
  gr->X = s;                                                // sum of efforts from commoners                                  
  gr->p = groupProduction(s, gr->lead[0].x[0], gr->nc);            // group production 
}

// updates payoff of commoners and leaders in group after production
void updatePayoffAfterProduction(int h, int j)
/*
 * h: polity
 * j: group
 */
{
  int i;
  double s;
  group *gr;
  gr = Polity[h].g+j;
  // commoner payoff : 1 + (1-theta)*b*P - c0*x
  s = 1.0 + (1.0 - gr->theta)*B0*gr->p;                                     // payoff benefit of each individual after getting share from production p after tax theta
  for(i = 0; i < gr->nc; i++){                                                     // through all commoners
    gr->com[i].pi = s - (C0*gr->com[i].x[0]);                           // payoff of commoner afer deduction of cost of effort made
  }
  // lead payoff : (1-eta)*theta*n*b*P - c1*y - c3*theta
  s = ( 1.0-Polity[h].eta ) * gr->theta * gr->nc * B0 * gr->p;                       // share leaders get after production
  
  for(i = 0; i < gr->nl; i++){                                                     // through all leaders in a group
    gr->lead[i].pi = s - C1*gr->lead[i].x[0] - C3*gr->theta;            // leaders payoff after deduction of cost of coordination effort they make
  }  
}

// updates polity strengths and payoff of chiefs and returns sum of polities strength
double updatePolityStrengthAndChiefPayoff()
{
  int h, j;
  polity *p;
  group *g;
  double s, t;
#if PS == 3
  double i_z; 
#endif
  // update polities' strength
  for(t = 0, h = 0; h < H; h++){                                       // through all polities in array    
    p = Polity+h;    
#if PS==3
    i_z = 1.0/p->chief[0].x[0];
    for( s = 0, j = 0; j < p->ng; j++){                                // through all groups in a polity
      g = p->g+j;      
      s += pow(g->theta * g->p * g->nc, i_z);                          // sum group production reaching to leaders in each group
    }    
    s = pow(p->eta*B0, i_z) * s;                                          // group production reaching to chief with 1/z exponent
#else    
    for( s = 0, j = 0; j < p->ng; j++){                                // through all groups in a polity
      g = p->g+j;      
      s += g->theta * g->p * g->nc;                                    // sum group production reaching to leaders in each group
    }    
    s = p->eta * s * B0;                                                    // group production reaching to chief
#endif    
    p->S = polityStrength(s, p->chief[0].x[0], p->ng);                        // polity strength 
    t += p->S;                                                         // sum of polities strength 
  }
  
  // update chiefs' payofffs  
  // chief payoff : b2*S/(sum of S / H) - c2*z - c4*eta
  if(t == 0.0){                                                        // if sum of polities strength is 0, then divide benefits
    s = B2;                                                            // shared benefit
    for(h = 0; h < H; h++){                                            // through all polities
      p = Polity+h;
      (p->chief)->pi = s - C2*(p->chief)->x[0] - C4*(p->eta);// payoff of chief after us vs them game       
    }
  }
  else{
    s = B2*H/t;
    for(h = 0; h < H; h++){                                              // through all polities
      p = Polity+h;
      (p->chief)->pi = s*p->S - C2*(p->chief)->x[0] - C4*(p->eta);// payoff of chief after us vs them game      
    }
  }
  
  return t;                                                            // return sum of polities strengths
}

// changes effort made by commoners
void changeStratCommoner(int v, int h, int j, int i, double theta, double y)
/*
 * v: type of strategy change method
 * h: polity
 * j: group
 * i: commoner index
 * theta: tax on commoners by leader before leader updates its tax
 * y: coordination effort made by leader before leader updates its effort
 */
{
  int k;
  double x;
  individual *cm;
  cm = Polity[h].g[j].com+i; 
  if(v == 1){                                                               // random mutation
    cm->x[0] = normal(cm->x[0], Sigma);
    cm->x[0] = MAX(cm->x[0], 0.0001);
    cm->x[0] = MIN(cm->x[0], 1.0/C0);
  }
  else if(v == 2){                                                          // copy from peer
    do{ k = rnd(Polity[h].g[j].nc); }while( k == i);                        // select another individual in group
    if(Polity[h].g[j].com[k].pi > cm->pi){                                  // if selected individual has more payoff copy strategies
      cm->x[0] = Polity[h].g[j].com[k].x[0];
    }
  }
  else if(v == 3){                                                          // myopic optimization
    double f, sx;
    x = normal(cm->x[0], Sigma);
    x = MAX(x, 0);
#if GP==3
    sx = Polity[h].g[j].X - pow(cm->x[0], 1.0/y) + pow(x,1.0/y);
#else
    sx = Polity[h].g[j].X - cm->x[0] + x;  
#endif
    // new sum of efforts    
    f = 1.0 + (1.0-theta)*B0*groupProduction(sx, y, Polity[h].g[j].nc) - C0*x;                  // new payoff; 
    
#if DEBUG
#if TestCommoner
    // printing expected and old payoffs for particular leader
    if(h == TestH && j == TestG && i == TestI)
      printf("Commoner H:%d g:%d i:%d \t expected payoff: %.3lf \told (real) payoff: %.3lf\t new x: %.3lf   old x: %.3lf \n", h, j, i, f, cm->pi, x, cm->x[0]);
#endif
#endif
    if(f > cm->pi){
      cm->x[0] = x;
#if DEBUG      
#if TestCommoner
      if(h == TestH && j == TestG && i == TestI)
	printf("new strategy chosen x: %.6lf expected payoff: %.3lf\n\n", cm->x[0], f);
#endif
#endif
      
    }
  }
}

// changes effort and taxes by leader
void changeStratLeader(int v, int h, int j, double eta)
/*
 * v: type of method to change
 * h: polity
 * j: group
 * eta: tax on leader by chief of polity
 */
{
  int k;
  double t, y;                                                           // stores mutated tax and effort
  individual *ld;
  group *gr;
  gr = Polity[h].g+j;
  ld = gr->lead+0;
  if(v == 1){                                                            // random mutation
#if UPDATE_LEAD_TAX
    gr->theta = normal(gr->theta, Sigma_t);                              // change tax level by normal distribution with deviation Sigma_t
    gr->theta = MAX(gr->theta, 0.0);
    gr->theta = MIN(gr->theta, 1.0);
#endif
#if UPDATE_LEAD_EFFORT
    ld->x[0] = normal(ld->x[0], Sigma);
    ld->x[0] = MAX(ld->x[0], 0.0);    
    //ld->x[0] = MIN((1.0-eta)*(B0*ld->x[0]*gr->X)*gr->theta/C1, ld->x[0]);  
#endif
  }
  else if(v == 2){                                                       // selective copying
    if(Polity[h].ng > 1){
      k = rnd(Polity[h].ng); //do{ k = rnd(Polity[h].ng); }while( k == j);                          // select another group in polity
      if(Polity[h].g[k].lead[0].pi > ld->pi){                              // if leader in selected group has higher payoff copy strategy
#if UPDATE_LEAD_TAX
	gr->theta = ((Polity+h)->g+k)->theta;                              // copy tax on commoners
#endif
#if UPDATE_LEAD_EFFORT
	ld->x[0] = ((Polity+h)->g+k)->lead->x[0];                          // copy coordination effort
#endif
      }
    }
  }
  else if(v == 3){                                                       // myopic optimization
    double f, np, s;    
#if UPDATE_LEAD_TAX
    t = normal(gr->theta, Sigma_t);                                      // change tax level by normal distribution with deviation Sigma_t
    t = MAX(t, 0.00);
    t = MIN(t, 1.0);
#else
    t = gr->theta;
#endif
#if UPDATE_LEAD_EFFORT
    y = normal(ld->x[0], Sigma);
    y = MAX(y, 0.0000);
    //y = MIN((1.0-eta)*(B0*y*gr->X)*t/C1, y);
#else
    y = ld->x[0];
#endif
#if GP==3
    int i;
    double i_y = 1.0/y;
    for(s = 0, i = 0; i < gr->nc; i++){
      s += pow(gr->com[i].x[0], i_y);
    }
#else
    s = gr->X;    
#endif
    np = groupProduction(s, y, gr->nc);                                                     // new group production
    f = (1.0-eta)*t*gr->nc*B0*np - C1*y - C3*t;                                            // new payoff
    //double pf, op;
    //op = groupProduction(s, ld->x[0], gr->nc);
    //pf = (1.0-eta)*gr->theta*gr->nc*B0*op - C1*ld->x[0] - C3*gr->theta;
    
#if TestLeader
#if DEBUG
    // printing expected and old payoffs for particular leader
    if(h == TestH && j == TestG)
      printf("Leader H:%d g:%d \t expected payoff: %.3lf \told (real) payoff: %.3lf \tnew theta: %.3lf   old theta: %.3lf   new y: %.3lf   old y: %.3lf \n", h, j, f, ld->pi, t, gr->theta, y, ld->x[0]);
#endif
#endif
    if(f > ld->pi){                                                       // if new payoff is higher than current payoff, then change to new strategies    
      gr->theta = t;
      ld->x[0] = y;   
#if TestLeader
#if DEBUG      
      if(h == TestH && j == TestG)
	printf("new strategy chosen theta:%.6lf y: %.6lf expected payoff: %.6lf\n\n", gr->theta, ld->x[0], f);
#endif
#endif
    }
    
  }
}

// change strategy of chief in a polity
void changeStratChief(int v, int h, double sps)
/*
 * v: method of change
 * h: polity
 * sps: sum of strength of polities in system
 */
{
  int k;
  double eta, z;                                                           // stores mutated tax and effort
  individual *chf;
  polity *p;
  p = Polity+h;    
  chf = p->chief;
  if(v == 1){                                                            // random mutation
#if UPDATE_CHIEF_TAX
    p->eta = normal(p->eta, Sigma_t);                                    // change tax level by normal distribution with deviation Sigma_t
    p->eta = MAX(p->eta, 0.0);
    p->eta = MIN(p->eta, 1.0);
#endif
#if UPDATE_CHIEF_EFFORT
    chf->x[0] = normal(chf->x[0], Sigma);
    chf->x[0] = MAX(chf->x[0], 0.0001);    
    chf->x[0] = MIN(B2/C2, chf->x[0]);  
    // chf->x[0] = MIN(chf->x[0], 2);
#endif
  }
  else if(v == 2){                                                       // selective copying
    if(H > 1){
      do{ k = rnd(H); }while( k == h);                                     // select another polity in system
      if(((Polity+k)->chief)->pi > chf->pi){                                 // if chief in selected polity has higher payoff copy strategy
#if UPDATE_CHIEF_TAX
	p->eta = (Polity+k)->eta;                                          // copy tax on leaders
#endif 
#if UPDATE_CHIEF_EFFORT
	chf->x[0] = (Polity+k)->chief->x[0];                               // copy coordination effort
#endif
      }
    }
  }
  else if(v == 3){                                                       // myopic optimization
    double f; 
    double sp, st;                                                       // stores sum of group production and strength
#if UPDATE_CHIEF_TAX
    eta = normal(p->eta, Sigma_t);                                       // generate candidate strategy tax
    eta = MAX(eta, 0.0);
    eta = MIN(eta, 1.0);
#else
    eta = p->eta;                                                         // no change in tax
#endif
#if UPDATE_CHIEF_EFFORT
    z = normal(chf->x[0], Sigma);                                         // generate candidate strategy effort
    z = MAX(z, 0.0001);  
    z = MIN(z, B2/C2);
#else
    z = chf->x[0];                                                         // no change
#endif
#if PS==3
    double i_z = 1.0/z;
    for( sp = 0.0, k = 0; k < p->ng; k++){                      // through all groups in a polity
      sp += pow((p->g+k)->theta * (p->g+k)->p * (p->g+k)->nc, i_z);                      // sum of group production and theta
    }  
    sp = pow(eta*B0, i_z)*sp;                                               // new group production reaching to chief
#else
    for( sp = 0.0, k = 0; k < p->ng; k++){                      // through all groups in a polity
      sp += (p->g+k)->theta * (p->g+k)->nc * (p->g+k)->p;        // sum of product of group production, theta and no. of commoners in group
    }  
    sp = eta*sp*B0;                                               // new group production reaching to chief
#endif
    st = polityStrength(sp, z, p->ng);                                             // new polity strength  
    
    sps = sps - p->S + st;                                            // new sum of polity strength
    if(sps == 0.0)                                                // if sum of polity of strength is 0, then share the benefit
      f = B2 - C2*z - C4*eta;
    else
      f = B2*H*st/sps - C2*z - C4*eta;                                         // new payoff
    if(f > chf->pi){                                              // if new payoff is higher than current payoff, then change to new strategies
      //printf("%d ", ++count);
      p->eta = eta;
      chf->x[0] = z;
    }
  }
}


// updates strategy of all chiefs, leaders and commoners
void updateStrategy(double s)
/*
 * s: sum of strength of polities
 */
{
  int h, i, j, k;
  polity *p;
  group *g;
  double eta, theta, y;
  for(h = 0; h < H; h++){
    p = Polity+h;
    eta = p->eta;                                                     // tax on leaders before update
    // update strategy of chief of polity
#if UPDATE_CHIEF_EFFORT || UPDATE_CHIEF_TAX
    k = drand(Mdist);                                                 // randomly select method of changing strategy
    if(k){                                                            // change only if k is not 0 
      changeStratChief(k, h, s);                                      // change strategy of chief of polity p[h]
    }
#endif
    for(j = 0; j < p->ng; j++){                                       // through all groups
      g = p->g+j;
      theta = g->theta;                                               // tax on commoner's before update
      y = g->lead->x[0];                                              // effort by leader before update
      // update strategy of leader of group*/
#if UPDATE_LEAD_EFFORT || UPDATE_LEAD_TAX
      k = drand(Mdist);                                               // randomly select method of changing strategy
      if(k){                                                          // change only if k is not 0 
	changeStratLeader(k, h, j, eta);                              // change strategy of leader in group g[j] of polity p[h]
      }
#endif
#if UPDATE_COM
      for(i = 0; i < g->nc; i++){                                     // through all commoners in group
	// update strategy of commoner 
	k = drand(Mdist);                                             // choose type of update method
	if(k){
	  changeStratCommoner(k, h, j, i, theta, y);
	}
      }     
#endif
    }    
  }
}

void playGame()
{
  int h, j;
  double s;  
  for(h = 0; h < H; h++){                                          // through all polities
    for(j = 0; j < (Polity+h)->ng; j++){                           // through all groups
      // us vs nature within groups
      updateGroupProduction(h, j);                                 // calculates and update group production p
      updatePayoffAfterProduction(h, j);                           // updates payoff of each commoners and leader after production             
    }    
  }  
  // us vs them game between polities
  s = updatePolityStrengthAndChiefPayoff();                        // updates polity strength and chief's payoff of all polities; gets sum of polities strengths   
  // update strategies /  efforts
  updateStrategy(s);  
}

/*****************************************************************************************************************************************/
/*****************************************************************************************************************************************/


void init()                                                          // initialize values of global variables
{
  int h, i, j;
  individual *cm;
  group *gr;
  for(h = 0; h < H; h++){
    (Polity+h)->ng = G;                                               // set no. of groups in a polity
    (Polity+h)->eta = Eta;                                            // set tax on leader
    (Polity+h)->S = 1;                                                // set polity strength to 1
    cm = (Polity+h)->chief+0;
    cm->pi = 1;                                                       // reset payoff of chief
    cm->t = 0;                                                        // set tax paid by chief to 0. chiefs don't pay tax at all for now
    cm->x[0] = INIT_CHIEF_EFFORT;                                             // set coordination efforts by chief
    for(j = 0; j < G; j++){
      gr = (Polity+h)->g+j;
      gr->nc = N;                                                     // no. of commoners in a group set to N
      gr->nl = 1;                                                     // no. of leaders in a group set to 1
      gr->p = 0;                                                      // reset group production
      gr->X = 0;                                                      // reset average group effort
	gr->theta = Theta;                                                // tax on commoners set
      cm = gr->lead+0;
      cm->t = 0;                                                      // reset tax paid by leader
      cm->pi = 1;                                                     // reset payoff of leader
      cm->x[0] = INIT_LEAD_EFFORT;                                          // initialize coordiation effort by leader
      for(i = 0; i < N; i++){
	cm = gr->com+i;
	cm->x[0] = INIT_COM_EFFORT;                                   // initialize effort
	cm->pi = 1;                                                   // reset payoff
	cm->t = 0;                                                    // reset tax
      }      
    }
  }
  // initialize mutation probability distribution
  Mdist->p[0] = 1-V1-V2-V3;
  Mdist->p[1] = V1; 
  Mdist->p[2] = V2;
  Mdist->p[3] = V3;
  initdist(Mdist, 1);
}


int main(int argc, char **argv)
{
#if DEBUG
  feenableexcept(FE_DIVBYZERO| FE_INVALID|FE_OVERFLOW); // enable exceptions
#endif
  if(argc ^ 2){
    printf("Usage: ./ml ml.config\n");
    exit(1);
  }
  if(read_config(argv[1])){           // read config
    printf("READDATA: Can't process %s \n", argv[1]);
    return 1;
  }  
  
  if(TestH >= H){
    printf("\n Program terminated! Error: TestH (test polity index) must be less than H. Note: [index starts from 0] \n");
    exit(1);
  }
  if(TestG >= G){
    printf("\n Program terminated! Error: TestG (test group index) must be less than G. Note: [index starts from 0] \n");
    exit(1);
  }
  
  //FILE *gp;
  initrand(Seed);
  allocStatVar();                                             // allocate memory for global statistic variables
  int r, i, h, j, k;
  double sf, su, tbp;
  group *gr;  
  unsigned long seed;  
  char xdata[200], str[30];
  time_t now;     
  // print headers for values to be displayed on std output
#if !CLUSTER
  printf("\nValues:\t  x\t  y\t  z\t  theta\t  eta\t  pi_0\t  pi_1\t  pi_2\t  P\t  S\t  Q\t  U\t  seed\n");
#endif
  
  // set B0 and B2
  B0 = K*N; 
  B1 = 1;
  
  E = e;   
  X0 = N*x_0;
  
  for( r = 0; r < Runs; r++){                                 // through all sets of runs    
    // initialize pseudorandom generator with seed
    if(Seed == 0){      
      now = time(0);
      seed = ((unsigned long)now) + r;  
    }
    else{
      seed = Seed;
    }
    //seed = 1461867884;
    Seed_i = seed;
    initrand(seed);      
    fgr = fopen("grtest.dat", "w");
    fprintf(fgr, "T x_avg y z theta eta cm_av lead chief P Q U\n");  
    //printf("\n run# %d seed: %lu\n",r+1, seed);
    // setup and initialize variables
    setup();                                                  // allocates memory for all polity system variables
    init();                                                   // intialize polity system state values
    calcStat(0, r);
    for( i = 0; i < T; i++){                                  // through all time points of simulation
      playGame();                                             // play us vs nature and us vs them game and update strategies
      if( (i+1) % SKIP == 0){                                 // every SKIP time, take snapshot of states of traits	
	calcStat((i+1)/SKIP, r);                              // calculate statistics and write individual runs data to file
      }   
      
#if DEBUG      
      
      h=TestH;
      for(su = 0, j = 0; j < (Polity+h)->ng; j++){                           // through all groups	  	  
	gr = (Polity+h)->g+j;	  
	tbp = gr->theta*gr->p;                // theta*B0*P
	su += tbp*Polity[h].eta;                 // sum of U's
	if(j == TestG){
	  for(sf = 0, k = 0; k < gr->nc; k++){
	    sf += gr->com[k].pi;
	  }
	  fprintf(fgr, "%d %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf %.4lf ", i, gr->X, gr->lead[0].x[0], (Polity[h].chief)->x[0], gr->theta, Polity[h].eta, sf/gr->nc, gr->lead[0].pi, Polity[h].chief[0].pi, gr->p, (1-Polity[h].eta)*tbp );
	}	   	    
      }  
      fprintf(fgr, "%.4lf\n", su/Polity[h].ng);
	
      
#endif
    }

#if !CLUSTER
    calcStat(-1, -1);                                        // free file pointers for individual run data files
#if !AVERAGE_GRAPH_ONLY
    if(Runs > 1)
      plotallIndividualRun(r, 0); 
#if !TURNOFF_TEST_PLOT
      plotTraitsDist(r);
#endif
#endif
    
#if GRAPHS
    if(Runs > 1)
      plotallIndividualRun(r, 1);
#endif
    // remove individual datafiles if ALLDATAFILE set to 0
#if !ALLDATAFILE        
      sprintf(str, "x%d.dat", r); prep_file(xdata, str); remove(xdata);
      sprintf(str, "p%d.dat", r); prep_file(xdata, str); remove(xdata);
      sprintf(str, "t%d.dat", r); prep_file(xdata, str); remove(xdata);
      sprintf(str, "gp%d.dat", r); prep_file(xdata, str); remove(xdata);
      sprintf(str, "s%d.dat", r); prep_file(xdata, str); remove(xdata);       
#endif
#endif
    
    cleanup();                                                // free all memory allocated for polity system
    fclose(fgr);
#if !CLUSTER
#if !AVERAGE_GRAPH_ONLY
#if !TURNOFF_TEST_PLOT
    //gp = popen("gnuplot -persistent", "w");
   // plotLines(gp, "grtest.dat", 8, "traits evolution");
    //fflush(gp);
    //pclose(gp);
    plotTest("grtest.dat", 8, "traits evolution");
#endif
#endif
#endif
  }
  writeDataToFile();  
  clearStatVar();                                             // free other memory allocated for statistics variables
#if !CLUSTER
  plotall(0);
#if GRAPHS
  plotall(1);
#endif
#endif  
  
  //printf("\n\nB0: %lf, B1: %lf, B2: %lf \n", B0, B1, B2);
  return 1;
}

/* 
gcc -Wall -O2 -march=native -pipe -o ml ml.c -lm
./ml ml.config
valgrind -v --track-origins=yes --leak-check=full --show-leak-kinds=all ./ml ml.config
gcc -g -o ml ml.c -lm
gdb ml
run ml.config

//profiling code
gcc -Wall -O2 -march=native -pipe -pg ml.c -o ml -lm
./ml ml.config
 gprof ml gmon.out > analysis.txt  
*/

