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


#define TestH 0                  // test polity interested in
#define TestG 4                  // test group interested in

#define CLUSTER 0            // is this simulation on cluster yes or no (1 or 0)
#define SKIP 1           // time interval between snapshot of states
#define STU 500             // summary time period for simulation

#define AVERAGE_GRAPH_ONLY 1  // if 1, generate only average graphs and supress individual run graphs
#define TURNOFF_TEST_PLOT 0 // if 1, does not produce frequency distribution and traits graphs for test group

#define ALLDATAFILE 1        // if 1, generate all data files for individual runs and summary too 
#define GRAPHS      0        // if 1, saves graphs as png files 

#define INIT_COM_EFFORT rnd(2)              // 0 or 1
#define INIT_LEAD_EFFORT 0
#define INIT_CHIEF_EFFORT 0

#define INIT_CHIEF_PUN_EFFORT 0
#define INIT_LEAD_PUN_EFFORT U01()
// change values to 0 if update strategy for any one of commoner or lead or chief is to be turned off
#define UPDATE_COM 1
#define UPDATE_LEAD_EFFORT 0
#define UPDATE_CHIEF_EFFORT 0

#define UPDATE_CHIEF_PUN_EFFORT 0
#define UPDATE_LEAD_PUN_EFFORT 1

#define TURNOFF_RUNS_DATA 1            // does not store individual runs dynamics data files in cluster

// strategy update options and values
// option_set = {random_mutation, selective_copy, optimization}
// strategy update option set for each level = {option1_set, option2_set, option_3_set, ...}
double Vc[][3] = { {0.01, 0.24, 0.0}, {0.00, 0.00, 0.25} }; // commoner strategy update method probability sets for different options (1, 2, 3)
double Vl[][3] = { {0.01, 0.24, 0.0}, {0.01, 0.24, 0.0} };  // leader strategy update method probability sets for different options (1, 2, 3)
double Vcf[][3] = { {0.01, 0.24, 0.0}, {0.01, 0.24, 0.0} };  // chief strategy update method probability sets for different options (1, 2, 3)

#define GP 1 // 1, 2, 3; uses different calculation for group production
#define PS 1 // 1, 2, 3; uses different calculation for polity strength

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
  double        pi;           // payoff   
  unsigned int  x;            // 1 or 0 
  double        k;            // punishment incurred from leader
} Commoner;

typedef struct
{
  double        pi;           // payoff 
  double        y;            // coordination effort
  double        p;            // punishment effort to punish commoners
  double        K;            // punishment incurred from chief
  double        d_l;          // small delta 
  
} Leader;

typedef struct
{
  double        pi;           // payoff 
  double        z;            // coordination effort
  double        q;            // punishment effort    
  double        d_c;          // capital delta
} Chief;

typedef struct
{
  double        P;            // group production
  double        X;            // sum of efforts from commoners in group
  double        t_u;          // theta_u: tax on commoners in a group 
  double        t_d;          // theta_d: rewards to commoner
  unsigned int  nc;           // number of commoners in a group
  Leader        *lead;        // leader
  Commoner      *com;         // group members / commoners; 
} group;

typedef struct
{
  unsigned int ng;            // no. of groups in a polity   
  double       e_u;           // tax  on leaders
  double       e_d;           // rewards to leaders
  double       SY;            // sum of production from all groups in a polity
  double       Q;             // 
  group        *g;            // polity members  
  Chief        *chief;        // chief;
  
} polity;

polity *Polity;               // "system structure"

// Global variables
/** configuration parameters **/
unsigned Seed;                    // Seed: Seed of psuedorandom number passed
unsigned long Seed_i;             // Seed_i: seed of each run
unsigned int Runs;                         // no. of runs
unsigned int G;                            // no. of groups
unsigned int n;                            // no. of commoners
unsigned int H;                            // no. of politites
unsigned int T;                            // max time of simulation
double Theta_u, Theta_d, Eta_u, Eta_d;     // taxes / rewards
double cx, cy, cz;                // cost parameter
double b, B;                      // expected benefit for a game
double Sigma;                     // standard deviation for distribution of mutation for contribution
double Sigma_t;                   // standard deviation for distribution of mutation for tax levels
double x_0, y_0, z_0, Y0, X0;     // x_0: half effort of commoner; y_0: half effort of leader; z_0: half effort of chief; Y0: half strength of polity
double e, E;                      // e: leader efficiency; E: chief efficiency
double k, K, delta, DELTA;        // punishments parameter
unsigned int L;                         // no. of least productive leaders to be punished by a chief 
double m;                         // probability of migration
unsigned int Vop;                       // strategy update set option (1, 2, 3)
double Rho;                       // probability that polity would go to us vs them game 
unsigned int *UvT;                // keep list of index if polities participate in usvsthem game or not

// statistical variables
double *xmean;          // average effort by individuals
double *ymean;            // average effort by leaders
double *zmean;            // average effort by chiefs
double *pi0mean;                  // average payoff of commoners
double *pi1mean;                  // average payoff of leaders
double *pi2mean;                  // average payoff of chiefs
double *pmean;                    // average punishment effort of leaders
double *qmean;                    // average punishment effort of chiefs
double *TUmean;                    // average tax on commoners
double *EUmean;                    // average tax on leaders
double *TDmean;                    // average rewards to commoners
double *EDmean;                    // average rewards to leaders
double *Pmean;                    // group production
double *Qmean;                    // production to leader per each commoner

dist *Vcdist, *Vldist, *Vcfdist;                      // 'dist' is defined in rand.c; Probability distribution of choosing method of changing efforts and tax levels 

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
  EXPECT(1, "int      H         = %d;", &H);
  
  EXPECT(1, "double   b         = %lf;", &b);  
  EXPECT(1, "double   B         = %lf;", &B);   
  EXPECT(1, "double   cx        = %lf;", &cx);
  EXPECT(1, "double   cy        = %lf;", &cy);
  EXPECT(1, "double   cz        = %lf;", &cz);     
  
  EXPECT(1, "unsigned L         = %u;", &L);
  EXPECT(1, "double   k         = %lf;", &k);
  EXPECT(1, "double   K         = %lf;", &K);
  EXPECT(1, "double   delta     = %lf;", &delta);
  EXPECT(1, "double   DELTA     = %lf;", &DELTA);
  
  EXPECT(1, "double   Theta_u   = %lf;", &Theta_u);
  EXPECT(1, "double   Theta_d   = %lf;", &Theta_d);
  EXPECT(1, "double   Eta_u     = %lf;", &Eta_u);
  EXPECT(1, "double   Eta_d     = %lf;", &Eta_d);  

  EXPECT(1, "unsigned Vop       = %d;", &Vop);
  EXPECT(1, "double   m         = %lf;", &m); 
  EXPECT(1, "double   Rho       = %lf;", &Rho)
  
  EXPECT(1, "double   x0        = %lf;", &x_0);
  EXPECT(1, "double   y0        = %lf;", &y_0);
  EXPECT(1, "double   z0        = %lf;", &z_0);
  EXPECT(1, "double   X0        = %lf;", &X0);
  EXPECT(1, "double   Y0        = %lf;", &Y0);  
  
  EXPECT(1, "double   e         = %lf;", &e);
  EXPECT(1, "double   E         = %lf;", &E);
  
  EXPECT(1, "double   Sigma     = %lf;", &Sigma);
  EXPECT(1, "double   Sigma_t   = %lf;", &Sigma_t);    

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
  sprintf(fname, "g%02dn%02db%0.2fB%.2fk%0.2fdl%.2ftu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fVop%dRho%.2fY0%.2fX0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s", G, n, b, B, k, delta, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, Vop, Rho, Y0, X0, e, E, x_0, y_0, z_0, apndStr);
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

// allocate memory for Polity system
void setup()
{
  int h, j;
  group *gr;
  // allocate memory for Polity and groups and commoners
  Polity = malloc(H*sizeof(polity));
  for(h = 0; h < H; h++){                                            // through all polities
    (Polity+h)->chief = malloc(sizeof(Chief));                 // allocate memory for chief
    (Polity+h)->g = malloc(G*sizeof(group));                          // allocate memory for groups in each polity    
    for(j = 0; j < G; j++){                                          // through all groups in a polity
      gr = (Polity+h)->g+j;                                             // get reference to group pointer
      gr->lead = malloc(sizeof(Leader));                       // allocate memory for leader in a group
      gr->com = malloc(n*sizeof(Commoner));                          // allocate memory for commoners in a group
    }
  }
  Vcdist = allocdist(4);                                              // allocate memory for distribution for update strategy method
  Vldist = allocdist(4);                                              // allocate memory for distribution for update strategy method
  Vcfdist = allocdist(4);                                              // allocate memory for distribution for update strategy method
  UvT = malloc(H*sizeof(unsigned int));  
}

void allocStatVar()
{
  // allocate memory to store statistics points of snapshot along time period of simulaltion
  xmean = calloc((int)(T/SKIP+1), sizeof(double));         // for effort by commoners
  ymean = calloc((int)(T/SKIP+1), sizeof(double));         // for effort by leaders
  zmean = calloc((int)(T/SKIP+1), sizeof(double));         // for effort by chiefs
  pi0mean = calloc((int)(T/SKIP+1), sizeof(double));                   // for payoff of commoners
  pi1mean = calloc((int)(T/SKIP+1), sizeof(double));                   // for payoff of leaders
  pi2mean = calloc((int)(T/SKIP+1), sizeof(double));                   // for payoff of chiefs
  pmean = calloc((int)(T/SKIP+1), sizeof(double));                     // for punishment effort by leaders 
  qmean = calloc((int)(T/SKIP+1), sizeof(double));                     // for punishment effort by chiefs
  TUmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for tax paid by commoners to leaders
  EUmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for tax paid by leaders to chiefs
  TDmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for tax paid by commoners to leaders
  EDmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for tax paid by leaders to chiefs
  Pmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for group production 
  Qmean = calloc((int)(T/SKIP+1), sizeof(double));                    // for polity production // singling out B from BQ    
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
  freedist(Vcdist); freedist(Vldist); freedist(Vcfdist); 
  free(UvT);
}

void clearStatVar()
{
  free(xmean); free(ymean); free(zmean); free(pi0mean); free(pi1mean); free(pi2mean); free(pmean); free(qmean);
  free(TUmean); free(EUmean); free(TDmean); free(EDmean);
  free(Pmean); free(Qmean);   
}

/***************************************Statistics calcualtions and writing to file*******************************************************/
/*****************************************************************************************************************************************/

// calculates all the stats
void calcStat(int k, int r) 
{
  double xm = 0, ym = 0, zm = 0, p0m = 0, p1m = 0, p2m = 0, tum = 0, tdm = 0, eum = 0, edm = 0, pm = 0, qm = 0, Pm = 0, Qm =  0;
  int h, j, i, ncm, nl, ncf;  
  polity *p;
  group *g;
  Commoner *com;
  Leader *ld;
  Chief *chief;
// #if !CLUSTER
  //FILE *fx, *fy, *fz, *fp0, *fp1, *fp2, *ft, *fe;
//   char fname[100];
//   
//   if(k== (T/SKIP)){                                // data files for freq. distribution
//     sprintf(fname, "x_%d.dat", r);
//     fx = fopen(fname, "w");   
//     sprintf(fname, "y_%d.dat", r);
//     fy = fopen(fname, "w");
//     sprintf(fname, "z_%d.dat", r);
//     fz = fopen(fname, "w");
//     sprintf(fname, "t_%d.dat", r);
//     ft = fopen(fname, "w");
//     sprintf(fname, "e_%d.dat", r);
//     fe = fopen(fname, "w");
//     sprintf(fname, "p0_%d.dat", r);
//     fp0 = fopen(fname, "w");
//     sprintf(fname, "p1_%d.dat", r);
//     fp1 = fopen(fname, "w");
//     sprintf(fname, "p2_%d.dat", r);
//     fp2 = fopen(fname, "w");
//     
//   }
// #endif

#if (!TURNOFF_RUNS_DATA || !CLUSTER )
  static FILE **fp = NULL;  // file pointers for individual run  
  if(k < 0){
    for (h = 5; h--; fclose(fp[h]));
    free(fp); 
    return;
  }     
  
  if(!k){   
    char ixdata[200], ipdata[200], itdata[200], igpdata[200], ipundata[200], tstr[100];
    fp = malloc(5*sizeof(FILE *));    
    sprintf(tstr, "x%d.dat", r); prep_file(ixdata, tstr);    
    sprintf(tstr, "p%d.dat", r); prep_file(ipdata, tstr);    
    sprintf(tstr, "t%d.dat", r); prep_file(itdata, tstr);    
    sprintf(tstr, "gp%d.dat", r); prep_file(igpdata, tstr); 
    sprintf(tstr, "pun%d.dat", r); prep_file(ipundata, tstr);
    fp[0] = fopen(ixdata, "w");
    fp[1] = fopen(ipdata, "w");
    fp[2] = fopen(itdata, "w");  
    fp[3] = fopen(igpdata, "w");
    fp[4] = fopen(ipundata, "w");
    
    // write headers       
    fprintf(fp[0], "0 \t com\tlead\tchief\n");
    fprintf(fp[1], "0 \t com\tlead\tchief\n");
    fprintf(fp[2], "0 \t theta_u \t theta_d \t eta_u \t eta_d\n");
    fprintf(fp[3], "0\tP\tQ\n");
    fprintf(fp[4], "0\tp\tq\n");
  }
#endif
  // calculate stat of traits
  for(h = 0; h < H; h++){                                                               // through all polities
    p = Polity+h;
    for(j = 0; j < p->ng; j++){                                                         // through all groups
      g = p->g+j;
      for(i = 0; i < g->nc; i++){                                                       // through all commoners
	com = g->com+i;
	xm += com->x;                                                                // sum of commoner's effort
	p0m += com->pi;                                                                 // sum of commoner's payoff
// #if !CLUSTER
// 	if(k== (T/SKIP)){
// 	  fprintf(fx, "%.6lf\n", (double)com->x);
// 	  fprintf(fp0, "%.6lf\n", com->pi);
// 	}
// #endif
      }      
      // leader stats
      ld = g->lead;                                                                // ind == leader
      ym += ld->y;                                                                // sum of leader's efforts
      p1m += ld->pi;                                                                 // sum of leader's payoff
      pm += ld->p;                                                                   // sum of leader's punishment effort
      tum += g->t_u;                                                                   // sum of tax factor theta_u on commoners
      tdm += g->t_d;                                                                   // sum of reward factor theta_d on commoners
/*#if !CLUSTER
      if(k== (T/SKIP)){
	fprintf(fy, "%.6lf\n", ld->y);
	fprintf(fp1, "%.6lf\n", ld->pi);     
	fprintf(ft, "%.6lf \t %.6lf\n", g->t_u, g->t_d);
      }
#endif   */  
      Pm += g->P;                                                                     // sum P over all groups        
    }
    // chief's stat
    chief = p->chief+0;
    zm += chief->z;                                                               // sum of chief's efforts
    p2m += chief->pi;
    qm += chief->q;                                                                         // sum of chief's punishment effort in each polity
    eum += p->e_u;                                                                       // sum of tax on leaders by chiefs in polity
    edm += p->e_d;                                                                       // sum of tax on leaders by chiefs in polity
    Qm += p->Q;    
    
// #if !CLUSTER
//     if(k== (T/SKIP)){
//       fprintf(fz, "%.6lf\n", (chief)->z);
//       fprintf(fp2, "%.6lf\n", (chief)->pi);
//       fprintf(fe, "%.6lf \t  %.6lf\n", p->e_u, p->e_d);
//     }
// #endif
  }
  // compute averages of traits
  ncm = H*G*n;         // total no. of commoners in system
  nl  = H*G;           // total no. of leaders in system
  ncf = H;             // total no. of chiefs in system
  
  xm /= ncm; 
  p0m /= ncm;
  
  ym /= nl;
  p1m /= nl;
  tum /= nl;
  tdm /= nl;
  pm  /= nl;  
  
  zm /= ncf;
  p2m /= ncf;  
  qm /= ncf;  
  eum /= ncf;  
  edm /= ncf; 
  
  Pm /= nl;        
  Qm /= ncf;    
  
  // store in averages in variable as sum to average it after multiple runs
  xmean[k] += xm;
  ymean[k] += ym;
  zmean[k] += zm;
  pi0mean[k] += p0m;
  pi1mean[k] += p1m;
  pi2mean[k] += p2m;
  pmean[k] += pm;
  qmean[k] += qm;
  TUmean[k] += tum;
  TDmean[k] += tdm;
  EUmean[k] += eum;
  EDmean[k] += edm;
  Pmean[k] += Pm;
  Qmean[k] += Qm;  
#if (!TURNOFF_RUNS_DATA || !CLUSTER )
  // write data for individual runs
  fprintf(fp[0], "%d  %.4lf  %.4lf  %.4lf\n", k, xm, ym, zm);
  fprintf(fp[1], "%d  %.4lf  %.4lf  %.4lf\n", k, p0m, p1m, p2m);
  fprintf(fp[2], "%d  %.4lf  %.4lf  %.4lf  %.4lf\n", k, tum, tdm, eum, edm);  
  fprintf(fp[3], "%d  %.4lf  %.4lf\n", k, Pm, Qm); 
  fprintf(fp[4], "%d  %.4lf  %.4lf\n", k, pm, qm); 
#endif
#if !CLUSTER
  // print final values for each run
  if(k == T/SKIP){
    printf("run#%d \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %lu\n", r, xm, ym, zm, pm, qm, tum, tdm, eum, edm, p0m, p1m, p2m, Pm, Qm, Seed_i);  
//     fclose(fx); fclose(fy); fclose(fz); fclose(ft); fclose(fe); fclose(fp0); fclose(fp1); fclose(fp2);
  }
  
#endif
}

void writeDataToFile()
{
  int j;
  char xdata[200], pdata[200], tdata[200], gpdata[200], pundata[200];
  FILE **fp = malloc(10*sizeof(FILE *));    
  prep_file(xdata, "x.dat");    
  prep_file(pdata, "p.dat");    
  prep_file(tdata, "t.dat");   
  prep_file(gpdata, "gp.dat");
  prep_file(pundata, "pun.dat");
  fp[0]= fopen(xdata, "w");
  fp[1]= fopen(pdata, "w");
  fp[2]= fopen(tdata, "w");
  fp[3]= fopen(gpdata, "w");
  fp[4]= fopen(pundata, "w");
#if ALLDATAFILE
  int stu = (T-STU)/SKIP;
  double xsum = 0, ysum = 0, zsum = 0, p0sum = 0, p1sum = 0, p2sum =0, t_usum = 0, e_usum = 0, t_dsum = 0, e_dsum = 0, psum = 0, qsum = 0, Psum = 0;  
#endif  
  // write headers
  fprintf(fp[0], "0 \t com\tlead\tchief\n");
  fprintf(fp[1], "0 \t com\tlead\tchief\n");
  fprintf(fp[2], "0 \t theta_u \t theta_d \t eta_u \t eta_d \n");
  fprintf(fp[3], "0 \t P \t Q\n");
  fprintf(fp[4], "0 \t p \t q\n");
  
  for(j = 0; j < (int)(T/SKIP)+1; j++){
    // average accumulated mean values for multiple runs
    xmean[j] /= (double)Runs;
    ymean[j] /= (double)Runs;
    zmean[j] /= (double)Runs;
    pi0mean[j] /= (double)Runs;
    pi1mean[j] /= (double)Runs;
    pi2mean[j] /= (double)Runs;
    pmean[j] /= (double)Runs;
    qmean[j] /= (double)Runs;
    TUmean[j] /= (double)Runs;
    TDmean[j] /= (double)Runs;
    EUmean[j] /= (double)Runs;
    EDmean[j] /= (double)Runs;
    Pmean[j] /= (double)Runs;
    Qmean[j] /= (double)Runs;    
    // write data to file
    fprintf(fp[0], "%d  %.4lf  %.4lf  %.4lf\n", j, xmean[j], ymean[j], zmean[j]);
    fprintf(fp[1], "%d  %.4lf  %.4lf  %.4lf\n", j, pi0mean[j], pi1mean[j], pi2mean[j]);
    fprintf(fp[2], "%d  %.4lf  %.4lf  %.4lf  %.4lf\n", j, TUmean[j], TDmean[j], EUmean[j], EDmean[j]); 
    fprintf(fp[3], "%d  %.4lf  %.4lf\n", j, Pmean[j], Qmean[j]);
    fprintf(fp[4], "%d  %.4lf  %.4lf\n", j, pmean[j], qmean[j]);    
    
#if ALLDATAFILE
    if(j < stu) continue;
    xsum += xmean[j];
    ysum += ymean[j];
    zsum += zmean[j];
    p0sum += pi0mean[j];
    p1sum += pi1mean[j];
    p2sum += pi2mean[j];
    t_usum += TUmean[j];
    t_dsum += TDmean[j];
    e_usum += EUmean[j];
    e_dsum += EDmean[j];
    psum += pmean[j];
    qsum += qmean[j];
    Psum += Pmean[j];
#endif
  }
  for(j = 0; j < 5; j++){
    fclose(fp[j]);
  } 
#if ALLDATAFILE
  double st = STU/SKIP + 1;
  xsum /= st; ysum /= st; zsum /= st; p0sum /= st; p1sum /= st; p2sum /= st; 
  t_usum /= st; t_dsum /= st; e_usum /= st; e_dsum /= st;                // computing average over summary period time
  psum /= st; qsum /= st;  // todo: Psum /= st; Qsum /= st; 
  // write data to file
  char  xsumdata[200], psumdata[200], tsumdata[200], punsumdata[200], Psumdata[200];
  prep_file(xsumdata, "xsum.dat");    
  prep_file(psumdata, "psum.dat");    
  prep_file(tsumdata, "tsum.dat");  
  prep_file(punsumdata, "punsum.dat"); 
  prep_file(Psumdata, "Psum.dat");
  fp[5]= fopen(xsumdata, "w");
  fp[6]= fopen(psumdata, "w");
  fp[7]= fopen(tsumdata, "w");
  fp[8]= fopen(punsumdata, "w");
  fp[9]= fopen(Psumdata, "w");
  // write headers  
  fprintf(fp[5], "com\tlead\tchief\n");
  fprintf(fp[6], "com\tlead\tchief\n");
  fprintf(fp[7], "theta_u \t theta_d \t eta_u \t eta_d\n");
  fprintf(fp[8], "p \t q \n"); 
  fprintf(fp[9], "P \n");
     
  fprintf(fp[5], "%.4lf  %.4lf  %.4lf\n", xsum, ysum, zsum);
  fprintf(fp[6], "%.4lf  %.4lf  %.4lf\n", p0sum, p1sum, p2sum);
  fprintf(fp[7], "%.4lf  %.4lf  %.4lf  %.4lf\n", t_usum, t_dsum, e_usum, e_dsum); 
  fprintf(fp[8], "%.4lf  %.4lf\n", psum, qsum); 
  fprintf(fp[9], "%.4lf\n", Psum);
  for(j = 5; j < 10; j++){
    fclose(fp[j]);
  }
#endif  
  free(fp);
#if !CLUSTER
  // print averaged final values
  printf("\nAvg: \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf \t%.4lf  %.4lf  %.4lf  %.4lf  %.4lf  %.4lf\n", xmean[T/SKIP], ymean[T/SKIP], zmean[T/SKIP],  pmean[T/SKIP], qmean[T/SKIP], TUmean[T/SKIP], TDmean[T/SKIP], EDmean[T/SKIP], EDmean[T/SKIP], pi0mean[T/SKIP], pi1mean[T/SKIP], pi2mean[T/SKIP], Pmean[T/SKIP], Qmean[T/SKIP]);  
#endif
}

void plotDynamicsGraphs(FILE *gp, char *title, int datacolumn, char *xdata, char *pdata, char *pundata, char *gpdata)
{
  fprintf(gp, "stats '%s' using 2:3 prefix 'A' nooutput\n", xdata);
  fprintf(gp, "stats '%s' using 4 prefix 'B' nooutput\n", xdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'C' nooutput\n", pdata);
  fprintf(gp, "stats '%s' using 4 prefix 'D' nooutput\n", pdata);
  fprintf(gp, "stats '%s' using 2:3 prefix 'E' nooutput\n", gpdata);  
  fprintf(gp, "stats '%s' using 2:3 prefix 'F' nooutput\n", pundata);

  fprintf(gp, "set xlabel 'Time' \n");
  fprintf(gp, "set title ''\n");
  fprintf(gp, "set multiplot layout 4,1 title '%s' \n", title);    // set subplots layout  
  
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
  fprintf(gp, "unset multiplot \n");
}

void plotAverage(int m)
/*
 * m: 0 or 1; to save graphs as image file or not
 */
{ 
  char xdata[200], pdata[200], tdata[200], gpdata[200], pundata[200];
  char title[200], xpng[200], str[100];
  int datacolumn = 3+1;
  sprintf(str, "x.dat"); prep_file(xdata, str);
  sprintf(str, "p.dat"); prep_file(pdata, str);
  sprintf(str, "t.dat"); prep_file(tdata, str); 
  sprintf(str, "gp.dat"); prep_file(gpdata, str); 
  sprintf(str, "pun.dat"); prep_file(pundata, str); 
  
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "b:%.1f, B:%.1f t_u:%.1f, t_d:%.1f, e_u:%.1f, e_d%.1f, c:%.1f,%.1f,%.1f, Vop:%d, Rho:%.2f, e:%.1f, E:%.1f, x0:%.2f, y0:%.2f, z0:%.2f, Y0:%.2f", b, B, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, Vop, Rho, e, E, x_0, y_0, z_0, Y0);
  fprintf(gp, "set key outside vertical  spacing 1 width 1 height -2\n");        
  if(m){ // save graphs as file
    sprintf(xpng, "g%02dn%02db%0.2fB%.2ftu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fRho%.2fY0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2fVop%d.png", G, n, b, B, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, Rho, Y0, e, E, x_0, y_0, z_0, Vop); 
    fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,8\" \n");
    fprintf(gp, "set output '%s' \n", xpng);
  }
  else{
    fprintf(gp, "set term x11 dashed %d \n", 0);
  }
  fprintf(gp, "set lmargin at screen 0.1 \n");
  fprintf(gp, "set rmargin at screen 0.8 \n");
  fprintf(gp, "set label 1 'Average' at screen 0.12,0.99 center font \",13\"\n");  // label to indicate average  
  plotDynamicsGraphs(gp, title, datacolumn, xdata, pdata, pundata, gpdata);  
  fprintf(gp, "unset label 1\n");    
  fflush(gp); 
  pclose(gp);  
  
#if !CLUSTER
#if !ALLDATAFILE
  remove(xdata);
  remove(pdata);
  remove(tdata);
  remove(gpdata);
  remove(pundata);
#endif
#endif
}

void plotIndividualRun(int r, int m)
/*
 * r: run
 * m: 0 or 1; to save graphs as image file or not
 */
{
    // write data to file
  char xdata[200], pdata[200], tdata[200], gpdata[200], pundata[200];
  char title[200], xpng[200], str[100];
  int datacolumn = 3+1;
  sprintf(str, "x%d.dat", r); prep_file(xdata, str);
  sprintf(str, "p%d.dat", r); prep_file(pdata, str);
  sprintf(str, "t%d.dat", r); prep_file(tdata, str); 
  sprintf(str, "gp%d.dat", r); prep_file(gpdata, str); 
  sprintf(str, "pun%d.dat", r); prep_file(pundata, str); 
  
  FILE * gp = popen ("gnuplot -persistent", "w"); // open gnuplot in persistent mode
  sprintf(title, "b:%.1f, B:%.1f t_u:%.1f, t_d:%.1f, e_u:%.1f, e_d%.1f, c:%.1f,%.1f,%.1f, Vop:%d, Rho:%.2f, e:%.1f, E:%.1f, x0:%.2f, y0:%.2f, z0:%.2f, Y0:%.2f", b, B, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, Vop, Rho, e, E, x_0, y_0, z_0, Y0);

  fprintf(gp, "set key outside vertical height -1\n");        
  if(m){ // save graphs as file
    sprintf(xpng, "g%02dn%02db%0.2fB%.2ftu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fRho%.2fY0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2fVop%d_%d.png", G, n, b, B, Theta_u, Theta_d, Eta_u, Eta_d, cx, cy, cz, Rho, Y0, e, E, x_0, y_0, z_0, Vop, r); 
    fprintf(gp, "set term pngcairo size 1024,768 enhanced color solid font \"Helvetica,8\" \n");
    fprintf(gp, "set output '%s' \n", xpng);
  }
  else{
    fprintf(gp, "set term x11 dashed %d \n", 0);
  }
  fprintf(gp, "set lmargin at screen 0.1 \n");
  fprintf(gp, "set rmargin at screen 0.8 \n");  
  fprintf(gp, "set label 1 'run: %d' at screen 0.13,0.98 center font \",13\"\n", r);  // label to indicate individual run index  
  plotDynamicsGraphs(gp, title, datacolumn, xdata, pdata, pundata, gpdata);  
  fprintf(gp, "unset label 1\n");  
  fflush(gp); 
  pclose(gp); 

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
  double r = ( 1.0 + (e*y)/(n*y_0+y) );
  return sx/( sx + (X0)/r );
}

double polityStrength(double sY, double z, int g)
/*
 * sx: sum of group production reaching to chief
 * z: effort by chief
 * g: polity size /  no. of groups in a polity
 */
{  
  double R = ( 1.0 + (E*z)/(z + g*z_0) );
  return sY/( sY + (g*Y0)/R );
}

// updates group production p
void updateGroupProduction(int h, int j)
{
  int i;
  double sx;
  group *gr;
  gr = Polity[h].g+j;
  for(sx = 0, i = 0; i < gr->nc; i++){                       // through all commoners
    sx += gr->com[i].x;                                   // sum production efforts from commoners
  }
  gr->X = sx;                                                // sum of efforts from commoners      
  gr->P = groupProduction(sx, gr->lead->y, gr->nc);            // group production 
}

// updates Q value in a polity
void updatePolityStrength(int h)
{
  int i;
  double sY;
  polity *p = Polity+h;  
  group *g;
  for(sY = 0, i = 0; i < p->ng; i++){                    // through all groups
    g = p->g+i;
    sY += g->t_u*g->P;                                   // sum of theta_u*P
  }
  sY *= p->e_u * n*b;                                      // sum of Y ; Y = eta_u*theta_u*n*b*P
  p->Q = polityStrength(sY, p->chief->z, p->ng);         // updates Q value
}

// returns payoff of leader before punishment
double pi_leader(double y, double eta_u, double theta_u, double P, double ppl)
/*
 * ppl: polity production share to each group leader
 */
{
  return (1.0-eta_u)*theta_u*n*b*P - cy*y + ppl;
}

// returns payoff of chief before punishment
double pi_chief(double z, double eta_d, double BQ, int uvt, int h1, double sQ)
/*
 * z : effort by chief
 * eta_d: downward reward tax by chief
 * BQ: benefit B * polity strength Q
 * uvt: us vs them game or not
 * h1: no. of polities in us vs them game pool
 * sQ: sum of polity strength in us vs them game pool
 */
{
  if(uvt == 1){    // us vs them    
    if(sQ > 0.0)
      return eta_d*n*G*h1*BQ/sQ - cz*z;
    else
      return eta_d*n*G*B - cz*z;      
  }
  else{            // us vs nature
    return eta_d*n*G*BQ - cz*z;
  }
}

// updates payoff of a whole polity
void updatePayoffBeforePunishment(int h, int h1, double sQ)
/*
 * h: polity index
 * h1: no. of polities in us vs them game pool for this time step
 * sQ: sum of polity strength of polities in us vs them game pool
 */
{
  int i, j;
  double ppg, ppc, gpc;
  polity *p = Polity+h;  
  group *g;
  Leader *ld;
  Commoner *com;
  double BQ = B*p->Q;
  ppg = (1.0-p->e_d)*n*BQ;                                 // polity production share to each group from chief 
  for(j = 0; j < p->ng; j++){
    g = p->g+j;
    ppc = (1.0 - g->t_d)*ppg/(g->nc);                              // polity production share to each commoner
    gpc = 1.0 + (1.0 - g->t_u)*b*g->P;                       // group production share to each commoner + 1
    // update payoff of commoners before punishment
    for(i = 0; i < g->nc; i++){
      com = g->com+i;
      com->pi = gpc - cx*com->x + ppc;                               // payoff due to group production - cost of x + polity production share to commoner
    }
    // update payoff of leader before punishment
    ld = g->lead;
    ld->pi = pi_leader(ld->y, p->e_u, g->t_u, g->P, ppg*g->t_d); // payoff due to group production - cost of y + polity production share to leader
  }
  p->chief->pi = pi_chief(p->chief->z, p->e_d, BQ, UvT[h], h1, sQ);     // payoff due to polity production - cost of z;
}

// returns index array 'l' of leader/ group sorted according to efforts made by leaders
int* sortedLeaderIndexArrayWithEffort(polity *p)
{
  int j;
  int *l = malloc(p->ng*sizeof(int));                               // leader/group index array
  double *y = malloc(p->ng*sizeof(double));                         // leader effort array
  for(j = 0; j < p->ng; j++){
    l[j] = j;
    y[j] = (p->g+j)->lead->y;
  }
  merge_sort_key_value(l, y, p->ng);                                // sort index array k in ascending according to values in effort array y  
  free(y);  
  return l;
}

// punish commoners and leaders not making enough efforts in a polity h
void punish(int h)
{
  int i, j, *l;
  polity *p = Polity+h;
  group *g;
  Leader *ld;
  Commoner *com;
  //punish commoners
  for(j = 0; j < p->ng; j++){
    g = p->g+j;
    ld = g->lead;
    for(i = 0; i < g->nc; i++){
      com = g->com+i;
      if(com->x == 0){                                              // if no effort, punish with probability p
	if(U01() < ld->p){
	  com->pi -= k;                                             // take off k from commoner's payoff
	  ld->pi -= delta;                                           // punishment is costly, so take off delta from leader's payoff
	}
      }
    }
  }
  //punish leaders
  //l = malloc(p->ng*sizeof(int));;
  l = sortedLeaderIndexArrayWithEffort(p);                          // sorted index array of leaders according to efforts made  
  for(j = 0; j < L; j++){ // make sure size of polity is larger than L here, if polity size varies         
    if(U01() < p->chief->q){                                          // punish with probability q
      ld = (p->g+l[j])->lead;
      ld->pi -= K;                                                    // take off K from leader's payoff
      p->chief->pi -= DELTA;                                          // punishment is costly, so take off DELTA from chief's payoff
    }
  }    
  free(l);
}

// returns new updated effort for a commoner with optimization
int updateStratCommoner_v3(int x, double X, double y, double pn, int n, double t_u, double ppc)
{  
    double p0, p1;                                       
    p0 = 1.0 + (1.0-t_u)*b*groupProduction(X-x+0, y, n) + ppc - k*pn;    // new payoff for 0 effort
    p1 = 1.0 + (1.0-t_u)*b*groupProduction(X-x+1, y, n) + ppc - cx - k*pn;    // new payoff for 0 effort    
    //return (U01() < 1.0/( 1.0+exp(p0 - p1) ) )? 1 : 0;
    return (p1 > p0)? 1: 0;  
}

// returns commoner's effort after selective copy
int sc_com(polity *p, group *g, Commoner *com, int i)
{
  int a;  
  if(g->nc > 1){
    do{ a = rnd(g->nc); }while( a == i);                        // select another individual in group
    if(g->com[a].pi > com->pi){                                  // if selected individual has more payoff copy strategies
      return g->com[a].x;
    }
  }
  return com->x; 
}
// changes effort made by commoners
void changeStratCommoner(int v, int h, int j, int i, double y, double pn, double ppc)
/*
 * v: type of strategy change method
 * h: polity
 * j: group
 * i: commoner index
 * y: coordination effort made by leader before leader updates its effort
 * pn: punishment probability
 * ppc: polity production share to each commoner
 */
{  
  polity *p;
  group *g;
  Commoner *com;
  p = Polity+h;
  g = p->g+j;
  com = g->com+i; 
  if(v == 1){                                                               // random mutation
    com->x = rnd(2);                                                        // 0 or 1    
  }
  else if(v == 2){                                                          // copy from peer    
    if(U01() > m){                                                // choose from same group 
      com->x = sc_com(p, g, com, i);      
    }
    else{                                                         // choose commoner from another group c in same polity
      int c;
      if(p->ng > 1){                                              // if more than one group in polity, choose different group
	do{ c = rnd(p->ng); }while(c == j);	
	g = p->g+c;                                               // change target group
      }
      com->x = sc_com(p, g, com, -1);
    }    
  }
  else if(v == 3){                                                          // myopic optimization
    com->x = updateStratCommoner_v3(com->x, g->X, y, pn, g->nc, g->t_u, ppc);
  }
}

// returns leader's effort after selective copy
void sc_lead(polity *p, Leader *ld, int j)
{
  int a;  
  if(p->ng > 1){
    do{ a = rnd(p->ng); }while( a == j);                        // select another individual in group
    if(p->g[a].lead->pi > ld->pi){                                  // if selected individual has more payoff copy strategies
#if UPDATE_LEAD_EFFORT
      ld->y = p->g[a].lead->y;
#endif
#if UPDATE_LEAD_PUN_EFFORT
      ld->p = p->g[a].lead->p;                          // copy punishment effort
#endif
    }
  }  
}

// updates strategy of leader with optimization
void  updateStratLeader_v3(polity *p, group *g, Leader *ld, int j, double ppl, int *l, double *ly, int n_x0)
{
  double y, f, nP, sx, pn;    
  int a, pl;
#if UPDATE_LEAD_EFFORT
    y = normal(ld->y, Sigma);
    y = MAX(y, 0.0);                           // new y    
#else
    y = ld->y;
#endif
#if UPDATE_LEAD_PUN_EFFORT
    pn = normal(ld->p, Sigma);                 // new punishment effort
    pn = MAX(pn, 0.0);
    pn = MIN(pn, 1.0);
#else
    pn = ld->p;
#endif
    sx = g->X;                                // sum of efforts from commoners
    nP = groupProduction(sx, y, g->nc);                                                     // new group production    
    f = pi_leader(y, p->e_u, g->t_u, nP, ppl);                                   // new payoff before punishment    
    // punish commoner if x = 0 with new probability p                      
    f -= pn*delta*n_x0;          // average cost due to punishment to commoner	      
    
    // punish leader if new y is among least L efforts made by leaders
    if(y <= ly[L-1]){                        // j is punishable if y is less than or equal to effort made by Lth leader, 
      f -= (p->chief->q*K);     
    }
    else if(y <= ly[L]){ //y < effort by Lth leader but less than L+1th leader; j could be still punishable if y <= effort made by L+1th leader       
      //check if j is leader among L leaders punished before
      pl = 0;
      for(a = 0; a < L; a++){
	if(l[a] == j){
	  pl = 1;
	  break;
	}
      }
      //if j is punished before, j is still punishable since new y is not greater than effort by Lth leader
      if(pl){
	  f -= (p->chief->q*K);     
      }      
    }      
    if(f > ld->pi){                                    // if new payoff is higher than current payoff, then change to new strategies    
      ld->y = y;
      ld->p = pn;             
    }    
}

// changes effort and taxes by leader
void changeStratLeader(int v, int h, int j, double ppl, int *l, double *ly, int n_x0)
/*
 * v: type of method to change
 * h: polity
 * j: group
 * ppl: polity production share to each leader
 * l: array to index of leaders punished by chief
 * ly: array of efforts by leaders punished by chief
 * n_x0: number of commoners in group with x = 0
 */
{
  int a;
  Leader *ld;
  group *g;    
  g = (Polity+h)->g+j;
  ld = g->lead+0;
  if(v == 1){                                                            // random mutation
#if UPDATE_LEAD_EFFORT
    ld->y = normal(ld->y, Sigma);
    ld->y = MAX(ld->y, 0.0);         
#endif
#if UPDATE_LEAD_PUN_EFFORT
    ld->p = normal(ld->p, Sigma);
    ld->p = MAX(ld->p, 0.0);
    ld->p = MIN(ld->p, 1.0);
#endif
  }  
  else if(v == 2){                                                       // selective copying
    if(U01() > m){                                                 // copy from leaders in same polity
      sc_lead(Polity+h, ld, j);
    }
    else{                                                                  // copy strategy from other leader in other polity
      if(H > 1){
	do{ a = rnd(H);} while( a == h);                                   // choose another polity
      }
      else{
	a = h;
      }
      sc_lead(Polity+a, ld, -1);
    }
  }
  else if(v == 3){                                                       // myopic optimization
    updateStratLeader_v3(Polity+h, g, ld, j, ppl, l, ly, n_x0);
  }  
}

// updates strategy of chief with optimization
void updateStrategyChief_v3(polity *p, Chief *chf, double sY, int uvt, int h1, double sQ)
/*
 * p: polity
 * chf: chief
 * sY: sum of group productions reaching to chief
 */
{
   double z, q, f, nQ;     
#if UPDATE_CHIEF_EFFORT
    z = normal(chf->z, Sigma);                                         // generate candidate strategy effort
    z = MAX(z, 0.0);       
#else
    z = chf->z;
#endif
#if UPDATE_CHIEF_PUN_EFFORT
    q = normal(chf->q, Sigma);
    q = MAX(q, 0.0);
    q = MIN(q, 1.0);
#else
    q = chf->q;
#endif
    nQ = polityStrength(sY, z, p->ng);                                 // new Q value        
    f = pi_chief(z, p->e_d, B*nQ, uvt, h1, sQ-p->Q+nQ);                // new payoff before punishment
    // payoff after punishing L leaders with probability q    
    f -= (L*q*DELTA);
       
    if(f > chf->pi){                                                   // if new payoff is higher than current payoff, then change to new strategies
      chf->z = z;
      chf->q = q;      
    }
}

// change strategy of chief in a polity
void changeStratChief(int v, int h, double sY, int uvt, int h1, double sQ )
/*
 * v: method of change
 * h: polity
 * sY: sum of group production share from leaders reaching to chief
 */
{
  Chief *chf;
  polity *p;
  p = Polity+h;    
  chf = p->chief;
  // begin region: v = 1
  if(v == 1){                                                            // random mutation
#if UPDATE_CHIEF_EFFORT
    chf->z = normal(chf->z, Sigma);
    chf->z = MAX(chf->z, 0.000);            
#endif
#if UPDATE_CHIEF_PUN_EFFORT
    chf->q = normal(chf->q, Sigma);
    chf->q = MAX(chf->q, 0.0);
    chf->q = MIN(chf->q, 1.0);
#endif
  }
  //end region: v = 1
  // begin region: v = 2
  else if(v == 2){                                                       // selective copying
    int k;  
    if(H > 1){
      do{ k = rnd(H); }while( k == h);                                     // select another polity in system
      if(((Polity+k)->chief)->pi > chf->pi){                                 // if chief in selected polity has higher payoff copy strategy
#if UPDATE_CHIEF_EFFORT
	chf->z = (Polity+k)->chief->z;                                 // copy coordination effort	
#endif
#if UPDATE_CHIEF_PUN_EFFORT
	chf->q = (Polity+k)->chief->q;                                 // copy punishment effort
#endif
      }
    }
  }  
  else if(v == 3){                                                       // myopic optimization       
    updateStrategyChief_v3(p, chf, sY, uvt, h1, sQ);
  } 
}


// updates strategy of all chiefs, leaders and commoners
void updateStrategy(int h1, double sQ)
/*
 * s: sum of strength of polities
 */
{
  int h, i, j, k, v, *l, n_x0;
  polity *p;
  group *g;
  double pn, ppc, y, ppl, sY, *ly;
  ly = malloc((L+1)*sizeof(double));
  for(h = 0; h < H; h++){
    p = Polity+h;    
    //l = malloc(p->ng*sizeof(int));
    l = sortedLeaderIndexArrayWithEffort(p);                      // list of indexes of leaders sorted according to effort made in a polity p
    for(k = 0; k < L+1; k++){
      ly[k] = (p->g+l[k])->lead->y;                                     // L+1 least efforts in ascending order made by leaders in polity p
    }
    // update strategy of chief of polity
    for(sY = 0, j = 0; j < p->ng; j++){                                 // through all groups
      g = p->g+j;
      y = g->lead->y;                                              // effort by leader before update
      pn = g->lead->p;
      ppl = (1.0-p->e_d)*n*B*p->Q;                                 // share chief sends to group
      ppc = (1.0 - g->t_d)*ppl/(g->nc);                            // polity production share to commoner
      ppl = g->t_d*ppl;                                            // polity production share that leader keeps 
      sY += g->t_u*g->P;                                           // sum of theta_u * P
      // update strategy of leader of group*/
      for(n_x0 = 0, i = 0; i < g->nc; i++){                                     // through all commoners in group
	if((g->com[i].x == 0)) n_x0++; 		
#if UPDATE_COM
	// update strategy of commoner 
	v = drand(Vcdist);                                             // choose type of update method
        if(v){
	  changeStratCommoner(v, h, j, i, y, pn, ppc);
	}
#endif
      }     

#if UPDATE_LEAD_EFFORT || UPDATE_LEAD_PUN_EFFORT
      v = drand(Vldist);                                               // randomly select method of changing strategy
      if(v){                                                          // change only if v is not 0 
	changeStratLeader(v, h, j, ppl, l, ly, n_x0);                 // change strategy of leader in group g[j] of polity p[h]	
      }
#endif
    }    
#if UPDATE_CHIEF_EFFORT || UPDATE_CHIEF_PUN_EFFORT
    sY *= p->e_u*n*b;                                                 // sum of eta_u*theta_u*b*P
    v = drand(Vcfdist);                                                 // randomly select method of changing strategy
    if(v){                                                            // change only if v is not 0 
      changeStratChief(v, h, sY, UvT[h], h1, sQ);                                     // change strategy of chief of polity p[h]
    }
#endif

    free(l);  
  }
  free(ly);
}

void playGame()
{
  int h, j;  
  double h1, sQ;
  for(h1 = 0, h = 0; h < H; h++){                                          // through all polities   
    // set pool of polities to play us vs them game or us vs nature game for every time step
    if(U01() < Rho){
      UvT[h] = 1;                                                  // plays us vs them game
      h1++;
    }
    else{
      UvT[h] = 0;
    }
  }
  for(sQ = 0, h = 0; h < H; h++){                                          // through all polities       
    // update group production and polity strength
    for(j = 0; j < (Polity+h)->ng; j++){                           // through all groups
      // us vs nature within groups
      updateGroupProduction(h, j);                                 // calculates and update group production p                   
    }       
    updatePolityStrength(h);                                       // updates polity strength of polity h; 
    if(UvT[h]){                                                    // if this polity participates in us vs them game
      sQ += (Polity+h)->Q;                                         // sum of polity strength of polities in us vs them game pool 
    }
  }
  for(h = 0; h < H; h++){
    updatePayoffBeforePunishment(h, h1, sQ);                               // updates payoff of each commoners and leader after production
    punish(h);                                                     // punish leaders and commoners in polity h
  }    
  // update strategies /  efforts
  updateStrategy(h1, sQ);
}

/*****************************************************************************************************************************************/
/*****************************************************************************************************************************************/


void init()                                                          // initialize values of global variables
{
  int h, i, j;
  Commoner *cm;
  Leader *ld;
  Chief *cf;
  polity *p;
  group *gr;  
  for(h = 0; h < H; h++){
    p = Polity + h;
    p->ng = G;                                               // set no. of groups in a polity
    p->e_u = Eta_u;                                            // set tax on leader
    p->e_d = Eta_d;
    p->Q = 0;                                                // set polity Q value to 0
    cf = p->chief;
    cf->pi = 0;                                                       // reset payoff of chief
    cf->z = INIT_CHIEF_EFFORT;                                        // set coordination efforts by chief
    cf->q = INIT_CHIEF_PUN_EFFORT;
    for(j = 0; j < G; j++){
      gr = p->g+j;
      gr->nc = n;                                                     // no. of commoners in a group set to N
      gr->P = 0;                                                      // reset group production
      gr->X = 0;                                                      // reset average group effort
      gr->t_u = Theta_u;                                              // tax on commoners set
      gr->t_d = Theta_d;
      ld = gr->lead;
      ld->p = INIT_LEAD_PUN_EFFORT;                                   // punishment effort by leader
      ld->pi = 0;                                                     // reset payoff of leader
      ld->y = INIT_LEAD_EFFORT;                                       // initialize coordiation effort by leader
      for(i = 0; i < n; i++){
	cm = gr->com+i;
	cm->x = INIT_COM_EFFORT;                                   // initialize effort
	cm->pi = 0;                                                   // reset payoff	
      }      
    }
  }  
  // set strategy update probability according to option (Vop) in config  
  // initialize strategy update method probability distribution
  j = Vop - 1;
  Vcdist->p[0] = 1.0-Vc[j][0]-Vc[j][1]-Vc[j][2];
  Vcdist->p[1] = Vc[j][0]; 
  Vcdist->p[2] = Vc[j][1];
  Vcdist->p[3] = Vc[j][2];
  initdist(Vcdist, 1.0);  
  Vldist->p[0] = 1.0-Vl[j][0]-Vl[j][1]-Vl[j][2];
  Vldist->p[1] = Vl[j][0];
  Vldist->p[2] = Vl[j][1];
  Vldist->p[3] = Vl[j][2];
  initdist(Vldist, 1.0); 
  Vcfdist->p[0] = 1.0 - Vcf[j][0]-Vcf[j][1]-Vcf[j][2];
  Vcfdist->p[1] = Vcf[j][0];
  Vcfdist->p[2] = Vcf[j][1];
  Vcfdist->p[3] = Vcf[j][2];
  initdist(Vcfdist, 1.0); 
}


int main(int argc, char **argv)
{
#if DEBUG
  feenableexcept(FE_DIVBYZERO| FE_INVALID|FE_OVERFLOW); // enable exceptions
#endif
  if(argc ^ 2){
    printf("Usage: ./mlp mlp.config\n");
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
  int r, i;
  unsigned long seed;  
#if !ALLDATAFILE
  char xdata[200], str[30];
#endif
  time_t now;     
  // print headers for values to be displayed on std output
#if !CLUSTER
  printf("\nValues:\t  x\t  y\t  z\t  p\t  q\t theta_u  theta_d  eta_u  eta_d\t  pi_0\t  pi_1\t  pi_2\t  P\t  Q\t  seed\n");
#endif
  
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
    }
    calcStat(-1, -1);                                        // free file pointers for individual run data files
#if !CLUSTER
    
#if !AVERAGE_GRAPH_ONLY
    if(Runs > 1)
      plotIndividualRun(r, 0); 
#if !TURNOFF_TEST_PLOT
      //plotTraitsDist(r);
#endif
#endif
    
#if GRAPHS
    if(Runs > 1)
      plotIndividualRun(r, 1);
#endif
    // remove individual datafiles if ALLDATAFILE set to 0
#if !ALLDATAFILE        
     sprintf(str, "x%d.dat", r); prep_file(xdata, str); remove(xdata);
      sprintf(str, "p%d.dat", r); prep_file(xdata, str); remove(xdata);
      sprintf(str, "t%d.dat", r); prep_file(xdata, str); remove(xdata);
      sprintf(str, "gp%d.dat", r); prep_file(xdata, str); remove(xdata);
      sprintf(str, "pun%d.dat", r); prep_file(xdata, str); remove(xdata);     
#endif
#endif
    cleanup();                                                // free all memory allocated for polity system
    
#if !CLUSTER
#if !AVERAGE_GRAPH_ONLY
#if !TURNOFF_TEST_PLOT
    //gp = popen("gnuplot -persistent", "w");
   // plotLines(gp, "grtest.dat", 8, "traits evolution");
    //fflush(gp);
    //pclose(gp);
    //plotTest("grtest.dat", 8, "traits evolution");
#endif
#endif
#endif
  }
  writeDataToFile();  
  clearStatVar();                                             // free other memory allocated for statistics variables
#if !CLUSTER
  plotAverage(0);
#if GRAPHS
  plotAverage(1);
#endif
#endif    
  return 1;
}

/* 
gcc -Wall -O2 -march=native -pipe -o mlp mlp.c -lm
./mlp mlp.config
valgrind -v --track-origins=yes --leak-check=full --show-leak-kinds=all ./mlp mlp.config
gcc -g -o mlp mlp.c -lm
gdb mlp
run mlp.config

//profiling code
gcc -Wall -O2 -march=native -pipe -pg mlp.c -o mlp -lm
./mlp mlp.config
 gprof mlp gmon.out > analysis.txt  
*/

