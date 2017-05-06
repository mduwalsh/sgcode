// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include<time.h>

#define FileName "images_1.html"
#define WIDTH 400
#define HEIGHT 225

/* sets of parameters to be simulated on */
int Vsim = 3;

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
double cy[]    = {0};
double cz[]    = {1};

double X0[]    = {2, 4, 8};    // X0
double z_0[]   = {1.0/16.0, 1./8., 1./4.};          // z0
double e[]     = {1, 2, 4};

/* end of sets of parameters to be simulated on */

/* other basic parameters */
unsigned Seed      = 0;     
int      Runs      = 10;  

int      G         = 500;   
int      H         = 1;
int      T         = 1500;   
double   Lambda    = 5; 
double   m         = 0.5;

double   Sigma     = 0.1; 
double   Sigma_t   = 0.1; 

/* end of other basic parameters */
int i[24], is[24];
int _n, _K; 
double _b, _delta, _k, _Theta, _cx, _cy, _cz, _z_0, _X0, _e;

char TITLE[300];


void prep_file_plot(char *fname, char *prpndStr, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "%sK%02dcx%.2fcy%.2fcz%.2fvc1%.2fvc2%.2fvc3%.2fvl1%.2fvl2%.2fvl3%.2fk%.2fdl%.2fe%.2fz0%.2f%s", prpndStr, _K, _cx, _cy, _cz, Vc[0], Vc[1], Vc[2], Vl[0], Vl[1], Vl[2], _k, _delta, _e, _z_0, apndStr);
}

void write_image_thumbnail_holder(FILE *fp, char *ifile)
{
  fprintf(fp, "<a href=\"%s\">\n", ifile);
  fprintf(fp, "<img src=\"thumbs/%s.jpg\" width=\"%d\" height=\"%d\" \n", ifile, WIDTH, HEIGHT);
  fprintf(fp, "alt=\"thumbs/%s.jpg\"/></a>\n", ifile);
  fprintf(fp, "<div>%s</div>\n", ifile);  
}

void write_image_table_format(FILE *fp, char *type)
{
  char ifile[300];
  // start table
  fprintf(fp, "\n<table>\n");
  //write first row of table / column headers
  fprintf(fp, "<tr class=\"rowsinfo\">\n");
  fprintf(fp, "<td></td>\n");
  for(i[10] = 0; i[10] < is[10]; i[10]++){  // z_0
    _z_0 = z_0[i[10]];
    fprintf(fp, "<td> z0 = %g </td>\n", _z_0);
  }
  fprintf(fp, "</tr>\n");  

  for(i[2] = 0; i[2] < is[2]; i[2]++){  // K
    _K = K[i[2]];       
    for(i[6] = 0; i[6] < is[6]; i[6]++){  // cx
      _cx = cx[i[6]];
      for(i[7] = 0; i[7] < is[7]; i[7]++){  // cy
	_cy = cy[i[7]];
	for(i[8] = 0; i[8] < is[8]; i[8]++){  // cz
	  _cz = cz[i[8]];			    		  
	  for(i[4] = 0; i[4] < is[4]; i[4]++){  // k
	  _k = k[i[4]];	
	    for(i[3] = 0; i[3] < is[3]; i[3]++){  // delta
	      _delta = delta[i[3]];
	      for(i[11] = 0; i[11] < is[11]; i[11]++){  // e
		_e = e[i[11]];  
		// create a row of images
		fprintf(fp, "<tr>\n");
		fprintf(fp, "<td class=\"rowsinfo\">\n");
		fprintf(fp, "k:%g, &delta;:%g, e:%g\n", _k, _delta, _e );
		fprintf(fp, "</td>\n");
		//major y axis columns in images.html layout file
		for(i[10] = 0; i[10] < is[10]; i[10]++){  // z_0
		  _z_0 = z_0[i[10]];
		  // get image file name
		  prep_file_plot(ifile, type, ".png");
		  fprintf(fp, "<td align='center'>\n");
		  write_image_thumbnail_holder(fp, ifile);
		  fprintf(fp, "</td>\n");
		}
		fprintf(fp, "</tr>\n");
		// end of a row
	      }
	    }
	  }
	}
      }
    }
  }   
  fprintf(fp, "</table>\n");  
}

void write_images_header(FILE *fp, char *title)
{
  fprintf(fp, "<?xml version=\"1.0\" encoding=\"UTF-8\" ?> \n");
  fprintf(fp, "<!DOCTYPE html PUBLIC \"-//W3C//DTD XHTML 1.1//EN\" \"http://www.w3.org/TR/xhtml11/DTD/xhtml11.dtd\"> \n");
  fprintf(fp, "<html xmlns=\"http://www.w3.org/1999/xhtml\"> \n");
  fprintf(fp, "<head> \n");
  fprintf(fp, "<title>%s</title>\n", title);
  fprintf(fp, "<meta http-equiv=\"content-type\" content=\"text/html; charset=UTF-8\"/> \n");
  fprintf(fp, "<meta name=\"GENERATOR\" content=\"KDE Konqueror KImgallery plugin version 4.13.3\"/> \n");
  fprintf(fp, "<style type='text/css'>\n");
  fprintf(fp, "BODY {color: #d0ffd0; background: #333333; \n");
  fprintf(fp, "font-family: Ubuntu, sans-serif; \n");
  fprintf(fp, "font-size: 8pt; margin: 8%%; } \n");
  fprintf(fp, "H1       {color: #d0ffd0;} \n");
  fprintf(fp, "TABLE    {text-align: center; margin-left: auto; margin-right: auto;} \n");
  fprintf(fp, "TD       { color: #d0ffd0; padding: 1em} \n");
  fprintf(fp, "IMG      { border: 1px solid #d0ffd0; } \n");
  fprintf(fp, ".layoutinfo {font-size: 11pt;} \n");
  fprintf(fp, ".rowsinfo {font-size: 12pt; font-weight: bold;} \n");
  fprintf(fp, "</style>\n");
  fprintf(fp, "</head> \n");
  fprintf(fp, "\n");
}

void write_images_html(char *title)
{
  FILE *fp;  
  fp = fopen(FileName, "w+");
  int ic1, ic2;

  // get time
  time_t t = time(NULL);
  struct tm *tm = localtime(&t);
  
  write_images_header(fp, title);
  fprintf(fp, "\n<body>\n");
  fprintf(fp, "<h1>%s</h1><p>\n", title);
  fprintf(fp, "<i>Created on</i>: %s</p> \n", asctime(tm));
  fprintf(fp, "<hr/>\n");
  write_image_table_format(fp, "ef");  
  write_image_table_format(fp, "py"); 
  fprintf(fp, "\n</body>\n");
  fprintf(fp, "</html>\n");
  
  fclose(fp);
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
  
  
  if(Vsim == 1){ Vc = Vc1; Vl = Vl1;}
  else if(Vsim == 2){ Vc = Vc2; Vl = Vl2;}
  else if(Vsim == 3){ Vc = Vc3; Vl = Vl3;}  
  sprintf(TITLE, "Vc_%.2f_%.2f_%.2f Vl_%.2f_%.2f_%.2f", Vc[0], Vc[1], Vc[2], Vl[0], Vl[1], Vl[2]);
  write_images_html(TITLE);
  
  return 0;
}




/** Usage:
    compile : gcc -o imagescript imagescript.c
    run     : ./imagescript
**/



 
