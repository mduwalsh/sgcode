// to run on cluster

#include<stdio.h>
#include<stdlib.h>
#include<unistd.h>
#include<string.h>
#include<time.h>

#define FileName "images_1.html"
#define TITLE "v1 0.01 v2 0.24 v3 0.0"
#define WIDTH 400
#define HEIGHT 225

/* sets of parameters to be simulated on */
int Vop[]      = {1, 2, 3};
int n[]        = {4};
double b[]     = {1, 2, 3};
double B[]     = {1};
double delta[] = {0};
double DELTA[] = {0};
double K[]     = {0};
double k[]     = {0};

double Theta_u[] = {0, 0.1, 0.2};
double Eta_u[]   = {0};
double Theta_d[] = {0};
double Eta_d[]   = {0};

double cx[]    = {1};
double cy[]    = {1};
double cz[]    = {0};

double x_0[]   = {1};    // x0
double y_0[]   = {0.125, 0.25, 0.5};          // y0
double z_0[]   = {1};          // z0
double X0[]    = {2, 4, 8, 16};
double Y0[]    = {1};
double e[]     = {1, 2, 4};
double E[]     = {1};
double Rho[]   = {0.4};


/* end of sets of parameters to be simulated on */

/* end of other basic parameters */
int i[24], is[24];
int _n, _Vop; 
double _b, _B, _delta, _DELTA, _K, _k, _Theta_u, _Theta_d, _Eta_u, _Eta_d, _cx, _cy, _cz, _Rho, _x_0, _y_0, _z_0, _X0, _Y0, _e, _E;

void prep_file_plot(char *fname, char *prpndStr, char *apndStr)
/*
 * fname: variable in which name of file to be stored
 * apndStr: string to be appended to file name string
 */
{
  sprintf(fname, "%stu%.2ftd%.2feu%.2fed%.2fcx%.2fcy%.2fcz%.2fVop%dRho%.2fY0%.2fe%.2fE%.2fx0%.2fy0%.2fz0%.2f%s", prpndStr, _Theta_u, _Theta_d, _Eta_u, _Eta_d, _cx, _cy, _cz, _Vop, _Rho, _Y0, _e, _E, _x_0, _y_0, _z_0, apndStr);
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
  //write first row of table / column headers  // outside x axis in html page
  fprintf(fp, "<tr class=\"rowsinfo\">\n");
  fprintf(fp, "<td></td>\n");
  for(i[3] = 0; i[3] < is[3]; i[3]++){  // delta
    _delta = delta[i[3]];
    fprintf(fp, "<td> &delta; = %g </td>\n", _delta);
  }
  fprintf(fp, "</tr>\n");  

  for(i[2] = 0; i[2] < is[2]; i[2]++){  // B
    //_B = B[i[2]];
    
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
					// create a row of images
					fprintf(fp, "<tr>\n");
					fprintf(fp, "<td class=\"rowsinfo\">\n");
					fprintf(fp, "&eta;:%g, z0:%g, E:%g, e:%g\n", _Eta_u, _z_0, _E, _e);
					fprintf(fp, "</td>\n");
					//major x axis columns in images.html layout file
					for(i[3] = 0; i[3] < is[3]; i[3]++){  // delta
					  _delta = delta[i[3]];
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
  //int ic1, ic2;

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
  
  write_images_html(TITLE);
  
  return 0;
}




/** Usage:
    compile : gcc -o imagescript imagescript.c
    run     : ./imagescript
**/



