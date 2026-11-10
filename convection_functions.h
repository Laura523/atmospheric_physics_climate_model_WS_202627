#ifndef FUNCTIONS_H
#define FUNCTIONS_H

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define G 9.80665
#define CP 1004           // heat capacity for constant pressure c_p [J/kg/K]
#define R 287.0           // specific gas constant of dry air [J/kg/K]
#define KAPPA (2.0 / 7.0) // kappa = R/c_p = 2/7 (previously 2.0)
#define SIGMA 5.670e-8    // Stefan-Boltzmann constant [W/m^2/K^4]

void swapdouble(double *a, double *b);
void bubblesort(int length, double *variable);
double TtoTheta(double p, double T);
double ThetatoT(double p, double theta);
void calcHeight(int nlyr, double dp, double *plyr, double *Tlyr, double *zlyr);
void calcRadiation(int nlyr, int mode, double alpha, double Eabs, double *Tlyr, double *Eup, double *Edn, double *dQ);
void schwarzschild(int nlyr, double *Tlyr,double alpha, double *Eup,double *Edn);

#endif

//baustellen: schwarzschild mit winkelintegral
// wellenlängenabhängigkeit (optische  dicke): tau fkt von höhe, auch von wellenlänge - 2d array 
// für sb mü =0.5 (alpha)
//double from [3] [nlayr] (3 wg 3 wellenlängenbereichen)
// siga t hoch 4s ersetzen durch planckk integral für wellenlängenberiehc 0 bis 8, 8 bis  12, 12 bis unendl. mikro
// integral: pi int (b(landa) dlanda)