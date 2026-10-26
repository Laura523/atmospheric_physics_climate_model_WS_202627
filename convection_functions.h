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

#endif