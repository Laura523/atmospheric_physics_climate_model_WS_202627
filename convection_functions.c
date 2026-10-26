#include "convection_functions.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>

// swap function
void swapdouble(double *a, double *b) {
  double temp = *a;
  *a = *b;
  *b = temp;
}

// bubble sort
void bubblesort(int length, double *variable) {
  for (int j = 0; j < length - 1; j++) {
    bool swapped = false;
    for (int i = 0; i < length - j - 1; i++) {
      if (variable[i] < variable[i + 1]) {
        swapdouble(&variable[i], &variable[i + 1]);
        swapped = true;
      }
    }
    if (!swapped)
      break;
  }
}

// convert T to potential temperature
double TtoTheta(double p, double T) {
  double p0 = 1000.0;
  double theta = T * pow(p0 / p, KAPPA);
  return theta;
}

// convert potential temperature to T
double ThetatoT(double p, double theta) {
  double p0 = 1000.0;
  double T = theta / pow(p0 / p, KAPPA);
  return T;
}

// hydrostatic height
void calcHeight(int nlyr, double dp, double *plyr, double *Tlyr,
                double *zlyr) {
  double zlev = 0.0;
  for (int i = nlyr - 1; i >= 0; i--) {
    double dz = R * Tlyr[i] * dp / (plyr[i] * G);
    zlyr[i] = zlev + dz / 2.0;
    zlev += dz;
  }
}

// * grey thermal radiation model (layer 0 = top, layer nlyr-1 = surface layer)
// * mode 0: only absorbed solar radiation heats the lowest layer
// * mode 1: additionally the surface cools by sigma*T_sur^4
// * mode 2: additionally every layer absorbs/emits thermal radiation
// * (alpha per layer) and the surface receives the downward thermal radiation
// * Eup[i], Edn[i]: upward/downward flux at the top of layer i [W/m^2]
// * dQ[i]: net energy gain per area of layer i [W/m^2]
void calcRadiation(int nlyr, int mode, double alpha, double Eabs,
                   double *Tlyr, double *Eup, double *Edn, double *dQ) { // *
  for (int i = 0; i < nlyr; i++) { // *
    Eup[i] = 0.0;                  // *
    Edn[i] = 0.0;                  // *
    dQ[i] = 0.0;                   // *
  }                                // *

  double Esur = SIGMA * pow(Tlyr[nlyr - 1], 4); // * surface emission

  if (mode == 0) {         // *
    dQ[nlyr - 1] = Eabs;   // *
    return;                // *
  }                        // *
  if (mode == 1) {                 // *
    Eup[nlyr - 1] = Esur;          // *
    dQ[nlyr - 1] = Eabs - Esur;    // *
    return;                        // *
  }                                // *

  // * downward flux: from the top of the atmosphere (no incoming thermal
  // * radiation) down to the surface
  Edn[0] = 0.0; // *
  for (int i = 0; i < nlyr - 1; i++)                                  // *
    Edn[i + 1] = Edn[i] * (1.0 - alpha) +
                 alpha * SIGMA * pow(Tlyr[i], 4);                     // *

  // * upward flux: from the surface up to the top of the atmosphere
  Eup[nlyr - 1] = Esur; // *
  for (int i = nlyr - 2; i >= 0; i--)                                 // *
    Eup[i] = Eup[i + 1] * (1.0 - alpha) +
             alpha * SIGMA * pow(Tlyr[i], 4);                         // *

  // * atmospheric layers: absorb from above and below, emit up and down
  for (int i = 0; i < nlyr - 1; i++)                                  // *
    dQ[i] = alpha * (Eup[i + 1] + Edn[i]) -
            2.0 * alpha * SIGMA * pow(Tlyr[i], 4);                    // *

  // * surface: solar - thermal emission + downward thermal radiation
  dQ[nlyr - 1] = Eabs - Eup[nlyr - 1] + Edn[nlyr - 1]; // *
}