// Goal week 1:
// * start with a temperature profile with p, T. The pressure profile should be
// equidistant in p. Assume an initial temperature profile. If you want, you can
// convert to altitude z(p) using the hydrostatic assumption, delta_p / delta_z
// = - rho * g, starting from z=0 at ground level.
// * convert to potential temperature and do convection
// * convert back to temperature

// Since this can become boring quickly, the next step would be to add some
// continuous heating to the layer next to the surface. The idea is to assume
// that the complete absorbed solar energy Eabs = 240 W/m2 is absorbed by the
// surface and then heats the lowermost layer:

// initialize temperature profile

// while (1==1) { // forever

// apply solar heating to lowermost layer assuming an appropriate time step

// do convection

// print profile
// }

// Does it converge? Do the results depend on your initial assumption of a
// temperature profile?

#include <math.h>
#include <stdbool.h>
#include <stdio.h>

#define G 9.80665
#define CP 1004             // heat capacity for constant pressure c_p [J/kg/K]
#define KAPPA (2.0 / 7.0)   //
#define alpha 0.15            // absorptivity, may be varied
#define sigma 5.670e-8         // Stefan Boltzmann Konstante

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

double ThetatoT(double p, double theta) {
  double p0 = 1000.0;
  double T = theta / pow(p0 / p, KAPPA);
  return T;
}

int main() {
  int nlyr = 10;                     // number of layers
  double dp = 1000.0 / (double)nlyr; // pressure increment

  // decleration of variables
  double plev[nlyr + 1];
  double plyr[nlyr];
  double Tlyr[nlyr];
  double TafterConvection[nlyr];
  double theta[nlyr];

  // solar heating parameters
  double Eabs = 240.0; // W/m^2
  double dt = 3600.0;  // time step: 1 hour in seconds
  double dT = Eabs * G * dt / (CP * dp);
  int nsteps = 100;

  // initialize pressure and temperature profile
  for (int i = 0; i < nlyr + 1; i++)
    plev[i] = i * dp;

  for (int i = 0; i < nlyr; i++) {
    plyr[i] = (plev[i] + plev[i + 1]) / 2.0;
    Tlyr[i] = 288; // asssume equilibrium
  }

  // calculate initial potential temperature
  for (int i = 0; i < nlyr; i++)
    theta[i] = TtoTheta(plyr[i], Tlyr[i]);

  // print inital layer values
  printf("initial layer values:\n");
  for (int i = 0; i < nlyr; i++)
    printf("i=%d p=%f T=%f theta=%f\n", i, plyr[i], Tlyr[i], theta[i]);

  // make plot to plot t after n timesteps, seemingly not possible in only c -
  // make file and move to python in python: ipywidgets daraus widgets dann
  // definde slider um nicht jeden timestep seperat zu plotten make results file

  // radiative cooling for ground layer
  double dQ = (Eabs - signgam * pow(Tlyr[9], 4)) * alpha * dT;
  double dTcool = dQ / CP;

  // timeloop for convection and heating
  for (int step = 0; step < nsteps; step++) {
    Tlyr[9] += dT - dTcool; // x += y entspricht x = x+y dh dT wird immer wieder
                            // zu T9 addiert
    // Convert temperature to potential temperature
    for (int i = 0; i < nlyr; i++)
      theta[i] = TtoTheta(plyr[i], Tlyr[i]);

    // Convective adjustment
    bubblesort(nlyr, theta);

    // Convert back to temperature
    for (int i = 0; i < nlyr; i++)
      TafterConvection[i] = ThetatoT(plyr[i], theta[i]);

    // Update temperature profile for the next time step
    for (int i = 0; i < nlyr; i++)
      Tlyr[i] = TafterConvection[i];
  }
  { printf("\nLayer values after %d timesteps:\n", nsteps); }
  // print profile after the timesteps for any number of timesteps
  int i;
  for (i = 0; i < nlyr; i++) {
    theta[i] = TtoTheta(plyr[i], Tlyr[i]);

    printf("i=%d p=%f T=%f theta=%f\n", i, plyr[i], Tlyr[i], theta[i]);
  }

}

// after doing things normally you can use : ./convection > convectionres.csv
// to get a file with the end results