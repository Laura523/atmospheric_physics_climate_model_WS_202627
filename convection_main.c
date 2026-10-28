#include "convection_functions.h"
#include <math.h>
#include <stdbool.h>
#include <stdio.h>


int main() {
  // declaration of variables
  int nlyr = 10;                     // number of layers
  double dp = 1000.0 / (double)nlyr; // pressure increment [hPa]
  double dpPa = dp * 100.0;          // pressure increment [Pa]

  double plev[nlyr + 1]; // pressure level
  double plyr[nlyr];     // pressure layer
  double zlyr[nlyr];     // height layer
  double Tlyr[nlyr];     // temp layer
  double theta[nlyr];    // potential temp layer

  double Told[nlyr];             // temp layer first timestep
  double TafterConvection[nlyr]; // temp layer after convection

  double Eup[nlyr]; // * upward thermal flux at top of each layer [W/m^2]
  double Edn[nlyr]; // * downward thermal flux at top of each layer [W/m^2]
  double dQ[nlyr];  // * net energy gain per layer [W/m^2]
  int nangles = 180; // angle for sun radiation dependence

  
  // initialize pressure profile
  for (int i = 0; i < nlyr + 1; i++)
    plev[i] = i * dp;

  for (int i = 0; i < nlyr; i++)
    plyr[i] = (plev[i] + plev[i + 1]) / 2.0;

  /* Task Week 1: heat lowest layer by absorbed solar energy and cool by
  emitting black body radiation. do convection. stop run when the single-column
  adjusted. */

  // solar heating parameters
  double Eabs = 240.0; // W/m^2
  double dt = 36000.0;  // time step: 10 hour in seconds

  double eps = 1e-6;      // convergence criterion
 int maxsteps = 1000000; // safety abort

  // * Task Week 2: grey thermal radiation model
  // * mode 0 = step (0) solar heating only, mode 1 = step (1) surface cooling,
  // * mode 2 = step (2) full grey radiation in the atmosphere
  int mode = 2; // *
  // * steps (3)/(4): same alpha for every layer. Tune it until the equilibrium
  // * surface temperature is 288 K (step 3), then increase it until the surface
  // * is 1.4 K warmer (step 4). Rerun after every change.
  double alpha = 0.15; // *

  // test different profiles
  double Tstart[3] = {200.0, 288.0, 350.0};
  int nruns = 3;

  for (int run = 0; run < nruns; run++) {

    // initialize temperature profile and potential temperature
    for (int i = 0; i < nlyr; i++) {
      Tlyr[i] = Tstart[run];
      theta[i] = TtoTheta(plyr[i], Tlyr[i]);
    }

    calcHeight(nlyr, dp, plyr, Tlyr, zlyr);

    printf("\n=== run %d: start temp %.0f K ===\n", run, Tstart[run]);
    printf("initial layer values:\n");
    for (int i = 0; i < nlyr; i++)
      printf("i=%d p=%f z=%f T=%f theta=%f\n", i, plyr[i], zlyr[i], Tlyr[i],
             theta[i]);

    // * open a new csv file for every run
    char filename[100]; // *
    snprintf(filename, sizeof(filename), "run_%d_mode%d_alpha%.3f.csv", run,
             mode, alpha);                  // *
    FILE *csv = fopen(filename, "w");       // *
    if (csv == NULL) {                      // *
      printf("Could not open %s\n", filename); // *
      return 1;                             // *
    }                                       // *

    // * csv header
    fprintf(csv, "step,Eup_sur,Edn_sur,dQ_sur"); // *
    for (int i = 0; i < nlyr; i++)               // *
      fprintf(csv, ",T%d", i);                   // *
    for (int i = 0; i < nlyr; i++)               // *
      fprintf(csv, ",theta%d", i);               // *
    fprintf(csv, "\n");                          // *

    int step = 0;
    bool converged = false;
    bool toolong = false;
    

    // Solar-angle dependence
    double totalCos = 0.0;

    for (int i = 0; i < nangles; i++) {
        double theta_angle =
            (i + 0.5) * (M_PI / 2.0) / nangles;

        totalCos += cos(theta_angle);
    }

    double meanCos = totalCos / nangles;

    // Normalize so Eabs remains the prescribed mean input
    double EabsAngle[nangles];

    for (int i = 0; i < nangles; i++) {
        double theta_angle =
            (i + 0.5) * (M_PI / 2.0) / nangles;

        EabsAngle[i] =
            Eabs * cos(theta_angle) / meanCos;
    }

    // heating, cooling and convection loop with directional dependence (yet to be completed)
    while (1 == 1) {

      for (int i = 0; i < nlyr; i++)
        Told[i] = Tlyr[i]; // start profile

      // * radiation: solar heating, thermal cooling of the surface and
      // * thermal heating/cooling of the atmosphere (depending on mode)
      double EabsNow = EabsAngle[step % nangles];

      calcRadiation(nlyr, mode, alpha, EabsNow,
              Tlyr, Eup, Edn, dQ);
      for (int i = 0; i < nlyr; i++)                              // *
        Tlyr[i] += dQ[i] * G * dt / (CP * dpPa);                  // *
      
        
      for (int i = 0; i < nlyr; i++)
        theta[i] = TtoTheta(plyr[i], Tlyr[i]); // T to pot temp

      bubblesort(nlyr, theta); // convection

      for (int i = 0; i < nlyr; i++)
        TafterConvection[i] = ThetatoT(plyr[i], theta[i]); // pot temp to T

      for (int i = 0; i < nlyr; i++)
        Tlyr[i] = TafterConvection[i]; // update T for next step

         
      step++;

      // * write all values of this time step to the csv file
      fprintf(csv, "%d,%f,%f,%f", step, Eup[nlyr - 1], Edn[nlyr - 1],
              dQ[nlyr - 1]);                    // *
      for (int i = 0; i < nlyr; i++)            // *
        fprintf(csv, ",%f", Tlyr[i]);           // *
      for (int i = 0; i < nlyr; i++)            // *
        fprintf(csv, ",%f", theta[i]);          // *
      fprintf(csv, "\n");                       // *

      // convergence test
      double maxdiff = 0.0;
      for (int i = 0; i < nlyr; i++) {
        double diff = fabs(Tlyr[i] - Told[i]);
        if (diff > maxdiff)
          maxdiff = diff;
      }

      converged = (maxdiff < eps);
      toolong = (step >= maxsteps);

      if (converged || toolong)
        break;
    }

    fclose(csv); // *

    calcHeight(nlyr, dp, plyr, Tlyr, zlyr);

    printf("\nStep %d:\n", step);
    if (converged)
      printf("Converged after %d steps\n", step); // *
    else
      printf("Not converged after %d steps\n", step); // *

    for (int i = 0; i < nlyr; i++) {
      theta[i] = TtoTheta(plyr[i], Tlyr[i]);
      printf("i=%d p=%f z=%f T=%f theta=%f\n", i, plyr[i], zlyr[i], Tlyr[i],
             theta[i]);
    }

    printf("Start %.0f K -> surface T at equilibrium: %f K\n", Tstart[run], // *
           Tlyr[nlyr - 1]);
    printf("Data written to %s\n", filename); // *
  }

  return 0;
}

// TODO: get a proper csv file for plotting (done: one csv file per run)  // *
// Write down what was learned from the task                              // *
// system converges to the same temperature indepent ent of the starting temp;alpha of 0.15 for 288 K  //*
// 1,5 grad ziel bei alpha 0.156 etwa//*
