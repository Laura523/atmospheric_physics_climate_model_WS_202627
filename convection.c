// Goal week 1:
// * start with a temperature profile with p, T. The pressure profile should be equidistant in p. Assume an initial temperature 
// profile. If you want, you can convert to altitude z(p) using the hydrostatic assumption, delta_p / delta_z = - rho * g, 
// starting from z=0 at ground level.
// * convert to potential temperature and do convection
// * convert back to temperature

// Since this can become boring quickly, the next step would be to add some continuous heating to the layer next to the surface. 
// The idea is to assume that the complete absorbed solar energy Eabs = 240 W/m2 is absorbed by the surface and then heats the 
// lowermost layer:

// initialize temperature profile

// while (1==1) { // forever

// apply solar heating to lowermost layer assuming an appropriate time step

// do convection

// print profile
// }

//Does it converge? Do the results depend on your initial assumption of a temperature profile?





#include <stdio.h>
#include <math.h>
#include <stdbool.h>

#define G 9.80665
#define CP 1004       // heat capacity for constant pressure c_p [J/kg/K]
#define KAPPA 2.0     // macro variable (2.0 or 7.0)

// swap function
void swapdouble(double *a, double *b)
{
    double *temp = a;
    a = b;
    b = temp;
}

// bubble sort
void bubblesort(int length, double *variable)
{
    for (int j=0; j<length-1; j++)
    {
        bool swapped = false;
        for (int i=0; i<length-j-1; i++)
        {
            if (variable[i]>variable[i+1])
            {
                swapdouble(&variable[i], &variable[i+1]);
                swapped = true;
            }
        }
        if (!swapped)
            break;
    }
}

// convert T to potential temperature
double TtoTheta (double p, double T) 
{
    double p0 = 1000.0;
    double theta = T * pow(p0/p, KAPPA);
    return theta;
}


// Bernhards version
int main ()
{
    int nlyr = 10;                          // number of layers
    double dp = 1000.0 / (double) nlyr;     // pressure increment

    // decleration of variables
    double plev[nlyr+1];
    double plyr [nlyr];
    double Tlyr[nlyr];
    double theta[nlyr];

    // initialize pressure and temperature profile
    for (int i=0; i<nlyr+1; i++)
        plev[i] = i*dp;

    for (int i=0; i<nlyr; i++)
    {
        plyr[i] = (plev[i]+plev[i+1])/2.0;
        Tlyr[i] = 288;                      // asssume equilibrium
    }

    // calculate initial potential temperature
    for (int i=0; i<nlyr; i++)
        theta[i] = TtoTheta(plyr[i], Tlyr[i]);

    // print inital layer values
    printf("initial layer values:\n");
    for (int i=0; i<nlyr; i++)
        printf("i=%d p=%f theta=%f\n", i, plyr[i], theta[i]);
    // now we do the convection part here...
}






// old version
/*
int main (int argc, char **argv)
{
    // Define variables
    int nlayer = 5;                 // number of levels (space inbetween called layers)
    double Eabs = 240;              // absorbed solar energy [W/m^2]
    double p[nlayer+1];             // pressure [hPa]
    double theta[nlayer];           // potential temperature
    double dp = 1000.0 / nlayer;    // pressure difference
    double dt = 10 * 60;            // [s]

    for (int i=0; i<nlayer+1; i++)  // give each level correct pressure value
        p[i] = i * dp;

    for (int i=0; i<nlayer; i++)    // give each layer correct potetial temperature value
        theta[i] = 288;

    for (int i=0; i<nlayer+1; i++)  // print all levels and their respective pressure 
        printf ("%d %f\n", i, p[i]);

    // heating step
    double dT = Eabs * dt * G / dp / CP;

    return 0;
}
*/