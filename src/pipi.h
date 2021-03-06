#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <string>
#include <sstream>
#include <stdlib.h>

using namespace std;
typedef complex<double> cd;

extern double mres;
const double pi = 3.1415926535897932384626433832795028841972;
const double conv = (pi / 180.);

//Masses
const double mPi = 0.1396;
const double mK = 0.496;
const double mEta = 0.54753;

const double mRho = .77545;
const double mF2 = 1.2754;

//Thresholds for pi, eta, and K
const double sthPi = 4.*mPi*mPi;
const double sthK = 4.*mK*mK;
const double sthEta = 4.*mEta*mEta;

//Unit imaginary and real
const cd xr(1., 0.);
const cd xi(0., 1.);

#define backN       6.*xr   //Background N, just needs to be "large enough"
#define maxN        3       //Truncated maximum n
#define DPOINTS     100    // Number of points to have in plotting functions
#define INTP        100     //Number of points for numerical integration

cd rtraj(double alph[], double s);
cd ctraj(double alph[], double s);
cd n_amp(int n, double alph[], double coupling[], double s, double t);


double phase_shift(int l, int iso, double s);
double inelasticity(int l, int iso, double s);
double conformal(double s, double s0);
double elastic_mom( double s, double sth);

complex<double> GKPRY_partial_wave(int l, int iso, double s);
complex<double> GKPRY_iso_amp(int iso, double s, double z);
complex<double> GKPRY_amplitude(double s, double z);
double GKPRY_cross_section(double s);

//partial-waves.cpp
double u_man(double s, double z);
double t_man(double s, double z);
double kallen(double s, double t, double u);
double legendre(int l, double x);
complex<double> VENEZ_iso_amp(int iso, double ** coup, double alph[], double s, double z);
complex<double> VENEZ_partial_wave(int l, int iso, double ** coup, double alph[], double s);
complex<double> VENEZ_amplitude(double ** coup, double alph[], double s, double z);
double VENEZ_cross_section(double ** coup, double alph[], double s);

//misc math stuff
cd cgamma(cd z);
void gauleg(double x1, double x2, double x[], double w[], int n);

//plotting
int WaveTranslate(string OPTION);
void plotGKPY(int CASE, int plot, string OPTION, string OUTPUT);
void getCOUPLING(string INPUT, double ** output);
void plotVENEZ(int MODE, int plot, string OPTION, string INPUT, string OUTPUT);
void plotFILE(string INPUT, string OUTPUT);

//fitting.cpp
int getDATA(string INPUT);
void fitVENEZ(int MODE, string OPTION, string INPUT, string OUTPUT);

extern double s_dat[], re_dat[], im_dat[];
extern int ISOCHOICE;
void get_Data();
