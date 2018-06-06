#include <iostream>
#include <fstream>
#include <math.h>
#include <cmath>
#include <complex>
#include <iomanip>
#include <cstdlib>
#include <stdio.h>
#include <string.h>
#include <sstream>
#include <stdlib.h>

using namespace std;
typedef complex<double> cd;

const double pi = 3.1415926535897932384626433832795028841972;
const double gamma_res = .145; //Width of resonance, this case the rho
const double m_res = .77545; //Mass of resonance (rho)
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

#define backN   20.*xr  //Background N, just needs to be "large enough"
#define maxN    6       //Truncated maximum n
#define INTP    500     //Number of points for numerical integration

//amp.cpp
cd rtraj(double a0, double ap, double s);
cd ctraj(double a0, double ap, double s);
cd n_amp(int n, double alph[], double coupling[], double s, double t);

//pipi-amp.cpp
cd isospin_amp(int iso, double coup[][maxN+1], double alph[], double s, double t, double u);
double phase_shift(int l, int iso, double s);
double inelasticity(int l, int iso, double s);
double conformal(double s, double s0);
double elastic_mom( double s, double sth);

//partial-waves.cpp
double u_man(double s, double z);
double t_man(double s, double z);
double kallen(double s, double t, double u);
double legendre(int l, double x);
cd partial_wave(int l, int iso, double coup[][maxN+1], double alph[], double s);
double cross_section(int iso, double coup[][maxN+1], double alph[], double s, double z);

//misc math stuff
cd cgamma(cd z);
void gauleg(double x1, double x2, double x[], double w[], int n);
