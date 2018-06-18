#include "pipi.h"

//Converts string of spectroscopic name of partial wave to an INT
int WaveTranslate(string OPTION)
{
        int wave;
        if (OPTION == "S0") {wave = 0;}
        else if (OPTION == "S2") {wave = 02;}
        else if (OPTION == "P1") {wave = 11;}
        else if (OPTION == "D0") {wave = 20;}
        else if (OPTION == "D2") {wave = 22;}
        else if (OPTION == "F1") {wave = 31;}
        else {cout << "Invalid Partial Wave " << OPTION << ". Quitting..." << endl; cout << endl; exit(1);}
        return wave;
}

// Mandelstam variables t and u as functions of s and s-channel scattering angle z
double t_man(double s, double z)
{
        double psqr = (s - 4.*pow(mPi,2.))/4.;
        double result = -2.*psqr* (1. - z);
        return result;
}

double u_man(double s, double z)
{
        double psqr = (s - 4.*pow(mPi,2.))/4.;
        double result = -2.*psqr* (1. + z);
        return result;
}

// Kallen triangle function
double kallen(double s, double t, double u)
{
        double result = pow(s,2.) + pow(t, 2.) + pow(u, 2.) + 2.*s*t + 2.*t*s + 2.*t*u;
        return result;
}

// Elastic momentum above threshold sth, as a function of s.
double elastic_mom( double s, double sth)
{
        return sqrt(s - sth) / 2.;
}

// Legendre Polynomials P_l(x) (up to l = 5)
double legendre(int l, double x)
{
        double PL;
        switch (l)
        {
        case 0: PL = 1.; break;
        case 1: PL = x; break;
        case 2: PL = .5 * (3.*pow(x, 2.) - 1.); break;
        case 3: PL = .5 * (5.*pow(x,3.) - 3.*x); break;
        case 4: PL = (35.*pow(x,4.) - 30.*pow(x,2.) + 3.)/8.; break;
        case 5: PL = (63.*pow(x,5.) - 70.*pow(x,3.) + 15.*x)/8.; break;
        default: cout << "Legendre Polynomials don't go that high!!!" << endl;
                exit(1);
        }
        return PL;
}

//  cgamma.cpp -- Complex gamma function.
//      Algorithms and coefficient values from "Computation of Special
//      Functions", Zhang and Jin, John Wiley and Sons, 1996.
//
//  (C) 2003, C. Bond. All rights reserved.
//
//  Returns gamma function for complex argument 'z'.
//
//  Returns (1e308,0) if the real part of the argument is a negative integer
//  or 0 or exceeds 171.
//

complex<double> cgamma(complex<double> z)
{
        complex<double> g,z0,z1;
        double x0,q1,q2,x,y,th,th1,th2,g0,gr,gi,gr1,gi1;
        double na,t,x1,y1,sr,si;
        int i,j,k;

        static double a[] = {
                8.333333333333333e-02,
                -2.777777777777778e-03,
                7.936507936507937e-04,
                -5.952380952380952e-04,
                8.417508417508418e-04,
                -1.917526917526918e-03,
                6.410256410256410e-03,
                -2.955065359477124e-02,
                1.796443723688307e-01,
                -1.39243221690590
        };

        x = real(z);
        y = imag(z);
        if (x > 171) return complex<double>(1e308,0);
        if ((y == 0.0) && (x == (int)x) && (x <= 0.0))
                return complex<double>(1e308,0);
        else if (x < 0.0) {
                x1 = x;
                y1 = y;
                x = -x;
                y = -y;
        }
        x0 = x;
        if (x <= 7.0) {
                na = (int)(7.0-x);
                x0 = x+na;
        }
        q1 = sqrt(x0*x0+y*y);
        th = atan(y/x0);
        gr = (x0-0.5)*log(q1)-th*y-x0+0.5*log(2.0*M_PI);
        gi = th*(x0-0.5)+y*log(q1)-y;
        for (k=0; k<10; k++) {
                t = pow(q1,-1.0-2.0*k);
                gr += (a[k]*t*cos((2.0*k+1.0)*th));
                gi -= (a[k]*t*sin((2.0*k+1.0)*th));
        }
        if (x <= 7.0) {
                gr1 = 0.0;
                gi1 = 0.0;
                for (j=0; j<na; j++) {
                        gr1 += (0.5*log((x+j)*(x+j)+y*y));
                        gi1 += atan(y/(x+j));
                }
                gr -= gr1;
                gi -= gi1;
        }
        // if (x1 <= 0.0) {
        //         q1 = sqrt(x*x+y*y);
        //         th1 = atan(y/x);
        //         sr = -sin(M_PI*x)*cosh(M_PI*y);
        //         si = -cos(M_PI*x)*sinh(M_PI*y);
        //         q2 = sqrt(sr*sr+si*si);
        //         th2 = atan(si/sr);
        //         if (sr < 0.0) th2 += M_PI;
        //         gr = log(M_PI/(q1*q2))-gr;
        //         gi = -th1-th2-gi;
        //         x = x1;
        //         y = y1;
        // }

        g0 = exp(gr);
        gr = g0*cos(gi);
        gi = g0*sin(gi);

        g = complex<double>(gr,gi);
        return g;
}

/*******************************************************************************
   Gauss-Legendre integration function, gauleg, from "Numerical Recipes in C"
   (Cambridge Univ. Press) by W.H. Press, S.A. Teukolsky, W.T. Vetterling, and
   B.P. Flannery
*******************************************************************************/


#define NR_END 1
#define EPS 3.0e-11 /* EPS is the relative precision. */

double *dvector(int nl, int nh)
/* allocate a double vector with subscript range v[nl..nh] */
{
        double *v;
        v=(double *)malloc((size_t) ((nh-nl+1+NR_END)*sizeof(double)));
        return v-nl+NR_END;
}

/******************************************************************************/
void gauleg(double x1, double x2, double x[], double w[], int n)
/*******************************************************************************
   Given the lower and upper limits of integration x1 and x2, and given n, this
   routine returns arrays x[1..n] and w[1..n] of length n, containing the abscissas
   and weights of the Gauss-Legendre n-point quadrature formula.
*******************************************************************************/
{
        int m,j,i;
        double z1,z,xm,xl,pp,p3,p2,p1;
        m=(n+1)/2; /* The roots are symmetric, so we only find half of them. */
        xm=0.5*(x2+x1);
        xl=0.5*(x2-x1);
        for (i=1; i<=m; i++) { /* Loop over the desired roots. */
                z=cos(3.141592654*(i-0.25)/(n+0.5));
                /* Starting with the above approximation to the ith root, we enter */
                /* the main loop of refinement by Newton's method.                 */
                do {
                        p1=1.0;
                        p2=0.0;
                        for (j=1; j<=n; j++) { /* Recurrence to get Legendre polynomial. */
                                p3=p2;
                                p2=p1;
                                p1=((2.0*j-1.0)*z*p2-(j-1.0)*p3)/j;
                        }
                        /* p1 is now the desired Legendre polynomial. We next compute */
                        /* pp, its derivative, by a standard relation involving also  */
                        /* p2, the polynomial of one lower order.                     */
                        pp=n*(z*p1-p2)/(z*z-1.0);
                        z1=z;
                        z=z1-p1/pp; /* Newton's method. */
                } while (fabs(z-z1) > EPS);
                x[i]=xm-xl*z; /* Scale the root to the desired interval, */
                x[n+1-i]=xm+xl*z; /* and put in its symmetric counterpart.   */
                w[i]=2.0*xl/((1.0-z*z)*pp*pp); /* Compute the weight             */
                w[n+1-i]=w[i];     /* and its symmetric counterpart. */
        }
}
