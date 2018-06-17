
// Contains functions related to VENEZ model applied to
// pi-pi scattering.
//
// Dependencies: venez-amp.cpp, cgamma.cpp
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

//TODO: Uncouple from big header.
#include  "pipi.h"

double mres;

//Real Linear Regge Trajectory
// a0 - intercept (~ .5 for rho)
// ap - slope parameter (~ .09 for rho)
cd rtraj(double alph[], double s)
{
        cd result = alph[0] + alph[1]*s;
        return result;
}

//TODO: Remove the global variables for gamma_res and m_res
//Complex Regge Trajectory (w/ Phase Space factor)
cd ctraj(double alph[], double s)
{
        cd phase_space = sqrt(cd(s - 4.*pow(mPi, 2.)));
        cd imagpart = xi * alph[2] * mres;
        cd result = alph[0] + alph[1] *(s - mres*mres)  + imagpart * phase_space;
        // cout << result << endl;
        return result;
}

//Veneziano-type function of two invariant variables, s and t.
// n - number of poles
// alph[] - array with alpha_0 and alpha^prime
// coupling[] - array of size n + 1 with couplings a_n,i for fixed n
cd n_amp(int n, double alph[], double coupling[], double s, double t)
{
        cd s_alpha = rtraj(alph, s);
        cd t_alpha = rtraj(alph, t);
        cd cs_alpha = ctraj(alph, s);
        cd ct_alpha = ctraj(alph, t);
        cd nn = double(n)*xr;

        //Prefactor containing the sum of two poles one in each variable
        cd coeff = (2.*nn - s_alpha - t_alpha)
                   / ((nn - cs_alpha) * (nn - ct_alpha));
        // cd coeff = 1.;

        //Background factor that introduces power-surpressed poles at high energy.
        //Reproduces correct asymptotic behavior.
        cd background = (cgamma( backN + 1. - s_alpha) * cgamma( backN + 1. - t_alpha)) / (cgamma(1. + backN - nn) * cgamma( 1. - cs_alpha - ct_alpha + backN + nn ));
        // cd background = 1.;

        //Sum over couplings to resonances on trajectory specified by n
        cd sum, ii;
        for (int i = 1; i < n + 1; i++)
        {
                ii = double(i)*xr- 1.;
                sum += coupling[i] * pow(-(s_alpha + t_alpha), ii);
        }

        cd result = coeff*background*sum;
        return result;
}

//TODO: Figure out how to get rid of global variables.
//TODO: Add Error Message if n > maxN.

// Outputs: Isospin projected Veneziano amplitude for pi-pi scattering.
// Inputs:  iso - isospin (0, 1, or 2)
//          coup - couplings matrix a_n,i
//          alph - Regge trajectory parameters
//          s, t, u - Mandelstam Variables
complex<double> VENEZ_iso_amp(int iso, double ** coup, double alph[], double s, double z)
{
        double n_coup[maxN+1];
        complex<double> amp = 0.;
        double t, u;

        t = t_man(s, z); u = u_man(s, z);
        for (int n = 1; n < maxN+1; n++)
        {
                for (int i = 1; i < maxN+1; i++)
                {
                        n_coup[i] = coup[n][i];
                }
                // cout << n_coup[1] << " " << n_coup[2] << " " << n_coup[3] << endl;
                switch (iso) {
                case 0: mres = .500; //f0(500)
                        amp += n_amp(n, alph, n_coup, t, u) - 3. *(n_amp(n, alph, n_coup, s, u) + n_amp(n, alph, n_coup, s, t));
                        break;
                case 1: mres = 0.770; //rho
                        amp += 2.* (n_amp(n, alph, n_coup, s, u) - n_amp(n, alph, n_coup, s, t));
                        break;
                case 2: mres = 0.770; //rho
                        amp += -2. * n_amp(n, alph, n_coup, t, u);
                        break;
                default: cout << "Isospin " << iso << " not allowed (only 0, 1, or 2). Quiting..." << endl;
                        exit(1);
                }
        }

        return amp;
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

//TODO: remove dependency on function???? (i.e. isospin_amp)
// Outputs the complex partial-wave amplitude for pi-pi scattering at fixed s.
// l - partial wave
// iso - isospin projection ( iso = 1, 2, 3)
// coup - a_n,i matrix of couplings
// alph - regge trajectory parameters
complex<double> VENEZ_partial_wave(int l, int iso, double ** coup, double alph[], double s)
{
        double Pl, z;
        double weights[INTP], abscissas[INTP];
        cd sum;

        //prepares arrays of weights and abscissas, using Legendre functions
        gauleg(-1., 1., abscissas, weights, INTP);

        //integrates by Gaussian quadrature
        for (int i = 0; i < INTP; i++)
        {
                z = abscissas[i];
                Pl = legendre(l, z);
                sum += weights[i] * Pl * VENEZ_iso_amp(iso, coup, alph, s, z);
        }
        return .5*sum;
}

complex<double> VENEZ_amplitude(double ** coup, double alph[], double s, double z)
{
        double clebcsh[3] = {.33333, .5, 1.6666};
        complex<double> amp;

        for (int i = 0; i < 3; i++)
        {
                amp += clebcsh[i] * VENEZ_iso_amp(i, coup, alph, s, z);
        }
        return amp;
}

double VENEZ_cross_section(double ** coup, double alph[], double s)
{
        double k = elastic_mom(s, sthPi);
        complex<double> amp = VENEZ_amplitude(coup, alph, s, 1.);

        double sigma = imag(amp) / (2. * sqrt(s) * k );
        return sigma;
}
