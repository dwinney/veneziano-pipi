
// Contains functions related to VENEZ model applied to
// pi-pi scattering.
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "pipi.h"

double mres;

//Real Linear Regge Trajectory
// a0 - intercept (~ .5 for rho)
// ap - slope parameter (~ .9 for rho)
complex<double> rtraj(double alph[], double s)
{
        complex<double> result = alph[0] + alph[1] * s;
        return result;
}

//TODO: Remove the global variables for gamma_res and m_res
//Complex Regge Trajectory (w/ Phase Space factor)
complex<double> ctraj(double alph[], double s)
{
        complex<double> phase_space = sqrt( (complex<double>)(s - sthPi));
        complex<double> imagpart = xi * 1.* alph[2] * mres * alph[1];
        // cout << phase_space << endl;
        complex<double> result = alph[0] + alph[1] *(s - mres*mres)  + imagpart * phase_space;
        // cout << result << endl;
        return result;
}

//Veneziano-type function of two invariant variables, s and t.
// n - number of poles
// alph[] - array with alpha_0 and alpha^prime
// coupling[] - array of size n + 1 with couplings a_n,i for fixed n
complex<double> n_amp(int n, double alph[], double coupling[], double s, double t)
{
        complex<double> s_alpha = rtraj(alph, s);
        complex<double> t_alpha = rtraj(alph, t);
        complex<double> cs_alpha = ctraj(alph, s);
        complex<double> ct_alpha = ctraj(alph, t);
        complex<double> nn = double(n)*xr;

        //Prefactor containing the sum of two poles one in each variable
        complex<double> coeff = (2.*nn - s_alpha - t_alpha)
                                / ((nn - cs_alpha) * (nn - ct_alpha));
        // complex<double> coeff = 1.;

        //Background factor that introduces power-surpressed poles at high energy.
        //Reproduces correct asymptotic behavior.
        complex<double> background = (cgamma( backN + 1. - s_alpha) * cgamma( backN + 1. - t_alpha)) / (cgamma(1. + backN - nn) * cgamma( 1. - cs_alpha - ct_alpha + backN + nn ));
        // complex<double> background = 1.;

        //Sum over couplings to resonances on trajectory specified by n
        complex<double> sum, ii;
        for (int i = 1; i < n + 1; i++)
        {
                ii = double(i)*xr- 1.;
                sum += coupling[i] * pow(-(s_alpha + t_alpha), ii);
        }

        complex<double> result = coeff*background*sum;
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



//TODO: remove dependency on function???? (i.e. isospin_amp)
// Outputs the complex partial-wave amplitude for pi-pi scattering at fixed s.
// l - partial wave
// iso - isospin projection ( iso = 0, 1, 2)
// coup - a_n,i matrix of couplings
// alph - regge trajectory parameters
complex<double> VENEZ_partial_wave(int l, int iso, double ** coup, double alph[], double s)
{
        double Pl, z, re, im;
        double weights[INTP], abscissas[INTP];
        complex<double> venez;

        //prepares arrays of weights and abscissas, using Legendre functions
        gauleg(-1., 1., abscissas, weights, INTP);

        //integrates by Gaussian quadrature
        re = 0.; im = 0.;
        for (int i = 0; i < INTP; i++)
        {
                z = abscissas[i];
                Pl = legendre(l, z);
                venez = VENEZ_iso_amp(iso, coup, alph, s, z);
                re += weights[i] * Pl * real(venez);
                im += weights[i] * Pl * imag(venez);
        }
        complex<double> amp(re, im);
        return .5*amp;
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
