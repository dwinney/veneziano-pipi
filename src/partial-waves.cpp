// Kinematic functions for extracting s-channel partial waves from an amplitude.
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@giu.edu
// ---------------------------------------------------------------------------

//TODO: Uncouple from big header.
#include "veneziano.h"

// Mandelstam variables t and u as functions of s and s-channel scattering angle z
double t_man(double s, double z)
{
        double psqr = (s - 4.*pow(xmpi,2.))/4.;
        double result = -2.*psqr* (1. - z);
        return result;
}

double u_man(double s, double z)
{
        double psqr = (s - 4.*pow(xmpi,2.))/4.;
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
cd partial_wave(int l, int iso, double coup[][maxN+1], double alph[], double s)
{
        double Pl, z, t, u;
        double weights[INTP], abscissas[INTP];
        cd sum;

        //prepares arrays of weights and abscissas, using Legendre functions
        gauleg(-1., 1., abscissas, weights, INTP);

        //integrates by Gaussian quadrature
        for (int i = 0; i < INTP; i++)
        {
                z = abscissas[i];
                t = t_man(s, z);
                u = u_man(s, z);

                Pl = legendre(l, z);
                sum += weights[i]*Pl*isospin_amp(iso, coup, alph, s,t,u);
        }
        return .5*sum;
}

double cross_section(int iso, double coup[][maxN+1], double alph[], double s, double z)
{
        double t, u;
        t = t_man(s, z);
        u =u_man(s, z);
        cd amp = isospin_amp(iso, coup, alph, s, t, u);
        cd result = amp*conj(amp);
        return real(result);
}
