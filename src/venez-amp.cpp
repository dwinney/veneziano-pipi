// Contains functions related to the a general amplitude in Veneziano form.
//
// Dependencies: cgamma.cpp
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

//TODO: Uncouple from big header.
#include "veneziano.h"

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
