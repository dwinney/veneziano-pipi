// Contains functions related specifically to pi-pi scattering.
//
// Dependencies: venez-amp.cpp, cgamma.cpp
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@giu.edu
// ---------------------------------------------------------------------------

//TODO: Uncouple from big header.
#include  "veneziano.h"



//TODO: Figure out how to get rid of global variables.
//TODO: Add Error Message if n > maxN.

// Outputs: Isospin projected Veneziano amplitude for pi-pi scattering.
// Inputs:  iso - isospin (0, 1, or 2)
//          coup - couplings matrix a_n,i
//          alph - Regge trajectory parameters
//          s, t, u - Mandelstam Variables
cd isospin_amp(int iso, double coup[][maxN+1], double alph[], double s, double t, double u)
{
        double couplings[6];
        cd amp = 0.;
        double n_coup[maxN+1];
        for (int n = 1; n < maxN + 1; n++)
        {
                for (int i = 1; i < maxN+1; i++)
                {
                        n_coup[i] = coup[n][i];
                }
                switch (iso) {
                case 0: amp += n_amp(n, alph, n_coup, t, u) - 3. *(n_amp(n, alph, n_coup, s, u) + n_amp(n, alph, n_coup, s, t));
                        break;
                case 1: amp += 2.* (n_amp(n, alph, n_coup, s, u) - n_amp(n, alph, n_coup, s, t));
                        break;
                case 2: amp += -2.*n_amp(n, alph, n_coup, t, u);
                        break;
                default: cout << "Isospin " << iso << " not allowed (only 0, 1, or 2). Quiting..." << endl;
                        exit(1);
                }
        }

        return amp;
}

//TODO: Documentation.

double phase_shift(int l, int iso, double s)
{
        int wave = (3*l - iso)/2;
        if (iso > 2 || iso < 0 || l < 0) {wave = 600;}

        double mPi, mK, mEta, sthPi, sthK, sthEta;
        //Masses
        mPi = 0.1396;
        mK = 0.496;
        mEta = 0.54753;

        //Thresholds for pi, eta, and K
        sthPi = 4.*pow(mPi, 2.) - pow(1.l, -8.);
        sthK = 4.*pow(mK, 2.);
        sthEta = 4.*pow(mEta, 2.);

        switch (wave) {
        //S2 - wave (l = 0, iso = 2)
        case -1:
                double cot_delta;
                double temp1, temp2, temp3;
                double B0, B1, z2, wl, k;
                double sh, sl, sm;

                sl = pow(1.05, 2.);
                sm = pow(.850, 2.);
                sh = pow(1.42, 2.);

                z2 = .1435;
                k = elastic_mom(s, sthPi);
                wl = conformal(s, sl);

                B0 = -79.4;
                B1 = -63.0;

                temp1 = sqrt(s) / (2.*k);
                temp2 = pow(mPi, 2.) / (s - 2.*pow(z2,2.));

                if ((s > sthPi) && (s <= sm)) //Low energy parameterization
                {
                        temp3 = B0 + B1 * wl;
                        cot_delta = temp1 * temp2 * temp3;
                }
                else if ((s > sm) && ( s <= sh)) //Intermediate energies
                {
                        double Bh0, Bh1, Bh2, wh, whsm, wlsm;
                        double temp4;

                        wh = conformal(s, sh);
                        whsm = conformal(sm, sh);
                        wlsm = conformal(sm, sl);

                        Bh2 = 32.;
                        Bh0 = B0 + B1 * wlsm; // see wl above
                        temp4 = pow((sqrt(sm) + sqrt(sh - sm)) / (sqrt(sm) + sqrt(sl - sm)), 2.);
                        Bh1 = B1* (sl / sh) * (sqrt(sh - sm) / sqrt(sl - sm)) * temp4;

                        temp3 = Bh0 + Bh1*(wh - whsm) + Bh2*pow(wh - whsm, 2.);

                        cot_delta = temp1 * temp2 * temp3;

                }
                else{break;}
                return atan(1. / cot_delta);
                break;

        case 0:
                cout << "S0 wave >:O" << endl;
                break;

        default:
                cout << "Invalid Phase Shift with l = " << l << " and isospin = " << iso << ". Quiting..." << endl;
                exit(1);

        }
        return 0.;
}

//Returns conformal variable given Mandelstam s and branchign point s0.
double conformal(double s, double s0)
{
        double numerator = sqrt(s) - sqrt(s0 - s);
        double denominator = sqrt(s) + sqrt(s0 - s);
        return numerator/denominator;
}

double elastic_mom( double s, double sth)
{
        return sqrt(s - sth) / 2.;
}
