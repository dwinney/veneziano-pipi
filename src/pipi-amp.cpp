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
        //Error Checks :)
        int wave = (3*l - iso)/2;
        if (iso > 2 || iso < 0) {wave = 600;}
        if (l % 2 != iso % 2) {wave = 600;}
        if (l > 4)
        {cout << "Contributions of l > 3 are completely negligable. Quitting..." << endl;;
         exit (1);}

        double cot_delta, delta, sh;
        double temp1, temp2, temp3;
        sh = pow(1.42, 2.);

        //Momenta
        double k, k2;
        k = elastic_mom(s, sthPi);
        k2 = elastic_mom(s, sthK);

//S2 wave (l = 0, iso = 2)
        if (wave == -1)
        {
                double B0, B1, z2, wl;
                double sl, sm;

                sl = pow(1.05, 2.);
                sm = pow(.850, 2.);

                wl = conformal(s, sl);

                z2 = .1435;
                B0 = -79.4;
                B1 = -63.0;

                temp1 = sqrt(s) / (2.*k);
                temp2 = mPi*mPi / (s - 2.*z2*z2);

                if ((s > sthPi) && (s <= sm)) //Low energy parameterization
                {
                        temp3 = B0 + B1 * wl;
                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else if ((s > sm) && ( s <= sh)) //Intermediate energies
                {
                        double Bh0, Bh1, Bh2, wh, whsm, wlsm;
                        double temp4;

                        wh = conformal(s, sh);
                        whsm = conformal(sm, sh);
                        wlsm = conformal(sm, sl);

                        Bh2 = 32.;
                        Bh0 = B0 + B1 * wlsm;
                        temp4 = pow((sqrt(sm) + sqrt(sh - sm)) / (sqrt(sm) + sqrt(sl - sm)), 2.);
                        Bh1 = B1 * (sl / sh) * (sqrt(sh - sm) / sqrt(sl - sm)) * temp4;

                        temp3 = Bh0 + Bh1*(wh - whsm) + Bh2*pow(wh - whsm, 2.);

                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else {delta = 0.;}

                return delta;
        }

//S0 wave (l = 0, iso = 0)
        if (wave == 0)
        {
                cout << "S0 wave >:O" << endl;
        }

//P1 wave (l = 1, iso = 1)
        else if (wave == 1)
        {
                double s0, w;

                s0 = pow(1.05, 2.);
                w = conformal(s, s0);

                if ((s > sthPi) && (s < sthK))
                {
                        double B0, B1;
                        B0 = 1.043;
                        B1 = .19;

                        temp1 = sqrt(s)*(mRho*mRho - s)/(2. * pow(k, 3.));
                        temp2 = 2.*pow(mPi, 3.)/(mRho*mRho*sqrt(s));

                        cot_delta = temp1*(temp2 + B0 + B1* w);
                        delta = atan2(1., cot_delta);
                }
                else if ((s >= sthK) && (s < sh))
                {
                        double lambda0, wK, KK;
                        double lambda1, lambda2;
                        lambda1 = 1.38;
                        lambda2 = -1.70;

                        //lambda0 = low energy (see above) at KK threshold
                        double B0, B1;
                        B0 = 1.043;
                        B1 = .19;
                        KK = elastic_mom(sthK, sthPi);
                        wK = conformal(sthK, s0);
                        temp1 = sqrt(sthK)*(mRho*mRho - s)/(2. * pow(KK, 3.));
                        temp2 = 2.*pow(mPi, 3.)/(mRho*mRho*sqrt(sthK));
                        lambda0 = temp1*(temp2 + B0 + B1* wK);

                        temp3 = (sqrt(s) / (2.*mK)) - 1.;

                        cot_delta = lambda0 + lambda1*temp3 + lambda2*temp3*temp3;
                        delta = atan2(1., cot_delta);
                }
                else {delta = 0.;}

                return delta;
        }

//D2 wave (l = 2, iso = 2)
        else if (wave == 2)
        {
                if ((s > sthPi) && (s < sh))
                {
                        double w, s0, Del;
                        double B0, B1, B2;

                        s0 = pow(1.45,2.);
                        B0 = 4.1e3;
                        B1 = 8.6e3;
                        B2 = 25.5e3;

                        w = conformal(s, s0);

                        temp1 = sqrt(2) / (2.* pow(k, 5.));
                        temp2 = (pow(mPi, 4.) * s) / (sthPi + 4.*Del*Del -s);
                        temp3 = B0 + B1 * w + B2 * w * w;

                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else {delta = 0.;}
                return delta;

        }

//D0 wave (l = 2, iso = 0)
        else if (wave == 3)
        {
                double B0, B1;
                double s0 = pow(1.05, 2.);

                B0 = 12.40;
                B1 = 10.06;

                temp1 = sqrt(s) / (2.*pow(k, 5.));
                temp2 = (mF2*mF2 - s)*mPi*mPi;

                if ((s > sthPi) && (s <= sthK))
                {
                        double w = conformal(s, s0);
                        temp3 = B0 + B1* w;
                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else if ((s > sthK) && (s < sh))
                {
                        double wh, s02, B0h, B1h;

                        B1h = 43.2;
                        s02 = pow(1.45, 2.);
                        wh = conformal(s, s02);
                        B0h = 18.69;

                        temp3 = B0h + B1h * wh;
                        cot_delta = temp1 * temp2 * temp3;
                        delta = atan2(1., cot_delta);
                }
                else {delta = 0.;}
        }

//Anything else
        else
        {
                cout << "Invalid Phase Shift with l = " << l << " and isospin = " << iso << ". Quiting..." << endl;
                exit(1);

        }
}

double inelasticity(int l, int iso, double s)
{
        double eta;
        double sh = pow(1.42, 2.);

        int wave = (3*l - iso)/2;
        if (iso > 2 || iso < 0 || l < 0) {wave = 600;}
        if (l % 2 != iso % 2) {wave = 600;}
        if (l > 3) {return eta = 1.;} //Inelasticities for high partial waves negligable

//S2 wave (l = 0, iso = 2)
        if (wave == -1)
        {
                double si = pow(1.05, 2.);

                if ((s > sthPi) && (s <= si))          //Low energy parameterization
                {
                        eta = 1;
                }
                else if ((s > si) && ( s <= sh))        //Intermediate energies
                {
                        double eps = .28;

                        eta = 1. - eps*pow(1 - (si/s), 1.5);
                }
                return eta;
        }

//S0 wave (l = 0, iso = 0)
        else if (wave == 0)
        {
                cout << "S0 wave >:O" << endl;
        }

//P1 wave (l = 1, iso = 1)
        else if (wave == 1)
        {
                if ((s > sthPi) && (s < sthK))
                {
                        eta = 1.;
                }
                else if ((s > sthK) && (s < sh))
                {
                        double temp;
                        double ep1, ep2;
                        ep1 = 0.;
                        ep2 = .07;

                        temp = 1. - (sthK/s);
                        eta = 1 - ep1 * sqrt(temp) - ep2 * temp;
                }
                return eta;
        }

//D2 wave (l = 2, iso = 2)
        else if (wave == 2)
        {
                double si = pow(1.05, 2.);
                if ((s > sthPi) && (s <= si))
                {
                        eta = 1.;
                }
                else if ((s > si) && (s < sh))
                {
                        double ep = 0.;
                        eta = 1 - ep* pow(1 - (si/s),3.);
                }
                return eta;
        }

//D0 wave (l = 1, iso = 0)
        else if (wave == 3)
        {
                if ((s > sthPi) && (s <= sthK))
                {
                        eta = 1.;
                }
                if ((s > sthK) && (s < sh))
                {
                        double eps, r, temp1, temp2;
                        eps = 0.254;
                        r = 2.29;

                        double k2 = elastic_mom(s, sthK);

                        temp1 = (1. - sthK/s) / (1. - sthK/(mF2*mF2));
                        temp2 = 1. - k2 / elastic_mom(mF2*mF2, sthK);

                        eta = 1. - eps*pow(temp1, 2.5)*(1. + r*temp2);
                }
                return eta;
        }


//Anything else
        else
        {
                cout << "Invalid Inelasticity with l = " << l << " and isospin = " << iso << ". Quiting..." << endl;
                exit(1);
        }
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
