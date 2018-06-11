#include "veneziano.h"

#include <TFitter.h>
#include <TRandom.h>
#include <TROOT.h>
#include <TMath.h>

// Apparently TMinuit needs a pointer here
double data_real[1938];

// Function to minimize
double chi_square(double coup[][maxN+1], double alph[])
{
        double chi2 = 0.;
        double GKPRY, VENEZ, error;
        double s, smax, step;
        int nstep;

        smax = pow(1.42, 2.);
        step = .1;
        nstep = (int) (smax - sthPi) / step;

        error = 10.;

        for (int i = 0; i < nstep; i++)
        {
                s = sthPi + 1e-5 + double(i)*step;

                GKPRY = data_real[i];
                VENEZ = real(VENEZ_isospin_amp(1, coup, alph, s, 1. ));
                chi2 += pow(((GKPRY - VENEZ) / error), 2.);
        }
        // cout << chi2 << endl;
        return chi2;

}

void WRAPPER(int& npar, double* g, double& result, double par[], int flag)
{
        double coup[4][4];
        coup[1][1] = par[3];
        // coup[2][1] = par[3];
        // coup[3][1] = par[4];
        // coup[3][3] = par[5];

        double alph [3];
        alph[1] = par[0];
        alph[2] = par[1];
        alph[3] = par[2];

        result = chi_square(coup, alph);
        // cout << par[0] << endl;
}

int main(int argc, char* argv[])
{
        double s, smax, step;
        int nstep;

        smax = pow(1.42, 2.);
        step = .1;
        nstep = (int) (smax - sthPi) / step;

        for (int i = 0; i < nstep; i++)
        {
                s = sthPi + 1e-5 + step * double(i);
                data_real[i] = real(GKPRY_iso_amp(1, s, 1.));
                cout << data_real[i];
        }


        // Fitting
        TFitter minuit(3);

        minuit.SetFCN(WRAPPER);

        double p1= 1;
        minuit.ExecuteCommand("SET PRINTOUT", &p1, 1);

        minuit.SetParameter(0, "alph_0", 0.3, .01, 0.1, 2.);
        minuit.SetParameter(1, "alph^p", 1., .01, .5, 1.5);
        minuit.SetParameter(2, "alph_im", 0.1, .01, .01, .3);
        minuit.SetParameter(3, "a11", -120., 0.01, 0., 0.);
        // minuit.SetParameter(3, "a21", 0., 1, 0., 0.);
        // minuit.SetParameter(4, "a31", 0., 1, 0., 0.);
        // minuit.SetParameter(5, "a33", 0., 1, 0., 0.);

        minuit.ExecuteCommand("SIMPLEX", 0, 0);
        // minuit.ExecuteCommand("MIGRAD", 0, 0);

        double a[4][4], alph[2];
        for (int c = 0; c < 4; c++)
        {
                for(int b = 0; b < 4; b++)
                {
                        a[c][b] = 0.;
                }
        }
        a[1][1] = minuit.GetParameter(3);
        // cout << a[1][1] << endl;
        // a[2][1] = minuit.GetParameter(3);
        // a[3][1] = minuit.GetParameter(4);
        // a[3][3] = minuit.GetParameter(5);
        //
        alph[0] = minuit.GetParameter(0);
        alph[1] = minuit.GetParameter(1);
        alph[2] = minuit.GetParameter(2);

        cd amp, gk;

        ofstream pwave;
        pwave.open("./output/rho.dat");
        pwave << left << setw(30) << "#s" << setw(30) << "Re[Amp]" << setw(30) << "Im[Amp]" << endl;
        for (int i = 0; i < 201; i++)
        {
                s = sthPi + .01*double(i);
                // amp = data_real[i];
                amp = VENEZ_isospin_amp(1, a, alph, s, 1. );
                gk = GKPRY_iso_amp(1, s, 1.);

                pwave << left << setw(30) << s << setw(30) << real(amp)  << setw(30) << imag(amp) << setw(30) << real(gk) << setw(30) << imag(gk) << endl;
        }
        pwave.close();
        system("gnuplot ./src/gnuplot/graph.gnu");
        system("okular ./output/rho.pdf");

        return 0;
}
