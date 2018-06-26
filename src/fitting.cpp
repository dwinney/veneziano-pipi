#include "pipi.h"
#include <TFitter.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TMath.h>

//TODO: COMMENTS

double s_dat[10000], re_dat[10000], im_dat[10000];
double regge[3];
int NPOINTS, wave, error;
int calls = 0;
int MODE = -600;

int getDATA(string INPUT)
{
        int NPOINTS = 0;
        ifstream data(INPUT.c_str());
        if (data.fail())
        {
                cout << " Could not open INPUT file. Quitting..." << endl;
                cout << endl;
                exit(1);
        }
        else
        {
                cout << " Reading data from " + INPUT << endl;
                cout << " NOTE: Energy should be in sqrt(s) (GeV)..." << endl;
                while (!data.eof())
                {
                        data >> s_dat[NPOINTS];
                        data >> re_dat[NPOINTS];
                        data >> im_dat[NPOINTS];
                        NPOINTS++;
                }
                cout << " Imported " << NPOINTS - 1 << " data points." << endl;
                cout << endl;
        }
        return NPOINTS;
}


double chi_square(double ** coup, double alph[])
{
        double chi = 0.;
        double s;
        complex<double> VENEZ;

        error = 10.;

        for (int i = 0; i < NPOINTS - 1; i++)
        {
                s = s_dat[i]*s_dat[i];

                if (MODE == 1)
                {
                        int iso = wave % 10;
                        int l = wave / 10;
                        VENEZ = VENEZ_partial_wave(l, iso, coup, alph, s);
                }
                else if (MODE == 2)
                {
                        VENEZ = VENEZ_iso_amp(wave, coup, alph, s, 1.);
                }
                else
                {
                        cout << " Invalid AMP selected. Qutting..." << endl;
                        cout << endl; exit(1);
                }


                chi += pow((im_dat[i] - imag(VENEZ)) / error, 2.);
                chi += pow((re_dat[i] - real(VENEZ)) / error, 2.);

        }

        calls++;
        if (calls % 50 == 0)
        {
                cout << "  ---"<< calls << " calls to chi-squared." << endl;
        }

        chi /= (NPOINTS-10);
        return chi;
}

static void WRAPPER(int& npar, double* g, double& result, double *par, int flag)
{
        double **a_matrix;
        a_matrix = new double *[maxN+1];
        for (int i = 0; i < maxN+1; i++) a_matrix[i] = new double[maxN+1];

        a_matrix[1][1] = par[0];
        a_matrix[2][1] = par[1];
        a_matrix[2][2] = par[2];
        a_matrix[3][1] = par[3];
        a_matrix[3][2] = par[4];
        a_matrix[3][3] = par[5];

        regge[0] = par[6];
        regge[1] = par[7];
        regge[2] = par[8];

        result = chi_square(a_matrix, regge);
        // cout << result << endl;
}

void fitVENEZ(int CASE, string OPTION, string INPUT, string OUTPUT)
{
        setprecision(17);
        NPOINTS = getDATA(INPUT);
        MODE = CASE;
        if (MODE == 1)
        {
                wave = WaveTranslate(OPTION);
        }
        else if (MODE == 2)
        {
                wave = atoi(OPTION.c_str());
        }

        cout << " Starting TMinuit..." << endl;
        cout << endl;
        TMinuit minuit(9);

        minuit.SetFCN(WRAPPER);

        double arglist[10];
        int ierflg = 0;
        minuit.mnparm(0, "a_1,1", 0, 0.01, 0., 0., ierflg);
        minuit.mnparm(1, "a_2,1", 0, 0.01, 0., 0., ierflg);
        minuit.mnparm(2, "a_2,2", 0, 0.01, 0., 0., ierflg);
        minuit.mnparm(3, "a_3,1", 0, 0.01, 0., 0., ierflg);
        minuit.mnparm(4, "a_3,2", 0, 0.01, 0., 0., ierflg);
        minuit.mnparm(5, "a_3,3", 0, 0.01, 0., 0., ierflg);

        minuit.mnparm(6, "intercept", 1., 0.1, 0, 0, ierflg);
        minuit.mnparm(7, "slope", .9, 0.1, 0, 0, ierflg);
        minuit.mnparm(8, "width", .18, 0.1, 0, 0, ierflg);

        arglist[0]=-1;
        minuit.mnexcm("SET PRINTOUT", arglist, 1, ierflg);
        minuit.mnexcm("SET NOW", arglist, 1, ierflg);
        arglist[0]=2;
        minuit.mnexcm("SET STR", arglist, 1, ierflg);
        minuit.SetErrorDef(1); //1 for chi square

        // if (CASE == 2) //ISOSPIN 1
        // {
        //         cout << endl;
        //         cout << " Fitting Isospin-1 amplidude." << endl;
        //         cout << " Using linear rho regge trajectory..." << endl;
        //
        //         arglist[0] = 100000;
        //         minuit.mnexcm("SIMPLEX", arglist,1,ierflg);
        //         minuit.mnexcm("MIGRAD", arglist, 1, ierflg);
        // }
        //
        if (OPTION == "P1")
        {
                cout << endl;
                cout << " For P-wave, only fitting spin-1 resonances." << endl;
                cout << " Fitting rho(770), rho(1450), and rho(1570)..." << endl;

                minuit.FixParameter(2);
                minuit.FixParameter(4);
                minuit.FixParameter(5);

                arglist[0] = 100000;
                minuit.mnexcm("SIMPLEX", arglist,1,ierflg);
                minuit.mnexcm("MIGRAD", arglist, 1, ierflg);

        }
        else if (OPTION == "F1")
        {
                cout << endl;
                cout << " For F-wave, only fitting spin-3 resonances." << endl;
                cout << " Fitting rho_3(1690)..." << endl;

                minuit.FixParameter(0);
                minuit.FixParameter(1);
                minuit.FixParameter(2);
                minuit.FixParameter(3);
                minuit.FixParameter(4);

                arglist[0] = 100000;
                minuit.mnexcm("SIMPLEX", arglist,1,ierflg);
                minuit.mnexcm("MIGRAD", arglist, 1, ierflg);
        }
        else if (OPTION == "S0")
        {
                cout << endl;
                cout << " For S-wave, fitting scalar resonances." << endl;
                cout << " Fitting f0(550), " << endl;

                minuit.FixParameter(0);
                minuit.FixParameter(1);
                minuit.FixParameter(3);
                minuit.FixParameter(5);
                arglist[0] = 100000;
                minuit.mnexcm("SIMPLEX", arglist,1,ierflg);
                minuit.mnexcm("MIGRAD", arglist, 1, ierflg);
        }

        cout << endl;
        cout << " Fit done." << endl;
        cout << " Fetching Best-fit values..." << endl;
        cout << endl;

        double coup[4][4], err;
        minuit.GetParameter(0, coup[1][1], err);
        minuit.GetParameter(1, coup[2][1], err);
        minuit.GetParameter(2, coup[2][2], err);
        minuit.GetParameter(3, coup[3][1], err);
        minuit.GetParameter(4, coup[3][2], err);
        minuit.GetParameter(5, coup[3][3], err);

        cout << " BEST TRAJECTORY COUPLINGS:" << endl;
        cout << setw(25) << "Intercept (a_0): " << setw(10) << regge[0] << endl;
        cout << setw(25) << "Slope (a^prime): " << setw(10) << regge[1] << endl;
        cout << setw(25) << "Rho width (Gamma): " << setw(10) << regge[2] << endl;
        cout << endl;

        cout << " BEST FIT COUPLINGS: " << endl;
        cout << setw(10) << "n" << setw(10) << "i" << setw(15) << "a_n,i" << endl;
        for (int i = 1; i < 4; i++)
        {
                for (int j = 1; j <= i; j++)
                {
                        cout << setw(10) << i << setw(10) << j << setw(15) << coup[i][j] << endl;
                }
        }
        cout << endl;

        OUTPUT += "-COUPLINGS.dat";
        cout << " Printing results.." << endl;

        ofstream couple;
        couple.open(OUTPUT.c_str());
        couple << left << setw(20) << regge[0] << setw(20) << regge[1] << setw(20) << regge[2] << endl;
        for (int i = 1; i < 4; i++)
        {
                for (int j = 1; j < i+1; j++)
                {
                        couple << left << setw(20) << i << setw(20) << j <<setw(20) << coup[i][j] << endl;
                }
        }
        couple.close();

        cout << " Output to " << OUTPUT << endl;
        cout << endl;

}
