#include "pipi.h"
#include <TFitter.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TMath.h>

double s_dat[10000], re_dat[10000], im_dat[10000];
int NPOINTS, wave;
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
        double s, error;
        complex<double> VENEZ;

        error = 1.0;

        for (int i = 0; i < NPOINTS - 1; i++)
        {
                s = s_dat[i];
                complex<double> data(re_dat[i], im_dat[i]);
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
                chi += real((data - VENEZ) * conj(data - VENEZ)) / (error*error);
        }
        // cout << setw(20) << chi << endl;
        return chi;
}

static void WRAPPER(int& npar, double* g, double& result, double *par, int flag)
{
        double **a_matrix;
        a_matrix = new double *[maxN+1];
        for (int i = 0; i < maxN+1; i++) a_matrix[i] = new double[maxN+1];

        a_matrix[1][1] = par[3];
        a_matrix[2][1] = par[4];
        a_matrix[2][2] = par[5];
        a_matrix[3][1] = par[6];
        a_matrix[3][2] = par[7];
        a_matrix[3][3] = par[8];

        double alph [3];
        alph[0] = par[0];
        alph[1] = par[1];
        alph[2] = par[2];

        result = chi_square(a_matrix, alph);
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

        double start[3], stepsiz[4], min[4], max[4];

        //alpha_0
        start[0] = 1.; stepsiz[0] = 0.01;
        min[0] = 0.; max[0] = 0.;
        //alpha^p
        start[1] = 0.928524; stepsiz[1] = 0.1;
        min[1] = 0.; max[1] = 0.;
        //width
        start[2] = .145; stepsiz[2] = .01;
        min[2] = 0.0; max[2] = 1.;
        //couplings
        start[3] = 0.; stepsiz[3] = .2;
        min[3] = 0.; max[3] = 0.;

        double arglist[10];
        int ierflg = 0;
        minuit.mnparm(0, "alpha_0", start[0], stepsiz[0], min[0], max[0], ierflg);
        minuit.mnparm(1, "alpha_p", start[1], stepsiz[1], min[1], max[1], ierflg);
        minuit.mnparm(2, "width", start[2], stepsiz[2], min[2], max[2], ierflg);
        minuit.mnparm(3, "a_1,1", 5., stepsiz[3], 0, 0, ierflg);
        minuit.mnparm(4, "a_2,1", 0., stepsiz[3], -1., 1., ierflg);
        minuit.mnparm(5, "a_2,2", start[3], stepsiz[3], min[3], max[3], ierflg);
        minuit.mnparm(6, "a_3,1", 0, stepsiz[3], 0, 0, ierflg);
        minuit.mnparm(7, "a_3,2", start[3], stepsiz[3], min[3], max[3], ierflg);
        minuit.mnparm(8, "a_3,3", start[3], stepsiz[3], min[3], max[3], ierflg);

        arglist[0]=-1;
        minuit.mnexcm("SET PRINTOUT", arglist, 1, ierflg);
        minuit.mnexcm("SET NOW", arglist, 1, ierflg);
        arglist[0]=2;
        minuit.mnexcm("SET STR", arglist, 1, ierflg);
        minuit.SetErrorDef(1); //1 for chi square

        // minuit.FixParameter(0);
        // minuit.FixParameter(1);
        // minuit.FixParameter(2);

        //TODO: Impose Bose symmetry
        if (OPTION == "P1")
        {
                cout << endl;
                cout << " Using rho(770) trajectory... " << endl;
                // minuit.FixParameter(0);
                // minuit.FixParameter(1);
                // minuit.FixParameter(2);
                cout << " Only fitting spin-1 resonances: a_1,1 a_2,1 and a3,1..." << endl;
                cout << endl;
                minuit.FixParameter(4);
                minuit.FixParameter(5);
                minuit.FixParameter(7);
                // minuit.FixParameter(6);
                minuit.FixParameter(8);
        }
        else if (OPTION == "F1")
        {
                cout << endl;
                cout << " Using rho(770) trajectory... " << endl;
                minuit.FixParameter(0);
                minuit.FixParameter(1);
                minuit.FixParameter(2);
                cout << "Only fitting spin-3 resonances: a3,3..." << endl;
                cout << endl;
                minuit.FixParameter(3);
                minuit.FixParameter(4);
                minuit.FixParameter(5);
                minuit.FixParameter(6);
                minuit.FixParameter(7);
        }
        else if (OPTION == "S0")
        {
                cout <<"I dunno wht im doing..." << endl;
        }
        // minuit.FixParameter(0);
        // minuit.FixParameter(1);
        // minuit.FixParameter(2);
        // minuit.FixParameter(3);
        // minuit.FixParameter(4);
        // minuit.FixParameter(5);
        // minuit.FixParameter(6);
        // minuit.FixParameter(7);
        // minuit.FixParameter(8);

        arglist[0] = 1500;
        minuit.mnexcm("SIMPLEX", arglist,1,ierflg);
        // minuit.mnexcm("MINIMIZE", arglist, 1, ierflg);
        minuit.mnexcm("MIGRAD", arglist, 1, ierflg);
        // minuit.mnexcm("MINOS", arglist, 1, ierflg);

        cout << endl;
        cout << " Fit done. Fetching Best-fit values..." << endl;
        cout << endl;

        double coup[4][4], err;
        minuit.GetParameter(3, coup[1][1], err);
        minuit.GetParameter(4, coup[2][1], err);
        minuit.GetParameter(5, coup[2][2], err);
        minuit.GetParameter(6, coup[3][1], err);
        minuit.GetParameter(7, coup[3][2], err);
        minuit.GetParameter(8, coup[3][3], err);

        double alph[3], err1, err2, err3;
        minuit.GetParameter(0, alph[0], err1);
        minuit.GetParameter(1, alph[1], err2);
        minuit.GetParameter(2, alph[2], err3);

        cout << " BEST FIT TRAJECTORY PARAMETERS: "<< endl;
        cout << setw(10) << "alpha_0" << setw(20) << alph[0] << endl;
        cout << setw(10) << "alpha_p" << setw(20) << alph[1] << endl;
        cout << setw(10) << "width" << setw(20) << alph[2] << endl;
        cout << endl;

        cout << " BEST FIT COUPLINGS: " << endl;
        cout << setw(10) << "n" << setw(10) << "i" << setw(10) << "a_n,i" << endl;
        for (int i = 1; i < 4; i++)
        {
                for (int j = 1; j <= i; j++)
                {
                        cout << setw(10) << i << setw(10) << j << setw(10) << coup[i][j] << endl;
                }
        }
        cout << endl;
        OUTPUT += "-COUPLINGS.dat";
        cout << " Printing results.." << endl;

        ofstream couple;
        couple.open(OUTPUT.c_str());
        couple << left << setw(20) << alph[0] << setw(20) << alph[1] << setw(20) << alph[2] << endl;
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
