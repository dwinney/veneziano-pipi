#include "veneziano.h"

#include <TFitter.h>
#include <TMinuit.h>
#include <TROOT.h>
#include <TMath.h>

double s_dat[1000], re_dat[1000], im_dat[1000];
int ISOCHOICE;

//Function to import amplitude data to be fit from file.
//Three columns s values, real part, imaginary part
//TODO: have general fit file name
void get_Data()
{
        string filename = "./output/data.dat";
        ifstream infile(filename.c_str());
        if (infile.fail())
        {
                cout << " Couldn't open data file. Quiting..." << endl;
                exit(1);
        }
        else{
                int n = 0;
                while (!infile.eof())
                {
                        infile >> s_dat[n];
                        infile >> re_dat[n];
                        infile >> im_dat[n];
                        n++;
                }
                cout << " Imported data from: " << filename << endl;

        }
}

//chi-square of the amplitude of given isospin to data.
//fits both real and imaginary parts simultaneously.
double chi_square(double coup[][maxN+1], double alph[])
{
        double chi = 0.;
        double s, error;
        complex<double> VENEZ;

        error = 10.0;

        for (int i = 0; i < 100; i++)
        {
                s = s_dat[i];
                VENEZ = VENEZ_iso_amp(ISOCHOICE, coup, alph, s, 1.);
                chi += pow(((re_dat[i] - real(VENEZ)) / error), 2.);
                chi += pow(((im_dat[i] - imag(VENEZ)) / error), 2.);
        }
        return chi;
}

//Wrapper function to
static void WRAPPER(int& npar, double* g, double& result, double *par, int flag)
{
        double coup[4][4];
        coup[1][1] = par[3];
        coup[2][1] = par[4];
        coup[2][2] = par[5];
        coup[3][1] = par[6];
        coup[3][2] = par[7];
        coup[3][3] = par[8];

        double alph [3];
        alph[0] = par[0];
        alph[1] = par[1];
        alph[2] = par[2];

        result = chi_square(coup, alph);
        // cout << result << endl;
}

int main(int argc, char* argv[])
{
        cout << endl;
        setprecision(10);
        int np = 100;
        int iso = 600;
        int plot = -1;
        for (int ii = 0; ii < argc; ii++)
        {
                if (strcmp(argv[ii],"-i")==0) iso = atof(argv[ii+1]);
                if (strcmp(argv[ii], "-plot")==0) plot = 1;
                // if (strcmp(argv[ii], "-n")==0) np = atoi(argv[ii+1]);
        }

        if (iso > 2)
        { cout << " Input valid Isospin to fit [-i isospin]. Quitting..." << endl;
          exit(1);}

        ISOCHOICE = iso;
        double s;
        double smax, step;
        complex<double> amp;

        smax = pow(1.42, 2.);
        step = (smax - sthPi) / double(np);         //get stepsize

        cout << "-------------------------------------------------------------" << endl;
        cout << " Printing Isospin-" << ISOCHOICE << " GKPY amplitude to ./output/data.dat" << endl;
        cout << " Using " << np << " data points. Step size = " << step << " GeV^2" << endl;

        ofstream data;
        data.open("./output/data.dat");
        for (int i = 0; i < np; i++)
        {
                s = sthPi + 1e-8 + step*double(i);
                amp = GKPRY_iso_amp(ISOCHOICE, s, 1.);
                data << left << setw(30) << s << setw(30) << real(amp) << setw(30) << imag(amp) << endl;
        }
        data.close();

        get_Data();
        cout << "-------------------------------------------------------------" << endl;
        cout << endl;
        cout << "-------------------------------------------------------------" << endl;
        TMinuit minuit(9);
        minuit.SetFCN(WRAPPER);

        double start[3], stepsiz[4], min[4], max[4];

        //alpha_0
        start[0] = .2; stepsiz[0] = 0.1;
        min[0] = 0.; max[0] = 1.5;
        //alpha^p
        start[1] = 0.75; stepsiz[1] = 0.1;
        min[1] = 0.5; max[1] = 1.5;
        //width
        start[2] = .164; stepsiz[2] = .05;
        min[2] = 0.0; max[2] = 1.;
        //couplings
        start[3] = 0.; stepsiz[3] = 2.;
        min[3] = -0.; max[3] = -0.;


        double arglist[10];
        int ierflg = 0;
        minuit.mnparm(0, "alpha_0", start[0], stepsiz[0], min[0], max[0], ierflg);
        minuit.mnparm(1, "alpha_p", start[1], stepsiz[1], min[1], max[1], ierflg);
        minuit.mnparm(2, "width", start[2], stepsiz[2], min[2], max[2], ierflg);
        minuit.mnparm(3, "a_1,1", start[3], stepsiz[3], min[3], max[3], ierflg);
        minuit.mnparm(4, "a_2,1", start[3], stepsiz[3], min[3], max[3], ierflg);
        minuit.mnparm(5, "a_2,2", start[3], stepsiz[3], min[3], max[3], ierflg);
        minuit.mnparm(6, "a_3,1", start[3], stepsiz[3], min[3], max[3], ierflg);
        minuit.mnparm(7, "a_3,2", start[3], stepsiz[3], min[3], max[3], ierflg);
        minuit.mnparm(8, "a_3,3", start[3], stepsiz[3], min[3], max[3], ierflg);

        arglist[0]=-1;
        minuit.mnexcm("SET PRINTOUT", arglist, 1, ierflg);
        minuit.mnexcm("SET NOW", arglist, 1, ierflg);
        arglist[0]=2;
        minuit.mnexcm("SET STR", arglist, 1, ierflg);
        minuit.SetErrorDef(1); //1 for chi square

        // PRESERVE BOSE SYMMETRY
        if (ISOCHOICE == 0 || ISOCHOICE == 2)
        {
                minuit.FixParameter(3);
                minuit.FixParameter(4);
                minuit.FixParameter(6);
                minuit.FixParameter(8);
        }
        else if (ISOCHOICE == 1)
        {
                minuit.FixParameter(5);
                minuit.FixParameter(7);
        }

        arglist[0] = 1000;
        minuit.mnexcm("SIMPLEX", arglist,1,ierflg);
        // minuit.mnexcm("MINIMIZE", arglist, 1, ierflg);
        minuit.mnexcm("MIGRAD", arglist, 1, ierflg);
        // minuit.mnexcm("MINOS", arglist, 1, ierflg);

        cout << "-------------------------------------------------------------" << endl;
        cout << endl;
        cout << "-------------------------------------------------------------" << endl;
        minuit.mnprin(1, arglist[0]);
        cout << "-------------------------------------------------------------" << endl;
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

        cout << "-------------------------------------------------------------" << endl;
        cout << " BEST FIT TRAJECTORY PARAMETERS: "<< endl;
        cout << setw(10) << "alpha_0" << setw(20) << alph[0] << endl;
        cout << setw(10) << "alpha_p" << setw(20) << alph[1] << endl;
        cout << setw(10) << "width" << setw(20) << alph[2] << endl;
        cout << "-------------------------------------------------------------" << endl;
        cout << endl;
        cout << "-------------------------------------------------------------" << endl;
        cout << " Printing coupling constant matrix to ./output/couplings.dat" << endl;
        cout << " Plotting best fit Isospin-" << ISOCHOICE << " amplitude to ./output/fit.data" << endl;
        cout << "-------------------------------------------------------------" << endl;

        ofstream fit;
        fit.open("./output/fit.dat");
        for (int i = 0; i < np; i++)
        {
                s = sthPi + step*double(i);
                cd amp = VENEZ_iso_amp(ISOCHOICE, coup, alph, s, 1.);
                fit << left << setw(30) << s << setw(30) << real(amp)  << setw(30) << imag(amp) << endl;
        }
        fit.close();

        ofstream couple;
        couple.open("./output/couplings.dat");
        couple << left << setw(20) << alph[0] << setw(20) << alph[1] << setw(20) << alph[2] << endl;
        for (int i = 1; i < 4; i++)
        {
                for (int j = 1; j < i+1; j++)
                {
                        couple << left << setw(20) << i << setw(20) << j <<setw(20) << coup[i][j] << endl;
                }
        }
        couple.close();

        if (plot > 0)
        {
                system("gnuplot ./src/gnuplot/graph.gnu");
                system("okular ./output/fit.pdf");
        }

        return 0;
}
