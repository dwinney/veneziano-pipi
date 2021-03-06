
// Contains Secondary functions used when plotting amplitudes, uses ROOT to plot
// (formerly gnuplot). Funtionality for GKPY and VENEZ models.
// Also can plot a generic amplitude (three columns: s, real, imag) from FILE.
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "pipi.h"

#include <TGraph.h>
#include <TMultiGraph.h>
#include <TCanvas.h>
#include <TApplication.h>
#include <TLegend.h>
#include <TError.h>

//If MODEL = GKPY
void plotGKPY(int MODE, int plot, string OPTION, string OUTPUT)
{
        setprecision(17);
        double s[10000], delta[10000], eta[10000], sigma[10000], amp_re[10000], amp_im[10000];
        double smax, step;

        complex<double> amp;

        smax = pow(1.42, 2.);


        step = ((smax - sthPi) / (double) DPOINTS );

        int wave = -1;
        switch (MODE) {
        case 1: {
                wave = atoi(OPTION.c_str());
                cout << " Printing forward GKPY isospin-" << OPTION << " amplitude..." << endl;
                break;
        }
        case 2: {
                cout << " Printing total forward amplitude from GKPY..." << endl;
                break;
        }
        case 3: {
                wave = WaveTranslate(OPTION);
                cout << " Printing phaseshift and inelasticities for " << OPTION << " wave..." << endl;
                break;
        }
        case 4: {
                wave = WaveTranslate(OPTION);
                cout << " Printing amplitude for GKPY " << OPTION << " partial-wave..." << endl;
                break;
        }
        }
        string DAT = OUTPUT + ".dat";
        ofstream gkpry;
        gkpry.open(DAT.c_str());
        for (int i = 0; i < DPOINTS; i++)
        {
                s[i] = sthPi + 1e-5 + step * double(i); //offset threshold
                gkpry << left << setw(30) << sqrt(s[i]);

                switch (MODE) {
                case 1:  //isospin amplitudes
                {
                        amp = GKPRY_iso_amp(wave, s[i], 1.);
                        amp_re[i] = real(amp); amp_im[i] = imag(amp);
                        gkpry << setw(30) << amp_re[i] << setw(30) << amp_im[i] << endl;
                        break;
                }
                case 2: //total amplitude
                {
                        amp = GKPRY_amplitude(s[i], 1.);
                        amp_re[i] = real(amp); amp_im[i] = imag(amp);
                        gkpry << setw(30) << amp_re[i] << setw(30) << amp_im[i] << endl;
                        break;
                }
                case 3: //phaseshifts and inelasticities
                {
                        int iso = wave % 10;
                        int l = wave / 10;
                        delta[i] = phase_shift(l, iso, s[i]);
                        eta[i] = inelasticity(l, iso, s[i]);
                        gkpry << setw(30) << delta[i] << setw(30) << eta[i] << endl;
                        break;
                }
                case 4: //partial waves
                {
                        int iso = wave % 10;
                        int l = wave / 10;
                        amp = GKPRY_partial_wave(l, iso, s[i]);
                        amp_re[i] = real(amp); amp_im[i] = imag(amp);
                        gkpry << setw(30) << amp_re[i] << setw(30) << amp_im[i] << endl;
                        break;
                }
                }
        }

        gkpry.close();
        cout << " Output to " << DAT << endl;
        cout << endl;

        if (plot > 0)
        {
                cout << " Plotting..." << endl;
                if (MODE != 3)
                {
                        gErrorIgnoreLevel = kWarning;
                        string PDF = OUTPUT + "-PLOT.pdf";
                        TCanvas * c1 = new TCanvas("c1",":)",200,10,700,500);
                        TMultiGraph *mg = new TMultiGraph();
                        mg->SetTitle("");

                        TGraph * real = new TGraph(DPOINTS, s, amp_re);
                        real->SetLineColor(2);
                        real->SetName("real");
                        real->SetLineWidth(4);

                        TGraph * imag = new TGraph(DPOINTS, s, amp_im);
                        imag->SetLineColor(4);
                        imag->SetName("imag");
                        imag->SetLineWidth(4);

                        mg->Add(real);
                        mg->Add(imag);
                        mg->Draw("AC  ");

                        TLegend* l = new TLegend(.1,0.77,0.27,0.9);
                        l->AddEntry("real","Real","l");
                        l->AddEntry("imag","Imaginary","l");
                        l->Draw();
                        c1->Print(PDF.c_str());

                        c1->Modified();
                        c1->Update();
                        cout << " Output to " << PDF << endl;
                        cout << endl;
                }
                else if (MODE == 3) //Special plotting instructions for inelastivities and phaseshifts
                {
                        gErrorIgnoreLevel = kWarning;
                        string PDF1 = "./output/GKPY-inelasticity-" + OPTION + "-PLOT.pdf";
                        TCanvas * c1 = new TCanvas("c1",":)",200,10,700,500);

                        TGraph *g1 = new TGraph(DPOINTS, s, eta);
                        g1->SetLineColor(2);
                        string title = "Inelasticity (" + OPTION + "-wave)";
                        g1->SetTitle(title.c_str());
                        g1->SetLineWidth(4);

                        g1->Draw("AC");
                        c1->Print(PDF1.c_str());

                        c1->Modified();
                        c1->Update();
                        cout << " Output to " << PDF1 << endl;

                        c1->Clear("g1");
                        string PDF2 = "./output/GKPY-phaseshift-" + OPTION + "-PLOT.pdf";
                        TGraph *g2 = new TGraph(DPOINTS, s, delta);
                        g2->SetLineColor(2);
                        title = "Phase Shift (" + OPTION + "-wave)";
                        g2->SetTitle(title.c_str());
                        g2->SetLineWidth(4);

                        g2->Draw("AC");
                        c1->Print(PDF2.c_str());

                        c1->Modified();
                        c1->Update();
                        cout << " and " << PDF2 << endl;
                        cout << endl;

                }
        }
}

// Imports couplings from a given string filename
void getCOUPLING(string INPUT, double ** output)
{
        double n[30], i[30];
        double aNI[30];
        ifstream coupling;
        coupling.open(INPUT.c_str());
        if (coupling.fail())
        {
                cout << "Could not open " << INPUT << ". Quitting..." << endl;
                cout << endl;
                exit(1);
        }
        else
        {
                int j = 0;
                int ii, nn;
                while(!coupling.eof())
                {
                        coupling >> n[j];
                        coupling >> i[j];
                        coupling >> aNI[j];

                        if (j > 0 && i[j] < maxN +  1 )
                        {
                                nn = (int) n[j]; ii = (int) i[j];         //switch from doubles to int for n and i
                                output[nn][ii] = aNI[j];
                        }
                        j++;
                }
                output[0][0] = n[0]; output[0][1] = i[0]; output [0][2] = aNI[0];
                // cout << output[0][0] << endl;
                cout << " Importing couplings and Regge Trajectory parameters from: " << endl;
                cout << " " << INPUT << endl;
                cout << endl;
        }
        coupling.close();

}

//VENEZIANO model
void plotVENEZ(int MODE, int plot, string OPTION, string INPUT, string OUTPUT)
{
        setprecision(17);
        double s[10000], amp_re[10000], amp_im[10000];
        double smax, step;
        complex<double> amp;

        smax = pow(1.42, 2.);
        step = .001;

        step = ((smax - sthPi) / double(DPOINTS));

        double **a_matrix;
        a_matrix = new double *[maxN+1];
        for (int i = 0; i < maxN+1; i++) a_matrix[i] = new double[maxN+1];

        double alph[3];
        getCOUPLING(INPUT, a_matrix);
        alph[0] = a_matrix[0][0]; alph[1] = a_matrix[0][1]; alph[2] = a_matrix[0][2];

        int wave = -1;
        switch (MODE) {
        case 1: {
                wave = atoi(OPTION.c_str());
                cout << " Printing forward VENEZ isospin-" << OPTION << " amplitude..." << endl;
                break;
        }
        case 2: {
                cout << " Printing total forward amplitude from VENEZ..." << endl;
                break;
        }
        case 3: {
                wave = WaveTranslate(OPTION);
                cout << " Printing amplitude for VENEZ " << OPTION << " partial-wave..." << endl;
                break;
        }
        }

        ofstream venez;
        string DAT = OUTPUT + ".dat";
        venez.open(DAT.c_str());
        for (int i = 0; i < DPOINTS; i++)
        {
                s[i] = sthPi + 1e-5 + step * double(i); //offset threshold
                venez << left << setw(30) << sqrt(s[i]);

                switch (MODE) {
                case 1: //isospin amplitudes
                {
                        amp = VENEZ_iso_amp(wave, a_matrix, alph, s[i], 1.);
                        amp_re[i] = real(amp); amp_im[i] = imag(amp);
                        venez << setw(30) << amp_re[i] << setw(30) << amp_im[i] << endl;
                        break;
                }
                case 2: //total amplitude
                {
                        amp = VENEZ_amplitude(a_matrix, alph, s[i], 1.);
                        amp_re[i] = real(amp); amp_im[i] = imag(amp);
                        venez << setw(30) << amp_re[i] << setw(30) << amp_im[i] << endl;
                        break;
                }
                case 3: //partial wave
                {
                        int iso = wave % 10;
                        int l = wave / 10;
                        amp = VENEZ_partial_wave(l, iso, a_matrix, alph, s[i]);
                        amp_re[i] = real(amp); amp_im[i] = imag(amp);
                        venez << setw(30) << amp_re[i] << setw(30) << amp_im[i] << endl;
                        break;
                }
                }
        }

        cout << " Output to " << DAT << endl;
        cout << endl;

        if (plot > 0)
        {
                cout << " Plotting..." << endl;

                gErrorIgnoreLevel = kWarning;
                string PDF = OUTPUT + "-PLOT.pdf";
                TCanvas * c1 = new TCanvas("c1",":)",200,10,700,500);
                TMultiGraph *mg = new TMultiGraph();
                mg->SetTitle("");

                TGraph * real = new TGraph(DPOINTS, s, amp_re);
                real->SetLineColor(2);
                real->SetName("real");
                real->SetLineWidth(4);

                TGraph * imag = new TGraph(DPOINTS, s, amp_im);
                imag->SetLineColor(4);
                imag->SetName("imag");
                imag->SetLineWidth(4);

                mg->Add(real);
                mg->Add(imag);
                mg->Draw("AC");

                TLegend* l = new TLegend(.1,0.77,0.27,0.9);
                l->AddEntry("real","Real","l");
                l->AddEntry("imag","Imaginary","l");
                l->Draw();
                c1->Print(PDF.c_str());

                c1->Modified();
                c1->Update();
                cout << " Output to " << PDF << endl;
                cout << endl;

        }
}

//Plot from FILE. Requries three columns for an amplitude.
void plotFILE(string INPUT, string OUTPUT)
{
        setprecision(17);
        double s[10000], amp_re[10000], amp_im[10000];
        int i = 0;

        ifstream infile;
        infile.open(INPUT.c_str());
        if (infile.fail())
        {
                cout << " Could not open " << INPUT << ". Quitting..." << endl;
                cout << endl;
                exit(1);
        }
        else
        {
                cout << " Reading data from file:" << endl;
                cout << " " << INPUT << endl;
                cout << endl;

                while(!infile.eof())
                {
                        infile >> s[i];
                        infile >> amp_re[i];
                        infile >> amp_im[i];
                        i++;
                }
                i -= 1;
        }
        infile.close();

        cout << " Plotting..." << endl;

        gErrorIgnoreLevel = kWarning;

        TCanvas * c1 = new TCanvas("c1",":)",200,10,700,500);
        TMultiGraph * mg = new TMultiGraph();

        TGraph * real = new TGraph(i, s, amp_re);
        real->SetLineColor(2);
        real->SetName("real");
        real->SetLineWidth(4);

        TGraph * imag = new TGraph(i, s, amp_im);
        imag->SetLineColor(4);
        imag->SetName("imag");
        imag->SetLineWidth(4);

        mg->Add(real);
        mg->Add(imag);
        mg->Draw("AL");

        // TLegend* l = new TLegend(.1,0.77,0.27,0.9);
        // l->AddEntry("real","Real","l");
        // l->AddEntry("imag","Imaginary","l");
        // l->Draw();
        c1->Print(OUTPUT.c_str());

        // c1->Modified();
        // c1->Update();

        cout << " Output to " << OUTPUT << endl;
        cout << endl;
}
