
// Main program, parses user inputs to fit and/or plot functions.
//
// Author:       Daniel Winney (2018)
// Affiliation:  Joint Physics Analysis Center (JPAC)
// Email:        dwinney@iu.edu
// ---------------------------------------------------------------------------

#include "pipi.h"

int main(int argc, char* argv[])
{
        cout << " " << endl; // Line break to make terminal look pretty :)
        int plot = -1; int fit = -1;
        string MODEL, AMP, OPTION, OPTION2;

        for (int ii = 0; ii < argc; ii++)
        {
                if (strcmp(argv[ii], "plot")==0)
                {
                        if (argc <  4)
                        {
                                cout << "Not enough arguments for plot, check README. Quitting..." << endl;
                                cout << endl;
                                exit(1);
                        }
                        plot = 1;
                        MODEL = argv[ii+1]; //GKPY or VENEZ
                        AMP = argv[ii+2]; //isospin, partial, total
                        if (argc > 4)
                        {
                                OPTION = argv[ii+3]; //wave
                        }
                        if (argc > 5)
                        {
                                OPTION2 = argv[ii+4]; //coupling inputfilename
                        }
                        break;
                }
                if (strcmp(argv[ii], "fit")==0)
                {
                        fit = 1;
                        if (argc < 5)
                        {
                                cout << " Not enough arguments for fit, check README. Quiting..." << endl;
                                cout << endl;
                                exit(1);
                        }
                        else
                        {
                                MODEL = argv[ii+1]; // VENEZ
                                AMP = argv[ii+2];
                                OPTION = argv[ii+3]; //wave
                                OPTION2 = argv[ii+4];   //file to fit
                        }
                }
        }

        //PLOTTING
        if (plot > 0)
        {
                string OUTPUT;
                OUTPUT = "./output/" + MODEL + "-" + AMP + "-" + OPTION;
                int CASE = -400;
                if (MODEL == "GKPY")
                {

                        if (AMP == "isospin")
                        {
                                CASE = 1;
                        }
                        else if (AMP == "total")
                        {
                                CASE = 2;
                        }
                        else if (AMP == "phaseshift")
                        {
                                CASE = 3;
                        }
                        else if (AMP == "partial")
                        {
                                CASE = 4;
                        }
                        else
                        {
                                cout << " Invalid option for GKPY parameterization. Quitting..." << endl;
                                cout << endl;
                                exit(1);
                        }
                        plotGKPY(CASE, plot, OPTION, OUTPUT);
                }
                else if (MODEL == "VENEZ")
                {
                        if (AMP == "isospin")
                        {
                                CASE = 1;
                        }
                        if (AMP == "total")
                        {
                                CASE = 2;
                                OPTION2 = OPTION;
                                OUTPUT = "./output/" + MODEL + "-" + AMP;
                        }
                        else if (AMP == "partial")
                        {
                                CASE = 3;
                        }
                        if (CASE < 0)
                        {
                                cout << " Invalid option for VENEZ parameterization. Quitting..." << endl;
                                cout << endl;
                                exit(1);
                        }
                        plotVENEZ(CASE, plot, OPTION, OPTION2, OUTPUT);

                }
                else if (MODEL == "FILE")
                {
                        if (OPTION == "") OPTION = "./output/FILE-DEFAULT-NAME-PLOT.pdf";
                        else OPTION = "./output/" + OPTION + ".pdf";
                        plotFILE(AMP, OPTION);
                }
                else
                {
                        cout << "Invalid MODEL chosen. Quiting..." << endl;
                        cout << endl;
                        exit(1);
                }
        }
        else if (fit > 0)
        {
                if (MODEL == "VENEZ")
                {
                        string OUTPUT;
                        OUTPUT = "./output/" + MODEL + "-" + AMP + "-" + OPTION;
                        int MODE = -600;
                        if (AMP == "partial")
                        {
                                MODE = 1;
                        }
                        else if (AMP == "isospin")
                        {
                                MODE = 2;
                        }
                        fitVENEZ(MODE, OPTION, OPTION2, OUTPUT);
                }
                else
                {
                        cout << "Invalid MODEL chosen. Quiting..." << endl;
                        cout << endl;
                        exit(1);
                }
        }
        else
        {
                cout << "Invalid commands entered. See README. Quitting..." << endl;
                cout << endl;
                exit(1);
        }
        return 6;
}
